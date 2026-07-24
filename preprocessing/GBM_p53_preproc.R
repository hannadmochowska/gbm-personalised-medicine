# GBM p53 ODE Preprocessing
#
# Purpose: Extract p53 pathway gene expression for TCGA-GBM patients,
# normalise to relative-expression units that the p53 ODE model (GBM_p53_model.ipynb) expects.
#
# Output:
#   GBM_patient_p53genes.csv — one row per TCGA-GBM patient, used by GBM_p53_model.ipynb
#
# Run this BEFORE GBM_p53_model.ipynb.
#
# Note on normalisation: TCGA RNA-seq data are z-scores; each value is shifted +1 before
# dividing by column mean, ensuring the relative factor is ≈ 1 at the cohort average
# (equivalent to the class approach of normalising absolute RPKM values by column mean).

# Packages
library(tidyverse)
library(stringr)

# 1. File paths
tcga_rna_file  <- "gbm_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
tcga_clin_file <- "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt"

# Eight p53 pathway genes used as ODE initial conditions
P53_GENES <- c("ATM", "CHEK2", "HIPK2", "MDM2", "PPM1D", "SIAH1", "TP53", "WSB1")

# Helper: convert expression to ODE-ready relative-expression values
# The ODE model multiplies reference initial conditions by patient/cell-line specific factors.
# We want factor ≈ 1 for an average sample, > 1 for over-expression, < 1 for under-expression, never exactly 0.
to_relative <- function(mat) {
  # mat: samples × genes matrix (numeric)
  # First clip negatives, then divide each gene by its column mean
  mat[mat < 0] <- 0
  col_means <- colMeans(mat, na.rm = TRUE)
  # avoid div-by-zero if a gene has mean 0 in this cohort
  col_means[col_means == 0] <- 1
  sweep(mat, 2, col_means, "/")
}

# 2. TCGA-GBM patients
cat("Loading TCGA-GBM RNA-seq z-scores …\n")

# Skip the 4-line cBioPortal header
clin_raw <- read.delim(tcga_clin_file, sep = "\t", skip = 4,
                       stringsAsFactors = FALSE)

rna_raw <- read.delim(tcga_rna_file, check.names = FALSE,
                      stringsAsFactors = FALSE)
rna_raw <- rna_raw[rna_raw$Hugo_Symbol != "", ]
rna_raw <- rna_raw[!duplicated(rna_raw$Hugo_Symbol), ]
rownames(rna_raw) <- rna_raw$Hugo_Symbol
rna_raw <- rna_raw[, !colnames(rna_raw) %in% c("Hugo_Symbol", "Entrez_Gene_Id")]

# Transpose: samples as rows
rna_t <- as.data.frame(t(rna_raw))
rna_t[is.na(rna_t)] <- 0

# TCGA z-scores: convert to relative expression
# z-score ≈ 0 → factor 1; z > 0 → above average; z < 0 → below average
# shift by +1 so the neutral value is 1 (then clip)
orig_rownames <- rownames(rna_t)
rna_t <- as.data.frame(
  lapply(rna_t, function(x) pmax(0.01, as.numeric(x) + 1))
)
rownames(rna_t) <- orig_rownames

# Trim patient IDs to 12 chars (TCGA-XX-XXXX)
patient_ids <- str_sub(rownames(rna_t), end = 12)

# Keep only genes present in P53_GENES
available_genes <- intersect(P53_GENES, colnames(rna_t))
missing_genes   <- setdiff(P53_GENES, colnames(rna_t))
if (length(missing_genes) > 0)
  cat("Warning: genes not found in TCGA data:", paste(missing_genes, collapse = ", "), "\n")

p53_patients <- rna_t[, available_genes, drop = FALSE]

# For any missing gene, fill with 1 (=average baseline)
for (g in missing_genes) p53_patients[[g]] <- 1.0
p53_patients <- p53_patients[, P53_GENES]

# Check against clinical data to keep only patients with survival info
clin_raw$PATIENT_ID <- trimws(clin_raw$PATIENT_ID)
clin_raw$OS_MONTHS  <- suppressWarnings(as.numeric(clin_raw$OS_MONTHS))

# Diagnostic: verify IDs are in expected format
cat(sprintf("  Sample RNA IDs (12-char): %s\n",
            paste(head(patient_ids, 3), collapse = ", ")))
cat(sprintf("  Sample clinical IDs:      %s\n",
            paste(head(clin_raw$PATIENT_ID, 3), collapse = ", ")))
cat(sprintf("  Clinical patients with OS data: %d\n",
            sum(!is.na(clin_raw$OS_MONTHS) & clin_raw$OS_MONTHS > 0)))

keep_ids <- patient_ids %in% clin_raw$PATIENT_ID[
  !is.na(clin_raw$OS_MONTHS) & clin_raw$OS_MONTHS > 0]

p53_patients <- p53_patients[keep_ids, ]
p53_patients$SAMPLE_ID  <- rownames(rna_t)[keep_ids]
p53_patients$PATIENT_ID <- patient_ids[keep_ids]

cat(sprintf("TCGA-GBM patients with RNA + clinical data: %d\n", nrow(p53_patients)))

# 3. Normalize (column-wise mean so reference factor = 1)
gene_cols <- P53_GENES[P53_GENES %in% colnames(p53_patients)]
p53_patients[, gene_cols] <- to_relative(as.matrix(p53_patients[, gene_cols]))

write.csv(p53_patients, "GBM_patient_p53genes.csv", row.names = FALSE)
cat("Saved: GBM_patient_p53genes.csv\n")
