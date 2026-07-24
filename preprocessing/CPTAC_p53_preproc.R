# CPTAC-GBM p53 ODE Preprocessing
#
# Purpose: Extract p53 pathway gene expression for CPTAC-GBM patients,
# normalise to relative-expression units that the p53 ODE model (CPTAC_p53_model.ipynb) expects.
#
# This is the brief's Q1 cohort: an independent GBM cohort (Wang et al., Cell 2021,
# gbm_cptac_2021 on cBioPortal) never used in Assignment 1 or 2.
#
# Output:
#   CPTAC_patient_p53genes.csv — one row per CPTAC-GBM patient, used by CPTAC_p53_model.ipynb
#
# Run this BEFORE CPTAC_p53_model.ipynb.

# Packages
library(tidyverse)
library(stringr)

# 1. File paths
cptac_rna_file  <- "gbm_cptac_2021/data_mrna_seq_fpkm_zscores_ref_all_samples.txt"
cptac_clin_file <- "gbm_cptac_2021/data_clinical_patient.txt"

# Eight p53 pathway genes used as ODE initial conditions (same as TCGA pipeline)
P53_GENES <- c("ATM", "CHEK2", "HIPK2", "MDM2", "PPM1D", "SIAH1", "TP53", "WSB1")

# Helper: convert expression to ODE-ready relative-expression values
# (identical logic to GBM_p53_preproc.R)
to_relative <- function(mat) {
  mat[mat < 0] <- 0
  col_means <- colMeans(mat, na.rm = TRUE)
  col_means[col_means == 0] <- 1
  sweep(mat, 2, col_means, "/")
}

# 2. CPTAC-GBM clinical data — derive survival columns
# CPTAC has no ready-made OS_MONTHS/OS_STATUS like TCGA; derive from
# VITAL_STATUS ("Deceased"/"Living") + PATH_DIAG_TO_DEATH_DAYS / PATH_DIAG_TO_LAST_CONTACT_DAYS.
cat("Loading CPTAC-GBM clinical data …\n")

clin_raw <- read.delim(cptac_clin_file, sep = "\t", skip = 4, stringsAsFactors = FALSE)
clin_raw$PATIENT_ID <- trimws(clin_raw$PATIENT_ID)

clin_raw$SURV_STATUS <- as.integer(clin_raw$VITAL_STATUS == "Deceased")
clin_raw$SURV_TIME   <- ifelse(
  clin_raw$SURV_STATUS == 1,
  suppressWarnings(as.numeric(clin_raw$PATH_DIAG_TO_DEATH_DAYS)) / 30.44,
  suppressWarnings(as.numeric(clin_raw$PATH_DIAG_TO_LAST_CONTACT_DAYS)) / 30.44
)
clin_raw <- clin_raw[!is.na(clin_raw$SURV_TIME) & clin_raw$SURV_TIME > 0, ]
cat(sprintf("CPTAC-GBM patients with derived survival data: %d\n", nrow(clin_raw)))

# 3. CPTAC-GBM RNA-seq (z-scores)
cat("Loading CPTAC-GBM RNA-seq z-scores …\n")

rna_raw <- read.delim(cptac_rna_file, check.names = FALSE, stringsAsFactors = FALSE)
rna_raw <- rna_raw[rna_raw$Hugo_Symbol != "", ]
rna_raw <- rna_raw[!duplicated(rna_raw$Hugo_Symbol), ]
rownames(rna_raw) <- rna_raw$Hugo_Symbol
rna_raw <- rna_raw[, !colnames(rna_raw) %in% c("Hugo_Symbol", "Entrez_Gene_Id")]

# Transpose: samples as rows
rna_t <- as.data.frame(t(rna_raw))
rna_t[is.na(rna_t)] <- 0

# z-scores -> relative expression: shift by +1 so neutral value is 1, then clip
orig_rownames <- rownames(rna_t)
rna_t <- as.data.frame(
  lapply(rna_t, function(x) pmax(0.01, as.numeric(x) + 1))
)
rownames(rna_t) <- orig_rownames

# CPTAC RNA column names ARE the patient IDs already (e.g. "C3L-00104") — no trimming needed
patient_ids <- rownames(rna_t)

# Keep only genes present in P53_GENES
available_genes <- intersect(P53_GENES, colnames(rna_t))
missing_genes   <- setdiff(P53_GENES, colnames(rna_t))
if (length(missing_genes) > 0)
  cat("Warning: genes not found in CPTAC data:", paste(missing_genes, collapse = ", "), "\n")

p53_patients <- rna_t[, available_genes, drop = FALSE]
for (g in missing_genes) p53_patients[[g]] <- 1.0
p53_patients <- p53_patients[, P53_GENES]

# Keep only patients with survival info
keep_ids <- patient_ids %in% clin_raw$PATIENT_ID
p53_patients <- p53_patients[keep_ids, ]
p53_patients$SAMPLE_ID  <- patient_ids[keep_ids]
p53_patients$PATIENT_ID <- patient_ids[keep_ids]

cat(sprintf("CPTAC-GBM patients with RNA + clinical data: %d\n", nrow(p53_patients)))

# 4. Normalize (column-wise mean so reference factor = 1)
gene_cols <- P53_GENES[P53_GENES %in% colnames(p53_patients)]
p53_patients[, gene_cols] <- to_relative(as.matrix(p53_patients[, gene_cols]))

write.csv(p53_patients, "CPTAC_patient_p53genes.csv", row.names = FALSE)
cat("Saved: CPTAC_patient_p53genes.csv\n")
