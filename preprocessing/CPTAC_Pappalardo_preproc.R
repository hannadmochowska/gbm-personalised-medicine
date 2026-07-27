# CPTAC-GBM Pappalardo (RTK/RAS/MAPK + PI3K/AKT/mTOR) ODE Preprocessing
#
# Purpose: Extract RTK/PI3K/MAPK pathway gene expression for CPTAC-GBM
# patients, in the same relative-expression units used for the TCGA-GBM run
# (GBM_Pappalardo_preproc.R). Like CPTAC_RbE2F_preproc.R, this is NOT for
# re-running the Pappalardo ODE model on CPTAC (TCGA discovery was a clean
# null) — it's to give the Q1 cross-dataset signature access to these 11
# genes' expression as candidate predictors in both cohorts.
#
# Genes (11, same as GBM_Pappalardo_preproc.R):
#   EGFR, PTEN, SOS1, KRAS, RAF1, MAP2K1, MAPK1, PIK3CA, AKT1, MTOR, RPS6KB1
#
# Output:
#   CPTAC_patient_Pappalardogenes.csv — one row per CPTAC-GBM patient
#
# Run this before Q1_build_signature_TCGA.R / Q1_validate_signature_CPTAC.R.

# Packages
library(tidyverse)
library(stringr)

# 1. File paths
cptac_rna_file  <- "gbm_cptac_2021/data_mrna_seq_fpkm_zscores_ref_all_samples.txt"
cptac_clin_file <- "gbm_cptac_2021/data_clinical_patient.txt"

PAPPALARDO_GENES <- c("EGFR", "PTEN", "SOS1", "KRAS", "RAF1", "MAP2K1",
                       "MAPK1", "PIK3CA", "AKT1", "MTOR", "RPS6KB1")

# Helper: identical logic to CPTAC_p53_preproc.R / CPTAC_RbE2F_preproc.R
to_relative <- function(mat) {
  mat[mat < 0] <- 0
  col_means <- colMeans(mat, na.rm = TRUE)
  col_means[col_means == 0] <- 1
  sweep(mat, 2, col_means, "/")
}

# 2. CPTAC-GBM clinical data — derive survival columns
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
rna_raw <- rna_raw[!is.na(rna_raw$Hugo_Symbol) & rna_raw$Hugo_Symbol != "", ]
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

# CPTAC RNA column names ARE the patient IDs already
patient_ids <- rownames(rna_t)

# Keep only genes present in PAPPALARDO_GENES
available_genes <- intersect(PAPPALARDO_GENES, colnames(rna_t))
missing_genes    <- setdiff(PAPPALARDO_GENES, colnames(rna_t))
if (length(missing_genes) > 0)
  cat("Warning: genes not found in CPTAC data:", paste(missing_genes, collapse = ", "), "\n")

patient_df <- rna_t[, available_genes, drop = FALSE]
for (g in missing_genes) patient_df[[g]] <- 1.0
patient_df <- patient_df[, PAPPALARDO_GENES]

# Keep only patients with survival info
keep_ids <- patient_ids %in% clin_raw$PATIENT_ID
patient_df <- patient_df[keep_ids, ]
patient_df$SAMPLE_ID  <- patient_ids[keep_ids]
patient_df$PATIENT_ID <- patient_ids[keep_ids]

cat(sprintf("CPTAC-GBM patients with RNA + clinical data: %d\n", nrow(patient_df)))

# 4. Normalize (column-wise mean so reference factor = 1)
gene_cols <- PAPPALARDO_GENES[PAPPALARDO_GENES %in% colnames(patient_df)]
patient_df[, gene_cols] <- to_relative(as.matrix(patient_df[, gene_cols]))

write.csv(patient_df, "CPTAC_patient_Pappalardogenes.csv", row.names = FALSE)
cat("Saved: CPTAC_patient_Pappalardogenes.csv\n")
