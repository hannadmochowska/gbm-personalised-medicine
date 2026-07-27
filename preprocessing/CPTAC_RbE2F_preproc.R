# CPTAC-GBM Rb-E2F ODE Preprocessing
#
# Purpose: Extract Rb-E2F restriction-point pathway gene expression for
# CPTAC-GBM patients, in the same relative-expression units used for the
# TCGA-GBM Rb-E2F run (GBM_RbE2F_preproc.R). This is NOT for re-running the
# Rb-E2F ODE model on CPTAC (discovery was a clean null in TCGA, so no
# validation was planned for that model) — it's purely to give the Q1
# cross-dataset signature (Q1_build_signature_TCGA.R /
# Q1_validate_signature_CPTAC.R) access to these 5 genes' expression as
# candidate predictors in both cohorts.
#
# Genes: MYC, CCND1, CCNE1, RB1, E2F1 (same 5 as GBM_RbE2F_preproc.R).
#
# Output:
#   CPTAC_patient_RbE2Fgenes.csv — one row per CPTAC-GBM patient
#
# Run this before Q1_build_signature_TCGA.R / Q1_validate_signature_CPTAC.R.

# Packages
library(tidyverse)
library(stringr)

# 1. File paths
cptac_rna_file  <- "gbm_cptac_2021/data_mrna_seq_fpkm_zscores_ref_all_samples.txt"
cptac_clin_file <- "gbm_cptac_2021/data_clinical_patient.txt"

RBE2F_GENES <- c("MYC", "CCND1", "CCNE1", "RB1", "E2F1")

# Helper: identical logic to CPTAC_p53_preproc.R, for consistency
to_relative <- function(mat) {
  mat[mat < 0] <- 0
  col_means <- colMeans(mat, na.rm = TRUE)
  col_means[col_means == 0] <- 1
  sweep(mat, 2, col_means, "/")
}

# 2. CPTAC-GBM clinical data — derive survival columns (identical to
# CPTAC_p53_preproc.R; CPTAC has no ready-made OS_MONTHS/OS_STATUS)
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
# CPTAC's file has 1 literal "NA" Hugo_Symbol row (not a blank), unlike TCGA's
# file — same guard as CPTAC_p53_preproc.R needs, applied here for consistency.
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

# CPTAC RNA column names ARE the patient IDs already (e.g. "C3L-00104") —
# no trimming needed, same as CPTAC_p53_preproc.R
patient_ids <- rownames(rna_t)

# Keep only genes present in RBE2F_GENES
available_genes <- intersect(RBE2F_GENES, colnames(rna_t))
missing_genes    <- setdiff(RBE2F_GENES, colnames(rna_t))
if (length(missing_genes) > 0)
  cat("Warning: genes not found in CPTAC data:", paste(missing_genes, collapse = ", "), "\n")

rbe2f_patients <- rna_t[, available_genes, drop = FALSE]
for (g in missing_genes) rbe2f_patients[[g]] <- 1.0
rbe2f_patients <- rbe2f_patients[, RBE2F_GENES]

# Keep only patients with survival info
keep_ids <- patient_ids %in% clin_raw$PATIENT_ID
rbe2f_patients <- rbe2f_patients[keep_ids, ]
rbe2f_patients$SAMPLE_ID  <- patient_ids[keep_ids]
rbe2f_patients$PATIENT_ID <- patient_ids[keep_ids]

cat(sprintf("CPTAC-GBM patients with RNA + clinical data: %d\n", nrow(rbe2f_patients)))

# 4. Normalize (column-wise mean so reference factor = 1)
gene_cols <- RBE2F_GENES[RBE2F_GENES %in% colnames(rbe2f_patients)]
rbe2f_patients[, gene_cols] <- to_relative(as.matrix(rbe2f_patients[, gene_cols]))

write.csv(rbe2f_patients, "CPTAC_patient_RbE2Fgenes.csv", row.names = FALSE)
cat("Saved: CPTAC_patient_RbE2Fgenes.csv\n")
