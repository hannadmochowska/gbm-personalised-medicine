# GBM Rb-E2F ODE Preprocessing
#
# Purpose: Extract Rb-E2F restriction-point pathway gene expression for TCGA-GBM
# patients, normalise to relative-expression units that the Rb-E2F ODE model
# (GBM_RbE2F_model.ipynb) expects.
#
# Model: Yao et al. 2008, "A bistable Rb-E2F switch underlies the restriction
# point" (Nat Cell Biol), BioModels BIOMD0000000318. 7 species: Myc, E2F, CycD,
# CycE, Rb, Rb-E2F complex, phosphorylated Rb. Verified directly from the SBML
# (species list, reactions, kinetic laws, parameters) before writing this script.
#
# Genes used as ODE synthesis-rate scalers (patient-specific): MYC, CCND1,
# CCNE1, RB1, E2F1 — the 5 species in the model that correspond to a single,
# specific human gene (the model lumps CyclinD/CDK4,6 and CyclinE/CDK2 into
# single catalytic species "CycD"/"CycE", so no separate CDK4/CDK6/CDK2 terms).
#
# Output:
#   GBM_patient_RbE2Fgenes.csv — one row per TCGA-GBM patient, used by GBM_RbE2F_model.ipynb
#
# Run this BEFORE GBM_RbE2F_model.ipynb.

# Packages
library(tidyverse)
library(stringr)

# 1. File paths
tcga_rna_file  <- "gbm_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
tcga_clin_file <- "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt"

# Five Rb-E2F pathway genes used as ODE synthesis-rate scalers
RBE2F_GENES <- c("MYC", "CCND1", "CCNE1", "RB1", "E2F1")

# Helper: convert expression to ODE-ready relative-expression values
# (identical logic to GBM_p53_preproc.R, for consistency across the project)
to_relative <- function(mat) {
  mat[mat < 0] <- 0
  col_means <- colMeans(mat, na.rm = TRUE)
  col_means[col_means == 0] <- 1
  sweep(mat, 2, col_means, "/")
}

# 2. TCGA-GBM patients
cat("Loading TCGA-GBM RNA-seq z-scores …\n")

clin_raw <- read.delim(tcga_clin_file, sep = "\t", skip = 4, stringsAsFactors = FALSE)

rna_raw <- read.delim(tcga_rna_file, check.names = FALSE, stringsAsFactors = FALSE)
# NB: TCGA's file has empty (not NA) unmapped Hugo_Symbol entries — verified
# directly against the raw file (13 blank rows, 0 literal "NA" rows) before
# reusing this exact filter from GBM_p53_preproc.R. CPTAC's file is different
# (1 literal "NA" row, 0 blanks) — see CPTAC_RbE2F_preproc.R for the
# corresponding !is.na() guard needed there.
rna_raw <- rna_raw[rna_raw$Hugo_Symbol != "", ]
rna_raw <- rna_raw[!duplicated(rna_raw$Hugo_Symbol), ]
rownames(rna_raw) <- rna_raw$Hugo_Symbol
rna_raw <- rna_raw[, !colnames(rna_raw) %in% c("Hugo_Symbol", "Entrez_Gene_Id")]

# Transpose: samples as rows
rna_t <- as.data.frame(t(rna_raw))
rna_t[is.na(rna_t)] <- 0

# TCGA z-scores: convert to relative expression (shift +1, clip at 0.01)
orig_rownames <- rownames(rna_t)
rna_t <- as.data.frame(
  lapply(rna_t, function(x) pmax(0.01, as.numeric(x) + 1))
)
rownames(rna_t) <- orig_rownames

# Trim patient IDs to 12 chars (TCGA-XX-XXXX)
patient_ids <- str_sub(rownames(rna_t), end = 12)

# Keep only genes present in RBE2F_GENES
available_genes <- intersect(RBE2F_GENES, colnames(rna_t))
missing_genes    <- setdiff(RBE2F_GENES, colnames(rna_t))
if (length(missing_genes) > 0)
  cat("Warning: genes not found in TCGA data:", paste(missing_genes, collapse = ", "), "\n")

rbe2f_patients <- rna_t[, available_genes, drop = FALSE]
for (g in missing_genes) rbe2f_patients[[g]] <- 1.0
rbe2f_patients <- rbe2f_patients[, RBE2F_GENES]

# Check against clinical data to keep only patients with survival info
clin_raw$PATIENT_ID <- trimws(clin_raw$PATIENT_ID)
clin_raw$OS_MONTHS  <- suppressWarnings(as.numeric(clin_raw$OS_MONTHS))

cat(sprintf("  Clinical patients with OS data: %d\n",
            sum(!is.na(clin_raw$OS_MONTHS) & clin_raw$OS_MONTHS > 0)))

keep_ids <- patient_ids %in% clin_raw$PATIENT_ID[
  !is.na(clin_raw$OS_MONTHS) & clin_raw$OS_MONTHS > 0]

rbe2f_patients <- rbe2f_patients[keep_ids, ]
rbe2f_patients$SAMPLE_ID  <- rownames(rna_t)[keep_ids]
rbe2f_patients$PATIENT_ID <- patient_ids[keep_ids]

# Deduplicate to one primary tumor sample per patient ID
rbe2f_patients <- rbe2f_patients[!duplicated(rbe2f_patients$PATIENT_ID), ]

cat(sprintf("TCGA-GBM patients with RNA + clinical data (deduplicated): %d\n", nrow(rbe2f_patients)))

# 3. Normalize (column-wise mean so reference factor = 1)
gene_cols <- RBE2F_GENES[RBE2F_GENES %in% colnames(rbe2f_patients)]
rbe2f_patients[, gene_cols] <- to_relative(as.matrix(rbe2f_patients[, gene_cols]))

write.csv(rbe2f_patients, "GBM_patient_RbE2Fgenes.csv", row.names = FALSE)
cat("Saved: GBM_patient_RbE2Fgenes.csv\n")
