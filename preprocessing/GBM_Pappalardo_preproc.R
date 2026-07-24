# GBM Pappalardo (RTK/RAS/MAPK + PI3K/AKT/mTOR) ODE Preprocessing
#
# Purpose: Extract RTK/PI3K/MAPK pathway gene expression for TCGA-GBM patients,
# normalise to relative-expression units that the Pappalardo ODE model
# (GBM_Pappalardo_model.ipynb) expects.
#
# Model: Pappalardo et al. 2016, "Computational Modeling of PI3K/AKT and MAPK
# Signaling Pathways in Melanoma Cancer" (PLoS ONE), BioModels BIOMD0000000666.
# 48 species / 48 reactions — verified directly from the SBML (via libsbml) before
# writing this script. The melanoma-specific BRAF-V600E/Dabrafenib sub-module was
# removed for GBM reuse (GBM is not BRAF-driven); see GBM_Pappalardo_model.ipynb.
#
# Genes used (11 — covers both the RAS/RAF/MEK/ERK and PI3K/AKT/mTOR arms
# downstream of the same receptor, i.e. both of GBM's RTK-pathway branches
# in one model):
#   EGFR, PTEN, SOS1, KRAS, RAF1, MAP2K1, MAPK1, PIK3CA, AKT1, MTOR, RPS6KB1
#
# Output:
#   GBM_patient_Pappalardogenes.csv — one row per TCGA-GBM patient, used by
#   GBM_Pappalardo_model.ipynb
#
# Run this BEFORE GBM_Pappalardo_model.ipynb.

# Packages
library(tidyverse)
library(stringr)

# 1. File paths
tcga_rna_file  <- "gbm_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
tcga_clin_file <- "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt"

# Eleven RTK/PI3K/MAPK pathway genes (all confirmed present in the TCGA file
# before writing this script)
PAPPALARDO_GENES <- c("EGFR", "PTEN", "SOS1", "KRAS", "RAF1", "MAP2K1",
                       "MAPK1", "PIK3CA", "AKT1", "MTOR", "RPS6KB1")

# Helper: convert expression to ODE-ready relative-expression values
# (identical logic to GBM_p53_preproc.R / GBM_RbE2F_preproc.R, for consistency
# across the project)
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
# TCGA's file has empty (not NA) unmapped Hugo_Symbol entries — verified directly
# against the raw file (13 blank rows, 0 literal "NA" rows), same as in
# GBM_p53_preproc.R and GBM_RbE2F_preproc.R, so the plain != "" filter is safe here.
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

# Keep only genes present in PAPPALARDO_GENES
available_genes <- intersect(PAPPALARDO_GENES, colnames(rna_t))
missing_genes    <- setdiff(PAPPALARDO_GENES, colnames(rna_t))
if (length(missing_genes) > 0)
  cat("Warning: genes not found in TCGA data:", paste(missing_genes, collapse = ", "), "\n")

patient_df <- rna_t[, available_genes, drop = FALSE]
for (g in missing_genes) patient_df[[g]] <- 1.0
patient_df <- patient_df[, PAPPALARDO_GENES]

# Check against clinical data to keep only patients with survival info
clin_raw$PATIENT_ID <- trimws(clin_raw$PATIENT_ID)
clin_raw$OS_MONTHS  <- suppressWarnings(as.numeric(clin_raw$OS_MONTHS))

cat(sprintf("  Clinical patients with OS data: %d\n",
            sum(!is.na(clin_raw$OS_MONTHS) & clin_raw$OS_MONTHS > 0)))

keep_ids <- patient_ids %in% clin_raw$PATIENT_ID[
  !is.na(clin_raw$OS_MONTHS) & clin_raw$OS_MONTHS > 0]

patient_df <- patient_df[keep_ids, ]
patient_df$SAMPLE_ID  <- rownames(rna_t)[keep_ids]
patient_df$PATIENT_ID <- patient_ids[keep_ids]

cat(sprintf("TCGA-GBM patients with RNA + clinical data: %d\n", nrow(patient_df)))

# 3. Normalize (column-wise mean so reference factor = 1)
gene_cols <- PAPPALARDO_GENES[PAPPALARDO_GENES %in% colnames(patient_df)]
patient_df[, gene_cols] <- to_relative(as.matrix(patient_df[, gene_cols]))

write.csv(patient_df, "GBM_patient_Pappalardogenes.csv", row.names = FALSE)
cat("Saved: GBM_patient_Pappalardogenes.csv\n")
