# Q1 — Extract the Assignment 2 DepMap/GDSC2 TMZ-Sensitivity Score
#
# Purpose: GBM_cellline_coefs.csv (from DepMap_signature.R, Assignment 2) is
# a 10-gene LASSO signature trained on GDSC2 temozolomide response in DepMap
# GBM cell lines. It already had the strongest single-feature discovery
# result of anything tried in this project (TCGA KM p = 0.013, per
# GBM_Assignment2_Report.docx) but was never added to the Q1 combined
# signature. This script projects it onto both patient cohorts using each
# patient's own z-scored expression (same approach GBM_Analysis.R used to
# project it onto TCGA originally), so it can be added as one more candidate
# feature in Q1_build_signature_TCGA.R / Q1_validate_signature_CPTAC.R.
#
# All 10 genes (SV2C, MTF1, PIH1D1, ABCB4, GTF2E1, TAF8, PRPF8, EIF3G, AZI2,
# TCEAL7) were confirmed present in both cohorts' RNA files before writing
# this script.
#
# Output:
#   Q1_TCGA_depmap_score.csv  — PATIENT_ID, DEPMAP_TMZ_SCORE
#   Q1_CPTAC_depmap_score.csv — PATIENT_ID, DEPMAP_TMZ_SCORE
#
# Run this before Q1_build_signature_TCGA.R / Q1_validate_signature_CPTAC.R.

library(tidyverse)
library(stringr)

coefs <- read.csv("GBM_cellline_coefs.csv", stringsAsFactors = FALSE)
cat("DepMap signature genes:", paste(coefs$Gene, collapse = ", "), "\n\n")

score_cohort <- function(rna_file, trim_to_12, out_file, cohort_label) {
  cat(sprintf("=== %s ===\n", cohort_label))

  rna_raw <- read.delim(rna_file, check.names = FALSE, stringsAsFactors = FALSE)
  rna_raw <- rna_raw[!is.na(rna_raw$Hugo_Symbol) & rna_raw$Hugo_Symbol != "", ]
  rna_raw <- rna_raw[!duplicated(rna_raw$Hugo_Symbol), ]
  rownames(rna_raw) <- rna_raw$Hugo_Symbol
  rna_raw <- rna_raw[, !colnames(rna_raw) %in% c("Hugo_Symbol", "Entrez_Gene_Id")]

  rna_t <- as.data.frame(t(rna_raw))
  rna_t[is.na(rna_t)] <- 0
  orig_rownames <- rownames(rna_t)
  rna_t <- as.data.frame(lapply(rna_t, as.numeric))
  rownames(rna_t) <- orig_rownames

  patient_ids <- if (trim_to_12) str_sub(rownames(rna_t), end = 12) else rownames(rna_t)

  available_genes <- intersect(coefs$Gene, colnames(rna_t))
  missing_genes    <- setdiff(coefs$Gene, colnames(rna_t))
  if (length(missing_genes) > 0)
    cat("Warning: signature genes not found:", paste(missing_genes, collapse = ", "), "\n")

  X <- as.matrix(rna_t[, available_genes, drop = FALSE])
  w <- coefs$Coefficient[match(available_genes, coefs$Gene)]
  score <- as.numeric(X %*% w)

  out <- data.frame(PATIENT_ID = patient_ids, DEPMAP_TMZ_SCORE = score)
  out <- out[!duplicated(out$PATIENT_ID), ]  # TCGA: keep first sample per patient if >1
  write.csv(out, out_file, row.names = FALSE)
  cat(sprintf("Saved: %s (%d patients)\n\n", out_file, nrow(out)))
  out
}

tcga_score <- score_cohort(
  rna_file     = "gbm_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",
  trim_to_12   = TRUE,
  out_file     = "Q1_TCGA_depmap_score.csv",
  cohort_label = "TCGA-GBM"
)

cptac_score <- score_cohort(
  rna_file     = "gbm_cptac_2021/data_mrna_seq_fpkm_zscores_ref_all_samples.txt",
  trim_to_12   = FALSE,
  out_file     = "Q1_CPTAC_depmap_score.csv",
  cohort_label = "CPTAC-GBM"
)
