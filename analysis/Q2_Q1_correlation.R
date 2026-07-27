# Q2-Q1 Correlation — Does the Drug-Viability Predictor Track the
# Patient-Response Predictor?
#
# Purpose: the brief's Q2 asks for a predictor of cell viability with
# respect to standard-of-care drug treatment (TMZ), and then explicitly
# asks: "Does it in any way correlate with p.1?" This has NOT been tested
# yet. Adding DEPMAP_TMZ_SCORE as one candidate feature inside the Q1
# LASSO model (the earlier "robustness check") is a different question —
# it asks whether the score helps predict survival jointly with other
# features. This script asks the question the brief actually poses:
# do the two predictors' OUTPUTS correlate with each other, patient by
# patient, independent of any shared model?
#
#   Q2 predictor = DEPMAP_TMZ_SCORE (Assignment 2's LASSO/GDSC2 signature,
#                  trained on GBM cell-line TMZ sensitivity, projected onto
#                  each patient's own expression)
#   Q1 predictor = the headline 9-feature signature risk score from
#                  Q1_build_signature_TCGA.R / Q1_validate_signature_CPTAC.R
#                  (MDM2, SIAH1, MAPK1, PTEN, p53-Ser15 @ DDR=0.12, MTOR,
#                  RB1, HIPK2, AGE), computed here from the coefficients
#                  frozen in the report — NOT re-read from
#                  Q1_signature_coefficients.csv, because that file on disk
#                  currently holds the later 11-feature DepMap-extended
#                  version, and correlating DEPMAP_TMZ_SCORE against a risk
#                  score that already contains DEPMAP_TMZ_SCORE as an input
#                  would be circular.
#
# Tested in both cohorts separately (TCGA-GBM = discovery, CPTAC-GBM =
# independent), since a correlation that only shows up in one cohort is a
# weaker claim than one that replicates in both.
#
# Prerequisites (run first, in this order):
#   GBM_p53_preproc.R, GBM_RbE2F_preproc.R, GBM_Pappalardo_preproc.R,
#   GBM_p53_model.ipynb, CPTAC_p53_preproc.R, CPTAC_RbE2F_preproc.R,
#   CPTAC_Pappalardo_preproc.R, CPTAC_p53_model.ipynb,
#   Q1_extract_depmap_signature.R

library(tidyverse)
library(stringr)

FIG_DIR <- getwd()

# Frozen v1 headline coefficients (from GBM_ThreeModel_ODE_Report.docx,
# "A Combined Cross-Dataset Signature" section). Hardcoded here deliberately
# so this script does not depend on which version is currently saved to
# Q1_signature_coefficients.csv on disk.
V1_COEFS <- c(
  MDM2            =  0.2296,
  SIAH1           = -0.1555,
  MAPK1           = -0.1539,
  PTEN            = -0.1379,
  s15_DDR_0.120   = -0.1269,
  MTOR            =  0.0620,
  RB1             = -0.0530,
  HIPK2           =  0.0284,
  AGE             =  0.0251
)

build_risk_score <- function(feat) {
  X <- as.matrix(feat[, names(V1_COEFS), drop = FALSE])
  storage.mode(X) <- "double"
  as.numeric(X %*% V1_COEFS)
}

run_cohort <- function(cohort_label, clin, p53_genes, rbe2f_genes, pappa_genes,
                        p53s15, depmap) {

  cat(sprintf("\n=== %s ===\n", cohort_label))

  # TCGA quirk (CPTAC is clean): 5 TCGA patients have two tumour samples
  # each, so PATIENT_ID is not unique in these per-sample files. Left
  # un-deduplicated, the chained inner_join() below compounds this into
  # severe pseudo-replication (confirmed via pandas: one patient ended up
  # represented 32x in an earlier run of the Q1 build script). Deduplicate
  # to one row per patient before joining, for both cohorts (a no-op for
  # CPTAC, which has no duplicate PATIENT_IDs to begin with).
  p53_genes   <- p53_genes[!duplicated(p53_genes$PATIENT_ID), ]
  rbe2f_genes <- rbe2f_genes[!duplicated(rbe2f_genes$PATIENT_ID), ]
  pappa_genes <- pappa_genes[!duplicated(pappa_genes$PATIENT_ID), ]
  p53s15      <- p53s15[!duplicated(p53s15$PATIENT_ID), ]

  ddr_cols <- grep("^DDR_", colnames(p53s15), value = TRUE)
  p53s15 <- p53s15[, c(ddr_cols, "PATIENT_ID")]
  colnames(p53s15)[colnames(p53s15) %in% ddr_cols] <- paste0("s15_", ddr_cols)

  feat <- clin %>%
    inner_join(p53_genes,   by = "PATIENT_ID") %>%
    inner_join(rbe2f_genes, by = "PATIENT_ID") %>%
    inner_join(pappa_genes, by = "PATIENT_ID") %>%
    inner_join(p53s15,      by = "PATIENT_ID") %>%
    inner_join(depmap,      by = "PATIENT_ID")

  cat(sprintf("Patients with both Q1 inputs and Q2 (DepMap) score present: %d\n", nrow(feat)))

  feat$Q1_risk_score <- build_risk_score(feat)
  feat$Q2_viability_score <- feat$DEPMAP_TMZ_SCORE

  pear <- cor.test(feat$Q1_risk_score, feat$Q2_viability_score, method = "pearson")
  spear <- cor.test(feat$Q1_risk_score, feat$Q2_viability_score, method = "spearman")

  cat(sprintf("Pearson  r = %.3f, p = %.4g\n", pear$estimate, pear$p.value))
  cat(sprintf("Spearman rho = %.3f, p = %.4g\n", spear$estimate, spear$p.value))

  p_scatter <- ggplot(feat, aes(x = Q2_viability_score, y = Q1_risk_score)) +
    geom_point(alpha = 0.6, color = "#4DBBD5") +
    geom_smooth(method = "lm", se = TRUE, color = "#E64B35") +
    labs(
      title = sprintf("Q2 drug-viability score vs Q1 risk score — %s", cohort_label),
      subtitle = sprintf("Pearson r = %.3f (p = %.3g); n = %d", pear$estimate, pear$p.value, nrow(feat)),
      x = "Q2: DepMap/GDSC2 TMZ-sensitivity score (Assignment 2 signature)",
      y = "Q1: combined signature risk score (this report)"
    ) +
    theme_bw(base_size = 13)

  fname <- sprintf("fig_Q2_Q1_correlation_%s.png", gsub("[^A-Za-z0-9]", "", cohort_label))
  png(file.path(FIG_DIR, fname), width = 1400, height = 1000, res = 150)
  print(p_scatter)
  dev.off()
  cat(sprintf("Saved: %s\n", fname))

  data.frame(
    cohort = cohort_label, n = nrow(feat),
    pearson_r = as.numeric(pear$estimate), pearson_p = pear$p.value,
    spearman_rho = as.numeric(spear$estimate), spearman_p = spear$p.value
  )
}

# ---------------------------------------------------------------- TCGA-GBM
tcga_clin_raw <- read.delim("gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
                             sep = "\t", skip = 4, stringsAsFactors = FALSE)
tcga_clin_raw$PATIENT_ID <- trimws(tcga_clin_raw$PATIENT_ID)
tcga_clin_raw$AGE <- suppressWarnings(as.numeric(tcga_clin_raw$AGE))
tcga_clin <- tcga_clin_raw[!is.na(tcga_clin_raw$AGE), c("PATIENT_ID", "AGE")]

tcga_p53   <- read.csv("GBM_patient_p53genes.csv", stringsAsFactors = FALSE)
tcga_rbe2f <- read.csv("GBM_patient_RbE2Fgenes.csv", stringsAsFactors = FALSE)
tcga_pappa <- read.csv("GBM_patient_Pappalardogenes.csv", stringsAsFactors = FALSE)
tcga_p53   <- tcga_p53[,   !colnames(tcga_p53)   %in% "SAMPLE_ID"]
tcga_rbe2f <- tcga_rbe2f[, !colnames(tcga_rbe2f) %in% "SAMPLE_ID"]
tcga_pappa <- tcga_pappa[, !colnames(tcga_pappa) %in% "SAMPLE_ID"]
tcga_p53s15 <- read.csv("GBM_p53s15DR_patients.csv", stringsAsFactors = FALSE)
tcga_depmap <- read.csv("Q1_TCGA_depmap_score.csv", stringsAsFactors = FALSE)

tcga_result <- run_cohort("TCGA-GBM", tcga_clin, tcga_p53, tcga_rbe2f, tcga_pappa,
                           tcga_p53s15, tcga_depmap)

# ---------------------------------------------------------------- CPTAC-GBM
cptac_clin_raw <- read.delim("gbm_cptac_2021/data_clinical_patient.txt",
                              sep = "\t", skip = 4, stringsAsFactors = FALSE)
cptac_clin_raw$PATIENT_ID <- trimws(cptac_clin_raw$PATIENT_ID)
cptac_clin_raw$AGE <- suppressWarnings(as.numeric(cptac_clin_raw$AGE))
cptac_clin <- cptac_clin_raw[!is.na(cptac_clin_raw$AGE), c("PATIENT_ID", "AGE")]

cptac_p53   <- read.csv("CPTAC_patient_p53genes.csv", stringsAsFactors = FALSE)
cptac_rbe2f <- read.csv("CPTAC_patient_RbE2Fgenes.csv", stringsAsFactors = FALSE)
cptac_pappa <- read.csv("CPTAC_patient_Pappalardogenes.csv", stringsAsFactors = FALSE)
cptac_p53   <- cptac_p53[,   !colnames(cptac_p53)   %in% "SAMPLE_ID"]
cptac_rbe2f <- cptac_rbe2f[, !colnames(cptac_rbe2f) %in% "SAMPLE_ID"]
cptac_pappa <- cptac_pappa[, !colnames(cptac_pappa) %in% "SAMPLE_ID"]
cptac_p53s15 <- read.csv("CPTAC_p53s15DR_patients.csv", stringsAsFactors = FALSE)
cptac_depmap <- read.csv("Q1_CPTAC_depmap_score.csv", stringsAsFactors = FALSE)

cptac_result <- run_cohort("CPTAC-GBM", cptac_clin, cptac_p53, cptac_rbe2f, cptac_pappa,
                            cptac_p53s15, cptac_depmap)

# ---------------------------------------------------------------- Summary
results <- bind_rows(tcga_result, cptac_result)
write.csv(results, "Q2_Q1_correlation_results.csv", row.names = FALSE)

cat("\n=== Q2-Q1 correlation summary ===\n")
print(results)
cat("\nSaved: Q2_Q1_correlation_results.csv, fig_Q2_Q1_correlation_TCGAGBM.png, fig_Q2_Q1_correlation_CPTACGBM.png\n")
cat("Interpretation guide: |r| < 0.2 ~ negligible, 0.2-0.4 ~ weak, 0.4-0.6 ~ moderate, > 0.6 ~ strong.\n")
cat("A significant NEGATIVE correlation would be the biologically coherent direction:\n")
cat("higher predicted TMZ sensitivity (better drug response) should track LOWER predicted\n")
cat("risk score (better survival) if the two predictors are picking up related biology.\n")
