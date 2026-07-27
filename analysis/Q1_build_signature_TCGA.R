# Q1 — Build a Cross-Dataset Survival Signature: TCGA-GBM Training
#
# (Duplicate of Q1_build_signature_TCGA.R, recreated under a new filename
# after the original became locked/inaccessible via a mount issue. Content
# is identical, including the TCGA duplicate-patient dedup fix. Once the
# original file is accessible again, this duplicate can be deleted.)
#
# Purpose: Answer the brief's Q1 ("build a predictor of patient response
# which works across different datasets") properly, by fixing the two
# things that made the earlier single-feature p53 test fail to replicate
# in CPTAC (see GBM_ThreeModel_ODE_Report.docx, Model 1):
#   1. A single feature hand-picked by scanning 20 candidates overfits.
#      Here, ALL candidate features are pooled and a penalised Cox model
#      (LASSO or Elastic Net, whichever cross-validates better — see below)
#      selects and weights them together, with both the penalty type and
#      strength chosen by cross-validation, not by picking whichever single
#      feature looks best.
#   2. Age, sex, and IDH mutation status (the strongest established GBM
#      prognostic markers this project's data can actually provide — see
#      Q1_extract_IDH_status.R for why MGMT isn't included) are added as
#      covariates, so any pathway signal has to earn its place alongside
#      known predictors, not be tested in isolation.
#
# Candidate features pooled here (48 total, n ~150 TCGA patients):
#   - 24 pathway genes' relative expression (p53: 8, Rb-E2F: 5, Pappalardo: 11)
#   - 20 p53 ODE steady-state outputs (10 DDR levels x 2 phosphosites) —
#     the only ODE-derived features used here, because they're the only ones
#     that also exist for CPTAC (Rb-E2F/Pappalardo ODE models were never run
#     on CPTAC, since their TCGA discovery was a clean null — see report).
#     This is the availability/completeness compromise: raw expression for
#     all three pathways in both cohorts, but ODE-simulated dynamics only
#     where both cohorts have them.
#   - AGE, SEX (dummy), IDH_MUT
#   - DEPMAP_TMZ_SCORE — the Assignment 2 LASSO/GDSC2 TMZ-sensitivity
#     signature (10 genes, all confirmed present in both cohorts), the
#     strongest single-feature discovery result in this project before now,
#     never previously tested for cross-dataset replication.
#
# This script TRAINS ONLY on TCGA-GBM (discovery). It saves the frozen,
# cross-validated coefficients to Q1_signature_coefficients.csv. The
# companion script, Q1_validate_signature_CPTAC.R, applies those coefficients
# unchanged to CPTAC-GBM — no re-fitting, no re-scanning, same discipline as
# the earlier p53 validation.
#
# Prerequisites (run first, in this order):
#   GBM_p53_preproc.R, GBM_RbE2F_preproc.R, GBM_Pappalardo_preproc.R,
#   GBM_p53_model.ipynb (produces GBM_p53s15DR_patients.csv / s46),
#   Q1_extract_IDH_status.R, Q1_extract_depmap_signature.R

library(tidyverse)
library(glmnet)
library(survival)
library(survminer)
library(stringr)

set.seed(42)  # cv.glmnet's fold assignment is random; fixed for reproducibility

FIG_DIR <- getwd()

cat("=== Q1: Building cross-dataset signature — TCGA-GBM training ===\n\n")

# 1. Clinical / survival data
clin_raw <- read.delim(
  "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
  sep = "\t", skip = 4, stringsAsFactors = FALSE
)
clin_raw$PATIENT_ID   <- trimws(clin_raw$PATIENT_ID)
clin_raw$SURV_TIME    <- suppressWarnings(as.numeric(clin_raw$OS_MONTHS))
clin_raw$SURV_STATUS  <- as.integer(clin_raw$OS_STATUS == "1:DECEASED")
clin_raw$AGE          <- suppressWarnings(as.numeric(clin_raw$AGE))
clin_raw$SEX_MALE     <- as.integer(clin_raw$SEX == "Male")
clin_raw <- clin_raw[!is.na(clin_raw$SURV_TIME) & clin_raw$SURV_TIME > 0 &
                      !is.na(clin_raw$AGE), ]
clin <- clin_raw[, c("PATIENT_ID", "SURV_TIME", "SURV_STATUS", "AGE", "SEX_MALE")]
cat(sprintf("TCGA-GBM patients with OS + age + sex: %d\n", nrow(clin)))

# 2. Pathway gene expression (3 files, 24 genes total)
p53_genes  <- read.csv("GBM_patient_p53genes.csv", stringsAsFactors = FALSE)
rbe2f_genes <- read.csv("GBM_patient_RbE2Fgenes.csv", stringsAsFactors = FALSE)
pappa_genes <- read.csv("GBM_patient_Pappalardogenes.csv", stringsAsFactors = FALSE)

p53_genes   <- p53_genes[,   !colnames(p53_genes)   %in% "SAMPLE_ID"]
rbe2f_genes <- rbe2f_genes[, !colnames(rbe2f_genes) %in% "SAMPLE_ID"]
pappa_genes <- pappa_genes[, !colnames(pappa_genes) %in% "SAMPLE_ID"]

# TCGA quirk: 5 patients (TCGA-06-0125, -0190, -0210, -0211, TCGA-14-1034)
# have TWO tumour samples each (a primary "-01" and a second "-02"), so
# PATIENT_ID is not unique in these per-sample gene files. Left as-is, the
# chained inner_join() calls below treat this as many-to-many and compound
# it: one patient ended up represented 32 times in an earlier run of this
# script (confirmed by replicating the join in pandas). Deduplicate to one
# row per patient before joining — first occurrence, which is the primary
# tumour ("-01") sample in every case checked.
p53_genes   <- p53_genes[!duplicated(p53_genes$PATIENT_ID), ]
rbe2f_genes <- rbe2f_genes[!duplicated(rbe2f_genes$PATIENT_ID), ]
pappa_genes <- pappa_genes[!duplicated(pappa_genes$PATIENT_ID), ]

# 3. p53 ODE steady-state outputs (20 features)
p53s15 <- read.csv("GBM_p53s15DR_patients.csv", stringsAsFactors = FALSE)
p53s46 <- read.csv("GBM_p53s46DR_patients.csv", stringsAsFactors = FALSE)
ddr_cols <- grep("^DDR_", colnames(p53s15), value = TRUE)
p53s15 <- p53s15[, c(ddr_cols, "PATIENT_ID")]
p53s46 <- p53s46[, c(ddr_cols, "PATIENT_ID")]
p53s15 <- p53s15[!duplicated(p53s15$PATIENT_ID), ]
p53s46 <- p53s46[!duplicated(p53s46$PATIENT_ID), ]
colnames(p53s15)[colnames(p53s15) %in% ddr_cols] <- paste0("s15_", ddr_cols)
colnames(p53s46)[colnames(p53s46) %in% ddr_cols] <- paste0("s46_", ddr_cols)

# 4. IDH mutation status
idh <- read.csv("Q1_TCGA_IDH_status.csv", stringsAsFactors = FALSE)

# 4b. Assignment 2 DepMap/GDSC2 TMZ-sensitivity score — added as one more
# candidate feature. This already had the strongest single-feature discovery
# result of anything tried in this project (TCGA KM p = 0.013) but was never
# folded into a cross-dataset test until now.
depmap <- read.csv("Q1_TCGA_depmap_score.csv", stringsAsFactors = FALSE)

# 5. Merge everything by PATIENT_ID (inner join — only patients present in
# every source go into training, so the frozen model never depends on a
# feature that's missing for some patients)
feat <- clin %>%
  inner_join(p53_genes,   by = "PATIENT_ID") %>%
  inner_join(rbe2f_genes, by = "PATIENT_ID") %>%
  inner_join(pappa_genes, by = "PATIENT_ID") %>%
  inner_join(p53s15,      by = "PATIENT_ID") %>%
  inner_join(p53s46,      by = "PATIENT_ID") %>%
  inner_join(idh,         by = "PATIENT_ID") %>%
  inner_join(depmap,      by = "PATIENT_ID")

cat(sprintf("TCGA-GBM patients with ALL feature sources present: %d\n", nrow(feat)))

feature_cols <- setdiff(colnames(feat), c("PATIENT_ID", "SURV_TIME", "SURV_STATUS"))
cat(sprintf("Candidate features: %d\n", length(feature_cols)))
cat(paste(feature_cols, collapse = ", "), "\n\n")

X <- as.matrix(feat[, feature_cols])
storage.mode(X) <- "double"
y <- Surv(feat$SURV_TIME, feat$SURV_STATUS)

# 6. Fit a penalised Cox model with cross-validation to choose both the
# penalty type (alpha) and strength (lambda). type.measure = "C" scores each
# candidate by cross-validated Harrell's C-index — directly the metric this
# signature will be judged on, rather than a proxy. (Falls back to
# "deviance" if your glmnet version is too old for type.measure = "C" —
# check with packageVersion("glmnet"); "C" needs glmnet >= 4.1.)
#
# alpha = 1 is pure LASSO (used previously); alpha = 0.5 is Elastic Net,
# which tends to keep correlated genes together instead of arbitrarily
# dropping one (e.g. MAPK1 and PTEN, both RTK/PI3K/MAPK pathway members that
# LASSO selected separately — Elastic Net may handle that correlation
# differently). Both are fit and the one with the better cross-validated
# C-index at its own lambda.1se is kept — chosen by the data, not assumed.
cat("Running cross-validated penalised Cox (LASSO vs Elastic Net; this can take a few seconds)...\n")

fit_one_alpha <- function(alpha_val) {
  tryCatch(
    cv.glmnet(X, y, family = "cox", type.measure = "C", alpha = alpha_val, nfolds = 10, standardize = TRUE),
    error = function(e) {
      cat(sprintf("type.measure = 'C' failed at alpha=%.1f (old glmnet version?) — falling back to 'deviance'.\n", alpha_val))
      cv.glmnet(X, y, family = "cox", type.measure = "deviance", alpha = alpha_val, nfolds = 10, standardize = TRUE)
    }
  )
}

cv_lasso  <- fit_one_alpha(1)    # LASSO
cv_enet   <- fit_one_alpha(0.5)  # Elastic Net

cindex_at_1se <- function(cv) cv$cvm[cv$lambda == cv$lambda.1se]
c_lasso <- cindex_at_1se(cv_lasso)
c_enet  <- cindex_at_1se(cv_enet)
cat(sprintf("\nCV score at lambda.1se — LASSO (alpha=1): %.4f | Elastic Net (alpha=0.5): %.4f\n",
            c_lasso, c_enet))

if (c_enet > c_lasso) {
  cat("Elastic Net wins on cross-validated score — using alpha = 0.5.\n")
  cv_fit <- cv_enet
} else {
  cat("LASSO wins (or ties) on cross-validated score — using alpha = 1.\n")
  cv_fit <- cv_lasso
}

png(file.path(FIG_DIR, "fig_Q1_signature_cv.png"), width = 1400, height = 1000, res = 150)
plot(cv_fit)
title("Q1 signature: cross-validated C-index vs log(lambda)", line = 2.5)
dev.off()

# lambda.1se (not lambda.min): the more conservative, sparser choice.
# With ~150 patients and now 48 candidate features, lambda.min tends to keep
# more coefficients and overfit; lambda.1se trades a little in-sample fit
# for a model more likely to actually generalise to CPTAC, which is the
# whole point of this exercise.
lambda_use <- cv_fit$lambda.1se
cat(sprintf("\nlambda.min = %.5f, lambda.1se = %.5f (using lambda.1se)\n",
            cv_fit$lambda.min, lambda_use))

coefs <- coef(cv_fit, s = lambda_use)
coef_df <- data.frame(feature = rownames(coefs), coefficient = as.numeric(coefs[, 1])) %>%
  filter(coefficient != 0) %>%
  arrange(desc(abs(coefficient)))

cat(sprintf("\nFeatures selected (non-zero coefficient) at lambda.1se: %d\n", nrow(coef_df)))
print(coef_df)

if (nrow(coef_df) == 0) {
  stop("LASSO selected zero features at lambda.1se — the penalty is too strong. ",
       "Re-run with s = cv_fit$lambda.min instead, or inspect fig_Q1_signature_cv.png.")
}

# 7. Save the frozen model: feature list + coefficients + the lambda used.
# Q1_validate_signature_CPTAC.R reads this file directly — it does NOT refit.
write.csv(coef_df, "Q1_signature_coefficients.csv", row.names = FALSE)
cat("\nSaved: Q1_signature_coefficients.csv (frozen model, apply unchanged to CPTAC)\n")

# 8. In-sample TCGA check (context only — NOT the generalisation test;
# CPTAC is). Risk score = linear predictor from the selected features.
risk_score <- as.numeric(X[, coef_df$feature, drop = FALSE] %*% coef_df$coefficient)
fit_check <- coxph(y ~ risk_score)
c_index_tcga <- summary(fit_check)$concordance[1]
cat(sprintf("\nTCGA in-sample check: C-index = %.3f, Cox p = %.4g\n",
            c_index_tcga, summary(fit_check)$coefficients[, "Pr(>|z|)"]))

# Median-split KM (illustrative only — median split, not threshold-optimised,
# to avoid the small-unstable-subgroup artefact seen in Models 2 and 3 of
# the ODE report)
med <- median(risk_score)
grp <- ifelse(risk_score > med, paste0("High risk score (n=", sum(risk_score > med), ")"),
                                 paste0("Low risk score (n=",  sum(risk_score <= med), ")"))
km_df <- data.frame(grp = grp)
km_fit <- survfit(y ~ grp, data = km_df)
p_km <- ggsurvplot(
  km_fit, data = km_df, pval = TRUE,
  xlab = "Overall Survival (months)", ylab = "Survival Probability",
  palette = c("#E64B35", "#4DBBD5"),
  title = "Q1 combined signature — TCGA-GBM (training set, median split)",
  legend.title = "", risk.table = TRUE, risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 13)
)
png(file.path(FIG_DIR, "fig_Q1_signature_TCGA_km.png"), width = 1400, height = 1100, res = 150)
print(p_km)
dev.off()

cat("\nDone. Next: run Q1_validate_signature_CPTAC.R to test this frozen signature on CPTAC-GBM.\n")
