# Q1 v2 — Build a Cross-Dataset Survival Signature Using Pathway-Level Scores
#
# Purpose: The v1 signature (Q1_build_signature_TCGA.R) pooled 47-48
# individual gene-expression and per-DDR-level ODE features. A follow-up
# test showed that ADDING more of those individual features made TCGA's
# in-sample fit better but CPTAC generalisation slightly worse (C-index
# 0.587 -> 0.576, Cox p 0.037 -> 0.054) — classic overfitting, since many of
# those features are correlated with each other (e.g. the 10 DDR levels of
# the p53 ODE model move together) and LASSO was picking among near-
# duplicates rather than genuinely independent signal.
#
# This version takes the opposite approach: collapse each pathway down to
# ONE summary score instead of many correlated ones, trading per-gene
# granularity (completeness) for a much healthier feature-to-patient ratio
# (availability of a stable signal) — the "compromise between data
# availability and data completeness" the brief asks for, applied directly.
#
# Candidate features (8 total, down from 48, n ~150 TCGA patients):
#   - P53_SCORE       = mean relative expression of the 8 p53 pathway genes
#   - RBE2F_SCORE     = mean relative expression of the 5 Rb-E2F genes
#   - PAPPALARDO_SCORE = mean relative expression of the 11 RTK/PI3K/MAPK genes
#   - P53_ODE_S15_MEAN, P53_ODE_S46_MEAN = each patient's p53 ODE output
#     averaged across all 10 DDR levels (rather than picking one "best"
#     level, which would just reintroduce single-feature cherry-picking)
#   - AGE, SEX_MALE, IDH_MUT
#
# Same discipline as v1: trains ONLY on TCGA-GBM, saves frozen coefficients,
# Q1_validate_signature_CPTAC_v2.R applies them unchanged to CPTAC-GBM.
# Uses "_v2" filenames throughout so it never overwrites the v1 result files
# (a real mix-up that happened during the v1 robustness check).
#
# Prerequisites (run first, in this order):
#   GBM_p53_preproc.R, GBM_RbE2F_preproc.R, GBM_Pappalardo_preproc.R,
#   GBM_p53_model.ipynb (produces GBM_p53s15DR_patients.csv / s46),
#   Q1_extract_IDH_status.R

library(tidyverse)
library(glmnet)
library(survival)
library(survminer)
library(stringr)

set.seed(42)
FIG_DIR <- getwd()

cat("=== Q1 v2: Building pathway-level cross-dataset signature — TCGA-GBM training ===\n\n")

# 1. Clinical / survival data
clin_raw <- read.delim(
  "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
  sep = "\t", skip = 4, stringsAsFactors = FALSE
)
clin_raw$PATIENT_ID  <- trimws(clin_raw$PATIENT_ID)
clin_raw$SURV_TIME   <- suppressWarnings(as.numeric(clin_raw$OS_MONTHS))
clin_raw$SURV_STATUS <- as.integer(clin_raw$OS_STATUS == "1:DECEASED")
clin_raw$AGE         <- suppressWarnings(as.numeric(clin_raw$AGE))
clin_raw$SEX_MALE    <- as.integer(clin_raw$SEX == "Male")
clin_raw <- clin_raw[!is.na(clin_raw$SURV_TIME) & clin_raw$SURV_TIME > 0 &
                      !is.na(clin_raw$AGE), ]
clin <- clin_raw[, c("PATIENT_ID", "SURV_TIME", "SURV_STATUS", "AGE", "SEX_MALE")]
cat(sprintf("TCGA-GBM patients with OS + age + sex: %d\n", nrow(clin)))

# 2. Pathway gene expression -> one mean score per pathway
p53_genes   <- read.csv("GBM_patient_p53genes.csv", stringsAsFactors = FALSE)
rbe2f_genes <- read.csv("GBM_patient_RbE2Fgenes.csv", stringsAsFactors = FALSE)
pappa_genes <- read.csv("GBM_patient_Pappalardogenes.csv", stringsAsFactors = FALSE)

P53_GENE_COLS  <- c("ATM", "CHEK2", "HIPK2", "MDM2", "PPM1D", "SIAH1", "TP53", "WSB1")
RBE2F_GENE_COLS <- c("MYC", "CCND1", "CCNE1", "RB1", "E2F1")
PAPPA_GENE_COLS <- c("EGFR", "PTEN", "SOS1", "KRAS", "RAF1", "MAP2K1", "MAPK1", "PIK3CA", "AKT1", "MTOR", "RPS6KB1")

# TCGA quirk: 5 patients have two tumour samples each, so PATIENT_ID is not
# unique in these per-sample gene files. Deduplicate to one row per patient
# (first occurrence = primary tumour "-01" sample in every case checked)
# before computing scores, otherwise those patients' scores would be
# duplicated rows going into the join chain below (see Q1_build_signature_TCGA.R
# for the full diagnosis of this bug, found via a pandas replication of the
# join that showed one patient represented 32x in an earlier run).
p53_genes   <- p53_genes[!duplicated(p53_genes$PATIENT_ID), ]
rbe2f_genes <- rbe2f_genes[!duplicated(rbe2f_genes$PATIENT_ID), ]
pappa_genes <- pappa_genes[!duplicated(pappa_genes$PATIENT_ID), ]

p53_score   <- data.frame(PATIENT_ID = p53_genes$PATIENT_ID,
                           P53_SCORE = rowMeans(p53_genes[, P53_GENE_COLS]))
rbe2f_score <- data.frame(PATIENT_ID = rbe2f_genes$PATIENT_ID,
                           RBE2F_SCORE = rowMeans(rbe2f_genes[, RBE2F_GENE_COLS]))
pappa_score <- data.frame(PATIENT_ID = pappa_genes$PATIENT_ID,
                           PAPPALARDO_SCORE = rowMeans(pappa_genes[, PAPPA_GENE_COLS]))

# 3. p53 ODE steady-state outputs -> one mean score per phosphosite
# (averaged across all 10 DDR levels, not one hand-picked level)
p53s15 <- read.csv("GBM_p53s15DR_patients.csv", stringsAsFactors = FALSE)
p53s46 <- read.csv("GBM_p53s46DR_patients.csv", stringsAsFactors = FALSE)
p53s15 <- p53s15[!duplicated(p53s15$PATIENT_ID), ]
p53s46 <- p53s46[!duplicated(p53s46$PATIENT_ID), ]
ddr_cols <- grep("^DDR_", colnames(p53s15), value = TRUE)
p53_ode_score <- data.frame(
  PATIENT_ID = p53s15$PATIENT_ID,
  P53_ODE_S15_MEAN = rowMeans(p53s15[, ddr_cols]),
  P53_ODE_S46_MEAN = rowMeans(p53s46[match(p53s15$PATIENT_ID, p53s46$PATIENT_ID), ddr_cols])
)

# 4. IDH mutation status
idh <- read.csv("Q1_TCGA_IDH_status.csv", stringsAsFactors = FALSE)

# 5. Merge everything by PATIENT_ID
feat <- clin %>%
  inner_join(p53_score,     by = "PATIENT_ID") %>%
  inner_join(rbe2f_score,   by = "PATIENT_ID") %>%
  inner_join(pappa_score,   by = "PATIENT_ID") %>%
  inner_join(p53_ode_score, by = "PATIENT_ID") %>%
  inner_join(idh,           by = "PATIENT_ID")

cat(sprintf("TCGA-GBM patients with ALL feature sources present: %d\n", nrow(feat)))

feature_cols <- setdiff(colnames(feat), c("PATIENT_ID", "SURV_TIME", "SURV_STATUS"))
cat(sprintf("Candidate features: %d\n", length(feature_cols)))
cat(paste(feature_cols, collapse = ", "), "\n\n")

X <- as.matrix(feat[, feature_cols])
storage.mode(X) <- "double"
y <- Surv(feat$SURV_TIME, feat$SURV_STATUS)

# 6. Fit penalised Cox with cross-validation, LASSO vs Elastic Net, same
# methodology as v1 for a fair comparison. With only 8 candidates and ~150
# patients, the penalty may end up doing very little — that's expected and
# fine; the point of v2 is a healthier feature-to-patient ratio, not
# necessarily heavier regularisation.
cat("Running cross-validated penalised Cox (LASSO vs Elastic Net)...\n")

fit_one_alpha <- function(alpha_val) {
  tryCatch(
    cv.glmnet(X, y, family = "cox", type.measure = "C", alpha = alpha_val, nfolds = 10, standardize = TRUE),
    error = function(e) {
      cat(sprintf("type.measure = 'C' failed at alpha=%.1f — falling back to 'deviance'.\n", alpha_val))
      cv.glmnet(X, y, family = "cox", type.measure = "deviance", alpha = alpha_val, nfolds = 10, standardize = TRUE)
    }
  )
}

cv_lasso <- fit_one_alpha(1)
cv_enet  <- fit_one_alpha(0.5)

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

png(file.path(FIG_DIR, "fig_Q1_signature_v2_cv.png"), width = 1400, height = 1000, res = 150)
plot(cv_fit)
title("Q1 v2 signature: cross-validated C-index vs log(lambda)", line = 2.5)
dev.off()

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
  stop("Zero features selected at lambda.1se — re-run with s = cv_fit$lambda.min, or inspect fig_Q1_signature_v2_cv.png.")
}

# 7. Save the frozen model — DIFFERENT filename from v1, so v1's headline
# result (Q1_signature_coefficients.csv, and its figures) is never touched.
write.csv(coef_df, "Q1_signature_coefficients_v2.csv", row.names = FALSE)
cat("\nSaved: Q1_signature_coefficients_v2.csv (frozen model, apply unchanged to CPTAC)\n")

# 8. In-sample TCGA check (context only)
risk_score <- as.numeric(X[, coef_df$feature, drop = FALSE] %*% coef_df$coefficient)
fit_check <- coxph(y ~ risk_score)
c_index_tcga <- summary(fit_check)$concordance[1]
cat(sprintf("\nTCGA in-sample check: C-index = %.3f, Cox p = %.4g\n",
            c_index_tcga, summary(fit_check)$coefficients[, "Pr(>|z|)"]))

med <- median(risk_score)
grp <- ifelse(risk_score > med, paste0("High risk score (n=", sum(risk_score > med), ")"),
                                 paste0("Low risk score (n=",  sum(risk_score <= med), ")"))
km_df <- data.frame(grp = grp)
km_fit <- survfit(y ~ grp, data = km_df)
p_km <- ggsurvplot(
  km_fit, data = km_df, pval = TRUE,
  xlab = "Overall Survival (months)", ylab = "Survival Probability",
  palette = c("#E64B35", "#4DBBD5"),
  title = "Q1 v2 pathway-level signature — TCGA-GBM (training set, median split)",
  legend.title = "", risk.table = TRUE, risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 13)
)
png(file.path(FIG_DIR, "fig_Q1_signature_v2_TCGA_km.png"), width = 1400, height = 1100, res = 150)
print(p_km)
dev.off()

cat("\nDone. Next: run Q1_validate_signature_CPTAC_v2.R to test this frozen signature on CPTAC-GBM.\n")
