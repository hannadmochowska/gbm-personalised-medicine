# Q1 — Validate the Cross-Dataset Signature: CPTAC-GBM (Independent Test)
#
# Purpose: Apply the frozen signature from Q1_build_signature_TCGA.R
# (Q1_signature_coefficients.csv) to CPTAC-GBM UNCHANGED — same features,
# same coefficients, no re-fitting, no re-scanning for a better threshold.
# This is the actual answer to Q1: does the combined signature generalise
# to an independent cohort? A single number (the C-index / Cox p-value
# below) answers it either way.
#
# Prerequisites (run first, in this order):
#   CPTAC_p53_preproc.R, CPTAC_RbE2F_preproc.R, CPTAC_Pappalardo_preproc.R,
#   CPTAC_p53_model.ipynb (produces CPTAC_p53s15DR_patients.csv / s46),
#   Q1_extract_IDH_status.R, Q1_extract_depmap_signature.R,
#   Q1_build_signature_TCGA.R (produces Q1_signature_coefficients.csv)

library(tidyverse)
library(survival)
library(survminer)
library(stringr)

FIG_DIR <- getwd()

cat("=== Q1: Validating cross-dataset signature — CPTAC-GBM (independent test) ===\n\n")

# 1. Clinical / survival data (CPTAC has no ready-made OS_MONTHS/OS_STATUS —
# same derivation as CPTAC_p53_preproc.R / CPTAC_Q1_validation.R)
clin_raw <- read.delim(
  "gbm_cptac_2021/data_clinical_patient.txt",
  sep = "\t", skip = 4, stringsAsFactors = FALSE
)
clin_raw$PATIENT_ID  <- trimws(clin_raw$PATIENT_ID)
clin_raw$SURV_STATUS <- as.integer(clin_raw$VITAL_STATUS == "Deceased")
clin_raw$SURV_TIME   <- ifelse(
  clin_raw$SURV_STATUS == 1,
  suppressWarnings(as.numeric(clin_raw$PATH_DIAG_TO_DEATH_DAYS)) / 30.44,
  suppressWarnings(as.numeric(clin_raw$PATH_DIAG_TO_LAST_CONTACT_DAYS)) / 30.44
)
clin_raw$AGE      <- suppressWarnings(as.numeric(clin_raw$AGE))
clin_raw$SEX_MALE <- as.integer(clin_raw$SEX == "Male")
clin_raw <- clin_raw[!is.na(clin_raw$SURV_TIME) & clin_raw$SURV_TIME > 0 &
                      !is.na(clin_raw$AGE), ]
clin <- clin_raw[, c("PATIENT_ID", "SURV_TIME", "SURV_STATUS", "AGE", "SEX_MALE")]
cat(sprintf("CPTAC-GBM patients with OS + age + sex: %d\n", nrow(clin)))

# 2. Pathway gene expression (3 files, 24 genes total)
p53_genes   <- read.csv("CPTAC_patient_p53genes.csv", stringsAsFactors = FALSE)
rbe2f_genes <- read.csv("CPTAC_patient_RbE2Fgenes.csv", stringsAsFactors = FALSE)
pappa_genes <- read.csv("CPTAC_patient_Pappalardogenes.csv", stringsAsFactors = FALSE)

p53_genes   <- p53_genes[,   !colnames(p53_genes)   %in% "SAMPLE_ID"]
rbe2f_genes <- rbe2f_genes[, !colnames(rbe2f_genes) %in% "SAMPLE_ID"]
pappa_genes <- pappa_genes[, !colnames(pappa_genes) %in% "SAMPLE_ID"]

# 3. p53 ODE steady-state outputs (20 features) — same DDR grid as TCGA
# (both notebooks use np.linspace(0.01, 1, num=10); verified when this
# feature was first validated in CPTAC_Q1_validation.R)
p53s15 <- read.csv("CPTAC_p53s15DR_patients.csv", stringsAsFactors = FALSE)
p53s46 <- read.csv("CPTAC_p53s46DR_patients.csv", stringsAsFactors = FALSE)
ddr_cols <- grep("^DDR_", colnames(p53s15), value = TRUE)
p53s15 <- p53s15[, c(ddr_cols, "PATIENT_ID")]
p53s46 <- p53s46[, c(ddr_cols, "PATIENT_ID")]
colnames(p53s15)[colnames(p53s15) %in% ddr_cols] <- paste0("s15_", ddr_cols)
colnames(p53s46)[colnames(p53s46) %in% ddr_cols] <- paste0("s46_", ddr_cols)

# 4. IDH mutation status
idh <- read.csv("Q1_CPTAC_IDH_status.csv", stringsAsFactors = FALSE)

# 4b. Assignment 2 DepMap/GDSC2 TMZ-sensitivity score (see
# Q1_extract_depmap_signature.R)
depmap <- read.csv("Q1_CPTAC_depmap_score.csv", stringsAsFactors = FALSE)

# 5. Merge (inner join, same logic as the TCGA training script)
feat <- clin %>%
  inner_join(p53_genes,   by = "PATIENT_ID") %>%
  inner_join(rbe2f_genes, by = "PATIENT_ID") %>%
  inner_join(pappa_genes, by = "PATIENT_ID") %>%
  inner_join(p53s15,      by = "PATIENT_ID") %>%
  inner_join(p53s46,      by = "PATIENT_ID") %>%
  inner_join(idh,         by = "PATIENT_ID") %>%
  inner_join(depmap,      by = "PATIENT_ID")

cat(sprintf("CPTAC-GBM patients with ALL feature sources present: %d\n", nrow(feat)))

# 6. Load the FROZEN signature from TCGA — no re-fitting past this point
coef_df <- read.csv("Q1_signature_coefficients.csv", stringsAsFactors = FALSE)
cat(sprintf("\nFrozen signature loaded: %d features (from Q1_build_signature_TCGA.R)\n", nrow(coef_df)))
print(coef_df)

missing_feats <- setdiff(coef_df$feature, colnames(feat))
if (length(missing_feats) > 0)
  stop("Frozen signature needs features not present in CPTAC data: ",
       paste(missing_feats, collapse = ", "))

X <- as.matrix(feat[, coef_df$feature, drop = FALSE])
storage.mode(X) <- "double"
y <- Surv(feat$SURV_TIME, feat$SURV_STATUS)

# 7. Apply the frozen coefficients — this is the actual cross-dataset test
risk_score <- as.numeric(X %*% coef_df$coefficient)

fit_val <- coxph(y ~ risk_score)
c_index <- summary(fit_val)$concordance[1]
c_index_se <- summary(fit_val)$concordance[2]
hr <- summary(fit_val)$coefficients[, "exp(coef)"]
p_val <- summary(fit_val)$coefficients[, "Pr(>|z|)"]

cat("\n--- Q1 result: frozen TCGA signature applied to CPTAC-GBM ---\n")
cat(sprintf("CPTAC patients tested:  %d\n", nrow(feat)))
cat(sprintf("Harrell's C-index:      %.3f (SE %.3f)\n", c_index, c_index_se))
cat(sprintf("Hazard ratio:           %.3f\n", hr))
cat(sprintf("Cox p-value:            %.4g\n", p_val))
cat(sprintf("Generalises (p<0.05):   %s\n", ifelse(p_val < 0.05, "YES", "NO")))
cat("(C-index of 0.5 = no better than chance; 1.0 = perfect ranking of who survives longer.)\n")

# 7b. Bootstrap confidence interval on the C-index. A single point estimate
# doesn't say how much it would move on a different sample of CPTAC
# patients; resampling patients with replacement and recomputing the C-index
# each time gives a proper interval instead.
set.seed(42)
n_boot <- 1000
boot_c <- numeric(n_boot)
n_pts  <- nrow(feat)
for (b in seq_len(n_boot)) {
  idx <- sample(seq_len(n_pts), n_pts, replace = TRUE)
  boot_fit <- tryCatch(
    # reverse = TRUE: risk_score is a hazard score (higher = worse survival,
    # same direction as HR above). Bare-formula concordance() defaults to the
    # OPPOSITE convention (higher predictor = better outcome) unless told
    # otherwise -- coxph's own $concordance auto-orients from the fitted
    # coefficient sign, but this direct formula call does not.
    survival::concordance(Surv(feat$SURV_TIME[idx], feat$SURV_STATUS[idx]) ~ risk_score[idx], reverse = TRUE)$concordance,
    error = function(e) NA
  )
  boot_c[b] <- boot_fit
}
boot_c <- boot_c[!is.na(boot_c)]
ci_lo <- quantile(boot_c, 0.025)
ci_hi <- quantile(boot_c, 0.975)
cat(sprintf("Bootstrap 95%% CI (n=%d resamples): C-index %.3f - %.3f\n", length(boot_c), ci_lo, ci_hi))
cat(ifelse(ci_lo > 0.5,
           "CI excludes 0.5 -- reasonably confident this beats chance, not just a lucky point estimate.\n",
           "CI includes 0.5 -- can't rule out this is within noise of chance-level performance.\n"))

p_forest <- ggforest(fit_val, data = data.frame(risk_score),
                     main = "Q1 signature (frozen from TCGA) -> CPTAC-GBM OS")
png(file.path(FIG_DIR, "fig_Q1_signature_CPTAC_cox.png"), width = 1400, height = 500, res = 150)
print(p_forest)
dev.off()

# Median-split KM for visualisation only (median of CPTAC's own risk score
# distribution — the Cox result above, on the continuous score, is the
# primary claim; this plot is not threshold-optimised for significance)
med <- median(risk_score)
grp <- ifelse(risk_score > med, paste0("High risk score (n=", sum(risk_score > med), ")"),
                                 paste0("Low risk score (n=",  sum(risk_score <= med), ")"))
km_df <- data.frame(grp = grp)
km_fit <- survfit(y ~ grp, data = km_df)
p_km <- ggsurvplot(
  km_fit, data = km_df, pval = TRUE,
  xlab = "Overall Survival (months)", ylab = "Survival Probability",
  palette = c("#E64B35", "#4DBBD5"),
  title = "Q1 combined signature (frozen from TCGA) — CPTAC-GBM validation",
  legend.title = "", risk.table = TRUE, risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 13)
)
png(file.path(FIG_DIR, "fig_Q1_signature_CPTAC_km.png"), width = 1400, height = 1100, res = 150)
print(p_km)
dev.off()

cat("\nDone. Figures saved: fig_Q1_signature_CPTAC_cox.png, fig_Q1_signature_CPTAC_km.png\n")
cat("This C-index/p-value is the answer to Q1 for this signature — whatever it says.\n")
