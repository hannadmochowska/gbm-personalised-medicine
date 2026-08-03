# CPTAC-GBM Independent Validation — Question 1
#
# Purpose: Test whether the p53 ODE feature identified as the strongest
# survival predictor in TCGA-GBM (GBM_Analysis.R, Question 1 block —
# p53-Ser15 phosphorylation at DDR = 0.12) replicates in CPTAC-GBM, an
# independent cohort never used to fit the ODE model or select that feature.
#
# Prerequisites: Run CPTAC_p53_preproc.R -> CPTAC_p53_model.ipynb first.

library(tidyverse)
library(survival)
library(survminer)
library(stringr)

FIG_DIR <- getwd()

# Pre-specified feature from the TCGA discovery analysis (GBM_Analysis.R, Q1)
BEST_DDR_TAG     <- "DDR_0.120"
BEST_PHOSPHOSITE <- "s15"        # p53-Ser15
best_name        <- "p53-Ser15"

cat("=== CPTAC-GBM: Independent validation of p53 ODE feature ===\n")

# 1. CPTAC clinical data — derive survival columns
# (identical logic to CPTAC_p53_preproc.R; CPTAC has no ready-made OS_MONTHS/
# OS_STATUS like TCGA, so this is re-derived here rather than assumed.)
clin_raw <- read.delim(
  "gbm_cptac_2021/data_clinical_patient.txt",
  sep = "\t", skip = 4, stringsAsFactors = FALSE
)
clin_raw$PATIENT_ID <- trimws(clin_raw$PATIENT_ID)
clin_raw$SURV_STATUS <- as.integer(clin_raw$VITAL_STATUS == "Deceased")
clin_raw$SURV_TIME   <- ifelse(
  clin_raw$SURV_STATUS == 1,
  suppressWarnings(as.numeric(clin_raw$PATH_DIAG_TO_DEATH_DAYS)) / 30.44,
  suppressWarnings(as.numeric(clin_raw$PATH_DIAG_TO_LAST_CONTACT_DAYS)) / 30.44
)
clin_raw <- clin_raw[!is.na(clin_raw$SURV_TIME) & clin_raw$SURV_TIME > 0, ]
rownames(clin_raw) <- clin_raw$PATIENT_ID
cat(sprintf("CPTAC-GBM patients with OS data: %d\n", nrow(clin_raw)))

# 2. Load p53 ODE outputs (from CPTAC_p53_model.ipynb)
p53s15 <- read.csv("CPTAC_p53s15DR_patients.csv", stringsAsFactors = FALSE)
p53s46 <- read.csv("CPTAC_p53s46DR_patients.csv", stringsAsFactors = FALSE)

ddr_cols <- grep("^DDR_", colnames(p53s15), value = TRUE)
cat("DDR levels modelled:", length(ddr_cols), "\n")

# CPTAC RNA/sample IDs ARE patient IDs already (e.g. "C3L-00104") —
# unlike TCGA, no 12-character trimming is needed here.
common_ids <- Reduce(intersect, list(
  rownames(clin_raw),
  p53s15$PATIENT_ID,
  p53s46$PATIENT_ID
))
cat("Patients with ODE output + clinical data:", length(common_ids), "\n")

clin <- clin_raw[common_ids, ]
s15  <- p53s15[match(common_ids, p53s15$PATIENT_ID), ddr_cols]
s46  <- p53s46[match(common_ids, p53s46$PATIENT_ID), ddr_cols]
colnames(s15) <- paste0("s15_", ddr_cols)
colnames(s46) <- paste0("s46_", ddr_cols)
all_features <- cbind(s15, s46)

surv_obj <- Surv(clin$SURV_TIME, clin$SURV_STATUS)

# 3. PRIMARY: validate the pre-specified TCGA-best feature in CPTAC

best_feat <- paste0(BEST_PHOSPHOSITE, "_", BEST_DDR_TAG)
if (!best_feat %in% colnames(all_features))
  stop("Feature '", best_feat, "' not found in CPTAC ODE output — check the DDR grid matches the TCGA run (both use np.linspace(0.01, 1, num=10)).")
best_vals <- all_features[[best_feat]]

best_df  <- data.frame(surv_obj, p53_activity = best_vals)
fit_best <- coxph(surv_obj ~ p53_activity, data = best_df)
cat(sprintf("\nValidating TCGA-best feature (%s) in CPTAC-GBM:\n", best_feat))
print(summary(fit_best)$coefficients)

p_forest <- ggforest(fit_best, data = best_df,
                     main = sprintf("Cox: %s (%s) -> CPTAC-GBM OS (independent validation)", best_name, BEST_DDR_TAG))
png(file.path(FIG_DIR, "fig_CPTAC_cox_p53_best.png"), width = 1400, height = 500, res = 150)
print(p_forest)
dev.off()

# Kaplan-Meier: threshold optimised within CPTAC purely for the plot's visual
# split — the headline Cox p-value above already uses the continuous score,
# so this KM p-value is a secondary/illustrative statistic, not the primary claim.
uni_vals   <- sort(unique(best_vals))
scan_pvals <- rep(1, length(uni_vals))
for (i in 2:(length(uni_vals) - 1)) {
  grp <- as.integer(best_vals > uni_vals[i])
  if (length(unique(grp)) > 1)
    scan_pvals[i] <- survdiff(surv_obj ~ grp)$pvalue
}
opt_thr  <- uni_vals[which.min(scan_pvals)]
n_high   <- sum(best_vals > opt_thr)
n_low    <- sum(best_vals <= opt_thr)
opt_pval <- min(scan_pvals)
cat(sprintf("KM optimal threshold: %.4f  (n_high=%d, n_low=%d, p=%.4f)\n",
            opt_thr, n_high, n_low, opt_pval))

km_df  <- data.frame(
  group = ifelse(best_vals > opt_thr,
                 paste0("High (n=", n_high, ")"),
                 paste0("Low (n=",  n_low,  ")"))
)
KM_fit <- survfit(surv_obj ~ group, data = km_df)
p_km   <- ggsurvplot(
  KM_fit, data = km_df,
  pval                  = TRUE,
  xlab                  = "Overall Survival (months)",
  ylab                  = "Survival Probability",
  palette               = c("#E64B35", "#4DBBD5"),
  title                 = sprintf("p53 ODE Feature Validated in CPTAC-GBM\n(%s at %s)", best_name, BEST_DDR_TAG),
  legend.title          = best_name,
  legend.labs           = c(paste0("High (n=", n_high, ")"), paste0("Low (n=", n_low, ")")),
  risk.table            = TRUE,
  risk.table.height     = 0.25,
  risk.table.y.text     = FALSE,
  risk.table.y.text.col = TRUE,
  ggtheme               = theme_bw(base_size = 13)
)
png(file.path(FIG_DIR, "fig_CPTAC_km_p53_best.png"), width = 1400, height = 1100, res = 150)
print(p_km)
dev.off()

# Multivariate Cox: both phosphosites at the same DDR level (mirrors TCGA analysis)
s15_col <- paste0("s15_", BEST_DDR_TAG)
s46_col <- paste0("s46_", BEST_DDR_TAG)

multi_df  <- data.frame(p53s15 = all_features[[s15_col]], p53s46 = all_features[[s46_col]])
fit_multi <- coxph(surv_obj ~ p53s15 + p53s46, data = multi_df)
cat("\nMultivariate Cox (p53s15 + p53s46 at best DDR) in CPTAC-GBM:\n")
print(summary(fit_multi)$coefficients)

p_multi <- ggforest(fit_multi, data = multi_df,
                    main = sprintf("Multivariate Cox: p53s15 + p53s46 (%s) — CPTAC-GBM", BEST_DDR_TAG))
png(file.path(FIG_DIR, "fig_CPTAC_cox_p53_multi.png"), width = 1400, height = 600, res = 150)
print(p_multi)
dev.off()

cat("\n--- CPTAC Validation Summary (primary, pre-specified feature) ---\n")
cat(sprintf("Feature tested (from TCGA):  %s\n", best_feat))
cat(sprintf("CPTAC patients:              %d\n", nrow(clin)))
cat(sprintf("Hazard Ratio:                %.3f\n", summary(fit_best)$coefficients[, "exp(coef)"]))
cat(sprintf("Cox p-value:                 %.4f\n", summary(fit_best)$coefficients[, "Pr(>|z|)"]))
cat(sprintf("KM log-rank p:               %.4f\n", opt_pval))
cat(sprintf("Significant (Cox, p<0.05):   %s\n",
            ifelse(summary(fit_best)$coefficients[, "Pr(>|z|)"] < 0.05, "YES", "NO")))

# 4. SECONDARY / exploratory only: full 20-feature scan within CPTAC, for
# discussion context (e.g. "did the same feature also look best when CPTAC is
# scanned on its own?"). This is NOT the headline validation result — treating
# CPTAC's own top hit as "the" finding would be circular, since it would mean
# selecting the feature using the same data used to test it.

cat("\n--- Secondary/exploratory: full feature scan within CPTAC (context only) ---\n")
feature_names <- colnames(all_features)
cox_pvals <- sapply(feature_names, function(feat) {
  df <- data.frame(y = surv_obj, x = all_features[[feat]])
  summary(coxph(y ~ x, data = df))$coefficients[, "Pr(>|z|)"]
})
cox_hr <- sapply(feature_names, function(feat) {
  df <- data.frame(y = surv_obj, x = all_features[[feat]])
  summary(coxph(y ~ x, data = df))$coefficients[, "exp(coef)"]
})
results_df <- data.frame(feature = feature_names, HR = cox_hr, p_value = cox_pvals) %>%
  arrange(p_value)

cat("Top 5 CPTAC features by Cox p-value (exploratory):\n")
print(head(results_df, 5))
cat(sprintf("TCGA-best feature (%s) ranks in CPTAC's own top 5: %s\n",
            best_feat, ifelse(best_feat %in% head(results_df$feature, 5), "YES", "NO")))

cat("\nDone. Figures saved: fig_CPTAC_cox_p53_best.png, fig_CPTAC_km_p53_best.png, fig_CPTAC_cox_p53_multi.png\n")
