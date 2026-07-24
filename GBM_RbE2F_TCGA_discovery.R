# TCGA-GBM Rb-E2F ODE Discovery — Question 1 (new pathway)
#
# Purpose: Discovery-only scan of the Rb-E2F restriction-point ODE model's
# candidate features against TCGA-GBM survival. Mirrors the Q1 block of
# GBM_Analysis.R (the p53 model discovery analysis), applied to the new
# Rb-E2F pathway instead.
#
# This is DISCOVERY, not validation: it scans all 20 candidate features
# (10 serum-stimulus levels x 2 readouts: E2F, phospho-Rb) to find the best
# one in TCGA. Per the agreed plan, only if this shows something worth
# reporting should CPTAC validation be attempted next (mirroring how the p53
# CPTAC validation only tested TCGA's pre-specified best feature, never
# re-scanned CPTAC's own data).
#
# Prerequisites: Run GBM_RbE2F_preproc.R -> GBM_RbE2F_model.ipynb first.

library(tidyverse)
library(survival)
library(survminer)
library(stringr)

FIG_DIR <- getwd()

cat("=== TCGA-GBM: Rb-E2F ODE discovery scan ===\n")

# 1. TCGA clinical survival data (same as GBM_Analysis.R)
clin_raw <- read.delim(
  "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
  sep = "\t", skip = 4, stringsAsFactors = FALSE
)
rownames(clin_raw)   <- clin_raw$PATIENT_ID
clin_raw$SURV_TIME   <- as.numeric(clin_raw$OS_MONTHS)
clin_raw$SURV_STATUS <- as.integer(clin_raw$OS_STATUS == "1:DECEASED")
clin_raw <- clin_raw[!is.na(clin_raw$SURV_TIME) & clin_raw$SURV_TIME > 0, ]
cat(sprintf("Clinical patients with OS data: %d\n", nrow(clin_raw)))

# 2. Load Rb-E2F ODE outputs
EF_df <- read.csv("GBM_RbE2F_EF_patients.csv", stringsAsFactors = FALSE)
RP_df <- read.csv("GBM_RbE2F_RP_patients.csv", stringsAsFactors = FALSE)

s_cols <- grep("^S_", colnames(EF_df), value = TRUE)
cat("Serum-stimulus levels modelled:", length(s_cols), "\n")

# Merge with clinical data (trim sample ID to 12-char patient ID, as in GBM_Analysis.R)
EF_df$PATIENT_ID_12 <- str_sub(EF_df$PATIENT_ID, end = 12)
RP_df$PATIENT_ID_12 <- str_sub(RP_df$PATIENT_ID, end = 12)

common_ids <- Reduce(intersect, list(
  rownames(clin_raw),
  EF_df$PATIENT_ID_12,
  RP_df$PATIENT_ID_12
))
cat("Patients with ODE output + clinical data:", length(common_ids), "\n")

clin <- clin_raw[common_ids, ]
ef   <- EF_df[match(common_ids, EF_df$PATIENT_ID_12), s_cols]
rp   <- RP_df[match(common_ids, RP_df$PATIENT_ID_12), s_cols]
colnames(ef) <- paste0("EF_", s_cols)
colnames(rp) <- paste0("RP_", s_cols)

surv_obj <- Surv(clin$SURV_TIME, clin$SURV_STATUS)

# 3. Univariate Cox: scan all 20 features (10 S levels x 2 readouts)
all_features  <- cbind(ef, rp)
feature_names <- colnames(all_features)

cox_pvals <- sapply(feature_names, function(feat) {
  df  <- data.frame(y = surv_obj, x = all_features[[feat]])
  fit <- coxph(y ~ x, data = df)
  summary(fit)$coefficients[, "Pr(>|z|)"]
})
cox_hr <- sapply(feature_names, function(feat) {
  df  <- data.frame(y = surv_obj, x = all_features[[feat]])
  fit <- coxph(y ~ x, data = df)
  summary(fit)$coefficients[, "exp(coef)"]
})

results_df <- data.frame(
  feature = feature_names,
  HR      = cox_hr,
  p_value = cox_pvals
) %>% arrange(p_value)

cat("\nTop 10 Rb-E2F ODE features by Cox p-value:\n")
print(head(results_df, 10))

# Best single-feature Cox + KM
best_feat <- results_df$feature[1]
cat(sprintf("\nBest feature: %s (HR=%.3f, p=%.4f)\n",
            best_feat, results_df$HR[1], results_df$p_value[1]))

best_vals <- all_features[[best_feat]]
best_name <- ifelse(grepl("^EF", best_feat), "E2F activity", "Phospho-Rb")
best_s    <- str_extract(best_feat, "S_[0-9.]+")

# Cox forest plot for best feature
best_df  <- data.frame(surv_obj, rbe2f_activity = best_vals)
fit_best <- coxph(surv_obj ~ rbe2f_activity, data = best_df)
p_forest <- ggforest(fit_best, data = best_df,
                     main = sprintf("Cox: %s (%s) -> GBM Overall Survival [DISCOVERY]", best_name, best_s))
png(file.path(FIG_DIR, "fig_RbE2F_TCGA_cox_best.png"), width = 1400, height = 500, res = 150)
print(p_forest)
dev.off()

# Kaplan-Meier: optimal threshold scan
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
cat(sprintf("Optimal threshold: %.4f  (n_high=%d, n_low=%d, p=%.4f)\n",
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
  title                 = sprintf("Rb-E2F ODE Model — GBM Survival [DISCOVERY]\n(%s at %s)", best_name, best_s),
  legend.title          = best_name,
  legend.labs           = c(paste0("High (n=", n_high, ")"), paste0("Low (n=", n_low, ")")),
  risk.table            = TRUE,
  risk.table.height     = 0.25,
  risk.table.y.text     = FALSE,
  risk.table.y.text.col = TRUE,
  ggtheme               = theme_bw(base_size = 13)
)
png(file.path(FIG_DIR, "fig_RbE2F_TCGA_km_best.png"), width = 1400, height = 1100, res = 150)
print(p_km)
dev.off()

# Multivariate Cox: E2F + phospho-Rb at the same S level as the best feature
best_s_tag <- str_extract(results_df$feature[1], "S_[0-9.]+")
ef_col <- paste0("EF_", best_s_tag)
rp_col <- paste0("RP_", best_s_tag)

if (ef_col %in% colnames(all_features) && rp_col %in% colnames(all_features)) {
  multi_df  <- data.frame(E2F = all_features[[ef_col]], phosphoRb = all_features[[rp_col]])
  fit_multi <- coxph(surv_obj ~ E2F + phosphoRb, data = multi_df)
  cat("\nMultivariate Cox (E2F + phospho-Rb at best S level):\n")
  print(summary(fit_multi)$coefficients)
  p_multi <- ggforest(fit_multi, data = multi_df,
                      main = sprintf("Multivariate Cox: E2F + phospho-Rb (%s) [DISCOVERY]", best_s_tag))
  png(file.path(FIG_DIR, "fig_RbE2F_TCGA_cox_multi.png"), width = 1400, height = 600, res = 150)
  print(p_multi)
  dev.off()
}

cat("\n--- Rb-E2F TCGA Discovery Summary ---\n")
cat(sprintf("Best predictor:  %s\n", best_feat))
cat(sprintf("Hazard Ratio:    %.3f\n", results_df$HR[1]))
cat(sprintf("Cox p-value:     %.4f\n", results_df$p_value[1]))
cat(sprintf("KM log-rank p:   %.4f\n", opt_pval))
cat(sprintf("Significant:     %s\n", ifelse(results_df$p_value[1] < 0.05, "YES", "NO")))
cat("\nNOTE: this is a discovery-only scan (feature selected and threshold\n")
cat("optimised on this same TCGA cohort). Do NOT report this p-value as\n")
cat("confirmatory. If the result looks promising, the next step is to test\n")
cat("ONLY this pre-specified feature (not a new scan) in CPTAC — mirroring\n")
cat("CPTAC_Q1_validation.R for the p53 model.\n")
