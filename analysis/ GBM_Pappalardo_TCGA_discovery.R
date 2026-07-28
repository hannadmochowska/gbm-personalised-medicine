# TCGA-GBM Pappalardo (RTK/RAS/MAPK + PI3K/AKT/mTOR) ODE Discovery — Question 1 (new pathway)
#
# Purpose: Discovery-only scan of the Pappalardo RTK/PI3K/MAPK ODE model's
# candidate features against TCGA-GBM survival. Mirrors GBM_RbE2F_TCGA_discovery.R
# and the Q1 block of GBM_Analysis.R (p53), applied to this third pathway.
#
# Prerequisites: Run GBM_Pappalardo_preproc.R -> GBM_Pappalardo_model.ipynb first.

library(tidyverse)
library(survival)
library(survminer)
library(stringr)

FIG_DIR <- getwd()

cat("=== TCGA-GBM: Pappalardo RTK/PI3K/MAPK ODE discovery scan ===\n")

# 1. TCGA clinical survival data (same as GBM_Analysis.R / GBM_RbE2F_TCGA_discovery.R)
clin_raw <- read.delim(
  "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
  sep = "\t", skip = 4, stringsAsFactors = FALSE
)
rownames(clin_raw)   <- clin_raw$PATIENT_ID
clin_raw$SURV_TIME   <- as.numeric(clin_raw$OS_MONTHS)
clin_raw$SURV_STATUS <- as.integer(clin_raw$OS_STATUS == "1:DECEASED")
clin_raw <- clin_raw[!is.na(clin_raw$SURV_TIME) & clin_raw$SURV_TIME > 0, ]
cat(sprintf("Clinical patients with OS data: %d\n", nrow(clin_raw)))

# 2. Load Pappalardo ODE outputs
ERK_df <- read.csv("GBM_Pappalardo_ERK_patients.csv", stringsAsFactors = FALSE)
AKT_df <- read.csv("GBM_Pappalardo_AKT_patients.csv", stringsAsFactors = FALSE)

gf_cols <- grep("^GF_", colnames(ERK_df), value = TRUE)
cat("Growth-factor levels modelled:", length(gf_cols), "\n")

# Merge with clinical data (trim sample ID to 12-char patient ID)
ERK_df$PATIENT_ID_12 <- str_sub(ERK_df$PATIENT_ID, end = 12)
AKT_df$PATIENT_ID_12 <- str_sub(AKT_df$PATIENT_ID, end = 12)

common_ids <- Reduce(intersect, list(
  rownames(clin_raw),
  ERK_df$PATIENT_ID_12,
  AKT_df$PATIENT_ID_12
))
cat("Patients with ODE output + clinical data:", length(common_ids), "\n")

clin <- clin_raw[common_ids, ]
erk  <- ERK_df[match(common_ids, ERK_df$PATIENT_ID_12), gf_cols]
akt  <- AKT_df[match(common_ids, AKT_df$PATIENT_ID_12), gf_cols]
colnames(erk) <- paste0("ERK_", gf_cols)
colnames(akt) <- paste0("AKT_", gf_cols)

surv_obj <- Surv(clin$SURV_TIME, clin$SURV_STATUS)

# 3. Univariate Cox: scan all 20 features (10 GF levels x 2 readouts)
all_features  <- cbind(erk, akt)
feature_names <- colnames(all_features)

# Guard against degenerate/near-constant features before fitting Cox models —
# the Rb-E2F discovery run showed that a feature with almost no variance across
# patients (most values piled at one extreme) produces a numerically unstable
# fit (HR collapsing to ~0 with a nonsensical p-value) that can spuriously "win"
# a p-value ranking. Flag any feature whose values aren't reasonably spread out.
feature_ok <- sapply(all_features, function(x) {
  spread <- diff(quantile(x, c(0.1, 0.9), na.rm = TRUE))
  spread > 1e-6 * max(abs(x), 1)
})
if (any(!feature_ok)) {
  cat("Excluding degenerate (near-constant) features from ranking:\n")
  cat(" ", paste(feature_names[!feature_ok], collapse = ", "), "\n")
}
feature_names <- feature_names[feature_ok]

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

cat("\nTop 10 Pappalardo ODE features by Cox p-value:\n")
print(head(results_df, 10))

# Best single-feature Cox + KM
best_feat <- results_df$feature[1]
cat(sprintf("\nBest feature: %s (HR=%.3f, p=%.4f)\n",
            best_feat, results_df$HR[1], results_df$p_value[1]))

best_vals <- all_features[[best_feat]]
best_name <- ifelse(grepl("^ERK", best_feat), "ERK activity", "AKT activity")
best_gf   <- str_extract(best_feat, "GF_[0-9.e+]+")

# Cox forest plot for best feature
best_df  <- data.frame(surv_obj, activity = best_vals)
fit_best <- coxph(surv_obj ~ activity, data = best_df)
p_forest <- ggforest(fit_best, data = best_df,
                     main = sprintf("Cox: %s (%s) -> GBM Overall Survival [DISCOVERY]", best_name, best_gf))
png(file.path(FIG_DIR, "fig_Pappalardo_TCGA_cox_best.png"), width = 1400, height = 500, res = 150)
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
if (min(n_high, n_low) < 15) {
  cat("CAUTION: optimal-threshold split is very lopsided (<15 in the smaller\n")
  cat("group) — treat this KM p-value as illustrative, not confirmatory, same\n")
  cat("issue flagged in the CPTAC p53 validation and the Rb-E2F discovery run.\n")
}

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
  title                 = sprintf("Pappalardo RTK/PI3K/MAPK ODE Model — GBM Survival [DISCOVERY]\n(%s at %s)", best_name, best_gf),
  legend.title          = best_name,
  legend.labs           = c(paste0("High (n=", n_high, ")"), paste0("Low (n=", n_low, ")")),
  risk.table            = TRUE,
  risk.table.height     = 0.25,
  risk.table.y.text     = FALSE,
  risk.table.y.text.col = TRUE,
  ggtheme               = theme_bw(base_size = 13)
)
png(file.path(FIG_DIR, "fig_Pappalardo_TCGA_km_best.png"), width = 1400, height = 1100, res = 150)
print(p_km)
dev.off()

# Multivariate Cox: ERK + AKT at the same GF level as the best feature
best_gf_tag <- str_extract(results_df$feature[1], "GF_[0-9.e+]+")
erk_col <- paste0("ERK_", best_gf_tag)
akt_col <- paste0("AKT_", best_gf_tag)

if (erk_col %in% colnames(all_features) && akt_col %in% colnames(all_features)) {
  multi_df  <- data.frame(ERK = all_features[[erk_col]], AKT = all_features[[akt_col]])
  fit_multi <- coxph(surv_obj ~ ERK + AKT, data = multi_df)
  cat("\nMultivariate Cox (ERK + AKT at best GF level):\n")
  print(summary(fit_multi)$coefficients)
  p_multi <- ggforest(fit_multi, data = multi_df,
                      main = sprintf("Multivariate Cox: ERK + AKT (%s) [DISCOVERY]", best_gf_tag))
  png(file.path(FIG_DIR, "fig_Pappalardo_TCGA_cox_multi.png"), width = 1400, height = 600, res = 150)
  print(p_multi)
  dev.off()
}

cat("\n--- Pappalardo TCGA Discovery Summary ---\n")
cat(sprintf("Best predictor:  %s\n", best_feat))
cat(sprintf("Hazard Ratio:    %.3f\n", results_df$HR[1]))
cat(sprintf("Cox p-value:     %.4f\n", results_df$p_value[1]))
cat(sprintf("KM log-rank p:   %.4f\n", opt_pval))
cat(sprintf("Significant:     %s\n", ifelse(results_df$p_value[1] < 0.05, "YES", "NO")))
cat("\nNOTE: this is a discovery-only scan (feature selected and threshold\n")
cat("optimised on this same TCGA cohort). Do NOT report this p-value as\n")
cat("confirmatory. If the result looks promising, the next step is to test\n")
cat("ONLY this pre-specified feature (not a new scan) in CPTAC.\n")
