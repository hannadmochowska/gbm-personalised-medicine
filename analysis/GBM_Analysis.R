# GBM p53 ODE Model — Survival Analysis
#
# Prerequisites: Run GBM_p53_preproc.R → GBM_p53_model.ipynb first.

# 0. Packages & shared settings
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
library(stringr)

FIG_DIR <- getwd()

# 1. Clinical survival data (TCGA-GBM)
cat("=== Loading clinical data ===\n")

clin_raw <- read.delim(
  "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
  sep = "\t", skip = 4, stringsAsFactors = FALSE
)
rownames(clin_raw)    <- clin_raw$PATIENT_ID
clin_raw$SURV_TIME    <- as.numeric(clin_raw$OS_MONTHS)
clin_raw$SURV_STATUS  <- as.integer(clin_raw$OS_STATUS == "1:DECEASED")
clin_raw              <- clin_raw[!is.na(clin_raw$SURV_TIME) & clin_raw$SURV_TIME > 0, ]
cat(sprintf("Clinical patients with OS data: %d\n", nrow(clin_raw)))

# 2. Does the p53 ODE model predict patient survival in TCGA-GBM?
cat("\n\n=== Does the p53 ODE model predict patient survival? ===\n")

# Load p53 ODE outputs for patients
p53s15 <- read.csv("GBM_p53s15DR_patients.csv", stringsAsFactors = FALSE)
p53s46 <- read.csv("GBM_p53s46DR_patients.csv", stringsAsFactors = FALSE)

ddr_cols <- grep("^DDR_", colnames(p53s15), value = TRUE)
cat("DDR levels modelled:", length(ddr_cols), "\n")

# Merge with clinical data (trim sample ID to 12-char patient ID)
p53s15$PATIENT_ID_12 <- str_sub(p53s15$PATIENT_ID, end = 12)
p53s46$PATIENT_ID_12 <- str_sub(p53s46$PATIENT_ID, end = 12)

# Deduplicate to one primary tumor sample per patient ID
p53s15 <- p53s15[!duplicated(p53s15$PATIENT_ID_12), ]
p53s46 <- p53s46[!duplicated(p53s46$PATIENT_ID_12), ]

common_ids <- Reduce(intersect, list(
  rownames(clin_raw),
  p53s15$PATIENT_ID_12,
  p53s46$PATIENT_ID_12
))
cat("Patients with ODE output + clinical data:", length(common_ids), "\n")

clin <- clin_raw[common_ids, ]
s15  <- p53s15[match(common_ids, p53s15$PATIENT_ID_12), ddr_cols]
s46  <- p53s46[match(common_ids, p53s46$PATIENT_ID_12), ddr_cols]
colnames(s15) <- paste0("s15_", ddr_cols)
colnames(s46) <- paste0("s46_", ddr_cols)

surv_obj <- Surv(clin$SURV_TIME, clin$SURV_STATUS)

# Univariate Cox: scan all DDR levels, both p53s15 and p53s46
all_features  <- cbind(s15, s46)
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

cat("\nTop 10 p53 ODE features by Cox p-value:\n")
print(head(results_df, 10))

# Best single-feature Cox + KM
best_feat <- results_df$feature[1]
cat(sprintf("\nBest feature: %s (HR=%.3f, p=%.4f)\n",
            best_feat, results_df$HR[1], results_df$p_value[1]))

best_vals <- all_features[[best_feat]]
best_name <- ifelse(grepl("^s15", best_feat), "p53-Ser15", "p53-Ser46")
best_ddr  <- str_extract(best_feat, "DDR_[0-9.]+")

# Cox forest plot for best feature
best_df  <- data.frame(surv_obj, p53_activity = best_vals)
fit_best <- coxph(surv_obj ~ p53_activity, data = best_df)
p_forest <- ggforest(fit_best, data = best_df,
                     main = sprintf("Cox: %s (%s) → GBM Overall Survival", best_name, best_ddr))
png(file.path(FIG_DIR, "fig_Q3_cox_p53_best.png"), width = 1400, height = 500, res = 150)
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
  title                 = sprintf("p53 ODE Model Predicts GBM Survival\n(%s at %s)", best_name, best_ddr),
  legend.title          = best_name,
  legend.labs           = c(paste0("High (n=", n_high, ")"), paste0("Low (n=", n_low, ")")),
  risk.table            = TRUE,
  risk.table.height     = 0.25,
  risk.table.y.text     = FALSE,
  risk.table.y.text.col = TRUE,
  ggtheme               = theme_bw(base_size = 13)
)
png(file.path(FIG_DIR, "fig_Q3_km_p53_best.png"), width = 1400, height = 1100, res = 150)
print(p_km)
dev.off()

# Multivariate Cox: p53s15 + p53s46 at best DDR level
best_ddr_tag <- str_extract(results_df$feature[1], "DDR_[0-9.]+")
s15_col      <- paste0("s15_", best_ddr_tag)
s46_col      <- paste0("s46_", best_ddr_tag)

if (s15_col %in% colnames(all_features) && s46_col %in% colnames(all_features)) {
  multi_df  <- data.frame(p53s15 = all_features[[s15_col]], p53s46 = all_features[[s46_col]])
  fit_multi <- coxph(surv_obj ~ p53s15 + p53s46, data = multi_df)
  cat("\nMultivariate Cox (p53s15 + p53s46 at best DDR):\n")
  print(summary(fit_multi)$coefficients)
  p_multi <- ggforest(fit_multi, data = multi_df,
                      main = sprintf("Multivariate Cox: p53s15 + p53s46 (%s)", best_ddr_tag))
  png(file.path(FIG_DIR, "fig_Q3_cox_p53_multi.png"), width = 1400, height = 600, res = 150)
  print(p_multi)
  dev.off()
}

cat("\n--- Summary ---\n")
cat(sprintf("Best predictor:  %s\n", best_feat))
cat(sprintf("Hazard Ratio:    %.3f\n", results_df$HR[1]))
cat(sprintf("Cox p-value:     %.4f\n", results_df$p_value[1]))
cat(sprintf("KM log-rank p:   %.4f\n", opt_pval))
cat(sprintf("Significant:     %s\n", ifelse(results_df$p_value[1] < 0.05, "YES", "NO")))
