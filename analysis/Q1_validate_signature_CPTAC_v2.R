# Q1 v2 — Validate the Pathway-Level Signature: CPTAC-GBM (Independent Test)
#
# Purpose: Apply the frozen v2 signature from Q1_build_signature_TCGA_v2.R
# (Q1_signature_coefficients_v2.csv) to CPTAC-GBM UNCHANGED. Same discipline
# as v1: no re-fitting, no re-scanning. Uses "_v2" filenames throughout so
# this never overwrites v1's headline result files.
#
# Prerequisites (run first, in this order):
#   CPTAC_p53_preproc.R, CPTAC_RbE2F_preproc.R, CPTAC_Pappalardo_preproc.R,
#   CPTAC_p53_model.ipynb (produces CPTAC_p53s15DR_patients.csv / s46),
#   Q1_extract_IDH_status.R, Q1_build_signature_TCGA_v2.R (produces
#   Q1_signature_coefficients_v2.csv)

library(tidyverse)
library(survival)
library(survminer)
library(stringr)

FIG_DIR <- getwd()

cat("=== Q1 v2: Validating pathway-level signature — CPTAC-GBM (independent test) ===\n\n")

# 1. Clinical / survival data
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

# 2. Pathway gene expression -> one mean score per pathway
p53_genes   <- read.csv("CPTAC_patient_p53genes.csv", stringsAsFactors = FALSE)
rbe2f_genes <- read.csv("CPTAC_patient_RbE2Fgenes.csv", stringsAsFactors = FALSE)
pappa_genes <- read.csv("CPTAC_patient_Pappalardogenes.csv", stringsAsFactors = FALSE)

P53_GENE_COLS   <- c("ATM", "CHEK2", "HIPK2", "MDM2", "PPM1D", "SIAH1", "TP53", "WSB1")
RBE2F_GENE_COLS <- c("MYC", "CCND1", "CCNE1", "RB1", "E2F1")
PAPPA_GENE_COLS <- c("EGFR", "PTEN", "SOS1", "KRAS", "RAF1", "MAP2K1", "MAPK1", "PIK3CA", "AKT1", "MTOR", "RPS6KB1")

p53_score   <- data.frame(PATIENT_ID = p53_genes$PATIENT_ID,
                           P53_SCORE = rowMeans(p53_genes[, P53_GENE_COLS]))
rbe2f_score <- data.frame(PATIENT_ID = rbe2f_genes$PATIENT_ID,
                           RBE2F_SCORE = rowMeans(rbe2f_genes[, RBE2F_GENE_COLS]))
pappa_score <- data.frame(PATIENT_ID = pappa_genes$PATIENT_ID,
                           PAPPALARDO_SCORE = rowMeans(pappa_genes[, PAPPA_GENE_COLS]))

# 3. p53 ODE steady-state outputs -> mean across all 10 DDR levels
p53s15 <- read.csv("CPTAC_p53s15DR_patients.csv", stringsAsFactors = FALSE)
p53s46 <- read.csv("CPTAC_p53s46DR_patients.csv", stringsAsFactors = FALSE)
ddr_cols <- grep("^DDR_", colnames(p53s15), value = TRUE)
p53_ode_score <- data.frame(
  PATIENT_ID = p53s15$PATIENT_ID,
  P53_ODE_S15_MEAN = rowMeans(p53s15[, ddr_cols]),
  P53_ODE_S46_MEAN = rowMeans(p53s46[match(p53s15$PATIENT_ID, p53s46$PATIENT_ID), ddr_cols])
)

# 4. IDH mutation status
idh <- read.csv("Q1_CPTAC_IDH_status.csv", stringsAsFactors = FALSE)

# 5. Merge
feat <- clin %>%
  inner_join(p53_score,     by = "PATIENT_ID") %>%
  inner_join(rbe2f_score,   by = "PATIENT_ID") %>%
  inner_join(pappa_score,   by = "PATIENT_ID") %>%
  inner_join(p53_ode_score, by = "PATIENT_ID") %>%
  inner_join(idh,           by = "PATIENT_ID")

cat(sprintf("CPTAC-GBM patients with ALL feature sources present: %d\n", nrow(feat)))

# 6. Load the FROZEN v2 signature — no re-fitting past this point
coef_df <- read.csv("Q1_signature_coefficients_v2.csv", stringsAsFactors = FALSE)
cat(sprintf("\nFrozen v2 signature loaded: %d features (from Q1_build_signature_TCGA_v2.R)\n", nrow(coef_df)))
print(coef_df)

missing_feats <- setdiff(coef_df$feature, colnames(feat))
if (length(missing_feats) > 0)
  stop("Frozen signature needs features not present in CPTAC data: ",
       paste(missing_feats, collapse = ", "))

X <- as.matrix(feat[, coef_df$feature, drop = FALSE])
storage.mode(X) <- "double"
y <- Surv(feat$SURV_TIME, feat$SURV_STATUS)

# 7. Apply the frozen coefficients
risk_score <- as.numeric(X %*% coef_df$coefficient)

fit_val <- coxph(y ~ risk_score)
c_index    <- summary(fit_val)$concordance[1]
c_index_se <- summary(fit_val)$concordance[2]
hr    <- summary(fit_val)$coefficients[, "exp(coef)"]
p_val <- summary(fit_val)$coefficients[, "Pr(>|z|)"]

cat("\n--- Q1 v2 result: frozen pathway-level TCGA signature applied to CPTAC-GBM ---\n")
cat(sprintf("CPTAC patients tested:  %d\n", nrow(feat)))
cat(sprintf("Harrell's C-index:      %.3f (SE %.3f)\n", c_index, c_index_se))
cat(sprintf("Hazard ratio:           %.3f\n", hr))
cat(sprintf("Cox p-value:            %.4g\n", p_val))
cat(sprintf("Generalises (p<0.05):   %s\n", ifelse(p_val < 0.05, "YES", "NO")))
cat("(C-index of 0.5 = no better than chance; 1.0 = perfect ranking of who survives longer.)\n")

# 7b. Bootstrap 95% CI (reverse = TRUE from the start this time: risk_score
# is a hazard score, higher = worse, matching the HR direction above —
# bare-formula concordance() needs to be told that explicitly).
set.seed(42)
n_boot <- 1000
boot_c <- numeric(n_boot)
n_pts  <- nrow(feat)
for (b in seq_len(n_boot)) {
  idx <- sample(seq_len(n_pts), n_pts, replace = TRUE)
  boot_c[b] <- tryCatch(
    survival::concordance(Surv(feat$SURV_TIME[idx], feat$SURV_STATUS[idx]) ~ risk_score[idx], reverse = TRUE)$concordance,
    error = function(e) NA
  )
}
boot_c <- boot_c[!is.na(boot_c)]
ci_lo <- quantile(boot_c, 0.025)
ci_hi <- quantile(boot_c, 0.975)
cat(sprintf("Bootstrap 95%% CI (n=%d resamples): C-index %.3f - %.3f\n", length(boot_c), ci_lo, ci_hi))
cat(ifelse(ci_lo > 0.5,
           "CI excludes 0.5 -- reasonably confident this beats chance.\n",
           "CI includes 0.5 -- can't rule out chance-level performance.\n"))

p_forest <- ggforest(fit_val, data = data.frame(risk_score),
                     main = "Q1 v2 pathway-level signature (frozen from TCGA) -> CPTAC-GBM OS")
png(file.path(FIG_DIR, "fig_Q1_signature_v2_CPTAC_cox.png"), width = 1400, height = 500, res = 150)
print(p_forest)
dev.off()

med <- median(risk_score)
grp <- ifelse(risk_score > med, paste0("High risk score (n=", sum(risk_score > med), ")"),
                                 paste0("Low risk score (n=",  sum(risk_score <= med), ")"))
km_df <- data.frame(grp = grp)
km_fit <- survfit(y ~ grp, data = km_df)
p_km <- ggsurvplot(
  km_fit, data = km_df, pval = TRUE,
  xlab = "Overall Survival (months)", ylab = "Survival Probability",
  palette = c("#E64B35", "#4DBBD5"),
  title = "Q1 v2 pathway-level signature (frozen from TCGA) — CPTAC-GBM validation",
  legend.title = "", risk.table = TRUE, risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 13)
)
png(file.path(FIG_DIR, "fig_Q1_signature_v2_CPTAC_km.png"), width = 1400, height = 1100, res = 150)
print(p_km)
dev.off()

cat("\nDone. Figures saved: fig_Q1_signature_v2_CPTAC_cox.png, fig_Q1_signature_v2_CPTAC_km.png\n")
cat("This C-index/p-value is the v2 answer to Q1 — compare against v1's C-index = 0.587, p = 0.037.\n")
