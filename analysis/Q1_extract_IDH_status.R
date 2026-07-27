# Q1 — Extract IDH1/IDH2 Mutation Status (TCGA-GBM + CPTAC-GBM)
#
# Purpose: IDH mutation status is the single strongest, most-established
# prognostic marker in glioma (alongside MGMT methylation, which isn't used
# here — see note at the bottom of this script). It isn't a clinical
# attribute in either cBioPortal download, but both studies do ship a
# data_mutations.txt MAF file, so it's derived from there: a patient is
# "IDH-mutant" if their tumour carries a missense mutation in IDH1 or IDH2
# (captures the canonical GBM hotspots, e.g. IDH1 R132H, without hardcoding
# a specific amino-acid change — verified against the actual file that all
# IDH1 rows in both cohorts are Missense_Mutation).
#
# Output:
#   Q1_TCGA_IDH_status.csv  — PATIENT_ID, IDH_MUT (0/1)
#   Q1_CPTAC_IDH_status.csv — PATIENT_ID, IDH_MUT (0/1)
#
# Run this before Q1_build_signature_TCGA.R / Q1_validate_signature_CPTAC.R.

library(tidyverse)
library(stringr)

IDH_GENES <- c("IDH1", "IDH2")

extract_idh <- function(mutations_file, clin_file, trim_to_12, out_file, cohort_label) {
  cat(sprintf("=== %s ===\n", cohort_label))

  clin_raw <- read.delim(clin_file, sep = "\t", skip = 4, stringsAsFactors = FALSE)
  clin_raw$PATIENT_ID <- trimws(clin_raw$PATIENT_ID)
  all_patients <- unique(clin_raw$PATIENT_ID)

  mut_raw <- read.delim(mutations_file, sep = "\t", stringsAsFactors = FALSE, quote = "")
  idh_rows <- mut_raw[mut_raw$Hugo_Symbol %in% IDH_GENES &
                       mut_raw$Variant_Classification == "Missense_Mutation", ]

  mut_ids <- idh_rows$Tumor_Sample_Barcode
  if (trim_to_12) mut_ids <- str_sub(mut_ids, end = 12)  # TCGA: barcode -> 12-char patient ID
  mut_ids <- unique(mut_ids)

  cat(sprintf("  IDH1/IDH2 missense mutations found in %d patients (of %d total)\n",
              length(intersect(mut_ids, all_patients)), length(all_patients)))

  out <- data.frame(
    PATIENT_ID = all_patients,
    IDH_MUT    = as.integer(all_patients %in% mut_ids)
  )
  write.csv(out, out_file, row.names = FALSE)
  cat(sprintf("  Saved: %s\n\n", out_file))
  out
}

tcga_idh <- extract_idh(
  mutations_file = "gbm_tcga_pan_can_atlas_2018/data_mutations.txt",
  clin_file      = "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
  trim_to_12     = TRUE,
  out_file       = "Q1_TCGA_IDH_status.csv",
  cohort_label   = "TCGA-GBM"
)

cptac_idh <- extract_idh(
  mutations_file = "gbm_cptac_2021/data_mutations.txt",
  clin_file      = "gbm_cptac_2021/data_clinical_patient.txt",
  trim_to_12     = FALSE,  # CPTAC barcodes are already patient-level IDs
  out_file       = "Q1_CPTAC_IDH_status.csv",
  cohort_label   = "CPTAC-GBM"
)
