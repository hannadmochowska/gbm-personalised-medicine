# GBM Personalised Medicine: Cross-Dataset Prediction of MGMT Promoter Methylation

Final project analysis for AI for Personalised Medicine. This work tests whether
a gene-expression classifier of an MGMT promoter methylation biomarker can be
developed in TCGA-GBM and transferred unchanged to CPTAC-GBM.

MGMT promoter methylation is associated with benefit from temozolomide, but it
is a molecular biomarker rather than a direct measurement of patient response.
The analysis therefore addresses the project’s cross-dataset prediction
question indirectly.

## Question addressed

- **Q1 — can a predictor of patient response work across datasets?** A LASSO
  gene-expression classifier was trained to predict a
  temozolomide-response-associated MGMT methylation outcome in TCGA and
  externally evaluated in CPTAC.

The precise study question was:

> Can a gene-expression predictor of an MGMT methylation biomarker associated
> with temozolomide response transfer from TCGA to CPTAC?

## Analysis design

1. Retain one tumour sample per patient.
2. Construct a comparable gene-level MGMT methylation score in both cohorts.
3. Derive the methylated/unmethylated threshold from TCGA only.
4. Retain genes measured in both cohorts.
5. Select the 500 most variable eligible genes using TCGA only.
6. Impute and standardise using TCGA parameters.
7. Generate held-out TCGA predictions through five-fold cross-validation.
8. Fit a final LASSO model in TCGA and freeze all coefficients and
   preprocessing.
9. Apply the frozen model once to CPTAC without refitting or selecting a new
   threshold.

## Result summary

| Evaluation | Patients | ROC AUC | Key interpretation |
|---|---:|---:|---|
| TCGA five-fold cross-validation | 51 | 0.610 | Weak internal discrimination |
| CPTAC external validation | 97 | 0.682 | Modest cross-dataset ranking ability |

The final TCGA LASSO retained **13 non-zero gene-expression features**.

At the prespecified probability threshold of 0.5:

- all 97 CPTAC patients were classified as methylated;
- sensitivity was 100%;
- specificity was 0%;
- apparent accuracy was 62.9%, reflecting methylated-class prevalence rather
  than effective binary classification.

The correct conclusion is that the model retains **modest cross-dataset
discrimination but does not generalise successfully as a calibrated binary
classifier**.

## Datasets

Both cohorts were downloaded from cBioPortal.

| Dataset | Role | Molecular inputs | Patients used |
|---|---|---|---:|
| TCGA-GBM PanCancer Atlas 2018 | Development and internal cross-validation | RNA-seq expression and HM450 methylation | 51 |
| CPTAC-GBM 2021 | Independent external validation | FPKM expression and EPIC methylation | 97 |

The matched sample sizes are smaller than the full clinical cohorts because
patients required both expression and MGMT-annotated methylation data.

## Important outcome qualification

The supplied TCGA methylation matrix was probe-level, whereas the supplied
CPTAC matrix was gene-indexed and did not retain CpG probe identifiers. The
analysis therefore averaged rows annotated to `MGMT` to construct a comparable
gene-level score.

This outcome is:

- study-specific;
- platform-dependent;
- not the clinical MGMT-STP27 classifier;
- not a direct measure of temozolomide response.

These qualifications must remain in any presentation or final submission.

## Report

The canonical corrected report is:

```text
/Users/saidatbello/Documents/AI in P - final project/
Q1_MGMT_cross_dataset_FIXED.qmd
```

Its rendered HTML is stored in the same directory.

## How to run

1. Open `Q1_MGMT_cross_dataset_FIXED.qmd` in RStudio.
2. Keep the TCGA and CPTAC clinical, expression and methylation files in the
   project directory using the filenames defined in the report.
3. Select **Render**.
4. The report checks that all required files and packages are present before
   loading the large molecular matrices.

## Dependencies

- R: `data.table`, `dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `janitor`,
  `glmnet`, `pROC`, `ggplot2`, `knitr`
- Quarto and RStudio for rendering

## Known strengths

- Independent external validation
- TCGA-only feature selection and preprocessing
- Frozen coefficients before CPTAC evaluation
- One tumour sample retained per patient
- Explicit compromise between feature completeness and cross-platform
  availability
- Honest reporting of discrimination and calibration failure

## Main limitations

- The target is an MGMT biomarker rather than observed response.
- Expression was quantified differently in TCGA and CPTAC.
- Methylation platforms and available annotations differ.
- The development cohort is small.
- The fixed threshold is poorly calibrated in CPTAC.
- Further validation against a clinically established MGMT assay is required.

## Interpretation

The analysis demonstrates why AUC alone is not enough. CPTAC AUC indicates that
the model contains some transferable ranking signal, but the fixed threshold
produces no useful separation into methylated and unmethylated groups. The
model should therefore be described as an exploratory cross-dataset biomarker
classifier, not a clinically deployable predictor of temozolomide response.

