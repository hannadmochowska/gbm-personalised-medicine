# GBM Personalised Medicine: Cross-Dataset Response Prediction and Longitudinal ODE Modelling

This analysis evaluates whether radiographic response in glioblastoma (GBM) can be predicted across independent datasets and whether a mechanistic ordinary differential equation (ODE) provides a useful representation of longitudinal tumour behaviour.

## Questions addressed

- **Q1 — can a predictor of patient response work across datasets?** A ridge
  logistic-regression model was developed in the Burdenko Glioblastoma
  Progression Dataset and transferred without refitting to the independent
  LUMIERE cohort. The model used age, sex, MGMT and IDH because these were the
  clinical variables available in both datasets.
- **Q3 — are there any useful ODE models?** A patient-specific exponential ODE
  was fitted to longitudinal LUMIERE tumour volumes. Its usefulness was assessed
  through model fit and the association between its estimated weekly growth rate
  and the first evaluable expert RANO response.

## Overall story

The two analyses reach a consistent conclusion. Static clinical variables did
not produce a transportable response classifier, and a constant exponential
growth rate did not capture the complexity of treated longitudinal tumour
trajectories. The analyses nevertheless establish reproducible baselines and
identify the additional information required: harmonised radiotherapy features,
reliable treatment timing and treatment-aware longitudinal models.

## Result summary

### Q1: cross-dataset clinical response model

| Evaluation | Patients | ROC AUC | Sensitivity | Specificity | Balanced accuracy |
|---|---:|---:|---:|---:|---:|
| Burdenko held-out cross-validation | 178 | 0.539 | 0.952 | 0.115 | 0.534 |
| LUMIERE external validation | 79 | 0.558 | 1.000 | 0.000 | 0.500 |

The external AUC 95% confidence interval was **0.418–0.699**. The
development-selected threshold classified every evaluable LUMIERE patient as
progression, producing zero specificity. An age-and-sex-only sensitivity model
also failed to provide useful external classification (AUC 0.594; balanced
accuracy 0.436).

The honest answer to Q1 is therefore that a cross-dataset pipeline can be built
and tested, but the available clinical variables do **not** support a useful
generalisable predictor.

### Q3: longitudinal exponential ODE

The fitted model was:

```text
dV/dt = rV
V(t) = V0 exp(rt)
```

where `V` is total enhancing plus non-enhancing segmented tumour volume and `r`
is the patient-specific weekly net growth rate.

| Result | Estimate |
|---|---:|
| Patients with an estimable ODE | 77 |
| Median patient-level R-squared | 0.266 |
| Growth-rate ROC AUC for RANO progression | 0.494 |
| AUC 95% confidence interval | 0.357–0.630 |
| Progression-group comparison p-value | 0.933 |

The exponential ODE offers a compact descriptive growth-rate parameter, but it
explains only a modest proportion of within-patient volume variation and does
not distinguish progression from non-progression. The answer to Q3 is that this
simple ODE has **limited descriptive value and is not useful as a response
predictor in its present form**.

## Datasets

| Dataset | Role | Data used | Patients used |
|---|---|---|---:|
| Burdenko-GBM-Progression | Q1 development | Clinical variables and first follow-up response | 178 |
| LUMIERE | Q1 external validation | Clinical variables and first evaluable expert RANO response | 79 |
| LUMIERE | Q3 longitudinal modelling | Automated enhancing/non-enhancing volumes and expert RANO response | 77 |

Burdenko RTPLAN and RTSTRUCT files were identified and filtered, but the DICOM
files require NIH controlled access. The downloaded manifest is retained in the
project; imaging-derived dose and target-volume features are described as future
work rather than completed results.

## Reports

```text
reports/Q1_radiotherapy_response_model.qmd  Q1 analysis
reports/Q3_LUMIERE_ODE_model.qmd             Q3 analysis
```

Both reports and rendered HTML files are stored in `reports/`.

## How to run

1. Open the relevant `.qmd` file in RStudio.
2. Confirm that the source CSV files remain in the project directory.
3. Select **Render**.
4. The report checks for required files and R packages before running.

Run Q1 first if the reports are being presented as one project story. Q3 is
standalone but explicitly explains how longitudinal modelling follows the Q1
transportability result.

## Dependencies

- R: `dplyr`, `tidyr`, `readr`, `stringr`, `janitor`, `ggplot2`, `knitr`,
  `glmnet`, `pROC`
- Quarto and RStudio for rendering

## Main limitations

- Burdenko and LUMIERE use different response constructions.
- Progression is substantially less common in the Burdenko development cohort.
- MGMT and IDH contain considerable missingness.
- The controlled Burdenko DICOM files were unavailable within the project
  timeframe.
- The ODE uses nominal study weeks and combines treatment, surgery and natural
  growth into one net rate.
- The ODE analysis is exploratory and does not constitute external validation.

## Interpretation

Negative external validation is still an informative result. The project does
not claim that cross-dataset prediction or mathematical modelling is impossible.
It shows that static clinical variables and a constant exponential ODE are
insufficient for this task. A stronger next model would combine clinical
variables, radiation dose and target volumes with a piecewise or
treatment-effect ODE informed by reliable treatment dates.

