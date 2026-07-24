# GBM Personalised Medicine: Are Mechanistic ODE Pathway Models Useful Predictors of Survival?

Final project for AI for Personalised Medicine. Tests whether three published, mechanistically grounded ODE models of glioblastoma (GBM) signalling pathways, parameterised from routine bulk RNA-seq, can predict patient overall survival.

Full write-up: `reports/GBM_ThreeModel_ODE_Report.docx`.

## Question addressed

- **Q3** (are there any useful ODE models?): three models tested, one each for DNA-damage response (p53), cell-cycle control (Rb-E2F), and growth-factor signalling (RTK/RAS/MAPK + PI3K/AKT/mTOR).
- **Q1** (does a predictor generalise across datasets?): the p53 model's TCGA discovery feature was tested, unchanged, in an independent cohort (CPTAC-GBM).

## Result summary

No model produced a validated, generalisable survival predictor.

| Model | Pathway | TCGA discovery | CPTAC validation |
|---|---|---|---|
| p53 | DNA-damage response | HR = 0.78, p = 0.065 (borderline) | HR = 0.98, p = 0.982 (did not replicate) |
| Rb-E2F | Cell-cycle restriction point | Clean null (all 20 features p > 0.3) | Not pursued, no discovery signal |
| RTK/PI3K/MAPK | Growth-factor signalling | Clean null (all 20 features HR ≈ 1) | Not pursued, no discovery signal |

See the report for the full discussion of why (MGMT/IDH status not modelled, bulk RNA-seq averaging over tumour heterogeneity, static snapshot vs. dynamic biology, small unstable Kaplan-Meier subgroups).

## Datasets

Both accessed via cBioPortal.

| Dataset | Study ID | Role | Patients used |
|---|---|---|---|
| TCGA-GBM (PanCancer Atlas 2018) | `gbm_tcga_pan_can_atlas_2018` | Discovery, all 3 models | 154 |
| CPTAC-GBM (Wang et al., Cancer Cell 2021) | `gbm_cptac_2021` | Validation, p53 model only | 96 |

Raw data files are not committed to this repo (too large). Download both studies from cBioPortal and place them in `data/gbm_tcga_pan_can_atlas_2018/` and `data/gbm_cptac_2021/` before running the pipeline.

## Repo structure

```
preprocessing/   R scripts: extract pathway genes per patient, normalise to relative expression
models/          Jupyter notebooks: solve each ODE model per patient, output candidate features
analysis/        R scripts: Cox regression / Kaplan-Meier discovery and validation
reports/         Final write-ups (.docx)
```

## How to run

Each pathway is a three-step pipeline: preprocess -> simulate -> analyse. Run in this order.

**p53 (TCGA discovery):**
1. `preprocessing/GBM_p53_preproc.R`
2. `models/GBM_p53_model.ipynb`
3. `analysis/GBM_Analysis.R`

**p53 (CPTAC validation, run after the TCGA step above):**
1. `preprocessing/CPTAC_p53_preproc.R`
2. `models/CPTAC_p53_model.ipynb`
3. `analysis/CPTAC_Q1_validation.R`

**Rb-E2F (TCGA discovery):**
1. `preprocessing/GBM_RbE2F_preproc.R`
2. `models/GBM_RbE2F_model.ipynb`
3. `analysis/GBM_RbE2F_TCGA_discovery.R`

**RTK/PI3K/MAPK (TCGA discovery):**
1. `preprocessing/GBM_Pappalardo_preproc.R`
2. `models/GBM_Pappalardo_model.ipynb`
3. `analysis/GBM_Pappalardo_TCGA_discovery.R`

## Dependencies

- R: `tidyverse`, `survival`, `survminer`, `stringr`
- Python: `numpy`, `pandas`, `scipy`, `matplotlib` (Jupyter notebook, no external ODE solver required)

## Model sources

- p53: Ma L, et al. PNAS. 2005;102(40):14266-14271.
- Rb-E2F: Yao G, et al. Nat Cell Biol. 2008;10(4):476-482.
- RTK/PI3K/MAPK: Pappalardo F, et al. PLoS ONE. 2016;11(3):e0152104.
