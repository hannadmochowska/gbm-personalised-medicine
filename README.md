# GBM Personalised Medicine: Are Mechanistic ODE Pathway Models Useful Predictors of Survival?

Final project for AI for Personalised Medicine. Tests whether three published, mechanistically grounded ODE models of glioblastoma (GBM) signalling pathways, parameterised from routine bulk RNA-seq, can predict patient overall survival — individually and combined into a single cross-dataset signature.

## Questions addressed

- **Q3** (are there any useful ODE models?): three models tested individually, one each for DNA-damage response (p53), cell-cycle control (Rb-E2F), and growth-factor signalling (RTK/RAS/MAPK + PI3K/AKT/mTOR). None predicted survival alone.
- **Q1** (does a predictor generalise across datasets?): two tests. (1) The p53 model's single TCGA discovery feature, tested unchanged in an independent cohort (CPTAC-GBM) — did not replicate. (2) A regularised signature pooling features from all three pathways plus age, cross-validated on TCGA-GBM and frozen, tested unchanged on CPTAC-GBM — **generalised**, with a real, non-trivial effect.
- **Q2 correlation check** (does the drug-viability predictor from Q2 correlate with Q1?): tested directly by correlating the two predictors' outputs patient by patient. No meaningful correlation in either cohort — the two predictors appear to capture different biology.

## Result summary

**Individual models (single best-scanned feature each):** no model produced a validated, generalisable survival predictor.

| Model | Pathway | TCGA discovery | CPTAC validation |
|---|---|---|---|
| p53 | DNA-damage response | HR = 0.78, p = 0.065 (borderline) | HR = 0.98, p = 0.982 (did not replicate) |
| Rb-E2F | Cell-cycle restriction point | Clean null (all 20 features p > 0.3) | Not pursued, no discovery signal |
| RTK/PI3K/MAPK | Growth-factor signalling | Clean null (all 20 features HR ≈ 1) | Not pursued, no discovery signal |

See the report for the full discussion of why (MGMT/IDH status not modelled, bulk RNA-seq averaging over tumour heterogeneity, static snapshot vs. dynamic biology, small unstable Kaplan-Meier subgroups).

**Combined cross-dataset signature (Q1):** pooling features across all three pathways, plus age, and letting a cross-validated penalised Cox model (LASSO vs. Elastic Net) select and weight them jointly — rather than hand-picking one feature — produces a signature that does generalise.

| Specification | Candidate features | CPTAC C-index | CPTAC Cox p | What survived regularisation |
|---|---|---|---|---|
| Individual genes + ODE outputs | 48 (24 pathway genes, 20 p53 ODE outputs, age/sex/IDH, DepMap TMZ score) | 0.589 | 0.022 | SIAH1, PTEN, MAPK1, one p53 ODE feature, age |
| Pathway-level mean scores | 8 (one mean score per pathway + age/sex/IDH) | 0.587 | 0.011 | age only |

Both specifications beat chance (bootstrap 95% CIs exclude 0.5), so the honest answer to Q1 is yes. But the two specifications reach it differently: individual-gene resolution finds real pathway biology (SIAH1, PTEN, MAPK1) alongside age; compressing each pathway into one mean score destroys that structure before the model ever sees it, leaving only age. This is the brief's "compromise between data availability and data completeness" playing out directly — full detail and interpretation in the report.

## Datasets

Both accessed via cBioPortal.

| Dataset | Study ID | Role | Patients used |
|---|---|---|---|
| TCGA-GBM (PanCancer Atlas 2018) | `gbm_tcga_pan_can_atlas_2018` | Discovery, all 3 models + Q1 signature | 154 (individual models) / 106 (Q1 signature, all feature sources + age required) |
| CPTAC-GBM (Wang et al., Cancer Cell 2021) | `gbm_cptac_2021` | Validation | 96 |

Patient counts differ between the individual-model analyses and the Q1 signature because the signature requires every candidate feature (all three pathways, ODE outputs, age, IDH, DepMap score) to be present for a patient, which is a stricter filter than any single model needs on its own.

Raw data files are not committed to this repo (too large). Download both studies from cBioPortal and place them in `data/gbm_tcga_pan_can_atlas_2018/` and `data/gbm_cptac_2021/` before running the pipeline.

## Known data-quality issue (fixed)

5 TCGA-GBM patients have two archived tumour samples each. An early version of the Q1 signature scripts merged per-gene feature files on patient ID without deduplicating first, which compounded across several joins and silently represented one patient's record 32 times in the training data — caught by independently re-deriving the merge in Python and finding the row count didn't match the number of unique patient IDs. All `Q1_*` and `Q2_Q1_correlation.R` scripts now deduplicate to one row per patient (the primary tumour sample) before joining. The individual-model discovery scripts (`GBM_Analysis.R`, `GBM_RbE2F_TCGA_discovery.R`, `GBM_Pappalardo_TCGA_discovery.R`) use the same raw per-gene files and merge pattern and have not yet been audited for the same issue — worth checking before treating their null results as final.

## Repo structure

```
preprocessing/   R scripts: extract pathway genes per patient, normalise to relative expression
models/          Jupyter notebooks: solve each ODE model per patient, output candidate features
analysis/        R scripts: Cox regression / Kaplan-Meier discovery, validation, and the Q1 combined signature
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

**Combined Q1 signature (run after all three pathways' TCGA preprocessing + model steps above, and the CPTAC steps for p53 and the CPTAC-only preprocessing for Rb-E2F/RTK-PI3K-MAPK):**
1. `analysis/Q1_extract_IDH_status.R` — derives IDH1/IDH2 mutation status for both cohorts
2. `analysis/Q1_extract_depmap_signature.R` — projects the Assignment 2 DepMap/GDSC2 TMZ-sensitivity signature onto both cohorts (needs `GBM_cellline_coefs.csv` from the Assignment 2 `DepMap_signature.R`, not included in this repo)
3. `analysis/Q1_build_signature_TCGA.R` — trains the individual-gene-level signature on TCGA-GBM, freezes coefficients
4. `analysis/Q1_validate_signature_CPTAC.R` — applies the frozen signature unchanged to CPTAC-GBM
5. `analysis/Q1_build_signature_TCGA_v2.R` / `analysis/Q1_validate_signature_CPTAC_v2.R` — same process using pathway-level mean scores instead of individual genes (the availability/completeness comparison)
6. `analysis/Q2_Q1_correlation.R` — correlates Q2's drug-viability score against the Q1 risk score directly, in both cohorts

## Dependencies

- R: `tidyverse`, `survival`, `survminer` (pulls in `ggpubr`), `glmnet`, `stringr`
- Python: `numpy`, `pandas`, `scipy`, `matplotlib` (Jupyter notebook, no external ODE solver required)

## Model sources

- p53: Ma L, et al. PNAS. 2005;102(40):14266-14271.
- Rb-E2F: Yao G, et al. Nat Cell Biol. 2008;10(4):476-482.
- RTK/PI3K/MAPK: Pappalardo F, et al. PLoS ONE. 2016;11(3):e0152104.
