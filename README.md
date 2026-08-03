# Glioblastoma Personalised Medicine: Mechanistic ODE Modeling, Cross-Dataset Validation & Target Discovery

[![Course Project](https://img.shields.io/badge/Course-AI_for_Personalised_Medicine-blue.svg)]()
[![Disease Focus](https://img.shields.io/badge/Disease-Glioblastoma_Multiforme_(GBM)-red.svg)]()

This repository contains the complete computational pipeline, ODE differential equation solvers, machine learning signatures, cross-dataset transportability evaluations, target discovery analyses, and Quarto reports for the **AI for Personalised Medicine** project on **Glioblastoma Multiforme (GBM)**.

---

## Project Executive Summary & Mapping to Course Questions

| Question | Core Topic | Computational / Analytical Strategy | Key Finding & Headline Result | Detailed Documentation |
| :---: | :--- | :--- | :--- | :--- |
| **Q1** | **Cross-Dataset Transportability** | Penalised Cox regression (LASSO) combining pathway features, gene expression, and clinical variables (TCGA $\to$ CPTAC). | **Signature Generalises ($C\text{-index} = 0.589, p = 0.022$)**. Individual gene resolution (*SIAH1*, *PTEN*, *MAPK1*) retains true pathway signal. | [`docs/Q1_Q3_pathway_signature.md`](docs/Q1_Q3_pathway_signature.md)<br>[`docs/Q1_MGMT_methylation.md`](docs/Q1_MGMT_methylation.md) |
| **Q2** | **Drug Viability Correlation Check** | Direct patient-by-patient correlation between Q2 in vitro DepMap TMZ score and Q1 patient overall survival risk. | **No Correlation ($r \approx 0.04, p > 0.50$)**. In vitro drug screens and in vivo survival capture orthogonal biology. | [`docs/Q1_Q3_pathway_signature.md`](docs/Q1_Q3_pathway_signature.md) |
| **Q3** | **Mechanistic ODE Pathway & Volume Models** | Parameterised 3 signaling ODE models (p53, Rb-E2F, RTK/PI3K/MAPK) from RNA-seq + 1 longitudinal MRI volume ODE ($dV/dt = rV$). | Standalone signaling ODE features failed to replicate alone ($HR = 0.98, p = 0.982$). Longitudinal growth rate $r$ failed RANO progression prediction ($AUC = 0.494$). | [`docs/Q1_Q3_pathway_signature.md`](docs/Q1_Q3_pathway_signature.md)<br>[`docs/Q1_Q3_radiotherapy_ODE.md`](docs/Q1_Q3_radiotherapy_ODE.md) |
| **Q4** | **Selective Target Discovery** | Genome-wide DepMap 23Q4 CRISPR knockout screen differential dependency scoring ($\Delta$) + DDR & LINCS L1000 mapping. | Nominated top novel candidates: **`CHMP4B`** ($\text{FDR} = 3.2\times 10^{-11}$), **`RPP25L`**, **`FERMT2`**, **`GPX4`**, **`ITGAV`**. Evaluated *SOX2*/*OLIG2* FDR limitations. | [`docs/Q4_target_discovery.md`](docs/Q4_target_discovery.md) |
| **Q5** | **Treatment Selection** | Clinical synthesis on prognostic vs. predictive biomarkers and platform calibration drift. | Delineated prognostic (Age, IDH, Q1 signature) vs. predictive (MGMT) markers. Highlighted threshold calibration shifts. | [`docs/Q5_treatment_selection.md`](docs/Q5_treatment_selection.md) |

---

## Datasets Overview

All datasets are publicly available from cBioPortal, TCGA, CPTAC, LUMIERE, or Broad DepMap:

| Dataset | Study ID / Source | Cohort Role | $N$ Patients | Data Modalities Used |
| :--- | :--- | :--- | :---: | :--- |
| **TCGA-GBM** | `gbm_tcga_pan_can_atlas_2018` | Discovery / Training | 154 / 106 | Bulk RNA-seq, HM450 Methylation, Survival, Age, Sex, IDH |
| **CPTAC-GBM** | `gbm_cptac_2021` | External Validation | 96 / 97 | Bulk RNA-seq (FPKM), EPIC Methylation, Survival, Age, Sex |
| **Burdenko-GBM** | Burdenko Progression Dataset | Radiotherapy Dev | 178 | Clinical variables (Age, Sex, MGMT, IDH) & Radiographic response |
| **LUMIERE** | LUMIERE MRI Cohort | Radiotherapy Val / ODE | 79 / 77 | Serial automated 3D MRI volume segmentations & Expert RANO |
| **DepMap 23Q4** | Broad Institute DepMap | CRISPR Target Discovery | 68 GBM / 1,032 Other | Genome-wide CRISPR-Cas9 Chronos knockout dependency scores |

---

## Repository Directory Structure

```
gbm-personalised-medicine/
├── README.md                        <- Primary entry point & project mapping table
│
├── preprocessing/                   <- R scripts: normalize bulk RNA-seq for ODE models
│   ├── GBM_p53_preproc.R            # Preprocess TCGA RNA-seq for p53 ODE model
│   ├── CPTAC_p53_preproc.R          # Preprocess CPTAC RNA-seq for p53 ODE model
│   ├── GBM_RbE2F_preproc.R          # Preprocess TCGA RNA-seq for Rb-E2F ODE model
│   ├── CPTAC_RbE2F_preproc.R        # Preprocess CPTAC RNA-seq for Rb-E2F ODE model
│   ├── GBM_Pappalardo_preproc.R     # Preprocess TCGA for RTK/PI3K/MAPK ODE model
│   └── CPTAC_Pappalardo_preproc.R   # Preprocess CPTAC for RTK/PI3K/MAPK ODE model
│
├── models/                          <- Jupyter notebooks: solve ODE systems per patient
│   ├── GBM_p53_model.ipynb          # Solve 20 p53 ODE state equations per TCGA patient
│   ├── CPTAC_p53_model.ipynb        # Solve p53 ODE state equations per CPTAC patient
│   ├── GBM_RbE2F_model.ipynb        # Solve 20 Rb-E2F ODE state equations per TCGA patient
│   ├── GBM_Pappalardo_model.ipynb   # Solve 20 RTK/PI3K/MAPK ODE equations per TCGA patient
│   └── Q4_GBM_DepMap_LINCS_Target_Discovery.ipynb # DepMap CRISPR target discovery notebook
│
├── analysis/                        <- R scripts: Cox regression, signatures, validation
│   ├── GBM_Analysis.R               # Discovery scan for p53 ODE features (TCGA)
│   ├── CPTAC_Q1_validation.R        # Independent validation of p53 ODE feature (CPTAC)
│   ├── GBM_RbE2F_TCGA_discovery.R   # Discovery scan for Rb-E2F ODE features
│   ├── GBM_Pappalardo_TCGA_discovery.R # Discovery scan for RTK/PI3K/MAPK ODE features
│   ├── Q1_extract_IDH_status.R      # Derive IDH1/IDH2 mutation status across cohorts
│   ├── Q1_extract_depmap_signature.R # Project Assignment 2 DepMap TMZ score onto cohorts
│   ├── Q1_build_signature_TCGA.R    # Train multi-pathway Cox LASSO signature on TCGA
│   ├── Q1_validate_signature_CPTAC.R # Validate frozen Q1 signature on CPTAC-GBM
│   ├── Q1_build_signature_TCGA_v2.R  # Train pathway-mean compressed signature
│   ├── Q1_validate_signature_CPTAC_v2.R # Validate pathway-mean signature
│   └── Q2_Q1_correlation.R          # Correlate Q2 drug score against Q1 survival risk
│
├── reports/                         <- Rendered Quarto reports and HTML deliverables
│   ├── Q1_MGMT_cross_dataset.qmd / .html
│   ├── Q1_radiotherapy_response_model.qmd / .html
│   └── Q3_LUMIERE_ODE_model.qmd / .html
│
├── figures/                         <- Output tables & publication figures for Q4
│   ├── gbm_top30_targets.csv
│   ├── gbm_selective_dependencies.csv
│   ├── gbm_ddr_signature_overlap.csv
│   ├── gbm_depmap_volcano.png
│   ├── gbm_depmap_waterfall.png
│   ├── gbm_depmap_target_boxplots.png
│   └── gbm_depmap_ddr_signature.png
│
├── docs/                            <- Detailed sub-analysis documentation
│   ├── Q1_Q3_pathway_signature.md   # Main signaling ODEs & Q1 signature documentation
│   ├── Q1_MGMT_methylation.md       # Cross-dataset MGMT promoter methylation classifier
│   ├── Q1_Q3_radiotherapy_ODE.md    # Radiotherapy response & LUMIERE volume ODE docs
│   ├── Q4_target_discovery.md       # DepMap CRISPR selective target discovery docs
│   └── Q5_treatment_selection.md    # Treatment selection & prognostic/predictive synthesis
│
└── presentation/                    <- Slides & final presentation deliverables
```

---

## Note on External Prerequisites & Assignment 2 Dependencies

- **Data Download**: Raw cBioPortal matrices (`gbm_tcga_pan_can_atlas_2018` and `gbm_cptac_2021`) are omitted due to file size. Place them in `data/` before running preprocessing.
- **`Q1_extract_depmap_signature.R`**: Relies on `GBM_cellline_coefs.csv` (the Assignment 2 GDSC2/DepMap 10-gene temozolomide sensitivity LASSO coefficients). If re-running from scratch, ensure `GBM_cellline_coefs.csv` is present in `analysis/`.

---

## How to Run the Pipeline

Each pathway follows a 3-step execution workflow: **Preprocessing $\to$ ODE Simulation $\to$ Survival Analysis**.

1. **p53 Pathway Pipeline**:
   ```bash
   Rscript preprocessing/GBM_p53_preproc.R
   # Execute models/GBM_p53_model.ipynb in Jupyter
   Rscript analysis/GBM_Analysis.R
   
   Rscript preprocessing/CPTAC_p53_preproc.R
   # Execute models/CPTAC_p53_model.ipynb in Jupyter
   Rscript analysis/CPTAC_Q1_validation.R
   ```

2. **Combined Q1 Multi-Pathway Signature**:
   ```bash
   Rscript analysis/Q1_extract_IDH_status.R
   Rscript analysis/Q1_extract_depmap_signature.R
   Rscript analysis/Q1_build_signature_TCGA.R
   Rscript analysis/Q1_validate_signature_CPTAC.R
   Rscript analysis/Q2_Q1_correlation.R
   ```

3. **Q4 DepMap CRISPR Target Discovery**:
   Execute `models/Q4_GBM_DepMap_LINCS_Target_Discovery.ipynb` directly in Jupyter.
