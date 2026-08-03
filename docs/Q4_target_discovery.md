# Question 4: Glioblastoma Target Discovery via DepMap CRISPR Screening & LINCS Integration

## Study Question

- **Q4 — What are the most selectively essential targets in Glioblastoma, and how do they map to DDR pathways and LINCS perturbation signatures?**

This analysis identifies lineage-specific functional dependencies in Glioblastoma (GBM) using genome-wide CRISPR-Cas9 knockout screens from the Broad Institute DepMap (23Q4) dataset, cross-references top hits against the p53-DDR pathway signature (connecting Q1/Q3 to Q4), and evaluates landmark gene coverage in LINCS L1000.

---

## Analytical Workflow

1. **Cohort Curation**:
   - **GBM / Glioma Lines**: 68 cell lines with matching CRISPR Chronos dependency data and Oncotree metadata.
   - **Pan-Cancer Background**: 1,032 non-GBM cell lines across 30+ tumor types.

2. **Differential Dependency Scoring ($\Delta$)**:
   $$\Delta_{\text{GBM}}(g) = \text{Mean}(\text{Chronos Score}_{g, \text{GBM}}) - \text{Mean}(\text{Chronos Score}_{g, \text{Other}})$$
   A negative $\Delta$ score indicates selective lethality in Glioblastoma cell lines relative to other cancer lineages.

3. **Statistical Evaluation**:
   - Welch's two-sample t-tests ($P$-value), Benjamini-Hochberg False Discovery Rate ($\text{FDR}$), and Cohen's $d$ effect sizes computed across all 18,443 genes.

---

## Key Results

### Top Selective Target Candidates in Glioblastoma

| Rank | Gene | $\Delta$ (Diff Score) | $P$-Value | FDR | Cohen's $d$ | Biological Function |
| :---: | :--- | :---: | :---: | :---: | :---: | :--- |
| **#1** | **`CHMP4B`** | **-0.578** | $5.76 \times 10^{-14}$ | $3.23 \times 10^{-11}$ | **-0.958** | Core ESCRT-III complex subunit; membrane scission and cytokinesis. |
| **#2** | **`RPP25L`** | **-0.538** | $3.78 \times 10^{-13}$ | $1.78 \times 10^{-10}$ | **-1.286** | Ribonuclease P/MRP subunit involved in RNA processing. |
| **#3** | **`FERMT2`** | **-0.510** | $8.35 \times 10^{-12}$ | $3.00 \times 10^{-9}$ | **-0.985** | Kindlin-2; integrin-binding focal adhesion protein driving migration/invasion. |
| **#4** | **`GPX4`** | **-0.460** | $4.88 \times 10^{-6}$ | $3.39 \times 10^{-4}$ | **-0.833** | Glutathione peroxidase 4, master regulator of ferroptosis susceptibility. |
| **#5** | **`KIF18A`** | **-0.421** | $2.40 \times 10^{-7}$ | $2.77 \times 10^{-5}$ | **-0.638** | Kinesin motor protein required for mitotic spindle assembly. |
| **#7** | **`ITGAV`** | **-0.384** | $4.29 \times 10^{-13}$ | $1.93 \times 10^{-10}$ | **-0.957** | Integrin $\alpha_v$ subunit, key driver of cell adhesion and angiogenesis. |

### Evaluation of Lineage Master Factors & Drivers

- **`SOX2` & `OLIG2`**: Uncorrected nominal p-values are $<0.05$ ($p = 0.034$ and $p = 0.048$), but differential effect sizes are modest ($\Delta = -0.114$ and $-0.028$). Neither survives genome-wide FDR correction ($\text{FDR} = 0.150$ and $0.188$). This reflects a limitation of 2D high-serum monolayer culture, where stem cell self-renewal programs maintained in 3D/in vivo niches are dispensable for 2D cell proliferation.
- **`CDK4`**: Exhibits $p = 0.562, \text{FDR} = 0.770, \Delta = +0.033$. It acts as a **pan-cancer dependency** across 2D cell lines (Mean score $\approx -0.74$ to $-0.77$) rather than a lineage-selective hit.

### p53-DDR Signature & LINCS Alignment

- **`BRCA1`** (Rank #322, $p = 0.0109, \text{FDR} = 0.074$): Shows significant selective essentiality in GBM.
- **LINCS L1000 Coverage**: 3 of top 30 candidates (`JUN`, `ITGB5`, `PTK2`) and 9 of 19 DDR signature genes (`BRCA1`, `POLB`, `RPA1`, `CHEK2`, `TP53`, `PARP1`, `CHEK1`, `TP53BP1`, `CDKN1A`) are present in the LINCS L1000 landmark gene set.

---

## Deliverables & Artifacts

- **Executable Notebook**: `models/Q4_GBM_DepMap_LINCS_Target_Discovery.ipynb`
- **Output Tables**: `figures/gbm_top30_targets.csv`, `figures/gbm_selective_dependencies.csv`, `figures/gbm_ddr_signature_overlap.csv`
- **Publication Figures**: `figures/gbm_depmap_volcano.png`, `figures/gbm_depmap_waterfall.png`, `figures/gbm_depmap_target_boxplots.png`, `figures/gbm_depmap_ddr_signature.png`
