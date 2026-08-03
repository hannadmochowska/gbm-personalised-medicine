# Question 5: Treatment Selection in Glioblastoma — Prognostic vs. Predictive Signatures

## Study Question

- **Q5 — How should molecular biomarkers and ODE-derived risk signatures be deployed for treatment selection in Glioblastoma? Are current signatures prognostic or predictive?**

---

## Key Synthesis & Clinical Takeaways

### 1. Prognostic vs. Predictive Distinction

- **Prognostic Biomarkers**: Indicate patient outcome (e.g., overall survival) regardless of specific therapy received.
  - *Example*: **Age** and **IDH mutation status** ($IDH1/2$) are strong, universal prognostic markers in Glioblastoma. Patients with $IDH$-mutant tumors have significantly longer survival regardless of the treatment regimen.
  - *Q1 Multi-Pathway Signature*: The 5-feature regularized signature (*SIAH1*, *PTEN*, *MAPK1*, p53 ODE state, age; $C\text{-index} = 0.589, p = 0.022$) is primarily **prognostic**, risk-stratifying patients into high- and low-risk survival tiers.
- **Predictive Biomarkers**: Specifically predict differential benefit or resistance from a particular treatment.
  - *Example*: **MGMT Promoter Methylation** is a predictive biomarker for benefit from alkylating chemotherapy (Temozolomide / TMZ). Methylation silences MGMT DNA repair, preventing repair of $O^6$-methylguanine lesions induced by TMZ.

---

### 2. Clinical Biomarker Deployment & Platform Limitations

1. **Cross-Platform Calibration Drift**:
   - The MGMT expression classifier developed in TCGA ($AUC = 0.682$ in CPTAC) demonstrated that while ranking discrimination was retained across cohorts, fixed probability thresholding failed ($100\%$ sensitivity, $0\%$ specificity at 0.50 cutoff).
   - *Clinical Implication*: Expression-based surrogate classifiers cannot replace standardized methylation assays (e.g., MGMT-STP27 or pyrosequencing) for patient selection in clinical trials without cohort-specific recalibration.

2. **Divergence of In Vitro Drug Screens and In Vivo Survival**:
   - The lack of correlation between Q2 in vitro DepMap TMZ drug-sensitivity scores and Q1 in vivo patient survival risk ($r \approx 0.04, p > 0.50$) highlights that in vitro viability assays lack key in vivo barriers:
     - The **blood-brain barrier (BBB)** restricting systemic drug penetration.
     - Tumor microenvironmental interactions and hypoxia.
     - Spatial tumor heterogeneity and residual stem cell niches.

---

### 3. Treatment-Aware Modeling Framework for Future Work

To transition from purely prognostic risk scoring to actionable predictive treatment selection models, future ODE and machine learning pipelines must incorporate:
- **Treatment-Aware ODE Formulations**: Explicit modeling of radiation dose fractioning ($d$, $D$) and temozolomide dosing schedules ($D_{\text{TMZ}}(t)$) rather than unadjusted natural growth rates.
- **Matched Multi-Omics + Radiomics**: Combining pre- and post-treatment 3D MRI radiomics (e.g., LUMIERE longitudinal volume segmentations) with molecular profiling to capture treatment response kinetics over time.
