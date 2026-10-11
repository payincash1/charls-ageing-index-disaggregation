# Multidimensional ageing index disaggregation — CHARLS (reproducibility archive)

Code accompanying: *What is lost when multidimensional ageing is reduced to one number? Evidence from the CHARLS cohort* (The Lancet Healthy Longevity, submitted).

## Overview
End-to-end pipeline: CHARLS 2011–2018 data preparation → multiple imputation (m=20) →
variational autoencoder latent extraction (20 imputed datasets, Procrustes-aligned) →
survival analysis (Cox, cluster-robust), discrimination/calibration/decision-curve metrics,
scalar-construction comparators, single-marker benchmarks, sensitivity analyses.

## Pipeline (run in order)

| Step | Script | Runtime |
|---|---|---|
| 1 Data + imputation | `scripts/01b_impute_m20.R` | 10–30 min |
| 2 VAE training + alignment | `scripts/02b_vae_m20.py` | ~2 h |
| 3 Main metrics | `scripts/03b_metrics.R` | ~25 min |
| 4 Calibration/AUC/Brier | `scripts/04_calibration.R` | ~5 min |
| 5 Single-marker benchmarks | `scripts/08_single_marker.R` | ~2 min |
| 6 Stability (10 seeds) | `scripts/vae_stability_4dim.py` | ~40 min |
| 7 Latent-dimensionality sensitivity | `scripts/15_dim_sensitivity.py` then `scripts/15b_dim_sensitivity_R.R` | ~25 min |
| 8 Held-out test evaluation | `scripts/18_test_set.py` then `scripts/18_test_set_R.R` | ~2 min |
| 9 CRP acute-phase sensitivity | `scripts/20_crp_acute_phase_sensitivity.R` | ~2 min |
| Aux | `scripts/10_table1_fill.R`, `scripts/11_tics_check.R`, `scripts/14_second_round_analyses.R` | minutes |

Steps 1→2→3→4 sequential; 5–9 independent. `base_path`/`base` at the top of each
script must point to the local CHARLS root folder (contains `chalrs2011` … `chalrs2018`).
Step 9 additionally requires the objects saved by step 3 (`output/metrics_final.RData`)
and the standardised indicator matrices written by step 1.

## Data
CHARLS data are publicly available upon application (https://charls.pku.edu.cn) and are
**not** included in this archive. Expected inputs: per-wave folders with the standard
CHARLS .dta files (Demographic_Background, Health_Status_and_Functioning, Blood,
weight, Exit_Interview, Sample_Infor).

## Dependencies
R ≥ 4.2: dplyr, haven, mice, survival, survminer, boot, survIDINRI, rms, riskRegression,
dcurves, ggplot2, lavaan (optional), readr.
Python 3.9+: torch, pandas, numpy, scipy (see `requirements.txt`).

## Key fixed settings
- ID convention: CHARLS 12-character standard IDs (conversion in `01b`).
- Imputation: mice PMM, m=20, maxit=20, seed 2024; event indicator + Nelson–Aalen included.
- VAE: encoder 8→64 (ReLU) → 4-dim Gaussian latent; decoder 4→64→8; KL weight 0.5;
  80/10/10 split; early stopping patience 20; Adam 1e-3; batch 256; seed 2024.
- Cox: community-clustered robust SEs; bootstrap B=1000 (cluster bootstrap).

## License
MIT (see LICENSE). Data availability per CHARLS terms.

## Citation
Archived release (concept DOI, resolves to the latest version):
https://doi.org/10.5281/zenodo.22745069
