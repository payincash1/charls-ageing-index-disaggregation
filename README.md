# charls-ageing-index-disaggregation

Reproducibility archive for:

> **What is lost when multidimensional ageing is reduced to one number? Evidence from the CHARLS cohort**
> (*The Lancet Healthy Longevity*, submitted)

All analysis code is under [`charls_mpai_reproducibility/`](charls_mpai_reproducibility/).

## What this code does

End-to-end pipeline for a national cohort of middle-aged and older Chinese adults
(CHARLS 2011–2018):

- Multiple imputation (m=20) of eight routine ageing indicators
- Variational autoencoder latent extraction, Procrustes-aligned across imputations
- Cox survival models comparing scalar composites against the four-dimensional representation
- Discrimination, calibration, decision-curve, NRI/IDI, and a battery of sensitivity analyses

The study quantifies how much prognostic signal is retained when multidimensional ageing
measurement is compressed into a single number, across alternative scalar constructions of the
same eight indicators.

## Quick start

See [`charls_mpai_reproducibility/README.md`](charls_mpai_reproducibility/README.md) for the
step-by-step pipeline, runtime estimates, and required inputs.

## Data

CHARLS data are publicly available upon application (https://charls.pku.edu.cn) and are **not**
included in this archive, in accordance with the CHARLS data use agreement.

## Citation

Analysis code for this study is archived at:

**https://doi.org/10.5281/zenodo.22745069** (concept DOI — resolves to the latest version)

Please also cite the accompanying manuscript once published.

## License

MIT — see [`charls_mpai_reproducibility/LICENSE`](charls_mpai_reproducibility/LICENSE).
