# Reduced-dilution ELISA titer prediction

This repository contains the analysis pipeline for estimating ELISA titers (EC50 / GMT) from a reduced number of dilution points, using shared-parameter dose–response models.

Predicted titers are compared with a gold standard obtained from individual five-parameter logistic (5PL) curves fitted to the full dilution series.

## Datasets

- Cobovax ancestral virus spike RBD
- Cobovax Omicron BA.2 full spike
- Cobovax N-CTD
- Pattinson et al. (2022)

Raw CSV files are not stored in this repository. If you have permission to use the data, place them in `data/raw/` as:

- `Cobovax_Ancestral.csv`
- `Cobovax_Omicron.csv`
- `Cobovax_N-CTD.csv`
- `pattinson_et_al_2022.csv`

## Methods (summary)

- Dilutions are log-transformed before model fitting.
- Candidate models: 3PL, 4PL, and 5PL. Shared parameters (slope and asymptotes; asymmetry for 5PL) are estimated across samples; only the titer parameter varies by sample.
- Model family is selected by cross-validation against the 5PL gold standard.
- For 2–5 dilution points, all combinations are evaluated; the combination with the lowest MSE versus the gold standard is retained.
- Additional analyses include a fixed assay budget (120 wells) and a train/test “real-world” split for combination selection.

## Repository layout

```text
data/raw/             local raw CSVs (not tracked)
R/                    helper functions
analysis/             main analysis scripts
figures/scripts/      figure code
output/               model objects and figures
docs/                 extra notes