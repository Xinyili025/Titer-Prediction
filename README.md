# Optimal Dilution Points for Reliable SARS-CoV-2 Antibody Titer Prediction
In this project, we evaluated whether titers could be predicted by fitting a shared-parameter logistic model with sample-specific inflection points to a smaller set of selected dilution points.

## Datasets
| Dataset | File |
| :--- | :--- |
| Cobovax - Ancestral virus spike RBD | `Cobovax_Ancestral.csv` |
| Cobovax - Omicron BA.2 full spike | `Cobovax_Omicron.csv` |
| Cobovax - N-CTD | `Cobovax_N-CTD.csv` |
| Pattinson et al. (2022) | `pattinson_et_al_2022.csv` |

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