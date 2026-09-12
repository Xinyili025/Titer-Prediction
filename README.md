# Optimal Dilution Points for Reliable SARS-CoV-2 Antibody Titer Prediction
In this project, we evaluated whether titers could be predicted by fitting a shared-parameter logistic model with sample-specific inflection points to a smaller set of selected dilution points.

## Datasets
Four datasets were analyzed in this study. Three of these were derived from the Cobovax study, an open-label, randomized trial of COVID-19 booster vaccinations in late 2021 and early 2022 in Hong Kong [1]. The fourth dataset was the Pattinson dataset, which comprised a community cohort study by the Marshfield Clinic Research Institute after the emergence of SARS-CoV-2 [2].

| Dataset | File |
| :--- | :--- |
| Cobovax - Ancestral virus spike RBD | `Cobovax_Ancestral.csv` |
| Cobovax - Omicron BA.2 full spike | `Cobovax_Omicron.csv` |
| Cobovax - N-CTD | `Cobovax_N-CTD.csv` |
| Pattinson | `pattinson_et_al_2022.csv` |

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
```
## References
[1] Leung NHL, Cheng SMS, Cohen CA, Martín-Sánchez M, Au NYM, Luk LLH, et al. Comparative antibody and cell-mediated immune responses, reactogenicity, and efficacy of homologous and heterologous boosting with CoronaVac and BNT162b2 (Cobovax): an open-label, randomised trial. Lancet Microbe. 2023;4: e670–e682.

[2] Pattinson D, Jester P, Guan L, Yamayoshi S, Chiba S, Presler R, et al. A Novel Method to Reduce ELISA Serial Dilution Assay Workload Applied to SARS-CoV-2 and Seasonal HCoVs. Viruses. 2022;14. doi:10.3390/v14030562

