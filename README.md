# Optimal Dilution Points for Reliable SARS-CoV-2 Antibody Titer Prediction
In this project, we evaluated whether titers could be predicted by fitting a shared-parameter logistic model with sample-specific inflection points to a smaller set of selected dilution points.

This repository contains the datasets and code to replicate all results and figures in the paper.

## Datasets
Four datasets were analyzed in this study. Three of these were derived from the Cobovax study, an open-label, randomized trial of COVID-19 booster vaccinations in late 2021 and early 2022 in Hong Kong [1]. The fourth dataset was the Pattinson dataset, which comprised a community cohort study by the Marshfield Clinic Research Institute after the emergence of SARS-CoV-2 [2].

| Dataset | File |
| :--- | :--- |
| Cobovax - Ancestral virus spike RBD | `Cobovax_Ancestral virus spike RBD.csv` |
| Cobovax - Omicron BA.2 full spike | `Cobovax_Omicron BA.2 full spike.csv` |
| Cobovax - N-CTD | `Cobovax_N-CTD.csv` |
| Pattinson | `Pattinson.csv` |

## Methods (summary)
- Cross-validation was used to select the optimal model for each dataset.
- All candidate reduced-point combinations were evaluated. The combination that yielded the lowest mean squared error (MSE) between predicted titers and the gold-standard titers was selected as optimal.
- Predictive accuracy for each optimal combination was further evaluated by calculating the ratio of predicted titers to gold-standard titers.
- Assay-capped scenario (for the Pattinson dataset): A limited total of 120 assays was allocated across different strategies (e.g., 60 samples × 2 points, 40 × 3, 30 × 4, 24 × 5). Predicted geometric mean titer (GMT) and its 95% confidence interval were calculated for each strategy for each iteration and compared to the gold-standard GMT.
- Real-world scenario (for the Pattinson dataset): A pilot subset of 50 samples was randomly selected as the training set, and the remaining samples were used as the test set. Predicted GMT and its 95% confidence interval were calculated for the 3-point strategy for each iteration and compared to the gold-standard GMT.

## How to run
1. Unzip or clone this repository.
2. Open `run_all.R` and set `setwd()` to your working directory.
3. Run `run_all.R` to reproduce the analyses and generate figures.

## Software

The analysis was conducted using R 4.4.1 and the following packages:

- tidyverse 2.0.0
- drc 3.0-1
- car 3.1-3
- scales 1.4.0
- writexl 1.5.1
- foreach 1.5.2
- doSNOW 1.0.20

## References
[1] Leung NHL, Cheng SMS, Cohen CA, Martín-Sánchez M, Au NYM, Luk LLH, et al. Comparative antibody and cell-mediated immune responses, reactogenicity, and efficacy of homologous and heterologous boosting with CoronaVac and BNT162b2 (Cobovax): an open-label, randomised trial. Lancet Microbe. 2023;4: e670–e682.

[2] Pattinson D, Jester P, Guan L, Yamayoshi S, Chiba S, Presler R, et al. A Novel Method to Reduce ELISA Serial Dilution Assay Workload Applied to SARS-CoV-2 and Seasonal HCoVs. Viruses. 2022;14. doi:10.3390/v14030562

## Contact
This repository was developed by Xinyi Li. For any enquiries, please contact Xinyi Li (u3012021@connect.hku.hk) or Prof. Ben Cowling (bcowling@hku.hk).

