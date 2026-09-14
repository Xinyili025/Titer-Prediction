# Users should set the working directory
# Replace this path with the location on your computer
setwd("/Users/lixinyi/Desktop/Titer-Prediction-main")

# R scripts are sourced in order
# Setup
source("R/00_Setup.R")

# Gold standard & shared models
source("R/01_Gold_standard_models.R")
source("R/02_Shared_parameter_models.R")

# OD curves, figure 1, appendix figure 1, appendix figure 2, and appendix table 1
source("R/03_Od_curves.R")
source("R/04_Fig01_od_curves.R")
source("R/05_Apptable01.R")

# Cross-validation for optimal models & appendix table 2
source("R/06_Cross_validation.R")

# Optimal combinations & appendix figure 3
source("R/07_Optimal_combinations.R")
source("R/08_Optimal_combinations_curves.R")
source("R/09_Appfig03_optimal_combinations_curves.R")

# Ratio & Figure 2
source("R/10_Ratio.R")
source("R/11_Fig02_ratio.R")

# Assay-capped scenario & Figure 3
source("R/12_Assay_capped_scenario.R")
source("R/13_Fig03_assay_capped_scenario.R")

# Real-world scenario & Appendix Figure 4
source("R/14_Real_world_scenario.R")
source("R/15_Appfig04_real_world_scenario.R")

cat("\nAnalysis complete! Check the \"Figures\" folder.\n")
