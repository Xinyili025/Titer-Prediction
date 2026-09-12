# Packages
library(tidyverse)
library(drc)
library(car)
library(scales)
library(writexl)
library(here)

# Set seed for reproducibility
set.seed(123)

# Import data
cobovax_ancestral <- read.csv(here("data", "raw", "Cobovax_Ancestral virus spike RBD.csv"), header = TRUE)
cobovax_omicron   <- read.csv(here("data", "raw", "Cobovax_Omicron BA.2 full spike.csv"), header = TRUE)
cobovax_nctd      <- read.csv(here("data", "raw", "Cobovax_N-CTD.csv"), header = TRUE)
pattinson         <- read.csv(here("data", "raw", "pattinson.csv"), header = TRUE)

# Names of datasets
cobovax_names <- c("cobovax_ancestral", "cobovax_omicron", "cobovax_nctd")
names_datasets <- c(cobovax_names, "pattinson")

# Reduced Dilution Points
reduced_points <- c(2:5)

# Colors for plotting
dp_colors <- c(n2="#C62828" ,n3="#1976D2",n4="#4CAF50",n5="#7A1FBD",full="black")
data_colors <- c(cobovax_ancestral="#E25C3C" ,cobovax_omicron="#388E3C",cobovax_nctd="#1976D2",pattinson="black")

# Functions
# Standard error
se <- function(x){
  sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))
}

# Mean squared error
mse <- function(x, true) {
  mean((x - true)^2, na.rm = TRUE)
}

# Log-transformation of dilutions
log_dilutions_cobovax <- function(x){
  log (x / 100, base = 2)
}

log_dilutions_pattinson <- function(x){
  log (x / 40, base = 4)
}

# Inverse transformation (back to original dilution)
original_dilutions_cobovax <- function(x){
  100 * 2 ^ x
}

original_dilutions_pattinson <- function(x){
  40 * 4 ^ x
}

# Combine four datasets
od <- list(
  cobovax_ancestral = cobovax_ancestral,
  cobovax_omicron = cobovax_omicron,
  cobovax_nctd = cobovax_nctd,
  pattinson = pattinson
)

# Add the log-transformed dilution column to each data frame
od[cobovax_names] <- lapply(od[cobovax_names], function(df) {
  df$logDilution <- log_dilutions_cobovax(df$Dilution)
  return(df)
})

od$pattinson["logDilution"] <- log_dilutions_pattinson(od$pattinson$Dilution)

# Dilutions and logDilutions
dilution_cobovax <- unique(cobovax_ancestral$Dilution)
dilution_pattinson <- unique(pattinson$Dilution)
logdilution_cobovax <- log_dilutions_cobovax(dilution_cobovax)
logdilution_pattinson <- log_dilutions_pattinson(dilution_pattinson)

# Function to generate all combinations of dilution points
create_dilution_datasets <- function(data,n_points, n_max){
  all_comb <- combn(1:n_max, n_points, simplify = FALSE)
  result <- lapply(all_comb, function(x) data[data$Dilution %in% data$Dilution[x], ])
  
  names(result) <- sapply(all_comb, function(x) paste0("c(", paste(x, collapse = ","), ")"))
  return(result)
}

# Datasets with reduced numbers of dilutions and the full dilution series
sets_od <- c(
  # Cobovax datasets
  lapply(od[cobovax_names], function(df) {
    setNames(
      lapply(c(reduced_points, length(dilution_cobovax)), 
             function(n) create_dilution_datasets(df, n, length(dilution_cobovax))), 
      paste0("n", c(reduced_points, length(dilution_cobovax)))
    )
  }),
  
  # Pattinson dataset
  list(pattinson = setNames(
    lapply(c(reduced_points, length(dilution_pattinson)), 
           function(n) create_dilution_datasets(od$pattinson, n, length(dilution_pattinson))),
    paste0("n", c(reduced_points, length(dilution_pattinson)))
  ))
)


