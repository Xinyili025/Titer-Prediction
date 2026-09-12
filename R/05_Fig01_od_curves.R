# Main Figure 1  = shared 5PL 
# Appendix Fig 1 = shared 4PL
# Appendix Fig 2 = shared 3PL

# Split each dataset by Sample
od_by_sample <- lapply(od, function(df) {
  split(df, df$Sample)
})

# Create datasets for prediction
prediction_data <- lapply(names_datasets, function(ds) {
  log_dil <- if(ds == "pattinson") logdilution_pattinson else logdilution_cobovax
  smooth_dil <- seq(log_dil[1], log_dil[length(log_dil)], length.out = 200)
  
  samples <- unique(od[[ds]]$Sample)
  expand.grid(logDilution = smooth_dil, Sample = samples)
})
names(prediction_data) <- names_datasets

# Prediction
create_predictions <- function(model,data) {
  results <- list()
  
  pred_values <- predict(model, newdata = data)
  
  results <- data.frame(
    Sample = data$Sample,
    logDilution = data$logDilution,
    OD_fit = pred_values
  )
  
  return(results)
}

# Predicted values by individual 5PL model
predicted_values_ind <- lapply(names_datasets, function(ds) {
  cat("\n=== Processing:", ds, "===\n")
  create_predictions(ec50_idv_models[[ds]], prediction_data[[ds]])
})

names(predicted_values_ind) <- names_datasets

# Fit shared full models
fit_shared_full <- function(data, model_pl, model_parameter, bounds) {
  n_samples <- length(unique(data$Sample))
  n_params <- length(model_parameter)
  
  if (n_params == 3) {
    lowerl <- c(bounds$b[1], bounds$d[1], rep(bounds$e[1], n_samples))
    upperl <- c(bounds$b[2], bounds$d[2], rep(bounds$e[2], n_samples))
  } else if (n_params == 4) {
    lowerl <- c(bounds$b[1], bounds$c[1], bounds$d[1], rep(bounds$e[1], n_samples))
    upperl <- c(bounds$b[2], bounds$c[2], bounds$d[2], rep(bounds$e[2], n_samples))
  } else if (n_params == 5) {
    lowerl <- c(bounds$b[1], bounds$c[1], bounds$d[1], rep(bounds$e[1], n_samples), bounds$f[1])
    upperl <- c(bounds$b[2], bounds$c[2], bounds$d[2], rep(bounds$e[2], n_samples), bounds$f[2])
  }
  
  drm(OD ~ logDilution, curveid = Sample,
      pmodels = model_parameter, fct = model_pl, data = data,
      lowerl = lowerl, upperl = upperl)
}

# 3PL
# shared_full <- list(
#   cobovax_ancestral = fit_shared_full(full_sets$cobovax_ancestral, L.3(), model_parameter_3pl, bounds_list$cobovax_ancestral),
#   cobovax_omicron   = fit_shared_full(full_sets$cobovax_omicron, L.3(), model_parameter_3pl, bounds_list$cobovax_omicron),
#   cobovax_nctd      = fit_shared_full(full_sets$cobovax_nctd, L.3(), model_parameter_3pl, bounds_list$cobovax_nctd),
#   pattinson         = fit_shared_full(full_sets$pattinson, L.3(), model_parameter_3pl, bounds_list$pattinson)
# )
# names(shared_full) <- names_datasets

# 4PL
# shared_full <- list(
#   cobovax_ancestral = fit_shared_full(full_sets$cobovax_ancestral, L.4(), model_parameter_4pl, bounds_list$cobovax_ancestral),
#   cobovax_omicron   = fit_shared_full(full_sets$cobovax_omicron, L.4(), model_parameter_4pl, bounds_list$cobovax_omicron),
#   cobovax_nctd      = fit_shared_full(full_sets$cobovax_nctd, L.4(), model_parameter_4pl, bounds_list$cobovax_nctd),
#   pattinson         = fit_shared_full(full_sets$pattinson, L.4(), model_parameter_4pl, bounds_list$pattinson)
# )
# names(shared_full) <- names_datasets

# 5PL
shared_full <- list(
  cobovax_ancestral = fit_shared_full(full_sets$cobovax_ancestral, L.5(), model_parameter_5pl, bounds_list$cobovax_ancestral),
  cobovax_omicron   = fit_shared_full(full_sets$cobovax_omicron, L.5(), model_parameter_5pl, bounds_list$cobovax_omicron),
  cobovax_nctd      = fit_shared_full(full_sets$cobovax_nctd, L.5(), model_parameter_5pl, bounds_list$cobovax_nctd),
  pattinson         = fit_shared_full(full_sets$pattinson, L.5(), model_parameter_5pl, bounds_list$pattinson)
)
names(shared_full) <- names_datasets

# Predicted values by shared model
predicted_values_shared <- lapply(names_datasets, function(ds) {
  cat("\n=== Processing:", ds, "===\n")
  create_predictions(shared_full[[ds]], prediction_data[[ds]])
})

names(predicted_values_shared) <- names_datasets
