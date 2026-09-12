fit_models <- function(data, best_df, model_pl, model_parameter, bounds) {
  
  dilutions_list <- lapply(best_df$dilutions, function(d) {
    as.numeric(strsplit(d, ", ")[[1]])
  })
  
  results <- lapply(dilutions_list, function(comb) {
    data_filtered <- data[data$Dilution %in% comb, ]
    n_samp <- length(unique(data_filtered$Sample))
    n_params <- length(model_parameter)
    
    if (n_params == 3) {
      lowerl <- c(bounds$b[1], bounds$d[1], rep(bounds$e[1], n_samp))
      upperl <- c(bounds$b[2], bounds$d[2], rep(bounds$e[2], n_samp))
    } else if (n_params == 4) {
      lowerl <- c(bounds$b[1], bounds$c[1], bounds$d[1], rep(bounds$e[1], n_samp))
      upperl <- c(bounds$b[2], bounds$c[2], bounds$d[2], rep(bounds$e[2], n_samp))
    } else if (n_params == 5) {
      lowerl <- c(bounds$b[1], bounds$c[1], bounds$d[1], rep(bounds$e[1], n_samp), bounds$f[1])
      upperl <- c(bounds$b[2], bounds$c[2], bounds$d[2], rep(bounds$e[2], n_samp), bounds$f[2])
    }
    
    drm(OD ~ logDilution, curveid = Sample, 
        pmodels = model_parameter, fct = model_pl, data = data_filtered,
        lowerl = lowerl, upperl = upperl)
  })
  
  names(results) <- c("n2", "n3", "n4", "n5")
  return(results)
}

fitted_optimal_models <- list(
  cobovax_ancestral = fit_models(od$cobovax_ancestral, optimal_models$cobovax_ancestral, 
                                 L.5(), model_parameter_5pl, bounds_list$cobovax_ancestral),
  cobovax_omicron   = fit_models(od$cobovax_omicron,   optimal_models$cobovax_omicron,   
                                 L.5(), model_parameter_5pl, bounds_list$cobovax_omicron),
  cobovax_nctd      = fit_models(od$cobovax_nctd,      optimal_models$cobovax_nctd,      
                                 L.3(), model_parameter_3pl, bounds_list$cobovax_nctd),
  pattinson         = fit_models(od$pattinson,         optimal_models$pattinson,         
                                 L.3(), model_parameter_3pl, bounds_list$pattinson)
)

names(fitted_optimal_models) <- names(optimal_models)

create_predictions_reduced_models <- function(models,data) {
  results <- list()
  
  for(name in names(models)) {
    cat("Prediction model", name, "...\n")
    
    pred_values <- predict(models[[name]], newdata = data)
    
    results[[name]] <- data.frame(
      Sample = data$Sample,
      logDilution = data$logDilution,
      OD_fit = pred_values
    )
  }
  
  return(results)
}

# Predicted values for models
predicted_values <- lapply(names(fitted_optimal_models), function(ds) {
  cat("\n=== Processing:", ds, "===\n")
  create_predictions_reduced_models(fitted_optimal_models[[ds]], prediction_data[[ds]])
})
names(predicted_values) <- names(fitted_optimal_models)
