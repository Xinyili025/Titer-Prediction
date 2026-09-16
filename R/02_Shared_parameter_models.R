set.seed(123)

# Fit shared-parameter models and extract EC50 values
shared_ec50_model <- function(data,model_pl,model_parameter,bounds) {
  
  n_samples <- length(unique(data$Sample))
  n_params <- length(model_parameter)
  
  if (n_params == 3) {
    # 3PL
    lowerl <- c(bounds$b[1], bounds$d[1], rep(bounds$e[1], n_samples))
    upperl <- c(bounds$b[2], bounds$d[2], rep(bounds$e[2], n_samples))
  } else if (n_params == 4) {
    # 4PL
    lowerl <- c(bounds$b[1], bounds$c[1], bounds$d[1], rep(bounds$e[1], n_samples))
    upperl <- c(bounds$b[2], bounds$c[2], bounds$d[2], rep(bounds$e[2], n_samples))
  } else if (n_params == 5) {
    # 5PL
    lowerl <- c(bounds$b[1], bounds$c[1], bounds$d[1], rep(bounds$e[1], n_samples), bounds$f[1])
    upperl <- c(bounds$b[2], bounds$c[2], bounds$d[2], rep(bounds$e[2], n_samples), bounds$f[2])
  }
  
  # Fit shared-parameter model
  model <- drm(OD ~ logDilution, 
               curveid = Sample, 
               pmodels = model_parameter, 
               fct = model_pl, 
               data = data,
               lowerl = lowerl,
               upperl = upperl)
  
  # Calculate EC50
  ec50 <- ED(model, 50, display = FALSE)[ , "Estimate"]
  
  # Clean names
  names(ec50) <- gsub("^e:|:50$", "", names(ec50))
  
  return(ec50)
  
}

# Loop over all combinations
batch_fit_shared_ec50 <- function(comb_list, antigen,model_pl,model_parameter,bounds) {
  result <- list()
  
  for(i in seq_along(comb_list)) {
    cat(antigen, "| Combination", i, "/", length(comb_list), "\n")
    
    # Return NA when model fails to converge
    result[[i]] <- tryCatch({
      
      shared_ec50_model(comb_list[[i]],model_pl,model_parameter,bounds)
      
    }, error = function(e) {
      cat("Fail", e$message, "\n")
      
      n_samples <- length(unique(comb_list[[i]]$Sample))
      
      return(rep(NA, n_samples))
    })
  }
  
  names(result) <- names(comb_list)
  return(result)
}

# Model syntax for 3PL model
model_pl_3pl <- L.3()
model_parameter_3pl <- list(~1, ~1, ~Sample-1)

# Model syntax for 4PL model
model_pl_4pl <- L.4()
model_parameter_4pl <- list(~1, ~1,  ~1,~Sample-1)

# Model syntax for 5PL model
model_pl_5pl <- L.5()
model_parameter_5pl <- list(~1, ~1, ~1, ~Sample-1, ~1)

ec50_shared <- list(
  cobovax_ancestral = lapply(sets_od$cobovax_ancestral, batch_fit_shared_ec50, 
                             antigen = "cobovax_ancestral", 
                             model_pl = model_pl_5pl, model_parameter = model_parameter_5pl,
                             bounds = bounds_list$cobovax_ancestral),
  
  cobovax_omicron = lapply(sets_od$cobovax_omicron, batch_fit_shared_ec50, 
                           antigen = "cobovax_omicron", 
                           model_pl = model_pl_5pl, model_parameter = model_parameter_5pl,
                           bounds = bounds_list$cobovax_omicron),
  
  cobovax_nctd = lapply(sets_od$cobovax_nctd, batch_fit_shared_ec50, 
                        antigen = "cobovax_nctd", 
                        model_pl = model_pl_3pl, model_parameter = model_parameter_3pl,
                        bounds = bounds_list$cobovax_nctd),
  
  pattinson = lapply(sets_od$pattinson, batch_fit_shared_ec50, 
                     antigen = "pattinson", 
                     model_pl = model_pl_3pl, model_parameter = model_parameter_3pl,
                     bounds = bounds_list$pattinson)
)

# Same sample order as the gold-standard EC50 values
for (ds in names(ec50_shared)) {
  ec50_shared[[ds]] <- lapply(ec50_shared[[ds]], function(n_list) {
    lapply(n_list, function(x) {
      x[names(ec50_idv[[ds]])]
    })
  })
}
