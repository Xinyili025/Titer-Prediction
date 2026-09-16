set.seed(123)

# Extract datasets with full dilution series from the sets_od list
full_sets <- list(cobovax_ancestral=sets_od$cobovax_ancestral$n12[[1]],
                  cobovax_omicron=sets_od$cobovax_omicron$n12[[1]],
                  cobovax_nctd=sets_od$cobovax_nctd$n12[[1]],
                  pattinson=sets_od$pattinson$n8[[1]])

# Gold standard: fit individual curves by 5PL model on the full dilution series
ec50_idv_models <- lapply(names(full_sets), function(ds) {
  n_samples <- length(unique(full_sets[[ds]]$Sample))
  bounds <- bounds_list[[ds]]
  
  lowerl <- c(rep(bounds$b[1], n_samples), rep(bounds$c[1], n_samples), 
              rep(bounds$d[1], n_samples), rep(bounds$e[1], n_samples), 
              rep(bounds$f[1], n_samples))
  
  upperl <- c(rep(bounds$b[2], n_samples), rep(bounds$c[2], n_samples), 
              rep(bounds$d[2], n_samples), rep(bounds$e[2], n_samples), 
              rep(bounds$f[2], n_samples))
  
  drm(OD ~ logDilution, curveid = Sample, fct = L.5(), data = full_sets[[ds]],
      lowerl = lowerl, upperl = upperl)
})

names(ec50_idv_models) <- names_datasets

# Extract EC50
ec50_idv <- lapply(ec50_idv_models, function(mod) {
  ec50 <- ED(mod, 50, display = FALSE)[ , "Estimate"]
  
  # Clean names: remove e: prefix and :50 suffix
  names(ec50) <- gsub("^e:|:50$", "", names(ec50))
  
  return(ec50)
})