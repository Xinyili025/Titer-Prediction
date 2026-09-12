# Find optimal dilution combination for each number of points (2-5)
find_optimal_model <- function(ed_shared_list,ec50_idv_full,dilutions_vector) {
  
  # Gold-standard titers from individual model on full dilution series
  gs <- ec50_idv_full
  
  results <- data.frame(n_points = integer(), optimal_comb = character(), MSE = numeric(), 
                        n_success = integer(), positions = character(), dilutions = character())
  
  # Calculate MSE for each dilution combination
  for(n in reduced_points) {
    
    # Skip combination with NAs
    mse <- sapply(ed_shared_list[[paste0("n", n)]], function(x) {
      if(all(is.na(x))) return(NA)
      mean((x - gs)^2)
    })
    
    # Remove failed combinations (all NA)
    mse_clean <- mse[!is.na(mse)]
    n_success <- length(mse_clean) 
    
    # Record results
    # No successful combinations
    if(n_success == 0) {
      results <- rbind(results, data.frame(
        n_points = n,
        optimal_comb = NA,
        MSE = NA,
        n_success = 0,
        positions = NA,
        dilutions = NA
      ))
    } else {
      optimal <- names(mse_clean)[which.min(mse_clean)]
      optimal_mse <- min(mse_clean)
      
      # Extract positions as numeric string and corresponding dilutions
      positions <- gsub("c\\(|\\)", "", optimal)
      idx <- as.numeric(strsplit(positions, ",")[[1]])
      dilutions <- paste(dilutions_vector[idx], collapse = ", ")
      
      results <- rbind(results, data.frame(
        n_points = n,
        optimal_comb = optimal,
        MSE = round(optimal_mse,5),
        n_success = n_success,
        positions = positions,
        dilutions = dilutions
      ))
    }
    
  }
  
  return(results)
}

# Optimal dilution combination selection
optimal_models <- lapply(cobovax_names, function(ds) {
  find_optimal_model(ec50_shared[[ds]], ec50_idv[[ds]], dilution_cobovax)
})
names(optimal_models) <- cobovax_names

optimal_models$pattinson <- find_optimal_model(ec50_shared$pattinson, ec50_idv$pattinson, dilution_pattinson)

# Extract optimal logtiters from optimal models
# Gold standard titers was extracted from individual model on full dilution series
extract_optimal_logtiter <- function(logtiter_list, optimal_df, ec50_idv_values) {
  data.frame(
    logtiter_n2 = logtiter_list$n2[[optimal_df$optimal_comb[optimal_df$n_points == 2]]],
    logtiter_n3 = logtiter_list$n3[[optimal_df$optimal_comb[optimal_df$n_points == 3]]],
    logtiter_n4 = logtiter_list$n4[[optimal_df$optimal_comb[optimal_df$n_points == 4]]],
    logtiter_n5 = logtiter_list$n5[[optimal_df$optimal_comb[optimal_df$n_points == 5]]],
    logtiter_gs = ec50_idv_values
  )
}

# Apply to all datasets using lapply
optimal_logtiters <- lapply(names(optimal_models), function(ds) {
  extract_optimal_logtiter(ec50_shared[[ds]], optimal_models[[ds]], ec50_idv[[ds]])
})
names(optimal_logtiters) <- names(optimal_models)

# Calculate ratio of predicted titer to gold-standard titer
optimal_logtiters[cobovax_names] <- lapply(optimal_logtiters[cobovax_names], function(x) {
  x[, paste0("ratio_n", 2:5)] <- x[, paste0("logtiter_n", 2:5)] - x$logtiter_gs
  return(x)
})
optimal_logtiters$pattinson[paste0("ratio_n", 2:5)] <- optimal_logtiters$pattinson[, paste0("logtiter_n", 2:5)] - optimal_logtiters$pattinson$logtiter_gs 
