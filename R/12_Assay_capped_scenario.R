# Fixed 120 assays on the Pattinson dataset
# Trade-off: More samples × fewer points  vs  fewer samples × more points
set.seed(123)

# Pattinson dataset
samples_pattinson <- unique(pattinson$Sample)

gs_it <- ec50_idv$pattinson
gs_logGMT <- mean(gs_it,na.rm = TRUE)
gs_GMT <- original_dilutions_pattinson(mean(gs_it,na.rm = TRUE))

# How many random subsets to draw per strategy
n_resample <- 100

# n2: 60×2, n3: 40×3, n4: 30×4, n5: 24×5
no_points <- 2:5
strategy_n_samples <- c(60, 40, 30, 24)
strategies <- paste0("n",2:5)
iter_names <- paste0("iter_", 1:n_resample)

# GMT and its SE
all_logGMT <- list(
  n2 = data.frame(logGMT = rep(NA, n_resample), log_se = rep(NA, n_resample)),
  n3 = data.frame(logGMT = rep(NA, n_resample), log_se = rep(NA, n_resample)),
  n4 = data.frame(logGMT = rep(NA, n_resample), log_se = rep(NA, n_resample)),
  n5 = data.frame(logGMT = rep(NA, n_resample), log_se = rep(NA, n_resample))
)

# For all dilution strategies
for(i in seq_along(strategies)) {
  s <- strategies[i]
  n_samples <- strategy_n_samples[i]
  
  cat("Strategy:", s, "(", i, "/", length(strategies), ") -", n_samples, "samples\n")
  
  dilutions_str <- optimal_models$pattinson$dilutions[[i]]
  selected_dilutions <- as.numeric(strsplit(dilutions_str, ", ")[[1]])
  
  for(r in 1:n_resample) {
    
    if(r %% 10 == 0) {
      cat("Iteration:", r, "/", n_resample, "\n")
    }
    
    selected_samples <- sample(samples_pattinson, n_samples, replace = FALSE)
    
    subset_data <- od$pattinson[
      od$pattinson$Sample %in% selected_samples & 
        od$pattinson$Dilution %in% selected_dilutions, ]
    
    result <- tryCatch({
      log_titers <- shared_ec50_model(data = subset_data,
                                      model_pl = model_pl_3pl,
                                      model_parameter = model_parameter_3pl,
                                      bounds = bounds_list$pattinson)
      list(
        logGMT = mean(log_titers, na.rm = TRUE),
        log_se = se(log_titers)
      )
      
    }, error = function(e) list(logGMT = NA, log_se = NA))
    
    all_logGMT[[s]]$logGMT[r] <- result$logGMT
    all_logGMT[[s]]$log_se[r] <- result$log_se
  }
}

# logGMT and its 95% CI
logGMT_fixed_assays_individual <- list(
  logGMT = data.frame(lapply(all_logGMT, function(x) x$logGMT)),
  ci_lower = data.frame(lapply(all_logGMT, function(x) x$logGMT - 1.96 * x$log_se)),
  ci_upper = data.frame(lapply(all_logGMT, function(x) x$logGMT + 1.96 * x$log_se))
)

# Rename columns
names(logGMT_fixed_assays_individual$logGMT) <- names(all_logGMT)
names(logGMT_fixed_assays_individual$ci_lower) <- names(all_logGMT)
names(logGMT_fixed_assays_individual$ci_upper) <- names(all_logGMT)

# Mse
mse_fixed_assays <- round(sapply(logGMT_fixed_assays_individual$logGMT, mse, true = gs_logGMT),4)
