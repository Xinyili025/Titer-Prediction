# Fraction of samples held out; number of random splits
holdout_pct = 0.2
n_iter = 100

# Manually calculate EC50 for 5PL model
# In 5PL, parameter 'e' is not the EC50
manual_EC50_5PL <- function(b, e, f) {
  e + (1/b) * log(2^(1/f) - 1)
}

cv_all_models <- function(data, holdout_pct, n_iter, seed = 123) {
  
  set.seed(seed)
  samples <- unique(data$Sample)
  n_holdout <- max(1, round(length(samples) * holdout_pct))
  
  # Store results 
  results_combined <- data.frame()
  
  for (iter in 1:n_iter) {
    if(iter %% 10 == 0) cat("Iteration:", iter, "\n")
    
    holdout_s <- sample(samples, n_holdout)
    train_s <- setdiff(samples, holdout_s)
    
    # Training models
    train_3pl <- tryCatch({
      drm(OD ~ logDilution, curveid = Sample,
          pmodels = list(~1, ~1, ~Sample-1),
          fct = L.3(), data = data[data$Sample %in% train_s, ])
    }, error = function(e) {
      cat("train_3pl failed:", e$message, "\n")
      NULL})
    
    train_4pl <- tryCatch({
      drm(OD ~ logDilution, curveid = Sample,
          pmodels = list(~1, ~1, ~1, ~Sample-1),
          fct = L.4(), data = data[data$Sample %in% train_s, ])
    }, error = function(e) {
      cat("train_4pl failed:", e$message, "\n")
      NULL})
    
    train_5pl <- tryCatch({
      drm(OD ~ logDilution, curveid = Sample,
          pmodels = list(~1, ~1, ~1, ~Sample-1, ~1),
          fct = L.5(), data = data[data$Sample %in% train_s, ])
    }, error = function(e) {
      cat("train_5pl failed:", e$message, "\n")
      NULL})
    
    # Training set failed for any model, skip this iteration
    if (any(sapply(list(train_3pl, train_4pl, train_5pl), is.null))) next
    
    # Fixed parameters by each model
    tc_3pl <- coef(train_3pl)
    b3 <- tc_3pl["b:(Intercept)"]
    d3 <- tc_3pl["d:(Intercept)"]
    
    tc_4pl <- coef(train_4pl)
    b4 <- tc_4pl["b:(Intercept)"]
    c4 <- tc_4pl["c:(Intercept)"]
    d4 <- tc_4pl["d:(Intercept)"]
    
    tc_5pl <- coef(train_5pl)
    b5 <- tc_5pl["b:(Intercept)"]
    c5 <- tc_5pl["c:(Intercept)"]
    d5 <- tc_5pl["d:(Intercept)"]
    f5 <- tc_5pl["f:(Intercept)"]
    
    # Testing samples
    for (s in holdout_s) {
      s_data <- data[data$Sample == s, ]
      
      # Calculate predicted and gold-standard titers
      
      # 3PL predicted titer
      pred_3pl <- tryCatch({
        m <- drm(OD ~ logDilution, fct = L.3(fixed = c(b3, d3, NA)), data = s_data)
        coef(m)["e:(Intercept)"]
      },error = function(e) {
        cat("pred_3pl failed:", e$message, "\n")
        NA})
      
      # 4PL predicted titer
      pred_4pl <- tryCatch({
        m <- drm(OD ~ logDilution, fct = L.4(fixed = c(b4, c4, d4, NA)), data = s_data)
        coef(m)["e:(Intercept)"]
      }, error = function(e) {
        cat("pred_4pl failed:", e$message, "\n")
        NA})
      
      # 5PL predicted titer
      pred_5pl <- tryCatch({
        m <- drm(OD ~ logDilution, fct = L.5(fixed = c(b5, c5, d5, NA, f5)), data = s_data)
        manual_EC50_5PL(b=b5,e=coef(m)["e:(Intercept)"],f=f5)
      }, error = function(e) {
        cat("pred_5pl failed:", e$message, "\n")
        NA} )
      
      # 5PL gold-standard titer
      gs_5pl <- tryCatch({
        m <- drm(OD ~ logDilution, fct = L.5(), data = s_data)
        
        manual_EC50_5PL(b=coef(m)["b:(Intercept)"],
                        e=coef(m)["e:(Intercept)"],
                        f=coef(m)["f:(Intercept)"])
        
      }, error = function(e) {
        cat("gs_5pl failed:", e$message, "\n")
        NA} )
      
      # For those all successfully calculated
      if (!any(is.na(c(pred_3pl,pred_4pl,pred_5pl,gs_5pl)))) {
        results_combined <- rbind(results_combined, data.frame(
          iter = iter, 
          sample = s,
          
          # 3PL
          pred_3pl = pred_3pl,
          # 4PL
          pred_4pl = pred_4pl,
          # 5PL
          pred_5pl = pred_5pl,
          gs_5pl = gs_5pl
          
        ))
      }
    }
  }
  
  results_3pl <- results_combined[, c("iter", "sample", "gs_5pl", "pred_3pl")]
  names(results_3pl) <- c("iter", "sample", "gs_e", "pred_e")
  
  results_4pl <- results_combined[, c("iter", "sample", "gs_5pl", "pred_4pl")]
  names(results_4pl) <- c("iter", "sample", "gs_e", "pred_e")
  
  results_5pl <- results_combined[, c("iter", "sample", "gs_5pl", "pred_5pl")]
  names(results_5pl) <- c("iter", "sample", "gs_e", "pred_e")
  
  return(list(
    `3PL` = results_3pl,
    `4PL` = results_4pl,
    `5PL` = results_5pl,
    `combined` = results_combined # combined results
  ))
}

# Run CV on each dataset
cv_results <- list()
for (ds in names(od)) {
  cat("\nProcessing:", ds, "\n")
  cv_results[[ds]] <- cv_all_models(od[[ds]], holdout_pct, n_iter)
}

# MSE for each iteration (mean and CI)
cv_mse <- data.frame()

for (ds in names(cv_results)) {
  for (model in c("3PL", "4PL", "5PL")) {
    res <- cv_results[[ds]][[model]]
    
    if (nrow(res) > 0) {
      mse_by_iter <- tapply(
        (res$pred_e - res$gs_e)^2,
        res$iter,
        FUN = mean,
        na.rm = TRUE
      )
      mse_by_iter <- mse_by_iter[!is.na(mse_by_iter)]
      
      # Mean of MSEs
      mean_mse <- mean(mse_by_iter)
      
      # Standard Error of MSEs
      se_mse <- se(mse_by_iter)
      
      # 95% Confidence Interval
      ci_lower <- mean_mse - 1.96 * se_mse
      ci_upper <- mean_mse + 1.96 * se_mse
      
      cv_mse <- rbind(cv_mse, data.frame(
        dataset = ds,
        model = model,
        mean_mse = mean_mse,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        n_iterations = length(unique(res$iter))
      ))
    }
  }
}

# Wide format with rounding
cv_mse$mse_ci <- paste0(
  round(cv_mse$mean_mse, 2),
  " (",
  round(cv_mse$ci_lower, 2),
  ", ",
  round(cv_mse$ci_upper, 2),
  ")"
)

cv_mse_table <- reshape(
  cv_mse[, c("dataset", "model", "mse_ci")],
  idvar = "dataset",
  timevar = "model",
  direction = "wide"
)

# > cv_mse_table
# dataset              mse_ci.3PL             mse_ci.4PL             mse_ci.5PL
# 1  cobovax_ancestral       5.64 (3.16, 8.11)      5.96 (3.34, 8.59)      2.37 (1.05, 3.69)
# 4    cobovax_omicron 407.04 (354.12, 459.97) 409.4 (356.28, 462.51) 201.69 (167.18, 236.2)
# 7       cobovax_nctd    77.18 (46.77, 107.6)  82.53 (50.05, 115.01) 129.66 (78.59, 180.73)
# 10         pattinson       4.11 (3.45, 4.77)      4.64 (3.95, 5.34)      4.21 (3.54, 4.87)

# Final model selection: cobovax_ancestral (5PL), cobovax_omicron (5PL), cobovax_nctd (3PL), pattinson (3PL)
