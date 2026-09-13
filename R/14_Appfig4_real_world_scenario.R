set.seed(123)

# Gold-standard log titers
ec50_idv_pt <- ec50_idv$pattinson

# Start 10 worker R processes for foreach
n_cores <- 10
cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# 3pl shared-parameter model
fit_shared_3pl <- function(data, bounds = bounds_list$pattinson) {
  
  samples <- unique(data$Sample)
  n_samples <- length(samples)
  
  lowerl <- c(bounds$b[1], 
              bounds$d[1], 
              rep(bounds$e[1], n_samples))
  
  upperl <- c(bounds$b[2],     
              bounds$d[2],     
              rep(bounds$e[2], n_samples))  
  
  fit <- tryCatch({
    drm(OD ~ logDilution, 
        curveid = Sample, 
        pmodels = list(~1, ~1, ~Sample-1), 
        fct = L.3(), 
        data = data,
        lowerl = lowerl,
        upperl = upperl)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(fit)) {
    return(rep(NA, n_samples))
  }
  
  ec50 <- ED(fit, 50, display = FALSE)[ , "Estimate"]
  names(ec50) <- gsub("^e:|:50$", "", names(ec50))
  
  return(ec50)
}

n_iter <- 100
samples_pattinson <- unique(pattinson$Sample)
n_train <- 50

# Pre-generate all 100 data splits (store both train and test samples)
all_splits <- list()
for (iter in 1:n_iter) {
  train_samples <- sample(samples_pattinson, n_train, replace = FALSE)
  test_samples <- setdiff(samples_pattinson, train_samples)
  all_splits[[paste0("iter_", iter)]] <- list(
    train = train_samples,
    test = test_samples
  )
}

# 3-Point combinations
n_points_list <- c("n3")
all_comb_names <- list(
  n3 = names(sets_od$pattinson$n3)
)

all_results <- list()

# Main Loop
total_start_time <- Sys.time()

for (np in n_points_list) {
  
  cat("Processing:", np, "\n")
  
  n_comb <- length(all_comb_names[[np]])
  combo_list <- sets_od$pattinson[[np]]
  
  pb <- txtProgressBar(max = n_iter, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  iter_results <- foreach(iter = 1:n_iter, 
                          .options.snow = opts,
                          .packages = "drc", 
                          .export = c("ec50_idv_pt", "fit_shared_3pl", "bounds_list", "od", "all_splits", "combo_list","se")) %dopar% {
                            
                            res <- list(best_comb = NA, 
                                        best_mse_train = NA, 
                                        best_mse_train_ci_lower = NA, 
                                        best_mse_train_ci_upper = NA,
                                        best_gmt_train = NA,
                                        best_gmt_train_ci_lower = NA,
                                        best_gmt_train_ci_upper = NA,
                                        
                                        # testing set
                                        mse_test = NA,
                                        mse_test_ci_lower = NA,
                                        mse_test_ci_upper = NA,
                                        gmt_test = NA,
                                        gmt_test_ci_lower = NA,
                                        gmt_test_ci_upper = NA)
                            
                            # ---- Use pre-generated split ----
                            train_samples <- all_splits[[paste0("iter_", iter)]]$train
                            test_samples <- all_splits[[paste0("iter_", iter)]]$test
                            
                            train_pt <- od$pattinson[od$pattinson$Sample %in% train_samples, ]
                            test_pt <- od$pattinson[od$pattinson$Sample %in% test_samples, ]
                            
                            # ---- Training set: individual 5PL (Gold Standard) ----
                            ec50_train_gs <- ec50_idv_pt[train_samples]
                            
                            if (all(is.na(ec50_train_gs))) {
                              return(res)
                            }
                            
                            # ---- Filter training set combinations ----
                            filter_training <- function(data) {
                              data[data$Sample %in% train_samples, ]
                            }
                            
                            combo_train <- lapply(combo_list, filter_training)
                            
                            # ---- Fit 3PL shared model ----
                            ec50_train_pred <- lapply(combo_train, fit_shared_3pl)
                            ec50_train_pred <- lapply(ec50_train_pred, function(x) x[names(ec50_train_gs)])
                            
                            # ---- Compute MSE and select the optimal combination ----
                            mse_all <- sapply(ec50_train_pred, function(x) {
                              if (all(is.na(x))) return(NA)
                              mean((x - ec50_train_gs)^2)
                            })
                            
                            best_idx <- which.min(mse_all)
                            res$best_comb <- names(mse_all)[best_idx]
                            
                            best_ec50 <- ec50_train_pred[[best_idx]]
                            res$best_mse_train <- mse_all[best_idx]
                            
                            sq_err_train <- (best_ec50 - ec50_train_gs)^2
                            res$best_mse_train_ci_lower <- res$best_mse_train - 1.96 * se(sq_err_train)
                            res$best_mse_train_ci_upper <- res$best_mse_train + 1.96 * se(sq_err_train)
                            
                            res$best_gmt_train <- mean(best_ec50, na.rm = TRUE)
                            res$best_gmt_train_ci_lower <- res$best_gmt_train - 1.96 * se(best_ec50)
                            res$best_gmt_train_ci_upper <- res$best_gmt_train + 1.96 * se(best_ec50)
                            
                            # ---- Testing set: individual 5PL (Gold Standard) ----
                            ec50_test_gs <- ec50_idv_pt[test_samples]
                            
                            if (all(is.na(ec50_test_gs))) {
                              return(res)
                            }
                            
                            # ---- Filter testing set and fit the best combination ----
                            filter_testing <- function(data) {
                              data[data$Sample %in% test_samples, ]
                            }
                            
                            combo_test <- lapply(combo_list, filter_testing)
                            ec50_test_pred <- fit_shared_3pl(combo_test[[best_idx]])
                            ec50_test_pred <- ec50_test_pred[names(ec50_test_gs)]
                            
                            # ---- Compute testing set MSE and GMT ----
                            if (all(is.na(ec50_test_pred))) {
                              res$mse_test <- NA
                              res$gmt_test <- NA
                            } else {
                              sq_err_test <- (ec50_test_pred - ec50_test_gs)^2
                              res$mse_test <- mean(sq_err_test, na.rm = TRUE)
                              res$mse_test_ci_lower <- res$mse_test - 1.96 * se(sq_err_test)
                              res$mse_test_ci_upper <- res$mse_test + 1.96 * se(sq_err_test)
                              
                              res$gmt_test <- mean(ec50_test_pred, na.rm = TRUE)
                              res$gmt_test_ci_lower <- res$gmt_test - 1.96 * se(ec50_test_pred)
                              res$gmt_test_ci_upper <- res$gmt_test + 1.96 * se(ec50_test_pred)
                            }
                            
                            return(res)
                            
                          }
  
  close(pb)
  
  all_results[[np]] <- list(
    # Training set
    best_comb               = sapply(iter_results, function(x) x$best_comb),
    best_mse_train          = sapply(iter_results, function(x) x$best_mse_train),
    best_mse_train_ci_lower = sapply(iter_results, function(x) x$best_mse_train_ci_lower),
    best_mse_train_ci_upper = sapply(iter_results, function(x) x$best_mse_train_ci_upper),
    best_gmt_train          = sapply(iter_results, function(x) x$best_gmt_train),
    best_gmt_train_ci_lower = sapply(iter_results, function(x) x$best_gmt_train_ci_lower),
    best_gmt_train_ci_upper = sapply(iter_results, function(x) x$best_gmt_train_ci_upper),
    
    # Testing set
    mse_test                = sapply(iter_results, function(x) x$mse_test),
    mse_test_ci_lower       = sapply(iter_results, function(x) x$mse_test_ci_lower),
    mse_test_ci_upper       = sapply(iter_results, function(x) x$mse_test_ci_upper),
    gmt_test                = sapply(iter_results, function(x) x$gmt_test),
    gmt_test_ci_lower       = sapply(iter_results, function(x) x$gmt_test_ci_lower),
    gmt_test_ci_upper       = sapply(iter_results, function(x) x$gmt_test_ci_upper)
  )
  
  # saveRDS(all_results, file = paste0("all_results_backup_", np, ".rds"))
}

stopCluster(cl)

# total_end_time <- Sys.time()
# print(total_end_time - total_start_time)
