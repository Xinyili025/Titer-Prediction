get_parameter_iqr_base <- function(model_list) {
  all_results <- data.frame()
  
  for (ds_name in names(model_list)) {
    mod <- model_list[[ds_name]]
    if (is.null(mod)) next
    
    # 1. all coefficients
    coeffs <- coef(mod)
    full_names <- names(coeffs)
    
    # 2. split by parameters
    split_names <- strsplit(full_names, ":")
    params <- sapply(split_names, function(x) x[1])
    
    # 3. dataframe
    df_coef <- data.frame(
      Parameter = params,
      Value = as.numeric(coeffs),
      stringsAsFactors = FALSE
    )
    
    # Median、Q25、Q75、IQR
    unique_params <- unique(df_coef$Parameter)
    
    for (p in unique_params) {
      vals <- df_coef$Value[df_coef$Parameter == p]
      vals_clean <- vals[!is.na(vals)]
      
      med <- median(vals_clean)
      q25 <- quantile(vals_clean, 0.25)
      q75 <- quantile(vals_clean, 0.75)
      iqr_val <- IQR(vals_clean)
      
      # combine results
      res_row <- data.frame(
        Dataset = ds_name,
        Parameter = p,
        Median = med,
        Q25 = q25,
        Q75 = q75,
        IQR = iqr_val,
        Formatted = sprintf("%.3f (%.3f–%.3f)", med, q25, q75),
        stringsAsFactors = FALSE
      )
      
      all_results <- rbind(all_results, res_row)
    }
  }
  
  return(all_results)
}

param_summary_table <- get_parameter_iqr_base(ec50_idv_models)
print(param_summary_table)
# write_xlsx(param_summary_table, path = "parameter_idv.xlsx")

param_summary_table_shared <- get_parameter_iqr_base(shared_full)
# write_xlsx(param_summary_table_shared, path = "parameter_shared_3pl.xlsx")
# write_xlsx(param_summary_table_shared, path = "parameter_shared_4pl.xlsx")
# write_xlsx(param_summary_table_shared, path = "parameter_shared_5pl.xlsx")
