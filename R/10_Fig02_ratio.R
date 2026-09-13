log4 <- function(x) log(x, base = 4)

# Gold-standard titers at or below -2 are plotted at -2
# Cobovax
for(ds in cobovax_names) {
  optimal_logtiters[[ds]]$logtiter_gs_limited <- ifelse(
    optimal_logtiters[[ds]]$logtiter_gs <= -2, -2, 
    optimal_logtiters[[ds]]$logtiter_gs
  )
}

# Pattinson
optimal_logtiters[["pattinson"]]$logtiter_gs_limited <- ifelse(optimal_logtiters[["pattinson"]]$logtiter_gs <= -2, -2, optimal_logtiters[["pattinson"]]$logtiter_gs)

# Cobovax: clip ratio to [0.5, 2] on the log2 scale
optimal_logtiters[cobovax_names] <- lapply(optimal_logtiters[cobovax_names], function(df) {
  for(n in 2:5) {
    ratio_col <- paste0("ratio_n", n)
    lim_col <- paste0("ratio_n", n, "_limited")
    df[[lim_col]] <- pmax(log2(0.5), pmin(df[[ratio_col]], log2(2)))
  }
  return(df)
})

# Pattinson: clip ratio to [0.25, 4] on the log4 scale
for(n in 2:5) {
  ratio_col <- paste0("ratio_n", n)
  lim_col <- paste0("ratio_n", n, "_limited")
  optimal_logtiters[["pattinson"]][[lim_col]] <- pmax(log4(0.25), 
                                                      pmin(optimal_logtiters[["pattinson"]][[ratio_col]], log4(4)))
}
