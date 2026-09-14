set.seed(123)

# Data Preparation for Panel A
comb_freq <- table(all_results$n3$best_comb)
comb_freq_sorted <- sort(comb_freq, decreasing = TRUE)

# Select combinations for plotting
comb_freq_top <- comb_freq_sorted[comb_freq_sorted >= 6]

comb_names <- names(comb_freq_top)
comb_counts <- as.numeric(comb_freq_top)
n_combs <- length(comb_names)

dilutions <- c("40", "160", "640", "2560", "10240", "40960", "163840", "655360")

right_labels <- paste0("N = ", comb_counts)

# Data preparation for Panel B
# Extract information in training set
rw_training <- list(
  logGMT = list(n3 = all_results$n3$best_gmt_train),
  ci_lower = list(n3 = all_results$n3$best_gmt_train_ci_lower),
  ci_upper = list(n3 = all_results$n3$best_gmt_train_ci_upper)
)
rw_training_mse <- round(mse(all_results$n3$best_gmt_train,gs_logGMT),4)

# Extract information in testing set
rw_testing <- list(
  logGMT = list(n3 = all_results$n3$gmt_test),
  ci_lower = list(n3 = all_results$n3$gmt_test_ci_lower),
  ci_upper = list(n3 = all_results$n3$gmt_test_ci_upper)
)
rw_testing_mse <- round(mse(all_results$n3$gmt_test,gs_logGMT),4)

# Appendix Figure 4
pdf(file = "Appfigure 4_Real_world_scenario.pdf", width=11, height=3.3)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.8, 1))
par(mar=c(2,2.5,2,2.5), oma=c(0,0,0,0), xpd=TRUE)

# Panel A
plot(1, 1, type = "n", 
     xlim = c(0.5, 8.5), 
     ylim = c(0.5, n_combs + 0.5),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "",
     main = "")

for (i in 1:n_combs) {
  
  comb_str <- comb_names[i]
  idx <- as.numeric(unlist(strsplit(gsub("c\\(|\\)", "", comb_str), ",")))
  y_pos <- n_combs - i + 1
  
  for (j in 1:8) {
    is_selected <- j %in% idx
    fill_col <- ifelse(is_selected, "#1976D2", "white")
    
    rect(xleft = j - 0.5, ybottom = y_pos - 0.5, 
         xright = j + 0.5, ytop = y_pos + 0.5, 
         col = fill_col, border = "black")
  }
}

axis(1, at = 1:8, labels = dilutions, las = 1, cex.axis = 1, tick = FALSE, line = -1.35)
mtext("Reciprocal dilution", side = 1, line = 0.8, font = 1)

# Y-axis
axis(2, at = n_combs:1, labels = right_labels, las = 1, cex.axis = 1, tick = FALSE, line = -1.78)
text(-0.42,8.55,labels = "A",font = 2,xpd=NA,cex = 1.2)

# Panel B
plot(c(0, 2), c(log_dilutions_pattinson(80), log_dilutions_pattinson(1280)),
     bty="n", col="transparent", xlab="", ylab="", yaxt="n", xaxt="n", main="")

# Gold-standard titer 
segments(x0 = 0, x1 = 2, 
         y0 = gs_logGMT, y1 = gs_logGMT, 
         col = "gray50", lwd = 1, lty = "dashed")

# Training set
x_pos <- jitter(rep(0.35, n_iter), amount = 0.1)

points(x_pos, rw_training$logGMT[[1]], pch = 16, col = alpha(dp_colors["n3"], 0.6), cex = 0.5)
arrows(x0 = x_pos, y0 = rw_training$ci_lower[[1]], y1 = rw_training$ci_upper[[1]],
       code = 0, col = alpha(dp_colors["n3"], 0.2), lwd = 0.5)

boxplot(rw_training$logGMT, at = 0.65,boxwex = 0.4,
        bty="n", xaxt="n", xlab="", yaxt="n", ylab="",
        frame=FALSE, add=TRUE, col=alpha(dp_colors["n3"], 0.5),
        outline=FALSE, whisklty=1, medlwd=1.75)

x_pos <- jitter(rep(1.35, n_iter), amount = 0.1)

points(x_pos, rw_testing$logGMT[[1]], pch = 16, col = alpha(dp_colors["n3"], 0.6), cex = 0.5)
arrows(x0 = x_pos, y0 = rw_testing$ci_lower[[1]], y1 = rw_testing$ci_upper[[1]],
       code = 0, col = alpha(dp_colors["n3"], 0.2), lwd = 0.5)

boxplot(rw_testing$logGMT, at = 1.65,boxwex = 0.4,
        bty="n", xaxt="n", xlab="", yaxt="n", ylab="",
        frame=FALSE, add=TRUE, col=alpha(dp_colors["n3"], 0.5),
        outline=FALSE, whisklty=1, medlwd=1.75)

axis(side=1, at=c(0, 2), labels=FALSE, lwd=1, lwd.ticks=0,line = -0.46)
axis(side=1, at=c(0,1,2), labels=FALSE, tick=TRUE, lwd.ticks=1, lwd=0, tck=-0.07,line=-0.46)
axis(side=1, at=c(0.5, 1.5), labels=c("Training sets", "Testing sets"), tick=FALSE, line= -1.5, cex.axis=1)
mtext("Data sets", side=1, line=0.8, font=1)

y_ticks <- log_dilutions_pattinson(c(80, 160, 320, 640, 1280))
axis(side=2, at=c(log_dilutions_pattinson(80), log_dilutions_pattinson(1280)), labels=FALSE, lwd=1, lwd.ticks=0,line=-0.54)
axis(side=2, at=y_ticks, labels=FALSE, tick=TRUE, lwd.ticks=1, lwd=0, tck=-0.05,line=-0.54)
axis(side=2, at=y_ticks, labels=c(80, 160, 320, 640, 1280), tick=FALSE, las=1, line=-0.8, cex.axis=1)
mtext("GMT", side=2, line=2, font=1)

text(x = 0.5, y = log_dilutions_pattinson(1200), 
     labels = rw_training_mse, cex = 1)
text(x = 1.5, y = log_dilutions_pattinson(1200), 
     labels = rw_testing_mse, cex = 1)

text(-0.65,log_dilutions_pattinson(1280)+0.3,labels = "B",font = 2,xpd=NA,cex = 1.2)

dev.off()
