figure1_panel <- function(data_source, logdilution, dilution_vector, 
                          lines_color, text_x, text_y, text_main, 
                          text_panel, x_panel, y_panel, y_max, 
                          type = c("raw", "pred")) {
  
  type <- match.arg(type)
  cex_fig <- 0.5
  lwd_fig <- 0.6
  
  plot(c(logdilution[1], logdilution[length(logdilution)]), c(0, y_max), 
       bty = "n", col = "white", xlab = "", ylab = "", yaxt = "n", xaxt = "n", main = "")
  
  if (type == "raw") {
    for (v in data_source) {
      lines(logdilution, v$OD, col = lines_color, lwd = lwd_fig)
    }
  } else {
    samples <- unique(data_source$Sample)
    for (i in seq_along(samples)) {
      sample_data <- subset(data_source, Sample == samples[i])
      lines(sample_data$logDilution, sample_data$OD_fit, col = lines_color, lwd = lwd_fig)
    }
  }
  
  # X-axis
  idx <- seq(1, length(dilution_vector), by = 2) - 1
  line_x <- -0.385
  axis(side = 1, at = range(logdilution), labels = FALSE, line = line_x, lwd = lwd_fig, lwd.ticks = 0)
  axis(side = 1, at = idx, labels = FALSE, tick = TRUE, line = line_x, tck = -0.045, lwd = 0, lwd.ticks = lwd_fig)
  axis(side = 1, at = idx + 1, labels = FALSE, tick = TRUE, line = line_x, tck = -0.025, lwd = 0, lwd.ticks = lwd_fig)
  axis(side = 1, at = logdilution[idx + 1], labels = dilution_vector[idx + 1], tick = FALSE, line = -1.11, cex.axis = 0.69)
  mtext(text_x, side = 1, line = 0.8, cex = cex_fig)
  
  # Y-axis
  line_y <- -0.46
  axis(side = 2, at = 0:y_max, labels = FALSE, tick = TRUE, line = line_y, lwd.ticks = lwd_fig, lwd = lwd_fig, tck = -0.045)
  axis(side = 2, at = 0:y_max, labels = 0:y_max, tick = FALSE, line = -0.8, las = 1, cex.axis = 0.69)

  mtext(text_y, side = 2, line = 0.75, cex = cex_fig)
  
  title(main = text_main, line = 0, cex.main = 0.825)
  text(labels = text_panel, x = x_panel, y = y_panel, cex = 1.1, font = 2, xpd = NA)
}

pdf(file = "Figure 1_Fitted_curves_5PL.pdf", width = 8.1, height = 5.75)
par(mfrow = c(3, 4))
par(mar = c(2.1, 2.5, 2, 0.4), oma = c(0, 0, 0, 0), xpd = TRUE)

# Row 1: Raw data
figure1_panel(od_by_sample$cobovax_ancestral, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_ancestral"], "", "Normalized OD",
              paste0("Ancestral virus spike RBD (n=", length(unique(od$cobovax_ancestral$Sample)), ")"), "A  Raw data", -0.2, 5.9, 5, type = "raw")

figure1_panel(od_by_sample$cobovax_omicron, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_omicron"], "", "",
              paste0("Omicron BA.2 full spike (n=", length(unique(od$cobovax_omicron$Sample)), ")"), "B", -1.5, 5.9, 5, type = "raw")

figure1_panel(od_by_sample$cobovax_nctd, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_nctd"], "", "",
              paste0("N-CTD (n=", length(unique(od$cobovax_nctd$Sample)), ")"), "C", -1.5, 5.9, 5, type = "raw")

figure1_panel(od_by_sample$pattinson, logdilution_pattinson, dilution_pattinson,
              data_colors["pattinson"], "", "",
              paste0("Pattinson (n=", length(unique(od$pattinson$Sample)), ")"), "D", -0.95, 2.36, 2, type = "raw")

# Row 2: Individual model curves
figure1_panel(predicted_values_ind$cobovax_ancestral, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_ancestral"], "", "Normalized OD",
              NULL, "E Individual model curves", 2.93, 5.9, 5, type = "pred")

figure1_panel(predicted_values_ind$cobovax_omicron, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_omicron"], "", "",
              NULL, "F", -1.5, 5.9, 5, type = "pred")

figure1_panel(predicted_values_ind$cobovax_nctd, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_nctd"], "", "",
              NULL, "G", -1.5, 5.9, 5, type = "pred")

figure1_panel(predicted_values_ind$pattinson, logdilution_pattinson, dilution_pattinson,
              data_colors["pattinson"], "", "",
              NULL, "H", -0.95, 2.36, 2, type = "pred")

# Row 3: Shared model curves
figure1_panel(predicted_values_shared$cobovax_ancestral, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_ancestral"], "Reciprocal dilution", "Normalized OD",
              NULL, "I  Shared model curves", 2.2, 5.9, 5, type = "pred")

figure1_panel(predicted_values_shared$cobovax_omicron, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_omicron"], "Reciprocal dilution", "",
              NULL, "J", -1.5, 5.9, 5, type = "pred")

figure1_panel(predicted_values_shared$cobovax_nctd, logdilution_cobovax, dilution_cobovax,
              data_colors["cobovax_nctd"], "Reciprocal dilution", "",
              NULL, "K", -1.5, 5.9, 5, type = "pred")

figure1_panel(predicted_values_shared$pattinson, logdilution_pattinson, dilution_pattinson,
              data_colors["pattinson"], "Reciprocal dilution", "",
              NULL, "L", -0.95, 2.36, 2, type = "pred")

dev.off()
