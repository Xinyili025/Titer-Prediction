# Figure 2: Predicted titer vs gold-standard titer (scatter + boxplots)
set.seed(123)

# Set up Y-axis for Figure 2 panels
setup_y_axis <- function(y_axis_at, y_labels, y_line, lwd_fig2, cex_axis_legend_fig2) {
  axis(side = 2, at = y_axis_at, labels = FALSE, line = y_line, 
       lwd = lwd_fig2, lwd.ticks = 0)
  axis(side = 2, at = y_axis_at, labels = FALSE, tick = TRUE, line = y_line,
       lwd.ticks = lwd_fig2, lwd = 0, tck = -0.045)
  axis(side = 2, at = y_axis_at, labels = y_labels, tick = FALSE, 
       line = -0.8, las = 1, cex.axis = cex_axis_legend_fig2)
}

figure2_function <- function(data,
                             logdilution_vector,
                             dilution_vector,
                             log_dilutions_func,
                             original_dilutions_func,
                             p1_xtext, p2_xtext,
                             optimal_combinations, 
                             x_axis_text2,
                             panel_label_text1,panel_label_text2,
                             main_title_text){
  # Figure 2 settings
  # Font sizes
  cex_panel_fig2 <- 1.5
  cex_title_fig2 <- 1.2
  cex_axis_legend_fig2 <- 1.05
  cex_pts_fig2 <- 0.65
  cex_axis_title_fig2 <- 0.8
  lwd_fig2 <- 0.6
  
  # positions
  # Axis line positions
  x_line <- -0.45
  y_line <- -0.545
  
  # Axis title positions
  x_title_line <- 1.6
  y_title_line <- 7.2
  
  # Panel label (A, B, C, ...) positions
  x_panel_pos1 <- ifelse(length(dilution_vector) == 12, -7.2, -6.05) #ACE,G
  x_panel_pos2 <- -0.555 # BDFH
  y_panel_pos <- 1.4
  
  # Boxplot positions (right panel)
  x_pos <- c(0.5,1.5,2.5,3.5)
  
  # X-axis range
  x_min_pos <- ifelse(length(logdilution_vector)==12,-2,-2)
  x_max_pos <- ifelse(length(logdilution_vector)==12,8,6)
  
  # Y-axis range
  y_min_pos <- ifelse(length(logdilution_vector)==12,0.5,0.25)
  y_max_pos <- ifelse(length(logdilution_vector)==12,2,4)
  
  # Scatter plot
  plot(c(x_min_pos, logdilution_vector[x_max_pos] ),
       c(log_dilutions_func(y_min_pos), log_dilutions_func(y_max_pos)),
       bty = "n", col = "white",
       xlab = "", ylab = "", yaxt = "n", xaxt = "n", main = "")
  
  segments(x0 = x_min_pos,
           y0 = log_dilutions_func(1),
           x1 = logdilution_vector[x_max_pos],
           y1 = log_dilutions_func(1),
           col = "#444444", lty = 2, lwd = lwd_fig2)
  
  for(i in 1:4) {
    points(data$logtiter_gs_limited,
           data[[paste0("ratio_n", i+1, "_limited")]],
           col = alpha(dp_colors[i], 0.55),
           cex = cex_pts_fig2, pch = 16)
  }
  
  # X-axis
  axis(side = 1,
       at = c(x_min_pos,logdilution_vector[x_max_pos]),
       labels = FALSE, line = x_line, lwd = lwd_fig2, lwd.ticks = 0)
  
  # X-axis ticks 
  all_ticks <- seq(-2, x_max_pos - 1, by = 1)
  long_idx  <- seq(1, length(all_ticks), by = 2)
  short_idx <- seq(2, length(all_ticks), by = 2)
  
  tick_at_long  <- all_ticks[long_idx]
  tick_at_short <- all_ticks[short_idx]
  
  val_neg2 <- original_dilutions_func(-2)
  val_neg1 <- original_dilutions_func(-1)
  
  all_labels <- c(as.expression(bquote(italic("") <= .(val_neg2))),
                  as.expression(bquote(italic("") <= .(val_neg1))),
                  dilution_vector[1:x_max_pos])
  
  label_values <- all_labels[long_idx]
  
  axis(side = 1,
       at = tick_at_long,
       labels = FALSE, tick = TRUE, line = x_line,
       lwd.ticks = lwd_fig2, lwd = 0, tck = -0.045)
  
  axis(side = 1,
       at = tick_at_short,
       labels = FALSE, tick = TRUE, line = x_line,
       lwd.ticks = lwd_fig2, lwd = 0, tck = -0.025)
  
  # Label
  axis(side = 1,
       at = tick_at_long,
       labels = label_values,
       tick = FALSE, line = -0.8, cex.axis = cex_axis_legend_fig2 )
  
  # Text
  mtext(p1_xtext, side = 1, line = x_title_line, cex = cex_axis_title_fig2, font = 1)
  
  # Y-axis
  if(length(dilution_vector) == 12) {
    y_axis_at <- log_dilutions_func(c(0.5,2^(-0.5), 1,2^(0.5), 2))
    y_labels <- c(expression("" <= "0.50"),
                  sprintf("%.2f", 2^(-0.5)), 
                  sprintf("%.2f", 1), 
                  sprintf("%.2f", 2^(0.5)),
                  expression("" >= "2.00"))
  } else {
    y_axis_at <- log_dilutions_func(c(0.25, 0.5, 1, 2, 4))
    y_labels <- c(expression("" <= "0.25"), 
                  sprintf("%.2f", 0.5),
                  sprintf("%.2f", 1),
                  sprintf("%.2f", 2),
                  expression("" >= "4.00"))
  }
  
  setup_y_axis(y_axis_at, y_labels, y_line, lwd_fig2, cex_axis_legend_fig2)
  
  mtext("Ratio of\npredicted titer\nto\ngold-standard titer",
        side = 2, line =y_title_line, cex =cex_axis_title_fig2, font = 1, las = 1, adj = 0.5)
  
  # Panel label
  text(labels = panel_label_text1,
       x = x_panel_pos1,
       y = y_panel_pos,
       cex = cex_panel_fig2, font = 2, xpd = NA)
  
  # Title
  title(main = main_title_text, line = 0.3, cex.main = cex_title_fig2, font = 2)
  
  # Boxplot
  plot(c(0,4),c(log_dilutions_func(y_min_pos),log_dilutions_func(y_max_pos)),bty="n",col="white",xlab = "",ylab = "",yaxt="n",xaxt="n",main="")
  
  segments(x0=0, y0=log_dilutions_func(1), x1=4, y1=log_dilutions_func(1),col = "#444444", lty = 2, lwd = lwd_fig2)
  
  # Points
  for(i in 1:4) {
    points(jitter(rep(x_pos[i], length(data[[paste0("ratio_n", i+1, "_limited")]])), amount = 0.2),
           data[[paste0("ratio_n", i+1, "_limited")]],
           pch = 16, cex=cex_pts_fig2,lwd=lwd_fig2,
           col = alpha(dp_colors[i], 0.65))
  }
  
  # Boxes
  for(i in 1:4) {
    boxplot(data[[paste0("ratio_n", i+1, "_limited")]],bty="n",xaxt="n",xlab = "",yaxt="n",ylab = "",
            frame=FALSE,add = TRUE,
            at=x_pos[i],col = alpha(dp_colors[i], 0.5),lwd=lwd_fig2,
            outline=FALSE,
            whisklty = 1)
  }
  
  # X-axis
  # Axis
  axis(side=1,at=c(0,4),labels=FALSE,line=x_line,lwd = lwd_fig2,lwd.ticks = 0)
  # tick
  axis(side=1,at=c(0:4),labels=FALSE,tick = TRUE,font = 2,line=x_line,lwd.ticks = lwd_fig2,lwd=0,tck=-0.045)
  # label
  axis(side=1,at=c(x_pos),labels=x_axis_text2,tick = FALSE,font = 1,line=-1.45,cex.axis =cex_axis_legend_fig2 )
  # text
  mtext(p2_xtext, side = 1, line = x_title_line, cex = cex_axis_title_fig2, font = 1)
  
  # Y-axis
  setup_y_axis(y_axis_at, y_labels, y_line, lwd_fig2, cex_axis_legend_fig2)
  
  # Panel label
  text(labels=panel_label_text2,x=x_panel_pos2,y=y_panel_pos, cex=cex_panel_fig2,font = 2,xpd=NA)
  
  title(main = main_title_text, line = 0.3, cex.main = cex_title_fig2, font = 2)
  
  # Legend
  if(!is.null(optimal_combinations)) {
    legend(x=4.175,
           y=1.08,
           legend=optimal_combinations,x.intersp = 0.5,
           col=dp_colors[1:4],pch=16, bty="n", title="Optimal dilution combination", cex=cex_axis_legend_fig2 ,xpd =NA)}
  
}

titles <- c("Ancestral virus spike RBD", "Omicron BA.2 full spike", "N-CTD")
panel_labels1 <- c("A", "C", "E")   
panel_labels2 <- c("B", "D", "F")   

pdf(file = "Figure 2_Ratio.pdf", width = 8, height = 9)
par(mfrow = c(4, 2),
    mar = c(2.5, 2.8, 2, 3),
    oma = c(0.5, 9, 1, 10.5),
    xpd = NA)

# Cobovax
for(i in seq_along(cobovax_names)) {
  ds_name <- cobovax_names[i]
  
  figure2_function(
    data = optimal_logtiters[[ds_name]],
    logdilution_vector = logdilution_cobovax,
    dilution_vector = dilution_cobovax,
    log_dilutions_func = log2,
    original_dilutions_func = original_dilutions_cobovax,
    p1_xtext = "",
    p2_xtext = "",
    optimal_combinations = optimal_models[[ds_name]]$dilutions,
    x_axis_text2 = 2:5,
    panel_label_text1 = panel_labels1[i],
    panel_label_text2 = panel_labels2[i],
    main_title_text = titles[i]
  )
}

# Pattinson
figure2_function(
  data = optimal_logtiters$pattinson,
  logdilution_vector = logdilution_pattinson,
  dilution_vector = dilution_pattinson,
  log_dilutions_func = log4,
  original_dilutions_func = original_dilutions_pattinson,
  p1_xtext = "Gold-standard titer",
  p2_xtext = "No. of dilution points",
  optimal_combinations = optimal_models$pattinson$dilutions,
  x_axis_text2 = 2:5,
  panel_label_text1 = "G",
  panel_label_text2 = "H",
  main_title_text = "Pattinson"
)

dev.off()
