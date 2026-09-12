figure_ap3_function <- function(pred,logdilution,dilution_vector,lines_color,text_x,text_y,text_title,x_panel,y_panel,text_panel,x_mse, y_mse,mse,y_max) {
  cex_fig2 <- 0.5
  lwd_fig2 <- 0.6
  
  plot(c(logdilution[1],logdilution[length(logdilution)]), c(0, y_max), bty="n",col="white",xlab = "",ylab = "",yaxt="n",xaxt="n",main="")
  
  samples <- unique(pred$Sample)
  for(i in seq_along(samples)) {
    sample_data <- subset(pred, Sample == samples[i])
    lines(sample_data$logDilution, sample_data$OD_fit,
          col=lines_color,lwd=lwd_fig2)
  }
  
  # X-axis
  # axis with ticks
  idx <- seq(1, length(dilution_vector), by = 2)-1
  line_x_fig2 <- -0.385
  axis(side = 1, at = range(logdilution), labels = FALSE, 
       line = line_x_fig2,lwd = lwd_fig2,lwd.ticks = 0)
  axis(side=1,at=idx,labels=FALSE,tick = TRUE,line=line_x_fig2,tck=-0.045,lwd = 0,lwd.ticks = lwd_fig2) # longer ticks
  axis(side=1,at=idx+1,labels=FALSE,tick = TRUE,line=line_x_fig2,tck=-0.025,lwd = 0,lwd.ticks = lwd_fig2) # shorter ticks
  
  # label
  axis(side=1,at=logdilution[idx+1],labels=dilution_vector[idx+1],tick = FALSE,line=-1.11,cex.axis = 0.69)
  # text
  mtext(text_x, side = 1, line = 0.8,cex=cex_fig2)
  
  # Y-axis
  # axis with ticks
  line_y_fig2 <- -0.46
  axis(side = 2, at = 0:y_max, labels = FALSE, tick = TRUE, 
       line = line_y_fig2, lwd.ticks = lwd_fig2,lwd = lwd_fig2, tck = -0.045)
  # label
  axis(side=2,at=0:y_max,labels=c(0:y_max),tick = FALSE,line=-0.8,las=1,cex.axis=0.69)
  # text
  mtext(text_y, side = 2, line =0.75,cex=cex_fig2)
  
  # title
  title(main = text_title, line = 0.3, font=2,cex.main=0.85)
  
  # panel label
  text(labels=text_panel,x=x_panel,y=y_panel, cex=1.1,font = 2,xpd=NA)
  
  # MSE
  if(!is.null(mse)){
    text(labels = bquote(MSE == .(format(round(mse, 4), scientific = FALSE,nsmall = 4))), 
         x = x_mse, y = y_mse, cex = 0.7 ,xpd = NA)
  }
}

panel_text_fig2 <- c("A Ancestral virus spike RBD", "B", "C", "D", 
                     "E Omicron BA.2 full spike", "F", "G", "H", 
                     "I N-CTD", "J", "K", "L", 
                     "M Pattinson", "N", "O", "P")

panel_position_fig2 <- c(4.3, -0.65, -0.65, -0.65, 
                         3.7, -0.65, -0.65, -0.65, 
                         0.1, -0.65, -0.65, -0.65, 
                         0.4, -0.65, -0.65, -0.65)

pdf(file = "Figure_ap3_Fitted_curves.pdf",width =7.4,height = 7.8)
par(mfrow=c(4,4))
par(mar=c(1.8,1,2.5,0.4),oma=c(0,1,0,0),xpd=TRUE)

# Cobovax (y_max = 5) 
idx <- 1
for(ds in cobovax_names) {
  for(i in 1:4) {
    figure_ap3_function(predicted_values[[ds]][[i]], 
                     logdilution_cobovax, dilution_cobovax,
                     data_colors[ds], 
                     "",
                     if(i == 1) "Normalized OD" else "",
                     c(optimal_models[[ds]]$dilutions[1:4], "Full dilution")[i], 
                     panel_position_fig2[idx],6.2,
                     panel_text_fig2[idx],
                     9.15,4.9,if(i <= 4) optimal_models[[ds]]$MSE[i] else NULL,
                     5)
    idx <- idx + 1
  }
}

# Pattinson (y_max = 2) 
for(i in 1:4) {
  figure_ap3_function(predicted_values$pattinson[[i]], 
                   logdilution_pattinson, dilution_pattinson,
                   data_colors["pattinson"], 
                   "Reciprocal dilution",
                   if(i == 1) "Normalized OD" else "",
                   c(optimal_models$pattinson$dilutions[1:4], "Full dilution")[i], 
                   panel_position_fig2[idx],2.5,panel_text_fig2[idx], 
                   5.82,1.958,if(i <= 4) optimal_models$pattinson$MSE[i] else NULL,
                   2)
  idx <- idx + 1
}

dev.off()
