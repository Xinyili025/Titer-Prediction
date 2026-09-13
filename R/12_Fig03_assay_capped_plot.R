set.seed(123)

figure3_function <- function(GMT_individual,jitter_amounts,x_text,y_text,x_label,panel_label_text,main_title_text,mse){
  
  # cex
  cex_panel_fig3 <- 1.2
  cex_title_fig3 <- 1
  cex_axis_title_fig3 <- 0.8
  cex_pts_fig3 <- 0.3
  cex_axis_label_fig3 <- 0.8
  
  # lwd
  lwd_fig3 <- 1
  
  # line
  x_line <- -0.51
  y_line <- -0.537
  
  x_title_line <- 0.45
  y_title_line <- 2
  
  title_line <- 0.5
  
  # X-axis positions
  x_pos_points <- seq(0,11,3)+0.5
  x_pos_box <- seq(0,11,3)+1.5
  x_pos_label <- seq(0,11,3)+1
  x_pos_ticks <- c(-0.5,2.5,5.5,8.5,11.5)
  
  plot(c(-0.5,11.5),c(log_dilutions_pattinson(40),log_dilutions_pattinson(1280)),bty="n",col="transparent",xlab = "",ylab = "",yaxt="n",xaxt="n",main="")
  
  # Gold-standard titer
  segments(x0 = -0.5, x1 = 11.5, 
           y0 = gs_logGMT, y1 = gs_logGMT, 
           col = "gray50", lwd = lwd_fig3, lty = "dashed")
  
  for(i in 1:4) {
    
    x_pos <- jitter(rep(x_pos_points[i], 100), amount = jitter_amounts[i])
    
    points(x_pos,
           GMT_individual$logGMT[[i]],
           pch = 16, col = alpha(dp_colors[i], 0.6), cex = cex_pts_fig3)
    
    arrows(x0 = x_pos,
           y0 = GMT_individual$ci_lower[[i]],
           y1 = GMT_individual$ci_upper[[i]],
           code = 0,         
           col = alpha(dp_colors[i], 0.2),
           lwd = 0.5)
  }
  
  boxplot(GMT_individual$logGMT,at = x_pos_box,
          bty="n",xaxt="n",xlab = "",yaxt="n",ylab = "",
          frame=FALSE,add = TRUE,
          col = alpha(dp_colors, 0.5),
          outline=FALSE,whisklty = 1,medlwd = 1.75)
  
  # X-axis
  # axis
  axis(side=1,at=c(-0.5,11.5),labels=FALSE,line=x_line,lwd = 1,lwd.ticks = 0)
  # tick
  axis(side=1,at=x_pos_ticks,labels=FALSE,tick = TRUE,line=x_line,lwd.ticks = 1,lwd=0,tck=-0.045)
  # label
  axis(side=1,at=x_pos_label,labels= x_text,tick = FALSE,font = 1,line=x_line-1.2,cex.axis=cex_axis_label_fig3)
  # text
  mtext(x_label, side = 1, line = x_title_line, cex = cex_axis_title_fig3, font = 1)
  
  # Y-axis
  # axis
  axis(side=2,at=c(log_dilutions_pattinson(40),log_dilutions_pattinson(1280)),labels=FALSE,line=y_line ,lwd = 1,lwd.ticks = 0)
  # tick
  axis(side=2,at=c(log_dilutions_pattinson(c(40,80,160,320,640,1280))),labels=FALSE,tick = TRUE,line=y_line,lwd.ticks = 1,lwd=0,tck=-0.045)
  # label
  axis(side=2,at=c(log_dilutions_pattinson(c(40,80,160,320,640,1280))),labels=c(c(40,80,160,320,640,1280)),tick = FALSE,font = 1,line=y_line-0.2,las=1,cex.axis=cex_axis_label_fig3)
  # text
  mtext(y_text, side = 2, line = y_title_line, cex = cex_axis_title_fig3, font = 1)

  title(main = main_title_text, line = title_line, cex.main = cex_title_fig3, font = 2)
  
  for(i in 1:4) {
    text(labels = mse[i], 
         x = x_pos_label[i], y = 2.44, cex = 0.9, xpd = NA)
  }  
}

pdf(file = "Figure 3_GMT.pdf",width =3.5,height =3.05)
par(mfrow=c(1,1))
par(mar=c(1.5,3,0,0),oma=c(0,0,0,0),xpd=TRUE)

figure3_function(
  GMT_individual = logGMT_fixed_assays_individual,
  jitter_amounts = c(0.4,0.4,0.4,0.4),
  x_text= expression(2 %*% 60, 3 %*% 40, 4 %*% 30, 5 %*% 24),
  y_text = "GMT",
  x_label="Dilution strategy",
  panel_label_text = "B",
  main_title_text = "",
  mse = mse_fixed_assays
)

dev.off()

