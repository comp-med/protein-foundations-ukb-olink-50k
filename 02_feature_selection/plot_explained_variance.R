#######################################
## function to plot explained variance
## for specific candidate proteins

plot.expl.var <- function(res, id, lab, single=T){
  
  ## 'res' -- results data frame (needs to be in long format)
  ## 'id'  -- id of the protein to plot
  ## 'lab' -- label for the categories
  
  ## subset to variable of interest
  res <- res[ var == id]
  
  ## sort by p.50
  res <- res[order(-p.50)]
  print(res)
  ## ?restrict to atmost 20 variables (sum the rest up)
  
  ## drop residuals
  res <- res[-which(variable == "Residuals")]
  
  ## divide the plot with a overall bar on top and detailed results underneath (excluding residuals)
  par(mar=c(1.5,1.5,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, mgp=c(.6,0,0), lwd=.5)
  if(single) layout(matrix(1:2,2,1), heights = c(.2,.8))
  
  #------------------------------------#
  ##--   bar overall contribution   --##
  #------------------------------------#
  
  ## empty plot
  plot(c(0,1), c(0,1), type="n", xaxt="n", yaxt="n", ylab="", xlab="Explained variance [%]")
  axis(1, at=seq(0,1,.1), labels=seq(0,100,10), lwd=.5)
  
  ## add explained variance
  for(j in nrow(res):1){
    rect(0, 0, sum(res$p.50[1:j]), 1, lwd=.1, border="white", col=sapply(res$variable[j], function(x) lab$cl[which(lab$short_name == x)]))
  }
  
  #------------------------------------#
  ##-- bar plot excluding residuals --##
  #------------------------------------#
  
  ## adopt graphical parameters
  par(mar=c(6,1.5,.5,.5), cex.axis=.5, cex.lab=.5, tck=-.01, mgp=c(.6,0,0))
  
  ## do the plotting
  plot(c(.5,nrow(res)+.5), c(0, max(res$p.975)*100), type="n", xlab="", ylab="Explained variance [%]", xaxt="n", yaxt="n")
  ## add axis
  axis(2, lwd=.5)

  ## add explained variance
  rect(1:nrow(res)-.4, 0, 1:nrow(res)+.4, res$p.50*100, lwd=.1, border="grey20", 
       col=sapply(res$variable, function(x){
         if(x %in% lab$short_name){
           lab$cl[which(lab$short_name == x)]
         }else{
           return("white")
         }
       }))
  
  ## add confidence intervals
  arrows(1:nrow(res), res$p.025*100, 1:nrow(res), res$p.975*100, length = 0, lwd=.2)
  
  ## add label (get plotting coordinates first)
  pm <- par("usr")
  text(1:nrow(res), pm[3]-(pm[4]-pm[3])*.05, srt=60, cex=.4, xpd=NA, pos=2, offset=0,
       labels=sapply(res$variable, function(x){
         if(x %in% lab$short_name){
           lab$label[which(lab$short_name == x)]
         }else{
           return(x)
         }
       }))
  
  
  
}
