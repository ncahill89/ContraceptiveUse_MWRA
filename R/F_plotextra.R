# plots for paper (not CIs over time not diag and not bar charts...)

PlotMetDemandinCountryBySubregion <- function(# Plot country estimates of met demand by sub-region.
  ### Plot used for Alkema et al:
  ### Plot met demand in 2010 against 1990 for all countries by sub-region and 
  ### color the plotting symbols by the PPPC.
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  fig.dir = NULL, ##<< Directory to store overview plots. If NULL, folder "fig" in current working directory.
  plot.tiff  = FALSE##<< TIFF-format? If FALSE, PDF-format is used.
  ){
  
  if (is.null(fig.dir)){
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  
  load(file = file.path(output.dir,"res.country.rda")) # change JR, 20140418
  #load(file = file.path(output.dir,"res.aggregate.rda")) # change JR, 20140418
  load(file = file.path(output.dir,"mcmc.meta.rda")) # change JR, 20140418
  country.info <- mcmc.meta$data.raw$country.info
  
  est.years <- dimnames(res.country$CIratio.Lg.Lcat.qt[[1]][[1]])[[2]]
  percentiles <- dimnames(res.country$CIratio.Lg.Lcat.qt[[1]][[1]])[[1]] # 0.5
  
  probBpos.c <- unlist(lapply(res.country$changeprop.Lg.Lcat.Ti, function(l) 
  l[["Met Demand"]]["1990-2010", "PPPC"]))
  x <- 100*unlist(lapply(res.country$CIratio.Lg.Lcat.qt, function(l) 
    l$"Met Demand"[percentiles =="0.5", est.years == "1990.5"]))
  y <- 100*unlist(lapply(res.country$CIratio.Lg.Lcat.qt, function(l) 
    l$"Met Demand"[percentiles =="0.5", est.years == "2010.5"]))
  group.c <- ifelse(country.info$dev.c=="Rich", "Developed", 
                    ifelse(country.info$namereg.c == "Oceania", "Oceania", 
                           paste(country.info$namesubreg.c)))
  
  #legendprobs <- c("Pr(increase) < 0.95","Pr(increase): 0.95-0.99", "Pr(increase) > 0.99")
  #colprobs <- c("red", "darkgrey", "green")
  #col.B <- ifelse(probBpos.c>=0.95,ifelse(probBpos.c>0.99, "green","grey"), "red")
  cutoffs <- c(0.95) #c(0.5, 0.8, 0.9, 0.95, 0.99)
  legendprobs <- c(paste("Pr(increase) =<", cutoffs[1]),
        #paste("Pr(increase): ", cutoffs[-length(cutoffs)],
        #  "-", cutoffs[-1],sep = ""),
        paste("Pr(increase) >", cutoffs[length(cutoffs)]))
  #colprobs <- heat.colors(n=length(cutoffs)+1)
  #length(cutoffs)
  colprobs <- c("red", "blue") #pink", "orange", "yellow", "green", "blue")
  pchprobs <- c(25,24)
  intervals <- c(0, cutoffs, 1)
#  probBpos.c <- runif(10,0,1)
  col.B <- pch.B <- rep(NA, length(probBpos.c))
  
  for (c in 1: length(probBpos.c)){
    col.B[c] <-  colprobs[sum(probBpos.c[c] >  intervals[-length(intervals)])]
    pch.B[c] <-  pchprobs[sum(probBpos.c[c] >  intervals[-length(intervals)])]
    #      probBpos.c[c] >= intervals[-length(intervals)])]
  }
  
  ##details<< Plot \code{MetDemand19902010} is added to \code{fig.dir}.
  fig.name <- file.path(fig.dir, paste0(run.name, "MetDemand19902010")) # change JR, 20140418
  if (plot.tiff){
    tiff(filename=paste0(fig.name,".tif"), 
         width = 14, height = 10,
         units="cm",
         pointsize=6,res=300, bg="white",compression="lzw")
  } else {
    pdf(file=paste0(fig.name,".pdf"), 
        width = 14, height = 10)
  }
  groupslist <- list()
  groupslist[[2]] <- c("Central Asia", "Eastern Asia","South-Eastern Asia", 
                       "Southern Asia", "Western Asia")
  groupslist[[1]] <- c( "Eastern Africa", "Middle Africa","Northern Africa",
                         "Southern Africa", "Western Africa")
  groupslist[[3]] <- c(  "Caribbean", "Central America", "South America","Oceania", #c( "Melanesia"  , "Micronesia","Polynesia"  )          
                         "Developed")
  nf <- layout(
    rbind(
      cbind(c(1,rep(1,2)), 
            matrix(c(seq(2, 16)),3,5, byrow = TRUE)
            ),
      c(0,17,18,18,19,19))
          , 
               widths = c(0.7, rep(1.5,5)), heights = c(rep(1.5,3),0.8), TRUE)
  #layout.show(nf)
  par( mar = c(1,1,2,2))
  plot(1, xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", bty = "n")        
  text(1,1, cex = 2,las = 3,srt = 90,
         labels = "Demand satisfied 2010")
  I <- length(groupslist)
  for (i in 1:I){
    for (j in 1:length(groupslist[[i]])){
    #par(mar = c(ifelse(j==1,5,1),ifelse(i==length(groupslist),5,1),1,0))
    select <- group.c==groupslist[[i]][j]
    ys <- y[select]
    xs <- x[select]
    cols <- col.B[select]
    pchs <- pch.B[select]
    
      plot(ys ~xs, xlim = c(0,100), ylim = c(0,100),
            cex = 1.5, cex.lab = 1.5, #col  = col.B, 
            cex.axis = 1.5, type = "n",
        xlab = ifelse(#i==I | 
            j==3&i==3, "Demand satisfied 1990",""),
         ylab = ifelse(j==1, "Demand satisfied 2010",""),
         xaxt = ifelse(i==length(groupslist),
                       #| i==3 & j==4 |
           #i==2&j==5,
                       "s", "n"
           ),
         yaxt = "n", #ifelse(j==1,"s", "n"),
         cex.main = 1.5, 
         main = groupslist[[i]][j])
    abline(0,1, lwd  =3)
    if (j==1)   axis(2, at = seq(0,100,20), labels = seq(0,100,20), las = 2, cex.axis = 1.5)
    
    points(ys ~xs, col = cols, lwd = 3, cex = 1.2, pch = pchs)     
    }
  }
  plot(1, xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", bty = "n")        
  plot(1, xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", bty = "n")        
  #plot(1, xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", bty = "n")        
  text(1,1, cex = 2,las = 1,#srt = 90,
       labels = "Demand satisfied 1990")
  
  par(mar = c(0,0,0,0))
   plot(1, xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", bty = "n")        
   legend("center", cex = 2, legend = paste(legendprobs), 
          lwd = 3,  col = colprobs, pch = pchprobs, lty = -1, bty = "n")
  #plot(1, xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", bty = "n")        
#   plot(1, xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", bty = "n")        
#   text(1,1, cex = 2,las = 3,#srt = 90,
#        labels = "Demand satisfied 1990")
  
  dev.off()
  cat("Met demand in country by subregion plotted.\n")
  ##value<< NULL
  return(invisible())
}

#--------------------------------------------------------------------
BarChartSubregion <- function(#Plot counts of MWRA with unmet need by subregion for 2010, broken down by country estimates.
  ### Plot used in Alkema et al for counts of MWRA with unmet need by subregion.
  ### Note: this function was used for other things before (which explains the commented code)
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  fig.dir = NULL, ##<< Directory to store overview plots. If NULL, folder "fig" in current working directory.
  plot.tiff  = FALSE##<< TIFF-format? If FALSE, PDF-format is used.
  ){
  if (is.null(fig.dir)){
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  
  load(file = file.path(output.dir,"mcmc.meta.rda")) # change JR, 20140418
  load(file = file.path(output.dir,"res.country.rda")) # change JR, 20140418
  load(file = file.path(output.dir,"res.aggregate.rda")) # change JR, 20140418
  country.info <- mcmc.meta$data.raw$country.info
  region.info <- mcmc.meta$data.raw$region.info
  name.g <- names(res.aggregate$CIprop.Lg.Lcat.qt)
  # combine regions in Oceania
  nameselect <-  unique(ifelse(is.element(region.info$name.subreg,
                                          c("Melanesia" , "Micronesia", "Polynesia")), "Oceania", 
                               region.info$name.subreg))
  select <- is.element(name.g, nameselect)
  name.s <- name.g[select]
  n.subreg <- length(name.s)
  
  # just assume 5 as usual...
  #percentiles <- dimnames(res.aggregate$CIprop.Lg.Lcat.qt[[1]][[1]])[[1]]
  est.years <- dimnames(res.aggregate$CIprop.Lg.Lcat.qt[[1]][[1]])[[2]]
  
  ##details<< Plot \code{totalcounts} is added to \code{fig.dir}.
  fig.name <- file.path(fig.dir, paste0(run.name, "totalcounts")) # change JR, 20140418
  if (plot.tiff){
    tiff(filename=paste0(fig.name,".tif"), 
         width = 7, height = 7,
         units="cm", pointsize=6,res=300)#, bg="white")#,compression="lzw")
  } else {
    pdf(paste0(fig.name,".pdf"), width = 7, height = 7)
  }
  par(mfrow = c(1,1), mar = c(5,10,3,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  # res is in 10^6!
  res2010.sq <- From.Lg.Lcat.qtTo.gq(res.Lg.Lcat.qt = res.aggregate$CIcount.Lg.Lcat.qt,
                                 cat = "Unmet", year = "2010.5")[
                                   is.element(name.g, name.s),]/1000
  #res2015.sq <- From.Lg.Lcat.qtTo.gq(res.Lg.Lcat.qt = res.aggregate$CIcount.Lg.Lcat.qt,
  #                                   cat = "Unmet", year = "2015.5")[
  #                                     is.element(name.g, name.s),]/1000
  
  order <- order(res2010.sq[,3])
  ## Regions with less than 1 million women with unmet need are left out.
  leaveout <- sum(res2010.sq[,3]< 1)
  #ymax = max(res2015.sq, res2010.sq)*1.05
  ymax = max(res2010.sq)*1.1
  plot(seq(1, n.subreg)~res2010.sq[order,3], type = "n", ylab = "",
       xlim = c(0,ymax), xaxt = "n",
       yaxt = "n", bty = "n", 
       ylim = c(leaveout+1,n.subreg),
#       xlab = "Number of women (million)"
       xlab = "Women aged 15-49, married/in a union (million)", cex.lab = 1.2)
  axis(2, at = seq(1, n.subreg), labels = name.s[order], las = 1, cex.axis = 1)
  #abline(v = c(0,70))
  for (s in seq(10,70,20)){
    polygon(c(s,s,s+10, s+10,s), c(0,30,30,0,0), border = "NA", 
            col = adjustcolor("lightgrey", alpha.f = 0.4))
  }
  axis(1, at = round(seq(0,100,length.out = 11),0), labels = round(seq(0,100,length.out = 11),0), las = 2,
       cex.axis = 1)
  palette("default")
  #mycols <- rep(adjustcolor(palette(), alpha.f = 0.4)[1],100) # #
  orderp <- list() # country order is determined by 2010
  res2010.c <- unlist(lapply(res.country$CIcount.Lg.Lcat.qt, function(l) l[["Unmet"]][3,est.years=="2010.5"]))
  # change LA, Feb 7, 2013
  #res2015.c <- unlist(lapply(res.country$CIcount.Lg.Lcat.qt, function(l) l[["Unmet"]][3,est.years=="2015.5"]))
  for (s in seq(1, n.subreg)){
    # props of country count in region by medians
    props <- res2010.c[is.element(paste(country.info$namesubreg.c), name.s[s])]/
      sum(res2010.c[is.element(paste(country.info$namesubreg.c), name.s[s])])
    orderp[[s]] <- order(props, decreasing = T)
  }
  width <- 0.15
  res.Lt.sq <- list(res2010.sq = res2010.sq)#, res2015.sq = res2015.sq)
  res.Lt.c <- list(res2010.c = res2010.c)#, res2015.c = res2015.c)
  # change LA
  for (t in 1:1){
  #for (t in 1:2){
    res.sq <- res.Lt.sq[[t]]
    res.c <- res.Lt.c[[t]]
    year <- c(2010.5, 2015.5)[t]
    #add <- c(0.15, -0.15)[t]
    add <- 0#c(0.15, -0.15)[t]
    for (j in seq(1, n.subreg)){
      s <- order[j] # get right region for plot index j
      # props of country count in region by medians
      # no longer used, not sure if it still work (July 2012)
#      props <- res.c[is.element(paste(country.info$namesubreg.c), name.s[s])]/
#        sum(res.c[is.element(paste(country.info$namesubreg.c), name.s[s])])
      # use country order from 2010
#      cprops <- c(0, cumsum(props[orderp[[s]]]))
#       for (i in 2:length(cprops)){
#         polygon(res.sq[s,3]*c(cprops[i-1],cprops[i],cprops[i],cprops[i-1],cprops[i-1]), 
#                 add + j + c(-width,-width,width, width, -width), 
#                 border = "NA", col = mycols[i])
#       }
#       # put box around
      polygon(c(0, res.sq[s,3],res.sq[s,3], 0,0), 
              add + j + c(-width,-width,+width, +width,-width), border = 1,
              col = adjustcolor("red", alpha.f = 0.4)[1])
    } # end s
    points(add+seq(1, n.subreg)~res.sq[order,3], pch = 20, lwd = ifelse(plot.tiff,1,2))
    segments(res.sq[order,1], add+seq(1, n.subreg), res.sq[order,5], 
           add+seq(1, n.subreg), lwd = ifelse(plot.tiff,1,2), col = 1)
  } # end t
  dev.off()
  cat("Counts of MWRA with unmet need by subregion broken down by country estimates plotted.\n")
  ##value<< NULL, 
  return(invisible())
}


#--------------------------------------------------------------------
CIPropChangesSubregions <- function( # Plot CIs for unmet and total for subregions in one plot.
  ### Plot used in Alkema et al:  Plot CIs for unmet and total for subregions in one plot.
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  fig.dir = NULL, ##<< Directory to store overview plots. If NULL, folder "fig" in current working directory.
  plot.tiff  = FALSE##<< TIFF-format? If FALSE, PDF-format is used.
  ){
  
  if (is.null(fig.dir)){
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  
  #load(file = file.path(output.dir,"res.country.rda")) # change JR, 20140418
  load(file = file.path(output.dir,"res.aggregate.rda")) # change JR, 20140418
  load(file = file.path(output.dir,"mcmc.meta.rda")) # change JR, 20140418
  percentiles <- dimnames(res.aggregate$CIprop.Lg.Lcat.qt[[1]][[1]])[[1]]
  est.years <- dimnames(res.aggregate$CIprop.Lg.Lcat.qt[[1]][[1]])[[2]]
  
  DoPlot <- function(){
    year <- "2010.5"
    cat <- "Total"
    med2010 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.5", est.years = year]))[
        is.element(name.g, name.s)]
    low2010 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.025", est.years = year]))[
        is.element(name.g, name.s)]
    up2010 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.975", est.years = year]))[
        is.element(name.g, name.s)]
    year <- "1990.5"
    med1990 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.5", est.years = year]))[
        is.element(name.g, name.s)]
    low1990 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.025", est.years = year]))[
        is.element(name.g, name.s)]
    up1990 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.975", est.years = year]))[
        is.element(name.g, name.s)]
    order <- order(med2010)
    for (s in seq(2,n.subreg,2)){
      polygon(c(0,0,1,1,0), c(s-0.5, s+0.5, s+0.5, s-0.5, s-0.5), border = "NA", 
              col = adjustcolor("lightgrey", alpha.f = 0.4))
    }
    segments(rep(0, n.subreg), seq(1, n.subreg)-0.5,
             rep(1, n.subreg),seq(1, n.subreg)-0.5)
    #  segments(rep(-1, n.subreg), seq(1, n.subreg)+0.5,rep(1, n.subreg),seq(1, n.subreg)+0.5)
    
    segments(low1990[order], seq(1+.25, n.subreg+0.25), up1990[order], 
             seq(1+0.25, n.subreg+0.25), lwd= 2*ifelse(plot.tiff,1,2), col = "turquoise")
    segments(low2010[order], seq(1+0.1, n.subreg+0.1), up2010[order], col = "blue", 
             seq(1+0.1, n.subreg+0.1), lwd= 2*ifelse(plot.tiff,1,2))
    segments(low1990[order], seq(1+.25, n.subreg+0.25), up1990[order], 
             seq(1+0.25, n.subreg+0.25), lwd= 1, col =1)
    segments(low2010[order], seq(1+0.1, n.subreg+0.1), up2010[order], col = 1, 
             seq(1+0.1, n.subreg+0.1), lwd= 1)
    
    med2010t <- med2010
    med1990t <- med1990
    
    # add unmet
    year <- "2010.5"
    cat <- "Unmet"
    med2010 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.5", est.years = year]))[
        is.element(name.g, name.s)]
    low2010 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.025", est.years = year]))[
        is.element(name.g, name.s)]
    up2010 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.975", est.years = year]))[
        is.element(name.g, name.s)]
    year <- "1990.5"
    med1990 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.5", est.years = year]))[
        is.element(name.g, name.s)]
    low1990 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.025", est.years = year]))[
        is.element(name.g, name.s)]
    up1990 <- unlist(lapply(res.aggregate$CIprop.Lg.Lcat.qt, function(l)
      l[[cat]][percentiles=="0.975", est.years = year]))[
        is.element(name.g, name.s)]
    segments(low1990[order], seq(1-.1, n.subreg-0.1), up1990[order], 
             seq(1-0.1, n.subreg-0.1), lwd= 2*ifelse(plot.tiff,1,2), col = "pink")
    segments(low2010[order], seq(1-0.25, n.subreg-0.25), up2010[order], col = "red", 
             seq(1-0.25, n.subreg-0.25), lwd=2*ifelse(plot.tiff,1,2))
    
    segments(low1990[order], seq(1-.1, n.subreg-0.1), up1990[order], 
             seq(1-0.1, n.subreg-0.1), lwd= 1, col = 1)
    segments(low2010[order], seq(1-0.25, n.subreg-0.25), up2010[order], col = 1, 
             seq(1-0.25, n.subreg-0.25), lwd=1)

    points(seq(1+0.25, n.subreg+0.25)~ med1990t[order], cex = 1.2, pch = 19, col = "turquoise")
    points(seq(1+0.25, n.subreg+0.25)~ med1990t[order],  pch = 1, cex = 1.2)
    points(seq(1+.1, n.subreg+0.1)~ med2010t[order], cex = 1.2, pch = 19, col = "blue")
    points(seq(1+.1, n.subreg+0.1)~ med2010t[order],  pch = 1,cex = 1.2)
    points(seq(1-.1, n.subreg-0.1)~ med1990[order], cex = 1.2, pch = 19, col = "pink")
    points(seq(1-.1, n.subreg-0.1)~ med1990[order],  pch = 1,cex = 1.2)
    points(seq(1-.25, n.subreg-0.25)~ med2010[order], cex = 1.2, pch = 19, col = "red")
    points(seq(1-.25, n.subreg-0.25)~ med2010[order],  pch = 1,cex = 1.2)
  
    axis(2, pos = 0, at = seq(1, n.subreg),  tick = T, 
         labels =name.s[order], las = 1, cex.axis = 1)
  }
  
  ##details<< Figure \code{CIspropsubregional} is added to \code{fig.dir}.
  fig.name <- file.path(fig.dir, paste0(run.name, "CIspropsubregional"))
  if (plot.tiff){
    tiff(filename=paste0(fig.name,".tif"),
         width = 8, height = 9,
         units="cm", pointsize=6,res=300, bg="white",compression="lzw")
  } else {
    pdf(paste0(fig.name,".pdf"), width = 8, height = 9)
  }
  
  # first UNDP subregional est's
  region.info <- mcmc.meta$data.raw$region.info
  name.g <- names(res.aggregate$CIprop.Lg.Lcat.qt)
  # combine regions in Oceania
  nameselect <-  unique(ifelse(is.element(region.info$name.subreg,
                                          c("Melanesia" , "Micronesia", "Polynesia")), 
                               "Mela-Micro-Polynesia", #"Oceania", 
                               region.info$name.subreg))
  select <- is.element(name.g, nameselect)
  name.s <- name.g[select]
  n.subreg <- length(name.s)
  
  n.plot <- ( # temp assignment to get axis for plot right
      length(nameselect)
      + 1.5 # world
      + 3.5 # dev etc
      )
  par(mfrow = c(1,1), mar = c(5,12,3,1), cex.main = 2, cex.axis = 2, cex.lab = 2)
  plot(1, ylab = "", bty = "n", 
       xlab = "Women aged 15-49, married/in a union (%)", cex.lab = 1.2,
       xlim = c(0,1),  xaxt= "n", yaxt = "n", ylim = c(0, n.plot))
  axis(1, at = seq(0,1,0.1),labels = seq(0,1,0.1)*100, las = 2, cex.axis = 1.2)
  

  DoPlot()
  center <- n.subreg+0.5
  segments(0,center, 1,center)
  i =0
  for (nameadd in c("Developing (excl. China)", "Developing countries",  "Developed countries","World")){
  i = i+1
  center <- center+1
  s <- 0
  if (i==4) center = center+0.5
  if (i!=2){
    polygon(c(0,0,1,1,0), center+c(s-0.5, s+0.5, s+0.5, s-0.5, s-0.5), border = "NA", col = "lightgrey")
    segments(0,center-0.5, 1,center-0.5)
  }
  axis(2, pos = 0, at = center, labels = nameadd, las = 1, cex.axis = 1)
  segments(0,center-0.5, 1,center-0.5)
  segments(0,center+0.5, 1,center+0.5)
  
  med9u <- res.aggregate$CIprop.Lg.Lcat.qt[[nameadd]][["Unmet"]][percentiles=="0.5", est.years == "1990.5"]
  CI9u.q2 <- res.aggregate$CIprop.Lg.Lcat.qt[[nameadd]][["Unmet"]][is.element(percentiles, c("0.025", "0.975")), 
                                                                   est.years == "1990.5"]
  
  med2010u <- res.aggregate$CIprop.Lg.Lcat.qt[[nameadd]][["Unmet"]][percentiles=="0.5", est.years == "2010.5"]
  
  CI2010u.q2 <- res.aggregate$CIprop.Lg.Lcat.qt[[nameadd]][["Unmet"]][is.element(percentiles, c("0.025", "0.975")), 
                                                                      est.years == "2010.5"]
  segments(CI9u.q2[1], center-0.1, CI9u.q2[2], center-0.1, col = "pink", lwd=2*ifelse(plot.tiff,1,2))
  segments(CI2010u.q2[1], center-0.25, CI2010u.q2[2], center-0.25, col = 2, lwd=2*ifelse(plot.tiff,1,2))
  segments(CI9u.q2[1], center-0.1, CI9u.q2[2], center-0.1, col = 1, lwd=1)
  segments(CI2010u.q2[1], center-0.25, CI2010u.q2[2], center-0.25, col = 1, lwd=1)
  points(center-0.1~ med9u, cex = 1.2, pch = 19, col = "pink")
  points(center-0.1~ med9u,pch = 1,cex = 1.2)
  points(center-0.25~ med2010u, cex = 1.2, pch = 19, col = 2)
  points(center-0.25~ med2010u, pch = 1,cex = 1.2)
  
  med9u <- res.aggregate$CIprop.Lg.Lcat.qt[[nameadd]][["Total"]][percentiles=="0.5", est.years == "1990.5"]
  CI9u.q2 <- res.aggregate$CIprop.Lg.Lcat.qt[[nameadd]][["Total"]][is.element(percentiles, c("0.025", "0.975")), est.years == "1990.5"]
  med2010u <- res.aggregate$CIprop.Lg.Lcat.qt[[nameadd]][["Total"]][percentiles=="0.5", est.years == "2010.5"]
  CI2010u.q2 <- res.aggregate$CIprop.Lg.Lcat.qt[[nameadd]][["Total"]][is.element(percentiles, c("0.025", "0.975")),est.years == "2010.5"]
  
  segments(CI9u.q2[1], center+0.25, CI9u.q2[2], center+0.25, col = "turquoise", lwd=2*ifelse(plot.tiff,1,2))
  segments(CI2010u.q2[1], center+0.1, CI2010u.q2[2], center+0.1, col = "blue", lwd = 2*ifelse(plot.tiff,1,2))
  
  segments(CI9u.q2[1], center+0.25, CI9u.q2[2], center+0.25, col =1, lwd=1)
  segments(CI2010u.q2[1], center+0.1, CI2010u.q2[2], center+0.1, col = 1, lwd = 1)
  points(center+0.25~ med9u, cex = 1.2, pch = 19, col = "turquoise")
  points(center+0.25~ med9u, pch = 1,cex = 1.2)
  points(center+0.1~ med2010u, cex = 1.2, pch = 19, col = "blue")
  points(center+0.1~ med2010u,pch = 1,cex = 1.2)
  }  
  
  
  abline(v = 0)
  abline(v=1)
  segments(center+1, 0, center+1, 1)
  #box()
  x = legend(#1,0, xjust = 1, yjust = 1, #
    "bottomright", plot = F,
         legend = c("Any method 1990", "Any method 2010",
                                   "Unmet need 1990", "Unmet need 2010"), bg = "white",
         col = c("turquoise","blue", "pink", 2), lwd =3, pch = 19, lty = 1, cex = 1.1)
  
  legend(x$rect$left, x$rect$top,
    legend = c("Any method 1990", "Any method 2010",
               "Unmet need 1990", "Unmet need 2010"), bg = "white",
    col = c("turquoise","blue", "pink", 2), lwd =3, pch = 19, lty = 1, cex = 1.1)
  
#   legend(x$rect$left, x$rect$top,
#          legend = rep("                                        0919",4), bty = "n",
#     col = 1, lwd =1, pch = 1, lty = 1, cex = 1.2)

  dev.off()
  cat("CIs for unmet and total for subregions plotted.\n")
  ##value<< NULL
  return(invisible())
}

#-----------------------------------------------------------------------------------------------
PlotLogisticParameters <- function (# Plot overview of country parameters of the logistic trends
  ### Plot overview of country parameters of the logistic trends: scatter plots of medians and CIs
  par.ciq, ##<< Output from \code{\link{GetParInfo}}
  country.info, region.info,
  fig.name = NULL) {
  
  ##details<< Note: pace parameter omega of logistic curve is expressed as 
  ## the no of years needed for an increase of 60% on the 0-asymptote scale
  #perc80 <- (logit(0.9) - logit(0.1))
  perc60 <- (logit(0.8) - logit(0.2))
  perc <- perc60
  yearsomega <- perc/par.ciq[,"omega.c",2]
  yearsomega.low <- perc/par.ciq[,"omega.c",3]
  yearsomega.up <- perc/par.ciq[,"omega.c",1]
  
  yearsRomega <- perc/par.ciq[,"Romega.c",2]
  yearsRomega.low <- perc/par.ciq[,"Romega.c",3]
  yearsRomega.up <- perc/par.ciq[,"Romega.c",1]
  
  t0low <- par.ciq[,"T.c",1]
  t0up <- par.ciq[,"T.c",3]
  Rt0low <- par.ciq[,"RT.c",1]
  Rt0up <- par.ciq[,"RT.c",3]
  Rt0.c <- par.ciq[,"RT.c" ,2]
  t0.c <- par.ciq[,"T.c",2]
  omega.c <- par.ciq[,"omega.c",2]
  Romega.c <-par.ciq[,"Romega.c",2]
  pmax.c <- par.ciq[,"pmax.c",2]
  Rmax.c <- par.ciq[,"Rmax.c",2]
  
  if (!is.null(fig.name)) {
    pdf(fig.name, width = 14, height = 8)
  }
  #----------------------------------------------------------------------------------
  # omega against Tmid  
  par(mfrow = c(1,2),  mar = c(5,5,5,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  plot(yearsomega~t0.c, type = "n", main = "Latent trend: Total prevalence",
       ylab = "# years needed for 60% increase", xlab = "Midpoint increase",
       xlim = c(min(t0.c, 1850),max(t0.c,2050)), ylim = c(0,max(yearsomega, 250)))
  text(perc/omega.c~t0.c, labels = country.info$code.c, col =1+as.numeric(country.info$reg.c))
  legend("topleft", legend = region.info$name.reg.short, col = seq(2, 1+region.info$n.reg), lwd = 3, 
         pch = 1, lty = -1)
  plot(yearsomega~t0.c, type = "n", main = "Latent trend: Total prevalence",
       ylab = "# years needed for 60% increase", xlab = "Midpoint increase",
       xlim = c(min(t0.c, 1850),max(t0.c,2050)), ylim = c(0,max(yearsomega, 250)))
  segments(t0.c, yearsomega.low, t0.c, yearsomega.up, col =1+as.numeric(country.info$reg.c))
  segments(t0low, yearsomega, t0up, yearsomega, col =1+as.numeric(country.info$reg.c))
  #----------------------------------------------------------------------------------
  par(mfrow = c(1,2),  mar = c(5,5,5,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  plot(yearsRomega~Rt0.c, type = "n", main = "Latent trend: Ratio",
       ylab = "# years needed for 60% increase", xlab = "Midpoint increase",
       xlim = c(min(Rt0.c, 1850),max(Rt0.c,2050)), ylim = c(0,max(yearsRomega, 250)))
  text(perc/Romega.c~Rt0.c, labels = country.info$code.c, col =1+as.numeric(country.info$reg.c))
  legend("topleft", legend = region.info$name.reg.short, col = seq(2, 1+region.info$n.reg), lwd = 3, 
         pch = 1, lty = -1)
  #----------------------------------------------------------------------------------
  plot(yearsRomega~Rt0.c, type = "n", main = "Latent trend: Ratio",
       ylab = "# years needed for 60% increase", xlab = "Midpoint increase",
       xlim = c(min(Rt0.c, 1850),max(Rt0.c,2050)), ylim = c(0,max(yearsRomega, 250)))
  segments(Rt0.c, yearsRomega.low, Rt0.c, yearsRomega.up, col =1+as.numeric(country.info$reg.c))
  segments(Rt0low, yearsRomega, t0up, yearsRomega, col =1+as.numeric(country.info$reg.c))
  #----------------------------------------------------------------------------------
  # asymptotes
  par(mfrow = c(1,2),  mar = c(5,5,5,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  plot(pmax.c~Rmax.c, type = "n", main = "Asymptotes",
       ylab = "Total prevalence", xlab = "Ratio",
       xlim = c(0.5,1), ylim = c(0.5,1))
  text(pmax.c~Rmax.c, labels = country.info$code.c, col =1+as.numeric(country.info$reg.c))
  legend("topleft", legend = region.info$name.reg.short, col = seq(2, 1+region.info$n.reg), lwd = 3, 
         pch = 1, lty = -1)
  plot(pmax.c~Rmax.c, type = "n", main = "Asymptotes",
       ylab = "Total prevalence", xlab = "Ratio",
       xlim = c(0.5,1), ylim = c(0.5,1))
  segments(Rmax.c, par.ciq[,"pmax.c",1],
           Rmax.c, par.ciq[,"pmax.c",3], col = 1+as.numeric(country.info$reg.c))
  segments( par.ciq[,"Rmax.c",1],pmax.c,
            par.ciq[,"Rmax.c",3], pmax.c, col =1+as.numeric(country.info$reg.c))
  #  dev.off()
  if (!is.null(fig.name)) {
    dev.off()
  }
  ##value<< NULL
  return(invisible())
}
#----------------------------------------------------------------------
PlotCountryEstimatesForAggregate <- function (# Create overview country estimates for aggregates
  ### Create overview plots of estimates of proportions/counts over time for
  ### countries within aggregates.
  CI.Lg.Lcat.qt, ##<< Object from class \code{CI.Lg.Lcat.qt}, either a proportion or a count (see next).
  fig.name = NULL, ## If NULL, plot appears in R, else it is saved as fig.name.
  cats = NULL, ## If NULL, all cats from the CI are plotted. Alternatively, subset of those cats
  country.info  ##Object of class \code{\link{country.info}}
  ){
  
  if (is.null(cats)){
    cats <- names(CI.Lg.Lcat.qt[[1]])
  } else {
    cats <- cats[is.element(cats, names(CI.Lg.Lcat.qt[[1]]))]
  } 
  nplot <- ceiling((length(cats)+1)/2)
  name.short.c <- InternalMakeCountryNamesShort(country.info$name.c)
  if (!is.null(fig.name)){
    pdf(fig.name, width = 7, height = 5)
  }
  for (subreg in 1:max(country.info$subreg.c)){
    select.c <- seq(1, length(country.info$name.c))[country.info$subreg.c == subreg]
    par(mfrow = c(2,nplot), mar = c(5,5,1,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
    for (cat in cats){
      InternalPlotCountryEstimatesForAggregate(
        CI.Lg.Lcat.qt = CI.Lg.Lcat.qt, cat = cat,
        select.c = select.c)
    }
    par( mar = c(0,0,5,0), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
    plot(1, type = "n", xaxt = "n", xlab = "", ylab = "", yaxt = "n",
         main = country.info$namesubreg.c[select.c[1]])
    legend("center", legend = name.short.c[country.info$subreg.c==subreg], 
           col = seq(1, sum(country.info$subreg.c==subreg)), 
           lty = seq(1, sum(country.info$subreg.c==subreg)),
           lwd = 3, cex = 0.8, bty = "n")
  } # end subregs
  if (!is.null(fig.name)){
    dev.off()  
  }
  return()
}

InternalPlotCountryEstimatesForAggregate <- function (# Create overview country estimates for aggregates
  ### Create overview plots of estimates of proportions/counts over time for
  ### countries within aggregates.
  CI.Lg.Lcat.qt, ##<< Object from class \code{CI.Lg.Lcat.qt}, either a proportion or a count (see next).
  cat, select.c
  ){
  est.years <- as.numeric(names(CI.Lg.Lcat.qt[[1]][[1]][1,]))
  name.g <- names(CI.Lg.Lcat.qt)
  plot(1, type = "n", xlim = c(min(est.years), max(est.years)), ylim = c(0,1), xlab = "Year", 
       ylab = paste(cat, "(%)"))
  col<-0
  lty = 0
  for (c in select.c){
    col <- col+1
    lty = lty+1
    lines(CI.Lg.Lcat.qt[[name.g[c]]][[cat]][3,]~est.years, lwd = 3, col = col, lty = lty)
  }    
}
#----------------------------------------------------------------------
# The End!
