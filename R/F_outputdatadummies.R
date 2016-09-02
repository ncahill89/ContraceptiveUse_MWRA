PlotSourceSigmas <- function (# Plot posteriors related to covariance matrices of (trad, modern)/total
  ### Plot posteriors of sigmas and rhos for (trad, modern)/total by source
  mcmc.array, 
  col_pxmr = c("red","purple", "green", "brown"),
  percentiles = c(0.025, 0.5, 0.975),
  return.res = FALSE
) {  
  I <- dim(mcmc.array)[1]*dim(mcmc.array)[2]
  # s refers to cov matrix for s = DHS, MICS, NS, Other
  psigma1.si <- matrix(NA,4, I)
  psigma2.si <- matrix(NA,4, I)
  prho.si <- matrix(NA,4, I)
  for (s in 1:4){ 
    parname <- paste0("T1.source.s[", s, "]")
    T11 <- mcmc.array[,,parname]
    parname <- paste0("T2.source.s[", s, "]")
    T22 <- mcmc.array[,,parname]
    parname <- paste0("T12.source.s[", s, "]")
    T12 <- mcmc.array[,,parname]
    for (i in 1:I){
      Sigma <- solve(cbind(c(T11[i], T12[i]), c(T12[i], T22[i])))
      psigma1.si[s,i] <- sqrt(Sigma[1,1])
      psigma2.si[s,i] <- sqrt(Sigma[2,2])
      prho.si[s,i] <- Sigma[1,2]/(psigma1.si[s,i]*psigma1.si[s,i])
    }
  }
  psigma1.sq <-  t(apply(psigma1.si, 1, quantile, percentiles))
  psigma2.sq <-  t(apply(psigma2.si, 1, quantile, percentiles))
  prho.sq <-  t(apply(prho.si, 1, quantile, percentiles))
  sourcenames  = c("DHS", "MICS", "National", "Other")
  rownames(psigma1.sq) =  rownames(psigma2.sq)=  rownames(prho.sq)= sourcenames
  par( mar = c(8,5,1,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  plot(psigma1.sq[, 2] ~ seq(0.9,3.9), ylim = c(0, max(psigma1.sq, psigma2.sq)*1.5),
       xlim = c(0.5,4.5),
       ylab = "sigma",pch = 19, lwd =8,
       xlab = "", xaxt = "n", col = col_pxmr[2])
  axis(1, las = 3, at = seq(1,4), labels = sourcenames)
  segments(seq(0.9,3.9), psigma1.sq[,1], seq(0.9,3.9), psigma1.sq[,3], col = col_pxmr[2], lwd = 3)
  points(psigma2.sq[,2]~seq(1.1,4.1), pch = 19 ,lwd =8, col = col_pxmr[3])
  segments(seq(1.1,4.1), psigma2.sq[,1], seq(1.1,4.1), psigma2.sq[,3],lwd =3, col = col_pxmr[3])
  legend("topright", legend = c("Traditional", "Modern"), col = col_pxmr[c(2,3)], pch = 19, lwd =3, cex = 1.2)
  
  plot(prho.sq[, 2] ~ seq(0.9,3.9), ylim = c(-1,1),
       xlim = c(0.5,4.5),
       ylab = "rho",pch = 19, lwd =8,
       xlab = "", xaxt = "n", col = 1)
  abline(h=0, lwd = 1)
  axis(1, las = 3, at = seq(1,4), labels = c("DHS", "MICS", "National", "Other")) # should correspond to source.name in GetBugsData
  segments(seq(0.9,3.9), prho.sq[,1], seq(0.9,3.9), prho.sq[,3], col = 1, lwd = 3)
  if (return.res) return(list(psigma1.sq = psigma1.sq, 
                        psigma2.sq = psigma2.sq,
                              prho.sq = prho.sq))
}
#-------------------------------------------------------------------------------------------------
PlotSourceSigmasUnmet <- function (# Plot posteriors of unmet/none by source
  ### Plot posteriors of unmet/none by source
  mcmc.array, percentiles = c(0.025, 0.5, 0.975),
  return.res  = FALSE) {
  parname <- "sigma.unmet.dhs"
  res.dhs.q <- quantile(c(mcmc.array[,,parname]), percentiles)
  parname <- "sigma.unmet.other"
  res.other.q <- quantile(c(mcmc.array[,,parname]), percentiles)
  par(mar = c(8,5,1,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  plot(1, res.dhs.q[2], ylim = c(0, max(res.dhs.q, res.other.q)*1.5),
       xlim = c(0.5,2.5),
       ylab = "sigma", pch = 19, lwd =8,
       xlab = "", xaxt = "n")
  axis(1, las = 3, at = seq(1,2), labels = c("DHS", "Other"))# should correspond to source.name.unmet in GetBugsData
  points(2, res.other.q[2], lwd = 8)
  segments(1, res.dhs.q[1], 1, res.dhs.q[3], lwd = 3, col = 1)
  segments(2, res.other.q[1], 2, res.other.q[3], lwd = 3, col = 1)
  legend("topright", legend = c("Unmet"), col = 1, pch = 19, lwd =3, cex = 1.5)
  if (return.res) return(list(res.other.q = res.other.q, res.dhs.q = res.dhs.q))
}
#---------------------------------------------------------------------------------
PlotBiases <- function(# Plot posteriors of bias parameters, and output the CIs
  ### Plot posteriors of bias parameters
  mcmc.array, percentiles = c(0.025,0.5, 0.975),
  return.res = FALSE){  
  biases.i3 <- NULL  
  # change Oct 20, 2012: change order of biases to correspond to order in text
  parnames.biases <- c("v.mneg", "v.mpos", "v.folk", "v.mics")
  parnames.biases.nice <- c("Sterilization excluded", "Sterilization included",
                            "Folk methods included", "Absence of probing questions")
  for (parname in (parnames.biases)){
    biases.i3 <- rbind(biases.i3, quantile(mcmc.array[,,parname],percentiles))
  }              
  rownames(biases.i3) <- parnames.biases.nice
  nbiases <- dim(biases.i3)[1]
  add <- 0
  par(mar = c(5,12.5,1,1))
  plot(seq(1-add, nbiases) ~ biases.i3[,2], yaxt = "n", ylab = "", 
       xlab = "Bias (proportion misclassified)", ylim = c(nbiases+1, 0),
       xlim = c(min(biases.i3), max(biases.i3)), type = "p",  pch = 20)
  axis(2, at = seq(1, nbiases), labels = parnames.biases.nice, las = 2, cex = 0.5)
  segments(biases.i3[,1], seq(1-add, nbiases), biases.i3[,3], seq(1-add, nbiases))
  abline(v = 0, col = 2)
  if (return.res) return(biases.i3)
}

#---------------------------------------------------------------------------------
PlotGeoEtc <- function(# Plot posteriors of perturbation multipliers
  ### Plot posteriors of perturbation multipliers
  par.V, 
  mcmc.array, 
  percentiles = c(0.025,0.5, 0.975)
){
  if (is.null(par.V$parnames.V.in.bugs)) { # change JR, 20131105
    return(invisible())
  }
  dummies.i3 <- NULL
  for (parname in par.V$parnames.V.in.bugs){
    dummies.i3 <- rbind(dummies.i3, quantile(mcmc.array[,,parname], percentiles))
  }
  ndummies.tot <- dim(dummies.i3)[1]
  add <- 0
  # max 50 dummies in each plot
  nplots <- ceiling(ndummies.tot/100)
  ndummiesperplot <- ceiling(ndummies.tot/nplots)
  # make one big plot...
  nf <- layout(t(seq(1, nplots)),
    widths = rep(8, nplots), heights = 15)
  #layout.show(nf)
  par(mar = c(2,5,2,1))
  for (j in 1:nplots){
    ndummies <- min( ndummiesperplot, ndummies.tot - (j-1)* ndummiesperplot)
    seq.select <- (j-1)* ndummiesperplot + seq(1, ndummies) 
    
    plot(seq(1-add, ndummies)~ dummies.i3[seq.select,2], yaxt = "n", ylab = "", 
         xlab = "Ratio multiplier (# obs)", cex.axis = 0.8, 
         xlim = c(0,  3), #max(dummies.i3)), 
         type = "n",
    #, #min(dummies.i3),
         # swap y-axis
         ylim = c(ndummies-1, 1.5), pch = 20)
    axis(3, at = seq(0,3,0.5), cex.axis = 0.8)
    axis(2, at = seq(1, ndummies), 
         labels = par.V$parnames.V.nice[seq.select],
         las = 2, cex = 0.001, cex.lab = 0.1, cex.axis = 0.4)
    abline(v = seq(1, ndummies), col = "grey")
    abline(v = seq(0,3,0.5), col = "grey", lwd = 2)
    abline(v = 1, col = 1, lwd= 2)
    segments(dummies.i3[seq.select,3], seq(1-add, ndummies), 
             dummies.i3[seq.select,1], seq(1-add, ndummies), 
             lwd = 1, col = 2)
    points(seq(1-add, ndummies)~ dummies.i3[seq.select,2], col = 2, pch = 20, lwd = 2)
  }
}
#---------------------------------------------------------------------------------------------
SummarizeBiases <- function( # Document number of bias parameters
  ### Documentation function to summarize how often biases were included
  winbugs.data ##<< Object of class \code{\link{winbugs.data}}
){
  ##value<<
  res <- list(nbias.folk = sum(winbugs.data$folk.ind1.j), ##<< folk
              nbias.MICS = sum(winbugs.data$source.MICS.ind1.j), ##<< MICS 
              nbias.modernneg = sum(winbugs.data$mneg.ind1.j), ##<< Modern[-]
              nbias.modernpos = sum(winbugs.data$mpos.ind1.j) ##<< Modern [+]
  )
  print(paste(names(res), res))
  return(res)
}
#----------------------------------------------------------------------------------------------
SummarizeMultipliers <- function(#Summarize multipliers
  ### Documentation function to summarize how many pertubation parameters were included
  winbugs.data,##<< Object of class winbugs.data
  data ##<< Object of class data
  ){
  multpliers.ncategories.inclbaseline <- 
    winbugs.data[c("ncat.age", "ncat.geo", "ncat.posbias",
                   "ncat.sa", "ncat.emal", "ncat.hw",
                   "ncat.posage", "ncat.negage")]  #data.frame(par.V$parnames.V.in.bugs,  par.V$parnames.V.nice)
  multpliers.ncategories.exclbaseline <-lapply(multpliers.ncategories.inclbaseline,"-", 1)
  print(paste(names(multpliers.ncategories.exclbaseline), multpliers.ncategories.exclbaseline))
  ##value<< List with number of multipliers used in each pertubation group
  return(multpliers.ncategories.exclbaseline)
}
#----------------------------------------------------------------------
# The End!
