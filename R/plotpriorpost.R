#--------------------------------------------------------------------------
# source_priorpost.R
# LA, Jan 2012
#--------------------------------------------------------------------------
PlotPriorPost <- function(# Plot posteriors and priors for hyper parameters
  ### Plot posteriors and priors for hyper parameters
  run.name = "test", ##<<,
  output.dir = NULL, ##<<, directory where main output are stored, 
  # defaults to output/run.name in current working directory
  fig.dir = NULL ##<< defaults to fig in current working directory
){
  if (is.null(output.dir))
    output.dir <- file.path(getwd(), "output", run.name, "/")
  if (is.null(fig.dir)) {
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  load(file.path(output.dir, "mcmc.meta.rda")) # change JR, 20140418
  load(file.path(output.dir, "mcmc.array.rda")) # change JR, 20140418
  
  if (!mcmc.meta$general$do.country.specific.run) { # for global run
    pdf(file.path(fig.dir, paste0(run.name, "priorpost_all.pdf")), width = 10, height = 10)
    #----------------------------------------------------------------------------------
    # prior and posteriors of kappas (variances in the BHMs for logistic curves) and two more gammas
    parnames <- c("sigma.wc", "sigma.Rwc", "sigma.lpc", "sigma.lrc", "sigma.Sc","sigma.RTc")
    rateforparnames <- unlist(mcmc.meta$winbugs.data[c(
      "halfnu0sigma2.wc0","halfnu0sigma2.Rwc0",
      "halfnu0sigma2.lpc0", "halfnu0sigma2.lrc0", 
      "halfnu0sigma2.RTc0" )])
    
    par(mfrow = c(1,1))
    if (mcmc.meta$general$change.priors.to.zerolower){ # unif priors
      # priors not implemented
      for (parname in c(parnames, "sigma.Tc", "sigma.earlierTc","sigma.unmetc")) {
        PlotPostOnly(post.samp =  c(mcmc.array[,,parname]), parname = parname)  
      }
    } else { # gamma priors
      p <- 0
      for (parname in parnames){
        p <- p+1
        PlotPostSDWithGammaPrior(post.samp =  c(mcmc.array[,,parname]), 
                                 priorshape = mcmc.meta$winbugs.data$halfnu0, 
                                 priorrate = rateforparnames[p],
                                 parname = parname)  
      }
      for (parname in "sigma.unmetc"){
        PlotPostSDWithGammaPrior(post.samp =  c(mcmc.array[,,parname]), 
                                 priorshape = 0.5, 
                                 priorrate = mcmc.meta$winbugs.data$halfsigma2.unmetc0,
                                 parname = parname)  
      }
      for (parname in "sigma.earlierTc"){
        PlotPostSDWithGammaPrior(post.samp =  c(mcmc.array[,,parname]), 
                                 priorshape =  mcmc.meta$winbugs.data$halfnu0_rich, 
                                 priorrate = mcmc.meta$winbugs.data$halfnu0_rich_sigma2.earlierTc0,
                                 parname = parname)  
      }
      for (parname in "sigma.Tc"){
        PlotPostSDWithGammaPrior(post.samp =  c(mcmc.array[,,parname]), 
                                 priorshape =  mcmc.meta$winbugs.data$halfnu0_poor, 
                                 priorrate = mcmc.meta$winbugs.data$halfnu0_poor_sigma2.Tc0,
                                 parname = parname)  
      }
    }
    #----------------------------------------------------------------------------------
    # other gammas, with shape 0.5
    parnames <- c("sigma.sourcetot")
    prrates <-  unlist(mcmc.meta$winbugs.data[c("halfsigma2.sourcetot0")])
    for (p in 1:length(parnames)){
      parname <- parnames[p]
      PlotPostSDWithGammaPrior(post.samp =  c(mcmc.array[,,parname]), 
                               priorshape = 0.5, priorrate = prrates[p],
                               parname = parname)
    }
    #--------------------------------------------------------------------------
    # wishart T.s
    # Simulate from Wishart prior
    # Note: bugs notation, thus R is on the variance scale
    # R <- cbind(c(0.1,0),c(0,0.1))
    R <- mcmc.meta$winbugs.data$R
    Vinv <- solve(R)
    k <- 3 # change JR, 20131113: from k <- 4 # change JR, 20131031: from k <- 2
    nsimu <- 1000
    sigma.1 <- rep(NA, nsimu)
    sigma.2 <- rep(NA,nsimu)
    rho <- rep(NA,nsimu)
    for (s in 1:nsimu){
      Sigma <- solve(rwish(v = k, S = Vinv)) # R specs: T ~ W(k,S) with E(T) = k*S, thus S = Vinv
      sigma.1[s] <- sqrt(Sigma[1,1])
      sigma.2[s] <- sqrt(Sigma[2,2])
      rho[s] <- Sigma[1,2]/(sigma.1[s]*sigma.2[s])
    }
    # hist(sigma.1, freq = F)
    # hist(1/sigma.1^2, freq = F)
    # curve(dgamma(x, 1/2, R[1,1]/2), add = T)
    # hist(sigma.2)
    # hist(rho)
    
    # posterior
    I <- dim(mcmc.array)[1]*dim(mcmc.array)[2]
    # s refers to cov matrix for s = DHS, MICS, NS, Other
    psigma1.si <- matrix(NA, 4, I)
    psigma2.si <- matrix(NA, 4, I)
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
    for (s in 1:4){
      par(mfrow = c(2,2))
      hist(psigma1.si[s,], freq = F,  col = "grey", main = "", xlab = "sigma11")
      lines(density(sigma.1, to =2), col= 2)
      hist(psigma2.si[s,], freq = F,  col =  "grey", xlim = c(0,0.5), main = "", xlab = "sigma22")
      lines(density(sigma.2, to =2), col= 2)
      hist(prho.si[s,], freq = F,  xlim = c(-1,1), col =  "grey", main = "", xlab = "rho")
      lines(density(rho), col= 2)
    }
    #--------------------------------------------------------------------------
    # AR parameters
    par(mfrow = c(2,2))
    for (parname in c("rho.tot", "rho.rat")){
      PlotPostWithUnifPrior(post.samp=  c(mcmc.array[,,parname]),
                            priorlow = 0, priorup = mcmc.meta$winbugs.data$rho.max, parname = parname)
    }
    for (parname in c("sigma.tot","sigma.rat")){
      PlotPostWithUnifPrior(post.samp=  c(mcmc.array[,,parname]),
                            priorlow = 0.01, priorup = mcmc.meta$winbugs.data$sigma.ar.max, parname = parname)
    }
    #--------------------------------------------------------------------------
    # hierarchical variances, non-country level
    par(mfrow = c(2,2))
    for (parname in  c("lp.world", "lr.world")){
      PlotPostWithNormalPrior(post.samp=  c(mcmc.array[,,parname]),
                              priormean = 0, priorsd = 1/sqrt(0.01), parname = parname) 
    }
    for (parname in  c("Rw.world", "w.world")){
      PlotPostWithNormalPrior(post.samp=  c(mcmc.array[,,parname]),
                              priormean = -1, priorsd = 1/sqrt(0.01), parname = parname) 
    }
    
    parnames <-c("RT.world", "T.world", "Tearlier")
    means <- unlist(mcmc.meta$winbugs.data[c("mean.RTworld", "mean.Tworld", "mean.Tearlier")])
    p <- 0
    for (parname in  parnames){
      p <- p+1
      PlotPostWithNormalPrior(post.samp=  c(mcmc.array[,,parname]),
                              priormean = means[p], priorsd = 1/sqrt(mcmc.meta$winbugs.data$tau0.T), parname = parname) 
    }
    par(mfrow = c(2,2))
    parnames <-c("sigma.RTreg", "sigma.Treg", "sigma.RTsubreg", "sigma.Tsubreg")
    for (parname in  parnames){
      PlotPostWithUnifPrior(post.samp=  c(mcmc.array[,,parname]),
                            priorlow = 0, priorup = mcmc.meta$winbugs.data$sigmaTregsubreg.upper, 
                            parname = parname) 
    }
    par(mfrow = c(2,2))
    parnames <-c("sigma.wreg", "sigma.wreg", "sigma.wsubreg", "sigma.wsubreg")
    for (parname in  parnames){
      PlotPostWithUnifPrior(post.samp=  c(mcmc.array[,,parname]),
                            priorlow = 0, priorup = mcmc.meta$winbugs.data$sigmawregsubreg.upper, 
                            parname = parname) 
    }
    #--------------------------------------------------------------------------
    # data dummies
    par(mfrow = c(2,2))
    parnames <-c("v.mics", "v.mneg", "v.folk", "v.mpos")
    for (parname in  parnames){
      PlotPostWithUnifPrior(post.samp=  c(mcmc.array[,,parname]),
                            priorlow = 0, priorup = 1, 
                            parname = parname) 
    }
    parnames <-c("sigma.geo.m[1]", "sigma.geo.m[2]")
    p <- 0
    for (parname in  parnames){
      p <- p+1
      PlotPostWithUnifPrior(post.samp = c(mcmc.array[,,parname]),
                            priorlow = 0, priorup = 2,#c(5,5,2,2)[p], 
                            parname = parname) 
    }
    parnames <-c("sigma.pos")
    p <- 0
    for (parname in  parnames){
      p <- p+1
      PlotPostWithUnifPrior(post.samp = c(mcmc.array[,,parname]),
                            priorlow = 0.01, priorup = 2,#c(5,5,2,2)[p], 
                            parname = parname) 
    }
    #parnames <-c("mu.pos.m[1]", "mu.pos.m[2]")
    par(mfrow = c(1,2))
    parnames <-c("mu.pos.m[1]", "mu.pos.m[2]")
    for (parname in  parnames){
      PlotPostWithNormalPrior(post.samp =c(mcmc.array[,,parname]),
                              priormean = -2, priorsd = 1/sqrt(0.64), parname = parname)
    }
    #--------------------------------------------------------------------------
    # unmet
    par(mfrow = c(2,2))
    parnames <- c("sigma.unmet.dhs", "sigma.unmet.other")
    for (parname in parnames){
      PlotPostWithUnifPrior(post.samp =  c(mcmc.array[,,parname]),
                            priorlow = 0.01, priorup = 2, parname = parname) 
    }
    parname <- "sigma.unmetworld"
    PlotPostWithUnifPrior(post.samp = c(mcmc.array[,,parname]),
                          priorlow = 0, priorup = 5, parname = parname)
    parname <- "c.unmet"
    PlotPostWithUnifPrior(post.samp = c(mcmc.array[,,parname]),
                          priorlow = -10, priorup = 0, parname = parname) 
    
    par(mfrow = c(2,2))
    parnames <- c("a.unmet", "b.unmet")
    sds <- 1/sqrt(unlist(mcmc.meta$winbugs.data[c("tau.a0", "tau.b0")]))
    means <- unlist(mcmc.meta$winbugs.data[c("a0.unmet", "b0.unmet")])
    p <- 0
    for (parname in parnames){
      p <- p+1
      PlotPostWithNormalPrior(post.samp=  c(mcmc.array[,,parname]),
                              priormean = means[p], priorsd = sds[p], parname = parname) 
    }
    dev.off()
  } else { # for country-specific run
    #--------------------------------------------------------------------------
    #parname <- "bias.modern"
    #if (parname %in% dimnames(mcmc.array)[[3]]) {
    #  pdf(file.path(fig.dir, paste0(run.name, "priorpost_all.pdf")), width = 10, height = 10)
    #  if (mcmc.meta$general$do.SS.run.first.pass) {
    #    PlotPostWithUnifPrior(post.samp = c(mcmc.array[,,parname]),
    #                          priorlow = 0, priorup = 30, parname = parname)
    #  } else {
    #    PlotPostWithTruncatedNormalPrior(post.samp = c(mcmc.array[,,parname]), 
    #                                     priormean = mcmc.meta$winbugs.data$median.bias.modern, 
    #                                     priorsd = mcmc.meta$winbugs.data$se.bias.modern, 
    #                                     priorlower = 0, parname = parname)
    #  }     
    #  dev.off()
    #}
  }
  #--------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------
  # Extra
  #----------------------------------------------------------------------------------
  # Joint posteriors for the hyper parameters
  # summ <- NULL
  # for (parname in c(parnames.list$parnames.h)){#, parnames.reg, parnames.c,parnames.subreg)){
  #   summ <- rbind(summ, quantile(mcmc.array[,,parname], c(0.025, 0.5, 0.975)))
  #  }
  # summ
  #parnames.list$parnames.h
  # stack.sp <- rbind(mcmc.array[,1,], mcmc.array[,2,])
  # res <- cor(stack.sp[,parnames.list$parnames.h])
  # res
  # library(arm)
  # corrplot(res, cutpts = c(0,0.25,0.5, 0.75,0.8,0.9,1))
  # #abline(v = seq(0,100), col = "grey")
  # #abline(h = seq(0,100), col = "grey")
  # 
  # parname1 <- "sigma.rat"
  # parname2 <- "rho.rat"
  # parname1 <- "mu.pos"
  # parname2 <- "sigma.pos"
  # parname1 <- "Tearlier"
  # parname2 <- "sigma.earlierTc"
  # par(mfrow = c(1,1), mar = c(5,5,1,1))
  # plot(mcmc.array[,,parname1] ~ mcmc.array[,,parname2], xlab = parname2, ylab = parname1)
  
  
  ##value<< NULL
  return(invisible(NULL))
}
#----------------------------------------------------------------------   
# The End!
