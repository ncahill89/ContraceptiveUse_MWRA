#----------------------------------------
# F_validation 
# Leontine Alkema
#----------------------------------------

# note: PlotValidationResults (at end) is the main function, that calls this function:
PlotValResults <- function(# Summarize and plot validation results
  ### Summarize and plot validation results for one type of left-out observation (e.g. unmet, total, trad or modern),
  ### Note: returns summary results
  include.j, ##<< vector with logicals, include observation or not? 
  P.j,##<< Percentiles to be summarized
  ## (element of P.jp data frame constructed in \code{\link{GetPercentilesLeftOut}})
  nameplot, ##<< Main for overview figure
  error.est.j, ##<< Errors to be summarized
  country.info, ##<< Object from class country.info
  winbugs.data, 
  data, 
  ## (element of error.est.jp data frame constructed in \code{\link{GetPercentilesLeftOut}})
  validation.list ##<< Info about validation exercise (from mcmc.meta)
  ){
  # specify subgroups of left-out observations:
  reg.j <-  country.info$reg.c[winbugs.data$getc.j]
  ssa.j <- country.info$ssa.c[winbugs.data$getc.j]
  mics.j <- ifelse(data$source.j=="MICS", T, F)
  dhs.j <- ifelse(data$source.j=="DHS", T, F)
  dev.j <- country.info$dev.c[winbugs.data$getc.j]
             
  select.ji <- cbind(rep(TRUE, length(ssa.j)), ssa.j=="No", ssa.j=="Yes", 
                     dhs.j, mics.j, !mics.j & !dhs.j, dev.j=="Poor", dev.j=="Rich")
  namesselect <- c("All", "OutsideSSA", "SSA", "DHS", "MICS", "Other", "Poor", "Rich")

  I <- dim(select.ji)[2] # number of subgroups
  table.ix <- NULL # x refers to table columns (outputs)
  table2.ix <- NULL # x refers to table columns (outputs)
  for (i in 1:I){
  #  nobs.i[i] <- sum(!is.na(Psselect[select.ji[,i]]))                   
    select <- select.ji[,i]& (include.j == TRUE)                
    table.ix <- rbind(table.ix, GetRow(P.c = P.j[select]))
    table2.ix <- rbind(table2.ix, c(mean(error.est.j[select], na.rm= T),
                                  median((error.est.j[select]), na.rm= T),
                                     median(abs(error.est.j[select]), na.rm= T),
                      mean(P.j[select] <0.5, na.rm = T) )) # if P.j is less than 0.5, the obs fell below median
  }
  nobs.i <- table.ix[, "# Obs"]
  rownames(table.ix) <- paste0(namesselect, " (", nobs.i, ")")
  rownames(table2.ix) <- paste0(namesselect, " (", nobs.i, ")")
  colnames(table2.ix) <- c("Mean error", "Median error", "Median abs error", "Prop below median")
    
  par(mfrow = c(2,4), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
  
  # 1. inside 80%
  par(mar = c(5,12,5,1))
  plot(1, type = "n", ylim = c(I, 0), main = "" , xlim = c(0.5,1), 
         xlab = "Coverage of 80% PI", ylab = "", yaxt = "n")
  # for inside 80
  for (i in 1:I){
    res1 <- qbinom(p = c(0.025, 0.975), size = nobs.i[i], prob = 0.8)/nobs.i[i]
    segments(res1[1], i, res1[2], i, lwd = 10, col = "grey")
    points(i~table.ix[i, "Within 80% CI"], col = 1, pch = 19, lwd = 3)
  }
  abline(v = 0.8)
  axis(2, las = 2, labels = rownames(table.ix), at  = seq(1, I))
  
  par(mar = c(5,5,5,1))
  # 2. outside 80% CI
  plot(1, type = "n", ylim = c(I,0), xlim = c(0,0.3), xlab = "Proportion below/above 80% PI", 
       ylab = "", yaxt = "n", main = nameplot)
  for (i in 1:I){
    res1 <- qbinom(p = c(0.025, 0.975), size = nobs.i[i], prob = 0.1)/nobs.i[i]
    segments(res1[1], i, res1[2], i, lwd = 10, col = "grey")
    points((i+0.1)~table.ix[i, "Below 80% CI"], col = 2, pch = 19, lwd = 3)
    points((i-0.1)~table.ix[i, "Above 80% CI"], col = 1, pch = 19, lwd = 3)
  }
  axis(1, las = 2, labels = NULL, at  = i)
  abline(v = 0.1)
  legend("bottomright", legend = c("Above", "Below"), col = c(1,2), lwd = 3, pch = 19, lty = -1, bty = "n")
  
  # 3. inside 95% CI
  plot(1, type = "n", ylim = c(I,0),
        main = "", xlim = c(0.7,1), 
         xlab = "Coverage of 95% PI", ylab = "", yaxt = "n")
  for (i in 1:I){
    res1 <- qbinom(p = c(0.025, 0.975), size = nobs.i[i], prob = 0.95)/nobs.i[i]
    segments(res1[1], i, res1[2], i, lwd = 10, col = "grey")
     points(i~table.ix[i, "Within 95% CI"], col = 1, pch = 19, lwd = 3)
  }
  abline(v = 0.95)
  
  # 4. outside 95%
  #par(mar = c(5,5,5,12))
  plot(1, type = "n", ylim = c(I,0), xlim = c(0,0.3), xlab = "Proportion below/above 95% PI", 
       ylab = "", yaxt = "n")
  for (i in 1:I){
    res1 <- qbinom(p = c(0.025, 0.975), size = nobs.i[i], prob = 0.025)/nobs.i[i]
    segments(res1[1], i, res1[2], i, lwd = 10, col = "grey")
    points((i+0.1)~table.ix[i, "Below 95% CI"],col = 2, pch = 19, lwd = 3)
    points((i-0.1)~table.ix[i, "Above 95% CI"],col = 1, pch = 19, lwd = 3)
  }
  axis(1, las = 2, labels = NULL, at  = i)
  abline(v = 0.025)
  legend("bottomright", legend = c("Above", "Below"), col = c(1,2), lwd = 3, 
         pch = 19, lty = -1, bty = "n")
  
  hist(P.j[include.j], freq = FALSE, main = "Percentiles", xlab = "P")
  abline(h=1, col = 2)  
  hist(P.j[ssa.j=="Yes"&include.j], freq = FALSE, breaks = 10, main = "Percentiles SSA", xlab = "P")
  abline(h=1, col = 2)  
  hist(error.est.j[include.j], freq = T, breaks = 20, main = "Errors", xlab = "Error")
  abline(v=0, col = 1)
  abline(v = mean(error.est.j[include.j], na.rm = T), col = 3)
  abline(v = median(error.est.j[include.j], na.rm = T), col = 4)
#  hist(st.error.j, freq = FALSE, main = "St. Errors", xlab = "St. Error")
#  abline(v=0, col = 2)
  hist(error.est.j[ssa.j=="Yes" & include.j], freq = T,  breaks = 20, main = "Errors SSA", xlab = "St. Error")
  abline(v=0, col = 1)
  abline(v = mean(error.est.j[ssa.j=="Yes" & include.j], na.rm = T), col = 3)
  abline(v = median(error.est.j[ssa.j=="Yes" & include.j], na.rm = T), col = 4)
  ##value<<
  return(list(table.ix = table.ix,##<< Summary table for P.j
              table2.ix = table2.ix ##<< Summary table for error.est.j
              ))
}

#----------------------------------------------------------------------------------
GetPercentilesLeftOut <- function(# Calculate where left-out obs falls in pred. posterior distribution
  ### Calculate where left-out obs falls in pred. posterior distribution:
  ### P.j  =(approx) P(Ynew <= yleftout) = percentile of left-out obs j in posterior predictive distribution
  ### and errors:
  ### error.est.j  = y.j - median(Ynew.j).

  data, mcmc.array, winbugs.data){
    Punmet.j  <- error.unmet.j <- error.est.unmet.j <- 
    Ptot.j <- Ptrad.j <- Pmodern.j <- 
        error.est.trad.j <-
        error.est.modern.j <- 
        error.est.tot.j <- rep(NA, winbugs.data$J)
  # if trad, mod and unmet were left out:
  if (!is.null(winbugs.data$getj.test.k)){
    for (j in winbugs.data$getj.test.k[1:winbugs.data$n.test.breakdown]){
      tradobs <- data$props.trad.j[j]
      modobs <- data$props.modern.j[j]
      totCPobs <-  tradobs + modobs
      parname.tr <- paste0("pred.logratio.ytrad.j[", j, "]")
      parname.mod <- paste0("pred.logratio.ymodern.j[", j, "]")
      a <- exp(c(mcmc.array[, , parname.tr])) + exp(c(mcmc.array[, , parname.mod]))
      totCP <- a/(a+1)
      trad <- (1-totCP)*exp(c(mcmc.array[, , parname.tr]))
      mod <- (1-totCP)*exp(c(mcmc.array[, , parname.mod]))
      Ptrad.j[j] <- mean(trad <= tradobs)
      Pmodern.j[j] <- mean(mod <= modobs)
      Ptot.j[j] <- mean(totCP <= totCPobs)
      error.est.trad.j[j] <- tradobs - median(trad)
      error.est.modern.j[j] <- modobs - median(mod)
      error.est.tot.j[j] <- totCPobs - median(totCP)
#      # from mcmc sample directly, using median of mean (expected to be slightly different)
#      parname.tr <- paste0("q.trad.j[", j, "]")
#      parname.mod <- paste0("q.modern.j[", j, "]")
#      error.trad.j[j] <- tradobs - median(c(mcmc.array[, , parname.tr]))
#      error.modern.j[j] <- modobs - median(c(mcmc.array[, , parname.mod]))
#      error.tot.j[j] <- totCPobs - median(c(mcmc.array[, , parname.tr])+c(mcmc.array[, , parname.mod]))
 
      if (is.element(j, winbugs.data$getj.test.unmet.k)){
          unmetobs <- data$props.unmet.j[j]
          parname <- paste0("pred.logitratio.yunmet.j[", j, "]")
          # use totCP to get sampled unmet 
          unmet <- (1-totCP)*invlogit(c(mcmc.array[, , parname]))
          Punmet.j[j] <- mean(unmet <= unmetobs) 
          error.est.unmet.j[j] <- unmetobs - median(unmet)
#          parname <- paste0("q.unmet.j[", j, "]")
#          error.unmet.j[j] <- unmetobs - median(c(mcmc.array[, , parname]))
      }
    } #end j-loop
  } else { # unmet only
      for (j in winbugs.data$getj.test.unmet.k){
          tradobs <- data$props.trad.j[j]
          modobs <- data$props.modern.j[j]
          totCPobs <-  tradobs + modobs
          parname.tr <- paste0("pred.logratio.ytrad.j[", j, "]")
          parname.mod <- paste0("pred.logratio.ymodern.j[", j, "]")
          a <- exp(c(mcmc.array[, , parname.tr])) + exp(c(mcmc.array[, , parname.mod]))
          totCP <- a/(a+1)

          # copied from above
          unmetobs <- data$props.unmet.j[j]
          parname <- paste0("pred.logitratio.yunmet.j[", j, "]")
          # use totCP to get sampled unmet 
          unmet <- (1-totCP)*invlogit(c(mcmc.array[, , parname]))
          Punmet.j[j] <- mean(unmet <= unmetobs)
          error.est.unmet.j[j] <- unmetobs - median(unmet)
#         parname <- paste0("q.unmet.j[", j, "]")
#          error.unmet.j[j] <- unmetobs - median(c(mcmc.array[, , parname]))
      } # end j-loop
  }# end else
  ##value<< List with data.frames P.jp and error.est.jp 
    ## with percentiles and errors respectively as explained in the description, 
    ## where j refers to the observation index, and p 
    ## refers to trad, modern, tot, unmet. 
  return(list(P.jp =   # ##<<P.jp is a data frame with 
    # ##describe<<
               data.frame(Ptrad.j, ##<<,
                      Pmodern.j, ##<<,
                      Ptot.j,##<<,
                      Punmet.j##<<,
                           ),
    # ##end<<
                           
#           error.jp = data.frame(error.trad.j,  ##<<,
#                        error.modern.j,  ##<<,
#                                 error.tot.j, ##<<, 
#                                 error.unmet.j ##<<,
#                        ), 
                          
              
          error.est.jp = data.frame(error.est.trad.j, 
                                   error.est.modern.j, 
                                    error.est.tot.j, 
                                    error.est.unmet.j
                                    )
              ))
}

#------------------------------------------------------------------------------
GetRow <- function(# Summarize proportion of left-out observations outside various PIs
  ###  Summarize proportion of left-out observations outside various PIs
  P.c ##<< Vector with percentiles (see \code{\link{GetPercentilesLeftOut}}) 
  ){
  out.in.CI95 <- round(OutsideAndInBounds(P.c = P.c),2)
  out.in.CI80 <- round(OutsideAndInBounds(P.c = P.c, CIlow = 0.1, CIup = 0.9),2)
  n <- sum(!is.na(P.c))
  res <- data.frame(n, out.in.CI80, out.in.CI95) 
  colnames(res) <- c("# Obs", 
        "Below 80% CI", "Within 80% CI", "Above 80% CI",
        "Below 95% CI", "Within 95% CI", "Above 95% CI")    
  ##value<< Data frame with one row and columns 
  ##c("# Obs",  "Below 80% CI", "Within 80% CI", "Above 80% CI",
  ##      "Below 95% CI", "Within 95% CI", "Above 95% CI")
  return(res)
}

#------------------------------------------------------------------------------
OutsideAndInBounds <- function(# Find proportion of left-out observations outside a specific PI
  ### Find proportion of left-out observations outside a specific PI
  P.c, ##<< Vector with percentiles (see \code{\link{GetPercentilesLeftOut}})
  CIlow = 0.025, ##<< Lower bound PI
  CIup = 0.975 ##<< Upper bound PI
  ){
  res <- c(mean(P.c < CIlow, na.rm = T),
          mean((P.c >= CIlow) & (P.c <= CIup), na.rm = T),
          mean(P.c > CIup, na.rm = T))
  ##value<<
  ## ("Below CI","In CI","Above CI")
  return(data.frame("Below CI" = res[1], "In CI" = res[2], "Above CI" = res[3]))
}

#----------------------------------------------------------------------------------
GetTraining <- function (# Construct training sets for validation exercise
  ### Construct training sets for validation exercise
  ### Choose one option from (at.random, at.end, exclude.unmet.only)
  data, winbugs.data = NULL, ##<< winbugs.data needs to be provided when leaving out obs at the end
  at.random = T, ##<< Logical: leave out 20% at random?
  at.end = FALSE, ##<< Logical: leave out observations after and IN \code{year.cutoff}?
  year.cutoff = 2005,##<< Used only when \code{at.end} is TRUE
  exclude.unmet.only = FALSE, ##<< Logical: leave out unmet data only?
  leave1 = FALSE, ##<< Not yet implemented
  seed = 12345##<< Set seed when sampling observation to leave out
  ) {
  
  # note: make sure training/test sets are not empty for breakdown/unmet
  # and that the training set is not empty for total
  set.seed(seed)
  if (at.random){
    # 20% at random (from break-down into modern/trad only)
    # (ok if some countries end up without observations in training, they're still in test)
    M <- sum(!is.na(data$props.modern.j))
    n.training.breakdown <- max(1, min(M-1, round(0.8*M)))
    # make sure unmet test or traning set is not empty (else loop does not work)
    unmet.empty <- T
    while (unmet.empty){
      getj.training.k <- sample(seq(1,M), n.training.breakdown)
      if (sum(!is.na(data$props.unmet.j[getj.training.k])) !=0 & sum(!is.na(data$props.unmet.j[-getj.training.k])) !=0){
        unmet.empty <- F
      }
    }
  }
  if (at.end){
    # leave out obs after AND IN year.cutoff
    getj.training.k <- seq(1, length(data$years.j))[ (data$years.j <= year.cutoff)]
    # make sure there's one obs in test/training for unmet, and one in training for total
    # (note: test set is not constructed for total)
    if( sum(!is.na(data$props.unmet.j[getj.training.k])) == 0){
          print("Warning: training set unmet is empty")
    }
    if( sum(!is.na(data$props.unmet.j[-getj.training.k])) == 0){
          print("Warning: test set unmet is empty")
    }
    if( sum(!is.na(data$props.tot.j[getj.training.k])) == 0){
          print("Warning: training set tot is empty")
    }
  }
  if (exclude.unmet.only){
    ##details<<  If \code{exclude.unmet.only}, leave out all data in round(20%  no countries with data)  
    iso.c <- unique(data$iso.j)
    C <- length(iso.c)
    n.unmet.c <- rep(NA, C)
    for (c in 1:C){
      n.unmet.c[c] <- sum(!is.na(data$props.unmet.j[data$iso.j == iso.c[c]]))
    }
    ncountrieswithdata <- sum(n.unmet.c !=0) # 108 countries
    nleaveout <- round(0.2*ncountrieswithdata)
    # sample countries to leave out
    indicesofcountriestoleaveout <-  sample(seq(1,C)[n.unmet.c>0], nleaveout)
    getj.training.k <- seq(1, length(data$years.j))[ is.element(data$iso.j, iso.c[-indicesofcountriestoleaveout])]
  }
  if (leave1) {
    # leave out all observations but one in set of countries with > 1 obs.
    # not yet implemented!
  }
  # save(getj.training.k, file = "getj.training.k.rda")
  # Note: reset seed after calling this function? 
  ##value<< vector with indices for training set (getj.training.k)   
  return(getj.training.k)
}

#----------------------------------------------------------
PlotValidationResults <- function(# Plot lots of results!
  ### Wrapper function to plot lots of results.
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  fig.dir = NULL, ##<< Directory to store overview plots. If NULL, folder "fig" in current working directory.
  #  plot.ind.country.results = FALSE ##<< Create zillion plots for all countries?
  ## If TRUE, plots are saved in subdirectory "country.plots" in fig.dir.
  plot.CIs  = TRUE, # create contr use overview plots for all countries?
  return.res.paper = FALSE
  ){
  
  if (is.null(fig.dir)){
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  # put separate?
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  
  load(file = file.path(output.dir,"res.country.rda")) # change JR, 20140418
  load(file = file.path(output.dir,"mcmc.meta.rda")) # change JR, 20140418
  validation <- !is.null(mcmc.meta$validation.list)
  validation.list <- mcmc.meta$validation.list
  
  data <- mcmc.meta$data.raw$data
  country.info <- mcmc.meta$data.raw$country.info
  region.info <- mcmc.meta$data.raw$region.info
  winbugs.data <- mcmc.meta$winbugs.data
  
  if (!validation){
    print("Not a validation exercise!")
    return()
  }
  load(file = file.path(output.dir,"Ps_validation.rda")) # change JR, 20140418
  if (validation.list$exclude.unmet.only){
    select.unmet.c <- unique(winbugs.data$getc.j[winbugs.data$getj.test.unmet.k])
  } else {
    select.unmet.c <- NULL
  }
  if (plot.CIs){
    PlotDataAndEstimates(data.raw = mcmc.meta$data.raw, #country.info = country.info, 
                         validation = validation, 
                         select.c = select.unmet.c,
                         getj.test.k = mcmc.meta$winbugs.data$getj.test.k,  
                         getj.test.unmet.k = mcmc.meta$winbugs.data$getj.test.unmet.k,
                         CI.Lg.Lcat.qt = res.country$CIprop.Lg.Lcat.qt,
                         CIstar.Lg.Lcat.qt = res.country$CIstar.Lg.Lcat.qt,
                         CIratio.Lg.Lcat.qt = res.country$CIratio.Lg.Lcat.qt,
                         fig.name = file.path(fig.dir, paste0(run.name, "CIs.pdf")) # change JR, 20140418
    )
  }
  # open text file to write results
  validation.res.file <- file.path(fig.dir, paste0(run.name, "validationres.html")) # change JR, 20140418
  cat("", file = validation.res.file, append = F)
  
  # when summarizing results in table: select only last observation year in each country
  # Note: can be different for unmet because less observations!
  j.include.c <- j.include.unmet.c <- rep(0, winbugs.data$C)
  J <- length(data$props.unmet.j)
  select.unmet.c <- unique(winbugs.data$getc.j[winbugs.data$getj.test.unmet.k])
  for (c in select.unmet.c){
    select <- seq(1, J)[winbugs.data$getc.j==c 
                        & is.element(seq(1, J), winbugs.data$getj.test.unmet.k)
                        & !is.na(data$props.unmet.j) ]
    kk <- which.max(data$years.j[select])
    j.include.unmet.c[c] <- seq(1, J)[select][kk]
  }
  # T/F indicator
  include.unmet.j <- is.element(seq(1, length(data$years.j)), j.include.unmet.c)

#   # which country has P.j missing?
#   length(select.unmet.c)
#   sum(include.unmet.j)
#   sum(!is.na(Ps$P.jp[include.unmet.j,3]))
#   sum(!is.na(winbugs.data$logitratio.yunmet.j[include.unmet.j]))
#   
#   sum(is.na(Ps$P.jp[include.unmet.j,3]))
#   which.max(is.na(Ps$P.jp[include.unmet.j,3]))
#   winbugs.data$getc.j[include.unmet.j]
#   data$name.j[include.unmet.j][40]
#  length(include.unmet.j)
  
  
  # now make the one for trad/modern
  if (!is.null(winbugs.data$getj.test.k)){
    select.c <- unique(winbugs.data$getc.j[winbugs.data$getj.test.k[1:winbugs.data$n.test.breakdown]])
    for (c in select.c){
      select <- seq(1, J)[winbugs.data$getc.j==c 
                          & is.element(seq(1, J), winbugs.data$getj.test.k[1:winbugs.data$n.test.breakdown])
                          & !is.na(data$props.trad.j)]
      # most recent when left out at end
      # random for 20%
      if (length(select)==1){
        kk <- select # kk different here from above!
      } else {
        kk <- ifelse(validation.list$at.end, select[which.max(data$years.j[select])],
                     #                  which.min(data$years.j[select]) )
                     #                  which.max(data$years.j[select]) )
                     sample(select, size = 1) )
        # watch out: sample (x,1) gives a sample from 1:x!
      }
      j.include.c[c] <- seq(1, J)[kk]
    }
    include.j <- is.element(seq(1, length(data$years.j)), j.include.c)  
  } else {
    include.j <- rep(FALSE, length(data$years.j))
  }
  # check:
  #data.frame(data$iso.j, data$years.j, data$props.modern.j, include.j)
  #data.frame(data$iso.j, data$years.j, data$props.unmet.j, include.unmet.j)
  
  nameplot2 <-ifelse(mcmc.meta$validation.list$exclude.unmet.only, "Unmet only",
                     ifelse(mcmc.meta$validation.list$at.random, "20%",
                            ifelse(mcmc.meta$validation.list$at.end,  paste("Up to (excl) year", 
                                                                            mcmc.meta$validation.list$year.cutoff),NA)))
  res.paper <- NULL
  pdf(file.path(fig.dir, paste0(run.name, "propoutside.pdf")), width = 16, height = 10) # change JR, 20140418
  for (p in 1:4){
    if (p==4){
      includeforp.j <- include.unmet.j
    } else {
      includeforp.j <- include.j 
    }
    P.j <- Ps$P.jp[,p]
    #error.j <- Ps$error.jp[,p]
    error.est.j <- Ps$error.est.jp[,p]
    nameplot1s <- c("Traditional", "Modern", "Total", "Unmet")
    if (sum(!is.na(P.j))>0){
      nameplot1 <- nameplot1s[p]
      res <- PlotValResults(include.j = includeforp.j, P.j = P.j, 
                            data = data, 
                            country.info = country.info,
                            winbugs.data = winbugs.data, 
                            nameplot = paste(nameplot1, nameplot2, sep = " - "),
                            error.est.j = error.est.j)
      # print table: html and tex
      print(xtable(res$table.ix, #digits = c(0,0,0,0,0), 
                   caption = paste(nameplot1, nameplot2)), 
            file = validation.res.file, append = T, type = "html")
      print(xtable(res$table2.ix, #digits = c(0,0,0,0,0), 
                   caption = paste(nameplot1, nameplot2)), 
            type="html", file = validation.res.file, append = T)
      addpaper <-unlist( c(
        res$table.ix[1, c("# Obs")],
        100*res$table.ix[1, c("Below 95% CI","Within 95% CI","Above 95% CI")] , # first one is all obs      
        100*unlist(res$table2.ix[1, c("Prop below median")]),
                           res$table2.ix[1, c("Median error","Median abs error")]
                           ))
      if (is.null(res.paper)){
        res.paper <- addpaper
      } else {
        res.paper <- rbind(res.paper, addpaper)
      }
    }
  }
  dev.off()
  if (prod(dim(res.paper))>1) rownames(res.paper) <- nameplot1s
  # # plot(exp(Ps$error.j$error.trad.j) ~ data$props.trad.j)
  # # plot(exp(Ps$error.j$error.trad.j[ssa.j=="Yes"]) ~ data$props.trad.j[ssa.j=="Yes"])
  # # plot(Ps$error.j$error.trad.j ~ data$props.trad.j)
  # # plot(Ps$error.j$error.trad.j ~ winbugs.data$ratios.trad.modern.jn[,1])
  # # plot(Ps$error.j$error.modern.j ~ winbugs.data$ratios.trad.modern.jn[,2])
  # # plot(Ps$error.j$error.unmet.j ~ winbugs.data$logitratio.yunmet.j)
  # 
  # #---------------------------------------------------------------------------------
  # # plot two sets of estimates for comparison
  # load(file = file.path(work.dir.save, paste0(name.dir.all,"CIs.rda")))
  # CIs.all <- CIs
  # load(file.path(work.dir.save, paste0(name.dir, "CIs.rda")))
  # CIs.tr <- CIs
  # 
  # pdf(file.path(work.dir.save, paste0(name.dir, "_comparison_", name.dir.all, "CIs.pdf")), width = 21, height = 12)
  #   PlotDataAndEstimates(data = data, country.info = country.info, 
  #                     # start.year = start.year, end.year = end.year, 
  #                            # par.ciq = par.ciq, plot.blue.line = FALSE,
  #                             validation = validation, 
  #                             getj.test.k = winbugs.data$getj.test.k,  
  #                             getj.test.unmet.k = winbugs.data$getj.test.unmet.k,
  #                             CIs = CIs.all, CIs2 = CIs.tr, 
  #                             name.dir = name.dir.all, name.dir2 = name.dir,
  #                             overview.country.plot = F, country.plot = F)
  # dev.off()
  
  #---------------------------------------------------------------------------------------------
  # Compare CIs in comp.year (e.g. 2008)
  # when leaving out data at the end for countries where data was left-out
  # or for unmet
  # a. get countries where data was left out
  # Note: for unmet we use the same countries as where trad/mod were left out
  # if (!is.null(winbugs.data$getj.test.k)){
  #   select.c <- unique(winbugs.data$getc.j[winbugs.data$getj.test.k[1:winbugs.data$n.test.breakdown]])
  # } else {
  #   select.c <- unique(winbugs.data$getc.j[winbugs.data$getj.test.unmet.k])
  # }
  # 
  # comp.year <- 2008.5
  # # b. results
  # # diff = all - training
  # # all (updated) CI outside old (training) 80% CI?
  # diff.cp <- above80.cp <- below80.cp <- matrix(NA, length(country.info$iso.c),4) # trad, mod, tot, unmet 
  # for (c in select.c){
  #   for (p in 1:4){
  #     CIall.median <- ifelse(p==1, CIs.all$CIs.trad.cqt[c,3,CIs.all$est.years==comp.year],
  #                     ifelse(p==2, CIs.all$CIs.modern.cqt[c,3,CIs.all$est.years==comp.year],
  #                     ifelse(p==3, CIs.all$CIs.tot.cqt[c,3,CIs.all$est.years==comp.year],
  #                       CIs.all$CIs.unmet.cqt[c,3,CIs.all$est.years==comp.year])))
  #     CItr.q <- ifelse(rep(p,5)==1, CIs.tr$CIs.trad.cqt[c,,CIs.tr$est.years==comp.year],
  #                     ifelse(rep(p,5)==2, CIs.tr$CIs.modern.cqt[c,,CIs.tr$est.years==comp.year],
  #                     ifelse(rep(p,5)==3, CIs.tr$CIs.tot.cqt[c,,CIs.tr$est.years==comp.year],
  #                       CIs.tr$CIs.unmet.cqt[c,,CIs.tr$est.years==comp.year])))
  #     diff.cp[c,p] <- CIall.median - CItr.q[3]
  #     below80.cp[c,p] <- ifelse ((CIall.median < CItr.q[2]), 1, 0)
  #     above80.cp[c,p] <- ifelse ((CIall.median > CItr.q[4]), 1, 0)
  # }}
  # 
  # # specify subgroups of left-out observations:
  # reg.c <-  country.info$reg.c
  # ssa.c <- country.info$ssa.c
  # dev.c <- country.info$dev.c
  # select.ci <- cbind(rep(TRUE, winbugs.data$C), ssa.c=="No", ssa.c=="Yes", 
  #                    dev.c=="Poor", dev.c=="Rich")
  # namesselect <- c("All", "OutsideSSA", "SSA", "Poor", "Rich")
  # I <- dim(select.ci)[2] # number of subgroups
  # # already defined
  # # nameplot2 <-ifelse(mcmc.meta$validation.list$exclude.unmet.only, "Unmet only",
  # #                    ifelse(mcmc.meta$validation.list$at.random, "20%",
  # #                     ifelse(mcmc.meta$validation.list$at.end,  paste("Up to (excl) year",
  # #                     mcmc.meta$validation.list$year.cutoff),NA)))
  # for (p in 1:4){
  #   table.ix <- NULL # x refers to table columns (outputs)
  #   nameplot1 <- c("Traditional", "Modern", "Total", "Unmet")[p]
  #   for (i in 1:I){
  #     select <- select.ci[,i]
  #     table.ix <- rbind(table.ix, 
  #                     c(sum(!is.na(diff.cp[select,p])), 
  #                       mean(diff.cp[select,p], na.rm = T), mean(abs(diff.cp[select,p]), na.rm = T),
  #                       mean(below80.cp[select,p], na.rm = T), mean(above80.cp[select,p], na.rm = T)))
  #   }
  #   nobs.i <- table.ix[, 1]
  #   rownames(table.ix) <- paste0(namesselect, " (", nobs.i, ")")
  #   colnames(table.ix) <- c("# Obs", "Mean diff", "Medan abs diff", "% Below 80CI", "% Above 80CI")
  #   # print table
  #   print(xtable(table.ix, #digits = c(0,0,0,0,0), 
  #                 caption = paste("Change in median/CI: Results for ", comp.year, ", ", nameplot1, " (", 
  #                                 nameplot2, ")")), 
  #                 file = validation.res.file, append = T, type = "html")
  # }
  # # plots
  # # hist(diff.cp[,1])
  if (return.res.paper){
    propleftoutunmet <- (sum(!is.na(winbugs.data$logitratio.yunmet.j[winbugs.data$getj.test.unmet.k]))/
          sum(!is.na(winbugs.data$logitratio.yunmet.j)      ))
    propleftoutmod <-  (sum(!is.na(winbugs.data$ratios.trad.modern.jn[,2][winbugs.data$getj.test.k]))/
           sum(!is.na(winbugs.data$ratios.trad.modern.jn[,2])      ))
        
    return(  list(res.paper = res.paper, propleftoutunmet = propleftoutunmet,
                  propleftoutmod = propleftoutmod, 
                  nunmetleftout = length(select.unmet.c ),
                  nunmetleftout2 = sum(include.unmet.j)))
                  #nleftout = length(select.c ))) # matches up for non-unmet exe, error for unmet exercise
    # select.unmet.c relevant for unmet val exericise only
  }
}# end non-validation results

#----------------------------------------------------------------------   
# The End!
