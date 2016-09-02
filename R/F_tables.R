#----------------------------------------------------------------------
# F_tables.R
#----------------------------------------------------------------------
GetPPCSymbols <- function(PPPC){# Gives vector with symbols for PPPC
  ### Gives vector witb symbols for PPPC:
  ### Probability that the 1990-2010 change was positive: *>0.9, **>0.95, ***>0.99.
  ### Probability that the 1990-2010 change was negative: '>0.9, ''>0.95, '''>0.99. 
  #PPPC <- seq(0,1,0.01)
  PPcats <- ifelse(PPPC <0.01, "\'\'\'",
                   ifelse(PPPC <0.05, "\'\'",
                          ifelse(PPPC <0.1, "\'",
                                 ifelse(PPPC >0.99, "***",
                                        ifelse(PPPC >0.95, "**",
                                               ifelse(PPPC >0.9, "*",""))))))
  return(PPcats)
}
#res <- data.frame(PPPC,PPcats)
#xtable(res)
#write.csv(res, file = "temp.csv")

#--------------------------------------------------------------------------------
GetSummaryCountries<- function (# Save a csv with summary info about countries
  ### Save a csv with summary info about countries
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  table.dir = NULL ##<< Directory to store tables. 
  ## If NULL, folder "tables" in current working directory.
  ){

  if (is.null(table.dir)){
    table.dir <- file.path(getwd(), "tables/")
    dir.create(table.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  load(file.path(output.dir, "mcmc.meta.rda"))
  load(file.path(output.dir, "par.ciq.rda"))
  
  country.info <- mcmc.meta$data.raw$country.info
  results <- data.frame(paste(mcmc.meta$data.raw$country.info$iso.c), 
                        paste(mcmc.meta$data.raw$country.info$code.c),
                        mcmc.meta$data.raw$country.info$N.c,  round(par.ciq[,,2],2))
  rownames(results) <- country.info$name.c
  ##details<< File \code{results.summ.csv} is written to \code{table.dir}, with
  ## country name, UN code, iso code, no of observations, and median estimates of the logistic parameters.
  colnames(results) <- c("num. code", "letter code", "#OBS", paste(colnames(par.ciq[,,2])))
  write.csv(results, file = file.path(table.dir, paste0(run.name, "results.summ.csv"))) # change JR, 20140418
  ##value<<NULL
  return(invisible())
}

#-------------------------------------------------------------------------------------
GetTablesRes <- function(# Save csv's with CIs for proportions and counts
  ### Save csv's with CIs for proportions and counts for countries or aggregates.
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  table.dir = NULL, ##<< Directory to store tables. 
  ## If NULL, folder "tables" in current working directory.
  res = NULL, ##<< If NULL, country/UNPD aggregates are summarized, specified with next argument.
  ## Alternatively, an object of class \code{\link{Results}} with \code{CIprop.Lg.Lcat.qt} and \code{CIcount.Lg.Lcat.qt}.
  name.res = "Country" ##<<Name used in csv file name AND to determine whether to
  ## save CIs for countries (Country) or UNPD aggregates (UNPDaggregate) if \code{res = NULL}.
  ){
  if (is.null(table.dir)) {
    table.dir <- file.path(getwd(), "tables/")
    dir.create(table.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)) {
    output.dir <- file.path(getwd(), "output", run.name, 
                            "/")
  }
  if (is.null(res)) {
    if (name.res == "Country") {
      load(file.path(output.dir, "res.country.rda")) # change JR, 20140418
      res <- res.country
    }
    if (name.res == "UNPDaggregate") {
      load(file.path(output.dir, "res.aggregate.rda")) # change JR, 20140418
      res <- res.aggregate
    }
  }
  est.years <- as.numeric(dimnames(res$CIprop.Lg.Lcat.qt[[1]][[1]])[[2]])
  percentiles <- as.numeric(dimnames(res$CIprop.Lg.Lcat.qt[[1]][[1]])[[1]])
  nperc <- length(percentiles)
  G <- length(res$CIprop.Lg.Lcat.qt)
  for (j in 1:3) {
    CI.Lg.Lcat.qt <- res[[c("CIprop.Lg.Lcat.qt", "CIcount.Lg.Lcat.qt", 
                            "CIratio.Lg.Lcat.qt")[j]]]
    for (i in 1:length(CI.Lg.Lcat.qt[[1]])) {
      CI.Lg.qt <- lapply(CI.Lg.Lcat.qt, function(l) l[[i]])
      estimates <- matrix(NA, nperc * G, length(est.years))
      for (g in 1:G) {
        indices <- seq((g - 1) * nperc + 1, g * nperc)
        estimates[indices, ] <- CI.Lg.qt[[g]]
      }
      results.all <- data.frame(paste(rep(names(CI.Lg.qt), 
                                          each = nperc)), paste(rep(res$iso.g, each = nperc)), 
                                paste(rep(percentiles, times = G)), estimates)
      dimnames(results.all)[[2]] <- c("Name", "Iso", "Percentile", 
                                      est.years)
      write.csv(results.all, file = file.path(table.dir, paste0(run.name, 
                                          "_", name.res, "_", c("perc", "count", "ratio")[j], 
                                          "_", paste(c(unlist(strsplit(names(CI.Lg.Lcat.qt[[1]])[i], 
                                                                       split = "/"))), collapse = "Over"), ".csv")), # change JR, 20140418
                row.names = FALSE)
    }
  }
  cat("Result tables written to", table.dir, "\n")
  return(invisible())
}
#----------------------------------------------------------------------------------
GetTablesChange <- function(# Save csv's with CIs for change in proportions and counts
  ### Save csv's with CIs for change in proportions and counts
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  table.dir = NULL, ##<< Directory to store tables. 
  ## If NULL, folder "tables" in current working directory.
  res = NULL, ##<< If NULL, country/UNPD aggregates are summarized, specified with next argument.
  ## Alternatively, an object of class \code{\link{Results}} with \code{changeprop.Lg.Lcat.qt} and \code{changecount.Lg.Lcat.qt}.
  name.res = "Country" ##<<Name used in csv file name AND to determine whether to
  ## save CIs for countries (Country) or UNPD aggregates (UNPDaggregate) if \code{res = NULL}.
){  
  if (is.null(table.dir)){
    table.dir <- file.path(getwd(), "tables/")
    dir.create(table.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name) # change JR, 20140418
  }
  if (is.null(res)){
    if (name.res == "Country"){
      load(file.path(output.dir, "res.country.rda"))
      res <- res.country
    }
    if (name.res == "UNPDaggregate"){
      load(file.path(output.dir, "res.aggregate.rda")) # change JR, 20140418
      res <- res.aggregate
      #print(res)
    }
  }  
  
  change.years.names <- dimnames(res$changeprop.Lg.Lcat.Ti[[1]][[1]])[[1]]
  infonames <- dimnames(res$changeprop.Lg.Lcat.Ti[[1]][[1]])[[2]]
  ninfo <- length(infonames)
  G <- length(res$changeprop.Lg.Lcat.Ti) # same for count and prop
  for (j in 1:2){
    CI.Lg.Lcat.Ti <- res[[c("changeprop.Lg.Lcat.Ti", "changecount.Lg.Lcat.Ti")[j]]]
    for (i in 1:length(CI.Lg.Lcat.Ti[[1]])){ # not the info i...
      CI.Lg.Ti <- lapply(CI.Lg.Lcat.Ti, function(l) l[[i]])
      estimates <- matrix(NA, length(change.years.names)*G, ninfo)
      for (g in 1:G){
        indices <- seq((g-1)*length(change.years.names) +1, g*length(change.years.names))
        estimates[indices,] <- CI.Lg.Ti[[g]]
      }
      results.all <- data.frame(
        paste(rep(names(CI.Lg.Ti), each = length(change.years.names))),
        paste(rep(res$iso.g, each = length(change.years.names))),
        paste(rep(change.years.names, times = G)),
        estimates)
      dimnames(results.all)[[2]] <- c("Name", "Iso", "Change", infonames)
      ##details<< Two csv-files (named \code{name.res_changes}) are saved: one with results for proportions, and one with results for counts.
      write.csv(results.all, file = file.path(table.dir, paste0(run.name,"_", name.res, "_changes_", c("perc", "count")[j], "_", names(CI.Lg.Lcat.Ti[[1]])[i],".csv")), # change JR, 20140418 
                row.names = FALSE)
    } # end i (cat) loop
  } # end j (prop or count) loop
  cat("Change tables written to", table.dir, "\n")
  return(invisible())
}
#----------------------------------------------------------------------------------
# SpitOutSomeRelevantNumbers <- function(# Spit out some relevant numbers
#   ## Info on MCMC chains and country trajectories # and world estimates
#   run.name = "test", ##<<,
#   output.dir = NULL ##<< directory where MCMC array and meta are stored, and new objects are added
#   ){
#   
#   
#   if (is.null(output.dir)){
#     output.dir <- file.path(getwd(), "output", run.name, "/")
#   }
#   
#   #load(file = file.path(output.dir,"par.ciq.rda")) # change JR, 20140418
#   load(file = file.path(output.dir,"res.country.rda")) # change JR, 20140418
#   #load(file = file.path(output.dir,"res.aggregate.rda")) # change JR, 20140418
#   load(file = file.path(output.dir,"mcmc.meta.rda")) # change JR, 20140418
#   
# #   poschange <- function(bla.i){
# #     return(paste(round(bla.i["50%"],2), 
# #                  ", 95% CI (", round(bla.i["2.5%"],2), ",",
# #                  round(bla.i["97.5%"],2), "), post. prob increase ", 
# #                  round(bla.i["PPPC"],3)
# #                  , sep = ""))
# #   }
# #   negchange <- function(bla.i){
# #     return(paste(round(-bla.i["50%"],2), 
# #                  ", 95% CI (", -round(bla.i["2.5%"],2), ",",
# #                  round(-bla.i["97.5%"],2), "), post. prob decrease ", 
# #                  round(1-bla.i["PPPC"],3)
# #                  , sep = ""))
# #   }
# #   levelcount <- function(bla.q){
# #     #NOTE: with .q no %!
# #     return(paste(round(bla.q["0.5"]), 
# #                  ", 95% CI (", round(bla.q["0.025"]), ",",
# #                  round(bla.q["0.975"]), ")", sep = ""))
# #   }
# #   levelprop <- function(bla.q){
# #     #NOTE: with .q no %!
# #     return(paste(round(bla.q["0.5"],2), 
# #                  ", 95% CI (", round(bla.q["0.025"],2), ",",
# #                  round(bla.q["0.975"],2), ")", sep = ""))
# #   }
#   
#   file.dir <-file.path(getwd(), "relevantnumbers.txt")
#   cat("Numbers to be written to", file.dir, "\n")
#   
#   cat("Info on MCMC algorithm from meta", "\n", file= file.dir, append = F)
#   cat("Number of iterations:", mcmc.meta$general$N.ITER, "\n", file=  file.dir, append = T)
#   cat("Number of chains:", length(mcmc.meta$general$ChainNums), "\n", file=  file.dir, append = T)
#   cat("Thinning:", mcmc.meta$general$N.THIN, "\n", file=  file.dir, append = T)
#   cat("Burnin:", mcmc.meta$general$N.BURNIN, "\n", file=  file.dir, append = T)
#   
#   cat("Info on samples used (from mcmc.array)", "\n", file= file.dir, append = T)
#   load(file = file.path(output.dir,"mcmc.array.rda")) # change JR, 20140418
#   cat("Number of samples:", dim(mcmc.array)[[1]], "\n", file=  file.dir, append = T)
#   cat("Number of chains:", dim(mcmc.array)[[2]], "\n", file=  file.dir, append = T)
#   
#   cat("Info on number of country trajectories used (and saved)", "\n", file= file.dir, append = T)
#   load(file = file.path(res.country$output.dir.countrytrajectories, "P.tp3s_country", 1, ".rda")) # change JR, 20140418
#   cat("Number of trajectories:", dim(P.tp3s)[[3]], "\n", file=  file.dir, append = T)
#   
#   cat("", "\n", file=  file.dir, append = T)
#   
# #   cat("Info on world estimates", "\n", file= file.dir, append = T)
# #   cat("INCREASE in total prevalence, 1990-2010:", poschange(res.aggregate$changeprop.Lg.Lcat.Ti[["World"]][["Total"]]["1990-2000",]),
# #       "\n", file = file.dir, append = T)      
# #   cat("DECREASE in % unmet need, 1990-2010:", negchange(res.aggregate$changeprop.Lg.Lcat.Ti[["World"]][["Unmet"]]["1990-2000",]),
# #       "\n", file = file.dir, append = T)      
# #   cat("Total MWRA with unmet need, 2010:", levelcount(
# #     res.aggregate$CIcount.Lg.Lcat.qt[["World"]][["Unmet"]][, "2010.5"]),
# #       "\n", file = file.dir, append = T)      
# # cat("", "\n", file= file.dir, append = T)
#   
#   
#   return(invisible())
# }
#----------------------------------------------------------------------
# The End!
