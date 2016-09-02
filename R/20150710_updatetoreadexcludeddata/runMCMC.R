#--------------------------------------------------------------------------
# F_runMCMC.R
# Leontine Alkema and Jin Rou New
#--------------------------------------------------------------------------
RunMCMC <- function(# Start MCMC sampling
  ### Start MCMC sampling for the CP model and save mcmc.meta and JAGS objects to \code{output.dir}.
  run.name = "test", ##<< Run name, used to create a directory \code{output/run.name}
  ## with JAGS output (after MCMC sampling), and estimates (in next steps).
  N.ITER = ifelse(!do.country.specific.run, 80000, 40000), ##<< Number of iterations, NOT including burn-in.
  N.STEPS = 4, ##<< For each N.ITER/N.STEPS iterations, the iterations will be saved. 
  N.THIN = ifelse(!do.country.specific.run, 30, 15), ##<< Thinning factor.
  N.BURNIN = 20000, ##<< Burnin (excluded samples at start of chain).
  ChainNums = seq(1,5), ##<< IDs of chains to run in series 
  ## (the IDs need to be numeric because they are used to set the seed).
  do.country.specific.run = FALSE, ##<< Logical: execute a country-specific (as opposed to global) run? # change JR, 20131104
  do.country.specific.targets.run = FALSE, ##<< Logical: execute a country-specific (as opposed to global) run for targets? # change JR, 20150301
  iso.select = NULL, ##<< (For country/subpopulation-specific run) Numeric or 3-character ISO country code for country/subpopulation to select 
  ## one country/subpopulation to run the model for, or if NULL, all countries/subpopulations in data are selected. Should be of length 1 
  ## if \code{do.country.specific.run} is \code{TRUE} # change JR, 20131104
  iso.country.select = NULL, ##<< (For subpopulation-specific run) Numeric or 3-character ISO country code for country that subpopulation
  ## belongs to. # change JR, 20140404
  run.name.global = NULL, ##<< (For country-specific run) Run name of global run # Change JR, 20131104
  change.priors.to.zerolower = FALSE, ##<< Logical: Change priors?
  ## FALSE: use gammas on kappa.c's, TRUE: use uniform priors.
  include.AR = TRUE , ##<< Logical: include AR(1)'s for total, modern/total and unmet?
  seed.MCMC = 1, ##<< seed for initializing MCMC, defaults to 1.
  output.dir = NULL, ##<< Directory where mcmc meta and raw MCMC output will be stored
  ##either an existing directory, or if NULL, directory \code{output/run.name} is created
  ## in current working directory
  data.csv = NULL, ##<< If \code{NULL}, contraceptive use data set included in package is used. 
  ## To use alternative data set, use \code{data.csv = .../dataCPmodel.csv}, where ... refers to the file path where file is located.
  regioninfo.csv = NULL, ##<< If \code{NULL}, region info included in package is used. 
  ## To use alternative csv file, use \code{regioninfo.csv = .../Country-and-area-classification.csv}, where ... refers to the file path where file is located.
  html.file = NULL,##<<If not NULL, summary results about the data set get written to this HTML file.
  exclude.unmet.only = FALSE, ##<< Logical: do validation for unmet need only?  (see description on validation exercises below)
  at.random = FALSE, ##<< Logical: do validation, leaving out obs at random?
  at.end = FALSE, ##<< Logical: do validation, leaving out obs at end? 
  year.cutoff = 2005,##<< Used only if \code{at.end} = \code{TRUE}):
  ## All data with observation year after and IN year.cutoff is excluded.
  seed.validation = 12345, ##<< For constructing training set in validation exercise (if applicable)
  generate.new.set = TRUE, ##<<Logical: generate a new training set in validation exercise?
  run.jags = TRUE, ##<< Logical: run JAGS?
  run.on.server = TRUE ##<< Logical: run on server (in parallel?)
){
  ChainNums <- unique(ChainNums)
  
  ##details<< All output is written to folder output.dir (you'll get a message if it already exists).
  ## JAGS objects are written to output.dir/temp.JAGSobjects.
  if (is.null(output.dir)){
    dir.create(file.path(getwd(), "output"), showWarnings = FALSE) 
    output.dir <-  file.path(getwd(), "output", run.name) # change JR, 20140414
  }
  if (file.exists(file.path(output.dir, "mcmc.meta.rda"))){
    cat(paste("The output directory", output.dir, "already contains an MCMC.meta (and probably an older MCMC run).\n"))
    cat(paste("Delete files in this directory or choose a different run.name/output directory,\n"))
    return(invisible())
  }
  dir.create(output.dir)
  if (run.jags) # change JR, 20131105
    dir.create(file.path(output.dir, "temp.JAGSobjects"))#, showWarnings = FALSE) 
  
  filename <- file.path(output.dir, "logfile.txt") # change JR, 20140418
  cat(paste("Seed and data files used are written to logfile ", filename), "\n")
  fileout <- file(filename, open = "wt")
  sink(fileout, split = T)
  
  #if (is.null(seed.MCMC)) seed.MCMC <- as.numeric(Sys.Date())
  # Note about Sys.Date, to reproduce it, use 
  # as.numeric(as.Date("2012-03-09"))
  cat(paste("seed.MCMC is", seed.MCMC), "\n")
  
  if (!is.null(iso.select)) { # change JR, 20131104
    if (!(all(grepl("[[:alpha:]]", iso.select)) | all(grepl("[[:digit:]]", iso.select)))) {
      error("iso.select: Must be NULL, else numeric or character ISO country code.")
      return(invisible())
    }
  }
  
  if (!is.null(iso.country.select)) { # change JR, 20140404
    if (!(all(grepl("[[:alpha:]]", iso.country.select)) | all(grepl("[[:digit:]]", iso.country.select)))) {
      error("iso.country.select: Must be NULL, else numeric or character ISO country code.")
      return(invisible())
    }
  }
  
  ##details<< Object \code{data.global} is loaded or created, which is NULL if this run is not country-specific.
  if (do.country.specific.run | do.country.specific.targets.run) { # change JR, 20150301
    if (is.null(run.name.global)) {
      run.name.global <- "Run20140520" # change JR, 2010612
      if (!file.exists(file.path("data/data.global.rda"))) {
        cat(paste0("Error: No default data.global file in data folder. Run global run first or specify run.name.global!\n"))
        return(invisible())
      } else {
        load(file.path("data/data.global.rda"))
        cat(paste0("Default global run loaded.\n"))
      }
      if (do.country.specific.targets.run) { # change JR, 20150301
        if (!file.exists(file.path("data/data.logratios.rda"))) {
          cat(paste0("Error: No default data.logratios file in data folder. Run global run first or specify run.name.global!\n"))
          return(invisible())
        } else {
          load(file.path("data/data.logratios.rda"))
          cat(paste0("Default data.logratios loaded.\n"))
        }
      }
    } else {
      if (!file.exists(file.path("output", run.name.global, "data.global.rda")))
        SummariseGlobalRun(run.name = run.name.global)
      load(file.path("output", run.name.global, "data.global.rda"))
      cat(paste0("Global run ", run.name.global, " loaded.\n"))
      if (do.country.specific.targets.run) { # change JR, 20150301
        load(file.path("output", run.name.global, "data.logratios.rda"))
        cat(paste0("Log ratio info list of ", run.name.global, " loaded.\n"))
      } else {
        data.logratios <- NULL # change JR, 20150301
      }
    }
  } else {
    data.global <- NULL
    data.logratios <- NULL # change JR, 20150301
  }
  
  ##details<< Object \code{validation.list} is created, which is NULL if this run is not a validation exercise.
  ## If several options in validation exercise were set to TRUE (in arguments of this function),
  ## observations are left out at random (first choice),
  ## or left out at the end (second choice).
  if (exclude.unmet.only|at.end|at.random) {
    if (exclude.unmet.only*at.random) exclude.unmet.only = FALSE
    if (at.end*at.random) at.end = FALSE
    if (at.end*exclude.unmet.only) exclude.unmet.only = FALSE
    do.validation <- TRUE
    #  ##details<< If it is a validation exercise, \code{validation.list} includes:
    #  ##describe<< 
    validation.list <- list(do.validation = do.validation,# ##<< TRUE
                            exclude.unmet.only = exclude.unmet.only,# ##<< From arguments
                            at.random = at.random,###<< From arguments
                            at.end = at.end,#  ##<< From arguments
                            year.cutoff = year.cutoff,# ##<< From arguments
                            seed = seed.validation,# ##<< From arguments
                            generate.new.set = generate.new.set # ##<< From arguments
    )
    # ##end<<
  } else {
    validation.list <- NULL
  }
  # Data 
  if (is.null(iso.country.select)) { # change JR, 20140409
    name.country.select <- NULL
  } else {
    name.country.select <- names(data.global$iso.c)[match(iso.country.select, data.global$iso.c)]  
  }
  
  data.preprocessed <- PreprocessData(data.csv = data.csv,
                                      iso.select = iso.select)
  
  select.ss <- grepl("Service statistic", data.preprocessed$Data.series.type[
    data.preprocessed$EXCLUDE1isyes == 0]) # change JR, 20150710
  if (!do.country.specific.run & any(select.ss))
    stop("SS data should not be used for global run!")
  do.SS.run <- do.country.specific.run & any(select.ss)
  if (do.SS.run & sum(select.ss) == 1) 
    stop("Only 1 observation of SS data is available. Not possible to use SS data for projection.")
  if (do.SS.run) {
    years.all <- (data.preprocessed$Start.year + data.preprocessed$End.year)[
      data.preprocessed$EXCLUDE1isyes == 0]/2 # change JR, 20150710
    max.survey.year <- max(years.all[!select.ss])
    diff.years.ss <- years.all[select.ss] - max.survey.year
    if (all(diff.years.ss <= 0))
      stop("Service statistics data do not extend beyond survey data.")
    if (length(diff.years.ss >= 0) == 1)
      stop("Only 1 observation of service statistics data at/beyond the latest survey data year is available. More SS data required!")
  }
  # for run with SS data, check if this is the first/second pass run
  do.SS.run.first.pass <- do.SS.run & !file.exists(file.path(output.dir, "res.country_pre.rda"))
  do.SS.run.second.pass <- do.SS.run & file.exists(file.path(output.dir, "res.country_pre.rda"))
  
  if (do.SS.run.second.pass) {
    cat("Starting the second pass run for model with service statistics data...\n")
    load(file.path(output.dir, "res.country_pre.rda"))
    # get SS data point, either the closest SS obs prior to/in most recent survey year if available, else the closest one after that year
    if (any(diff.years.ss <= 0)) { # if SS data overlaps with non-SS data
      year.ss <- max(years.all[select.ss & years.all <= max.survey.year])
    } else {
      year.ss <- min(years.all[select.ss])
    }
    modern.CP.ss <- data.preprocessed$Contraceptive.use.MODERN[
      data.preprocessed$EXCLUDE1isyes == 0][select.ss & years.all == year.ss]/100 # change JR, 20150710
    cat(paste0("The most recent observation year of non-SS data is ", max.survey.year, ".\n"))
    cat(paste0("The service statistics data observation in ", year.ss, " is used to calculate the relative SS bias.\n"))
    # get bias.modern
    modern.CP.est <- res.country$CIprop.Lg.Lcat.qt[[1]][["Modern"]]["0.5", ]
    years.est <- as.numeric(names(modern.CP.est))
    modern.CP.ssyear <- approx(x = years.est, y = modern.CP.est, xout = year.ss)$y
    bias.modern <- (modern.CP.ss*modern.CP.ssyear - modern.CP.ss)/(modern.CP.ss*modern.CP.ssyear - modern.CP.ssyear)
    cat(paste0("The relative service statistics bias is ", bias.modern, ".\n"))
    data.SS <- list(bias.modern = bias.modern)
    rm(res.country)
  } else {
    if (do.SS.run.first.pass)
      cat("Starting the first pass run for model with service statistics data...\n")
    data.SS <- NULL
  }
  
  filename.append <- ifelse(do.SS.run.first.pass, "_pre", "")
  # save a copy of preprocessed data file for Shiny
  write.csv(data.preprocessed, file = file.path(output.dir, paste0("dataCPmodel_input", filename.append, ".csv")), row.names = F)
  cat(paste0("Pre-processed data saved to ", file.path(output.dir, paste0("dataCPmodel_input", filename.append, ".csv")), "\n"))
  
  data.raw <- ReadDataAll(data.csv = data.csv,
                          regioninfo.csv = regioninfo.csv,
                          iso.select = iso.select,
                          iso.country.select = iso.country.select,
                          name.country.select = name.country.select,
                          do.SS.run.first.pass = do.SS.run.first.pass,
                          html.file = html.file)
  # change JR, 20150710
  data.raw.excluded <- ReadDataAll(data.csv = data.csv,
                                   regioninfo.csv = regioninfo.csv,
                                   iso.select = iso.select,
                                   iso.country.select = iso.country.select,
                                   name.country.select = name.country.select,
                                   do.SS.run.first.pass = do.SS.run.first.pass,
                                   get.excluded.data = TRUE,
                                   html.file = html.file)
  winbugs.data <- GetBugsData(data = data.raw$data, 
                              country.info = data.raw$country.info,
                              data.global = data.global,
                              data.SS = data.SS, # change JR, 20140414
                              data.logratios = data.logratios, # change JR, 20150301
                              output.dir = output.dir, 
                              validation.list  = validation.list,
                              do.country.specific.run = do.country.specific.run,
                              do.SS.run.second.pass = do.SS.run.second.pass, # change JR, 20140414
                              do.country.specific.targets.run = do.country.specific.targets.run, # change JR, 20150301
                              change.priors.to.zerolower = change.priors.to.zerolower)
  # save names of V parameters separately
  # this function is used to find out which parameter was assigned to which country
  # (through indices of observations)
  if (!do.country.specific.targets.run) {
    par.V <- InternalGetParnamesV(winbugs.data = winbugs.data,
                                  name.short.j = InternalMakeCountryNamesShort(data.raw$data$name.j))
  } else {
    par.V <- NULL
  }
  parnames.list <- GetParNames(winbugs.data = winbugs.data, 
                               validation.list = validation.list,
                               do.country.specific.run = do.country.specific.run, # change JR, 20131104
                               do.country.specific.targets.run = do.country.specific.targets.run) # change JR, 20150301
  
  ##details<<
  ##describe<< Object mcmc.meta is saved, which is a list with
  mcmc.meta <- list(
    general = ##<< General info:
      ##describe<<
      list(N.ITER = N.ITER,##<< From arguments
           ChainNums = ChainNums, ##<< ##<< From arguments (updated later if chains are added)
           N.STEPS = N.STEPS, ##<< From arguments
           N.THIN = N.THIN, ##<< From arguments
           N.BURNIN = N.BURNIN, ##<< From arguments
           seed.MCMC = seed.MCMC, ##<< From arguments
           do.country.specific.run = do.country.specific.run, ##<< From arguments # change JR, 20131104
           do.SS.run.first.pass = do.SS.run.first.pass, # change JR, 20140414
           do.country.specific.targets.run = do.country.specific.targets.run, # change JR, 20150301
           run.name.global = run.name.global, ##<< From arguments # change JR, 20131104
           change.priors.to.zerolower = change.priors.to.zerolower, ##<< From arguments
           output.dir = output.dir##<< From arguments
      ),
    ##end<<
    parnames.list = parnames.list, ##<< Output from \code{GetParNames}, list with all parnames in BUGS
    par.V = par.V, ##<< Output from \code{InternalInternalGetParnamesV}, BUGS and ``nice'' parameter names for data multipliers 
    validation.list = validation.list, ##<< Details about the validation exercise
    include.AR = include.AR,##<< From arguments
    winbugs.data = winbugs.data, ##<< Object from \code{\link{GetBugsData}}
    data.raw = data.raw, ##<< Object from \code{\link{ReadDataAll}}
    data.global = data.global, ##<< Object from \code{\link{SummariseGlobalRun}} # change JR, 20131104
    data.SS = data.SS, ##<< Object summarising information from first pass of run with SS data # change JR, 20140414
    data.logratios = data.logratios ##<< Object summarising information about log ratios # change JR, 20150301
  )
  ##end<<
  save(mcmc.meta, file = file.path(output.dir, paste0("mcmc.meta", filename.append, ".rda"))) # change JR, 20140414
  ##details<< See \code{\link{AddMCMCChain}} for adding an additional chain (with same number of iterations etc).
  #Details<< BUGS model is stored in \code{output.dir} using \code{WriteModel} or \code{WriteCountryModel} 
  ### if do.country.specific.run is \code{TRUE}.
  if (do.country.specific.run & !do.country.specific.targets.run) { # change JR, 20150301
    WriteCountryModel(mcmc.meta = mcmc.meta)
  } else if (do.country.specific.run & do.country.specific.targets.run) { # change JR, 20150301
    WriteCountryModelForTargets(mcmc.meta = mcmc.meta)
  } else {
    WriteModel(mcmc.meta = mcmc.meta)
  }
  #file.show(file.path(output.dir, "model.txt"))
  
  sink()
  closeAllConnections()
  
  if (run.jags) {
    if (run.on.server) {
      foreach(chainNum=ChainNums) %dopar% {
        cat(paste("Start chain ID ", chainNum), "\n")
        #tryCatch(
        InternalRunOneChain(chainNum = chainNum, mcmc.meta = mcmc.meta)
        #                   , warning=function(w, output.dir = mcmc.meta$general$output.dir){
        #                   file.w = file.path(output.dir, "warnings.txt");
        #                   cat(paste(w), file=file.w, append=T); return()},#;  close(file.w)},
        #                   error=function(e){#},  output.dir = mcmc.meta$general$output.dir){
        #                   file.e=file.path(output.dir, "errors.rda");
        #                   save(e, file = "file.e"); # does not give anything.... 
        #                   return()},                
        #                   #cat(paste(unlist(e)), file = file.e,append=T); return()},#; close(file.e)},                
        #                   #sink(file=file.path(mcmc.meta$general$output.dir, "errors.txt"),append=T);print("Errors:");print(e);sink();return(NULL)},
        #                   finally=function(j){ return()}) # does not do anything either...
      } # end chainNums
    } else {
      for (chainNum in ChainNums){
        cat(paste("Start chain ID ", chainNum), "\n")
        #tryCatch(
        InternalRunOneChain(chainNum = chainNum, mcmc.meta = mcmc.meta)
        #                   , warning=function(w, output.dir = mcmc.meta$general$output.dir){
        #                   file.w = file.path(output.dir, "warnings.txt");
        #                   cat(paste(w), file=file.w, append=T); return()},#;  close(file.w)},
        #                   error=function(e){#},  output.dir = mcmc.meta$general$output.dir){
        #                   file.e=file.path(output.dir, "errors.rda");
        #                   save(e, file = "file.e"); # does not give anything.... 
        #                   return()},                
        #                   #cat(paste(unlist(e)), file = file.e,append=T); return()},#; close(file.e)},                
        #                   #sink(file=file.path(mcmc.meta$general$output.dir, "errors.txt"),append=T);print("Errors:");print(e);sink();return(NULL)},
        #                   finally=function(j){ return()}) # does not do anything either...
      }
    }
    #  cat("Any JAGS related errors/warnings are written to files in output.dir, but can usually be ignored.. (if there is something wrong, you'll find out soon enough in the next steps :))", "\n")
    cat("All chains have finished!\n")
  }
  ##value<< NULL (mcmc.meta is saved to output.dir, and JAGS objects are saved in their own directory).
}# end function

#-----------------------------------------------------
InternalRunOneChain <- function(#Do MCMC sampling
  ###Do MCMC sampling for one chain
  #, use with tryCatch function to avoid zillion errors from JAGS...
  chainNum, ##<< Chain ID
  mcmc.meta ##<< List, described in \code{\link{RunMCMC}}
){
  # set seed before sampling the initial values
  set.seed.chain <- chainNum*mcmc.meta$general$seed.MCMC*
    min(as.numeric(mcmc.meta$data.raw$country.info$iso.c[1]), na.rm = T) # change JR, 20140617: added min # change JR, 20140310: added seed
  # in Jags version (July 21, 2012) jags.seed doesn't work, and inits need to be provided as a function
  # note that even with same Jags seed, as long as inits from R are different, any non-initialized pars will haev different starting values
  # and seed in Jags is consistent
  mcmc.info <- list(set.seed.chain = set.seed.chain, chainNum = chainNum)
  filename.append <- ifelse(mcmc.meta$general$do.SS.run.first.pass, "_pre", "") # change JR, 20140414
  mcmc.info.file <- file.path(mcmc.meta$general$output.dir, paste0("mcmc.info", filename.append, ".", chainNum, ".rda")) # change JR, 20140414 # change JR, 20140418
  if (file.exists(mcmc.info.file)){
    cat(paste("The output directory", mcmc.meta$general$output.dir, "already contains info on chain", chainNum))
    cat(paste("No new samples are added"))
    return(invisible())
  }
  save(mcmc.info, file = mcmc.info.file)
  cat("JAGS is called to obtain posterior samples, and there will be some info about the model and steps written to file.", "\n")
  #  cat("And lots of error messages (potentially), which I do not know how to get rid off", "\n")
  cat("Just wait for statement that MCMC run has finished", "\n")    
  jags.dir <- file.path(mcmc.meta$general$output.dir, "temp.JAGSobjects/")
  #   mod<-tryCatch(jags(data=mcmc.meta$winbugs.data, 
  #               inits=fix.init,
  #               parameters.to.save=unlist(mcmc.meta$parnames.list), 
  #               model.file=file.path(mcmc.meta$general$output.dir, "model.txt"),
  #               n.chains=1, 
  #               n.iter=mcmc.meta$general$N.BURNIN+mcmc.meta$general$N.ITER/mcmc.meta$general$N.STEPS, 
  #               n.burnin=mcmc.meta$general$N.BURNIN,
  #               n.thin=mcmc.meta$general$N.THIN,
  #               DIC=FALSE, 
  #               working.directory=mcmc.meta$general$output.dir),
  # 
  #               warning=function(w, output.dir = mcmc.meta$general$output.dir){
  #                 file.w = file.path(output.dir, "warnings.txt");
  #                 cat(paste(unlist(w)), file=file.w, append=T); return()},#;  close(file.w)},
  #               error=function(e){#},  output.dir = mcmc.meta$general$output.dir){
  #                 #file.e=file.path(output.dir, "errors.rda");
  #                 save(e, file = "bla.rda");#file.e); 
  #                 return()},#; close(file.e)},                
  #                 #cat(paste(unlist(e)), file = file.e,append=T); return()},#; close(file.e)},                
  #                 #sink(file=file.path(mcmc.meta$general$output.dir, "errors.txt"),append=T);print("Errors:");print(e);sink();return(NULL)},
  #               finally=function(mod.upd, i.temp =1, chainNum.temp = chainNum, jags.dir.temp = jags.dir){
  #                 save(mod.upd,file=paste(jags.dir.temp, "/jags_mod", chainNum.temp, "update_", i.temp, ".Rdata", sep = ""));
  #                 return(mod.upd)} ) # return doesn't help much...
  set.seed(set.seed.chain) # note: seed only useful if inits function didn't change!
  # need to sample something in R first
  temp <- rnorm(1)
  mod<-jags(data=mcmc.meta$winbugs.data, 
            # inits=fix.init,
            # note this functionality does not work as of July 21, 2012
            # instead: need a function without input arguments, that samples 1!
            inits= InternalMCMCinits(
              winbugs.data = mcmc.meta$winbugs.data,
              do.country.specific.run = mcmc.meta$general$do.country.specific.run, # change JR, 20131104
              do.country.specific.targets.run = mcmc.meta$general$do.country.specific.targets.run, # change JR, 20150301
              change.priors.to.zerolower = mcmc.meta$general$change.priors.to.zerolower),
            parameters.to.save=unlist(mcmc.meta$parnames.list), 
            model.file=file.path(mcmc.meta$general$output.dir, paste0("model", filename.append, ".txt")), # change JR, 20140414
            n.chains=1, 
            n.iter=mcmc.meta$general$N.BURNIN+mcmc.meta$general$N.ITER/mcmc.meta$general$N.STEPS, 
            n.burnin=mcmc.meta$general$N.BURNIN,
            n.thin=mcmc.meta$general$N.THIN,
            DIC=FALSE, 
            # jags.seed = set.seed.chain, 
            # note this functionality does not work as of July 21, 2012
            jags.seed = set.seed.chain, # change JR, 20140617: changed from 123
            working.directory=mcmc.meta$general$output.dir)
  # in theory... to update in the future: load("jags_mod.Rdata");recompile(mod)
  # in practice.... that never worked for me!
  
  i = 1 # index for which update
  mod.upd <- mod
  save(mod.upd, file=file.path(mcmc.meta$general$output.dir, "temp.JAGSobjects", paste0("jags_mod", filename.append, chainNum, "update_", i, ".Rdata"))) # change JR, 20140414 # change JR, 20140418
  #load(file=file.path(mcmc.meta$general$output.dir, "temp.JAGSobjects", paste0("jags_mod", filename.append, chainNum, "update_", i, ".Rdata")) # change JR, 20140418
  cat(paste("MCMC results step", 1, " for chain ", chainNum, " written to folder temp.JAGSobjects in ", mcmc.meta$general$output.dir), "\n")
  
  #--- update MCMC ----------
  if (mcmc.meta$general$N.STEPS >1){
    for (i in 2:(mcmc.meta$general$N.STEPS)){
      mod.upd <-update(mod.upd, parameters.to.save=unlist(mcmc.meta$parnames.list), # change JR, 20131104
                       n.iter=mcmc.meta$general$N.ITER/mcmc.meta$general$N.STEPS,
                       n.thin=mcmc.meta$general$N.THIN)
      save(mod.upd, file = file.path(mcmc.meta$general$output.dir, "temp.JAGSobjects", paste0("jags_mod", filename.append, chainNum, "update_", i, ".Rdata"))) # change JR, 20140414 # change JR, 20140418
      #load(file = file.path(mcmc.meta$general$output.dir, "temp.JAGSobjects", paste0("jags_mod", filename.append, chainNum, "update_", i, ".Rdata"))) # change JR, 20140414 # change JR, 20140418
      cat(paste("MCMC results step", i, " for chain ", chainNum, " written to folder temp.JAGSobjects in ", mcmc.meta$general$output.dir), "\n")
    }
  }  
  #cat("Ignore all errors above....")
  cat(paste("Hooraah, Chain", chainNum, "has finished!"), "\n")
  ##note<< Called from \code{\link{RunMCMC}} and \code{\link{AddMCMCChain}}.
  ## This function can give errors and warnings when JAGS output is read into R,
  ## no worries about that here, convergence will be checked later.
  ##value<< NULL
  return(invisible())
}
#----------------------------------------------------------------------------------
AddMCMCChain <- function(# Add additional MCMC chain to existing run.
  ### Add additional MCMC chain with same specs as chains from meta from \code{run.name}.
  run.name = "test", ##<< run name, usually Frunnumber, from run where chain(s) should be added
  ChainNums = 6, ##<< IDs of additional chains to add, mcmc.meta is updated with additional chain info.
  ## If one or more of the IDs were already there, a message will be printed and no chains will be added.
  do.SS.run.first.pass = FALSE, ##<< is this for the first pass run with SS data?
  output.dir = NULL ##<< Directory where mcmc meta of run name was stored.
  ##  If NULL, it's \code{output/run.name/} in the current working directory.
){
  ##details<< See \code{\link{RunMCMC}} for initial run.
  ## This function will crash if you specified a run for which mcmc.meta has not yet been constructed
  output.dir <- file.path(getwd(), "output", run.name) # change JR, 20140418
  filename.append <- ifelse(do.SS.run.first.pass, "_pre", "")
  load(file.path(output.dir, paste0("mcmc.meta", filename.append, ".rda"))) # change JR, 20140418
  if (sum(is.element(ChainNums, mcmc.meta$general$ChainNums))>0){
    ChainNums <- setdiff(ChainNums, mcmc.meta$general$ChainNums)
    if (sum(ChainNums)==0){
      cat("MCMC run(s) for ChainNum(s) and run.name already exist(s)!", "\n")
      return(invisible())
    }
  }
  # add chain info to mcmc.meta
  
  mcmc.meta$general$ChainNums <- unique(c(mcmc.meta$general$ChainNums, ChainNums))
  save(mcmc.meta, file = file.path(output.dir, paste0("mcmc.meta", filename.append, ".rda"))) # change JR, 20140418
  for (chainNum in ChainNums){
    #tryCatch(
    InternalRunOneChain(chainNum = chainNum, mcmc.meta = mcmc.meta)
    #              , warning=function(w, output.dir = mcmc.meta$general$output.dir){
    #                file.w = file.path(output.dir, "warnings.txt");
    #                cat(paste(w), file=file.w, append=T); return()},#;  close(file.w)},
    #              error=function(e){#},  output.dir = mcmc.meta$general$output.dir){
    #                file.e=file.path(output.dir, "errors.rda");
    #                save(e, file = "file.e"); # does not give anything.... 
    #                return()},                
    #              #cat(paste(unlist(e)), file = file.e,append=T); return()},#; close(file.e)},                
    #              #sink(file=file.path(mcmc.meta$general$output.dir, "errors.txt"),append=T);print("Errors:");print(e);sink();return(NULL)},
    #              finally=function(j){ return()}) # does not do anything either...
  } # end chainNums
  cat("Any JAGS related errors/warnings are written to txt files in output.dir, but can usually be ignored.. (if there is something wrong, you'll find out soon enough in the next steps :))", "\n")
  cat("All additional chains have finished!")
  ##value<< NULL
}# end function
#----------------------------------------------------------------------------------
InternalMCMCinits <- function(# Initialize MCMC run in BUGS
  ### Initialize MCMC run in Bugs
  winbugs.data, ## Object of \code{\link{winbugs.data}}
  do.country.specific.run, ##<< Logical # change JR, 20131104
  do.country.specific.targets.run, ##<< Logical # change JR, 20150301
  change.priors.to.zerolower ##<< Logical
){
  list.varpar <- NULL
  inits.list1 <- list(
    T.c = c(rtnorm(winbugs.data$C, 1980, 80, lower=1800), NA), # change JR, 2013119: changed from winbugs.data$C+1
    RT.c = c(rtnorm(winbugs.data$C, 1980, 80, lower=1800), NA), # change JR, 2013119: changed from winbugs.data$C+1
    logitRomega.c = logit(runif(winbugs.data$C, 0.02, 0.25)),
    logitomega.c = logit(runif(winbugs.data$C, 0.02, 0.25)),
    logitRmax.c = LogitMinMax(runif(winbugs.data$C, 0.55, 0.95), xmin = 0.5, xmax = 1),
    logitpmax.c = LogitMinMax(runif(winbugs.data$C, 0.55, 0.95), xmin = 0.5, xmax = 1)
  )
  if (!do.country.specific.run & !do.country.specific.targets.run) { # change JR, 20150301
    list.varpar <- NULL # change JR, 20150301
    inits.list1 <- c(inits.list1, list(
      # "T.s[1,,]" =  rwish(v = 2, S = solve(diag(0.1,2))),
      # "T.s[2,,]" =  rwish(v = 2, S = solve(diag(0.1,2))),
      # "T.s[3,,]" =  rwish(v = 2, S = solve(diag(0.1,2))),
      # "T.s[4,,]" =  rwish(v = 2, S = solve(diag(0.1,2))),
      # # R specs: T ~ W(k,S) where T is prec. matrix, with E(T) = k*S, thus S = Vinv 
      tau.sourcetot = 1/runif(1,0.05, 1)^2,
      sigma.unmetworld = runif(1,0.05,1),
      mu.pos.m = rnorm(2, -2, 2), 
      sigma.pos = runif(1, 0.01, 2),  
      # tau.pos.m = 1/runif(2, 0.01, 1)^2,  
      sigma.geo.m = runif(2, 0.01, 2),  
      # note: first 3 refer to non-rich countries
      w.world = rtnorm(1,logit(0.07), 1,  lower=-4.5, upper = 0),
      T.world = rtnorm(1, 1980, 40,  lower=1800),
      Tearlier = rtnorm(1, 1900, 50,  lower=1800),
      RT.world = rtnorm(1, 1980, 40,  lower=1800),
      Rw.world = rtnorm(1,logit(0.07), 1, lower=-4.5, upper = 0),
      w.reg = rtnorm(winbugs.data$n.reg, log(0.07), sd = 1.5,  lower=-4.5, upper = 0), 
      T.reg = rtnorm(winbugs.data$n.reg, 1980, 50,  lower=1800),
      RT.reg = rtnorm(winbugs.data$n.reg, 1980, 50,  lower=1800),
      Rw.reg = rtnorm(winbugs.data$n.reg, log(0.07), sd = 1.5,  lower=-4.5, upper = 0),
      sigma.Treg = runif(1, 2, 80),
      sigma.RTreg = runif(1, 2, 80),
      sigma.wreg = runif(1, 0.01, 2),
      sigma.Rwreg = runif(1, 0.01, 2),
      sigma.Tsubreg = runif(1, 2, 80),
      sigma.RTsubreg = runif(1, 2, 80),
      sigma.wsubreg = runif(1, 0.01, 2),
      sigma.Rwsubreg = runif(1, 0.01, 2),
      lp.world = rnorm(1, -0.5, 1), 
      lr.world  = rnorm(1, 2, 2),
      sigma.unmet.dhs = runif(1, 0.01, 1),
      sigma.unmet.other = runif(1, 0.01, 1),
      sigma.ar.unmet  = runif(1,0.01,1),
      rho.unmet  = runif(1,0,1),
      sigma.tot  = runif(1,0.01,1),
      sigma.rat  = runif(1,0.01,1),
      rho.tot  = runif(1,0,1),
      rho.rat  = runif(1,0,1),
      a.unmet = rnorm(1,0,1),
      b.unmet = rnorm(1,0,1),
      c.unmet = runif(1,-7,0),
      v.mics = runif(1),v.folk = runif(1),v.mneg = runif(1),v.mpos = runif(1)
    ))
    if (!do.country.specific.run) { # change JR, 20150301
      if (change.priors.to.zerolower){
        list.varpar <- c(list.varpar, list(
          sigma.lpc = runif(1,0,5),
          sigma.lrc = runif(1,0,5),
          sigma.wc = runif(1, 0,2),
          sigma.Rwc = runif(10,2),
          sigma.Tc = runif(1,0,30),
          sigma.RTc = runif(1,0,30),
          sigma.earlierTc = runif(1,0,70),
          sigma.unmetc = runif(1,0,5) ))
      } else {
        list.varpar <- c(list.varpar, list(
          # eps.ci  = matrix(rnorm(winbugs.data$C*, mean = 0, sd = ),..,..),
          # eta.ci  = rnorm(winbugs.data$C*, mean = 0, sd = ),
          # theta.ci  = rnorm(winbugs.data$C*, mean = 0, sd = ),
          tau.lpc = 1/runif(1,0,5)^2,
          tau.lrc = 1/runif(1,0,5)^2,
          tau.wc = 1/runif(1,0,2)^2,
          tau.Rwc = 1/runif(1,0,2)^2, # change JR 20131101
          tau.Tc= 1/runif(1,0,30)^2,
          tau.RTc = 1/runif(1,0,30)^2,
          tau.earlierTc = 1/runif(1,0,70)^2,
          tau.unmetc = 1/runif(1,0,5)^2 ))
      }
      inits.list1 <- c(inits.list1, list(
        T.subreg = rtnorm(winbugs.data$n.subreg, 1980, 60,  lower=1800),
        w.subreg = rtnorm(winbugs.data$n.subreg, logit(0.07), sd = 1.5,  lower=-4.5, upper = 0), 
        RT.subreg = rtnorm(winbugs.data$n.subreg, 1980, 60,  lower=1800),
        Rw.subreg = rtnorm(winbugs.data$n.subreg, logit(0.07), sd = 1.5,  lower=-4.5, upper = 0)
      ))
    }
  }
  inits.list <- c(inits.list1, list.varpar)
  ##value<<List with initial values for model parameters 
  return(list(inits.list))
} # end inits
#----------------------------------------------------------------------
GetBugsData <- function( # Construct winbugs.data object
  ### Construct winbugs.data object
  data, ##<< class data
  country.info, ##<< class country.info
  data.global = NULL, ##<< class data.global
  data.SS = NULL, ##<< class data.SS
  data.logratios = NULL, ##<< class data.logratios
  output.dir = NULL, ##<< used in validation exercise only (to read in training set or save training set).
  verbose = TRUE, ##<< logical: print info?
  validation.list = NULL, ##<< NULL or of class validation.list (see \code{\link{RunMCMC}})
  do.country.specific.run = FALSE, ##<< logical, see \code{\link{RunMCMC}}
  do.SS.run.second.pass = FALSE, ##<< logical, see \code{\link{RunMCMC}}
  do.country.specific.targets.run = FALSE, ##<< logical, see \code{\link{RunMCMC}} # change JR, 20140301
  change.priors.to.zerolower = FALSE, ##<< logical, see \code{\link{RunMCMC}}
  names.sources = c("DHS", "MICS", "NS", "Other", "SS"), ##<< defines numeric IDs of data sources for contraceptive use.
  #, based on \code{unique(data$source.j)}
  names.sources.unmet = c("DHS", "Other") ##<< defines ordering of 
  ## data sources for unmet need.
  #, based on \code{unique(source.unmetj)}
) {
  if (do.SS.run.second.pass & is.null(data.SS))
    stop("data.SS cannot be NULL if this is the second pass of a run with SS data!")
  
  if (do.country.specific.targets.run) {
    #------------------------------------------------------------------------------------
    # 1. Extract information about last available year of estimates
    logratio.y.jn <- matrix(NA, 1, 3)
    c.select <- which(data.global$iso.c == gsub(" ", "", country.info$code.c[1]))
    logratio.y.jn[1, ] <- c(data.logratios$median.log.trad.noneed.c[c.select],
                            data.logratios$median.log.mod.noneed.c[c.select],
                            data.logratios$median.log.unmet.noneed.c[c.select])
    # Tau.logratios <- solve(data.logratios$Sigma.medianmin.33)
    Tau.logratios <- solve(data.logratios$Sigma.c33[c.select, , ]) # country-specific
    #------------------------------------------------------------------------------------
    # 2. Find all sorts of indices of the observations (year index, country, i)
    # round years at half-year
    round.years.j <- floor(data.logratios$year) + 0.5 # change JR, 20150301
    # for obs j, find the country c and year index t
    gett.j <- round.years.j  # so t refers to midpoint of calendar year
    J <- length(gett.j)
    C <- length(country.info$iso.c)
    
    # getc.J: which UNIQUE obs years are there in country c?
    N.unique.c <- rep(1, C) # change JR, 20150301
    gett.ci <- matrix(gett.j, C, 1) # change JR, 20150301
    geti.j <- getc.j <- rep(1, J)
    # For AR: find the indices of the countries with more than 1 observation:  
    getc.z <- seq(1, C)[N.unique.c>1]
    n.countriesmorethan1obs <- length(getc.z)
    #------------------------------------------------------------------------------------
    # 3. Regional information
    # Note: reg.c and subreg.c are categorical variables/factors, so as.numeric gives their level
    n.reg <- length(unique(country.info$reg.c))
    n.subreg <- length(unique(country.info$subreg.c))
    subreg.c <- country.info$subreg.c
    reg.subreg <- rep(NA, n.subreg)
    for (subreg in 1:n.subreg){
      c <- which.max(country.info$subreg.c==subreg)
      reg.subreg[subreg] <- country.info$reg.c[c]
    }
    crich.index <- seq(1, C)[country.info$dev.c=="Rich"] # dev.c is factor
    n.rich <- length(crich.index)
    cnotrich.index <- seq(1, C)[country.info$dev.c!="Rich"]
    n.notrich <- length(cnotrich.index)
    #----------------------------------------------------------------------------------
    ##details<< list.no.validation contains
    ##describe<< 
    list.no.validation <- list( 
      logratio.y.jn = logratio.y.jn, 
      Tau.logratios = Tau.logratios,
      pmid.for.unmet = 0.4, ##<< constant in model unmet
      getc.z = getc.z, ##<< indices of countries with more than 1 observation (for AR loop)
      n.countriesmorethan1obs = n.countriesmorethan1obs, ##<< count
      getc.j = getc.j, ##<< country index for each obs index j
      N.unique.c = c(N.unique.c), ##<< number of unique obs years per country
      gett.ci = gett.ci, ##<< year index for obs i in country c
      geti.j = geti.j, ##<< index of obs year in country c
      C = C, ##<< no of countries
      J = J, ##<< no of observations
      crich.index = crich.index,##<< indices of developED (rich) countries
      cnotrich.index = cnotrich.index, ##<< indices of developING (notrich) countries
      n.rich = n.rich, ##<< no of developed countries
      n.notrich = n.notrich,##<< no of developing countries
      reg.subreg = reg.subreg, ##<< region index for each subregion
      n.subreg = n.subreg,##<< no of subregions
      n.reg = n.reg, ##<< no of regions
      subreg.c = subreg.c##<< subregion index for each country
    )
    list.validation <- NULL
    ##end<<
  } else {
    # note: ordering of obs is not changed (else training wrong!)
    #------------------------------------------------------------------------------------
    # 0. If do.SS.run.first.pass, all obs of SS already removed from data in ReadDataAll.
    #    If do.SS.run.second.pass, remove all observations of SS prior to and including SS obs used to estimate SS bias for the JAGS model.
    if (do.SS.run.second.pass) {
      max.survey.year <- max(data$years.j[data$source.j != "SS"])
      diff.years.ss <- data$years.j[data$source.j == "SS"] - max.survey.year
      # get SS data point, either the closest SS obs prior to/in most recent survey year if available, else the closest one after that year
      if (any(diff.years.ss <= 0)) { # if SS data overlaps with non-SS data
        year.ss <- max(data$years.j[data$source.j == "SS" & data$years.j <= max.survey.year])
      } else {
        year.ss <- min(data$years.j[data$source.j == "SS"])
      }
      remove <- data$source.j == "SS" & data$years.j <= year.ss
      data <- data[!remove, ]
    }
    #------------------------------------------------------------------------------------
    # 1. Set min. observed trad/modern/no need/unmet to 1% and reduce other categories
    # print("Observed proportions less than 1% are set to 1%")
    props.tot.j <- data$props.tot.j
    props.modern.j <- data$props.modern.j
    props.trad.j <- data$props.trad.j
    props.unmet.j <- data$props.unmet.j
    
    select <- ifelse(props.tot.j < 0.01 & !is.na(props.modern.j) & !is.na(props.tot.j), T,F) # change JR, 20131120
    props.modern.j[select] <- 0.01
    props.trad.j[select] <- 0.01
    props.tot.j[select] <- 0.02
    select <- ifelse(props.tot.j < 0.01 & is.na(props.modern.j) & !is.na(props.tot.j), T,F) # change JR, 20131120
    props.tot.j[select] <- 0.01
    
    select <- ifelse(props.modern.j < 0.01 & !is.na(props.modern.j) & !is.na(props.tot.j), T,F) # change JR, 20131120
    props.tot.j[select] <- props.tot.j[select] + (0.01 - props.modern.j[select])
    props.modern.j[select] <- 0.01
    select <- ifelse(props.trad.j < 0.01 & !is.na(props.modern.j) & !is.na(props.tot.j), T,F) # change JR, 20131120
    props.tot.j[select] <- props.tot.j[select] + (0.01 - props.trad.j[select])
    props.trad.j[select] <- 0.01
    
    select <- ifelse(props.unmet.j < 0.01  & !is.na(props.unmet.j) & !is.na(props.tot.j), T,F) # change JR, 20131120
    props.tot.j[select] <- props.tot.j[select] - (0.01 - props.unmet.j[select])
    props.unmet.j[select] <- 0.01
    
    logratio.ymodern.j <- log(props.modern.j/(1-props.tot.j))
    logratio.ytrad.j <- log(props.trad.j/(1-props.tot.j))
    logit.ytot.j <- log(props.tot.j/(1-props.tot.j))
    logitratio.yunmet.j <- logit(props.unmet.j/(1-props.tot.j))
    y.modern.j <- props.modern.j # change JR, 20131120 # SS obs selected for in model using getj.training.modern.k
    se.modern.j <- rep(0.025, length(y.modern.j)) # change JR, 20140806
    #------------------------------------------------------------------------------------
    # 2. Find all sorts of indices of the observations (year index, country, i)
    # round years at half-year
    cat("Observation years are grouped into their respective calendar years.\n")
    round.years.j <- floor(data$years.j) + 0.5
    # for obs j, find the country c and year index t
    gett.j <- round.years.j  # so t refers to midpoint of calendar year
    J <- length(gett.j)
    C <- length(country.info$iso.c)
    
    # getc.J: which UNIQUE obs years are there in country c?
    # to get getc.j,  watch out with as.numeric(data$name.j), gives order alphabetically!
    # so just do a loop!
    N.unique.c <- rep(NA, C)
    # gett.ci: the SORTED indices of the obs years for i = 1, .., N.unique.c[c] (for AR(1))
    # Note: i does not relate to the observations anymore, just the unique years!
    gett.ci <- matrix(NA, C, max(country.info$N.c))
    for (c in 1:C) {
      select <- seq(1, J)[data$iso.j == country.info$iso.c[c]]
      N.unique.c[c] <- length(unique(gett.j[select]))
      gett.ci[c,1:N.unique.c[c]] <- sort(unique(gett.j[select]))  
    }
    geti.j <- getc.j <- rep(NA, J)
    for (c in 1:C){
      select <- seq(1, J)[data$iso.j == country.info$iso.c[c]]
      getc.j[select] <- c
      # to get i for obs j, find year and see which index it has...
      for (jc in 1:length(select)){
        geti.j[select[jc]] <- which.max(gett.j[select][jc]==sort(unique(gett.j[select])))
      }
    }
    # For AR: find the indices of the countries with more than 1 observation:  
    getc.z <- seq(1, C)[N.unique.c>1]
    n.countriesmorethan1obs <- length(getc.z)
    
    # For unmet: find # obs and indices to skip the missing observation j's
    N.unmet <- sum(!is.na(props.unmet.j))
    getj.unmet.k <- seq(1,J)[!is.na(props.unmet.j)]
    #------------------------------------------------------------------------------------
    # 3. Regional information
    # Note: reg.c and subreg.c are categorical variables/factors, so as.numeric gives their level
    n.reg <- length(unique(country.info$reg.c))
    n.subreg <- length(unique(country.info$subreg.c))
    subreg.c <- country.info$subreg.c
    reg.subreg <- rep(NA, n.subreg)
    for (subreg in 1:n.subreg){
      c <- which.max(country.info$subreg.c==subreg)
      reg.subreg[subreg] <- country.info$reg.c[c]
    }
    #country.info$name.c[country.info$dev.c=="Rich"]
    crich.index <- seq(1, C)[country.info$dev.c=="Rich"] # dev.c is factor
    n.rich <- length(crich.index)
    cnotrich.index <- seq(1, C)[country.info$dev.c!="Rich"]
    n.notrich <- length(cnotrich.index)
    #----------------------------------------------------------------------------------
    # 4. Data dummies
    source.ind.j <- ifelse(data$source.j==names.sources[1],1,
                           ifelse(data$source.j==names.sources[2],2,
                                  ifelse(data$source.j==names.sources[3],3,
                                         ifelse(data$source.j==names.sources[4],4,5)))) # change JR, 20131120
    # for unmet, only DHS versus non-DHS  
    source.ind.unmet.j <- ifelse(data$source.unmet.j==names.sources.unmet[1], 1, 2)
    
    # ind for SA, EMAL, HW, age refers to a 0,1,2,3 etc indicator, (O = not applicable)
    # note that within countries, only one multiplier is assigned
    # approach for those categories is to use "" for the first category
    # followed by a new level for each additional country
    # use as.numeric for factors, then 1 corresponds to first level 
    # (alphabetically ordered, thus "" first) 
    # levels(as.factor(c("", "   ")))
    # levels(as.factor(c("", "0")))
    
    # note: list of multipliers is verified using 
    # par.V <- GetParnamesV(winbugs.data = winbugs.data, name.short.j = MakeCountryNamesShort(data$name.j))
    # par.V$parnames.V.in.bugs
    # par.V$parnames.V.nice
    # data.frame(par.V$parnames.V.in.bugs,  par.V$parnames.V.nice)
    
    # note: each of these vectors are augmented with a (J+1)th element containing the baseline category
    # such that there are no errors when the baseline category is not present! # change JR, 2013111
    
    select <- !is.na(props.tot.j) # change JR, 20131120
    
    sa.factor <- as.factor(c(ifelse(data$poptype.j=="SA" & select, data$name.j, ""), "")) # change JR, 20131120
    sa.ind.j <- as.numeric(sa.factor)
    ncat.sa <- nlevels(sa.factor)
    
    emal.factor <- as.factor(c(ifelse( (data$poptype.j=="EM"|data$poptype.j=="AL") & select, data$name.j, ""), "")) # change JR, 20131120
    emal.ind.j <- as.numeric(emal.factor)
    ncat.emal <- nlevels(emal.factor)
    
    hw.factor <- as.factor(c(ifelse(data$poptype.j=="HW" & select, data$name.j, ""), "")) # change JR, 20131120
    hw.ind.j <- as.numeric(hw.factor)
    ncat.hw <- nlevels(hw.factor)
    
    age.factor <- as.factor(c(ifelse(data$age.cat.j=="?" & select, data$name.j, ""), "")) # change JR, 20131120
    age.ind.j <- as.numeric(age.factor)
    ncat.age <- nlevels(age.factor)
    
    posage.factor <- as.factor(c(ifelse(data$age.cat.j=="+" & select, data$name.j, ""), "")) # change JR, 20131120
    posage.ind.j <- as.numeric(posage.factor)
    ncat.posage <- nlevels(posage.factor)
    
    negage.factor <- as.factor(c(ifelse(data$age.cat.j=="-" & select, data$name.j, ""), "")) # change JR, 20131120
    negage.ind.j <- as.numeric(negage.factor)
    ncat.negage <- nlevels(negage.factor)
    
    ##<< Combine iso code/name with the explanation of the subgroup for geo and posbias 
    ## (such that different multipliers are added, e.g. if different geo regions used within one country)
    geo.factor <- as.factor(c(ifelse(data$geo.j!="" & select, 
                                     paste(data$name.j, data$geo.j, sep = ": "), ""), "")) # change JR, 20131120
    geo.ind.j <- as.numeric(geo.factor)
    ncat.geo <- nlevels(geo.factor)
    
    posbias.factor <- as.factor(c(ifelse(data$posbias.j!="" & select, 
                                         paste(data$name.j, data$posbias.j, sep = ": "), ""), "")) # change JR, 20131120
    posbias.ind.j <- as.numeric(posbias.factor)
    ncat.posbias <- nlevels(posbias.factor)
    
    # ind1 refers to a 0/1 indicator, (1 = Yes)
    # For biases:
    source.MICS.ind1.j <- ifelse(data$source.j == "MICS" & select, 1,0) # change JR, 20131120
    folk.ind1.j <- ifelse(data$folkbias.j != "" & select, 1,0) # change JR, 20131120
    mneg.ind1.j <- ifelse(data$mod.bias.j == "-" & select, 1,0) # change JR, 20131120
    mpos.ind1.j <- ifelse(data$mod.bias.j == "+" & select, 1,0) # change JR, 20131120
    
    #  if (return.Vinfo){
    #     return(list(posbias.factor = unique(posbias.factor),
    #                 geo.factor = unique(geo.factor),
    #                 negage.factor = unique(negage.factor), 
    #                 posage.factor = unique(posage.factor),  
    #                 age.factor = unique(age.factor),
    #                 hw.factor = unique(hw.factor),  
    #                 sa.factor = unique(sa.factor),
    #                 emal.factor = unique(emal.factor)
    #                 ))
    #   }
    ##details<< list.no.validation contains
    ##describe<< 
    list.no.validation <- list( 
      N.unmet = N.unmet, ##<<count unmet
      getj.unmet.k = getj.unmet.k, ##<< indices unmet
      pmid.for.unmet = 0.4, ##<< constant in model unmet
      logitratio.yunmet.j = logitratio.yunmet.j,  ##<< logit(unmet/1-total) observed
      ratios.trad.modern.jn = as.matrix(cbind(logratio.ytrad.j, logratio.ymodern.j)),##<< matrix, observed logit(trad/tot, mod/tot))
      logit.ytot.j = logit.ytot.j,##<< logit(tot/none) observed
      y.modern.j = y.modern.j, ##<< observations of modern from SS # change JR, 20131121
      se.modern.j = se.modern.j, ##<< SE of observations of modern from SS (set to 2.5 percent) # change JR, 20131120
      getc.z = getc.z, ##<< indices of countries with more than 1 observation (for AR loop)
      n.countriesmorethan1obs = n.countriesmorethan1obs, ##<< count
      getc.j = getc.j, ##<< country index for each obs index j
      N.unique.c = c(N.unique.c), ##<< number of unique obs years per country
      gett.ci = gett.ci, ##<< year index for obs i in country c
      geti.j = geti.j, ##<< index of obs year in country c
      C = C, ##<< no of countries
      J = J, ##<< no of observations (total non-missing observations on total contraceptive use)
      crich.index = crich.index,##<< indices of developED (rich) countries
      cnotrich.index = cnotrich.index, ##<< indices of developING (notrich) countries
      n.rich = n.rich, ##<< no of developed countries
      n.notrich = n.notrich,##<< no of developing countries
      reg.subreg = reg.subreg, ##<< region index for each subregion
      n.subreg = n.subreg,##<< no of subregions
      n.reg = n.reg, ##<< no of regions
      subreg.c = subreg.c,##<< subregion index for each country
      source.ind.unmet.j = source.ind.unmet.j, ##<< indicator for unmet source (1,2)
      source.ind.j = source.ind.j,##<< indicator for source (1,2,3,4)
      sa.ind.j = sa.ind.j, ##<< indicator for SA women (1 = NA, 2+ gives unique pertubation multiplier)
      ncat.sa = ncat.sa,##<< 1 + no of permutation parameters 
      posbias.ind.j = posbias.ind.j, ##<< indicator for pos bias (1 = NA, 2+ gives unique pertubation multiplier)
      ncat.posbias = ncat.posbias,##<< 1 + no of permutation parameters
      age.ind.j = age.ind.j, ##<< indicator for age (1 = NA, 2+ gives unique pertubation multiplier)
      ncat.age = ncat.age,##<<1 + no of permutation parameters
      hw.ind.j = hw.ind.j, ##<< indicator for HW (1 = NA, 2+ gives unique pertubation multiplier)
      ncat.hw = ncat.hw,##<<1 + no of permutation parameters
      emal.ind.j = emal.ind.j,##<<  indicator for EM/AL (1 = NA, 2+ gives unique pertubation multiplier)
      ncat.emal = ncat.emal,##<<1 + no of permutation parameters
      geo.ind.j = geo.ind.j, ##<< indicator for geo bias (1 = NA, 2+ gives unique pertubation multiplier)
      ncat.geo = ncat.geo,##<<1 + no of permutation parameters
      posage.ind.j = posage.ind.j,##<< indicator for pos age bias (1 = NA, 2+ gives unique pertubation multiplier)
      ncat.posage = ncat.posage,##<<1 + no of permutation parameters
      negage.ind.j = negage.ind.j, ##<< indicator for neg age bias (1 = NA, 2+ gives unique pertubation multiplier)
      ncat.negage = ncat.negage,##<<1 + no of permutation parameters  
      source.MICS.ind1.j = source.MICS.ind1.j, ##<< indicator for MICS (1 = yes, 0 = no)
      folk.ind1.j = folk.ind1.j,##<<  indicator for folk(1 = yes, 0 = no)
      mpos.ind1.j = mpos.ind1.j, ##<<  indicator for modern[+] (1 = yes, 0 = no)
      mneg.ind1.j = mneg.ind1.j ##<<  indicator for modern[-] (1 = yes, 0 = no)
    )
    if (do.SS.run.second.pass) { # change JR, 20140414
      list.no.validation <- c(list.no.validation, 
                              list(bias.modern = data.SS$bias.modern ##<< calculated SS modern bias
                              ))                               
    }
    ##end<<
    
    if (!is.null(validation.list)){
      ##details<< If \code{!is.null(validation.list)}, \code{getj.training.k} is constructed
      ## using \code{\link{GetTraining}}
      if (validation.list$generate.new.set){
        getj.training.k <- GetTraining(data, winbugs.data = NULL, 
                                       at.random = validation.list$at.random, 
                                       at.end = validation.list$at.end, 
                                       year.cutoff = validation.list$year.cutoff,
                                       exclude.unmet.only = validation.list$exclude.unmet.only, 
                                       seed = validation.list$seed)
        save(getj.training.k, file = file.path(output.dir, "getj.training.k.rda")) # change JR, 20140418
      } else { # not used!
        load(file = file.path(output.dir, "getj.training.k.rda")) # change JR, 20140418
      }
      # when leaving out obs at random, all tot are included in training
      # when leaving out obs at end, some tot might be excluded (but tot in training is not zero)
      # and test set should not include tot obs! 
      
      # setdiff means first set is the baseline, remove any elements that are in the 2nd set
      getj.test.unmet.k <- setdiff(seq(1, J)[!is.na(props.unmet.j)], getj.training.k)
      getj.training.unmet.k <- setdiff(seq(1, J)[!is.na(props.unmet.j)],getj.test.unmet.k)
      n.training.unmet <- length(getj.training.unmet.k)
      n.test.unmet <- length(getj.test.unmet.k)
      if (validation.list$exclude.unmet.only){
        # include all total/breakdown in training!!!
        getj.training.tot.k <- seq(1, J)[is.na(logratio.ymodern.j) & !is.na(props.tot.j)] # change JR, 20131120
        getj.training.k <- seq(1, J)[!is.na(logratio.ymodern.j) & !is.na(props.tot.j)] # change JR, 20131120
        n.training.tot <- length(getj.training.tot.k)
        n.training.breakdown <- length(getj.training.k)
        getj.training.modern.k <- seq(1, J)[data$source.j == "SS"] # change JR, 20140612
        n.training.modern <- length(getj.training.modern.k) # change JR, 20131120
        list.validation <- list(
          getj.training.k  = getj.training.k, 
          n.training.breakdown = n.training.breakdown,
          getj.training.tot.k  = getj.training.tot.k,
          n.training.tot = n.training.tot,
          getj.training.modern.k = getj.training.modern.k, # change JR, 20131120
          n.training.modern = n.training.modern, # change JR, 20131120
          getj.training.unmet.k = getj.training.unmet.k, 
          getj.test.unmet.k = getj.test.unmet.k,
          n.training.unmet = n.training.unmet,
          n.test.unmet = n.test.unmet
        )
      } else { # at end or at random
        getj.test.k <- setdiff(seq(1, J)[!is.na(logratio.ymodern.j) & !is.na(props.tot.j)], getj.training.k) # change JR, 20131120
        n.training.breakdown <- length(getj.training.k)
        n.test.breakdown <- length(getj.test.k)
        getj.training.tot.k <- seq(1, J)[is.na(logratio.ymodern.j) & !is.na(props.tot.j)] # change JR, 20131120
        n.training.tot <- length(getj.training.tot.k)
        getj.training.modern.k <- seq(1, J)[data$source.j == "SS"] # change JR, 20140612
        n.training.modern <- length(getj.training.modern.k) # change JR, 20131120
        ##details<< Then \code{list.validation} is given by:
        ##describe<<
        list.validation <- list(
          getj.training.k  = getj.training.k, ##<< indices j for training obs k=1,...
          getj.test.k = getj.test.k,##<< indices (left-out if unmet left-out only)
          n.training.breakdown = n.training.breakdown,##<< count of training obs
          n.test.breakdown = n.test.breakdown,##<< count of test obs (left-out if unmet left-out only)
          getj.training.tot.k  = getj.training.tot.k,##<< indices
          n.training.tot = n.training.tot,##<< count of training obs
          getj.training.modern.k = getj.training.modern.k, # change JR, 20131120
          n.training.modern = n.training.modern, # change JR, 20131120
          getj.training.unmet.k  = getj.training.unmet.k, ##<< indices
          getj.test.unmet.k = getj.test.unmet.k,##<< indices
          n.training.unmet = n.training.unmet,##<< count of training obs
          n.test.unmet = n.test.unmet##<< count of test obs
        ) 
        ##end<<
        ## or smaller list if only unmet data were left out.
      } # end else for unmet only or not
    } else { #all, no validation
      getj.training.tot.k <- seq(1, J)[is.na(logratio.ymodern.j) & !is.na(props.tot.j)] # change JR, 20131120 
      getj.training.k <- seq(1, J)[!is.na(logratio.ymodern.j) & !is.na(props.tot.j)] # change JR, 20131120
      getj.test.k <- NULL
      getj.training.modern.k <- seq(1, J)[data$source.j == "SS"] # change JR, 20140612
      getj.training.unmet.k <- seq(1, J)[!is.na(props.unmet.j)]
      getj.test.unmet.k <- NULL
      n.training.breakdown <- length(getj.training.k)
      n.training.tot <- length(getj.training.tot.k)
      n.training.modern <- length(getj.training.modern.k) # change JR, 20131120
      n.training.unmet <- length(getj.training.unmet.k)
      ##details<< If it is not a validation exercise, list.validation contains
      ##describe<<
      list.validation <- list(
        getj.training.k  = getj.training.k, ##<< indices
        n.training.breakdown = n.training.breakdown, ##<< count
        getj.training.tot.k  = getj.training.tot.k, ##<<indices
        n.training.tot = n.training.tot, ##<< count
        getj.training.modern.k = getj.training.modern.k, ##<< indices # change JR, 20131120
        n.training.modern = n.training.modern, ##<< count # change JR, 20131120
        getj.training.unmet.k  = getj.training.unmet.k, ##<< indices 
        n.training.unmet = n.training.unmet ##<< count
      )
      ##end<<
    }
  }
  priorspecs <- GetBugsPriorSpecs(change.priors.to.zerolower = change.priors.to.zerolower,
                                  do.country.specific.run = do.country.specific.run, # change JR, 20131104
                                  do.country.specific.targets.run = do.country.specific.targets.run, # change JR, 20131104
                                  name.reg = as.character(country.info$namereg.c)[1], # change JR, 20150301
                                  name.subreg = as.character(country.info$namesubreg.c)[1], # change JR, 20131104
                                  iso.country.select = country.info$iso.country.select[1], # change JR, 20140404
                                  data.global = data.global) # change JR, 20131104
  ##value<< One combined list that includes elements from 
  #  ##describe<<
  winbugs.data <- c(list.no.validation, ##<< See details.
                    list.validation, ##<< See details.
                    priorspecs) ##<< Prior specs from \code{GetBugsPriorSpecs}.
  #  ##end<<
  return(winbugs.data) 
}

#----------------------------------------------------------------------------------
GetBugsPriorSpecs <- function( #Set priors parameters
  ### Set prior parameters for CP model. The names for the prior parameters that are not explained below
  ### follow from
  ### the names used in R for the model parameters and
  ### the specification of the prior distributions 
  ### (see"Prior distributions" and the table with names used in \code{R} in the 
  ### web appendix of Alkema et al).
  change.priors.to.zerolower = FALSE, ##<< logical indicating if gammas are used for kappa.c's
  do.country.specific.run = FALSE, ##<< logical indicating if run is country-specific
  do.country.specific.targets.run = FALSE, ##<< logical indicating if run is country-specific for targets # change JR, 20140301
  name.reg = NULL, ##<< (For country-specific run for targets) character giving the name of the region
  ## to which country belongs to # change JR, 20150301
  name.subreg = NULL, ##<< (For country-specific run) character giving the name of the subregion
  ## to which country belongs to # change JR, 20131104
  iso.country.select = NULL, ##<< (For subpopulation-specific run) character giving the 3-character ISO country code of the country the
  ## subpopulation belongs to # change JR, 20140404
  data.global = NULL, ##<< (For country-specific run) object from \code{\link{SummariseGlobalRun}} # change JR, 20131104
  rho.max = 1, ##<< upper bound rho for AR with total/ratio
  sigma.ar.max =1, ##<< upper bound sd sigma for AR with total/ratio
  rho.max.unmet = 1, ##<< upper bound rho for AR with unmet
  sigma.ar.max.unmet = 1, ##<< upper bound sd sigma for AR with unmet
  sigma2.sourcetot0 = 0.0225, 
  R = as.matrix(cbind(c(0.1,0), c(0, 0.1))),##<< prior precision matrix for logit(trad/total, mod/total)
  nu0 = 10, ##<< prior sample size for kappa.c's
  # upto June 6, 2012
  #   sigma2.lpc0 = 0.43,##<<
  #   sigma2.lrc0 = 0.87, ##<<
  #   sigma2.wc0 = 0.17,##<<
  #   sigma2.Rwc0 =  0.18,##<<
  #   sigma2.Tc0 = 253.87,##<<
  #   sigma2.RTc0 =  151.66,##<<
  #   sigma2.earlierTc0 = 339.94, ##<<
  #   sigma2.unmetc0 = 0.4^2,##<<
  
  # based on run20120608_noAR
  sigma2.lpc0 = 0.63,##<<
  sigma2.lrc0 = 1.97, ##<<
  sigma2.wc0 =  0.35,##<<
  sigma2.Rwc0 =   0.65,##<<
  sigma2.Tc0 = 300,##<<
  sigma2.RTc0 =  197,##<<
  sigma2.earlierTc0 = 263, ##<<
  sigma2.unmetc0 = 0.08,##<<
  
  mean.Tworld = 1980,##<<
  mean.RTworld = 1980,##<<
  mean.Tearlier = 1920,##<<
  sigmaTregsubreg.upper = 80, ##<<,
  tau0.T = 1/50^2, ##<< corresponding to SD of 50 years
  sigmawregsubreg.upper = 3, ##<< (refers to logit scale)
  a0.unmet = -0.38, ##<< From LS-fit
  b0.unmet = 0.12, ##<< From LS-fit
  tau.a0 = 1, ##<< 
  tau.b0 = 1 ##<<             
){
  if (!do.country.specific.run & !do.country.specific.targets.run) { # change JR, 20150301
    ##details<< \code{prior.list1} has elements that do not relate to the priors for the kappa.c's.
    prior.list1 <- list(
      rho.max = rho.max, 
      sigma.ar.max = sigma.ar.max, 
      rho.max.unmet = rho.max.unmet,  
      sigma.ar.max.unmet = sigma.ar.max.unmet, 
      halfsigma2.sourcetot0 = 0.5*sigma2.sourcetot0, 
      R = R, 
      mean.Tworld = mean.Tworld, 
      mean.RTworld = mean.RTworld,
      mean.Tearlier = mean.Tearlier,
      sigmaTregsubreg.upper = sigmaTregsubreg.upper, 
      tau0.T = tau0.T, 
      sigmawregsubreg.upper = sigmawregsubreg.upper, 
      a0.unmet = a0.unmet, 
      b0.unmet = b0.unmet, 
      tau.a0 = tau.a0, 
      tau.b0 = tau.b0)
  } else { # change JR, 20131104
    mcmc.post <- data.global$mcmc.post
    prior.list1 <- list(
      rho.tot0 = mcmc.post$rho.tot,
      sigma.tot0 = mcmc.post$sigma.tot,
      rho.rat0 = mcmc.post$rho.rat,
      sigma.rat0 = mcmc.post$sigma.rat,      
      lp.world0 = mcmc.post$lp.world,
      lr.world0 = mcmc.post$lr.world,
      Tearlier0 = mcmc.post$Tearlier,
      rho.unmet0 = mcmc.post$rho.unmet,
      sigma.ar.unmet0 = mcmc.post$sigma.ar.unmet,
      a.unmet0 = mcmc.post$a.unmet,
      b.unmet0 = mcmc.post$b.unmet,
      c.unmet0 = mcmc.post$c.unmet) # change JR, 20150301
    if (!do.country.specific.targets.run) { # change JR, 20150301
      prior.list1 <- c(prior.list1, list(
        v.mics0 = mcmc.post$v.mics,
        v.mneg0 = mcmc.post$v.mneg,
        v.folk0 = mcmc.post$v.folk, 
        v.mpos0 = mcmc.post$v.mpos,
        sigma.pos0 = mcmc.post$sigma.pos,
        mu.pos.m0 = c(mcmc.post[['mu.pos.m[1]']], mcmc.post[['mu.pos.m[2]']]),
        tau.sourcetot0 = 1/(mcmc.post$sigma.sourcetot^2),
        T1.source.s0 = c(mcmc.post[['T1.source.s[1]']], mcmc.post[['T1.source.s[2]']],
                         mcmc.post[['T1.source.s[3]']], mcmc.post[['T1.source.s[4]']]),
        T2.source.s0 = c(mcmc.post[['T2.source.s[1]']], mcmc.post[['T2.source.s[2]']],
                         mcmc.post[['T2.source.s[3]']], mcmc.post[['T2.source.s[4]']]),
        T12.source.s0 = c(mcmc.post[['T12.source.s[1]']], mcmc.post[['T12.source.s[2]']],
                          mcmc.post[['T12.source.s[3]']], mcmc.post[['T12.source.s[4]']]),
        sigma.unmet.other0 = mcmc.post$sigma.unmet.other,
        sigma.unmet.dhs0 = mcmc.post$sigma.unmet.dhs,
        sigma.unmetworld0 = mcmc.post$sigma.unmetworld,
        sigma.geo.m0 = c(mcmc.post[['sigma.geo.m[1]']], mcmc.post[['sigma.geo.m[2]']])))
      if (is.null(iso.country.select)) { # for country-specific run
        subreg.global <- which(data.global$name.subreg == name.subreg)
        prior.list1 <- c(prior.list1, list(
          unmet.subreg0 = mcmc.post[[paste0("unmet.subreg[", subreg.global, "]")]],
          w.subreg0 = mcmc.post[[paste0("w.subreg[", subreg.global, "]")]],
          T.subreg0 = mcmc.post[[paste0("T.subreg[", subreg.global, "]")]],
          Rw.subreg0 = mcmc.post[[paste0("Rw.subreg[", subreg.global, "]")]],
          RT.subreg0 = mcmc.post[[paste0("RT.subreg[", subreg.global, "]")]]
        ))
      } else { # for subpopulation-specific run
        country.global <- which(data.global$iso.c == iso.country.select) 
        prior.list1 <- c(prior.list1, list(
          unmet.subreg0 = mcmc.post[[paste0("unmet.intercept.c[", country.global, "]")]],
          w.subreg0 = LogitMinMax(mcmc.post[[paste0("omega.c[", country.global, "]")]], 0.01, 0.5),
          T.subreg0 = mcmc.post[[paste0("T.c[", country.global, "]")]],
          Rw.subreg0 = LogitMinMax(mcmc.post[[paste0("Romega.c[", country.global, "]")]], 0.01, 0.5),
          RT.subreg0 = mcmc.post[[paste0("RT.c[", country.global, "]")]]
        ))
      }
      if (change.priors.to.zerolower) {
        prior.list1 <- c(prior.list1, list(
          sigma.lpc0 = mcmc.post$sigma.lpc,
          sigma.lrc0 = mcmc.post$sigma.lrc,
          sigma.wc0 = mcmc.post$sigma.wc,
          sigma.Rwc0 = mcmc.post$sigma.Rwc,
          sigma.Tc0 = mcmc.post$sigma.Tc,
          sigma.RTc0 = mcmc.post$sigma.RTc,
          sigma.earlierTc0 = mcmc.post$sigma.earlierTc,
          sigma.unmetc0 = mcmc.post$sigma.unmetc))
      } else {
        prior.list1 <- c(prior.list1, list(
          tau.unmetc0 = 1/(mcmc.post$sigma.unmetc^2),
          tau.earlierTc0 = 1/(mcmc.post$sigma.earlierTc^2),
          tau.Tc0 = 1/(mcmc.post$sigma.Tc^2),
          tau.lpc0 = 1/(mcmc.post$sigma.lpc^2),
          tau.lrc0 = 1/(mcmc.post$sigma.lrc^2),
          tau.wc0 = 1/(mcmc.post$sigma.wc^2),
          tau.Rwc0 = 1/(mcmc.post$sigma.Rwc^2),
          tau.RTc0 = 1/(mcmc.post$sigma.RTc^2)))
      }
    } else { # change JR, 20150301
      reg.global <- which(data.global$name.reg == name.reg)
      prior.list1 <- c(prior.list1, list(
        sigma.unmetworld0 = mcmc.post[["sigma.unmetworld"]],
        w.reg0 = mcmc.post[[paste0("w.reg[", reg.global, "]")]],
        T.reg0 = mcmc.post[[paste0("T.reg[", reg.global, "]")]],
        Rw.reg0 = mcmc.post[[paste0("Rw.reg[", reg.global, "]")]],
        RT.reg0 = mcmc.post[[paste0("RT.reg[", reg.global, "]")]],
        sigma.wsubreg0 = mcmc.post[["sigma.wsubreg"]],
        sigma.Rwsubreg0 = mcmc.post[["sigma.Rwsubreg"]],
        sigma.RTsubreg0 = mcmc.post[["sigma.RTsubreg"]],
        sigma.Tsubreg0 = mcmc.post[["sigma.Tsubreg"]]
      ))
    }
  }
  
  if (!do.country.specific.targets.run | (do.country.specific.targets.run & 
                                            do.country.specific.targets.run)) {
    if (!change.priors.to.zerolower) {
      halfnu0_rich =2/2
      halfnu0_poor =8/2
      # tau.earlierTc ~ dgamma(halfnu0_rich,halfnu0_rich_sigma2.earlierTc0)
      # tau.Tc ~ dgamma(halfnu0_poor,halfnu0_poor_sigma2.Tc0)
      
      ##details<<
      ## If \code{!change.priors.to.zerolower}, \code{prior.list2} contains prior gamma parameters,
      ## in the form of \code{halfnu0sigma2...}.
      prior.list2 <- list(
        halfnu0 = nu0/2, 
        halfnu0_rich = halfnu0_rich,
        halfnu0_poor = halfnu0_poor,
        halfnu0_rich_sigma2.earlierTc0 = halfnu0_rich*sigma2.earlierTc0, 
        halfnu0_poor_sigma2.Tc0 =  halfnu0_poor*sigma2.Tc0,
        #halfnu0sigma2.earlierTc0 = nu0/2*sigma2.earlierTc0, 
        #halfnu0sigma2.Tc0 =  nu0/2*sigma2.Tc0,
        halfnu0sigma2.lpc0 = nu0/2*sigma2.lpc0,
        halfnu0sigma2.lrc0 =  nu0/2*sigma2.lrc0, 
        halfnu0sigma2.wc0 = nu0/2*sigma2.wc0,
        halfnu0sigma2.Rwc0 =  nu0/2*sigma2.Rwc0,
        halfnu0sigma2.RTc0 = nu0/2*sigma2.RTc0,
        halfsigma2.unmetc0 = 0.5*sigma2.unmetc0
      )
    }
  } else {
    prior.list2 <- NULL
  }
  
  return(c(prior.list1, prior.list2))
  ### List with \code{priorlist1} and \code{priorlist2}, described above.
}
#----------------------------------------------------------------------------------
GetParNames <- function(# Get list of parnames
  ### Get list of parnames
  winbugs.data, ##<< Object from \code{\link{winbugs.data}}
  validation.list = NULL, ##<< Validation parameters are included if this is of object validation.list.
  do.country.specific.run = FALSE, ##<< Logical: do country-specific run?
  do.country.specific.targets.run = FALSE ##<< Logical: do country-specific targets run? # change JR, 20150301
){
  JAGS = TRUE 
  # Note: I used to have the option to run it in Bugs (and not specify the elements individually),
  # but that option is no longer needed
  # Note that for running in JAGS, 
  # you get an error when simply including parameter vectors with NAs.
  # So just specify all of the NAs ones individually (does mean they allllllll show up in Winbugs...)
  
  parnames.c <- c("omega.c", "T.c", "pmax.c", "Romega.c" ,"RT.c", "Rmax.c", "unmet.intercept.c")
  parnames.V <- parnames.h <- parnames.subreg <- parnames.reg <- 
    T1.source.s <- T2.source.s <- T12.source.s <- 
    parnames.ss <- NULL # change JR, 20140612
  
  if (is.null(validation.list)){
    parnames.validation <- NULL
  } else {
    if (validation.list$exclude.unmet.only){
      #parnames.validation <- c("pred.logitratio.yunmet.j", "q.unmet.j")
      #if (JAGS){
      parnames.validation <- c(
        paste0("pred.logitratio.yunmet.j[", winbugs.data$getj.test.unmet.k, "]"),
        paste0("q.unmet.j[", winbugs.data$getj.test.unmet.k, "]"),
        # save modern and trad as well to get the sampled total
        # (maybe more efficient to do in Bugs)
        paste0("pred.logratio.ytrad.j[", winbugs.data$getj.test.unmet.k, "]"),
        paste0("pred.logratio.ymodern.j[", winbugs.data$getj.test.unmet.k, "]")
      )
      #} # end JAGS
    }  else { # validation for unmet and modern/trad
      #parnames.validation <- c("pred.logratio.ytrad.j", "pred.logratio.ymodern.j", "pred.logit.ytotal.j",
      #                            "pred.logitratio.yunmet.j")
      #if (JAGS){
      parnames.validation <- c(
        paste0("pred.logratio.ytrad.j[", 
               winbugs.data$getj.test.k[1:winbugs.data$n.test.breakdown], "]"),
        #        paste0("q.trad.j[", winbugs.data$getj.test.k[1:winbugs.data$n.test.breakdown], "]"),
        paste0("pred.logratio.ymodern.j[", 
               winbugs.data$getj.test.k[1:winbugs.data$n.test.breakdown], "]"),
        #        paste0("q.modern.j[", winbugs.data$getj.test.k[1:winbugs.data$n.test.breakdown], "]"),
        paste0("pred.logitratio.yunmet.j[", winbugs.data$getj.test.unmet.k, "]")
        #        paste0("q.unmet.j[", winbugs.data$getj.test.unmet.k, "]")
      )    
      #} # end JAGS
    } # end validation unmet/trad/modern
  } # end validation 
  parnames.theta <- "theta.ci"
  parnames.eps <- "eps.ci"
  parnames.eta <- "eta.ci"
  parnames.eps <- parnames.theta <- parnames.eta <- NULL
  for (c in 1:winbugs.data$C){
    for (i in 1:winbugs.data$N.unique.c[c]){
      parnames.theta <- c(parnames.theta, paste0("theta.ci[", c, ",", i, "]"))
      parnames.eps <- c(parnames.eps, paste0("eps.ci[", c, ",", i, "]"))
      parnames.eta <- c(parnames.eta, paste0("eta.ci[", c, ",", i, "]"))
    }
  }
  if (!do.country.specific.targets.run) # change JR, 20150619
    parnames.V <- c("V.geo.12i", "V.age.12i", "V.hw.12i", "V.emal.12i",
                    "V.sa.12i", "V.posbias.12i", "V.posage.12i", "V.negage.12i")
  if (!do.country.specific.run & !do.country.specific.targets.run) { # for global run # change JR, 20150301
    # Note that there's a function for parnames.V for plotting
    parnames.h <- c(parnames.h,
                    "sigma.sourcetot", "sigma.unmet.dhs", "sigma.unmet.other",
                    "sigma.ar.unmet", "rho.unmet", 
                    "a.unmet", "b.unmet", "c.unmet",
                    "sigma.unmetc", 
                    "lp.world", "lr.world", "sigma.lrc", "sigma.lpc",  
                    "rho.tot", "sigma.tot", "rho.rat", "sigma.rat",
                    "sigma.geo.m[1]", "sigma.geo.m[2]", 
                    "sigma.pos",
                    "mu.pos.m[1]", "mu.pos.m[2]",
                    "v.mics", "v.folk", "v.mneg", "v.mpos",
                    "Tearlier", "sigma.earlierTc",
                    "sigma.wc", "sigma.Rwc", "sigma.Tc", "sigma.RTc",
                    # the following parameters are not found in WriteCountryModel at all
                    "sigma.unmetworld", 
                    "w.world", "Rw.world", "T.world", "RT.world",
                    "sigma.wreg", "sigma.Rwreg", "sigma.Treg", "sigma.RTreg",
                    "sigma.wsubreg", "sigma.Rwsubreg", "sigma.Tsubreg", "sigma.RTsubreg")
    parnames.reg <- c(parnames.reg, "w.reg","Rw.reg","T.reg","RT.reg")
    parnames.subreg <- c(parnames.subreg, "w.subreg","Rw.subreg","T.subreg","RT.subreg", "unmet.subreg")
    T1.source.s <- c(T1.source.s, "T1.source.s")
    T2.source.s <- c(T2.source.s, "T2.source.s")
    T12.source.s <- c(T12.source.s, "T12.source.s")
  } else if (do.country.specific.run & do.country.specific.targets.run) { # for country-specific targets run # change JR, 20150301
    parnames.h <- c(parnames.h,
                    "sigma.unmetc", "sigma.lrc", "sigma.lpc", "sigma.earlierTc",
                    "sigma.wc", "sigma.Rwc", "sigma.Tc", "sigma.RTc")
    parnames.subreg <- c(parnames.subreg, "w.subreg","Rw.subreg","T.subreg","RT.subreg", "unmet.subreg")
  }
  
  ##value<< List with 
  parnames <- list(parnames.V = parnames.V, ##<< Multipliers, e.g. V.geo.12i etc (no indices)
                   parnames.reg = parnames.reg, ##<< Regional mean parameters, e.g. w.reg (no indices)
                   parnames.subreg = parnames.subreg, ##<< Subregional mean parameters, e.g. w.subreg (no indices)
                   T1.source.s = T1.source.s, ##<< T11 of precision matrix trad/tot, modern/tot
                   T2.source.s = T2.source.s, ##<< T22 of precision matrix trad/tot, modern/tot
                   T12.source.s = T12.source.s, ##<< T12 of precision matrix trad/tot, modern/tot
                   parnames.h = parnames.h, ##<< Hyper parameters
                   parnames.c = parnames.c, ##<< Country parameters, e.g. omega.c (no indices)
                   parnames.theta = parnames.theta,##<< AR for unmet, e.g. theta.ci[193,3] (with indices)
                   parnames.eta = parnames.eta,##<< AR for ration (with indices)
                   parnames.eps = parnames.eps, ##<< AR for total (with indices)
                   parnames.ss = parnames.ss, ##<< parameters related to service statistics
                   parnames.validation = parnames.validation##<< NULL if no validation or not\code{JAGS}, 
                   ##else the q.j[j]'s and pred.logratio.ytrad.j[j]'s (with indices)
  ) 
  return(parnames)
}

#----------------------------------------------------------------------------------
InternalGetParnamesV <- function( # Get nice parameter names for the multipliers
  ### Get nice parameter names for the multipliers
  winbugs.data, ##<< Object
  name.short.j ##<< country names, e.g. data$name.j
){
  parnames.V.in.bugs <- NULL
  parnames.V.nice <- NULL
  if (winbugs.data$ncat.geo > 1) { # change JR, 20131104
    catindex.temp <- GetCatIndex(ncat = winbugs.data$ncat.geo, ind.j = winbugs.data$geo.ind.j, 
                                 name.short.j = name.short.j) # change JR, 20131112
    for (i in 2:winbugs.data$ncat.geo){       
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.geo.12i[1,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$geo.ind.j==i][1], " geo, trad",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$geo.ind.j==i),")"))
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.geo.12i[2,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$geo.ind.j==i][1], " geo, mod",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$geo.ind.j==i),")"))
    }
  }
  if (winbugs.data$ncat.age > 1) { # change JR, 20131104
    catindex.temp <- GetCatIndex(ncat = winbugs.data$ncat.age, ind.j = winbugs.data$age.ind.j, 
                                 name.short.j = name.short.j) # change JR, 20131112
    for (i in 2:winbugs.data$ncat.age){ 
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.age.12i[1,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$age.ind.j==i][1], " age, trad",
                                                   catindex.temp[i], # change JR, 20131112  
                                                   " (", sum(winbugs.data$age.ind.j==i),")"))
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.age.12i[2,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$age.ind.j==i][1], " age, mod",
                                                   catindex.temp[i], # change JR, 20131112  
                                                   " (", sum(winbugs.data$age.ind.j==i),")"))
    }
  }
  if (winbugs.data$ncat.hw > 1) { # change JR, 20131104
    catindex.temp <- GetCatIndex(ncat = winbugs.data$ncat.hw, ind.j = winbugs.data$hw.ind.j, 
                                 name.short.j = name.short.j) # change JR, 20131112
    for (i in 2:winbugs.data$ncat.hw){ 
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.hw.12i[1,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$hw.ind.j==i][1], " HW, trad",
                                                   catindex.temp[i], # change JR, 20131112  
                                                   " (", sum(winbugs.data$hw.ind.j==i),")"))
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.hw.12i[2,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$hw.ind.j==i][1], " HW, mod",
                                                   catindex.temp[i], # change JR, 20131112  
                                                   " (", sum(winbugs.data$hw.ind.j==i),")"))
    }
  }
  if (winbugs.data$ncat.emal > 1) { # change JR, 20131104
    catindex.temp <- GetCatIndex(ncat = winbugs.data$ncat.emal, ind.j = winbugs.data$emal.ind.j, 
                                 name.short.j = name.short.j) # change JR, 20131112
    for (i in 2:winbugs.data$ncat.emal){ 
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.emal.12i[1,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$emal.ind.j==i][1], " EM/AL, trad",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$emal.ind.j==i),")"))
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.emal.12i[2,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$emal.ind.j==i][1], " EM/AL, mod",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$emal.ind.j==i),")"))
    }
  }
  if (winbugs.data$ncat.sa > 1) { # change JR, 20131104
    catindex.temp <- GetCatIndex(ncat = winbugs.data$ncat.sa, ind.j = winbugs.data$sa.ind.j, 
                                 name.short.j = name.short.j) # change JR, 20131112
    for (i in 2:winbugs.data$ncat.sa){ 
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.sa.12i[1,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$sa.ind.j==i][1], " SA, trad",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$sa.ind.j==i),")"))
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.sa.12i[2,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$sa.ind.j==i][1], " SA, mod",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$sa.ind.j==i),")"))
    }
  }
  if (winbugs.data$ncat.posbias > 1) { # change JR, 20131104
    catindex.temp <- GetCatIndex(ncat = winbugs.data$ncat.posbias, ind.j = winbugs.data$posbias.ind.j, 
                                 name.short.j = name.short.j) # change JR, 20131112
    for (i in 2:winbugs.data$ncat.posbias){ 
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.posbias.12i[1,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$posbias.ind.j==i][1], " +, trad",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$posbias.ind.j==i),")"))
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.posbias.12i[2,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$posbias.ind.j==i][1], " +, mod",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$posbias.ind.j==i),")"))
    }
  }
  if (winbugs.data$ncat.posage > 1) { # change JR, 20131104
    catindex.temp <- GetCatIndex(ncat = winbugs.data$ncat.posage, ind.j = winbugs.data$posage.ind.j, 
                                 name.short.j = name.short.j) # change JR, 20131112
    for (i in 2:winbugs.data$ncat.posage){ 
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.posage.12i[1,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$posage.ind.j==i][1], " +Age, trad",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$posage.ind.j==i),")"))
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.posage.12i[2,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$posage.ind.j==i][1], " +Age, mod",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$posage.ind.j==i),")"))
    }
  }
  if (winbugs.data$ncat.negage > 1) { # change JR, 20131104
    catindex.temp <- GetCatIndex(ncat = winbugs.data$ncat.negage, ind.j = winbugs.data$negage.ind.j, 
                                 name.short.j = name.short.j) # change JR, 20131112
    for (i in 2:winbugs.data$ncat.negage){ 
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.negage.12i[1,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$negage.ind.j==i][1], " -Age, trad",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$negage.ind.j==i),")"))
      parnames.V.in.bugs <- c(parnames.V.in.bugs, paste0("V.negage.12i[2,", i, "]"))
      parnames.V.nice <- c(parnames.V.nice, paste0(name.short.j[winbugs.data$negage.ind.j==i][1], " -Age, mod",
                                                   catindex.temp[i], # change JR, 20131112 
                                                   " (", sum(winbugs.data$negage.ind.j==i),")"))
    }
  }
  ##value<< 
  return(list(parnames.V.in.bugs = parnames.V.in.bugs, ##<< vector with parnames used in Bugs
              parnames.V.nice = parnames.V.nice ##<< vector with longer description 
              ## for each parname, including country, pop group, mod/trad and number of observations
  ))
}

GetCatIndex <- function(# Get category index
  ### Internal function to get category index
  ncat,
  ind.j,
  name.short.j
) {
  catindex <- NULL
  name.previous <- " "
  for (i in 2:ncat){
    name.current <- name.short.j[ind.j==i][1]
    if (name.current == name.previous) {
      cat.by.country <- cat.by.country + 1
    } else {
      cat.by.country <- 1
    }
    catindex <- c(catindex, cat.by.country)
    name.previous <- name.current
  }
  catindex <- ifelse(diff(c(catindex, 1)) == 0, "", paste0(" ", catindex))
  # add dummy first element
  catindex <- c(0, catindex)
  return(catindex)
}
#-----------------------------------------------------------------
# The End!
