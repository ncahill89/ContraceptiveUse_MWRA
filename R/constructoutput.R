#--------------------------------------------------------------------------
# F_constructoutput.R
# Leontine Alkema, 2011
#--------------------------------------------------------------------------
ConstructOutput <- function(# Construct output for MCMC run 
  ### Construct output for MCMC run: country trajectories, and results for countries and UNDP aggregates.
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored, and new objects are added
  ## (if NULL, it's output/run.name, default from \code{runMCMC}).
  MWRA.csv = NULL,  ##<< If \code{NULL},
  ## estimates of the number of MWRA that are included in package are used. 
  ## To use an alternative csv file, use \\
  ## \code{MWRA.csv = .../Number-of-women-married-in union_15-49.csv}, 
  ## where ... refers to the file path where file is located.
  do.SS.run.first.pass = FALSE, ##<< do first pass of run with SS data? # change JR, 20140414
  seedAR = 1234, ##<< Seed for sampling the AR(1) country trajectories.
  start.year = 1990.5, ##<< First year of estimation period (will be centered at half-year)
  ## If given user-input, it will use min of 1990 and user input
  end.year = 2015.5, ##<< First year of estimation period (will be centered at half-year)
  ## If given user-input, it will use max of 2015 and user input
  years.change = matrix(c(1990.5, 2000.5,
                          2000.5, 2010.5,
                          1990.5, 2010.5), 
                        3, 2, byrow = TRUE), ##<< Matrix with 2 columns, with column 1 
  ## containing yyyy1 and column 2 containing yyyy2 for calculating change yyyy1-yyyy2 
  years.change2 = matrix(c(2005.5, 2010.5, 2015.5,
                           2000.5, 2005.5, 2010.5,
                           1995.5, 2000.5, 2005.5,
                           1990.5, 1995.5, 2000.5,
                           1990.5, 2000.5, 2010.5),                         
                         5, 3, byrow = TRUE) ##<< Matrix with 3 columns, with column 1 
  ## containing yyyy1, column 2 containing yyyy2 and column 3 containing yyyy3 for 
  ## calculating change (yyyy2-yyyy3) - (yyyy1-yyyy2) 
  ## The years 1990, 2000 and 2010 are always included.
  ## Years outside (\code{start.year}, \code{end.year}) are ignored. 
  ## Mid-point years closest to the given \code{years.change} are used for calculations.
  ){
  nrepeatARsampling = 1 # Removed from set of input (July 2012), not to be changed.
  if (is.null(output.dir))
    output.dir <- file.path(getwd(), "output", run.name, "/")
  filename.append <- ifelse(do.SS.run.first.pass, "_pre", "")
  load(file.path(output.dir, paste0("mcmc.meta", filename.append, ".rda"))) # change JR, 20140418
  load(file.path(output.dir, paste0("mcmc.array", filename.append, ".rda"))) # change JR, 20140418
  parnames.list <- GetParNames(winbugs.data = mcmc.meta$winbugs.data, 
                               validation.list  = mcmc.meta$validation.list,
                               do.country.specific.run = mcmc.meta$general$do.country.specific.run) # change JR, 20131104
  #------------------------------------------------------------------------------------------
  if (!file.exists(file.path(output.dir, paste0("res.country", filename.append, ".rda")))) {
    res.country <- GetCIs(mcmc.meta  = mcmc.meta,
                          mcmc.array = mcmc.array, 
                          MWRA.csv = MWRA.csv,
                          # thin = 1000, 
                          seed = seedAR,
                          include.AR = mcmc.meta$include.AR,
                          output.dir= output.dir,
                          nrepeatARsampling = nrepeatARsampling,
                          start.year = start.year,
                          end.year = end.year,
                          years.change = years.change,
                          years.change2 = years.change2) # change JR, 20140317
    save(res.country, file = file.path(output.dir, paste0("res.country", filename.append, ".rda"))) # change JR, 20140418
  } else {
    stop(paste0("res.country", filename.append, ".rda already exists!\n"))
  }
  #------------------------------------------------------------------------------------------
  if (!do.SS.run.first.pass) {
  # Validation
  if (!is.null(mcmc.meta$validation.list)){
    Ps <- GetPercentilesLeftOut(data = mcmc.meta$data.raw$data, 
                                mcmc.array = mcmc.array, 
                                winbugs.data = mcmc.meta$winbugs.data) 
    save(Ps, file = file.path(output.dir, "Ps_validation.rda")) # change JR, 20140418
    return(invisible()) # no aggregates etc constructed for validation run
  }
  #------------------------------------------------------------------------------------------
  # Posteriors of model parameters of logistic curves
  par.ciq <- GetParInfo(mcmc.array = mcmc.array, parnames.list = parnames.list, 
                        winbugs.data = mcmc.meta$winbugs.data, country.info = mcmc.meta$data.raw$country.info)
  save(par.ciq, file = file.path(output.dir, paste0("par.ciq", filename.append, ".rda"))) # change JR, 20140418
  #------------------------------------------------------------------------------------------
  do.country.specific.run <- mcmc.meta$general$do.country.specific.run
  if (!do.country.specific.run) {
    res.aggregate <- GetAggregates(run.name = run.name,
                                   years.change = years.change,
                                   years.change2 = years.change2) # change JR, 20140317
    save(res.aggregate, file = file.path(output.dir, "res.aggregate.rda")) # change JR, 20140418
  }
  cat(paste0("Results constructed and saved to ", output.dir, "\n"))
  }
  #------------------------------------------------------------------------
  ##value<< NULL, output objects from class
  ## \code{\link{Results}} for all countries (\code{res.country.rda})
  ## and UNDP aggregates (\code{res.aggregate.rda}) are written to \code{output.dir}.
  ## Additional output written to the same directory: \code{par.ciq} with CIs for the parameters
  ## of the logistic curves for total and modern/total.
  ## For validation runs, output object \code{Ps} from class \code{\link{Validation.Results}}
  ## are added.
  return(invisible(NULL))
}
#--------------------------------------------------------------------------
# The End!
