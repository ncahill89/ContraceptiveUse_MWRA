#--------------------------------------------------------------------------
# plotresults.R
# Leontine Alkema, 2011
#--------------------------------------------------------------------------
PlotResults <- function(# Plot lots of results!
  ### Wrapper function to plot lots of results.
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  fig.dir = NULL, ##<< Directory to store overview plots. If NULL, folder "fig" in current working directory.
  do.SS.run.first.pass = FALSE, ##<< do first pass of run with SS data? # change JR, 20140414
  plot.ind.country.results = FALSE, ##<< Create zillion plots for all countries?
  ## If TRUE, plots are saved in subdirectory "country.plots" in fig.dir.
  plot.prior.post = TRUE, ##<< If TRUE, plots of priors and posteriors are generated.
  plot.estimates = TRUE, ##<< If TRUE, plots of estimates are generated.
  plot.parameters = TRUE, ##<< If TRUE, plots of parameters are generated.
  start.year = NULL, ##<< Which years are plotted? Defaults to estimation years used in \code{CI.Lg.Lcat.qt}, or 1990 otherwise.
  end.year = NULL, ##<< Which years are plotted? Defaults to estimation years used in \code{CI.Lg.Lcat.qt}, or 1990 otherwise.
  ...
){
  if (is.null(fig.dir)){
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  if (plot.ind.country.results){
    fig.dir.countries <- file.path(fig.dir, "country.plots/")
  }
  # put separate?
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  filename.append <- ifelse(do.SS.run.first.pass, "_pre", "")
  if (do.SS.run.first.pass & !file.exists(file.path(output.dir, paste0("par.ciq", filename.append, ".rda")))) {
    par.ciq <- NULL
  } else {
    load(file.path(output.dir, paste0("par.ciq", filename.append, ".rda")))
  }
  load(file.path(output.dir, paste0("mcmc.meta", filename.append, ".rda")))
  do.country.specific.run <- mcmc.meta$general$do.country.specific.run
  load(file.path(output.dir, paste0("res.country", filename.append, ".rda")))
  if (!do.country.specific.run)
    load(file.path(output.dir, "res.aggregate.rda"))
  validation <- !is.null(mcmc.meta$validation.list)
  data <- mcmc.meta$data.raw$data
  country.info <- mcmc.meta$data.raw$country.info
  region.info <- mcmc.meta$data.raw$region.info
  winbugs.data <- mcmc.meta$winbugs.data

  ##details<< For validation run you have to use different function \code{\link{PlotValidationResults}}. Otherwise the following plots are made:
  if (validation){
    print("Use PlotValidationResults!")
    return(invisible())
  }# end validation results
  
  ##details<< Plot priors and posteriors using \code{\link{PlotPriorPost}}.
  if (plot.prior.post)
    PlotPriorPost(run.name = run.name)
  #------------------------------------------------------------------------------------------
  if (plot.estimates) {
  ##details<< Plot country overview plots for proportions with and without details using
  ##\code{\link{PlotDataAndEstimates}}.
  fig.name.years <- ifelse(!is.null(start.year) | !is.null(end.year),
                           paste0("_",
                                  ifelse(!is.null(start.year), paste0("from", floor(start.year)), ""),
                                  ifelse(!is.null(end.year), paste0("to", floor(end.year)), "")),
                           "")
  PlotDataAndEstimates(data.raw = mcmc.meta$data.raw, 
                       par.ciq = par.ciq, 
                       CI.Lg.Lcat.qt = res.country$CIprop.Lg.Lcat.qt,
                       CIstar.Lg.Lcat.qt = res.country$CIstar.Lg.Lcat.qt,
                       CIratio.Lg.Lcat.qt = res.country$CIratio.Lg.Lcat.qt,
                       start.year = start.year,
                       end.year = end.year,
                       fig.name = file.path(fig.dir, paste0(run.name, filename.append, "CIs", fig.name.years, ".pdf")),
                       ...)
  
  PlotDataAndEstimates(data.raw = mcmc.meta$data.raw, 
                       CI.Lg.Lcat.qt = res.country$CIprop.Lg.Lcat.qt,
                       CIratio.Lg.Lcat.qt = res.country$CIratio.Lg.Lcat.qt,
                       start.year = start.year,
                       end.year = end.year,
                       fig.name = file.path(fig.dir, paste0(run.name, filename.append, "CIs_nopar", fig.name.years, ".pdf")),
                       ...)
  
  # to plot individual country results
  if (plot.ind.country.results) {
    PlotDataAndEstimates(data.raw = mcmc.meta$data.raw, 
                         CI.Lg.Lcat.qt = res.country$CIprop.Lg.Lcat.qt,
                         CIratio.Lg.Lcat.qt = res.country$CIratio.Lg.Lcat.qt,
                         ind.country.overviewplot = TRUE,
                         start.year = start.year,
                         end.year = end.year,
                         run.name = run.name,
                         ...)
  }
  
  ##details<< Plot country overview plots for counts using
  ##\code{\link{PlotDataAndEstimates}}.
  PlotDataAndEstimates(
    CI.Lg.Lcat.qt = res.country$CIcount.Lg.Lcat.qt,
    plot.prop = FALSE, 
    start.year = start.year,
    end.year = end.year,
    fig.name = file.path(fig.dir, paste0(run.name, filename.append, "CIscountcountry", fig.name.years, ".pdf")),
    ...)
  
  ##details<< Plot UNDP aggregates overview plots for proportions and counts using
  ##\code{\link{PlotDataAndEstimates}}.
  if (!do.country.specific.run) {
    PlotDataAndEstimates(CI.Lg.Lcat.qt = res.aggregate$CIprop.Lg.Lcat.qt,
                         CIratio.Lg.Lcat.qt = res.aggregate$CIratio.Lg.Lcat.qt,
                         start.year = start.year,
                         end.year = end.year,
                         fig.name = file.path(fig.dir, paste0(run.name, filename.append, "CIsaggregate", fig.name.years, ".pdf")),
                         ...)
    PlotDataAndEstimates(
      CI.Lg.Lcat.qt = res.aggregate$CIcount.Lg.Lcat.qt,
      plot.prop = FALSE, 
      start.year = start.year,
      end.year = end.year,
      fig.name = file.path(fig.dir, paste0(run.name, "CIscountaggregate", fig.name.years, ".pdf")),
      ...) 
    if (plot.ind.country.results){
      PlotDataAndEstimates(CI.Lg.Lcat.qt = res.aggregate$CIprop.Lg.Lcat.qt,
                           CIratio.Lg.Lcat.qt = res.aggregate$CIratio.Lg.Lcat.qt,
                           ind.country.overviewplot = TRUE,
                           start.year = start.year,
                           end.year = end.year,
                           run.name = run.name,
                           ...)
    }
    PlotCountryEstimatesForAggregate(
      CI.Lg.Lcat.qt = res.country$CIprop.Lg.Lcat.qt,
      country.info = country.info,
      fig.name = file.path(fig.dir, paste0(run.name, "CIscountryestsinaggregates", fig.name.years, ".pdf")))
  }
  }
  #------------------------------------------------------------------------------------------
  if (plot.parameters) {
  ##details<< Plot model parameters using \code{\link{PlotLogisticParameters}}.
  if (is.null(par.ciq)) {
    stop("par.ciq is NULL. Model parameters cannot be plotted.")
  } else {
  PlotLogisticParameters(par.ciq = par.ciq, country.info = country.info,
                         region.info = region.info,
                         fig.name = file.path(fig.dir, paste0(run.name, filename.append, "modelpar.pdf"))  ) # change JR, 20140418
  }
  #------------------------------------------------------------------------------------------
  ##details<< Plot info about data parameters (variance by source, biases and perturbation multipliers). 
  load(file.path(output.dir, paste0("mcmc.array", filename.append, ".rda"))) # change JR, 20140418
  
  if (!do.country.specific.run) { # change JR, 20131105
    # variance by source
    pdf(file.path(fig.dir, paste0(run.name, "sigmasource.pdf")), width = 7, height = 5)
    PlotSourceSigmas(mcmc.array = mcmc.array)
    PlotSourceSigmasUnmet(mcmc.array = mcmc.array)
    dev.off()
    # biases,
    pdf(file.path(fig.dir, paste0(run.name, "Biases.pdf")), width = 7, height = 5, useDingbats = FALSE)
    PlotBiases(mcmc.array = mcmc.array)
    dev.off()
  }
  
  # (long!) list of names of unique multipliers
  par.V <- mcmc.meta$par.V
  #InternalGetParnamesV(winbugs.data = winbugs.data, 
  #                      name.short.j = InternalMakeCountryNamesShort(mcmc.meta$data.raw$data$name.j))
  #par.V$parnames.V.in.bugs
  #par.V$parnames.V.nice
  #data.frame(par.V$parnames.V.in.bugs,  par.V$parnames.V.nice)
  pdf(file.path(fig.dir, paste0(run.name, filename.append, "Multipliers.pdf")))#, height = 15, width = 20)
    PlotGeoEtc(par.V = par.V, mcmc.array = mcmc.array)
  dev.off()
  
  # # To check data dummies:
  # # number of observations with bias (folk/mics/M+/M-): 
  # summary.biases <-SummarizeBiases(winbugs.data) # 237 biases
  # # number of unique multipliers included (per composition):
  # summary.multipliers <- SummarizeMultipliers(winbugs.data, data) # 214 sets
  }
  cat("All results plotted and saved to ", fig.dir, "\n")
  ##value<< NULL
  return(invisible(NULL))
} # end function
#----------------------------------------------------------------------   
# The End!
