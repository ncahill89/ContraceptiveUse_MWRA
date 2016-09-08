#----------------------------------------------------------------------
# summariseglobalrun.R
# Jin Rou New, Nov 2013
# Modified, NC 2016
#----------------------------------------------------------------------
SummariseGlobalRun <- function(# Summarise info from global run required for country-specific runs
  ### Summarise info from global run required for country-specific runs
  run.name = "test", ##<< Run name
  output.dir = NULL ##<< Directory where MCMC array and meta are stored, and new objects are added
  ## (if NULL, it's output/run.name, default from \code{runMCMC}).
) {
  if (is.null(output.dir))
    output.dir <- file.path(getwd(), "output", run.name)
  
  load(file.path(output.dir, "mcmc.meta.rda"))
  load(file.path(output.dir, "mcmc.array.rda"))
  
  parnames.list <- mcmc.meta$parnames.list
  parnames.to.save <- c(parnames.list$parnames.h,
                        c(sapply(parnames.list$parnames.subreg, paste0, "[", 1:mcmc.meta$winbugs.data$n.subreg, "]")),
                        c(sapply(parnames.list$parnames.reg, paste0, "[", 1:mcmc.meta$winbugs.data$n.reg, "]")), # change JR, 20150301
                        paste0(parnames.list$T1.source.s, "[", 1:4, "]"),
                        paste0(parnames.list$T2.source.s, "[", 1:4, "]"),
                        paste0(parnames.list$T12.source.s, "[", 1:4, "]"),
                        c(sapply(c("unmet.intercept.c", "setlevel.c", "RT.c",
                                 "omega.c", "Romega.c"), paste0, "[", 1:mcmc.meta$winbugs.data$C, "]")))                        
  parnames.not.to.save <- c("w.world", "Rw.world", "S.world", "RT.world",
                            "sigma.wreg", "sigma.Rwreg", "sigma.Sreg", "sigma.RTreg"
                            # , "sigma.wsubreg", "sigma.Rwsubreg", "sigma.Tsubreg", "sigma.RTsubreg" # change JR, 20150301
                            )
  parnames.to.save <- parnames.to.save[!is.element(parnames.to.save, parnames.not.to.save)]   
  
  # check that all required parameters are found in the mcmc.array
  if (sum(!is.element(parnames.to.save, dimnames(mcmc.array)[[3]])) > 0) {
    cat("Warning: The following parameters cannot be found in the mcmc.array:\n")
    cat(setdiff(parnames.to.save, dimnames(mcmc.array)[[3]]))
  }
  parnames.to.save <- intersect(parnames.to.save, dimnames(mcmc.array)[[3]])
  
  # get posterior median of all parameters to save
  mcmc.post <- lapply(parnames.to.save, function(parname) median(c(mcmc.array[, , parname]), na.rm = T))
  names(mcmc.post) <- parnames.to.save
  iso.c <- gsub(" ", "", mcmc.meta$data.raw$country.info$code.c)
  names(iso.c) <- mcmc.meta$data.raw$country.info$name.c
  data.global <- list(run.name.global = run.name, 
                      name.subreg = as.character(mcmc.meta$data.raw$region.info$name.subreg),
                      name.reg = as.character(mcmc.meta$data.raw$region.info$name.reg), # change JR, 20150301
                      iso.c = iso.c,
                      cnotrich.index = mcmc.meta$winbugs.data$cnotrich.index,
                      mcmc.post = mcmc.post)
  save(data.global, file = file.path(output.dir, "data.global.rda"))
  ##value<< \code{NULL}; Saves \code{data.global} to \code{output.dir}.
  return(invisible())
}
#----------------------------------------------------------------------   
# The End!
