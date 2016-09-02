#--------------------------------------------------------------------------------
# Leontine Alkema and Jin Rou New, 2012-2015
# CUmain_cs.R
#--------------------------------------------------------------------------------
# See manual for additional information

# OVERVIEW:
# 0. Install packages (do this once!)
# 1. Change settings
# 2. Run MCMC
# 3. Construct output and make tables/plots

# Note:
# When running this script in R-studio
# Use the commands [ctr] + Enter to run selected lines of code on a pc
# When running this script in R
# Use the commands [ctr] + r to run selected lines of code on a pc

if (FALSE) {
  #--------------------------------------------------------------------------------
  # 0. Install packages
  #--------------------------------------------------------------------------------
  # Install packages from CRAN:
  install.packages(pkgs = c("xtable", "foreach", "MCMCpack", "rjags", "R2jags", 
                            "lattice","abind", "proto", "msm", "reshape"))
  
  # add on Windows PC:
  install.packages(pkgs = c("ggplot2",  "plyr"))
  
  # Add on a Linux machine (to run chains in parallel):
  #install.packages(pkgs = c("doMC"))
  
  # Install ContraceptiveUse from local zip-file, two options for that:
  # 1. Go to ``Tools'', "Install packages'', and select \code{ContraceptiveUse.zip}
  # Or
  # 2. Use
  ## install.packages(pkgs = "...\ContraceptiveUse.zip"),
  ##                 repos = NULL)
  # where ... refers to the directory where ContraceptiveUse.zip is stored.
  # If you have saved the zip-file in your current working directory, you can use this command:
  #install.packages(pkgs = file.path(getwd(), "ContraceptiveUse.zip"),
  #                 repos = NULL)
  
  # Note:
  # Next time you use this R script (after installing the packages)
  # You can skip the install steps.
}
#--------------------------------------------------------------------------------
# 1. Settings and load libraries and codes
#--------------------------------------------------------------------------------
rm(list = ls()); gc()
run.on.server <- FALSE

if (!run.on.server) {
  ### CHANGE THIS:
  # Set your working directory, choose one!
  # Note: Reverse the "/" when copying the name of the directory from Windows explorer
  mydir <- "/Users/admin/Dropbox/*PostDoc2016/ContraceptiveUse_Original"
  setwd(mydir)
}
# Load libraries
# library(ContraceptiveUse)
#library(MCMCpack)
library(rjags)
library(R2jags)
library(lattice)
library(abind)
library(msm)
library(proto)
library(plyr)
library(reshape2)
library(ggplot2)

if (run.on.server) { # On Linux machine:
  library(foreach)
  library(doMC)
  registerDoMC()
}

# library(RPushbullet)
# options(error = function() { # Be notified when there is an error
#   pbPost(type = "note", 
#          title = "Error!", 
#          body = geterrmessage(),
#          recipients = c(1, 2))
# })

# pbPost(type = "note", 
#        title = paste0("CUmain_cs.R"), 
#        body = paste0("Job started!"),
#        recipients = c(1, 2))


Rfiles <- list.files(file.path(getwd(), "R"))
Rfiles <- Rfiles[grepl(".R", Rfiles)]
sapply(paste0("R/", Rfiles), source)

# User settings
monitor.time.taken <- TRUE
# specify run name
run.name.global <- "RUN20160523"
#GlobalSum<-SummariseGlobalRun(run.name = run.name.global)
#load(file="output/Test/data.global.rda")

iso.country.select <- NULL
regioninfo.csv <- "data/Country-and-area-classification.csv"
MWRA.csv <- "data/Number-of-women-married-in union_15-49.csv"
data.csv <- "data/dataCPmodel.csv"
data <- read.csv(data.csv, stringsAsFactors = FALSE)
if ("Country..letter.code" %in% names(data)) {
  isos.all <- unique(data$Country..letter.code[is.na(data$EXCLUDE1isyes) | data$EXCLUDE1isyes != 1])
} else if ("Country.letter.code" %in% names(data)) {
  isos.all <- unique(data$Country.letter.code[is.na(data$EXCLUDE1isyes) | data$EXCLUDE1isyes != 1])
} else {
  stop("Neither of the columns Country..letter.code nor Country.letter.code was found in data.csv.")
}
isos.all <- gsub(" ", "", isos.all)
# Run for FP2020 countries only
info <- read.csv("data/Country-and-area-classification-inclFP2020.csv", stringsAsFactors = FALSE)
info.codes <- read.csv("data/Country-names-and-codes.csv", stringsAsFactors = FALSE)[, c("Country.letter.code", "ISO.code")]
names(info.codes)[names(info.codes) == "ISO.code"] <- "ISO.Code"
info <- join(info, info.codes, by = "ISO.Code")
isos.fp2020 <- info$Country.letter.code[info$FP2020 == "Yes"]
isos.all <- intersect(isos.all, isos.fp2020) # Only ESH, Western Sahara dropped because no data available
print(isos.all)
print(length(isos.all))
### END CHANGE
#------------------------------------------------------------------------
# 2. Start MCMC run
#------------------------------------------------------------------------
iso.select<-"IND"
isos.all<-isos.all[which(isos.all==iso.select)]

for (iso.select in isos.all) {
  run.name <- paste0(run.name.global, "_FP2020_", iso.select)
  
   #Start a file to output run times for each country
  if (monitor.time.taken)
    cat(paste("ISO", paste(names(proc.time()), collapse = ","), sep = ","), 
        file = file.path("output", run.name.global, paste0("onecountryruntimes_noSS.txt")), 
        fill = T, append = F)
  
  if (file.exists(file.path("output", run.name)))
    next()
  
  # check if do.SS.run
  data.preprocessed <- PreprocessData(data.csv = data.csv,
                                      iso.select = iso.select)
  select.ss <- grepl("Service statistic", data.preprocessed$Data.series.type)
  do.SS.run <- any(select.ss)
  
  if (!do.SS.run) {
    RunMCMC(run.name = run.name, 
            do.country.specific.run = TRUE,
            iso.select = iso.select, 
            iso.country.select = iso.country.select,
            run.name.global = run.name.global,
            run.on.server = run.on.server,
            data.csv = data.csv, 
            regioninfo.csv = regioninfo.csv
            #,ChainNums = seq(1,2)
            #include.AR = FALSE
    ) # for default settings
    # see the help file of RunMCMC to change settings:
    #?RunMCMC
  } else {
    message("Service statistics data is available.")
    # First pass
    RunMCMC(run.name = run.name, 
            do.country.specific.run = TRUE,
            iso.select = iso.select, 
            iso.country.select = iso.country.select,
            run.name.global = run.name.global,
            run.on.server = run.on.server,
            data.csv = data.csv, 
            regioninfo.csv = regioninfo.csv,
            ChainNums = seq(1,2)
            
    )
    closeAllConnections()
    ConstructMCMCArray(run.name = run.name, do.SS.run.first.pass = TRUE)
    ConstructOutput(run.name = run.name,
                    MWRA.csv = MWRA.csv
                    , start.year = 1970.5, end.year = 2020.5, # change to shorter period including year required to estimate bias.modern to speed up? 
                    do.SS.run.first.pass = TRUE)
    #----------------------------------------------------------------------
    # Second pass
    RunMCMC(run.name = run.name, 
            do.country.specific.run = TRUE,
            iso.select = iso.select, 
            iso.country.select = iso.country.select,
            run.name.global = run.name.global,
            run.on.server = run.on.server,
            data.csv = data.csv, 
            regioninfo.csv = regioninfo.csv,
            ChainNums = seq(1,2)
            
    )
  }
  if (monitor.time.taken)
    cat(paste(iso.select, paste(proc.time(), collapse = ","), sep = ","), 
       file = file.path("output", run.name.global, paste0("onecountryruntimes_noSS.txt")),
        fill = T, append = T)
  closeAllConnections()  
  ConstructMCMCArray(run.name = run.name)
  # CheckConvergence(run.name = run.name, 
  #                  plot.trace = TRUE, 
  #                  check.convergence = TRUE,
  #                  use.all.parameters = TRUE) 
  ConstructOutput(run.name = run.name,
                  MWRA.csv = MWRA.csv
                  # , start.year = 1990.5, end.year = 2020.5
                  , start.year = 1970.5, end.year = 2035.5
  )
  
  if (FALSE) {
    # Plots of trends over time in props/counts for countries and aggregates
    PlotResults(run.name = run.name,
                # Do you want to save individual country/aggregate plots in tiff?
                # If yes, set this to TRUE, then TIFFS are saved in subfolder in "fig"
                start.year = 1970.5, end.year = 2035.5,
                plot.ind.country.results = FALSE)
    PlotResults(run.name = run.name,
                # Do you want to save individual country/aggregate plots in tiff?
                # If yes, set this to TRUE, then TIFFS are saved in subfolder in "fig"
                start.year = 1990.5, end.year = 2020.5,
                plot.ind.country.results = FALSE)
    
    # Tables 
    # Get tables (csv's) for estimates:
    GetTablesRes(run.name = run.name, name.res = "Country")
    # Get tables (csv's) for changes:
    GetTablesChange(run.name = run.name, name.res = "Country")
  } 
}
# pbPost(type = "note", 
#        title = paste0("CUmain_cs.R"), 
#        body = paste0("Job done!"),
#        recipients = c(1, 2))
#----------------------------------------------------------------
# The End!
