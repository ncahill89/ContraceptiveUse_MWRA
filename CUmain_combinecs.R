#--------------------------------------------------------------------------------
# Leontine Alkema and Jin Rou New, 2012-2015
# CUmain_combinecs.R
#--------------------------------------------------------------------------------
# Run after CUmain_cs.R

#--------------------------------------------------------------------------------
# OVERVIEW of steps:
# 1. Settings
# 2. Read MCMC info (construct MCMC array)
# 3. Diagnostics
# 4. Construct output (e.g. CIs)
# 5. Results (plots and tables)

#--------------------------------------------------------------------------------
# 1. Settings
#--------------------------------------------------------------------------------
rm(list = ls())# clear work space
run.on.server <- FALSE

if (!run.on.server) {
  ### CHANGE THIS:
  # Set your working directory
  # Use the same working directory as in CUmain1.R if you want to use default function calls.
  # Note: Reverse the "/" when copying the name of the directory from Windows explorer
  mydir <- "/Users/admin/Dropbox/*PostDoc2016/ContraceptiveUse"
  setwd(dir = mydir)
}


# load libraries
library(MCMCpack)
library(rjags)
library(R2jags)
library(lattice)
library(abind)
library(plyr)
library(proto)
library(ggplot2)
library(reshape2)
# library(ContraceptiveUse)

# library(RPushbullet)
# options(error = function() { # Be notified when there is an error
#   pbPost(type = "note", 
#          title = "Error!", 
#          body = geterrmessage(),
#          recipients = c(1, 2))
# })
# 
# pbPost(type = "note", 
#        title = paste0("CUmain_combinecs.R"), 
#        body = paste0("Job started!"),
#        recipients = c(1, 2))

Rfiles <- list.files(file.path(getwd(), "R"))
Rfiles <- Rfiles[grepl(".R", Rfiles)]
sapply(paste0("R/", Rfiles), source)

# User settings
run.name.global <- "RUN20160523_FP2020"
run.name.all <- paste0(run.name.global, "_all")
paste0(run.name.global, "_all")
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
#--------------------------------------------------------------------------------
# 2. Combine country-specific runs
#--------------------------------------------------------------------------------
if (!file.exists(file.path("output", run.name.all, "mcmc.array.rda")))
  CombineCountrySpecificRuns(run.name.global = run.name.global,
                             iso.all = isos.all,
                             iso.country.select = iso.country.select,
                             data.csv = data.csv,
                             regioninfo.csv = regioninfo.csv)
#--------------------------------------------------------------------------------
# 3. Diagnostics
#--------------------------------------------------------------------------------
CheckConvergence(run.name = run.name.all, 
                 plot.trace = TRUE, 
                 check.convergence = TRUE,
                 use.all.parameters = TRUE) 
#--------------------------------------------------------------------------------
# 4. Results (plots and tables)
#--------------------------------------------------------------------------------
# Plots of trends over time in props/counts for countries and aggregates
PlotResults(run.name = run.name.all,
            # Do you want to save individual country/aggregate plots in tiff?
            # If yes, set this to TRUE, then TIFFS are saved in subfolder in "fig"
            start.year = 1970.5, end.year = 2035.5,
            plot.ind.country.results = FALSE)

# Tables 
# Get tables (csv's) for estimates:
GetTablesRes(run.name = run.name.all, name.res = "Country")
# Get tables (csv's) for changes:
GetTablesChange(run.name = run.name.all, name.res = "Country")

# Plot comparison
PlotComparison(run.name = "Run20150423_FP2020pre2012_all", 
               run.name2 = "Run20150423_FP2020_all", 
               legend = "FP2020_2015 (data before 2012)",
               legend2 = "FP2020_2015")

pbPost(type = "note", 
       title = paste0("CUmain_combinecs.R"), 
       body = paste0("Job done!"),
       recipients = c(1, 2))
#-------------------------------------------------------------------
# The End!
