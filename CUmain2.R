#--------------------------------------------------------------------------------
# Leontine Alkema and Jin Rou New
# CUmain2.R
#--------------------------------------------------------------------------------
# Run after CUmain1.R

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
rm(list = ls())
run.on.server <- FALSE

if (!run.on.server) {
  ### CHANGE THIS:
  # Set your working directory, choose one!
  # Note: Reverse the "/" when copying the name of the directory from Windows explorer
  mydir <- "/Users/admin/ContraceptiveUse_MWRA"
  #mydir <- "V:/FertilitySection/R/R code" # Fertility section
  #mydir <- getwd() # LA on server
  #mydir <- "V:/FertilitySection/Methodological issues/Alkema_Joint project on contraceptive use trends/R code/"
  setwd(mydir)
}

# specify run name
run.name <- "Run20160822"
### END CHANGE

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

Rfiles <- list.files(file.path(getwd(), "R"))
Rfiles <- Rfiles[grepl(".R", Rfiles)]
sapply(paste0("R/", Rfiles), source)

#--------------------------------------------------------------------------------
# 2. Read MCMC info (e.g. construct MCMC array) 
#--------------------------------------------------------------------------------
ConstructMCMCArray(run.name = run.name)
#load(file="output/Test4/mcmc.array.rda")
#attributes(mcmc.array)
# In default call, output is stored in "output/run.name/.."

#--------------------------------------------------------------------------------
# 3. Diagnostics
#--------------------------------------------------------------------------------
# Note: you need at least two chains and a large number of iterations (e.g 1000 per chain)
# Default is to check convergence for a subset of countries/AR-distortions/multipliers only.

# This can take a while!
# CheckConvergence(run.name = run.name, 
#                  plot.trace = TRUE, 
#                  check.convergence = TRUE) 
# Note: if this function crashes, a text file or pdf might still be open,
# just restart R or use following comments to close things:
#sink()
#closeAllConnections()
#dev.off()

#--------------------------------------------------------------------------------
# 4. Construct output (CIs etc)
#--------------------------------------------------------------------------------

ConstructOutput(run.name = run.name,
                MWRA.csv = "data/Number-of-women-married-in union_15-49.csv",
                start.year = 1970.5,
                end.year = 2035.5) 
#?ConstructOutput

# If you want results for user-specified years, e.g. from 1970 to 2020, use this:
# ConstructOutput(run.name = run.name, 
#                 start.year = 1970.5, end.year = 2020.5)
# But note that you'll need MWRA numbers for those years as well!

# If you want change results for non-default years, use this:
# years.change <- matrix(c(yyyy1, yyyy2,
#                          yyyy1, yyyy2),
#                        nrow = 2, ncol = 2, byrow = TRUE) # change nrow accordingly
# years.change2 <- matrix(c(yyyy1, yyyy2, yyyy3
#                           yyyy1, yyyy2), yyyy3
#                         nrow = 2, ncol = 3, byrow = TRUE) # change nrow accordingly
# ConstructOutput(run.name = run.name,
#                 years.change = years.change, 
#                 years.change2 = years.change2)

# If you want to use an updated data set for MWRA, use this:
# ConstructOutput(run.name = run.name,  end.year = 2020.5,
#                 MWRA.csv = "data/Married-or-in-union-women_WPP2010_1990-2020_South-Sudan&Sudan.csv")
#--------------------------------------------------------------------------------
# 5. Results (plots and tables)
#--------------------------------------------------------------------------------
# By using the default function calls,
# all figures and csv-files are saved in directories
# "fig" and "tables" in your current working directory.
# The run.name is added to the name of the figure/csv.

# TIFFS note: by adding the plot.tiff  = TRUE statement, tiffs are produced.
##Summarize Global Run for FPET
run.name.global <- run.name
GlobalSum<-SummariseGlobalRun(run.name = run.name.global)

# Standard output:
GetAllBarCharts(run.name = run.name)#, plot.tiff  = TRUE)
 #GetAllPolar(run.name = run.name)

# Plots of trends over time in props/counts for countries and aggregates
PlotResults(run.name = run.name,
            # Do you want to save individual country/aggregate plots in tiff?
            # If yes, set this to TRUE, then TIFFS are saved in subfolder in "fig"
            # start.year = 1990.5,
            # end.year = 2020.5,
            plot.ind.country.results = FALSE,
            plot.prior.post = FALSE,
            plot.parameters = TRUE)

# Plots in paper:
 BarChartSubregion(run.name = run.name)#, plot.tiff  = TRUE)
 CIPropChangesSubregions(run.name = run.name)#, plot.tiff  = TRUE)
# PlotMetDemandinCountryBySubregion(run.name = run.name)

# Note: if a function crashes, the pdf might still be open,
# use "dev.off()" several times until you get an error
#dev.off()

# Tables 
# Get tables (csv's) for estimates:
GetTablesRes(run.name = run.name, name.res = "Country")
GetTablesRes(run.name = run.name, name.res = "UNPDaggregate")
# Get tables (csv's) for changes:
GetTablesChange(run.name = run.name, name.res = "UNPDaggregate")
GetTablesChange(run.name = run.name, name.res = "Country")

# If you want to construct estimates for different aggregates,
# Go to "Main_ExtraAggregates.R"
#-------------------------------------------------------------------
# The End!