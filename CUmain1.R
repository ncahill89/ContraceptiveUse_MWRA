#--------------------------------------------------------------------------------
# Leontine Alkema, 2012
# CU_main1.R
#--------------------------------------------------------------------------------

# OVERVIEW:
# 0. Install packages (do this once!)
# 1. Change settings
# 2. Run MCMC 
# To work with MCMC output: go to CUmain2.R

# 3. (Optional) Additional runs
# 4. (Optional) Extra R code to print info about data, and plot data only, get table with data info

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

#--------------------------------------------------------------------------------
# 1.Settings
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

# Specify a run name
run.name <- "Run20160822"
### END CHANGE
#------------------------------------------------------------------------
# 2. Start MCMC run
#------------------------------------------------------------------------
# load libraries
# library(ContraceptiveUse)
library(MCMCpack)
library(rjags)
library(R2jags)
library(lattice)
library(abind)
library(msm)
library(foreach)
library(proto)
library(plyr)
library(reshape)
library(ggplot2)
library(xtable)

if (run.on.server) {
  # On Linux machine:
  library(doMC)
  registerDoMC()
}

Rfiles <- list.files(file.path(getwd(), "R"))
Rfiles <- Rfiles[grepl(".R", Rfiles)]
sapply(paste0("R/", Rfiles), source)


RunMCMC(run.name = run.name, 
        data.csv = "data/dataCPmodel.csv", regioninfo.csv = "data/Country-and-area-classification.csv"
        ,ChainNums = seq(1,2)
) # for default settings
# see the help file of RunMCMC to change settings:
#?RunMCMC

#------------------------------------------------------------------------
# 3. (Optional) Additional runs
#------------------------------------------------------------------------

# To add an MCMC chain to run.name:
# AddMCMCChain(run.name = run.name, ChainNums = 6)

# # To do a very short test run, called "test":
#  RunMCMC(data.csv = "data/dataCPmodel.csv", regioninfo.csv = "data/Country-and-area-classification.csv",
#     run.name = "test", 
#    N.ITER = 5,  
#    N.STEPS = 1,     
#    seed.MCMC = 1, 
#    N.THIN = 1,      
#    N.BURNIN = 1, 
#    ChainNums = c(1,2))

# To do a run with updated data sets for CP use and/or region:
# RunMCMC(run.name = "UpdateRun20120215", 
#         data.csv = "data/dataUpdated.csv",
#         regioninfo.csv = "data/regionUpdated.csv")
#}

# if (FALSE) {
# #-------------------------------------------------------
# # 4. (Optional) 
# # Extra R code to print info about data, and plot data only
# #-------------------------------------------------------
# 
# # Create a directory "fig" in current working directory to store figures
# dir.create("fig", showWarnings = FALSE)
# 
# # To print info about data to file:
# data.raw <- ReadDataAll(data.csv = "data/dataCPmodel.csv", regioninfo.csv = "data/Country-and-area-classification.csv",html.file="test.html") 
# #?ReadDataAll
# 
# 
# # To plot the data
# pdf(paste("fig/data_all.pdf", sep = ""), width = 24, height = 12)
#  PlotDataAndEstimates(data = data.raw)
# dev.off()
# 
# # To get overview figure of data availability
# pdf(paste("fig/", "datainfo_total.pdf", sep = ""), width = 6, height = 12)
#  PlotDataAvailability(data = data.raw$data, country.info = data.raw$country.info, 
#                        summarize.unmet = FALSE)
# dev.off()     
# 
# pdf(paste("fig/", "datainfo_unmet.pdf", sep = ""), width = 6, height = 12)
#   PlotDataAvailability(data = data.raw$data, country.info = data.raw$country.info, 
#                        summarize.unmet = TRUE)
# dev.off()     
# }
# # # extra
# # # to get number of multipliers
# ###Find in F_outputdatadummies
#  output.dir <- file.path(getwd(), "output", run.name, "/")
#  load(paste(output.dir, "mcmc.meta.rda", sep = ""))
#  winbugs.data = mcmc.meta$winbugs.data
#  res <- SummarizeMultipliers(winbugs.data = winbugs.data,
#                              data = data.raw)
#  data.frame(res)
#  
#  ###Example using SA to show where summaries come from.
#  data<-data.raw$data
#  select <- !is.na(data$props.tot.j) # change JR, 20131120
#  
#  sa.factor <- as.factor(c(ifelse(data$poptype.j=="SA" & select, data$name.j, ""), "")) # change JR, 20131120
#  sa.ind.j <- as.numeric(sa.factor)
#  ncat.sa <- nlevels(sa.factor)
#  
# # # number of biases (or get this from info.html)
#  res <- SummarizeBiases( winbugs.data = winbugs.data)
#  res
#----------------------------------------------------------------
# The End!