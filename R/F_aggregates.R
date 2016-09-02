InternalGetAggregates <- function(#  Find aggregates for set of countries
  ### Find aggregates for set of countries
  W.Lc.t, ##<< Estimates of number of MWRA
  select.c, ##<< Selected countries to be included for the aggregate estimate
  dir.traj, ##<< Where are the country trajectories saved?
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
  # percentiles hard-coded
  percentiles <- c(0.025, 0.1, 0.5, 0.9, 0.975) 
  #as.numeric(names(res.country$CIprop.Lg.Lcat.qt[[1]][[1]][,1]))
  # find number of samples n.s:
  c <- select.c[1]
  load(file = file.path(dir.traj, paste0("P.tp3s_country", c, ".rda"))) # change JR, 20140830
  n.s <- dim(P.tp3s)[3]
  nyears <- length(W.Lc.t[[1]])
  #print(nyears)
  est.years <- as.numeric(dimnames(P.tp3s)[[1]] )#res.country$CIprop.Lg.Lcat.qt[[1]][[1]][1,]))
  #print(length(est.years))
  # change JR, 20140317
  # make sure 1990, 2000 and 2010 are included: 
  years.change <- unique(rbind(floor(years.change)+0.5,
                               matrix(c(1990.5, 2000.5,
                                        2000.5, 2010.5,
                                        1990.5, 2010.5), 
                                      3, 2, byrow = TRUE)))
  years.change2 <- unique(rbind(floor(years.change2)+0.5,
                                matrix(c(2005.5, 2010.5, 2015.5,
                                         2000.5, 2005.5, 2010.5,
                                         1995.5, 2000.5, 2005.5,
                                         1990.5, 1995.5, 2000.5,
                                         1990.5, 2000.5, 2010.5),                         
                                       5, 3, byrow = TRUE)))
  check1 <- !(c(floor(years.change)) %in% floor(est.years))
  check2 <- !(c(floor(years.change2)) %in% floor(est.years))
  if (any(check1))
    stop(paste0(c(years.change)[check1], " is not found in estimation years."))
  if (any(check2))
    stop(paste0(c(years.change2)[check2], " is not found in estimation years."))
  years.change.unique <- unique(c(c(years.change), c(years.change2)))
 
  cumsum.trad.ts  <-  cumsum.modern.ts <- cumsum.unmet.ts <- matrix(0, nyears, n.s)
  cumsumW.t <- rep(0, nyears)
  for (c in select.c){
    load(file = file.path(dir.traj, paste0("P.tp3s_country", c, ".rda"))) # change JR, 20140418
    for (t in 1:nyears){
      cumsum.trad.ts[t,] <- cumsum.trad.ts[t,] + W.Lc.t[[c]][t]*P.tp3s[t,1,]
      cumsum.modern.ts[t,] <- cumsum.modern.ts[t,] + W.Lc.t[[c]][t]*P.tp3s[t,2,]
      cumsum.unmet.ts[t,] <- cumsum.unmet.ts[t,] + W.Lc.t[[c]][t]*(P.tp3s[t,3,])
      cumsumW.t[t] <- cumsumW.t[t] + W.Lc.t[[c]][t]
    }    
  } # end country loop
  # divide by cumsumW.t (to have results same as for country results)
  CIprop.Lcat.qt <- CIratio.Lcat.qt <- list()
  
  temp <- matrix(0, length(percentiles), nyears)
  colnames(temp) <- est.years
  rownames(temp) <- percentiles
  # note: these names are the same as in the country CIs! see F_output.R
  CIprop.Lcat.qt[["Traditional"]] <- CIprop.Lcat.qt[["TotalPlusUnmet"]]<- CIprop.Lcat.qt[["TradPlusUnmet"]] <- 
    CIprop.Lcat.qt[["Modern"]] <- CIprop.Lcat.qt[["Total"]] <- CIprop.Lcat.qt[["Unmet"]] <- 
    CIratio.Lcat.qt[["Met Demand"]] <- CIratio.Lcat.qt[["Met Demand with Modern Methods"]] <- # change JR, 20140830: added demand met with modern methods
    CIratio.Lcat.qt[["Modern/Total"]]<- temp
  cumsum.tot.ts <- cumsum.trad.ts + cumsum.modern.ts
  for (t in 1:nyears){
    CIprop.Lcat.qt[["Traditional"]][,t] <- quantile(cumsum.trad.ts[t,]/cumsumW.t[t],  percentiles)
    CIprop.Lcat.qt[["TotalPlusUnmet"]][,t] <- quantile((cumsum.unmet.ts[t,]+cumsum.tot.ts[t,])/cumsumW.t[t],  percentiles)
    CIprop.Lcat.qt[["TradPlusUnmet"]][,t] <- quantile((cumsum.trad.ts[t,]+cumsum.unmet.ts[t,])/cumsumW.t[t],  percentiles)
    CIprop.Lcat.qt[["Modern"]][,t] <- quantile(cumsum.modern.ts[t,]/cumsumW.t[t],  percentiles)
    CIprop.Lcat.qt[["Total"]][,t] <- quantile(cumsum.tot.ts[t,]/cumsumW.t[t],  percentiles)
    CIprop.Lcat.qt[["Unmet"]][,t] <- quantile(cumsum.unmet.ts[t,]/cumsumW.t[t], percentiles)
    CIratio.Lcat.qt[["Met Demand"]][,t]<- quantile(cumsum.tot.ts[t,]/
                                (cumsum.tot.ts[t,]+cumsum.unmet.ts[t,]), percentiles)
    CIratio.Lcat.qt[["Met Demand with Modern Methods"]][,t]<- quantile(cumsum.modern.ts[t,]/
                                                     (cumsum.tot.ts[t,]+cumsum.unmet.ts[t,]), percentiles) # change JR, 20140830: added demand met with modern methods
    CIratio.Lcat.qt[["Modern/Total"]][,t]<- quantile(cumsum.modern.ts[t,]/
      (cumsum.tot.ts[t,]), percentiles)
  }

  # change JR, 20140317
  # find changes based on posterior samples
  P.yp3s  <- array(NA, c(length(years.change.unique), 3, dim(cumsum.trad.ts)[2]))
  dimnames(P.yp3s) <- list(years.change.unique, c("Traditional", "Modern", "Unmet"), NULL)
  # for the props:
  for (t in 1:length(years.change.unique)) {
    select <- est.years==years.change.unique[t]
    P.yp3s[t,1,] <- cumsum.trad.ts[select,]/cumsumW.t[select]
    P.yp3s[t,2,] <- cumsum.modern.ts[select,]/cumsumW.t[select]
    P.yp3s[t,3,] <- cumsum.unmet.ts[select,]/cumsumW.t[select]
  }
  changeprop.Lcat.Ti <-  GetInfoChange(P.yp3s = P.yp3s, years.change = years.change, years.change2 = years.change2)

  # for the counts:  
  for (t in 1:length(years.change.unique)) {
    select <- est.years==years.change.unique[t]
    P.yp3s[t,1,] <- cumsum.trad.ts[select,]
    P.yp3s[t,2,] <- cumsum.modern.ts[select,]
    P.yp3s[t,3,] <- cumsum.unmet.ts[select,]
  }
  changecount.Lcat.Ti <- GetInfoChange(P.yp3s = P.yp3s, type.is.prop = FALSE, years.change = years.change, years.change2 = years.change2)

  ##value<< Object of class \code{\link{Results}} with 
  return(list(
    changecount.Lcat.Ti = changecount.Lcat.Ti, ##<< Info on changes in counts
    changeprop.Lcat.Ti = changeprop.Lcat.Ti, ##<< Info on changes in proportions
    W.t = cumsumW.t,##<< Number of MWRA for the aggregate
    CIprop.Lcat.qt = CIprop.Lcat.qt, ##<< CIs for props
    CIratio.Lcat.qt = CIratio.Lcat.qt ##<< CIs for counts
    ))
}

#--------------------------------------------------------------------
GetAggregates <- function(# Construct aggregate estimates
  ### Estimate the proportion/number of MWRA in various categories for aggregates.
  run.name,##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored, 
  ## as well as folder with country trajectories
  ## If NULL, it's \code{output/run.name}, the default from \code{runMCMC}.
  file.aggregates = NULL, ##<< If NULL (default), UNDP aggregates are constructed. 
  ##Alternatively, file path of alternative grouping should be given, e.g.
  ##\code{file.aggregates = "data/MDGgroupings.csv")}. Such data file needs to contain 
  ## columns \code{iso.country}, \code{groupname} and \code{iso.group} (which may contain missing values).
  ## Each country can only be included once (can only be part of one grouping).
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
                         5, 3, byrow = TRUE), ##<< Matrix with 3 columns, with column 1 
  ## containing yyyy1, column 2 containing yyyy2 and column 3 containing yyyy3 for 
  ## calculating change (yyyy2-yyyy3) - (yyyy1-yyyy2) 
  ## The years 1990, 2000 and 2010 are always included.
  ## Years outside (\code{start.year}, \code{end.year}) are ignored. 
  ## Mid-point years closest to the given \code{years.change} are used for calculations.
  winbugs.data = NULL, ##<< Object of class \code{winbugs.data}, needed only for UNDP aggregates
  region.info = NULL##<< Object of class \code{region.info}, needed only for UNDP aggregates.
  ){
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  load(file.path(output.dir, "res.country.rda")) # change JR, 20140418
  W.Lc.t <- res.country$W.Lg.t
#  if (is.null(output.dir.countrytrajectories) | is.null(W.Lc.t)){
#    print("First save country trajectories and read in number of MWRA!")
#    return()
#  }   
  output.dir.countrytrajectories <- file.path(output.dir, "/countrytrajectories/")
  nyears <-  length(W.Lc.t[[1]])
  load(file.path(output.dir, "mcmc.meta.rda")) # change JR, 20140418
  country.info <- mcmc.meta$data.raw$country.info
  C <- length(country.info$name.c)
  res.aggregate <- list()
  if (is.null(file.aggregates)){
    # construct UNPD aggregates
    cat("Overview: Constructing aggregates for UNPD regions, and dev/dev-ing countries (excl China)", "\n")
    region.info <- mcmc.meta$data.raw$region.info
    G <-  (3+ # dev. dev-ing, dev-ing excl china 
             1+ # world
             1+ # Oceania, deving only
             region.info$n.subreg
           +region.info$n.reg -1) # -1 because northern america is subregion and region
    iso.g <- rep(NA, G) # maybe to include later!
    #G <- length(iso.g)
    W.Lg.t <- list()
    # add dev, dev-ing, dev-ing excl china:
    cat("Constructing aggregates for the developED countries","\n")
    select.c <- seq(1, C)[country.info$dev.c =="Rich"]    
    nameg <- "Developed countries"
    res.aggregate[[nameg]] <- InternalGetAggregates(W.Lc.t = W.Lc.t, select.c = select.c, 
                                                    dir.traj = output.dir.countrytrajectories, 
                                                    years.change = years.change,
                                                    years.change2 = years.change2)
    cat("Constructing aggregates for the developING countries","\n")
    select.c <- seq(1, C)[country.info$dev.c !="Rich"]    
    nameg <- "Developing countries"
    res.aggregate[[nameg]] <- InternalGetAggregates(W.Lc.t = W.Lc.t, select.c = select.c, 
                                                    dir.traj = output.dir.countrytrajectories,
                                                    years.change = years.change,
                                                    years.change2 = years.change2)
    cat("Constructing aggregates for the developing countries, excl. China","\n")
    select.c <- seq(1, C)[country.info$dev.c !="Rich" & country.info$name.c!="China"]    
    nameg <- "Developing (excl. China)"
    res.aggregate[[nameg]] <- InternalGetAggregates(W.Lc.t = W.Lc.t, select.c = select.c, 
                                                    dir.traj = output.dir.countrytrajectories,
                                                    years.change = years.change,
                                                    years.change2 = years.change2)
    cat("Constructing aggregates for Mela-Micro-Polynesia","\n")
    select.c <- seq(1, C)[is.element(country.info$namesubreg.c, 
                          c("Melanesia" , "Micronesia", "Polynesia"))]
    nameg <- "Mela-Micro-Polynesia"
    res.aggregate[[nameg]] <- InternalGetAggregates(W.Lc.t = W.Lc.t, select.c = select.c, 
                                                    dir.traj = output.dir.countrytrajectories,
                                                    years.change = years.change,
                                                    years.change2 = years.change2)
    # Add subregs
    for (subreg in 1:region.info$n.subreg){
      cat("Constructing aggregates for subregion", subreg, "(out of", region.info$n.subreg, "subregions)","\n")
      select.c <- seq(1, C)[country.info$subreg.c==subreg]
      nameg <- region.info$name.subreg[subreg]
      res.aggregate[[nameg]] <- InternalGetAggregates(W.Lc.t = W.Lc.t, select.c = select.c, 
                                             dir.traj = output.dir.countrytrajectories,
                                             years.change = years.change,
                                             years.change2 = years.change2)
    }  
    # For regs
    for (reg in 1:region.info$n.reg){
      cat("Constructing aggregates for region", reg, "(out of", region.info$n.reg, "regions)","\n")
      select.c <- seq(1, C)[country.info$reg.c==reg]
      nameg <- region.info$name.reg[reg]
      res.aggregate[[nameg]] <- InternalGetAggregates(W.Lc.t = W.Lc.t, select.c = select.c, 
                                                      dir.traj = output.dir.countrytrajectories,
                                                      years.change = years.change,
                                                      years.change2 = years.change2)
    }
    # World
    cat("Constructing aggregates for the world","\n")
    select.c <- seq(1, C)
    nameg <- "World"
    res.aggregate[[nameg]] <- InternalGetAggregates(W.Lc.t = W.Lc.t, select.c = select.c, 
                                                    dir.traj = output.dir.countrytrajectories,
                                                    years.change = years.change,
                                                    years.change2 = years.change2)    
  } else {
    # read alternative grouping (length m, need iso and groupname)
    group.data <- read.csv(file = file.aggregates, stringsAsFactors = FALSE) 
    # this file needs two columns (any other columns are ignored)
    # iso and groupname
    if (is.null(group.data$iso.country) | is.null(group.data$groupname) | is.null(group.data$iso.group)){
      cat("The csv file provided does not contain column(s) iso.country and/or groupname and/or iso.group!", "\n")
      cat("Fix that (note that missing values in iso.group column are allowed).", "\n")
    } else {
      groupname.m <- group.data$groupname
      iso.m <- group.data$iso.country
      groupnames <- unique(groupname.m)
      G <- length(groupnames)
      iso.g  <- rep(NA, G)
      W.gt <- matrix(NA, G, length(W.Lc.t))
      for (g in 1:G){
        nameg <- paste(groupnames[g])
        cat("Constructing aggregates for", nameg, "(", G, "groups in total)","\n")
        select.c <- seq(1, C)[is.element(country.info$iso.c, 
                                                      iso.m[groupname.m==nameg])]
        iso.g[g] <- group.data$iso.group[select.c[1]]
        res.aggregate[[nameg]] <- InternalGetAggregates(W.Lc.t = W.Lc.t, select.c = select.c, 
                                                        dir.traj = output.dir.countrytrajectories,
                                                        years.change = years.change,
                                                        years.change2 = years.change2)
      }
    }
  }

  W.Lg.t <- lapply(res.aggregate, function(l) l$W.t)
  CIprop.Lg.Lcat.qt <- lapply(res.aggregate, function(l) l$CIprop.Lcat.qt)
  CIratio.Lg.Lcat.qt <- lapply(res.aggregate, function(l) l$CIratio.Lcat.qt)
  changeprop.Lg.Lcat.Ti <- lapply(res.aggregate, function(l) l$changeprop.Lcat.Ti)
  changecount.Lg.Lcat.Ti <- lapply(res.aggregate, function(l) l$changecount.Lcat.Ti)

  # add counts for levels:
  CIcount.Lg.Lcat.qt <-
    FromCIpropToCIcount(CIprop.Lg.Lcat.qt = CIprop.Lg.Lcat.qt, W.Lg.t = W.Lg.t)
  
  ##value<< Object of class \code{\link{Results}}, 
  ## either for all subregions, regions and the world (UNDP aggregates),
  ## or for alternative groupings.
  return(list(iso.g = iso.g, ##<< Iso codes, but note that aggregate names (never missing) are used to name the results lists.
              CIprop.Lg.Lcat.qt = CIprop.Lg.Lcat.qt, ##<< Proportions for Total, Traditional, Modern, Unmet, TotalPlusUnmet, TradPlusUnmet. 
              ## Percentages are included as names in the \code{.qt}-matrix without percentage signs.
              CIratio.Lg.Lcat.qt = CIratio.Lg.Lcat.qt,##<<Ratios Met Demand, Z (unmet/none) and modern/total (R). R and Z might not be included for aggregates.
              CIcount.Lg.Lcat.qt = CIcount.Lg.Lcat.qt,##<< Counts for same categories as under proportions.
              changeprop.Lg.Lcat.Ti = changeprop.Lg.Lcat.Ti,##<< Changes in proportions for \code{T = years.change}, \code{i} refers to percentiles and PPPC (the posterior probability of a positive change).
              changecount.Lg.Lcat.Ti = changecount.Lg.Lcat.Ti,##<< Changes in counts, same notation as for proportions.
              W.Lg.t = W.Lg.t##<< Number of MWRA
              ))
}
#-------------------------------------------------------------------------
SelectCountriesInRegion <- function(# Select countries in a particular region
  ### Select countries in a particular region
  region.name, 
  country.info,
  region.info
) {
  if (region.name %in% region.info$name.subreg) {
    select.c <- country.info$namesubreg.c == region.name
  } else if (region.name %in% region.info$name.reg) {
    select.c <- country.info$namereg.c == region.name
  } else if (region.name == "World") {
    select.c <- rep(TRUE, length(country.info$name.c))
  } else if (region.name == "Developed countries") {
    select.c <- country.info$dev.c == "Rich"
  } else if (region.name == "Developing countries") {
    select.c <- country.info$dev.c != "Rich"
  } else if (region.name == "Developing (excl. China)") {
    select.c <- country.info$dev.c != "Rich" & country.info$name.c != "China"
  } else if (region.name == "Mela-Micro-Polynesia") {
    select.c <- is.element(country.info$namesubreg.c, 
                           c("Melanesia" , "Micronesia", "Polynesia"))
  } else {
    stop(paste0(region.name, " does not exist!"))
  }
  return(select.c)
}
#-------------------------------------------------------------------------
# The End!
