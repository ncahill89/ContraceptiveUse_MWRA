#----------------------------------------------------------------------------------
# Leontine Alkema and Jin Rou New
# F_plotCIstime.R
# Functions to plot CIs only for any object with CIprop and CIcounts
# and country data
#----------------------------------------------------------------------------------
PlotDataAndEstimates <- function (# Create overview country/aggregate plots 
  ### Create overview plots of estimates of proportions/counts over time 
  ### for either 
  ### all countries, or a set of aggregates.
  ### Note that all input arguments with default setting NULL are optional 
  ### (e.g. data are not plotted if data.raw=NULL).
  data.raw = NULL, ##<< Non-NULL only for country plots: List with \code{data},
  ## \code{country.info} and \code{region.info}.
  start.year = NULL, ##<< Which years are plotted? Defaults to estimation years used in \code{CI.Lg.Lcat.qt}, or 1990 otherwise.
  end.year = NULL, ##<< Which years are plotted? Defaults to estimation years used in \code{CI.Lg.Lcat.qt}, or 1990 otherwise.
  CI.Lg.Lcat.qt = NULL, ##<< Object from class \code{CI.Lg.Lcat.qt} (optional), either a proportion or a count (see next).
  plot.prop = TRUE, ##<< Are the CIs proportions or counts?
  ymin.at.0 = FALSE, ##<< Set lower bound on y-axis at 0?
  ymax.at.100 = FALSE, ##<< Set upper bound on y-axis at percent = 100%? Only applies if plot.prop is TRUE.
  add.info = TRUE, ##<< Add information related to misclassification or different populations.
  par.ciq = NULL, ##<< Add info on logistic parameters (optional, relevant for country proportion plots only)
  CIstar.Lg.Lcat.qt = NULL, ##<< Systematic curves, object from class \code{CI.Lg.Lcat.qt} (optional, relevant for country proportion plots only)
  CIratio.Lg.Lcat.qt = NULL, ##<< Systematic curves, object from class \code{CI.Lg.Lcat.qt} (optional, relevant for proportion plots only)
  meanpropratio.Lg.Lcat.t = NULL, ##<< Mean (as opposed to median) proportions of indicators (optional, relevant for country proportion plots only).
  validation = FALSE, ##<< Logical, if TRUE, plot observation from test sets in grey.
  getj.test.k = NULL, ##<< Used if \code{validation}, indices of test sets for observations on
  ## traditional and modern use.
  getj.test.unmet.k = NULL,##<< Used if \code{validation}, indices of test set for observations
  ## on unmet need.
  fig.name = NULL, ##<< Used for overview plots only:
  ## If NULL, plot appears in R, else it is saved as fig.name.
  ## Not used if individual country plots are requested using 
  ## \code{ind.country.overviewplot} or 
  ## \code{ind.country.indplot} below.
  select.c = NULL, ##<< For country plots only, indices of countries to plot
  ind.country.overviewplot = FALSE, ##<< Logical: for each country, overview plot 
  ## saved to \code{figdir.indcountries}? (in TIFF format)
  ind.country.indplot = FALSE, ##<< Logical: for each country, 
  ## individual plots saved to \code{figdir.indcountries}? (NOT USED YET)
  figdir.indcountries = NULL, ##<< Directory to store country plots, if NULL, directory \code{fig/ind.country.plots} 
  ## is created in current working directory.
  categories.to.plot = NULL, ##<< Selected names of categories to plot. If NULL, all available categories are plotted.
  run.name = NULL, ##<< Run name, used only to add to directory name when plotting results for individual countries.
  CI2.Lg.Lcat.qt = NULL, ##<< Add a second set of CIs (optional)
  CIratio2.Lg.Lcat.qt = NULL, ##<< Add a second set of CIs (optional)
  CI3.Lg.Lcat.qt = NULL, ##<< Add a third set of CIs (optional)
  CIratio3.Lg.Lcat.qt = NULL, ##<< Add a third set of CIs (optional)
  CI4.Lg.Lcat.qt = NULL, ##<< Add a fourth set of CIs (optional)
  CIratio4.Lg.Lcat.qt = NULL, ##<< Add a fourth set of CIs (optional)
  name.dir1 = NULL, ##<< Used if CIs2 are added, name used for CIs
  name.dir2 = NULL, ##<< Used if CIs2 are added, name used for CIs2
  name.dir3 = NULL, ##<< Used if CIs3 are added, name used for CIs3
  name.dir4 = NULL, ##<< Used if CIs4 are added, name used for CIs4
  cex.adj.factor = 1, ## Factor to adjust size of plots by
  #  plot.blue.line = FALSE, ##<< add main trends? only if next one not null
  shiny = FALSE, ##<< \code{TRUE} if plot function is used for Shiny FPET app
  turn.legend.off=FALSE,
  use.smaller.points=FALSE
){
  if (!is.null(data.raw)){ # for country plots only
    data <- data.raw$data
    country.info <- data.raw$country.info
    region.info <- data.raw$region.info
    mainlong.g <- paste0(country.info$name.c #, " (", country.info$namesubreg.c, ")"
                         )
    # for validation exercise unmet need, only plot countries where data were left out
  } else {
    mainlong.g <- names(CI.Lg.Lcat.qt)
  }
  name.g <- mainlong.g
  G <- length(name.g)
  
  fig.name.years <- ifelse(!is.null(start.year) | !is.null(end.year),
                           paste0("_",
                                  ifelse(!is.null(start.year), paste0("from", floor(start.year)), ""),
                                  ifelse(!is.null(end.year), paste0("to", floor(end.year)), "")),
                           "")
  if (is.null(start.year)){
    if (!is.null(CI.Lg.Lcat.qt)){
      est.years <- as.numeric(names(CI.Lg.Lcat.qt[[1]][[1]][1,]))
      start.year <- est.years[1]
    } else {
      start.year <- 1990.5
    } 
  }
  if (is.null(end.year)){
    if (!is.null(CI.Lg.Lcat.qt)){
      # LAchange20140610: add estyears here again (in case start year was not NULL)
      est.years <- as.numeric(names(CI.Lg.Lcat.qt[[1]][[1]][1,]))
      end.year <- est.years[length(est.years)]
    } else {
      end.year <- 2015.5      
    } 
  }
  xmin <- start.year #min(years.i, start.year, na.rm = T)
  xmax <- end.year
  
  cats.all <- c("Total", "Modern", "Traditional", "Modern/Total", 
                "TradPlusUnmet", "TotalPlusUnmet", "Met Demand with Modern Methods",
                "Unmet", "TotalPlusUnmet", "Met Demand")
  # ADHOC CHANGE JR
  #   axisnamecats.all <- c("CP (any)", "CP (modern)", "CP (traditional)", "Ratio (modern/any method)", 
  #                         "Unmet need for\nmodern methods", "Total demand", "Demand satisfied\nwith modern methods",
  #                         "Unmet need", "Total demand", "Demand satisfied")
  axisnamecats.all <- c("Total\ncontraceptive prevalence", "Modern\ncontraceptive prevalence", 
                        "Traditional\ncontraceptive prevalence", 
                        "Ratio (modern CP/total CP)", 
                        "Unmet need for\nmodern methods", "Total demand", 
                        "Demand satisfied\nwith modern methods",
                        "Unmet need", "Total demand", "Demand satisfied")
  if (!is.null(CIratio.Lg.Lcat.qt)){
    select.cats <- seq_along(cats.all)
    cats.from.ratio <- c(F,F,F,T,F,F,T,F,F,T)
  } else {
    select.cats <- c(1, 2, 3, 5, 8, 9)
    cats.from.ratio <- c(F,F,F,F,F,F)
  }
  cats <- cats.all[select.cats]
  axisnamecats <- axisnamecats.all[select.cats]
  
  if (is.null(categories.to.plot))
    categories.to.plot <- cats
  if (is.null(CIratio.Lg.Lcat.qt))
    categories.to.plot <- categories.to.plot[
      !(categories.to.plot %in% c("Modern/Total", "Met Demand with Modern Methods", "Met Demand"))]
  cats.select <- match(categories.to.plot, cats)
  if(turn.legend.off){
  nplots <- length(categories.to.plot) + ifelse(is.null(data.raw), 0, 0) # + 1 for legend if plotting data
  }
  if(!turn.legend.off){
    nplots <- length(categories.to.plot) + ifelse(is.null(data.raw), 0, 1) # + 1 for legend if plotting data
  }
  ncols <- min(ifelse(plot.prop, 6, 3), nplots)
  nrows <- max(ceiling(nplots/ncols))
  
  if (ind.country.overviewplot | ind.country.indplot){
    if (is.null(figdir.indcountries)){
      figdir.indcountries <- file.path(getwd(), "fig", paste0("ind.country.plots_", run.name))
      dir.create(file.path(getwd(), "fig"), showWarnings = FALSE) 
    }
    dir.create(figdir.indcountries, showWarnings = FALSE) 
    if (ind.country.overviewplot & ind.country.indplot){
      ## details<< If \code{ind.country.indplot} and \code{ind.country.indplot} are both TRUE, 
      ## currently, \code{ind.country.indplot} is used.
      ind.country.overviewplot <- FALSE
    }
  } else {
    if (!is.null(fig.name)) {
      pdf(fig.name, width = ncols*5.25, height = nrows*5.5 + 1)
      # pdf(fig.name, width = 25.5, height = 21)
    }
  }
  
  # cats.from.ratio refers to taking it from CI or CIratio
  if (!is.null(select.c)) {
    gseq <- select.c
  } else {
    gseq <- 1:G
  }
  
  for (g in gseq) {
    # ##details<< If \code{ind.country.overviewplot}, plots saved in \code{figdir.indcountries}
    # if (ind.country.overviewplot) 
    #   pdf(file.path(figdir.indcountries, paste0(country.info$name.c[g], fig.name.years, ".pdf")), width = 21, height = 12)
    if (ind.country.overviewplot) 
      tiff(file.path(figdir.indcountries, paste0(name.g[g], fig.name.years, ".tif")), width = 21, height = 12,
           pointsize=6,res=300, bg="white",compression="lzw")
    if (!is.null(data.raw)) { # for country plots only
      select <- seq(1, length(data[,1]))[is.element(data$iso.j, country.info$iso.c[g])]
      years.i <- data$years.j[select]
      include.SS <- ifelse(any(data$source.j[select]=="SS"), TRUE, FALSE)
      
      check <- c("DHS", "MICS", "National survey", 
                 "Other survey", "Service statistics",
                 "Subpopulation", "", 
                 "+:  Higher contraceptive use", # + for age pos and pos bias because of nonpreg etc
                 "-:  Lower contraceptive use",
                 "A:  Other age group", 
                 "F:  Folk methods included", 
                 "S-: Sterilization included", 
                 "S+: Sterilization excluded","",
                 "Married women", 
                 "Sexually active women", "Ever married/All women",
                 "Both sexes and husband/wives")
      
      data.type.included <- c(sapply(c("DHS", "MICS", "NS", "Other", "SS"),
                                     CheckIfAny, data$source.j[select]),
                              any(data$geo.j[select] != ""),
                              any((data$age.cat.j[select] == "+") | (data$posbias.j[select] != "")),
                              CheckIfAny("-", data$age.cat.j[select]),
                              CheckIfAny("?", data$age.cat.j[select]),
                              any(data$folkbias.j[select] != ""),
                              CheckIfAny("-", data$mod.bias.j[select]),
                              CheckIfAny("+", data$mod.bias.j[select]),
                              CheckIfAny("MW", data$poptype.j[select]),
                              CheckIfAny("SA", data$poptype.j[select]),
                              CheckIfAny(c("EM", "AL"), data$poptype.j[select]),
                              CheckIfAny(c("BS", "HW"), data$poptype.j[select]))
      if (add.info) {
        pch.i <- ifelse(data$source.j[select] == "SS", 19,
                        ifelse(data$poptype.j[select]=="SA", 24,
                               ifelse(data$poptype.j[select]=="EM" | data$poptype.j[select]=="AL", 25,
                                      ifelse(data$poptype.j[select]=="HW" | data$poptype.j[select]=="BS", 22, 21))))
      } else {
        pch.i <- ifelse(data$source.j[select] == "SS", 19,
                        ifelse(data$poptype.j[select] == "MW" & data$geo.j[select] == "" &
                                 data$age.cat.j[select] == "0" & data$posbias.j[select] == "" &
                                 data$mod.bias.j[select] == "" & data$folkbias.j[select] == "", 21, # standard 
                               22))
      }
      cex.i <- ifelse(data$source.j[select]=="SS", 1.5,
                      ifelse(add.info, 6, 4))
      col.i <- ifelse(data$source.j[select]=="DHS","red",
                      ifelse(data$source.j[select]=="MICS","green",
                             ifelse(data$source.j[select]=="NS","blue",
                                    ifelse(data$source.j[select]=="SS","saddlebrown",1))))
      bg.i <- bg.unmet.i <- rep("white", length(select)) 
      # do unmet seperately, because validation can be different
      if (!is.null(getj.test.k )){
        # validation exercise, use grey for obs that were left out  
        bg.i <- ifelse(is.element(select, getj.test.k), "darkgrey", bg.i)
      }
      if (!is.null(getj.test.unmet.k )){ # for unmet only exercise
        # note that other props show up as if there were left out as well
        bg.i <- ifelse(is.element(select, getj.test.unmet.k), "darkgrey", bg.unmet.i)
      } 
      # end data.raw loop
    } else {
      include.SS <- FALSE
    }
    
    if (shiny & length(cats.select) == 10) {
      nf <- layout(rbind(rep(1, ncols),
                         (1:ncols) + 1,
                         c(rep(ncols + 2, ncols - 1), ncols + 3),
                         (1:ncols) + ncols + 3,
                         c(rep(2*ncols + 4, ncols -1), 2*ncols + 5),
                         (1:ncols) + 2*ncols + 5),
                   widths = rep(1.5, ncols), 
                   heights = c(0.25, 1.5, 0.2, 1.5, 0.2, 1.5), TRUE)
    } else {
      nf <- layout(rbind(rep(1, ncols),
                         matrix(seq(2, nrows*ncols+1), nrows, ncols, byrow = TRUE)),
                   widths = rep(1.5, ncols), 
                   heights = c(0.25, rep(1.5, nrows)), TRUE)
    }
    
    # ADHOC CHANGE JR ######## single plot
    #     nplots <- length(categories.to.plot) + 1
    #     ncols <- nplots
    #     nrows <- 1
    #     nf <- layout(rbind(rep(1, ncols),
    #                        matrix(seq(2, nrows*ncols+1), nrows, ncols, byrow = TRUE)),
    #                  widths = rep(1.5, ncols), 
    #                  heights = c(0.25, rep(1.5, nrows)), TRUE)
    
    #     # ADHOC CHANGE JR ######## all plots
    #     nplots <- length(categories.to.plot) + 1
    #     ncols <- min(nplots, ifelse(plot.prop, 4, 3))
    #     nrows <- max(ceiling(nplots/ncols))
    #     nf <- layout(rbind(rep(1, ncols),
    #                        matrix(seq(2, nrows*ncols+1), nrows, ncols, byrow = TRUE)),
    #                  widths = rep(1.5, ncols), 
    #                  heights = c(0.25, rep(1.5, nrows)), TRUE)
    
    # layout.show(nf)
    if(turn.legend.off){
    InternalPlotTitle(title = " ", position = "left", cex = 3.2*cex.adj.factor)
    legend.select.comparison <- c(TRUE, !is.null(CI2.Lg.Lcat.qt[[g]]), 
                                  !is.null(CI3.Lg.Lcat.qt[[g]]), !is.null(CI4.Lg.Lcat.qt[[g]]))
    if (sum(legend.select.comparison) > 1)
      legend("center", legend = c(name.dir1, name.dir2, name.dir3, name.dir4)[legend.select.comparison], 
             col = c("red","blue","green","#984EA3")[legend.select.comparison], lwd = 3, cex = 1.5*cex.adj.factor,bty="n")
    }
    
    if(!turn.legend.off){
      InternalPlotTitle(title = paste(mainlong.g[g], "   "), position = "center", cex = 3.2*cex.adj.factor)
      legend.select.comparison <- c(TRUE, !is.null(CI2.Lg.Lcat.qt[[g]]), 
                                    !is.null(CI3.Lg.Lcat.qt[[g]]), !is.null(CI4.Lg.Lcat.qt[[g]]))
      if (sum(legend.select.comparison) > 1)
        legend("right", legend = c(name.dir1, name.dir2, name.dir3, name.dir4)[legend.select.comparison], 
               col = c("red","blue","green","#984EA3")[legend.select.comparison], lwd = 3, cex = 1.5*cex.adj.factor)
    }
    
    # par(mar = c(7,6,6,1), cex.main = 2, cex.axis = 2, cex.lab = 2)
    for (cat in (1:length(cats))[cats.select]) {   
      if (!is.null(CI.Lg.Lcat.qt)) {
        if (!is.null(CIstar.Lg.Lcat.qt[[g]][[cats[cat]]])) {
          CIstar.qt <- CIstar.Lg.Lcat.qt[[g]][[cats[cat]]]
        } else {
          CIstar.qt <- NULL
        }
        if (!cats.from.ratio[cat]) {
          est.years <- as.numeric(names(CI.Lg.Lcat.qt[[gseq[1]]][[1]][1,]))
          CI.qt <- CI.Lg.Lcat.qt[[g]][[cats[cat]]]
          CI.qt[, est.years < start.year | est.years > end.year] <- NA
          if (!is.null(CI2.Lg.Lcat.qt)) {
            est.years2 <- as.numeric(names(CI2.Lg.Lcat.qt[[gseq[1]]][[1]][1, ]))
            CI2.qt <- CI2.Lg.Lcat.qt[[g]][[cats[cat]]]
            CI2.qt[, est.years2 < start.year | est.years2 > end.year] <- NA
          } else {
            CI2.qt <- NULL
          }
          if (!is.null(CI3.Lg.Lcat.qt)) {
            est.years3 <- as.numeric(names(CI3.Lg.Lcat.qt[[gseq[1]]][[1]][1, ]))
            CI3.qt <- CI3.Lg.Lcat.qt[[g]][[cats[cat]]]
            CI3.qt[, est.years3 < start.year | est.years3 > end.year] <- NA
          } else {
            CI3.qt <- NULL
          }
          if (!is.null(CI4.Lg.Lcat.qt)) {
            est.years4 <- as.numeric(names(CI4.Lg.Lcat.qt[[gseq[1]]][[1]][1,]))
            CI4.qt <- CI4.Lg.Lcat.qt[[g]][[cats[cat]]]
            CI4.qt[, est.years4 < start.year | est.years4 > end.year] <- NA
          } else {
            CI4.qt <- NULL
          }
        } else {
          est.years <- as.numeric(names(CIratio.Lg.Lcat.qt[[gseq[1]]][[1]][1,]))
          CI.qt <- CIratio.Lg.Lcat.qt[[g]][[cats[cat]]]
          CI.qt[, est.years < start.year | est.years > end.year] <- NA
          if (!is.null(CIratio2.Lg.Lcat.qt)) {
            est.years2 <- as.numeric(names(CIratio2.Lg.Lcat.qt[[gseq[1]]][[1]][1,]))
            CI2.qt <- CIratio2.Lg.Lcat.qt[[g]][[cats[cat]]]
            CI2.qt[, est.years2 < start.year | est.years2 > end.year] <- NA
          } else {
            CI2.qt <- NULL
          }
          if (!is.null(CIratio3.Lg.Lcat.qt)) {
            est.years3 <- as.numeric(names(CIratio3.Lg.Lcat.qt[[gseq[1]]][[1]][1,]))
            CI3.qt <- CIratio3.Lg.Lcat.qt[[g]][[cats[cat]]]
            CI3.qt[, est.years3 < start.year | est.years3 > end.year] <- NA
          } else {
            CI3.qt <- NULL
          }
          if (!is.null(CIratio4.Lg.Lcat.qt)) {
            est.years4 <- as.numeric(names(CIratio4.Lg.Lcat.qt[[gseq[1]]][[1]][1,]))
            CI4.qt <- CIratio4.Lg.Lcat.qt[[g]][[cats[cat]]]
            CI4.qt[, est.years4 < start.year | est.years4 > end.year] <- NA
          } else {
            CI4.qt <- NULL
          }
        }
        if (!is.null(meanpropratio.Lg.Lcat.t)) {
          mean.qt <- matrix(NA, 5, length(est.years))
          if (!is.null(meanpropratio.Lg.Lcat.t[[g]][[cats[cat]]])) {
            mean.qt[3, ] <- meanpropratio.Lg.Lcat.t[[g]][[cats[cat]]]
            mean.qt[3, est.years < start.year | est.years > end.year] <- NA
          }
        } else {
          mean.qt <- NULL
        }
      }
      name.cat <- cats[cat]
      if (ind.country.indplot) {
        par(mfrow = c(1,1))
        pdf(file.path(figdir.indcountries, 
                      paste0(name.g[g], "_total", fig.name.years, ".pdf")), width = 12, height = 12)
        main.extra <- name.g[g]
      } else {
        main.extra <- ifelse(cat==2, mainlong.g[g], "")
      }
      
      # actual plotting for that category
      par(mar = c(6, ifelse(cex.adj.factor == 1, 6, 4), 5, 1), 
          cex.main = 2.5*cex.adj.factor, cex.axis = 2*cex.adj.factor, cex.lab = 2*cex.adj.factor)
      if (!is.null(data.raw)) {
        props.j <- NULL
        trad <- modern <- FALSE
        if (name.cat=="Total") props.j <- data$props.tot.j
        if (name.cat=="Modern/Total") props.j <- data$props.modern.j/data$props.tot.j
        if (name.cat=="TradPlusUnmet") props.j <- data$props.trad.j + data$props.unmet.j
        if (name.cat=="TotalPlusUnmet") props.j <- data$props.tot.j + data$props.unmet.j
        if (name.cat=="Met Demand with Modern Methods") props.j <- data$props.modern.j/(data$props.tot.j + data$props.unmet.j)
        if (name.cat=="Unmet") props.j <- data$props.unmet.j
        if (name.cat=="Met Demand") props.j <- data$props.tot.j/(data$props.tot.j + data$props.unmet.j)
        if (name.cat=="Modern") {
          props.j <- data$props.modern.j
          modern <- TRUE
        }
        if (name.cat=="Traditional") {
          props.j <- data$props.trad.j
          trad <- TRUE
        }
        if (sum(!is.na(props.j[select])) > 0) {
          data.props <- props.j[select]
          if (xmin > 1989)
            data.props <- data.props[years.i > 1989]
          data.props <- data.props[!is.na(data.props)]
        } else {
          data.props <- NULL
        }
      } else {
        data.props <- NULL
      }
      yall <- c(CI.qt, CI2.qt, CI3.qt, CI4.qt, data.props)
      ymin <- ifelse(ymin.at.0, 0, ifelse(min(yall, na.rm = T) < 0.05, 0, 
                                          min(yall, na.rm = T)-0.2*diff(range(yall, na.rm = T))))
      ymax <- ifelse(plot.prop & ymax.at.100, 1, 
                     max(yall, na.rm = T)+0.2*diff(range(yall, na.rm = T)))
      InternalPlotEmpty(ylab = ifelse(plot.prop,"", "Count"),
                        main = axisnamecats[cat], 
                        xlim = c(xmin, xmax),
                        ylim = c(ymin, ymax),
                        plot.prop = plot.prop) 
      # NOTE: might want to change xlim for incl earlier data
      if (!is.null(CI2.Lg.Lcat.qt) | !is.null(CI3.Lg.Lcat.qt) | !is.null(CI4.Lg.Lcat.qt)) {
        if (!is.null(CI4.qt))
          InternalPlotCIs(CIs.qt = CI4.qt, col.median = "#984EA3", seq.years = est.years4,
                          col95 = "#984EA332") #purple 
        if (!is.null(CI3.qt))
          InternalPlotCIs(CIs.qt = CI3.qt, col.median = "green", seq.years = est.years3,
                          col95 = "#00FF0032")
        if (!is.null(CI2.qt))
          InternalPlotCIs(CIs.qt = CI2.qt, col.median = "blue", seq.years = est.years2,
                          col95 = "#0000FF32")
        InternalPlotCIs(CIs.qt = CI.qt, col.median = "red", seq.years = est.years,
                        col95 = "#FF000020")       
      } else {
        if (!is.null(CI.Lg.Lcat.qt)){
          InternalPlotCIs(CIs.qt = CI.qt, col.median = 1, #col_pxmr[1],
                          seq.years =  est.years, col95 = "#0000FF30", 
                          CIs.star.qt = CIstar.qt)
        }
      }
      if (!is.null(meanpropratio.Lg.Lcat.t)) {
        InternalPlotCIs(CIs.qt = mean.qt, col.median = "red", lty = 2, seq.years =  est.years)
      }
      if (!is.null(data.raw)) {
        if (sum(!is.na(props.j[select])) > 0) {
          if (ind.country.overviewplot) {
            cexuse.i <- cex.i*0.6
          } else {
            cexuse.i <- cex.i
          }
          
          if(!use.smaller.points){
           InternalPlotData(props.i = props.j[select], years.i, 
                            col.i = col.i, pch.i = pch.i, cex.i = cexuse.i,
                            add.info = add.info, data = data, select = select, bg.i = bg.i,
                            trad = trad, modern = modern, cex.adj.factor = cex.adj.factor)
          }
          if(use.smaller.points){
          InternalPlotData(props.i = props.j[select], years.i, 
                           col.i = col.i, pch.i = 19, cex.i = 1.5,
                           add.info = FALSE, data = data, select = select, bg.i = bg.i,
                           trad = trad, modern = modern, cex.adj.factor = cex.adj.factor)
          }
          if (min(years.i[!is.na(props.j[select])]) < 1990 & xmin > 1989) {
            text(xmin-1, ymax, "*Data incl. before 1990", cex = 1.7*cex.adj.factor, pos = 4)
          }
        }
      }
      if (!is.null(par.ciq) & name.cat=="Total") 
        InternalPlotParInfoTot(par.ciq[g,,], cex.adj.factor = cex.adj.factor)
      if (!is.null(par.ciq) & name.cat=="Modern/Total") 
        InternalPlotParInfoRat(par.ciq[g,,], cex.adj.factor = cex.adj.factor)
      if (ind.country.indplot) 
        dev.off()
      if (shiny & length(cats.select) == 10 & name.cat=="Modern/Total") {
        InternalPlotTitle(title = "Unmet need and demand for modern methods", 
                          position = "center", cex = 2.8*cex.adj.factor)
        InternalPlotNull()
      }
      if (length(cats.select) == 10 & name.cat == "Met Demand with Modern Methods")
        InternalPlotNull()
      if (shiny & length(cats.select) == 10 & name.cat == "Met Demand with Modern Methods") {
        InternalPlotTitle(title = "Unmet need and demand for any method", 
                          position = "center", cex = 2.8*cex.adj.factor)
        InternalPlotNull()
      }
    }
    if (!ind.country.indplot & !is.null(data.raw)){
      if(!turn.legend.off){
       InternalPlotDataLegend(TIFF = ind.country.overviewplot, 
                              add.info = add.info, data.type.included = data.type.included,
                              cex.adj.factor = cex.adj.factor)
      }
  }
    if (ind.country.overviewplot) 
      dev.off()
  } # end country loop
  
  if (ind.country.indplot) {
    pdf(file.path(figdir.indcountries, paste0("legend", fig.name.years, ".pdf")), width = 12, height = 12) # change JR, 20140418
    InternalPlotDataLegend(add.info = add.info, include.SS = include.SS, cex.adj.factor = cex.adj.factor)
    dev.off()
  }
  if (!is.null(fig.name) & !(ind.country.overviewplot | ind.country.indplot))
    dev.off()
  return(invisible())
}
#----------------------------------------------------------------------
InternalPlotTitle <- function(
  title, ##<< Title to plot
  position = "center", ##<< Relative position of title
  cex = 1 ##<< Size adjustment factor
) {
  par(mar = c(0,0,0,0))
  plot(1~1, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
  legend(position, title, text.font = 2, cex = cex, bty = "n")
}
#----------------------------------------------------------------------
InternalPlotNull <- function() {# Make a completely empty plot
  ### Make a completely empty plot
  plot(1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
}
#----------------------------------------------------------------------------------
InternalPlotEmpty <- function(# Start empty plot with labeled axes
  ###F Start empty plot with labeled axes
  main = "", xlab = "", ylab = "", 
  xlim = c(1990, 2015), ylim = c(0, 1),
  plot.prop = FALSE
){
  plot(1, type = "n", xlim = xlim, ylim = ylim,
       main = "", xlab = xlab, ylab = ylab,
       xaxt = "n", yaxt = "n")
  mtext(main, line = 1, cex = 1.6, font = 2)
  abline(h = 0)
  #   if (ymax > 1.1){ # counts
  #     axis(2, at = round(seq(0,ymax,length.out = 10)), labels = round(seq(0,ymax/1000,length.out = 10)), las = 2)
  # #  } else {
  # #    if (ymax > 1.1){ # counts
  # #      axis(2, at = round(seq(0,ymax,length.out = 10)), labels = round(seq(0,ymax,length.out = 10)), las = 2)
  #       } else {
  #         if (ymax > 0.2) axis(2, at = seq(0,1,0.1), labels = seq(0,1,0.1)*100, las = 2)
  #         if (ymax <= 0.2 & ymax >0.1) axis(2, at = seq(0,ymax,0.05), labels = seq(0,ymax,0.05)*100, las = 2)
  #         if (ymax < 0.1) axis(2, at = seq(0,ymax,0.02), labels = seq(0,ymax,0.02)*100, las = 2)
  #       }
  #  # }
  # LAchange20140610
  axis(1, at = seq(5*floor(min(xlim)/5), 5*ceiling(max(xlim)/5),5), las = 3)
  x <- axis(2, labels = FALSE)
  if (plot.prop){
    axis(2, at = x, labels = paste(100*x,"%", sep = ""), las = 2)
  } else {
    axis(2)
  }
}
#----------------------------------------------------------------------------------
InternalPlotCIs <- function(# Add CIs to a plot
  ### Add CIs to a plot
  CIs.qt, ##<< q has to refer to (2.5, 10, 50, 80, 97.5)th percentiles
  CIs.star.qt = NULL, ##<< Add main trend if non-NULL
  col.median = 1, ##<< Color for line to represent median, lines not added if \code{NULL}
  col80 = col.median, ##<< Color for lines to represent 80% CI, lines not added if \code{NULL}
  lty = 1,
  seq.years, ##<< Estimation years
  col95 = "#0000FF30" ##<< Color for shaded area to represent 95% CI, area not added if \code{NULL}
){
  # adds areas with col95 for 95% CI
  # Note: c(1,5) from CI are used!!!
  nyears <- length(seq.years)
  CI.low.t <- CIs.qt[1,]
  CI.up.t <- CIs.qt[5,]
  if (!is.null(col95)) {
    for (t in 2:nyears) {
      polygon(c(seq.years[t-1], seq.years[t-1], seq.years[t], seq.years[t],seq.years[t-1]),
              c(CI.low.t[t-1], CI.up.t[t-1], CI.up.t[t], CI.low.t[t], CI.low.t[t-1]),
              col=col95, border = NA)
    }
  }
  if (!is.null(col80)) {
    for (q in c(2,4))
      lines(CIs.qt[q, ] ~ seq.years, type = "l",  lwd = 2, lty = lty, col = col80)
  }
  if (!is.null(col.median)) {
    lines(CIs.qt[3, ] ~ seq.years, type = "l",  lwd = 2, lty = lty, col = col.median)
  }
  if (!is.null(CIs.star.qt))
    lines(CIs.star.qt[3,] ~ seq.years, type = "l", col = "blue", lwd = 3)
}
#----------------------------------------------------------------------------------
InternalPlotData <- function(# Add data to a plot
  ### Add data to a plot
  props.i, ##<< i refers to observation index (selected in indices in vector of length j)
  years.i,
  lwd.i = 3, bg.i = NULL, col.i =1, pch.i = 19, cex.i = 1,
  add.info = FALSE, ##<< Logical: do you want to plot the different symbols to represent source, biases etc?
  data = NULL, ##<< used if \code{add.info}, include object from class \code{data} 
  select = NULL, ##<< used if \code{add.info}, info only added for indices in select
  trad = FALSE, ##<< used if \code{add.info}; logical, are you plotting traditional data? 
  modern = FALSE, ##<< used if \code{add.info}; logical, are you plotting modern data?
  cex.adj.factor = 1
) {
  if (sum(!is.na(props.i)) != 0) {
    points(props.i ~ years.i, bg = bg.i, 
           cex = cex.i*cex.adj.factor, col = col.i, lwd = lwd.i*cex.adj.factor, pch = pch.i)
    
    if (add.info) {
      points(props.i[data$geo.j[select] !=""] ~ years.i[data$geo.j[select] !=""], 
             col = "orange", cex = 1.5*cex.adj.factor, lwd = 5*cex.adj.factor, 
             pch = pch.i[data$geo.j[select] !=""]) 
      # age stuff
      if (sum(data$age.cat.j[select] == "?", na.rm = TRUE) !=0){
        text.label = ifelse(data$age.cat.j[select] == "?", "A", "")
        text(props.i ~ years.i, labels = text.label, col = 1, cex = 2*cex.adj.factor)
      }
      if (sum(data$age.cat.j[select] == "+", na.rm = TRUE) !=0){
        text.label = ifelse(data$age.cat.j[select] == "+", "+", "")
        text(props.i ~ years.i, labels = text.label, col = 1, cex = 2*cex.adj.factor)
      }
      if(sum(data$age.cat.j[select] == "-", na.rm = TRUE)!=0){
        text.label = ifelse( data$age.cat.j[select]=="-", "  -", "")
        text(props.i ~ years.i, labels = text.label, col = 1, cex = 2.5*cex.adj.factor)
      }  
      if (sum(data$posbias.j[select]!="", na.rm = TRUE) !=0){
        text.label = ifelse(data$posbias.j[select]!="","  +","")
        text(props.i ~ years.i, labels = text.label,col = 1, cex =2*cex.adj.factor)
      }
      if (modern) {
        if (sum(data$mod.bias.j[select]=="+", na.rm = TRUE) != 0) {
          text.label = ifelse(data$mod.bias.j[select]=="+","S+","")
          text(props.i ~ years.i, labels = text.label,col = 1, cex =2*cex.adj.factor)
        }
        if (sum(data$mod.bias.j[select]=="-", na.rm = TRUE) != 0) {
          text.label = ifelse(data$mod.bias.j[select]=="-","S-","")
          text(props.i ~ years.i, labels = text.label,col = 1, cex =2*cex.adj.factor)
        }
      }
      if (trad){
        if (sum(data$folkbias.j[select]!="", na.rm = TRUE) !=0) {
          text.label = ifelse(data$folkbias.j[select]!="", "F", "")
          text(props.i ~ years.i, labels = text.label, col = 1, cex =2*cex.adj.factor)
        }
      }
    }
  }
}
#----------------------------------------------------------------------------------
InternalPlotParInfoTot <- function(#Add information about logistic parameters to plot
  ### Add information about logistic parameters of total to plot
  par.iq, ##<< matrix with CIs from par.ciq (constructed using \code{\link{GetParInfo}} for one country
  cex.adj.factor = 1
){
  abline(h = par.iq["pmax.c", 2])
  abline(h = par.iq["pmax.c", c(1,3)], lty = 2)
  legend("bottomright", cex = 1.8*cex.adj.factor, legend = c(
    paste0("omega", " = ", round(par.iq["omega.c", 2],2), 
           " (", round(par.iq["omega.c", 3] - par.iq["omega.c", 1],2), ")"),
    paste0("p[max]", " = ", round(par.iq["pmax.c", 2],2), 
           " (", round(par.iq["pmax.c", 3] - par.iq["pmax.c", 1],2), ")"))
  )
}

#----------------------------------------------------------------------------------
InternalPlotParInfoRat <- function(#Add information about logistic parameters to plot
  ### Add information about logistic parameters of ratio modern/total to plot
  par.iq, ##<< matrix with CIs from par.ciq (constructed using \code{\link{GetParInfo}} for one country
  cex.adj.factor = 1
){
  abline(v = par.iq["RT.c", 2])
  abline(v = par.iq["RT.c", c(1,3)], lty = 2)
  abline(h = par.iq["Rmax.c", 2])
  abline(h = par.iq["Rmax.c", c(1,3)], lty = 2)
  location <- "bottomright"
  legend(location, cex = 1.8*cex.adj.factor, legend = c(
    paste0("R",expression(omega), " = ", round(par.iq["Romega.c", 2],2), 
           " (", round(par.iq["Romega.c", 3] - par.iq["Romega.c", 1],2), ")"),
    paste0("RT = ", round(par.iq["RT.c", 2]), 
           " (", round(par.iq["RT.c", 3] - par.iq["RT.c", 1]), ")"),
    paste0(expression(R[max]), " = ", round(par.iq["Rmax.c", 2],2), 
           " (", round(par.iq["Rmax.c", 3] - par.iq["Rmax.c", 1],2), ")"))
  )
}

#----------------------------------------------------------------------------------
InternalPlotDataLegend <- function(# Plot legend with details on biases etc.
  ### Plot legend with details on biases etc.
  TIFF = FALSE, ##<< TIFF (if false, it's PDF)
  add.info = TRUE, 
  data.type.included,
  cex.adj.factor = 1
){
  if (TIFF){
    par(mar = c(0,ifelse(cex.adj.factor == 1, 6, 4),0,1), 
        cex.main = 1.5*cex.adj.factor, cex.axis = 1.5*cex.adj.factor, cex.lab = 1.5*cex.adj.factor)
  } else {
    par(mar = c(0,ifelse(cex.adj.factor == 1, 6, 4),0,1), 
        cex.main = 1.5*cex.adj.factor, cex.axis = 1.5*cex.adj.factor, cex.lab = 1.5*cex.adj.factor)
  }
  if (add.info) {
    data.type.included <- c(data.type.included[1:6], TRUE, data.type.included[7:12], TRUE,
                            data.type.included[13:length(data.type.included)])
    plot(1, type = "n", xaxt = "n", xlab = "", ylab = "", yaxt = "n", bty = "n")
    legend("left", legend = c("DHS", "MICS", "National survey", 
                              "Other survey", "Service statistics",
                              "Subpopulation", "", 
                              "+:  Higher contraceptive use", # + for age pos and pos bias because of nonpreg etc
                              "-:  Lower contraceptive use",
                              "A:  Other age group", 
                              "F:  Folk methods included", 
                              "S-: Sterilization included", 
                              "S+: Sterilization excluded", "",
                              "Married women", 
                              "Sexually active women", "Ever married/All women",
                              "Both sexes and husband/wives")[data.type.included],
           col = c("red","green","blue",1,ifelse(data.type.included[5], "saddlebrown", "white"),
                   "orange","white", rep("grey",6),"white", rep("grey", 4))[data.type.included],
           bty = "o", cex = ifelse(TIFF, 1.6, 2.1)*cex.adj.factor,
           lwd = 2.8*cex.adj.factor, lty = -1,
           pch = c(rep(21, 4),19,21,21,-1,-1,-1,-1,-1,-1,-1,21,24,25,22)[data.type.included])
  } else {
    data.type.included <- c(data.type.included[1:5], 
                            any(data.type.included[6:length(data.type.included)]))
    plot(1, type = "n", bty = "n", xaxt = "n", xlab = "", ylab = "", yaxt = "n")
    legend("left", legend = c("DHS", "MICS", "National survey", 
                              "Other survey", "Service statistics",
                              "Non-standard observation")[data.type.included],
           col = c("red","green","blue",1,ifelse(data.type.included[5], "saddlebrown", "white"),
                   "grey")[data.type.included],
           bty = "o", cex = ifelse(TIFF, 1.6, 2.1)*cex.adj.factor,
           lwd = 2.8*cex.adj.factor, lty = -1,
           pch = c(rep(21, 4),19,22)[data.type.included])
  }
}
#----------------------------------------------------------------------------------
BreakLongStrings <- function(# Breaks long strings with newline
  ### Breaks long strings with newline, used for plotting
  vector, ##<< Vector of strings
  max.nchar = 20 ##<< Maximum length of string on one line
) {
  vector.output <- sapply(vector, InternalBreakLongString, max.nchar = max.nchar)
  return(vector.output)
}
#----------------------------------------------------------------------------------
InternalBreakLongString <- function(# Breaks a long string with newline
  string, ##<< String
  max.nchar = 20 ##<< Maximum length of string on one line
) {
  if (nchar(string) > max.nchar) {
    niters <- floor(nchar(string)/max.nchar)
    for (iter in niters:1) {
      max.nchar.temp <- max.nchar*iter
      if (nchar(string) > max.nchar.temp) {
        indices.space <- which(strsplit(string, split = "")[[1]] == " ")
        if (length(indices.space) == 0) break()
        indices <- c(indices.space, nchar(string))
        index.replace <- indices[sum(indices <= max.nchar.temp)]
        substr(string, index.replace, index.replace) <- "\n"      
      }
    }
  }
  return(string)
}
#----------------------------------------------------------------------------------
CheckIfAny <- function(values, data.vector) {
  return(any(data.vector %in% values))
}
#----------------------------------------------------------------------------------
# Not used
# InternalPlotModelLegend <- function(# Plot legend with details on model estimates
#   ){
#   par( mar = c(0,0,0,0), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
#   plot(1, type = "n", xaxt = "n", xlab = "", ylab = "", yaxt = "n")
#      legend("center", legend = c("Legend for Estimates:", "Lines: Median and 80% CIs", 
#      "Grey area represents 95% CI", 
#      "Blue line: Excl. AR(1)"),
#      col = c("white", "white", "white", "blue"),
#       #bty = "n", 
#       cex = 2.5,
#       lwd = 3, lty = 1,
#       pch = c(-1,-1,-1))
#       
#   # put back default plot settings
#   par(mar = c(5,5,3,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
# }
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------
# The End!
