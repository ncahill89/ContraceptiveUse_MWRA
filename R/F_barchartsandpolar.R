
GetPolarExtraAggr <- function(# Construct polar histograms for new aggregates
  ### Construct polar histograms for new aggregates
  ### Adapted from:
  ### C. Ladroue, after Kettunen, J. et al. Genome-wide association study identifies multiple loci influencing human serum metabolite levels. Nat Genet advance online publication (2012). URL http://dx.doi.org/10.1038/ng.1073
  extra.aggregate, ##<< Object of class \code{\link{Results}} with \code{CIprop.Lg.Lcat.qt}
  ## for the aggregates.
  add.to.plot.name  = "extra", ## Name to include in name of figures.
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  fig.dir = NULL ##<< Directory to store overview plots. If NULL, folder "fig" in current working directory.
#  plot.tiff  = FALSE##<< TIFF-format? If FALSE, PDF-format is used.
  ){
  
  if (is.null(fig.dir)){
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  resall <- extra.aggregate
  est.years <- dimnames(resall$CIratio.Lg.Lcat.qt[[1]][[1]])[[2]]
  years <- c(1990, 2000, 2010)
  item <- names(resall$CIprop.Lg.Lcat.qt)
  for (year in years){
    unmet.c <- unlist(lapply(resall$CIprop.Lg.Lcat.qt,function(l) 
      l$Unmet[3, est.years == paste(year+0.5)]))
    total.c <- unlist(lapply(resall$CIprop.Lg.Lcat.qt,function(l) 
      l$Total[3, est.years == paste(year+0.5)]))
    family <- item
    df.all <-data.frame(family=family,  
                        item= item,
                        Total  = total.c, Unmet = unmet.c, None = 1-total.c-unmet.c)
    df<-melt(df.all,c("family","item"),variable_name="score") # from wide to long
    p<-polarHistogram(df,familyLabels=FALSE, direction="outwards",
                      guides=c(20,40,60, 80))
    ##details<< Polar histograms for 1990, 2000 and 2010 are saved in \code{fig.dir}.
#    if (plot.tiff){
      
#    } else {
    pdf(file.path(fig.dir, paste0(run.name, "polarhist", add.to.plot.name, year, ".pdf")), width = 10, height = 8) # change JR, 20140418
    print(p)
    dev.off()
  }
  ##value<< NULL
}

#----------------------------------------------------------------------
GetAllPolar <- function(# Construct polar histograms
  ### Construct polar histograms. 
  ### Adapted from:
  ### C. Ladroue, after Kettunen, J. et al. Genome-wide association study identifies multiple loci influencing human serum metabolite levels. Nat Genet advance online publication (2012). URL http://dx.doi.org/10.1038/ng.1073
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  fig.dir = NULL ##<< Directory to store overview plots. If NULL, folder "fig" in current working directory.
  # plot.tiff  = FALSE##<< TIFF-format? If FALSE, PDF-format is used.
){
  if (is.null(fig.dir)){
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  
  load(file = file.path(output.dir,"mcmc.meta.rda")) # change JR, 20140418
  do.country.specific.run <- mcmc.meta$general$do.country.specific.run # change JR, 20131105
  load(file = file.path(output.dir,"res.country.rda")) # change JR, 20140418
  if (!do.country.specific.run)
    load(file = file.path(output.dir,"res.aggregate.rda")) # change JR, 20140418
  country.info <- mcmc.meta$data.raw$country.info
  region.info <- mcmc.meta$data.raw$region.info
  
  est.years <- dimnames(res.country$CIratio.Lg.Lcat.qt[[1]][[1]])[[2]]
  years <- c(1990, 2000, 2010)
  if (!do.country.specific.run) { # change JR, 20131105
    for (res in 1:2) {
      resall <- list(res.aggregate, res.country)[[res]]
      item <- names(resall$CIprop.Lg.Lcat.qt)
      for (year in years){
        unmet.c <- unlist(lapply(resall$CIprop.Lg.Lcat.qt,function(l) 
          l$Unmet[3, est.years == paste(year+0.5)]))
        total.c <- unlist(lapply(resall$CIprop.Lg.Lcat.qt,function(l) 
          l$Total[3, est.years == paste(year+0.5)]))
        if (res==2){
          family <- country.info$namesubreg.c
          item <- InternalMakeCountryNamesShort(country.info$name.c)
        } else {  
          family <- item
        }  
        df.all <-data.frame(family = family,  
                            item = item,
                            Total = total.c, Unmet = unmet.c, None = 1-total.c-unmet.c)
        if (res==2){
          for (reg in region.info$name.reg){
            df <- df.all[country.info$namereg.c==reg,]
            df<-melt(df,c("family","item"),variable_name="score") # from wide to long
            p<-polarHistogram(df,familyLabels=FALSE, direction="outwards",
                              guides=c(20,40,60, 80))
            #,circleProportion = 0.9)
            pdf(file.path(fig.dir, paste0(run.name, "polarhist", reg, year, ".pdf")), width = 10, height = 8) # change JR, 20140418
            print(p)
            dev.off()
          }
        } else {  
          df<-melt(df.all,c("family","item"),variable_name="score") # from wide to long
          p<-polarHistogram(df,familyLabels=FALSE, direction="outwards",
                            guides=c(20,40,60, 80))
          pdf(file.path(fig.dir, paste0(run.name, "polarhist", "ALL", year, ".pdf")), width = 10, height = 8) # change JR, 20140418
          print(p)
          dev.off()
        }
      }
    }
  } else { # change JR, 20131105
    resall <- res.country
    item <- names(resall$CIprop.Lg.Lcat.qt)
    for (year in years){
      unmet.c <- unlist(lapply(resall$CIprop.Lg.Lcat.qt,function(l) 
        l$Unmet[3, est.years == paste(year+0.5)]))
      total.c <- unlist(lapply(resall$CIprop.Lg.Lcat.qt,function(l) 
        l$Total[3, est.years == paste(year+0.5)]))
      family <- country.info$namesubreg.c
      item <- InternalMakeCountryNamesShort(country.info$name.c)
      df.all <-data.frame(family = family,  
                          item = item,
                          Total = total.c, Unmet = unmet.c, None = 1-total.c-unmet.c)
      for (reg in region.info$name.reg){
        df <- df.all[country.info$namereg.c==reg,]
        df<-melt(df,c("family","item"),variable_name="score") # from wide to long
        p<-polarHistogram(df,familyLabels=FALSE, direction="outwards",
                          guides=c(20,40,60, 80))
        #,circleProportion = 0.9)
        pdf(file.path(fig.dir, paste0(run.name, "polarhist", reg, year, ".pdf")), width = 10, height = 8)
        print(p)
        dev.off()
      }
    }
  }
  cat(paste0("Polar histograms constructed and saved in ", fig.dir, "\n"))
  ##details<< Polar histograms for 1990, 2000 and 2010 are saved in \code{fig.dir} for all subregions, and for all countries grouped by region.
  ##value<< NULL 
}
#-------------------------------------------------------------------
GetAllBarCharts <- function(# Bar charts for props and counts for 1990, 2000 and 2010
  ### Bar charts for props and counts for 1990, 2000 and 2010
  run.name = "test", ##<< Run name
  output.dir = NULL, ##<< Directory where MCMC array and meta are stored.
  ## If NULL, it's \code{output/run.name}, default from \code{runMCMC}.
  fig.dir = NULL, ##<< Directory to store overview plots. If NULL, folder "fig" in current working directory.
  extra.aggregate = NULL, ##<< If NULL, UNPD aggregates are used. 
  ## Alternatively, an object of class \code{\link{Results}} with \code{CIprop.Lg.Lcat.qt}
  ## for plotting alternative aggregates.
  add.to.plot.name  = "extra", ## Name to include in name of figures, used only if alternative aggregates were provided.
  plot.tiff  = FALSE##<< TIFF-format? If FALSE, PDF-format is used.
){
  if (is.null(fig.dir)){
    fig.dir <- file.path(getwd(), "fig/")
    dir.create(fig.dir, showWarnings = FALSE)
  }
  if (is.null(output.dir)){
    output.dir <- file.path(getwd(), "output", run.name, "/")
  }
  years <- c(1990, 2000, 2010)
  if (is.null(extra.aggregate)){
    load(file = file.path(output.dir,"mcmc.meta.rda")) # change JR, 20140418
    #country.info <- mcmc.meta$data.raw$country.info
    region.info <- mcmc.meta$data.raw$region.info
    # find order of subregions in regions...
    order.subreg.g <- rev(order(mcmc.meta$winbugs.data$reg.subreg))
    load(file = file.path(output.dir,"res.aggregate.rda")) # change JR, 20140418
    select.g <- is.element(names(res.aggregate$CIprop.Lg.Lcat.qt), region.info$name.subreg)
  } else {
    res.aggregate <- extra.aggregate
    select.g <- rep(TRUE, length(res.aggregate$CIprop.Lg.Lcat.qt))
    order.subreg.g <- seq(1, length(select.g))
  }  
  est.years <- dimnames(res.aggregate$CIratio.Lg.Lcat.qt[[1]][[1]])[[2]]
  tempinternal <- function(){
    prop.gq3t <- NULL
    for (year in years){
      prop.gq3 <- From.Lg.Lcat.qtTo.gq(res.Lg.Lcat.qt = res.Lg.Lcat.qt,
                                       cat = cat, year = paste(year+0.5))[, c(1,3,5)]
      if (!is.null(prop.gq3t)){ 
        prop.gq3t <- abind(prop.gq3t,prop.gq3, along = 3 )
      } else {
        prop.gq3t <- prop.gq3
      }
      fig.name <- file.path(fig.dir, paste0(run.name, "BarChart_", nameadd, "_", floor(year))) # change JR, 20140418
      if (plot.tiff){
        tiff(filename=paste0(fig.name,".tif"),
             width = 8, height = 12,
             units="cm", pointsize=6,res=300, bg="white",compression="lzw")
      } else {
        pdf(paste0(fig.name,".pdf"), width = 8, height = 12)
      }
      BarChart(res.gq3t = prop.gq3[select.g,],
               main = paste(cat, " (", floor(year), ")", sep = ""), 
               order.g = order.subreg.g
               )
      dev.off()
    } # end years
    fig.name <- file.path(fig.dir, paste0(run.name, "BarChart_", nameadd,"_allyears")) # change JR, 20140418
    if (plot.tiff){
      tiff(filename=paste0(fig.name,".tif"),
           width = 8, height = 12,
           units="cm", pointsize=6,res=300, bg="white",compression="lzw")
    } else {
      pdf(paste0(fig.name,".pdf"), width = 8, height = 12)
    }
       BarChart(res.gq3t = prop.gq3t[select.g,,], main = paste(cat),
               yearsnames  = years, order.g = order.subreg.g, plot.tiff = plot.tiff)
    dev.off()
  } # end function
  #percentiles <- dimnames(res.country$CIratio.Lg.Lcat.qt[[1]][[1]])[[1]] # 0.5
  # assume it's as usual...
  
  list.Lg.Lcat.qt <- list(res.aggregate$CIprop.Lg.Lcat.qt, res.aggregate$CIcount.Lg.Lcat.qt)
  for (cat in c("Total", "Unmet")){
    for (res in 1:2){
      res.Lg.Lcat.qt <- list.Lg.Lcat.qt[[res]]
      nameadd <- paste0(add.to.plot.name, "_", cat, c("prop", "count")[res])
      tempinternal()
    } # end res   
  } # end cat
  cat <- "Met Demand"
  res.Lg.Lcat.qt <- res.aggregate$CIratio.Lg.Lcat.qt
  nameadd <- paste0(add.to.plot.name, cat)
  tempinternal()
  ##details<< Bar charts for 1990, 2000 and 2010 (individually, as well as for all years combined) 
  ## are saved in \code{fig.dir} for either all subregions (in the default function call)
  ## or for alternative groupings (if provided by user).
  ##value<< NULL 
  return(invisible())  
}
#--------------------------------------------------------------------
GetNiceSequence <- function(# Get a nice sequence for plotting
  ### Get a nice sequence for plotting
  xmin =0, 
  xmax){
  # assume that xmin is a nice number already...
  # get outcome between 0 and 1 times multiplier...
  # should be an easier way!!
  l <-  length(strsplit(as.character(round(xmax)), split = NULL)[[1]])  
  multiplier <- 10^ifelse(l==1,0, l)
  xmin.new <- xmin/multiplier
  xmax.new <- xmax/multiplier
  if ((xmax.new - xmin.new) <= 0.1){
    seq.nice <- multiplier*seq(xmin.new, 0.02+xmax.new, by = 0.02)
    return(seq.nice)
  }
  if ((xmax.new - xmin.new) <= 0.4){
    seq.nice <-multiplier*seq(xmin.new, 0.05+xmax.new, by = 0.05)
    return(seq.nice)
  }
  if ((xmax.new - xmin.new) <= 1){
    seq.nice <- multiplier*seq(xmin.new, min(1, 0.1+xmax.new), by = 0.1)
    return(seq.nice)
  }
}

#GetNiceSequence(xmin =0, xmax=1)
#GetNiceSequence(xmin =0, xmax=56)
#GetNiceSequence(xmin =0, xmax=981)
#--------------------------------------------------------------------
BarChart <- function(# Plot bar chart
  ### for various years/outcomes/...
  res.gq3t, ##<< rownames are the names to plot
  ## q is CI and median, t is optional for various years
  #xlab = NULL,
  main = NULL, 
  yearsnames = NULL, 
  order.g  = NULL, ##<< If NULL, plotted in order of res.g
  prop.break.down.Lg.ct  = NULL, ##<< Optional, 
  ## for break down of counts into prop by country
  plot.tiff = FALSE
  ){
  
  G <- dim(res.gq3t)[1]
  if (is.null(order.g)) order.g <- seq(G, 1,-1)
  
  if (length(dim(res.gq3t))==3){
    nyears <- dim(res.gq3t)[3]
    if (nyears >4){
      print("Stop, plot max 4 years!")
      return()
    }
    final.gq3t <- res.gq3t[order.g,,]
  } else {
    nyears <- 1
    name.all.g <- dimnames(res.gq3t)[[1]]
    # take care: 3rd dimension seems to be lost when remnaming, reordering etc
    final.gq3t <- res.gq3t[order.g,]
    dim(final.gq3t) <- c(dim(final.gq3t),1)
    dimnames(final.gq3t)[[1]] <- name.all.g[order.g]
  }
  name.g <- dimnames(final.gq3t)[[1]]
  name.g <- ifelse(name.g=="Australia and New Zealand", "Australia and NZ", paste(name.g))
  name.g <- ifelse(name.g=="Latin America and the Caribbean", "LAC", paste(name.g))
  if (max(final.gq3t)<1){
    xlab = "Percentage"
    final.gq3t = 100*final.gq3t
    xmax <- min(100, max(final.gq3t)*1.05)
  } else {
    xmax <- max(final.gq3t)*1.05
    xlab = "Number of MWRA (*1,000)"
  }
  # grey vertical bounds in background... at...
  vertbounds <- GetNiceSequence(xmin = 0, xmax = xmax)
  
  par(mfrow = c(1,1), mar = c(5,10,3,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  plot(seq(1, G)~rep(1,G), type = "n", 
       ylab = "",  yaxt = "n",
       xlim = c(0,xmax), xaxt = "n",
       bty = "n", 
       ylim = c(1,G),
       xlab = xlab, main = main)
  axis(2, at = seq(1, G), labels = name.g, las = 1, cex.axis = 1)
  for (s in seq(1, length(vertbounds)-1,2)){
    polygon(vertbounds[c(s,s,s+1,s+1,s)],
            c(0,G+1,G+1,0,0), border = "NA", 
            col = adjustcolor("lightgrey", alpha.f = 0.4))
  }
  axis(1, at = vertbounds, labels = vertbounds, las = 1, cex.axis = 1)
  palette("default")
  mycols <- adjustcolor(palette(), alpha.f = 0.4)
  
  # depending on number of years, halfwidths and add need to adjusted
  # note that ylab is increasing as usual,
  # so if you want years in decreasing order, add needs to decrease
  if (nyears==1){
    halfwidth <- 0.3 
    add.t <- c(0, NA)
  } 
  if(nyears==2){
    halfwidth <- 0.15
    add.t <- c(0.15, -0.15)
  }
  if(nyears==3){
    halfwidth <- 0.1
    add.t <- c(0.2,0, -0.2)
  }
  if(nyears==4){
    halfwidth <- 0.075
    add.t <- c(0.3, 0.15, 0, -0.15, -0.3)
  }
  
  #res.Lt.c <- list(res2010.c = res2010.c, res2015.c = res2015.c)
  for (t in 1:nyears){
    res.gq <- final.gq3t[,,t]
    #res.c <- res.Lt.c[[t]]
    add <- add.t[t]
# TO CLEAN UP!
#     for (j in seq(1, G)){
#       s <- order[j] # get right region for plot index j
#       # props of country count in region by medians
#       props <- res.c[is.element(paste(country.info$subreg.c), name.s[s])]/
#         sum(res.c[is.element(paste(country.info$subreg.c), name.s[s])])
#       # use country order from 2010
#       cprops <- c(0, cumsum(props[orderp[[s]]]))
#       for (i in 2:length(cprops)){
#         polygon(res.sq[s,3]*c(cprops[i-1],cprops[i],cprops[i],cprops[i-1],cprops[i-1]), 
#                 add + j + c(-width,-width,width, halfwidth, -width), 
#                 border = "NA", col = mycols[i])
#       }
#       # put box around
#       polygon(c(0, res.sq[s,3],res.sq[s,3], 0,0), 
#               add + j + c(-width,-width,+width, +width,-width), border = 1)
#     } # end s
    if (!is.null(prop.break.down.Lg.ct)){
    
    } else {
      for (g in 1:G){
        polygon(c(0, res.gq[g,2],res.gq[g,2], 0,0), 
              add + g + c(-halfwidth,-halfwidth,+halfwidth, +halfwidth,-halfwidth), border = 1,
                col = mycols[t] )      
      }
    }
    for (g in 1:G){
      points(add+seq(1, G)~res.gq[,2], pch = 20, lwd = ifelse(plot.tiff,1,2))
      segments(res.gq[,1], add+seq(1, G), res.gq[,3], 
             add+seq(1, G), lwd = ifelse(plot.tiff,1,2), col = 1)
    } # end g
  } # end t
  if (!is.null(yearsnames)){
    legend("bottomright", legend = yearsnames, col = mycols[seq(1, nyears)], pch = -1, lty = 1, 
           cex = 0.8, lwd = 5)
  }
  return(invisible())
}
#----------------------------------------------------------------------
# The End!
