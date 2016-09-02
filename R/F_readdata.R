#----------------------------------------------------------------------------------
# Leontine Alkema
# F_readdata.R
#----------------------------------------------------------------------------------
# Contains:
# general math functions
# function to read in data, and make summary tables (plots)
# functions to read/summarize regional info
# function to change country names

#----------------------------
logit <- function(# Logit transform
  ### Logit transform
  x##<< between 0 and 1
  ){
  log(x/(1-x))
}

#----------------------------
invlogit <- function(# Inverse-Logit transform
  ### Inverse-Logit transform
  y){
 x <- 1/(1+exp(-y)) 
 return(x)
}

#----------------------------
LogitMinMax <- function(# Logit-minmax transform
  ### y = log((x-xmin)/(xmax-x))
  x, ##<< between xmin and xmax
  xmin, xmax){
   y <- log((x-xmin)/(xmax-x))
   return(y)
}

#----------------------------------------------------------------------------------
InvLogitMinMax <- function(#  Inverse Logit-minmax transform
  ### x = (ymax*exp(y) + ymin)/(1+exp(y)), ends up between ymin and ymax
  y, ymin=0, ymax = 1){
  x <- (ymax*exp(y) + ymin)/(1+exp(y))
  return(x)
}
#----------------------------------------------------------------------------------
PreprocessData <- function(# Pre-process contraceptive prevalence data
  ### Pre-process contraceptive prevalence data
  data.csv = NULL, ##<< If \code{NULL}, data set included in package is used. 
  ## To use alternative data set, use \code{data.csv = .../dataCPmodel.csv}, where ... refers to the file path where file is located.
  iso.select = NULL ##<< If \code{NULL}, data of all countries/subpopulations in data set are read in, else
  ## only data of countries/subpopulations with selected ISO code is read in. # change JR, 20131104
) {
  data.raw <- read.csv(file = data.csv, header = TRUE, as.is = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
  names(data.raw)[grepl("Country.letter.code|Country..letter.code", names(data.raw))] <- "Country..letter.code"
  names(data.raw)[names(data.raw) == "Age.range"] <- "Age..range"
  names(data.raw)[grepl("exclude", names(data.raw), ignore.case = TRUE)] <- "EXCLUDE1isyes"
  names(data.raw)[grepl("GEO.biases", names(data.raw), ignore.case = TRUE)] <- "GEO.biases..unknown.direction."
  if (is.null(data.raw$Note.on.methods)) 
    data.raw$Note.on.methods <- rep(NA, nrow(data.raw))
  # Select data only for countries in iso.select and order data by country order given in iso.select
  if (!is.null(iso.select)) {
    if (all(grepl("[[:digit:]]", iso.select))) {
      data.raw <- data.raw[as.numeric(data.raw$ISO.code) %in% iso.select, ]
      data.raw <- data.raw[order(factor(data.raw$ISO.code, levels = iso.select)), ]
    } else if (all(grepl("[[:alpha:]]", iso.select))) {
      data.raw <- data.raw[gsub(" ", "", as.character(data.raw$Country..letter.code)) %in% iso.select, ]
      data.raw <- data.raw[order(factor(gsub(" ", "", as.character(data.raw$Country..letter.code)), 
                                        levels = iso.select)), ]
    }
    cat(paste0("Note: Only data for ", paste(iso.select, collapse = ", "), " is read in.\n"))
  }
  
  ##details<<
  ## Observations are excluded if column \code{EXCLUDE1isyes} == 1.
  ## Observations are excluded if \code{Note.on.methods} is "Data pertain to methods used since the last pregnancy."
  ## or "Data pertain to past or current use.".
  ## Observations are excluded if total use is missing, unless they are from service statistics.
  # (one survey in Bhutan with modern only)
  if (is.null(data.raw$Note.on.methods))
    data.raw$Note.on.methods <- rep(NA, nrow(data.raw))
  remove <- (!is.na(data.raw$EXCLUDE1isyes) & data.raw$EXCLUDE1isyes == 1) |
    # is.na(prop.raw) | 
    (is.na(data.raw$Contraceptive.use.ANY) & !grepl("Service statistic", data.raw$Data.series.type)) |
    (!is.na(data.raw$Note.on.methods) & (data.raw$Note.on.methods == "Data pertain to methods used since the last pregnancy." |
                                           data.raw$Note.on.methods == "Data pertain to past or current use."))
  if (any(remove)) {
    data.preprocessed <- data.raw[!remove, ]
    cat(paste0(sum(remove), " observations removed.\n"))
  } else {
    data.preprocessed <- data.raw
  }
  return(data.preprocessed = data.preprocessed)
}
#----------------------------------------------------------------------------------
ReadDataAll <- function(# Read contraceptive prevalence data
  ### Read contraceptive prevalence data
  data.csv = NULL, ##<< If \code{NULL}, data set included in package is used. 
  ## To use alternative data set, use \code{data.csv = .../dataCPmodel.csv}, where ... refers to the file path where file is located.
  regioninfo.csv = NULL, ##<< If \code{NULL}, region info included in package is used. 
  ## To use alternative csv file, use \code{regioninfo.csv = .../Country-and-area-classification.csv}, where ... refers to the file path where file is located.
  iso.select = NULL, ##<< If \code{NULL}, data of all countries/subpopulations in data set are read in, else
  ## only data of countries/subpopulations with selected ISO code is read in.
  iso.country.select = NULL, ##<< If not \code{NULL}, data is treated as data from subpopulation of country with the ISO code iso.country.select.
  name.country.select = NULL, ##<< Country name corresponding to \code{iso.country.select}. Cannot be \code{NULL} if \code{iso.country.select} is not \code{NULL}.
  do.SS.run.first.pass = do.SS.run.first.pass, ##<< Do first pass run of run with SS data?
  countrycodes.csv = "data/Country-names-and-codes.csv", ##<< csv file with ISO 3-character and 3-digit country codes, only used to convert iso.country.select
  ## into 3-digit country code if given in 3-character country code.
  html.file = NULL ##<<If not NULL, summary results are written to this HTML file.
){
  if (!is.null(iso.country.select) & is.null(name.country.select)) # change JR, 20140409
    stop("name.country.select cannot be NULL if iso.country.select is non-NULL.")
  # Create new html file
  if (!is.null(html.file)){
    cat("", file = html.file, append = F)
    print(paste("Summary stats written to", html.file))
  }
  if (is.null(data.csv)){
    data.csv <- file.path(find.package("ContraceptiveUse"), "data", "dataCPmodel.csv")
  }  
  cat(paste("CP data read from", data.csv), "\n")
  if (is.null(regioninfo.csv)){
    regioninfo.csv <- file.path(find.package("ContraceptiveUse"), "data", "Country-and-area-classification.csv")
    
  }
  cat(paste("Country/region info read from", regioninfo.csv), "\n")
  #----------------------------------------------------------------------
  data.unsorted <- PreprocessData(data.csv = data.csv,
                                  iso.select = iso.select)
  if (do.SS.run.first.pass) {
    remove <- grepl("Service statistic", data.unsorted$Data.series.type)
    if (any(remove)) {
      data.unsorted <- data.unsorted[!remove, ]
      cat(paste0(sum(remove), " service statistics data observations removed.\n"))
    }
  }
  #----------------------------------------------------------------------
  J <- nrow(data.unsorted)
  # For some obs where there is no break-down of total into modern and trad CP
  # Re-order the rows such that all modern obs are first, followed by obs with total only
  # (makes things easier when constructing the data set for BUGS)
  # sum(is.na(data.unsorted$Contraceptive.use.MODERN))
  order <- c(seq(1, J)[!is.na(data.unsorted$Contraceptive.use.MODERN)], seq(1, J)[is.na(data.unsorted$Contraceptive.use.MODERN)])
  # data.unsorted$Contraceptive.use.MODERN[order]
  years.j <- ((data.unsorted$Start.year+data.unsorted$End.year)/2)[order]
  props.tot.j <- data.unsorted$Contraceptive.use.ANY[order]/100
  props.modern.j <- data.unsorted$Contraceptive.use.MODERN[order]/100
  props.trad.j <- props.tot.j - props.modern.j
  data <- data.unsorted[order,]
  if (is.null(iso.country.select)) { # change JR, 20140404
    name.j <- data$Country
    name.unsorted <- data.unsorted$Country
  } else {
    name.j <- data$New.population
    name.unsorted <- data.unsorted$New.population
  }
  ##details<< Age categorization:
  # 1. Fix age groups that excel made into dates
  # change JR, 20140409
  data$Age..range <- gsub(" ", "", data$Age..range)
  Age <- InternalFixRange(data$Age..range)
  # Age <- ifelse(data$Age..range == "Dec-49", "12-49", 
  #                 ifelse(data$Age..range == "Oct-49", "10-49", 
  #                 paste(data$Age..range)) )
  ##details<< - If group starts at 13-17 and ends at 47-51: base line (0).
  ##details<< - If group start at 13-17 but ends after 51: negative bias ("neg").
  ##details<< - "See notes" is positive.
  ##details<< - Other groups get flagged with "?" (other).
  age.cat.j <- rep(NA,J)
  options(warn=-1) # next part will give warnings, suppress those
  for (j in 1:J){
    split <- as.integer(strsplit(Age[j], split = "")[[1]]) # gives a warning
    age.cat.j[j] <- ifelse(length(split)!=5 ,"?",
            ifelse(split[1]==1 & split[2]<8 & split[2]>2 & ((split[4]==4 & split[5]>6)|(split[4]==5 & split[5]<2)) , 0,
            ifelse(split[1]==1 & split[2]<8 & split[2]>2 & split[4]==5 & split[5]>1 , "-","?")))
  }
  options(warn=0) # back to default setting
  age.cat.j <- ifelse(Age == "See notes", "+", age.cat.j)
  if (!is.null(html.file)){
    print(xtable(xtabs(~Age+age.cat.j), digits = c(0,0,0,0), type = "html", 
          caption = "Age"), type="html", file = html.file, append = T)
    print(xtable(table(age.cat.j), digits = c(0,0), type = "html", 
          caption = "Age"), type="html", file = html.file, append = T)
    # print(xtabs(~Age+age.cat.j))
  }
  # get other columns in shape so that they're meaningful
  # then have separate function to get the input for bugs
  source <- data$Data.series.type
  source.j <- ifelse(source == "DHS", "DHS", 
      ifelse(is.element(source, c("MICS", "MICS1", "MICS2", "MICS3", "MICS4")), "MICS", # change JR, 20140310: added MICS and MICS4
      ifelse(source == "National survey", "NS",
      ifelse(grepl("Service statistic", source), "SS", # change JR, 20131120 # change JR, 20140418
             "Other"))))
  if (!is.null(html.file)){
    print(xtable(xtabs( ~ source + source.j), digits = c(0,0,0,0,0), type = "html", 
          caption = "Source"), type="html", file = html.file, append = T)
    print(xtable(table(source.j), digits = c(0,0), type = "html", 
                 caption = "Source"), type="html", file = html.file, append = T)
  }
  geo.j <- data$GEO.biases..unknown.direction.
  if (all(is.na(geo.j))) # change JR, 20131121: convert NA's to ""
    geo.j[is.na(geo.j)] <- ""  
  geo.short.j <- ifelse(geo.j != "", 1, 0)
  if (!is.null(html.file)){
    print(xtable(table(geo.j), digits = c(0,0), type = "html", 
                 caption = "Geo"), type="html", file = html.file, append = T)
    print(xtable(table(geo.j!=""), digits = c(0,0), type = "html", 
                 caption = "Geo"), type="html", file = html.file, append = T)
  }
  
  poptype.j <- data$Population.type
  poptype.short.j <- ifelse(is.element(poptype.j, c("BS", "HW")), "BSHW", paste(poptype.j))
  if (!is.null(html.file)){
    print(xtable(xtabs(~ poptype.j + poptype.short.j), digits = rep(0,6), type = "html", 
          caption = "Pop type"), type="html", file = html.file, append = T)
    print(xtable(table(poptype.short.j), digits = c(0,0), type = "html", 
                 caption = "Pop type"), type="html", file = html.file, append = T)
  }
  
  folkbias.j <- data$Folk.method.positive.bias
  if (all(is.na(folkbias.j))) # change JR, 20131121: convert NA's to ""
    folkbias.j[is.na(folkbias.j)] <- ""  
  #  negbias.j <- data$Negative.bias
  #   if (all(is.na(negbias.j))) # change JR, 20131121: convert NA's to ""
  #     negbias.j[is.na(negbias.j)] <- ""  
  posbias.j <- data$Non.pregnant.and.other.positive.biases
  if (all(is.na(posbias.j))) # change JR, 20131121: convert NA's to ""
    posbias.j[is.na(posbias.j)] <- ""  
  mod.bias.j <- data$Modern.method.bias
  if (all(is.na(mod.bias.j))) # change JR, 20131121: convert NA's to ""
    mod.bias.j[is.na(mod.bias.j)] <- ""  
  
  if (!is.null(html.file)){
    print(xtable(table(folkbias.j), digits = rep(0,2), type = "html", 
                 caption = "Folk"), type="html", file = html.file, append = T)
    print(xtable(table(folkbias.j!=""), digits = rep(0,2), type = "html", 
          caption = "Folk"), type="html", file = html.file, append = T)
    print(xtable(table(posbias.j), digits = rep(0,2), type = "html", 
          caption = "+bias"), type="html", file = html.file, append = T)
    print(xtable(table(posbias.j!=""), digits = rep(0,2), type = "html", 
                 caption = "+bias"), type="html", file = html.file, append = T)
    print(xtable(table(mod.bias.j), digits = rep(0,2), type = "html", 
          caption = "Bias modern"), type="html", file = html.file, append = T)
  }
  source.unmet.j <- ifelse(source.j=="DHS", "DHS", "Other")
  props.unmet.j <- data$Unmet/100
  if (!is.null(html.file)){
    print(xtable(table(source.unmet.j[!is.na(props.unmet.j)]), digits = rep(0,2), type = "html", 
                 caption = "Source Unmet need"), type="html", file = html.file, append = T)
  }
  ##details<<
  ## \code{data} is a data frame with
  ## (iso.j, name.j,
  ##                  years.j, props.tot.j, props.modern.j, props.trad.j, 
  ##                  age.cat.j, source.j, poptype.j, geo.j, 
  ##                  posbias.j, folkbias.j, mod.bias.j,
  ##                  props.unmet.j, source.unmet.j).
  data <- data.frame(iso.j = data$ISO.code, name.j = name.j, # change JR, 20140404
                     years.j, props.tot.j, props.modern.j, props.trad.j, 
                     age.cat.j, source.j, poptype.j, geo.j, 
                     posbias.j, folkbias.j, mod.bias.j,
                     props.unmet.j = props.unmet.j, source.unmet.j = source.unmet.j,
                     stringsAsFactors=FALSE) 
  country.info <- unique(data.frame(name.unsorted, # change JR, 20140404 
                                    gsub(" ", "", data.unsorted$ISO.code), # change JR, 20131104 
                                    gsub(" ", "", data.unsorted$Country..letter.code), # change JR, 20140404
                                    stringsAsFactors=FALSE))
  names(country.info) <- c("name.c", "iso.c", "code.c")
  # name.short.c <- InternalMakeCountryNamesShort(country.info$name.c)
  # note: when using table, the output is in alphabetical order, thus China and China, Hong Kong are swapped
  # swap back!
  # sum(names(table(data$Country)[name.c])!=  paste(name.c))
  # sum((table(data$Country)[name.c])!=  table(data$Country))
  N.c <- table(name.unsorted)[country.info$name.c]
  if (is.null(iso.country.select)) { # change JR, 20140404
    reg.country.info.c <- InternalGetRegionInfoForCountry(iso.c = country.info$iso.c, 
                                                          regioninfo.csv = regioninfo.csv,
                                                          countrycodes.csv = countrycodes.csv) # change JR, 20140407
  } else {
    reg.country.info.c <- InternalGetRegionInfoForCountry(iso.c = iso.country.select, 
                                                          regioninfo.csv = regioninfo.csv,
                                                          countrycodes.csv = countrycodes.csv) # change JR, 20140407
  }
  # print(reg.country.info.c)
  country.info <- data.frame(country.info, N.c, 
      # note: not storage efficient to store a name and number for (sub)-regions,
      # but having some issues with factors in earlier versions, so this is safest option :)
      reg.country.info.c[, c("namereg.c", "reg.c" ,  
                             "namesubreg.c", "subreg.c", "ssa.c", "dev.c")],
                             stringsAsFactors=FALSE)
  if (!is.null(iso.country.select)) # change JR, 20140404
    country.info <- data.frame(country.info, 
                               iso.country.select = iso.country.select, 
                               name.country.select = name.country.select, # change JR, 20140409
                               stringsAsFactors=FALSE)
  if (!is.null(html.file)){
    print(xtable(country.info, type = "html", 
          caption = "Country info"), type="html", file = html.file, append = T)
  }
  region.info <- InternalGetRegionInfoGeneral(country.info = country.info)

  ##value<<               
  return(list(data = data,##<< Object of class \code{data}, data frame, described here in Details
              country.info = country.info,##<< Object of class \code{country.info}
              region.info = region.info##<< Object of class \code{region.info}
  ))
}    
#----------------------------------------------------------------------------------
InternalFixRange <- function(# Fix range character vector that got converted to dates in Excel
  ### Replaces month names with numbers
  x
) {
  x <- as.character(x)
  from <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  to <- as.character(1:12)
  sub.table <- data.frame(from = from, to = to, stringsAsFactors = FALSE)
  for (i in 1:nrow(sub.table))
    x <- gsub(sub.table$from[i], sub.table$to[i], x)
  return(x)
}
#----------------------------------------------------------------------------------
InternalGetRegionInfoGeneral <- function(# Summarize regional info
  ### Summarize regional info to construct region.info
  country.info
){
  n.subreg <- length(unique(country.info$subreg.c))
  n.reg <- length(unique(country.info$reg.c))
  name.subreg <- rep("", n.subreg)
  name.reg <- rep("", n.reg)
  # make sure that name.subreg[1] corresponds to same 1 in country.info
  #name.subreg <- levels(as.factor(country.info$subreg.c))
  #name.reg <- levels(as.factor(country.info$reg.c))
  for (subreg in 1:n.subreg){
    name.subreg[subreg] <- country.info$namesubreg.c[country.info$subreg.c==subreg][1]
  }
  for (reg in 1:n.reg){
    name.reg[reg] <- country.info$namereg.c[country.info$reg.c==reg][1]
  }
  name.reg.short <- name.reg
  name.reg.short <- ifelse(name.reg=="Latin America and the Caribbean","LAC", paste(name.reg.short))
  name.reg.short <- ifelse(name.reg=="Northern America","NA", paste(name.reg.short))
  
  region.info <- list(name.subreg = as.character(name.subreg), 
                      name.reg = as.character(name.reg), 
                      name.reg.short = as.character(name.reg.short), 
                      n.subreg = n.subreg, n.reg = n.reg)
  ##value<< 
  return(region.info##<< Object of class \code{\link{region.info}}
         ## where name.subreg[1] refers to index 1 in country.info$subreg.c etc
         )
}

#---------------------------------------------------------------------------------------
InternalGetRegionInfoForCountry <- function(# Find (sub)region info for country vector
  ### Find (sub)region info for country vector
  iso.c,
  regioninfo.csv, ##<< csv with region info
  countrycodes.csv ##<< csv with country codes info
){
  regions <- read.csv(file = regioninfo.csv, header = TRUE, stringsAsFactors = FALSE)
  
  # change JR, 20140404
  if (!all(grepl("[[:digit:]]", iso.c))) {
    if (all(grepl("[[:alpha:]]", iso.c))) {
      iso.c <- InternalGetCountryCodes(iso.c = iso.c, countrycodes.csv = countrycodes.csv)
    } else {
      stop("Elements of iso.c must be all alphabet letters or all digits.")
    }
  }
  
  C <- length(iso.c)
  dev.c <-  ldev.c <- namesubreg.c <- namereg.c <- ssa.c <- rep("", C)
  for (c in 1:C){
    namesubreg.c[c] <- paste(regions$Region[which.max(regions$ISO.Code==iso.c[c])])
    namereg.c[c] <- paste(regions$Major.area[which.max(regions$ISO.Code==iso.c[c])])
    ssa.c[c] <- paste(regions$Sub.Saharan..Africa[which.max(regions$ISO.Code==iso.c[c])])
    dev.c[c] <- ifelse(paste(regions$Least.developed.country[which.max(regions$ISO.Code==iso.c[c])])=="Yes","Poor",
                ifelse(paste(regions$Developed..region[which.max(regions$ISO.Code==iso.c[c])])=="Yes", "Rich", "Med"))
  }
  subreg.c <- as.numeric(as.factor(namesubreg.c))
  reg.c <- as.numeric(as.factor(namereg.c))
  reg.country.info.c <- data.frame(as.character(iso.c), 
                                   reg.c, namereg.c, 
                                   subreg.c, namesubreg.c, 
                                   ssa.c, dev.c,
                                   stringsAsFactors=FALSE)
  names(reg.country.info.c) <- c("iso.c", "reg.c", "namereg.c",
        "subreg.c", "namesubreg.c", "ssa.c", "dev.c")
  
  ##value<< Data frame with (iso.c, reg.c, namereg.c, subreg.c, namesubreg.c, ssa.c, dev.c)
  ## where ssa.c and dev.c are factors, 
  ## and reg.c and subreg.c are just integers
  ## part of Object \code{\link{country.info}}.
  return(reg.country.info.c)
}
#----------------------------------------------------------------------------------
InternalGetCountryCodes <- function(# Get ISO 3-digit country codes from character codes and vice versa.
  ### Get vector of 3-digit country codes from character codes and vice versa.
  iso.c,
  countrycodes.csv ##<< csv file with ISO 3-character and 3-digit country codes 
) {
  # change JR, 20140404
  countrycodes <- read.csv(countrycodes.csv, header = T, stringsAsFactors = F)
  countrycodes$ISO.code <- as.character(countrycodes$ISO.code)
  iso.c <- as.character(iso.c)  
  if (all(grepl("[[:digit:]]", iso.c))) {
    iso.c.output <- join(data.frame(ISO.code = iso.c), countrycodes)$Country.letter.code
  } else if (all(grepl("[[:alpha:]]", iso.c))) {
    iso.c.output <- join(data.frame(Country.letter.code = iso.c), countrycodes)$ISO.code
  } else {
    stop("Elements of iso.c must be all alphabet letters or all digits.")
  }
  if (any(is.na(iso.c.output)))
    warning(paste0("The country code(s) ", paste(iso.c[is.na(iso.c.output)], collapse = ", "), 
                   " cannot be found in countrycodes.csv!"))
  return(iso.c.output)
}
#----------------------------------------------------------------------------------
InternalMakeCountryNamesShort <- function(# Shorten country names (and make consistent)
  ### Shorten country names (and make consistent)
  name.c){
    name.c = ifelse(name.c=="Cote d Ivoire", paste("Cote d'Ivoire"), paste(name.c))
    name.c = ifelse(name.c=="Sao Tome and Principe" |name.c=="Sao Tome & Principe", paste("Sao Tome Pr"), paste(name.c))
  name.c = ifelse(name.c=="Vietnam", paste("Viet Nam"), paste(name.c))
  name.c = ifelse(name.c=="Gambia The", paste("Gambia"), paste(name.c))
  name.c = ifelse(name.c=="Brunei Darussalam", paste("Brunei"), paste(name.c))
  name.c = ifelse(name.c=="Saint Kitts and Nevis", paste("Saint Kitts & Nevis"), paste(name.c))
  name.c = ifelse(name.c=="Timor Leste", paste("Timor-Leste"), paste(name.c))
  name.c = ifelse(name.c=="Dominican Rep.", paste("Dominican Republic"), paste(name.c))
  name.c = ifelse(name.c=="Syrian Arab Republic", paste("Syria"), paste(name.c))
  name.c = ifelse(name.c=="GuineaBissau", paste("Guinea-Bissau"), paste(name.c))
  name.c = ifelse(name.c=="Libyan Arab Jamahiriya", paste("Libya"), paste(name.c))
  name.c = ifelse(name.c=="Ukraine ", paste("Ukraine"), paste(name.c))
  name.c = ifelse(name.c=="Republic of Moldova"|name.c=="Moldova, Rep. of", paste("Moldova"), paste(name.c))
  name.c = ifelse(name.c=="Federated States of Micronesia"|name.c=="Micronesia (Federated States of )"|
  name.c=="Micronesia, Federated States of"|name.c=="Micronesia, Fed. States of"|
      name.c=="Micronesia (Fed. States of)", paste("Micronesia"), paste(name.c))
  name.c = ifelse(name.c=="United Kingdom", paste("U.K."), paste(name.c))
  name.c = ifelse(name.c=="United States of America"|name.c=="United States", paste("U.S."), paste(name.c))
  name.c = ifelse(name.c=="Congo, Dem. Rep."|name.c=="Democratic Republic of the Congo"|name.c=="Congo DR", paste("DRC"), paste(name.c))
  name.c = ifelse(name.c=="The former Yugoslav Republic of Macedonia"|name.c=="TFYR Macedonia", paste("Macedonia"), paste(name.c))
  name.c = ifelse(name.c=="Bosnia and Herzegovina"|name.c=="Bosnia & Herzegovina", paste("Bosn&Herze"), paste(name.c))
  name.c = ifelse(name.c=="Trinidad and Tobago"|name.c=="Trinidad & Tobago", paste("Trinidad&T"), paste(name.c))
  name.c = ifelse(name.c=="China, Hong Kong SAR", paste("Hong Kong"), paste(name.c))
  name.c = ifelse(name.c=="United Repulic of Tanzania",   paste("Tanzania"), paste(name.c))
  name.c = ifelse(name.c=="United States Virgin Islands",   paste("US Virgin Isl."), paste(name.c))
  name.c = ifelse(name.c=="United Arab Emirates",   paste("Arab Emirates"), paste(name.c))
  name.c = ifelse(name.c=="Lao People's Democratic Republic"|name.c=="Lao People's Dem. Rep."|
        name.c =="Lao PDR",   paste("Laos"), paste(name.c))
  name.c = ifelse(name.c=="Republic of Korea"|name.c == "Republic of Korea "|
        name.c == "Korea Rep"|name.c=="Korea, Rep. of", paste("South Korea"), paste(name.c))
  name.c = ifelse(name.c=="Democratic People's Republic of Korea"|name.c=="Korea DPR"|
        name.c=="Dem. People's Republic of Korea"|name.c=="Korea, Dem. People's Rep.", paste("North Korea"), paste(name.c))
  name.c = ifelse(name.c=="Central African Republic"|name.c=="Central African Rep.", paste("CAR"), paste(name.c))
  name.c = ifelse(name.c=="Iran (Islamic Republic of)"| name.c == "Iran, Islamic Republic of", paste("Iran"), paste(name.c))
  name.c = ifelse(name.c=="United Republic of Tanzania"| name.c=="Tanzania, United Republic of",   paste("Tanzania"), paste(name.c))
  name.c = ifelse(name.c=="Venezuela (Bolivarian Republic of)",   paste("Venezuela"), paste(name.c))
  name.c = ifelse(name.c=="Bolivia (Plurinational State of)",   paste("Bolivia"), paste(name.c))
  name.c = ifelse(name.c=="Antigua and Barbuda"|name.c=="Antigua & Barbuda", paste("Antigua and B."), paste(name.c))
  name.c = ifelse(name.c=="Northern Mariana Islands", paste("N. Mariana Isl."), paste(name.c))
  name.c = ifelse(name.c=="Occupied Palestinian Territory"|name.c=="OPT", 
                  #paste("Occ. Palestinian Terr."), 
                  paste("OPT"), 
                  paste(name.c))
  name.c = ifelse(name.c=="Saint Vincent and the Grenadines" | 
      name.c=="Saint Vincent & the Grenadines" |name.c=="Saint Vincent/Grenadines"|
      name.c == "St Vincent & the Grenadines", paste("St. Vincent & Gren."), paste(name.c))  
  ##value<< Vector with length of \code{name.c}, but some names replaced
  return(name.c)
}

#----------------------------------------------------------------------------------
InternalWhichReg <- function(# Find region for vector of country iso codes
  ### Find region for vector of country iso codes
  iso.j, ##<< vector of country iso codes
  iso.c, ##<< vector with unique country iso codes
  reg.c##<< vector with region integers, corresponding to iso.c
  ){
  reg.j <- rep(NA, length(iso.j))
  for (iso in unique(iso.j)){
    reg <- reg.c[which.max(iso.c == iso)]
    reg.j[iso.j==iso] <- reg    
  }
  ##value<< vector of length \code{iso.j} with region codes
  return(reg.j)
}

#----------------------------------------------------------------------------------
InternalInternalWhichSubreg <- function(# Find region for vector of country iso codes
  ### Find region for vector of country iso codes
  iso.j, ##<< vector of country iso codes
  iso.c, ##<< vector with unique country iso codes
  subreg.c##<< vector with sub-region integers, corresponding to iso.c
  ){
  subreg.j <- rep(NA, length(iso.j))
  for (iso in unique(iso.j)){
    subreg <- subreg.c[which.max(iso.c == iso)]
    subreg.j[iso.j==iso] <- subreg    
  }
  ##value<< vector of length \code{iso.j} with subregion codes
  return(subreg.j)
}
#----------------------------------------------------------------------------------
GetObsCounts <- function(# Find proportion of countries with 0, 1, nres observations
  ### Find proportion of countries with 0, 1, nres observations
  data.j, ##<< observations, can include NAs
  iso.j, #<< ID for country
  select.j = NULL, ##<< subset of observations to select (NULL = all)
  selected.isos  = NULL,##<< set of iso codes to select (NULL = all)
  nres = 5 ##<< nres gives max obs, e.g. nres = 4 means 4+ is the last open-ended cat
  ## set to max(nres, 5) or min(4, nres)
  ){
  if (is.null(select.j)){
    select.j <- rep(T, length(data.j))
  } 
  if (is.null(selected.isos)){
    selected.isos <- unique(iso.j)
  } 
  # isos for countries with observations in the period
  isos.select <- iso.j[!is.na(data.j) & select.j==TRUE & is.element(iso.j, selected.isos)]
  if (nres >5) nres = 5
  if (nres <4) nres = 4
  if (length(isos.select)==0){
    res <- c(1, rep(0, nres))
    names(res) <- c("0", "1", "2", "3", "4", "5+")[1:(1+nres)]
    return(res)
  } 
  # all countries, with or without observations in that particular period
  isos.all <- iso.j[ is.element(iso.j, selected.isos)]#select.j==TRUE &
  bla <- table(table(isos.select))
  if (nres==4){
    res <- c(length(unique(isos.all)) -  length(unique(isos.select)), # countries without data
             sum(bla["1"], na.rm =T), sum(bla["2"], na.rm =T), sum(bla["3"], na.rm =T),  
             sum(bla[as.numeric(names(bla))>=4], na.rm =T)
             )
    names(res) <- c("0", "1", "2", "3", "4+")
  } else {
    res <- c(length(unique(isos.all)) -  length(unique(isos.select)),
             sum(bla["1"], na.rm =T), sum(bla["2"], na.rm =T), sum(bla["3"], na.rm =T), 
             sum(bla["4"], na.rm =T),
             sum(bla[as.numeric(names(bla))>=5], na.rm =T)
             )
    names(res) <- c("0", "1", "2", "3", "4", "5+")
  }
  ##value<< vector with number of countries with 0,..,nres observations
  return(res/sum(res, na.rm = T))
}
#----------------------------------------------------------------------------------
InternalPlotPropsDataAvailability <- function(# Plot proportions of data availability
  ### Internal function for plotting proportions of data availability
  ymin, ##<< which height should stuff be plotted?
  res, ##<< proportions to be plotted
  col = 2, ##<< color to visualize proportion
  include.zero ##<< should first box with zero be included?
  ){
  yvalues <- ymin + res
  istart <- ifelse(include.zero,1,2)
  for (i in istart:length(res)){
    polygon(c(i-1,i-1, i,i,i-1), 
            c(ymin,yvalues[i], yvalues[i], ymin, ymin),
            col = col, border = NA)
    # get box back
    polygon(c(i-1,i-1, i,i,i-1), 
            c(ymin,ymin+1, ymin+1, ymin, ymin),
            bg = NA, border = 1)
    text(i-0.5,ymin+0.55, labels = round(100*res[i]), cex =1, col = 1)
  }
}

#----------------------------------------------------------------------------------
PlotDataAvailability <- function(# Create a plot to visualize data availability.
  ### Create a plot to visualize data availability.
  data, ##<< Object of \code{data} 
  country.info, ##<< Object of \code{country.info}
  nres = 5, ##<< Determines the number of columns (nres+ is the last column, max is 5).
  summarize.unmet = FALSE ##<< Logical: total prevalence is summarized when FALSE, unmet need when TRUE.
  ){
  
  if (!summarize.unmet){
    data.j <- data$props.tot.j
  } else {
    data.j <- data$props.unmet.j
  }
  iso.j  <- data$iso.j
  
  # get regional info
  isos.Lr <- list()
  isos.Lr[["All"]] = paste(country.info$iso.c)
  isos.Lr[["Developed"]] = paste(country.info$iso.c[country.info$dev.c=="Rich"])
  isos.Lr[["Africa"]] = paste(country.info$iso.c[country.info$dev.c!="Rich"&
    country.info$namereg.c=="Africa"])
  isos.Lr[["Asia"]] =paste(country.info$iso.c[country.info$dev.c!="Rich"&
    country.info$namereg.c=="Asia"])
  isos.Lr[["LAC"]] =paste(country.info$iso.c[country.info$dev.c!="Rich"&
    country.info$namereg.c=="Latin America and the Caribbean"])
  isos.Lr[["Oceania"]] =paste(country.info$iso.c[country.info$dev.c!="Rich"&
    country.info$namereg.c=="Oceania"])
  regs <- names(isos.Lr)
  nc.r <- list() # number of countries in each region
  for (reg in regs){
    nc.r[[reg]] <- length(unique(isos.Lr[[reg]]))
  }
  
  include.zero = TRUE
  #res <- GetObsCounts(data.j = data.j, iso.j = iso.j)
  # could be that break down by region+period introduces zeroes
  #if (res[1]>0 ) include.zero = TRUE
  
  par(mfrow = c(1,1), mar = c(5,10,1,1), cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
  plot(1, type = "n", ylim = c(0, -1+length(regs)*6), bty = "n", 
       xlim = c(as.numeric(!include.zero),nres+1), xaxt = "n", yaxt = "n", 
       ylab = "", xlab = "Number of observations",
       main = "")
  axis(1, at = 0.5+seq(0,nres), labels = c(seq(0,nres-1), paste(nres, "+", sep = "")))
  # year is the cutoff moment (1990.0 is in second row)
  years <- c(0,1990, 2000, 2005, Inf)
  yearnow <- max(floor(data$years.j))+1 # note: 2011 should give 2012 as the end point
  labelsperiods <- c("Before 1990", "1990-1999", "2000 - 2004", paste("2005 -", yearnow))
  ymin <- ystart <- -1 + length(regs)*6
  namesregs <- paste(regs, " (", unlist(nc.r), ")", sep = "")
  for (reg in 1:length(regs)){
    axis(2, at = -1+ymin-seq(1,5), labels = rep("",5), las = 1, tick=TRUE)
    axis(2, at = -1.5+ymin-seq(1,4), labels = labelsperiods, las = 1, tick=FALSE)
    #  axis(2, at = -1+ymin-seq(1,5), labels = c("", years[-1]), las = 1)
    axis(2, at = -1+ymin-0.5, labels = paste(namesregs[reg], "", sep = ""), las = 1)
    ymin <- ymin-2
    region <- regs[reg]
    iso.select <- isos.Lr[[region]] 
    res <- GetObsCounts(data.j = data.j, iso.j = iso.j, 
                        selected.isos = iso.select)
    InternalPlotPropsDataAvailability(col = 3,ymin = ymin, res = res, include.zero = include.zero)
    for (t in 2:length(years)){
      ymin <- ymin-1
      res <- GetObsCounts(data.j = data.j, iso.j = iso.j, 
                          select.j = (data$years.j>=years[t-1] & data$years.j<years[t]),
                          selected.isos = iso.select)
      InternalPlotPropsDataAvailability(col = "lightblue",ymin = ymin, res = res, include.zero = include.zero)
    }
  }# end plot
  
  if (include.zero) segments(1, -1, 1, ystart-1, lwd = 5)
  ##value<< NULL. Plot appears in R-console.
  return(invisible())
}

#----------------------------------------------------------------------------------
# The End.
