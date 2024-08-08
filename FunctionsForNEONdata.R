###############################################################################
## Functions for estimating small mammal species-specific capture rates
###############################################################################

library(tidyverse)

## Logit function ----------------------------------------------------------- ##


logit <- function(x) {
  log(x / (1 - x))
}


# Reverse logit funciton ---------------------------------------------------- ##


rev.logit <- function(x) {
  exp(x) / (1 + exp(x))
}


## Function to get trapping info including date & number traps set ---------- ##

trap.status.function <- function(captures = NEONsmammalcaptures){
  
  statuslong <- NEONsmammalcaptures %>%
    # since trapStatus=="4 - more than 1 capture in one trap" will have multiple rows per trap, 
    # I need to also include trapCoordinate so can make sure don't count same trap multiple times
    mutate(trapStatus = ifelse(trapStatus=="4 - more than 1 capture in one trap", 
                                paste(trapStatus, trapCoordinate, sep="_"), trapStatus)) %>%
    group_by(siteID, plotID, collectDate) %>%
    count(trapStatus) %>%
    mutate(trapStatus2 = as.numeric(str_trunc(trapStatus, width=1, ellipsis = "")))
    
  traps.avail <- statuslong %>%
    group_by(siteID, plotID, collectDate) %>%
    summarise(traps_available = sum(n[trapStatus2 %in% c(2,3,5,6)]) + 
                length(which(trapStatus2 == 4))) %>% 
    # group trapStatus2 == 2-6 are available, [1 is not set] but trapStatus2 == 4 is more than 1 
    # capture per trap - from those, calculate number of unique trapCoordinates, which are rows
    filter(traps_available != 0) # remove all the rows with 0 traps set

  return(traps.avail)
  
  }




## Function to look for duplicates ------------------------------------------ ##

# does not show first occurrence
dupes.fun <- function(data, ...) {
  data %>%
    group_by(.dots = lazyeval::lazy_dots(...)) %>%
    mutate(n = n()) %>%
    filter(n() > 1)
}



## Multisite primary session function ----------------------------------- ##

# as long data
NEON.session.function <- function(
    trap.availability.longdata = trapping.info, # output from trap.status.function
    int.break = 10){  #if interval is greater than 10 days it's a different primary session
  
  data <- trap.availability.longdata %>%
    group_by(siteID, plotID) %>%
    mutate(ints = c( diff.Date(collectDate), int.break + 1 )) %>%
    mutate(prim.session = rep(1:length(which(ints>int.break)), times=diff(c(0,which(ints>int.break))))) %>%
    group_by(plotID, prim.session) %>%
    mutate(Sec = row_number()) %>%
    ungroup()
  
   return(data)
}
# returns a long dataframe that contain dates trapped, number of traps set,
# time until next trapping day (with 11 on the last occasion because that's 1 more 
# than the default break. will need to be removed if do anything with intervals)
# Also primary occasion starting at 1 for first primary sessoin per plot
# Sec is secondary occasion or day
#   siteID plotID   collectDate         traps_available ints     prim.session   Sec
#   <chr>  <chr>    <dttm>                        <int> <drtn>          <int> <int>
# 1 ABBY   ABBY_004 2021-07-06 00:00:00             100   1 days           17     1
# 2 ABBY   ABBY_004 2021-07-07 00:00:00             100   1 days           17     2
# 3 ABBY   ABBY_004 2021-07-08 00:00:00             100  25 days           17     3
# 4 ABBY   ABBY_004 2021-08-02 00:00:00              98   1 days           18     1
# 5 ABBY   ABBY_004 2021-08-03 00:00:00              98   1 days           18     2
# 6 ABBY   ABBY_004 2021-08-04 00:00:00             100  28 days           18     3
# 7 ABBY   ABBY_004 2021-09-01 00:00:00              99   1 days           19     1
# 8 ABBY   ABBY_004 2021-09-02 00:00:00              99   1 days           19     2
# 9 ABBY   ABBY_004 2021-09-03 00:00:00             100 290 days           19     3





## function to create capture history as long data ------------------------- ##

# Function to get longdata format capture history for estimating recapture
# rates, only using sessions when an animal was caught at least once.
# This keeps capture info as longdata but adds rows to data for when an individual
# was caught within that primary session but not a secondary occasion
# add a State=0 for not caught 
# This is for estimating recapture rates only, with grouping columns assumed
# to be in the data (e.g., column called species).
# This will not work for survival estimation (because only have sessions animals were caught).
# Because of the messy data - sometimes same tag same date, same plot, multiple species
# take the first duplicate row per plot,date,ID as covariate info (not all that many)
NEON.pRDcapture.history.long.fun <- function(
    cleaned.data = reduced.smammal.captures, # 'cleaned' needs column 'plotID', 
    ID.col = "plot_taxon_tag", # column name of unique individual ID
    session.info = reduced.trapping.session.info, # output from NEON.session.function
    #prim.session.list = reduced.prim.session.list, # list with plots, dates, and prim.session
    cols = c("genus", "species"), # which column names to keep as individual covariates (in addition to site, plot, taxon)
    taxon = "all") { # taxonID (4_letter code) or "all" 
  
  names(cleaned.data)[which(names(cleaned.data)==ID.col)] <- "site_tag"
  
  plots <- sort(unique(session.info$plotID))

  cleaned.data <- cleaned.data %>%
    filter(plotID %in% plots)
  
  # separate out species wanted, if applicable  
  if(taxon != "all"){
    sp.data <- cleaned.data %>%
      filter(taxonID == species) 
  } else {
    sp.data <- cleaned.data
  }
  
  sp.data$State <- 1 # State == 1  for caught
  
  long.dat <- left_join(sp.data, session.info)
  
  # keep individual-level data we care about for each prim.session
  prim.ses.dat <- long.dat %>%
    group_by(siteID, plotID, site_tag, prim.session, across(all_of(cols))) %>%
    summarise(n=n()) %>%
    select(!n)

  
  eg <- expand.grid(site_tag = unique(sp.data$site_tag[which(sp.data$plotID==plots[1])]), 
           collectDate = session.info$collectDate[which(session.info$plotID==plots[1])])
  eg$plotID <-  plots[1]
  for(p in 2:length(plots)){
    egn <- expand.grid(site_tag = unique(sp.data$site_tag[which(sp.data$plotID==plots[p])]), 
                       collectDate = session.info$collectDate[which(session.info$plotID==plots[p])])
    egn$plotID <- plots[p]
    eg <- rbind(eg, egn)
  }

  full.session.dat <- left_join(eg, session.info) 
  
  full.session.dat <- full_join(full.session.dat, prim.ses.dat, relationship="many-to-one")

  long.dat <- long.dat %>%
    select(siteID, plotID, site_tag, all_of(cols), collectDate, prim.session, 
           Sec, State)

  long.dat <- full_join(long.dat, full.session.dat)
  long.dat$State[which(is.na(long.dat$State))] <- 0

  keepIDses <- long.dat %>%
    group_by(plotID, site_tag, prim.session) %>%
    summarise(Here = ifelse(sum(State) == 0, 0, 1))

  long.dat <- left_join(long.dat, keepIDses) %>%
    filter(Here==1) ### could leave in Here == 0 later if want to do survival estimation
  
  return(long.dat)
}
  
  
  
## function to calculate number of captures per any taxa (group) ----------- ##
# per primary occasion, and trapnights


taxa.capture.fun <- function(
    captures = smammal.captures.clean, 
    session.info = trapping.session.info,
    group = "species", # "species" or "genus" or other group identified in captures
    ID.col = "plot_taxon_tag", # individual identifier as character vector, e.g., "tagID"
    taxa.to.include = NULL, #  if NULL, then use all. Otherwise a vector of 4-letter species codes (or other levels of group) to include
    include.hanta = TRUE, # if TRUE then need output to be longdata to see it. 
    wide.or.long = "long"){ # output as wide (with columns for taxa) or long?
  
  names(captures)[which(names(captures)==ID.col)] <- "ID"
  names(captures)[which(names(captures)==group)] <- "group.tmp"
  
  if(length(taxa.to.include)==0){
    taxa.to.include <- sort(as.character(unique(as.data.frame(captures)[ , "group.tmp"])))
  } 
  
  captures$n <- 1
  
  out.data <- left_join(captures, session.info)
  
  # need to add option for hantadata or other covariates
  out.data <- out.data %>%
    group_by(siteID, plotID, prim.session, group.tmp) %>% 
    count(ID) %>% 
    count(group.tmp) 
  
  if(include.hanta == TRUE){
    hanta.data <- left_join(captures, session.info) %>%
      group_by(siteID, plotID, prim.session, group.tmp) %>%
      summarise(hantaPos = length(which(hantaResult=="Positive")),
                hantaTested = length(which(hantaResult == "Positive" | hantaResult == "Negative")))
    out.data <- left_join(out.data, hanta.data)
  }
  
  prim.ses <- session.info %>%
    group_by(plotID, prim.session) %>%
    summarise(first.date = min(collectDate),
              n.sec.occ = length(unique(collectDate)),
              trapnights = sum(traps_available))
  
  out.data <- left_join(out.data, prim.ses)
  
  if(wide.or.long == "wide"){
    out.wide <- out.data %>%
      pivot_wider(
        names_from = group.tmp,
        values_from = n,
        values_fn = ~ sum(.x, na.rm=TRUE),
        values_fill = 0,
        names_sort = TRUE)
    
    #gr <- sort(unique(out.data$group.tmp))
    #inds <- match(gr, names(out.wide))
    #out.wide <- out.wide[ , c(1:(min(inds)-1), inds)]
    return(out.wide)
  }
  if(wide.or.long == "long"){
    names(out.data)[which(names(out.data)=="group.tmp")] <-  group
    return(out.data)
  }
  
}


## Functions to estimate population abundance per primary session ----------- ##
# given captures and species-specific (and trapnight adjusted) capture rates 


# this only works for species model, that assumes number of traps don't matter.
# p.species.estimates needs to have column "species", and "intercept", 
NEON.N.estimates.fromcaptures.speciesmodel <- function(
    num.captures = smammal.species.num.captures.longdata, #output from taxa.capture.fun function in long format
    p.param.estimates = p.species.estimates, # summary of Bayesian estimates as data.frame
    group = "species", # same label as column in p.param.estimates
    #model = "species-trap", # options are "species-trap" or "species" model. see model run code for model specification
#    calc.CI = TRUE, # calculate upper and lower CIs based on CIs of estimates
    model = expression(), # using column names in p.param.estimates (and covariates in num.captures)
    wide.or.long = "long"){ 
  
  names(num.captures)[which(names(num.captures)=="n")] <- "n_caps"
  
  data <- left_join(num.captures, p.estimates)
  
  data <- data %>%
      mutate(p_night = intercept,
             p_eff = 1 - (1 - p_night)^n.sec.occ) 
  
  
  data <- data %>%
    mutate(N.est = n_caps / p_eff)
 
  if(wide.or.long == "wide"){
    out.wide <- data %>%
      pivot_wider(
        names_from = all_of(group),
        values_from = N.est,
        values_fn = ~ sum(.x, na.rm=TRUE),
        values_fill = 0,
        names_sort = TRUE)
    #inds <- match(p.param.estimates[group][,1], names(out.wide))
    #out.wide <- out.wide[ , c(1:(min(inds)-1), inds)]
    return(out.wide)
  }
  if(wide.or.long == "long"){
    return(data)
  }
  
}


# try specifying model by expression
# and eval(expr, envir) # expression and named data


# Michaelis-Menten model where estimate is half saturation constant of 
# the effect of traps
# columns in p.param.estimates are 'species', 'trap_k' (estimate),
# 'trap_kLCI' and 'trap_kUCI' for lower and upper CI of estimate 
# See "ModelRunSpeciesTrapMMCaptureRateEstimations.R"
NEON.N.estimates.fromcaptures.MMmodel <- function(
    num.captures = smammal.species.num.captures.longdata, #output from taxa.capture.fun function in long format
    p.param.estimates = p.speciestrapMM.model.estimates, # summary of Bayesian estimates as data.frame
    group = "species", # same label as column in p.param.estimates
    wide.or.long = "long"){ 
  
  names(num.captures)[which(names(num.captures)=="n")] <- "n_caps"
  
  data <- left_join(num.captures, p.param.estimates)
  

  data <- data %>%
     mutate(
       trapnights_transformed = trapnights/100, # reverse of transformation before model
       p_night = (trapnights_transformed/n.sec.occ) / (trap_k + (trapnights_transformed/n.sec.occ)), # assuming trapnights evenly among secondary nights of primary occasion
       p_night_upper = (trapnights_transformed/n.sec.occ) / (trap_kUCI + (trapnights_transformed/n.sec.occ)),
       p_night_lower = (trapnights_transformed/n.sec.occ) / (trap_kLCI + (trapnights_transformed/n.sec.occ)),
       p_eff = 1 - (1 - p_night)^n.sec.occ,
       p_eff_upper = 1 - (1 - p_night_upper)^n.sec.occ,
       p_eff_lower = 1 - (1 - p_night_lower)^n.sec.occ,
       N.est = n_caps / p_eff,
       N.est_upper = n_caps / p_eff_upper,
       N.est_lower = n_caps / p_eff_lower)
  
  if(wide.or.long == "wide"){ # doesn't include confidence intervals
    out.wide <- data %>%
      pivot_wider(
        names_from = all_of(group),
        values_from = N.est,
        values_fn = ~ sum(.x, na.rm=TRUE),
        values_fill = 0,
        names_sort = TRUE)
    #inds <- match(p.param.estimates[group][,1], names(out.wide))
    #out.wide <- out.wide[ , c(1:(min(inds)-1), inds)]
    return(out.wide)
  }
  if(wide.or.long == "long"){
    return(data)
  }
  
}

