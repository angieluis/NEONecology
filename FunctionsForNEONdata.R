###############################################################################
## Functions for estimating small mammal species-specific capture rates
###############################################################################

library(tidyverse)


## To do: 
# Could add "yearmon" option to CH instead of prim.session
# or update "yearmon" so adjusts based on prim.session?

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

  }

site.date.function <- function()


## Function to look for duplicates ------------------------------------------ ##

# does not show first occurrence
dupes.fun <- function(data, ...) {
  data %>%
    group_by(.dots = lazyeval::lazy_dots(...)) %>%
    mutate(n = n()) %>%
    filter(n() > 1)
}



## Multisite primary session function ----------------------------------- ##

NEON.primary.session.function <- function(site.dates = smammal.sitedates,# a list within a list 
                                          int.break = 10){  #if interval is greater than 10 it's a different session
  primary.session.list <- list()
  counter <- 1
  for(i in 1:length(site.dates)){
    datesite <- site.dates[[i]]
    for(p in 1:length(site.dates[[i]])){
      ints <- diff(datesite[[p]]) #gives intervals
      ints <- c(ints, int.break+1) #add extra one on end so counts last interval
      breaks <- which(ints>int.break)
      num <- length(breaks)
      session.num <- rep(1:num, times=diff(c(0,breaks)))
      sessiondf <- data.frame(date = datesite[[p]], prim.session = session.num)
      yearmon <- character()
      for(m in 1:length(unique(sessiondf$prim.session))){
        d <- gsub("-", "", as.character(sessiondf$date[max(which(sessiondf$prim.session==m))]))
        yearmon[which(sessiondf$prim.session==m)] <- unique(substr(d,start=1,stop=6))
      }
      sessiondf$yearmon <- yearmon
      primary.session.list[[counter]] <- sessiondf
      names(primary.session.list)[counter] <- names(datesite)[p]
      counter <- counter + 1
    }
  }
  return(primary.session.list)
}
# returns a list of dataframes (for each plot) that contain dates trapped, primary session (1:m), and yearmon
# yearmon is still a problem, when there was trapping early in the month and late in the month
# prim.session is better. e.g., 
#          date prim.session yearmon
# 4  2019-07-30            2  201908
# 5  2019-07-31            2  201908
# 6  2019-08-01            2  201908
# 7  2019-08-27            3  201908
# 8  2019-08-28            3  201908
# 9  2019-08-29            3  201908
# But then prim.session numbers won't line up across sites over time.


## function to create capture history as long data ------------------------- ##

# Function to get longdata format capture history for estimating recapture
# rates, only using sessions when an animal was caught at least once.
# Keeping as longdata but need to add rows to data for when an individual
# was caught within that primary session but not a secondary occasion
# add a State=0 for not caught 
# This is for estimating recapture rates only, with grouping columns assumed
# to be in the data (e.g., column called species).
# This will not work for survival estimation (because only have sessions animals were caught).
# This is using prim.session from the primary session list from above function
# not yearmon. Could update code so could use either.
# Because of the messy data - sometimes same tag same date, same plot, multiple species
# take the first duplicate row per plot,date,ID as covariate info (not all that many)
NEON.pRDcapture.history.long.fun <- function(
    cleaned.data = reduced.smammal.captures, # 'cleaned' needs column 'plotID', 
    ID.col, # column name of unique individual ID
    prim.session.list = reduced.prim.session.list, # list with plots, dates, and prim.session
    #session.type = "prim.session", # alternative is "yearmon". see above function --- not implemented
    cols = c("genus", "species"), # which column names to keep as individual covariates (in addition to site, plot, taxon)
    taxon = "all") { # taxonID (4_letter code) or "all" 
  
  names(cleaned.data)[which(names(cleaned.data)==ID.col)] <- "site_tag"
  
  plots <- names(prim.session.list)
  prim.session.long <- data.frame(plotID = character(0),
                                  date = Date(0),
                                  prim.session = numeric(0),
                                  Sec = numeric(0))
  
  for(i in 1:length(prim.session.list)){
    dat <- data.frame(plotID = names(prim.session.list)[i],
                      prim.session.list[[i]])
    sec <- numeric()
    for(s in unique(dat$prim.session)){
      sec <- c(sec, 1:length(which(dat$prim.session==s)))
    }
    dat$Sec <- sec
    
    prim.session.long <- rbind(prim.session.long, dat)
  }
  
  if(length(which(names(cleaned.data)=="date"))==0){
    cleaned.data$date <- cleaned.data$collectDate
  }
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
  
  long.dat <- left_join(sp.data, prim.session.long)
  
  # keep individual-level data we care about for each prim.session
  prim.ses.dat <- long.dat %>%
    group_by(siteID, plotID, site_tag, prim.session, across(all_of(cols))) %>%
    summarise(n=n())

  
  eg <- expand.grid(site_tag = unique(sp.data$site_tag[which(sp.data$plotID==plots[1])]), 
           date = prim.session.long$date[which(prim.session.long$plotID==plots[1])])
  eg$plotID <-  plots[1]
  for(p in 2:length(plots)){
    egn <- expand.grid(site_tag = unique(sp.data$site_tag[which(sp.data$plotID==plots[p])]), 
                       date = prim.session.long$date[which(prim.session.long$plotID==plots[p])])
    egn$plotID <- plots[p]
    eg <- rbind(eg, egn)
  }
  # eg$State <- 0 # State == 0 for not caught # Wait - is this overwriting those that were caught?
  
  full.session.dat <- left_join(eg, prim.session.long) 
  
  full.session.dat <- full_join(full.session.dat, prim.ses.dat, relationship="many-to-one")

  #long.dat$State[which(is.na(long.dat$State))] <- 0
  
  long.dat <- long.dat %>%
    select(siteID, plotID, site_tag, all_of(cols), date, prim.session, 
           #yearmon, 
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

taxa.capture.fun <- function(
  captures = reduced.smammal.captures, 
  prim.session.list = reduced.prim.session.list,
  group = "species", # "species" or "genus" or other group identified in captures
  taxa.to.include = NULL, #  if NULL, then use all. Otherwise a vector of 4-letter species codes
  wide.or.long = "wide"){ # output as wide (with columns for taxa) or long?
  # (or other levels of group) to include
  
  
  if(length(taxa.to.include)==0){
    taxa.to.include <- sort(as.character(unique(as.data.frame(reduced.smammal.captures)[ , group])))
  } 
  n.taxa <- length(taxa.to.include)
  
  if(length(which(names(captures)=="date"))==0){
    captures$date <- captures$collectDate
  }
  
  prim.session.long <- data.frame(plotID = character(0),
                                  date = Date(0),
                                  prim.session = numeric(0),
                                  Sec = numeric(0))
  
  for(i in 1:length(prim.session.list)){
    dat <- data.frame(plotID = names(prim.session.list)[i],
                      prim.session.list[[i]])
    sec <- numeric()
    for(s in unique(dat$prim.session)){
      sec <- c(sec, 1:length(which(dat$prim.session==s)))
    }
    dat$Sec <- sec
    
    prim.session.long <- rbind(prim.session.long, dat)
  }
  
  captures$n <- 1
  
  out.data <- left_join(captures, prim.session.long)
  
  out.data <- out.data %>%
    group_by(siteID, plotID, prim.session, species) %>%
    count(plot_taxon_tag) %>%
    count(species)

  prim.ses <- prim.session.long %>%
    group_by(plotID, prim.session) %>%
    summarise(first.dat = min(date),
              n.sec.occ = length(unique(date)))
  
  out.data <- left_join(out.data, prim.ses)
  
  if(wide.or.long == "wide"){
    out.wide <- out.data %>%
       pivot_wider(
         names_from = all_of(group),
         values_from = n,
         values_fn = ~ sum(.x, na.rm=TRUE))
    out.wide <- replace(out.wide, is.na(out.wide), 0)
    out.wide <- out.wide[, c(1:5, order(names(out.wide)[6:dim(out.wide)[2]])+5)]
    return(out.wide)
  }
  if(wide.or.long == "long"){
    return(out.data)
  }
  
}


#<------------------------------------------------------------------------------------- Stopped here
## Function to estimate population abundance -------------------------------- ##
# given captures and species-specific capture rates


NEON.Nestimates.fromcaptures <- function(num.captures = smammal.species.num.captures.widedata, #output from taxa.capture.fun function
                                         p.estimates = p.species.estimates, # summary of Bayesian estimates as data.frame
                                         group = "species", # same label as column in p.estimates
                                         prim.session.list = reduced.prim.session.list,
                                         trapStatuslist = trapStatuslist){ # do I need to reduce this?
  
  mat.ind <- which(names(num.captures)=="n.sec.occ") + 1
  sp.caps <- num.captures[, mat.ind:(dim(num.captures)[2])]
  restof.data <- num.captures[, 1:(mat.ind-1)]
  
  trap.nights <- numeric()
  for(j in 1:dim(restof.data)[1]){
    site <- str_split_1(restof.data$plot[j], "_")[1]
    prim.sessions <- prim.session.list[[restof.data$plot[j]]]
    session.dates <- prim.sessions$date[which(prim.sessions$prim.session == restof.data$prim.session[j])]
    tsl <- trapStatuslist[[site]][[restof.data$plot[j]]]
    trap.nights[j] <- sum(tsl$traps_available[match(session.dates, tsl$date)])
  }
  restof.data$trapnights <- trap.nights
  
  n.species <- dim(sp.caps)[2]
  p.estimates <- p.estimates[match(names(sp.caps), p.estimates[,group] ),] ###########################
  p <- p.estimates$mean  
  p.effmat <- matrix(NA, nrow=dim(sp.caps)[1], ncol=dim(sp.caps)[2])
  for(i in 1:n.species){
    for(t in 1:dim(p.effmat)[1]){
      p.effmat[t, i] <- 1 - (1 - p[i])^num.captures$n.sec.occ[t]
    }
  }
  
  N.estmat <- sp.caps / p.effmat / (restof.data$trapnights / restof.data$n.sec.occ / 100)
  Nestdf <- data.frame(restof.data, N.estmat, check.names = FALSE)
  
  
  return(Nestdf)
}

