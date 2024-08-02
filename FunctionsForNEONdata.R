###############################################################################
## Functions for estimating small mammal species-specific capture rates
###############################################################################

library(tidyverse)


## To do: 
# need tagID in CHlong state==0
# Move adding the Secondary occasion to the primary.session.function
#    instead of within the capture history function.
# Could add "yearmon" option to CH instead of prim.session
# or update "yearmon" so adjusts based on prim.session?



## Function to look for duplicates ------------------------------------------ ##


# does not show first occurrence
dupes.fun <- function(data, ...) {
  data %>%
    group_by(.dots = lazyeval::lazy_dots(...)) %>%
    mutate(n = n()) %>%
    filter(n() > 1)
}



## Multisite primary session function ----------------------------------- ##

NEON.primary.session.function <- function(site.dates = reducedSWsitedates,# a list within a list 
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
  
  
  