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

  return(traps.avail)
  
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

# as long data
NEON.session.function <- function(trap.availability.longdata = trapping.info, # output from trap.status.function
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
# Keeping as longdata but need to add rows to data for when an individual
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

taxa.capture.fun <- function(
  captures = reduced.smammal.captures, 
  #prim.session.list = reduced.prim.session.list,
  session.info = reduced.trapping.session.info,
  group = "species", # "species" or "genus" or other group identified in captures
  taxa.to.include = NULL, #  if NULL, then use all. Otherwise a vector of 4-letter species codes (or other levels of group) to include
  wide.or.long = "wide"){ # output as wide (with columns for taxa) or long?
  
  
  
  if(length(taxa.to.include)==0){
    taxa.to.include <- sort(as.character(unique(as.data.frame(reduced.smammal.captures)[ , group])))
  } 
  n.taxa <- length(taxa.to.include)
  
  # if(length(which(names(captures)=="date"))==0){
  #   captures$date <- captures$collectDate
  # }
  # 
  # prim.session.long <- data.frame(plotID = character(0),
  #                                 date = Date(0),
  #                                 prim.session = numeric(0),
  #                                 Sec = numeric(0))
  # 
  # for(i in 1:length(prim.session.list)){
  #   dat <- data.frame(plotID = names(prim.session.list)[i],
  #                     prim.session.list[[i]])
  #   sec <- numeric()
  #   for(s in unique(dat$prim.session)){
  #     sec <- c(sec, 1:length(which(dat$prim.session==s)))
  #   }
  #   dat$Sec <- sec
  #   
  #   prim.session.long <- rbind(prim.session.long, dat)
  # }
  
  captures$n <- 1
  
  out.data <- left_join(captures, session.info)
  
  out.data <- out.data %>%
    group_by(siteID, plotID, prim.session, species) %>%
    count(plot_taxon_tag) %>%
    count(species)

  prim.ses <- session.info %>%
    group_by(plotID, prim.session) %>%
    summarise(first.dat = min(collectDate),
              n.sec.occ = length(unique(collectDate)),
              trapnights = sum(traps_available))
  
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

