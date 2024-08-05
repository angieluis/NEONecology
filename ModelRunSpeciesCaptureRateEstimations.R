
############################################################################
# Model specification and run code for estimating species specific capture 
# rates using closed Bayesian mark-recapture models from capture data 

# Only covariate is species. These all ignore differences in number of 
# traps set per night --------------------------------------------------- ##
# assumes all nights all plots are the same ----------------------------- ##
# So not a good model


# see "FunctionsForNEONdata.R" for data functions
# see "CleanOrganizeSmammalCaptureData.R" for cleaning and formatting data
############################################################################

# source in functions
source("FunctionsForNEONdata.R")

# load data formatted for these models
load("NEONchForPmodels.RData")


# Need to add group numbers for the taxon level I'm interested in:
NEON.CHlong.all$species.num <- taxon.cap.sum$species.num[match(NEON.CHlong.all$species,
                                                               taxon.cap.sum$species)]
NEON.CHlong.all$genus.num <- taxon.cap.sum$genus.num[match(NEON.CHlong.all$genus,
                                                               taxon.cap.sum$genus)]



################################################################################
## Specify the Robust Design monthly Closed Model---------------------------- ##
################################################################################


sink("RDp_species.bug")

cat("
     
     model {
       
      ##### PRIORS AND CONSTRAINTS #####
      
      
      ##### PRIORS #####
      for(s in 1:n.group){
        p[s] ~ dunif(0, 1)    #  capture rate per group
      }    

        
        ##### LIKELIHOOD #####
        
        ##### OBSERVATION PROCESS
        for(obs in 1:n.obs) {
          y[obs] ~ dbern(p[group[obs]]) # assumes that everything in obs only includes sessions that animal was caught
          } # obs
        
        } # model
    ", fill = TRUE)

sink()


## -------------------------------------------------------------------------- ##



##############################################################################
## Run the Model for species
##############################################################################



# packages
library(R2jags)


## Bundle data -------------------------------------------------------------- ##

# list of bugs data
bugs.data <- list(
  n.obs                 = length(NEON.CHlong.all$State),
  y                     = NEON.CHlong.all$State, 
  n.group               = max(as.numeric(NEON.CHlong.all$species.num)), # group is species for this model
  group                 = as.numeric(NEON.CHlong.all$species.num) 
)


# supply initial values
inits <- function() {
  list(
    p     = runif(bugs.data$n.group, 0, 1)
  )
}

# parameters to monitor
parameters <- c("p")

# run the model
date()
NEONAllsites.species.p.mod <- jags(data = bugs.data,
                                   inits, 
                                   parameters, 
                                   "RDp_species.bug", 
                                   n.chains = 3, 
                                   n.thin   = 5, 
                                   n.iter   = 10000, 
                                   n.burnin = 2000)
date() # 1.5 minutes

save(NEONAllsites.species.p.mod, file="SmammalSpeciesPestModOutput.RData")

p.species.estimates <- data.frame(species=taxon.cap.sum$species, 
  NEONAllsites.species.p.mod$BUGSoutput$summary[2:(length(unique(NEON.CHlong.all$species))+1), ])

summary(p.species.estimates)


# Make some plots -----------------------------------------------------------##

library(mcmcplots)
library(MCMCvis)

mcmcplot(NEONAllsites.species.p.mod,
         parms = parameters)

pdf(file="SmammalSpeciesCaptureRates.pdf", width=8, height=20)
MCMCplot(NEONAllsites.species.p.mod,
         xlim = c(0, 1),
         sz_labels = 0.7,
         ref_ovl = TRUE,
         params = parameters,
         labels = taxon.cap.sum$species,
         main="Species-specific capture rates for NEON data")
dev.off()




##############################################################################
## Run the Model for Genus
##############################################################################

# packages
library(R2jags)


## Bundle data -------------------------------------------------------------- ##

# list of bugs data
bugs.data <- list(
  n.obs                 = length(NEON.CHlong.all$State),
  y                     = NEON.CHlong.all$State, 
  n.group               = max(as.numeric(NEON.CHlong.all$genus.num)), # group is genus for this model
  group                 = as.numeric(NEON.CHlong.all$genus.num) 
)


# supply initial values
inits <- function() {
  list(
    p     = runif(bugs.data$n.group, 0, 1)
  )
}

# parameters to monitor
parameters <- c("p")

# run the model
date()
NEONAllsites.genus.p.mod <- jags(data = bugs.data,
                                   inits, 
                                   parameters, 
                                   "RDp_species.bug", 
                                   n.chains = 3, 
                                   n.thin   = 5, 
                                   n.iter   = 10000, 
                                   n.burnin = 2000)
date() # 1.5 minutes

save(NEONAllsites.genus.p.mod, file="SmammalGenusPestModOutput.RData")

# labels are not alphabetical but by genus.number so need to match back up with info
p.genus.estimates <- data.frame(genus=taxon.cap.sum$genus[match(1:bugs.data$n.group, taxon.cap.sum$genus.num)], 
                                NEONAllsites.genus.p.mod$BUGSoutput$summary[2:(length(unique(NEON.CHlong.all$genus))+1), ])

summary(p.genus.estimates)


# Make some plots -----------------------------------------------------------##

library(mcmcplots)
library(MCMCvis)

mcmcplot(NEONAllsites.genus.p.mod,
         parms = parameters)

pdf(file="SmammalGenusCaptureRates.pdf", width=8, height=10)
MCMCplot(NEONAllsites.genus.p.mod,
         xlim = c(0, 1),
         sz_labels = 0.9,
         ref_ovl = TRUE,
         params = parameters,
         labels = taxon.cap.sum$genus[match(1:bugs.data$n.group, taxon.cap.sum$genus.num)],
         main="Genus-specific capture rates for NEON data")
dev.off()



###############################################################################
## Use estimates to adjust captures for abundance estimates
###############################################################################

# should do this for all captures not just the reduced ones 
# ie. also include those only trapped 1 night per session

smammal.species.num.captures.widedata <- taxa.capture.fun(
    captures = reduced.smammal.captures, 
    session.info = reduced.trapping.session.info,
    group = "species", # "species" or "genus" or other group identified in captures
    taxa.to.include = NULL, #  if NULL, then use all. Otherwise specify level of group
    wide.or.long = "wide") # output as wide (with columns for taxa) or long?
  
