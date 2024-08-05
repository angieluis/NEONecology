
############################################################################
# Model specification and run code for estimating species specific capture 
# rates using closed Bayesian mark-recapture models from capture data 

# Covariates are species and number of traps ---------------------------- ##
# All species have an intercept of 0 and separate slopes for how -------- ##
# number of traps affect daily capture rates ---------------------------- ##


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


# transforming traps by dividing by 100 below

################################################################################
## Specify the Robust Design monthly Closed Model---------------------------- ##
################################################################################


sink("RDp_species_slopetrap.bug")

cat("
      model{
       
      
  
      ##### PRIORS #####
      for(s in 1:n.group){
        p.trap[s] ~ dunif(0, 1.2)    #  slope of how number of traps affect capture rate per group
      }    
      
      

      
        ##### LIKELIHOOD #####
        
        ##### OBSERVATION PROCESS
        for(obs in 1:n.obs) {
      
          #### CONSTRAINTS ####
          p.eff = min(p.trap[group[obs]] * traps_available[obs], 1) # not transforming so make so can't be >1
      
          y[obs] ~ dbern(p.eff)

          
          } # obs
        
        ##### Monitoring derived parameters for 1 species ####
        for(tr in 1:n.tr){
         p.monitor[tr] =  min(p.trap[sp.monitor] * trap_range[tr], 1)
        }
        
        } # model
    ", fill = TRUE)

sink()


## -------------------------------------------------------------------------- ##


##############################################################################
## Run the Model for species + traps
##############################################################################



# packages
library(R2jags)

## To test output and CI, I'll monitor one species with 75 captures:
# Chaetodipus californicus
CHCA.num <- which(taxon.cap.sum$species=="Chaetodipus californicus")
# over this standardized trap range
trap_range <- seq(0,1.20, length=21)

## Bundle data -------------------------------------------------------------- ##

# list of bugs data
bugs.data <- list(
  n.obs                 = length(NEON.CHlong.all$State),
  y                     = NEON.CHlong.all$State, 
  n.group               = max(as.numeric(NEON.CHlong.all$species.num)), # group is species for this model
  group                 = as.numeric(NEON.CHlong.all$species.num),
  traps_available       = NEON.CHlong.all$traps_available/100,
  sp.monitor            = CHCA.num, # just for monitoring derived parameters of 1 species
  trap_range            = trap_range, # just for monitoring derived parameters
  n.tr                  = length(trap_range) # just for monitoring derived parameters
)


# supply initial values
inits <- function() {
  list(
    p.trap = runif(n.group, 0, 0.5)
  )
}

# parameters to monitor
parameters <- c("p.trap", "p.monitor")

# run the model
date()
NEONAllsites.speciestrapslope.p.mod <- jags.parallel(data = bugs.data,
                                   inits, 
                                   parameters, 
                                   "RDp_species_slopetrap.bug", 
                                   n.chains = 3, 
                                   n.thin   = 5, 
                                   n.iter   = 10000, 
                                   n.burnin = 2000)
date() #  10 minutes with parallel

save(NEONAllsites.speciestrapslope.p.mod, file="SmammalSpeciesTrapSlopePestModOutput.RData")


p.speciestrap.slope.chains <- NEONAllsites.speciestrapslope.p.mod$BUGSoutput$sims.list$p.trap # 100 traps is 1 when transformed so don't need to multiply
p.speciestrapslope.estimates.per100traps <- data.frame(species=taxon.cap.sum$species, 
      p.est = apply(p.speciestrap.slope.chains, 2, mean),
      lowerCI = apply(p.speciestrap.slope.chains, 2, quantile, probs = 0.025),
      upperCI = apply(p.speciestrap.slope.chains, 2, quantile, probs = 0.975)) 
p.speciestrapslope.estimates.per100traps$upperCI <- replace(p.speciestrapslope.estimates.per100traps$upperCI, p.speciestrapslope.estimates.per100traps$upperCI>1,1)
summary(p.speciestrapslope.estimates.per100traps)

p.speciestrapslope.model.estimates <- data.frame(
  species = taxon.cap.sum$species,
  # intercept = apply(NEONAllsites.speciestrap.p.mod$BUGSoutput$sims.list$p.0, 2, mean),
  # interceptLowerCI = apply(NEONAllsites.speciestrap.p.mod$BUGSoutput$sims.list$p.0, 2, quantile, probs = 0.025),
  # interceptUpperCI = apply(NEONAllsites.speciestrap.p.mod$BUGSoutput$sims.list$p.0, 2, quantile, probs = 0.975),
  trap_slope = apply(NEONAllsites.speciestrapslope.p.mod$BUGSoutput$sims.list$p.trap, 2, mean),
  trap_slopeLCI = apply(NEONAllsites.speciestrapslope.p.mod$BUGSoutput$sims.list$p.trap, 2, quantile, probs=0.025),
  trap_slopeUCI = apply(NEONAllsites.speciestrapslope.p.mod$BUGSoutput$sims.list$p.trap, 2, quantile, probs=0.925))
                                            

# Make some plots -----------------------------------------------------------##

library(mcmcplots)
library(MCMCvis)

mcmcplot(NEONAllsites.speciestrap.p.mod,
         parms = parameters)

pdf(file="SmammalSpeciesTrapCaptureRatesSlope.pdf", width=8, height=20)
ggplot(data=p.speciestrapslope.estimates.per100traps, aes(x = species, y = p.est)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI)) +
  ylim(0,1) +
  theme_minimal() +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  ggtitle("Species - Trap Slope model: 
     Daily capture rates per 100 traps") 
dev.off()

plot(p.speciestrapslope.estimates.per100traps$p.est, p.species.estimates$mean)
# compare the species model with the 

## Compare CI estimates using monitored ------------------------------------ ##

# monitored # "Chaetodipus californicus"

# how does using the mean, upper and lower 95%CI compare to using the full
# posterior and then calculating the CIs?


# First - precalculated CIs:
CHCA.mean.p.fromsum_slope = if_else(p.speciestrapslope.model.estimates$trap_slope[CHCA.num] * trap_range>1, 1,
                                    p.speciestrapslope.model.estimates$trap_slope[CHCA.num] * trap_range)
CHCA.lowerCI.p.fromsum_slope = if_else(p.speciestrapslope.model.estimates$trap_slopeLCI[CHCA.num] * trap_range>1, 1,
                                       p.speciestrapslope.model.estimates$trap_slopeLCI[CHCA.num] * trap_range)
CHCA.upperCI.p.fromsum_slope = if_else(p.speciestrapslope.model.estimates$trap_slopeUCI[CHCA.num] * trap_range>1, 1,
                                       p.speciestrapslope.model.estimates$trap_slopeUCI[CHCA.num] * trap_range)

plot((trap_range)*100, CHCA.mean.p.fromsum_slope, pch=19, ylim=c(0,1), 
     xlab="number of traps per night", ylab="daily capture probability",
     main="Chaetodipus californicus")
lines((trap_range)*100, CHCA.lowerCI.p.fromsum_slope, lty=2)
lines((trap_range)*100, CHCA.upperCI.p.fromsum_slope, lty=2)

# full posterior then get CI
d <- dim(NEONAllsites.speciestrapslope.p.mod$BUGSoutput$sims.list$p.trap)[1]
rangemat <- matrix(trap_range, nrow=d,  ncol=length(trap_range), byrow = T)
post <- NEONAllsites.speciestrapslope.p.mod$BUGSoutput$sims.list$p.trap[,CHCA.num] * rangemat # nothing >1 
points((trap_range)*100, apply(post, 2, quantile, probs=0.5), pch=19, ylim=c(0,1), col="red")
lines((trap_range)*100, apply(post, 2, quantile, probs=0.025), lty=2, col="red")
lines((trap_range)*100, apply(post, 2, quantile, probs=0.975), lty=2, col="red")

# CIs a little wider on the top end when using full posterior.
# with - traps set, still have 40% capture prob. That's not right. shoudl I force intercept to 0?
# if so, then would have to calculate different slope per species instead of intercept.
NEONAllsites.speciestrapslope.p.mod$BUGSoutput$DIC
# 296060.8

loo::waic(-NEONAllsites.speciestrapslope.p.mod$BUGSoutput$sims.list$deviance)
# 592315.9


### Also look at deer mice ------------------------------------------------- ##
PEMA.num <- 85
PEMA.mean.p.fromsum_slope = if_else(p.speciestrapslope.model.estimates$trap_slope[PEMA.num] * trap_range>1, 1,
                                    p.speciestrapslope.model.estimates$trap_slope[PEMA.num] * trap_range)
PEMA.lowerCI.p.fromsum_slope = if_else(p.speciestrapslope.model.estimates$trap_slopeLCI[PEMA.num] * trap_range>1, 1,
                                       p.speciestrapslope.model.estimates$trap_slopeLCI[PEMA.num] * trap_range)
PEMA.upperCI.p.fromsum_slope = if_else(p.speciestrapslope.model.estimates$trap_slopeUCI[PEMA.num] * trap_range>1, 1,
                                       p.speciestrapslope.model.estimates$trap_slopeUCI[PEMA.num] * trap_range)

plot((trap_range)*100, PEMA.mean.p.fromsum_slope, pch=19, ylim=c(0,1), 
     xlab="number of traps per night", ylab="daily capture probability",
     main="Peromyscus maniculatus")
lines((trap_range)*100, PEMA.lowerCI.p.fromsum_slope, lty=2)
lines((trap_range)*100, PEMA.upperCI.p.fromsum_slope, lty=2)

