
############################################################################
# Model specification and run code for estimating species specific capture 
# rates using closed Bayesian mark-recapture models from capture data 

# Covariates are species and number of traps ---------------------------- ##
# Using Michaelis-Menten equation with asymptote of 1 
# All species have an intercept of 0 and different half saturation------- ##
# constants that describe how number of traps affect daily capture rates  ##


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

# Michaelis-Menten equation:

# n_traps / (k[species] + n_traps)

################################################################################
## Specify the Robust Design monthly Closed Model---------------------------- ##
################################################################################


sink("RDp_species_MMtrap.bug")

cat("
      model{
       
      
  
      ##### PRIORS #####
      for(s in 1:n.group){
        k.trap[s] ~ dunif(0, 5)    #  slope of how number of traps affect capture rate per group
      }    
      
      

      
        ##### LIKELIHOOD #####
        
        ##### OBSERVATION PROCESS
        for(obs in 1:n.obs) {
      
          #### CONSTRAINTS ####
          p.eff[obs] =  traps_available[obs] / (k.trap[group[obs]] + traps_available[obs])
          y[obs] ~ dbern(p.eff[obs])

          
          } # obs
        
        ##### Monitoring derived parameters for 1 species ####
        for(tr in 1:n.tr){
         p.monitor[tr] =  trap_range[tr] / (k.trap[sp.monitor] + trap_range[tr])
         
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
    k.trap = runif(n.group, 0, 5)
  )
}

# parameters to monitor
parameters <- c("k.trap", "p.monitor")

# run the model
date()
NEONAllsites.speciestrapMM.p.mod <- jags.parallel(data = bugs.data,
                                   inits, 
                                   parameters, 
                                   "RDp_species_MMtrap.bug", 
                                   n.chains = 3, 
                                   n.thin   = 5, 
                                   n.iter   = 10000, 
                                   n.burnin = 2000)
date() #  10 minutes with parallel

save(NEONAllsites.speciestrapMM.p.mod, file="SmammalSpeciesTrapMMPestModOutput.RData")

# chains per 100 traps
p.speciestrap.MM.chains <- 1/ (NEONAllsites.speciestrapMM.p.mod$BUGSoutput$sims.list$k.trap +1)
p.speciestrapMM.estimates.per100traps <- data.frame(species=taxon.cap.sum$species, 
      p.est = apply(p.speciestrap.MM.chains, 2, mean),
      lowerCI = apply(p.speciestrap.MM.chains, 2, quantile, probs = 0.025),
      upperCI = apply(p.speciestrap.MM.chains, 2, quantile, probs = 0.975)) 
summary(p.speciestrapMM.estimates.per100traps)

p.speciestrapMM.model.estimates <- data.frame(
  species = taxon.cap.sum$species,
  trap_k = apply(NEONAllsites.speciestrapMM.p.mod$BUGSoutput$sims.list$k.trap, 2, mean),
  trap_kLCI = apply(NEONAllsites.speciestrapMM.p.mod$BUGSoutput$sims.list$k.trap, 2, quantile, probs=0.025),
  trap_kUCI = apply(NEONAllsites.speciestrapMM.p.mod$BUGSoutput$sims.list$k.trap, 2, quantile, probs=0.925))
                                            

# Make some plots -----------------------------------------------------------##

library(mcmcplots)
library(MCMCvis)

mcmcplot(NEONAllsites.speciestrapMM.p.mod,
         parms = parameters)

MCMCplot(NEONAllsites.speciestrapMM.p.mod,
         params = "k.trap")


pdf(file="SmammalSpeciesTrapCaptureRatesMM.pdf", width=8, height=20)
ggplot(data=p.speciestrapMM.estimates.per100traps, aes(x = species, y = p.est)) +
  geom_pointrange(aes(ymin = lowerCI, ymax = upperCI)) +
  ylim(0,1) +
  theme_minimal() +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  ggtitle("Species - Trap MM model: 
     Daily capture rates per 100 traps") 
dev.off()


## Compare CI estimates using monitored ------------------------------------ ##

# monitored # "Chaetodipus californicus"

# how does using the mean, upper and lower 95%CI compare to using the full
# posterior and then calculating the CIs?

# First - precalculated CIs:
CHCA.mean.p.fromsum_MM = trap_range / (p.speciestrapMM.model.estimates$trap_k[CHCA.num] + trap_range)
CHCA.lowerCI.p.fromsum_MM = trap_range / (p.speciestrapMM.model.estimates$trap_kLCI[CHCA.num] + trap_range)
CHCA.upperCI.p.fromsum_MM = trap_range / (p.speciestrapMM.model.estimates$trap_kUCI[CHCA.num] + trap_range)

plot((trap_range)*100, CHCA.mean.p.fromsum_MM, pch=19, ylim=c(0,1), 
     xlab="number of traps per night", ylab="daily capture probability",
     main="Chaetodipus californicus")
lines((trap_range)*100, CHCA.lowerCI.p.fromsum_MM, lty=2)
lines((trap_range)*100, CHCA.upperCI.p.fromsum_MM, lty=2)

# full posterior then get CI
d <- dim(NEONAllsites.speciestrapMM.p.mod$BUGSoutput$sims.list$k.trap)[1]
rangemat <- matrix(trap_range, nrow=d,  ncol=length(trap_range), byrow = T)
post_MM <- rangemat / (NEONAllsites.speciestrapMM.p.mod$BUGSoutput$sims.list$k.trap[,CHCA.num] + rangemat)
points((trap_range)*100, apply(post, 2, quantile, probs=0.5), pch=19, ylim=c(0,1), col="red")
lines((trap_range)*100, apply(post, 2, quantile, probs=0.025), lty=2, col="red")
lines((trap_range)*100, apply(post, 2, quantile, probs=0.975), lty=2, col="red")

# CIs a little wide on the high end when using full posterior.

### compare to linear model with species-slope:
plot((trap_range)*100, CHCA.mean.p.fromsum_MM, type="l",lwd=3, ylim=c(0,1), 
     xlab="number of traps per night", ylab="daily capture probability",
     main="Chaetodipus californicus")
lines((trap_range)*100, CHCA.lowerCI.p.fromsum_MM, lty=2)
lines((trap_range)*100, CHCA.upperCI.p.fromsum_MM, lty=2)
lines((trap_range)*100, CHCA.mean.p.fromsum_slope, lwd=2, col="blue")
lines((trap_range)*100, CHCA.lowerCI.p.fromsum_slope, lty=2, col="blue")
lines((trap_range)*100, CHCA.upperCI.p.fromsum_slope, lty=2, col="blue")
legend("topleft", c("Michaelis-Menten model", "linear species-slope model"),
       col=c("black", "blue"), lwd=3, bty="n")

NEONAllsites.speciestrapMM.p.mod$BUGSoutput$DIC
# 293693.8
loo::waic(-NEONAllsites.speciestrapMM.p.mod$BUGSoutput$sims.list$deviance)
# 587549.6

# by both DIC and WAIC, MM is better than species slope model

### Also look at deer mice ------------------------------------------------- ##
PEMA.num <- 85

PEMA.mean.p.fromsum_MM = trap_range / (p.speciestrapMM.model.estimates$trap_k[PEMA.num] + trap_range)
PEMA.lowerCI.p.fromsum_MM = trap_range / (p.speciestrapMM.model.estimates$trap_kLCI[PEMA.num] + trap_range)
PEMA.upperCI.p.fromsum_MM = trap_range / (p.speciestrapMM.model.estimates$trap_kUCI[PEMA.num] + trap_range)

plot((trap_range)*100, PEMA.mean.p.fromsum_MM, pch=19, ylim=c(0,1), 
     xlab="number of traps per night", ylab="daily capture probability",
     main="Peromyscus maniculatus")
lines((trap_range)*100, PEMA.lowerCI.p.fromsum_MM, lty=2)
lines((trap_range)*100, PEMA.upperCI.p.fromsum_MM, lty=2)

###############################################################################
## Use estimates to adjust captures for abundance estimates
###############################################################################

# should do this for all captures not just the reduced ones 
# ie. also include those only trapped 1 night per session

smammal.species.num.captures.longdata <- taxa.capture.fun(
  captures = smammal.captures.clean, 
  session.info = trapping.session.info,
  group = "species", # "species" or "genus" or other group identified in captures
  ID.col = "plot_taxon_tag", # individual identifier as character vector, e.g., "tagID"
  taxa.to.include = NULL, #  if NULL, then use all. Otherwise specify level of group
  wide.or.long = "long") # output as wide (with columns for taxa) or long?


smammal.species.estimates.longdata <- NEON.N.estimates.fromcaptures.MMmodel(
    num.captures = smammal.species.num.captures.longdata, #output from taxa.capture.fun function in long format
    p.param.estimates = p.speciestrapMM.model.estimates, # summary of Bayesian estimates as data.frame
    group = "species", # same label as column in p.param.estimates
    wide.or.long = "long")
# 5 rows with NAs, which are species that were not in the reduced data used to estimate. 
# "Tamiasciurus douglasii" "Sylvilagus nuttallii"   "Sorex nanus"            "Neotoma cinerea"        "Sylvilagus nuttallii"  
# each of these were caught once in the whole dataset. So just use captures for that instead of N. 
smammal.species.estimates.longdata$N.est[which(is.na(smammal.species.estimates.longdata$N.est))] <-
  smammal.species.estimates.longdata$n_caps[which(is.na(smammal.species.estimates.longdata$N.est))]

# for sessions with very low trapnights, it's going to overestimate N 
# because prob of capture is so low
# (with unrealistically skinny CI based on the constraints of the model)

dat <- trapping.session.info %>%
  group_by(plotID, prim.session) %>%
  summarise(trapnights = sum(traps_available))
summary(dat)
hist(dat$trapnights, breaks=30)
hist(dat$trapnights[which(dat$trapnights<100)], breaks=30)
quantile(dat$trapnights, probs= c(0.01, 0.1, 0.25, 0.5))

# remove sessions with fewer than 51 trapnights [1% of the plot-sessions] ----##

smammal.species.estimates.longdata <- smammal.species.estimates.longdata %>%
  filter(trapnights > 50)
summary(smammal.species.estimates.longdata)

ggplot(smammal.species.estimates.longdata, aes(x=n_caps, y=N.est, col=trapnights)) +
  geom_point() +
  geom_line(aes(x=n_caps, y=n_caps), color = "gray") +
  ggtitle("Comparison of captures and N estimates")
 

       

save(smammal.species.estimates.longdata, )