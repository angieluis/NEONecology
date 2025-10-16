############################################################################
## Format NEON data from all sites all species so I can run Bayesian
## Closed Robust design models to estimate a species-specific capture
## rate, p. I will use this to estimate abundance for each species
## for each month for each site.
##
## see "FunctionsForAllSpeciesEstimations.R" for data functions
## See "ModelRunAllSpeciesEstimation.R" for model and run code
############################################################################

library(tidyverse)

setwd("~/Documents/NEON data/NEON_count-small-mammals/NEONsmammals/")

# load capture data
load("reducedSWNEONsmammalcaptures.RData")

# source in functions
source("~/Documents/JAGS/BayesianMarkRecapSNV/RealData/Code/FunctionsForAllSpeciesEstimations.R")

reducedSWplotsummary <- reducedSWcaptures.clean %>%
  group_by(plotID) %>%
  summarise(n.months = length(unique(paste(year(ymd(date)),month(ymd(date)), sep="-"))), 
            n.species = length(unique(taxonID)))

# remove some species  like "leam" # Lepus americanus, "othe" # other mammal, "rosp" # rodentia
sp.remove <- read.csv("NEONSpecies_remove.csv")
rm.sp <- sp.remove$taxonID[which(sp.remove$remove.other==1)]

reducedSWcaptures.clean <- reducedSWcaptures.clean %>%
  filter(!taxonID %in% rm.sp)

# do the whole thing twice. first time removing all those not id to species
# second time group together all to genus (so incl all the 'sp' ones along with all the others)

reducedSWcaptures.clean$genus <- str_split_i(reducedSWcaptures.clean$scientificName, " ", 1)
# genera <- sort(unique(reducedSWcaptures.clean$genus))


##############################################################################
## Use functions to format data for all years
##############################################################################


NEON.CHlist.all <- NEON.pRDcapture.history.list.fun( 
  cleaned.data = reducedSWcaptures.clean,  
  prim.session.list = reducedSW.prim.session.list,  
  plotID = reducedplots, 
  species = "all")  

length(NEON.CHlist.all$Ch.list) # 47
dim(NEON.CHlist.all$Ch.list[[1]]) # 12516     3
sum(unlist(lapply(NEON.CHlist.all$Ch.list, sum, na.rm=TRUE))) 
# total number of captures = 31497 reducedSWcaptures.clean has 31634 rows. these should be the same, 
# missing about 130 captures. In the cleaned data there are no duplicates by site_tag and date.
# missing plots or dates?

NEON.individual.covariates <- NEON.p.individual.covariate.fun(
  cleaned.data = reducedSWcaptures.clean,
  CH.secondary = NEON.CHlist.all$Ch.list,  
  tags = rownames(NEON.CHlist.all$Ch.list[[1]]),   
  plots = reducedplots
  
)


plot.sp.sum <- NEON.individual.covariates %>%
  group_by(plot, species) %>%
  summarise(nind=length(unique(site_tag)))
print(plot.sp.sum,n=125)

n.species <- length(unique(NEON.individual.covariates$species))
#levels(NEON.individual.covariates$species)
sp.sum <- NEON.individual.covariates %>%
  group_by(species) %>%
  summarise(nind=length(unique(site_tag)))
print(sp.sum, n=40)
# several of these species there is only 1 individual. 
# and those not id to species are still in here
# I will leave it in for estimating p but not use those to estimate N.


# turn species and genus into a number
species.info <- reducedSWcaptures.clean %>%
  group_by(taxonID, scientificName, genus) %>%
  summarise(n_caps = n())
species.info$species.num <- 1:dim(species.info)[1]
genera.df <- data.frame(genus = sort(unique(species.info$genus)), genus.num = 1:length(sort(unique(species.info$genus))))
species.info$genus.num <- match(species.info$genus, genera.df$genus)
write_excel_csv(species.info, file="NEONspecies.csv")

NEON.individual.covariates$species.num <- match(NEON.individual.covariates$species, species.info$taxonID)
NEON.individual.covariates$genus.num <- species.info$genus.num[match(NEON.individual.covariates$genus, species.info$genus)]

obs.dat <- NEON.monthly.longdata.CH.for.p.fun(
  CH.secondary = NEON.CHlist.all$Ch.list,
  individual.covariates = NEON.individual.covariates)
  

save(obs.dat, species.info, NEON.individual.covariates, NEON.CHlist.all, file="NEONDataforPmodels.RData")


# num.captures.all.species.longdata <- species.capture.fun(CH.secondary = CHlist.all$Ch.list,
#        temporal.covariates = covariate.data.all$temporal.covariates,
#        individual.covariates = NEON.individual.covariates,
#        n.sec.occ.list = CHlist.all$n.sec.occ
#      )
# head(num.captures.all.species.longdata)
# 
# save(num.captures.all.species.longdata, obs.dat, covariate.data.all, CHlist.all, file="NEONDataforPmodels.RData")


## See "ModelRunAllSpeciesEstimationsNEON.R" to Run the Model 
