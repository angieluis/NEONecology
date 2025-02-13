###############################################################################
### Means of datasets over site - habitat - year within the 4 month window:
###############################################################################

# describe how I calculated the 4-month window
# and why not just matching by plot (don't match up)

setwd("~/Documents/NEON data/NEONecology")
load("ReducedMergedPlantSoilData.RData")
load("SmammalSpeciesEstimates.RData")

library(tidyverse)
library(kableExtra)
library(neonDivData)
library(sf)
library(scales)
library(tigris)
library(data.table)
library(iNEXT)

source("~/Documents/R Files/Numerical Ecology with R-2ed/NEwR2-Functions/panelutils.R")


###############################################################################
# Plant Species Richness per site-habitat-year -------------------------------#
# Using rarefaction based on incidence in the iNEXT package (see Hsieh et al. 2016 appendix)
# using all subplot-dates within the 4-month peak window
###############################################################################

# Do rarefaction using the 4 subplots per plot, each time sampled within 
# 4-month peak window for the year, as the sampling unit

# need data with species as rows and site, habitat and year as columns, with the first row
# containing the number of subplots. Each entry gives the number of subplots
# where that species was present (so no more than the first row)

# I created a separate column for species that truncates to first 2 words

plant.inc_bysubplot <- reduced.plant.cover %>%
  filter(peak.window.4mo == "Y", taxon_rank == "species" | 
           taxon_rank == "subspecies" | 
           taxon_rank == "variety") %>% ## only using those IDed to species
  select(siteID, plotID, nlcdClass, subplot_id, observation_datetime, species, presence_absence) %>%
  group_by(siteID, nlcdClass, plotID, subplot_id, year=year(observation_datetime), observation_datetime) %>%
  pivot_wider(names_from = species,
             values_from = presence_absence,
             values_fn = ~ ifelse(sum(.x) > 0, 1, 0),
             values_fill = 0,
             names_sort = TRUE) 

cols <- names(plant.inc_bysubplot)[7:dim(plant.inc_bysubplot)[2]]
plant.inc <- plant.inc_bysubplot %>% 
  group_by(siteID, nlcdClass, year) %>%
  summarise(n_subplots = n(),
            across(all_of(cols), sum))

plant_inc_t <- t(plant.inc[, 4:dim(plant.inc)[2]])

plantRarefaction_sitehabitatyear <- iNEXT(plant_inc_t, q=0, datatype = "incidence_freq") #q = 0 for species richness

ggiNEXT(plantRarefaction_sitehabitatyear, type=1) +
  theme_bw(base_size = 18) + theme(legend.position="none")

plantS_sitehabitatyear <- estimateD(plant_inc_t, q=0, datatype = "incidence_freq", base="size")
 # default with base="size" and level=NULL, computes the diversity estimates for the minimum among all the coverage values
plantS_sitehabitatyear <- data.frame(plant.inc[,1:4], plantS_sitehabitatyear)
# it rarefied to 8 subplots (I think that's what t is)

# compare to observed species richness (total number of unique species per site-habitat-year)
obsS <- reduced.plant.cover %>%
  filter(peak.window.4mo == "Y", taxon_rank == "species") %>% 
  group_by(siteID, nlcdClass, year = year(observation_datetime)) %>%
  summarise(plant.richness = length(unique(taxon_name))) 
plot(obsS$plant.richness,plantS_sitehabitatyear$qD, xlab="observed richness",
     ylab="Richness rarefied to 8 subplots")




# could also do the abundance version instead of incidence if want shannon or simpson diversity


###############################################################################
# Plant Species Percent Cover per site-habitat-year --------------------------#
# 
###############################################################################

# plant cover is assessed at typically 6 sub-subplots (1m2) per plot
# I will average % cover for each species over all subsubplots within the 
# 4 month window, per site-habitat-year


# For the microbial analyses, Ylva and Dean talked about taking one sample 
# closest to peak date or closest to sampling event,
# but there can be several plots sampled over several days (could maybe use bout?)
# but I need to talk with them some more about that


plantcover.nplotwindows <- reduced.plant.cover %>%
  group_by(siteID, plotID, nlcdClass, year = year(observation_datetime), 
           observation_datetime) %>%
  summarise(n_plotdates_4mo = length(unique(observation_datetime[which(peak.window.4mo=="Y")])),
            n_plotdates_3mo = length(unique(observation_datetime[which(peak.window.3mo=="Y")])),
            n_plotdates_2mo = length(unique(observation_datetime[which(peak.window.2mo=="Y")]))) %>%
  group_by(siteID,nlcdClass,year) %>%
  summarise(n_plotdates4mo = sum(n_plotdates_4mo),
            n_plots4m = length(unique(plotID[which(n_plotdates_4mo==1)])), 
            n_plotdates3mo = sum(n_plotdates_3mo),
            n_plots3m = length(unique(plotID[which(n_plotdates_3mo==1)])),
            n_plotdates2mo = sum(n_plotdates_2mo),
            n_plots2m = length(unique(plotID[which(n_plotdates_2mo==1)])))
            
 
x <- reduced.plant.cover %>%
  group_by(plotID, year = year(observation_datetime)) %>%
  summarise(n = length(unique(observation_datetime)))
summary(factor(x$n))
#    1    2    3    4 
# 2883  610   42    3 
# most plots have one date per year but some have up to 4.
# for those that are more than one are they considered different bouts?


plant.cover.subsubplotsitehabitat.date <- reduced.plant.cover %>%
  filter(peak.window.4mo == "Y", 
         sample_area_m2 == 1,            # only use the subsubplots since that's where they quantified percent cover, rest is presence/absence
         taxon_rank == "species" | 
           taxon_rank == "subspecies" | 
           taxon_rank == "variety") %>% ## only using those IDed to species
  select(siteID, plotID, subplotID, nlcdClass, observation_datetime, taxon_rank, species, percent_cover) %>% 
  ### need to first sum percent cover over each species because some subspecies, varieties were separated
  group_by(siteID, plotID, subplotID, nlcdClass, observation_datetime, species) %>% 
  summarise(percent_cover = sum(percent_cover)) %>%
  ## now pivot wider, and put in percent cover per subplot date if present, if not put in a 0 
  pivot_wider(names_from = species,
              values_from = percent_cover,
              values_fill = 0,
              names_sort = TRUE) 


cols <- names(plant.cover.subsubplotsitehabitat.date)[6:dim(plant.cover.subsubplotsitehabitat.date)[2]]

plant.cover.sitehabitatyear <- plant.cover.subsubplotsitehabitat.date %>% 
  group_by(siteID, nlcdClass, year=year(observation_datetime)) %>%
  summarise(across(all_of(cols), mean))

### mean per site
cols <- names(plant.cover.sitehabitatyear)[4:dim(plant.cover.sitehabitatyear)[2]]
plant.cover.site <- plant.cover.sitehabitatyear %>%
  group_by(siteID) %>%
  summarise(across(all_of(cols), mean))

save(plant.cover.site, file="MeanPlantPercentCoverPerSite.RData")

###############################################################################
# Plant Biomass per site-habitat-year
# taking mean of all samples within the 4 month peak season
###############################################################################

# Talk to others, may want to take one sample closest to peak date


plant.biomass.sitehabitatyear <- productivity.persample %>%
  filter(peak.window.4mo == "Y" ) %>%
  group_by(siteID, nlcdClass, year) %>%
  summarise(plant.biomass = mean(sum.dryMass),
            uCI.biomass = plant.biomass + 1.96*sd(sum.dryMass)/sqrt(length((sum.dryMass))),
            lCI.biomass = plant.biomass - 1.96*sd(sum.dryMass)/sqrt(length((sum.dryMass))))






###############################################################################
# Small Mammal Data ###########################################################
# from raw data, estimated species specific probability of capture
# using closed Bayesian mark-recapture models 
# Covariates were species and number of traps --------------------------- ##
# Using Michaelis-Menten equation with asymptote of 1 
# All species have an intercept of 0 and different half saturation------- ##
# constants that describe how number of traps affect daily capture rates  ##
###############################################################################

# population estimates per primary capture occasion (usually over 3 days)
smammal.species.estimates.longdata

# remove those not IDed to species
sm.species <- sort(unique(smammal.species.estimates.longdata$species))
rm1 <- grep("/",sm.species)
rm2 <- grep("sp[.]", sm.species)
sm.species <- sm.species[-c(rm1,rm2)]

## add productivity window information
# first add habitat column
plot.info <- neonDivData::neon_location
sm.plot.info <- plot.info %>%
  filter(str_detect(location_id,"mammalGrid"))
smammal.species.estimates.longdata <- left_join(smammal.species.estimates.longdata,
                                                sm.plot.info[,c("siteID","plotID","nlcdClass")] )

# add columns for different productivity windows
smammal.species.estimates.longdata <- left_join(left_join(smammal.species.estimates.longdata %>%
                                   mutate(year = year(first.date))
                                 , site.habitat.peakDate),
                       site.peakDate) %>%
  mutate(peak.date = mdy(paste(if_else(is.na(peakdate), site.peakdate, peakdate), year, sep="-")), #if it's present, use the site-habitat peak date, otherwise use the site peak date
         # is this within 1 month on either side of peak date (so 2 month window?)
         peak.window.2mo = factor(if_else(abs(difftime(ymd(first.date), 
                                                       peak.date, units="days"))<32 , "Y", "N")),
         # is this within 1.5 month on either side of peak date (so 3 month window?)
         peak.window.3mo = factor(if_else(abs(difftime(ymd(first.date), 
                                                       peak.date, units="days"))<46 , "Y", "N")),
         peak.window.4mo = factor(if_else(abs(difftime(ymd(first.date), 
                                                       peak.date, units="days"))<62 , "Y", "N")))
# sites which are only wetlands or ag don't have productivity data, so NA



smammal.plot.date <- smammal.species.estimates.longdata %>%
  filter(peak.window.4mo == "Y", 
         species %in% sm.species) %>% ## only using those IDed to species
  select(siteID, plotID, nlcdClass, first.date, species, N.est) %>%
  group_by(siteID, plotID, nlcdClass, first.date, species) %>% 
  pivot_wider(names_from = species,
              values_from = N.est,
              values_fill = 0,
              names_sort = TRUE) 


cols <- names(smammal.plot.date)[6:dim(smammal.plot.date)[2]]

smammal.sitehabitatyear <- smammal.plot.date %>% 
  group_by(siteID, nlcdClass, year=year(first.date)) %>%
  summarise(across(all_of(cols), mean))








###############################################################################
# Beetles (from DivData)

###############################################################################

beetlesDD <- neonDivData::data_beetle

names(beetlesDD)[which(names(beetlesDD)=='value')] <- "abundance"    #"count per trap day"

beetlesDD <- beetlesDD %>% 
  mutate_at(vars(taxon_id, taxon_name, taxon_rank, nativeStatusCode,
                 nlcdClass), factor)

# create new column called species 
beetlesDD$species <- rep(NA, dim(beetlesDD)[1])
beetlesDD$species[beetlesDD$taxon_rank=="species"|beetlesDD$taxon_rank=="subspecies"] <- 
  paste(str_split_i(beetlesDD$taxon_name[beetlesDD$taxon_rank=="species"|beetlesDD$taxon_rank=="subspecies"]," ",1), 
          str_split_i(beetlesDD$taxon_name[beetlesDD$taxon_rank=="species"|beetlesDD$taxon_rank=="subspecies"]," ",2))


# add columns for different productivity windows
beetlesDD <- left_join(left_join(beetlesDD %>%
                                             mutate(year = year(observation_datetime))
                                           , site.habitat.peakDate),
                                 site.peakDate) %>%
  mutate(peak.date = mdy(paste(if_else(is.na(peakdate), site.peakdate, peakdate), year, sep="-")), #if it's present, use the site-habitat peak date, otherwise use the site peak date
         # is this within 1 month on either side of peak date (so 2 month window?)
         peak.window.2mo = factor(if_else(abs(difftime(ymd(observation_datetime), 
                                                       peak.date, units="days"))<32 , "Y", "N")),
         # is this within 1.5 month on either side of peak date (so 3 month window?)
         peak.window.3mo = factor(if_else(abs(difftime(ymd(observation_datetime), 
                                                       peak.date, units="days"))<46 , "Y", "N")),
         peak.window.4mo = factor(if_else(abs(difftime(ymd(observation_datetime), 
                                                       peak.date, units="days"))<62 , "Y", "N")))
# sites which are only wetlands or ag don't have productivity data, so NA


# site-habitat-year average abundances (count per trap day) during 4-month window
beetle.plot.date <- beetlesDD %>%
  filter(peak.window.4mo == "Y", !is.na(species)) %>% ## only using those IDed to species & within window
  select(siteID, plotID, nlcdClass, trapID, observation_datetime, species, abundance) %>%
  group_by(siteID, plotID, nlcdClass, observation_datetime, species) %>% # do I want to use bout instead of day?
  pivot_wider(names_from = species,
              values_from = abundance, 
              values_fn = sum, # over trapIDs per plot
              values_fill = 0,
              names_sort = TRUE) 

cols <- names(beetle.plot.date)[6:dim(beetle.plot.date)[2]]

beetles.sitehabitatyear <- beetle.plot.date %>% 
  group_by(siteID, nlcdClass, year=year(observation_datetime)) %>%
  summarise(across(all_of(cols), mean))



### Rarefaction by plot (incidence) -------------------------------------------#
### this doesn't look good - re-evaluate - by count?


# beetles.inc_byplot <- beetle.plot.date #replace ifelse(abundance>0,1,0)

# cols <- names(beetles.inc_byplot)[29:dim(beetles.inc_byplot)[2]]
# beetles.inc <- beetles.inc_byplot %>% 
#   group_by(siteID, nlcdClass, year) %>%
#   summarise(n_subplots = n(),
#             across(all_of(cols), sum))
# 
# beetles_inc_t <- t(beetles.inc[, 4:dim(beetles.inc)[2]])
# 
# #beetleRarefaction_sitehabitatyear <- iNEXT(beetles_inc_t, q=0, datatype = "incidence_freq") #q = 0 for species richness
# 
# #ggiNEXT(beetleRarefaction_sitehabitatyear, type=1) +
# #  theme_bw(base_size = 18) + theme(legend.position="none")
# 
# beetleS_sitehabitatyear <- estimateD(beetles_inc_t, q=0, datatype = "incidence_freq", base="size")
# # default with base="size" and level=NULL, computes the diversity estimates for the minimum among all the coverage values
# beetleS_sitehabitatyear <- data.frame(beetles.inc[,1:4], beetleS_sitehabitatyear)
# # several sites have only one species. Estimation is not robust for those.
# 
# # compare to observed species richness (total number of unique species per site-habitat-year)
# obsS <- beetlesDD %>%
#   filter(peak.window.4mo == "Y", !is.na(species)) %>% 
#   group_by(siteID, nlcdClass, year = year(observation_datetime)) %>%
#   summarise(beetle.richness = length(unique(species))) 
# plot(obsS$beetle.richness,beetleS_sitehabitatyear$qD, xlab="observed richness",
#      ylab="Richness rarefied to 2 plots")
# ### this is pretty crap. Max richness is 2. max observed richness >40. 
# # trapping days is usually 14 but not always. think about that.

#------------------------------------------------------------------------------#

  

############################################################################
## below is from before - need to update




site.habitat.year.means.4mo <- full_join(site.habitat.year.means.4mo,
                                     beetle.richness.site.habitat.4mo)

bird.richness.site.habitat.4mo <- birds %>%
  filter(peak.window.4mo == "Y") %>%
  group_by(siteID, nlcdClass, plotID, year=year(observation_datetime), observation_datetime) %>%
  summarise(bird.richness = length(unique(taxon_id[which(taxon_rank=="species")]))) %>%
  group_by(siteID,nlcdClass, year) %>%
  summarise(birdS = mean(bird.richness))

site.habitat.year.means.4mo <- full_join(site.habitat.year.means.4mo,
                                     bird.richness.site.habitat.4mo)


smammal.richness.site.habitat.4mo <- smammals %>%
  filter(peak.window.4mo == "Y") %>%
  group_by(siteID, nlcdClass, plotID, year=year(observation_datetime), observation_datetime) %>%
  summarise(smammal.richness = length(unique(taxon_id[which(taxon_rank=="species")]))) %>%
  group_by(siteID,nlcdClass, year) %>%
  summarise(smammalS = mean(smammal.richness))

site.habitat.year.means.4mo <- full_join(site.habitat.year.means.4mo,
                                         smammal.richness.site.habitat.4mo)

  
soil.p.samples.sitehabitat4 <- soil.periodic.merge %>%
  filter(peak.window.4mo == "Y" ) %>%
  group_by(siteID, nlcdClass, year = year(collectDate)) %>%
  summarise(Npercent = mean(nitrogenPercent, na.rm = T),
            Cpercent = mean(organicCPercent, na.rm = T),
            CNratio = mean(CNratio, na.rm = T),
            soilMoisture = mean(soilMoisture, na.rm=T),
            soilInWaterpH = mean(soilInWaterpH, na.rm=T))

site.habitat.year.means.4mo <- full_join(site.habitat.year.means.4mo,
                                     soil.p.samples.sitehabitat4)

  
## means across years for each site-habitat -------------------------------- #

cs <- names(site.habitat.year.means.4mo)[4:dim(site.habitat.year.means.4mo)[2]]

site.habitat.means.4mo <- site.habitat.year.means.4mo %>%
  group_by(siteID, nlcdClass) %>%
  summarise(across(all_of(cs), ~ mean(.x, na.rm=TRUE)))

## add in the initial soil data to this... 

# probably not appropriate- because averaging across horizons
soil.initial.sitehabitat <- soil.initial %>%
  group_by(siteID, nlcdClass, plotID) %>%
  summarise(ctonRatio = mean(ctonRatio, na.rm=T),
            MehlichIIITotP = mean(MehlichIIITotP, na.rm = T),
            bulkDensCenterDepth = mean(bulkDensCenterDepth, na.rm = T),
            sandTotal = mean(sandTotal, na.rm = T),
            siltTotal = mean(siltTotal, na.rm = T),
            clayTotal = mean(clayTotal, na.rm = T)) %>%
  group_by(siteID, nlcdClass) %>%
  summarise(iCNratio = mean(ctonRatio, na.rm=T),
            iMehlichP = mean(MehlichIIITotP, na.rm=T),
            ibulkDensCenter = mean(bulkDensCenterDepth, na.rm=T),
            iSand = mean(sandTotal, na.rm=T),
            iSilt = mean(siltTotal, na.rm=T),
            iClay = mean(clayTotal, na.rm=T)
            )
  
site.habitat.means.4mo <- full_join(site.habitat.means.4mo,
                                    soil.initial.sitehabitat) 



pairs(
  drop_na(site.habitat.means.4mo[,-c(1:2)]),
  lower.panel = panel.smooth,
  upper.panel = panel.cor,
  diag.panel = panel.hist,
  main = "site-habitat correlations"
)

write_csv(site.habitat.means.4mo, file="SiteHabitatMeans.csv")



#############################################################################
## Bring in smammal community
## see "ModelRunSpeciesTrapMMCaptureRateEstimations.R" for how estimate N
#############################################################################


load("SmammalSpeciesEstimates.RData")
smammal.species.estimates.longdata
# add in habitat type for each plot
smammal.species.estimates.longdata <- left_join(smammal.species.estimates.longdata ,
                                                smammal.plot.data[,1:7])

#choose 'growing season' equivalent to plants - 4 month window
# add a column to this for the window like the other datasets.
smammal.species.estimates.longdata$first.date <- ymd(smammal.species.estimates.longdata$first.date)
smammal.species.estimates.longdata <- left_join(left_join(smammal.species.estimates.longdata %>%
                                             mutate(year = year(first.date))
                                           , site.habitat.peakDate),
                                 site.peakDate) %>%
  mutate(peak.date = mdy(paste(if_else(is.na(peakdate), site.peakdate, peakdate), year, sep="-")), #if it's present, use the site-habitat peak date, otherwise use the site peak date
         # is this within 1 month on either side of peak date (so 2 month window?)
         peak.window.2mo = factor(if_else(abs(difftime(ymd(first.date), 
                                                       peak.date, units="days"))<32 , "Y", "N")),
         # is this within 1.5 month on either side of peak date (so 3 month window?)
         peak.window.3mo = factor(if_else(abs(difftime(ymd(first.date), 
                                                       peak.date, units="days"))<46 , "Y", "N")),
         peak.window.4mo = factor(if_else(abs(difftime(ymd(first.date), 
                                                       peak.date, units="days"))<62 , "Y", "N")))

#NA if before 2019.



#### get mean of each species abundances within the 4-month window for each site-habitat 
smammal.site.habitat.means <- smammal.species.estimates.longdata %>%
  filter(peak.window.4mo == "Y") %>%
  group_by(siteID, nlcdClass, species ) %>%
  summarise(meanN = mean(N.est, na.rm=T)) %>%
  pivot_wider(names_from = species,
              values_from = meanN,
              values_fill = 0)

# remove ones that aren't species
sp <- sort(colnames(smammal.site.habitat.means[,-c(1:2)]))
rm1 <- grep("/", sp)
rm2 <- grep("sp[.]",sp)
sp <- sp [-c(rm1,rm2)]
match(sp, colnames(smammal.site.habitat.means))

smammal.site.habitat.means <- smammal.site.habitat.means[,c(1,2,
                                                         match(sp, colnames(smammal.site.habitat.means)))]

sort(colSums(smammal.site.habitat.means[,3:131]))
# 17 species with sum <2
# 43 species with sum <5
# 62 species <10
# need to look at correlations within site among habitats - not likely independent, probably need site to be the replicate

smammal.site.means <- smammal.species.estimates.longdata %>%
  filter(peak.window.4mo == "Y") %>%
  group_by(siteID, species ) %>%
  summarise(meanN = mean(N.est, na.rm=T)) %>%
  pivot_wider(names_from = species,
              values_from = meanN,
              values_fill = 0)
smammal.site.means <- smammal.site.means[,
                              c(1,match(sp, colnames(smammal.site.means)))]
sort(colSums(smammal.site.means[,2:130]))
# 18 <2
# 78 <10

smammal.site.data <- smammal.plot.data %>%
  group_by(siteID) %>%
  summarise(latitude = mean(latitude),
            longitude = mean(longitude),
            elevation = mean(elevation),
            sum.dates = sum(n.dates),
            sum.trapnights = sum(trapnights),
            nlcdClasses = str_c(unique(nlcdClass), collapse = " | "))


plant.site.means  <- full_join(reduced.plant.cover %>%
    filter(peak.window.4mo == "Y") %>%
    group_by(siteID, nlcdClass, plotID, year=year(observation_datetime), observation_datetime) %>%
    summarise(plant.richness = length(unique(species))) %>%
    group_by(siteID) %>%
    summarise(plantS = mean(plant.richness))
,
  plant.biomass.sitehabitat.4mo <- productivity.persample %>%
    filter(peak.window.4mo == "Y" ) %>%
    group_by(siteID) %>%
    summarise(plant.biomass = mean(sum.dryMass))
)
# MOAB is missing from the plant cover data - was never sampled during peak window - which doesn't make sense

save(plant.site.means, smammal.site.data, smammal.site.habitat.means, smammal.site.means, file="SmammalMeans.RData")

