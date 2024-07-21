
library(tidyverse)
library(neonUtilities)
library(neonDivData)
library(maps)
library(scales)




###############################################################################
## Look at diversity in plants using neonDivData
###############################################################################

# reloaded the neonDivData package to get updated data
# install.packages("neonDivData", repos = 'https://daijiang.r-universe.dev')

#neon_plants <- data_plant #before update
plantsDD <- neonDivData::data_plant

plant.sites <- sort(unique(plantsDD$siteID))
plant.plots <- sort(unique(plantsDD$plotID))
length(plant.plots)
#1577

print(plantsDD %>%
  group_by(siteID) %>%
  summarise(n.plots = length(unique(plotID)))
  ,n=length(plant.sites))


names(plantsDD)[which(names(plantsDD)=='value')] <- "percent_cover"


plantsDD <- plantsDD %>% 
  mutate_at(vars(taxon_id, taxon_name, taxon_rank, nativeStatusCode, heightPlantOver300cm,
                 nlcdClass), factor)

### some are IDed to subspecies
# [981] Calystegia silvatica (Kit.) Griseb.                                                               
# [982] Calystegia silvatica (Kit.) Griseb. ssp. fraterniflora (Mack. & Bush) Brummitt                    
# when IDed to subspecies the 4 letter code is longer than 4 letters, so if 
# truncate to first 4 letters, it removes subspecies, but it's a problem see below

# # add column that truncates taxonID to 4 letters, so removes subspecies
# plantsDD$taxonID4 <- factor(str_trunc(as.character(plantsDD$taxon_id), 4, 
#                                          "right", ellipsis = ""))
# # this doesn't work. see below


# add columms of genus names that truncates to first 1 word & species names
# truncates to first 2 words
# some are IDed to 'species group' so these might be misleading?
plantsDD$genus <- str_split_i(plantsDD$taxon_name," ",1)
plantsDD$species <- paste(plantsDD$genus , str_split_i(plantsDD$taxon_name," ",2))

# # tomake sure species lines up to taxonID4:
# dat <- plantsDD %>% 
#   group_by(taxonID4) %>%
#   summarise(n.sp=length(unique(species)))
# summary(dat)
# # nope
#
# dat2 <- plantsDD %>% 
#     group_by(species) %>%
#     summarise(n.sp=length(unique(taxonID4)))
# summary(dat2)
# 
# # each species has only 1 4-letter code, but some 4-letter codes have >1 species, up to 14.
# # look at see which
# dat3 <- plantsDD %>% 
#   group_by(taxonID4) %>%
#   summarise(n.sp=length(unique(species)), taxon=paste(unique(plantsDD$taxon_name)),
#             species=paste(unique(plantsDD$species)))
# 
# #AAAAAAAAGGH! 
# # "ABFR" and "ABFR2" are 2 different species
# # Abies fraseri"    "Abronia fragrans"
# # so truncating to first 4 does not make them unique species
# # hopefully using first 2 words of taxon_name will work for species.
# # remove taxonID4

# plantsDD <- plantsDD[,-(which(names(plantsDD)=="taxonID4"))]

species.codes <- plantsDD %>%
  group_by(species, genus, taxon_rank) %>%
  summarise(codes=str_flatten(unique(taxon_id),collapse=" | "))

write.csv(species.codes, file="PlantSpecies.csv", row.names = F)

# <--------------------------- still need to clean this up for diversity estimation
### Also some where are not ID to species
#  [983] Calystegia sp. 
# also some say "sp." and some say "spp."
# will probably wnat to make sure that something from this genus is already 
# counted in the diversity estimate or make it's own species
# I checked, and all those with 'sp.' and 'spp.', say 'genus' under taxon_rank. 

## <------------------------------------- Need to redo all this because 
# If subplotID ends in .1 then they measured percent cover at 1 square-meter plot. 
# If ends in .10 (or no decimal) then the 10 sq meter plot (or 100 sqm plot)
# then they just wrote down all the species they saw – so NA under percent cover
# I need to make sure those NAs aren't removed and are counted in the species 
# richness estimates

# this is getting species richness per plot per day - should I sum for each plot
# over any day? all unique species per plot over any time sampled? number of 
# sampling occasions and plots and times of year will affect it
plantdiversity_byplotdate <- plantsDD %>%
  group_by(siteID, plotID, observation_datetime) %>%
  summarise(subplots = length(unique(subplot_id)), 
            sp_richness = length(unique(taxon_name)), 
            any_invasives = any(nativeStatusCode=="I"),
            total_cover = sum(percent_cover, na.rm=TRUE))

# how should I account for different number of subplots?

# mean over all dates sampled for a plot
plantdiversity_byplot <- plantdiversity_byplotdate %>%
  group_by(siteID, plotID) %>%
  summarise(mean_sp_richness = mean(sp_richness), n_surveys = n(), 
            any_invasives = any(any_invasives))

# even within just Abby there is a lot of variation in mean species richness
# mean and SE of mean species richness for each plot within a site

plantdiversity_bysite <- plantdiversity_byplot %>%
  group_by(siteID) %>%
  summarise(mean_mean_sp_richness=mean(mean_sp_richness),
            SE_mean_sp_richness=sd(mean_sp_richness)/sqrt(n()))

# joined with total species richness - all taxa ever recorded for the site 
plantdiversity_bysite <- full_join(plantdiversity_bysite,
                                   plantsDD %>%
                                       group_by(siteID) %>%
                                       summarise(total_sp_richness=length(unique(taxon_name))),
                                     by="siteID")

summary(plantdiversity_bysite)
# mean and total pretty correlated

### right now, some are IDed to subspecies
# [981] Calystegia silvatica (Kit.) Griseb.                                                               
# [982] Calystegia silvatica (Kit.) Griseb. ssp. fraterniflora (Mack. & Bush) Brummitt                    
# when IDed to subspecies the 4 letter code is longer than 4 letters, so if 
# truncate to first 4 letters, it removes subspecies
### Also some where are not ID to species
#  [983] Calystegia sp. 
# will probably wnat to make sure that something from this genus is already 
# counted in the diversity estimate or make it's own species

# 6/24/2013 - 10/13/2021

###############################################################################
## Download raw data for plant cover and see what DivData did to clean it 
# remove bare ground, etc..
###############################################################################

# plants.raw <- loadByProduct("DP1.10058.001",
#                        token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
# )
                            
# endDate 6/24/2013 - 12/1/2022


###############################################################################
## Look at diversity in beetles using neonDivData
###############################################################################

beetles <- data_beetle

names(beetles)[which(names(beetles)=='value')] <- "abundance"    #"count per trap day"

beetles <- beetles %>% 
  mutate_at(vars(taxon_id, taxon_name, taxon_rank, nativeStatusCode,
                 nlcdClass), factor)

beetles$trappingDays <- as.numeric(beetles$trappingDays)

# again this is getting species richness per plot per day
beetlesdiversity_byplotdate <- beetles %>%
  group_by(siteID, plotID, observation_datetime) %>%
  summarise(sp_richness = length(unique(taxon_name)), 
            any_invasives = any(nativeStatusCode=="I"),
            total_abundance = sum(abundance,
                                  na.rm=TRUE))

# mean over all dates sampled for a plot
beetlesdiversity_byplot <- beetlesdiversity_byplotdate %>%
  group_by(siteID, plotID) %>%
  summarise(mean_sp_richness = mean(sp_richness), n_surveys = n(), 
            any_invasives = any(any_invasives),)


beetlesdiversity_bysite <- beetlesdiversity_byplot %>%
  group_by(siteID) %>%
  summarise(mean_mean_sp_richness=mean(mean_sp_richness),
            SE_mean_sp_richness=sd(mean_sp_richness)/sqrt(n()))

# joined with total species richness - all taxa ever recorded for the site 
beetlesdiversity_bysite <- full_join(beetlesdiversity_bysite,
  beetles %>%
    group_by(siteID) %>%
    summarise(total_sp_richness=length(unique(taxon_name))),
  by="siteID")

summary(beetlesdiversity_bysite)
plot(beetlesdiversity_bysite$mean_mean_sp_richness, beetlesdiversity_bysite$total_sp_richness)
# mean and total not very well correlated


#############################################################################
## Bird data
#############################################################################


birds <- data_bird

names(birds)[which(names(birds)=='value')] <- "abundance"    #"count of individuals"

birds <- birds %>% 
  mutate_at(vars(siteID, plotID, taxon_id, taxon_name, taxon_rank, targetTaxaPresent, detectionMethod,
                 visualConfirmation, sexOrAge, observedHabitat, nativeStatusCode,
                 nlcdClass), factor)

birds <- birds %>% 
  mutate_at(vars(pointCountMinute, observerDistance, startCloudCoverPercentage, endCloudCoverPercentage, startRH,
                 endRH, observedAirTemp, kmPerHourObservedWindSpeed), as.numeric)



birdsdiversity_byplotdate <- birds %>%
  group_by(siteID, plotID, observation_datetime) %>%
  summarise(sp_richness = length(unique(taxon_name)), 
            total_abundance = sum(abundance, na.rm=TRUE))

birdsdiversity_byplot <- birdsdiversity_byplotdate %>%
  group_by(siteID, plotID) %>%
  summarise(mean_sp_richness = mean(sp_richness), n_surveys = n())


birdsdiversity_bysite <- birdsdiversity_byplot %>%
  group_by(siteID) %>%
  summarise(mean_mean_sp_richness=mean(mean_sp_richness),
            SE_mean_sp_richness=sd(mean_sp_richness)/sqrt(n()))

# joined with total species richness - all taxa ever recorded for the site 
birdsdiversity_bysite <- full_join(birdsdiversity_bysite,
                                   birds %>%
                                    group_by(siteID) %>%
                                    summarise(total_sp_richness=length(unique(taxon_name))),
                                     by="siteID")

summary(birdsdiversity_bysite)

plot(birdsdiversity_bysite$mean_mean_sp_richness, birdsdiversity_bysite$total_sp_richness)
# not bad


#############################################################################
## Smammals
#############################################################################

# unit is
# "unique individuals per 100 trap nights per plot per month"

smammals <- neonDivData::data_small_mammal

names(smammals)[which(names(smammals)=='value')] <- "abundance"    

smammals <- smammals %>% 
  mutate_at(vars(siteID, plotID, taxon_id, taxon_name, taxon_rank, variable_name,
                 unit, nativeStatusCode, nlcdClass), factor)


smammalsdiversity_byplotdate <- smammals %>%
  group_by(siteID, plotID, observation_datetime) %>%
  summarise(sp_richness = length(unique(taxon_name)), 
            total_abundance = sum(abundance, na.rm=TRUE))

smammalsdiversity_byplot <- smammalsdiversity_byplotdate %>%
  group_by(siteID, plotID) %>%
  summarise(mean_sp_richness = mean(sp_richness), n_surveys = n())


smammalsdiversity_bysite <- smammalsdiversity_byplot %>%
  group_by(siteID) %>%
  summarise(mean_mean_sp_richness=mean(mean_sp_richness),
            SE_mean_sp_richness=sd(mean_sp_richness)/sqrt(n()))

# joined with total species richness - all taxa ever recorded for the site 
smammalsdiversity_bysite <- full_join(smammalsdiversity_bysite,
                                      smammals %>%
                                     group_by(siteID) %>%
                                     summarise(total_sp_richness=length(unique(taxon_name))),
                                   by="siteID")

summary(smammalsdiversity_bysite)

plot(smammalsdiversity_bysite$mean_mean_sp_richness, smammalsdiversity_bysite$total_sp_richness)
# not bad

###############################################################################
## join
###############################################################################

dat <- rbind(
  data.frame(group=rep("plant",dim(plantdiversity_bysite)[1]),plantdiversity_bysite),
  data.frame(group=rep("beetle",dim(beetlesdiversity_bysite)[1]),beetlesdiversity_bysite),
  data.frame(group=rep("bird",dim(birdsdiversity_bysite)[1]),birdsdiversity_bysite),
  data.frame(group=rep("smammal",dim(smammalsdiversity_bysite)[1]),smammalsdiversity_bysite))

  
total.species.richness.bysite <- pivot_wider(dat[,c(1,2,5)], names_from = group,
            values_from = total_sp_richness)

summary(total.species.richness.bysite)

#violin plot
ggplot(dat, aes(x = group, y = total_sp_richness, col=group)) +
  geom_violin(aes(fill = group), show.legend = F) +
  theme_classic(base_size = 18) +
  ylab("total species richness over all sampling per site") +
  ggtitle("Distributions of total species richness across sites")  



##############################################################################
## Productivity
##############################################################################

# prod <- loadByProduct("DP1.10023.001", 
#                       token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
# )

head(prod$hbp_massdata)

biomass <- prod$hbp_massdata
sort(unique(biomass$siteID))
# all 47 sites are here

sort(unique(biomass$herbGroup))
sort(unique(biomass$siteID[which(biomass$herbGroup=="All herbaceous plants")]))
# only 15 of 47 so need to look at groups separately
sort(unique(biomass$siteID[which(biomass$herbGroup=="Annual and Perennial Forbs")]))
# all 47
sort(unique(biomass$siteID[which(biomass$herbGroup=="N-fixing Plants")]))
#all 47
sort(unique(biomass$siteID[which(biomass$herbGroup=="Leguminous Forbs")]))
# missing 2
sort(unique(biomass$siteID[which(biomass$herbGroup=="Woody-stemmed Plants")]))
# all 47

biomass.byplot <- biomass %>%
  group_by(siteID, plotID) %>%
  summarise(nplots = n(), 
            mean.forbs=mean(dryMass[herbGroup=="Annual and Perennial Forbs"], na.rm=TRUE),
            mean.woodyplants=mean(dryMass[herbGroup=="Woody-stemmed Plants"], na.rm=TRUE),
            mean.Nfixers=mean(dryMass[herbGroup=="N-fixing Plants"], na.rm=TRUE),
            mean.legumes=mean(dryMass[herbGroup=="Leguminous Forbs"], na.rm=TRUE)) 

# subplots not listed here, so I'm assuming just one per plot?
# lots of variability across plots within one site 
# probably should separate based on habitat type of plot (within site)

biomass.bysite <- biomass.byplot %>%
  group_by(siteID) %>%
  summarise(n.plots=n(), mean_mean_forbs = mean(mean.forbs, na.rm=TRUE),
            mean_mean_woodyplants = mean(mean.woodyplants, na.rm=TRUE),
            mean_mean_Nfixers = mean(mean.Nfixers, na.rm=TRUE),
            mean_mean_legumes = mean(mean.legumes, na.rm=TRUE))

            

dat <- pivot_longer(biomass.bysite, cols= c(mean_mean_forbs, mean_mean_woodyplants,
                                            mean_mean_Nfixers, mean_mean_legumes),
                    names_to="type", values_to="biomass")
dat$ type <- gsub("mean_mean_", "",dat$type )

#violin plot
ggplot(dat, aes(x = type, y = biomass, col=type)) +
  geom_violin(aes(fill = type), show.legend = F) +
  theme_classic(base_size = 18) +
  ylab("biomass means (of plot means) per site") +
  ggtitle("Distributions of biomass across sites")  





summary(factor(prod$hbp_perbout$exclosure))
# some plots have an exclosure?


##############################################################################
## Site information
##############################################################################

plotsummary <- plantsDD %>%
  group_by(plotID, siteID, latitude, longitude, elevation) %>%
  summarise(num.days=length(unique(observation_datetime)), 
            num.plant.occ=length(unique(paste(year(observation_datetime),month(observation_datetime)))))

sitesummary <-  plantsDD %>%
  group_by(siteID) %>%
  summarise(num.plots= length(unique(plotID)),
            latitude = mean(latitude), 
            longitude = mean(longitude), 
            elevation = mean(elevation),
            num.days = length(unique(observation_datetime)), 
            num.months = length(unique(paste(year(observation_datetime),month(observation_datetime)))))

sitesummary$habitat <- rep(NA, dim(sitesummary)[1])
for(i in 1:dim(sitesummary)[1]){
  sitesummary$habitat[i] <- paste(unique(plantsDD$nlcdClass[which(plantsDD$siteID==sitesummary$siteID[i])]),collapse=",")
}

sitesummary$color <- hue_pal()(length(unique(sitesummary$siteID))) 

maps::map("usa")
map("state",lwd=2, add=TRUE)#mar=c(1, 1, 1, 1)) #c(bottom, left, top, right)
points(sitesummary$longitude,sitesummary$latitude,pch=19, col=sitesummary$color)  
text(sitesummary$longitude,sitesummary$latitude-0.8, sitesummary$siteID,cex=0.8, col=sitesummary$color)
title("NEON terrestrial sites")

# export data 
write.csv(plantsDD, file="AllPlantData.csv", row.names = F)
write.csv(sitesummary ,file="SiteSummaries.csv",  row.names = F)
write.csv(beetles ,file="AllBeetleData.csv",  row.names = F)
write.csv(birds ,file="AllBirdData.csv",  row.names = F)
write.csv(beetles ,file="AllBeetleData.csv",  row.names = F)
write.csv(smammals ,file="AllSmammalData.csv",  row.names = F)
write.csv(biomass ,file="AllProductivityData.csv",  row.names = F)
write.csv(prod$hbp_perbout, file="AllProductivitySampleData.csv", row.names=F)



##############################################################################
### Summary weather statistics
##############################################################################
# Present summary statistics for biometeorological variables for NEON weather 
# stations at core TIS sites. Statistics will include means, standard deviations, 
# maxima, and minima for periods of days, months, and years. Engineering-grade 
# product only.

weather <- loadByProduct("DP4.00001.001",
                      token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

sort(unique(weather$wss_daily_precip$siteID))
sort(unique(weather$wss_daily_humid$siteID))
sort(unique(weather$wss_daily_temp$siteID))
# only 20 sites

##############################################################################
# Soil geochem and properties
##############################################################################

# this is just initial characterization:
soil <- loadByProduct("DP1.10047.001",
                      token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

soil_perplot <- soil$spc_perplot 

soil_perplot %>%
  group_by(nlcdClass) %>%
  summarise(n.site = length(unique(siteID)))


# basic data
save(plantsDD, beetles, birds, smammals, prod, soil, file="BasicNEONdata.RData")
# eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA

names(soil$spc_biogeochem)

biogeochem <- soil$spc_biogeochem %>%
  select(siteID, plotID, collectDate, horizonID, horizonName,  caNh4d, kNh4d, 
         mgNh4d, naNh4d, cecdNh4, phCacl2)

# prob only want to look at the horizon where the microbial data comes from
sort( unique(biogeochem$horizonName)) #>300 horizon names. might have to match up each sample?

biogeo_n_plot_dates <- biogeochem %>%
  group_by(siteID, plotID) %>%
  summarise(n.dates=length(unique(collectDate)))

summary(biogeo_n_plot_dates)
# each plot sampled once

biogeo_n_plot_dates <- biogeochem %>%
  group_by(siteID, plotID) %>%
  summarise(dates=unique(collectDate))
summary(biogeo_n_plot_dates)
# dates range from 9/21/2015 to 10/20/2021

write.csv(soil$spc_biogeochem, file="SoilBiogeochem_initial.csv")

soil_periodic <- loadByProduct("DP1.10086.001",
                      token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

save(soil, soil_periodic, file="SoilData.RData")

# I think this will tell about sampling
soil_persamp <- soil_periodic$sls_bgcSubsampling %>%
  group_by(siteID, plotID) %>%
  summarise(n.dates=length(unique(collectDate)))
# 517 plots

soilchem_persamp <- soil_periodic$sls_soilChemistry %>%
  group_by(siteID, plotID) %>%
  summarise(n.dates=length(unique(collectDate)))
# 513 plots

soilcore_persamp <- soil_periodic$sls_soilCoreCollection %>%
  group_by(siteID, plotID, nlcdClass) %>%
  summarise(n.dates=length(unique(collectDate)))


write.csv(soil_periodic$sls_soilChemistry, file="SoilChemPeriodic.csv")
write.csv(soil_periodic$sls_soilMoisture, file="SoilMoisturePeriodic.csv")
write.csv(soil_periodic$sls_soilpH, file="SoilpHPeriodic.csv")

##############################################################################
# Microbial Data
##############################################################################



microbial <- loadByProduct("DP1.10081.001",
                      token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)






##############################################################################
# Precipitation
##############################################################################

# Bulk precipitation collected using up to three methods - primary, secondary, 
# and throughfall. Bulk precipitation is determined at five- and thirty-minute 
# intervals for primary precipitation and at one- and thirty-minute intervals for 
# secondary and throughfall precipitation.


precip <- loadByProduct("DP1.00006.001", startdate = "2020-01", enddate = "2022-12",
                         token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)
## 26 GB. need to subset.

# tried "2020-01", enddate = "2022-12",
# still 20 GB. Tried to download and "vector memory exhausted (limit reached?)"



##############################################################################
# Temperature
##############################################################################
# Triple aspirated air temperature
# Air temperature, available as one- and thirty-minute averages derived from 
# triplicate 1 Hz temperature observations. Observations are made by sensors 
# located at the top of the tower infrastructure. Temperature observations are 
# made by three platinum resistance thermometers, which are housed together in a 
# fan aspirated shield to reduce radiative biases.

temp <- loadByProduct("DP1.00003.001", startdate = "2019-01", enddate = "2022-12",
                      token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)
 



################# Additional things to download:

# bundled atmospheric data:
# zipsByProduct("DP4.00200.001",
#                       token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
# )
## 405 GB!!! need to subset somehow. what do I want?





# 

# Vegetation indices - spectrometer - mosaic
# DP3.30026.001
# Bundled vegetation indices derived from surface directional reflectance 
# (DP1.30006.001) containing commonly used indices that are used as proxies for 
# vegetation health. Indices in this data product include the Normalized Difference 
# Vegetation Index (NDVI), Enhanced Vegetation Index (EVI), Atmospherically 
# Resistant Vegetation Index (ARVI), Photochemical Reflectance Index (PRI), and 
# Soil Adjusted Vegetation Index (SAVI). 

# fPAR - spectrometer - flightline
# DP2.30014.001 
# The fraction of incident photosynthetically active radiation (400-700 nm) 
# absorbed by the green elements of a vegetation canopy


# IR biological temperature
# DP1.00005.001
# Infrared temperature, available as one- and thirty-minute averages of 1 Hz 
# observations. Biological temperature (i.e. surface temperature) is measured via 
# IR temperature sensors located in the soil array and at multiple heights on the 
# tower infrastructure.

# Leaf surface area?

# Photosynthetically active radiation (PAR)
# DP1.00024.001
# Photosynthetically Active Radiation (PAR) observations represent the radiation 
# flux at wavelengths between 400-700 nm, which constitute the wavelengths that 
# drive photosynthesis. This data product is available as one- and thirty-minute 
# averages of 1 Hz observations. Observations are made by sensors located at multiple 
# heights on the tower infrastructure and by sensors located on the aquatic meteorological station.


# Relative humidity
# DP1.00098.001
# Relative humidity, temperature, and dew or frost point temperature, available 
# as one- and thirty-minute averages of 1 Hz observations. Observations are made 
# by sensors located at the top of the tower infrastructure, in the soil array, 
# and on the aquatic meteorologic station.

# longwave and shortwave radiation?

# Single aspirated air temperature
# DP1.00002.001
# Air temperature, available as one- and thirty-minute averages of 1 Hz 
# observations. Observations are made by sensors located at multiple heights on 
# the tower infrastructure

# Slope and Aspect - LiDAR
# DP3.30025.001
# Slope is a ratio of rise over run (height over distance) of the bare earth 
# elevation product given in degrees; aspect is the direction of the steepest 
# slope of the bare earth elevation product

# Soil CO2 concentration
# DP1.00095.001
# CO2 concentration in soil air at various depth below the soil surface starting 
# at approximately 2 cm. Data are from all five Instrumented Soil Plots per site 
# and presented as 1-minute and 30-minute averages.

# Soil heat flux plate
# DP1.00040.001
# The amount of thermal energy moving by conduction across an area of soil in a 
# unit of time. Measured at three of the five sensor-based soil plots at each 
# terrestrial site..

# Soil temperature
# DP1.00041.001
# Temperature of the soil at various depth below the soil surface from 2 cm up to 
# 200 cm at non-permafrost sites (up to 300 cm at Alaskan sites). Data are from 
# all five Instrumented Soil Plots per site and presented as 1-minute and 30-minute averages.

# Soil water content and water salinity
# DP1.00094.001
# Soil volumetric water content and an index of salinity at various depths below 
# the soil surface from 6 cm down to 200 cm at non-permafrost sites (down to 300 
# cm at Alaskan sites). Data are from all five Instrumented soil plots per site 
# and presented as 1-minute and 30-minute averages. Measurement depths are not 
# currently being reported correctly in the sensor positions file in the Soil 
# water content and water salinity data product (DP1.00094.001). A readme file 
# (readme_swc_depths.txt) and an associated file containing the sensor installation 
# depths (swc_depthsV2.csv) have been added to the data product download package 
# until the data product algorithms are corrected.



##############################################################################
# Daymet Data
##############################################################################

library(daymetr)
daymet.info <- sitesummary[,c(1,3,4)]


# https://daymet.ornl.gov/overview 
write.csv(daymet.info,file="NEONdaymetinfo.csv",row.names = FALSE)
df_batch <- download_daymet_batch(file_location = "NEONdaymetinfo.csv", 
                                  start = 2013, end = 2021, simplify = TRUE)
head(df_batch)

unique(df_batch$measurement)

df_batch$date <- as.Date(df_batch$yday-1, origin=paste(df_batch$year,"01-01",sep="-"))
df_batch$month <- month(df_batch$date)


daymetsummaries.long <- df_batch %>%
  group_by(site,measurement) %>%
  summarise(mean = mean(value))

daymetmeans <- pivot_wider(daymetsummaries.long, 
            names_from = measurement,
            values_from = mean)

write.csv(daymetmeans, file="DaymetEnvDataMeans.csv", row.names = F)

# Day length	dayl	s/day	Duration of the daylight period in seconds per day. 
# This calculation is based on the period of the day during which the sun is above 
# a hypothetical flat horizon
# Precipitation	prcp	mm/day	Daily total precipitation in millimeters per day, 
# sum of all forms converted to water-equivalent. Precipitation occurrence on any 
# given day may be ascertained.
# Shortwave radiation	srad	W/m2	Incident shortwave radiation flux density in 
# watts per square meter, taken as an average over the daylight period of the day. 
# NOTE: Daily total radiation (MJ/m2/day) can be calculated as follows: 
#   ((srad (W/m2) * dayl (s/day)) / l,000,000)
# Snow water equivalent	swe	kg/m2	Snow water equivalent in kilograms per square meter. The amount of water contained within the snowpack.
# Maximum air temperature	tmax	degrees C	Daily maximum 2-meter air temperature in degrees Celsius.
# Minimum air temperature	tmin	degrees C	Daily minimum 2-meter air temperature in degrees Celsius.
# Water vapor pressure	vp	Pa	Water vapor pressure in pascals. Daily average partial pressure of water vapor.

save(daymetsummaries.long, daymetmeans, file="DaymetMeans.RData")




###############################################################################
## Summarize plot / habitat info
###############################################################################

# which plots (within sites) overlap for the different data types?
plantplots <- sort(unique(plantsDD$plotID))
length(plantplots)
#1577
beetleplots <- sort(unique(neon_beetles$plotID))
length(beetleplots)
#511
birdplots <- sort(unique(birds$plotID))
length(birdplots)
#644
smammalplots <- sort(unique(smammals$plotID))
length(smammalplots)
#326
biomassplots <- sort(unique(biomass$plotID))
length(biomassplots)
#1628

# habitat types
sort(unique(plantsDD$nlcdClass))

plant.habitat.plots <- plantsDD %>%
  group_by(nlcdClass) %>%
  summarise(plant.plots = length(unique(plotID)))

habitat.plots <- full_join(plant.habitat.plots,
                           neon_beetles %>%
                             group_by(nlcdClass) %>%
                             summarise(beetle.plots = length(unique(plotID))))

habitat.plots <- full_join(habitat.plots,
                           birds %>%
                             group_by(nlcdClass) %>%
                             summarise(bird.plots = length(unique(plotID))))

habitat.plots <- full_join(habitat.plots,
                           smammals %>%
                             group_by(nlcdClass) %>%
                             summarise(smammal.plots = length(unique(plotID))))

habitat.plots <- full_join(habitat.plots,
                           soil_perplot %>%
                             group_by(nlcdClass) %>%
                             summarise(soil.plots = length(unique(plotID))))



# How many of the smammal plots are also in the other datasets?
summary(smammalplots %in% birdplots)
summary(smammalplots %in% beetleplots)
summary(smammalplots %in% plantplots)



plant.site.habitats <- plantsDD %>%
  group_by(siteID, nlcdClass) %>%
  summarise(plant.n = length(unique(plotID)))

site.habitat.plots <- full_join(plant.site.habitats,
  birds %>%
  group_by(siteID, nlcdClass) %>%
  summarise(bird.n = length(unique(plotID)))
)

site.habitat.plots <- full_join(site.habitat.plots,
                                neon_beetles %>%
                                  group_by(siteID, nlcdClass) %>%
                                  summarise(beetle.n = length(unique(plotID)))
)

site.habitat.plots <- full_join(site.habitat.plots,
                                smammals %>%
                                  group_by(siteID, nlcdClass) %>%
                                  summarise(smammal.n = length(unique(plotID)))
)

site.habitat.plots <- full_join(site.habitat.plots,
                                soil_perplot %>%
                                  group_by(siteID, nlcdClass) %>%
                                  summarise(soil.n = length(unique(plotID)))
)

print(site.habitat.plots,n=113)
# number of unique plot IDs - could be sampled for differing amounts of time.

# per habitat type, how many sites have at least one plot with data on all 4 taxa?
# (not necessarily same plots but must be within the same habitat type)

site.habitat.plots %>%
  drop_na() %>%
  group_by(nlcdClass) %>%
  summarise(n.sites = n())
  



#############################################################################
# Save raw data separately from neonDivData

prod, 
  