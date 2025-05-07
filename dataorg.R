###############################################################################
## Downloading, organizing and taking first look at NEON data
###############################################################################

# plant cover, beetles, birds, and small mammals are from the neonDivData package
# plant productivity, and soil data are directly from NEON using neonUtilities package

library(tidyverse)
library(neonUtilities)
library(neonDivData)
library(maps)
library(scales)




###############################################################################
## Look at plant cover and diversity using neonDivData
###############################################################################

# reloaded the neonDivData package to get updated data
# install.packages("neonDivData", repos = 'https://daijiang.r-universe.dev')
# these data go through 2022-10-26

plantsDD <- neonDivData::data_plant

plant.sites <- sort(unique(plantsDD$siteID))
plant.plots <- sort(unique(plantsDD$plotID))
length(plant.plots)
#1579

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
# truncate to first 4 letters, it removes subspecies, but it's a problem 
# because 4-letter codes are not unique - "ABFR" and "ABFR2" are 2 different species
# Abies fraseri"    "Abronia fragrans"
# instead to get unique species will truncate species to first 2 words


# add columms of genus names that truncates to first 1 word & species names
# truncates to first 2 words
# some are IDed to 'species group' so these might be misleading?
plantsDD$genus <- str_split_i(plantsDD$taxon_name," ",1)
plantsDD$species <- paste(plantsDD$genus , str_split_i(plantsDD$taxon_name," ",2))

species.codes <- plantsDD %>%
  group_by(species, genus, taxon_rank) %>%
  summarise(codes=str_flatten(unique(taxon_id),collapse=" | "))

# write file of species names so Jed can look up plant traits:
write.csv(species.codes, file="PlantSpecies.csv", row.names = F)

# <--------------------------- need to deal with this up for diversity estimation
### Some were not ID'ed to species
#  [983] Calystegia sp. 
# also some say "sp." and some say "spp."
# will probably want to make sure that something from this genus is already 
# counted in the diversity estimate or make it's own species
# I checked, and all those with 'sp.' and 'spp.', say 'genus' under taxon_rank. 
# I'm leaving those for now and will deal with them when do diversity estimation
# for each plot


# Also downloaded raw data but currently using the DivData
# plants.raw <- loadByProduct("DP1.10058.001",
#                        token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
# )
                            
# endDate 6/24/2013 - 12/1/2022






###############################################################################
## Look at diversity in beetles using neonDivData
###############################################################################

beetlesDD <- neonDivData::data_beetle

names(beetlesDD)[which(names(beetlesDD)=='value')] <- "abundance"    #"count per trap day"

beetlesDD <- beetlesDD %>% 
  mutate_at(vars(taxon_id, taxon_name, taxon_rank, nativeStatusCode,
                 nlcdClass), factor)

beetlesDD$trappingDays <- as.numeric(beetlesDD$trappingDays)
# what does this mean? How should it be accounted for?

length(unique(beetlesDD$taxon_name))
# [1] 793 unique taxa

dat <- beetlesDD %>%
  group_by(taxon_name,taxon_rank) %>%
  summarise(unique(taxon_id))
summary(dat$taxon_rank)
# family      genus        order      species speciesgroup speciesGroup    subfamily     subgenus 
# 4           39            1          665            1            2            1           21 
# subspecies  tribe 
# 54            5 

# Clinidium apertum apertum ID'ed to subspecies but says genus, so change to subspecies
beetlesDD$taxon_rank[which(beetlesDD$taxon_name=="Clinidium apertum apertum")] <- "subspecies"
# all rest of those ID'ed to genus says genus sp.

# will have to think more about all those not ID'ed to species - the subgenus adds extra complication

### First initial look at diversity (counting each unique taxon_id as separate - 
# which it isn't since genus and subgenus and subspecies etc)
# this is getting species richness per plot per day
# Also need to think about date is appropriate because sampling occasion could be over several days?
beetlesdiversity_byplotdate <- beetlesDD %>%
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
  beetlesDD %>%
    group_by(siteID) %>%
    summarise(total_sp_richness=length(unique(taxon_name))),
  by="siteID")

summary(beetlesdiversity_bysite)
plot(beetlesdiversity_bysite$mean_mean_sp_richness, beetlesdiversity_bysite$total_sp_richness)
# mean and total not very well correlated


#############################################################################
## Bird data
#############################################################################


birdsDD <- neonDivData::data_bird

names(birdsDD)[which(names(birdsDD)=='value')] <- "abundance"    #"count of individuals"

birdsDD <- birdsDD %>% 
  mutate_at(vars(siteID, plotID, taxon_id, taxon_name, taxon_rank, targetTaxaPresent, detectionMethod,
                 visualConfirmation, sexOrAge, observedHabitat, nativeStatusCode,
                 nlcdClass), factor)

birdsDD <- birdsDD %>% 
  mutate_at(vars(pointCountMinute, observerDistance, startCloudCoverPercentage, endCloudCoverPercentage, startRH,
                 endRH, observedAirTemp, kmPerHourObservedWindSpeed), as.numeric)


length(unique(birdsDD$taxon_name))
# [1] 584 unique taxa

dat <- birdsDD %>%
  group_by(taxon_name,taxon_rank) %>%
  summarise(unique(taxon_id))
summary(dat$taxon_rank)
# class       family      genus      species speciesGroup    subfamily   subspecies 
# 1           16           16          539            3            1            8 


#### First look, but will need to deal with sampling and ID issues:
birdsdiversity_byplotdate <- birdsDD %>%
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
                                   birdsDD %>%
                                    group_by(siteID) %>%
                                    summarise(total_sp_richness=length(unique(taxon_name))),
                                     by="siteID")

summary(birdsdiversity_bysite)

plot(birdsdiversity_bysite$mean_mean_sp_richness, birdsdiversity_bysite$total_sp_richness)



#############################################################################
## Smammals from DivData
## Also downloaded raw data below
#############################################################################

# unit is
# "unique individuals per 100 trap nights per plot per month"

smammalsDD <- neonDivData::data_small_mammal

names(smammalsDD)[which(names(smammalsDD)=='value')] <- "abundance"    

smammalsDD <- smammalsDD %>% 
  mutate_at(vars(siteID, plotID, taxon_id, taxon_name, taxon_rank, variable_name,
                 unit, nativeStatusCode, nlcdClass), factor)

length(unique(smammalsDD$taxon_name))
# [1] 154 unique taxa

dat <- smammalsDD %>%
  group_by(taxon_name,taxon_rank) %>%
  summarise(unique(taxon_id))
summary(dat$taxon_rank)
# genus    species subspecies 
# 18        135          1 

### First look at species richness
smammalsdiversity_byplotdate <- smammalsDD %>%
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
                                      smammalsDD %>%
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
## Raw data on Smammal trapping - each capture, etc 
##############################################################################

# I had already downloaded and cleaned and saved the captures
# see  "~/Documents/NEON data/NEON_count-small-mammals/NEONsmammals/NEONdataorg.R"
# called captures in "/Documents/NEON data/NEON_count-small-mammals/NEONsmammals/NEONsmammalcaptures.RData"
# Downloaded it last year, so it only goes through 2021
# that's past when they stopped hanta testing, but doens't match up to rest of diversity data

smammal.capturedata <- loadByProduct("DP1.10072.001", 
                       token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

names(smammal.capturedata)

save(smammal.capturedata, file="rawSmammalCaptureData.RData")

#see "CleanOrganizeSmammalCaptureData.R" for cleaning and organization to estimate abundances


##############################################################################
## Rodent Pathogen Status, hantavirus 
##############################################################################

# hantavirus data only goes through 2019 or 2020. then they switched to 
# tick borne diseases
# if want to match up will need to make sure that things didn't differ a lot
# over time at these sites

hanta <- loadByProduct("DP1.10064.001", 
                       token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

names(hanta)

head(hanta$rpt_bloodtesting)
hanta.bloodtesting <- as.tibble(hanta$rpt_bloodtesting) %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  mutate_at(vars(siteID, plotID, sampleCondition, testPathogenName, testResult), factor)

save(hanta.bloodtesting, file="HantavirusTestResults.RData")



##############################################################################
# Rodent pathogen status, tick-borne
# DP1.10064.002
##############################################################################

# Presence/absence of tick-borne diseases (or antibodies to same) in each single
# rodent sample from 2020-onward. Prior to 2020, the protocol was used to detect
# hantavirus; these data are available in the Rodent-borne pathogen status data 
# product (DP1.10064.001).

TBD <- loadByProduct("DP1.10064.002", 
                       token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

names(TBD)

head(TBD$rpt2_pathogentesting)
TBD.bloodtesting <- as_tibble(TBD$rpt2_pathogentesting) %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  mutate_at(vars(siteID, plotID, sampleCondition, testPathogenName, testResult,
                  sampleCode, testingID), factor)

save(TBD.bloodtesting, file="TickBornePathogenTestResults.RData")


##############################################################################
# Ticks sampled using drag cloths
# DP1.10093.001
##############################################################################
# Abundance and density of ticks collected by drag and/or flag sampling 
# (by species and/or lifestage)

ticks <- loadByProduct("DP1.10093.001", 
                     token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)
names(ticks)
# [1] "categoricalCodes_10093"      "citation_10093_RELEASE-2025" "issueLog_10093"             
# [4] "readme_10093"                "tck_fielddata"               "tck_taxonomyProcessed"      
# [7] "validation_10093"            "variables_10093" 

# I will use tck_fielddata which lists how many ticks were caught but not IDed
# some were IDed and that's in tck_taxonomyProcessed but it wasn't all of them

tick_sampling <- ticks$tck_fielddata
save(tick_sampling, file="TickSampling.RData")

##############################################################################
# Tick pathogen status
# DP1.10092.001
##############################################################################

# Presence/absence of a pathogen in each single tick sample
tickPathogenStatus <- loadByProduct("DP1.10092.001", 
                       token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

head(tickPathogenStatus$tck_pathogen)

save(tickPathogenStatus, file="TickPathogenStatus.RData")




##############################################################################
## Productivity
## Raw data downloaded using neonUtilities
##############################################################################

prod <- loadByProduct("DP1.10023.001", 
                       token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

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
# some plots have an exclusores - to exclude cattle and were measured multiple times per year


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
write.csv(beetlesDD ,file="AllBeetleData.csv",  row.names = F)
write.csv(birdsDD ,file="AllBirdData.csv",  row.names = F)
write.csv(smammalsDD ,file="AllSmammalData.csv",  row.names = F)
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

# won't use this. instead got Daymet data. See below.





##############################################################################
# Soil geochem and properties
# Only done once - initial characterization
##############################################################################

soil <- loadByProduct("DP1.10047.001",
                      token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

soil_perplot <- soil$spc_perplot 

soil_perplot %>%
  group_by(nlcdClass) %>%
  summarise(n.site = length(unique(siteID)))



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





##############################################################################
# Soil Data done periodically along with microbial sampling
##############################################################################

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
# Soil microbe sequence metadata
##############################################################################

# they are redoing all their sequencing
# release doesn't go through end of 2019
# will download provisional 
microbial <- loadByProduct("DP1.10108.001", include.provisional = TRUE,
                      token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)




##############################################################################
# Soil microbe biomass
##############################################################################

microbial.biomass.raw <- loadByProduct("DP1.10104.001", 
                           token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

names(microbial.biomass.raw)
# scaled and not scaled dataframes. not sure what we want

# are they the same samples? some.
summary(microbial.biomass.raw$sme_microbialBiomass$sampleID %in% 
          microbial.biomass.raw$sme_scaledMicrobialBiomass$sampleID)
#    Mode   FALSE    TRUE 
# logical    1135    4748 
summary(microbial.biomass.raw$sme_scaledMicrobialBiomass$sampleID %in%
          microbial.biomass.raw$sme_microbialBiomass$sampleID)
#    Mode   FALSE    TRUE 
# logical    2326    4741 

write_csv(microbial.biomass.raw$sme_microbialBiomass, file="MicrobialBiomass.csv")
write_csv(microbial.biomass.raw$sme_scaledMicrobialBiomass, file="ScaledMicrobialBiomass.csv")

save(microbial.biomass.raw, file="MicrobialBiomassRaw.RData")

###############################################################################
## Root biomass and chemistry, Megapit
###############################################################################


root.biomass.raw <- loadByProduct("DP1.10066.001", 
                                       token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

save(root.biomass.raw, file="RootBiomassRaw.Rdata")

write_csv(root.biomass.raw$mpr_perrootsample, file="RootBiomass_perRootSample.csv")
write_csv(root.biomass.raw$mpr_perdepthincrement, file= "RootBiomass_perDepth.csv")
write_csv(root.biomass.raw$mpr_perpitprofile, file="RootBiomass_perPitProfile.csv")
write_csv(root.biomass.raw$mpr_carbonNitrogen, file="RootBiomass_CarbonNitrogen.csv")



###############################################################################
## Root biomass and chemistry, periodic
###############################################################################


root.biomass.periodic.raw <- loadByProduct("DP1.10067.001", 
                                  token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE4NzMwMzE4NjksImlhdCI6MTcxNTM1MTg2OSwiZW1haWwiOiJhbmdlbGEubHVpc0B1bW9udGFuYS5lZHUifQ.YLxLG3mCbxvV8RTI2amQFiOum--sxt5q5PgL4UIWaOnsILZTCu1kBRbAkoroYJEs5vNeDOI_6Tgk2913yV7NiA"
)

save(root.biomass.periodic.raw, file="RootBiomassPeriodicRaw.Rdata")



###############################################################################
# Save basic data for further Analyses
###############################################################################

save(plantsDD, beetlesDD, birdsDD, smammalsDD, prod, soil, 
     soil_periodic, microbial, hanta.bloodtesting,
     file="BasicNEONdata.RData")




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

# Mosquito sampled from CO2 traps
# DP1.10043.001
# Taxonomically identified mosquitoes and the plots and times from which they were collected

# Mosquito pathogen status
# DP1.10041.001
# Presence/absence of a pathogen in a single mosquito sample (pool)

# Tick pathogen status
# DP1.10092.001	
# Presence/absence of a pathogen in each single tick sample

# Ticks sampled using drag cloths
# DP1.10093.001	
# Abundance and density of ticks collected by drag and/or flag sampling 
# (by species and/or lifestage)


##############################################################################
# Daymet Data
##############################################################################
# initially did this by site, but this is by plot 

library(daymetr)

daymet.info.plot <- plotsummary[,c(1,3,4)]


# https://daymet.ornl.gov/overview 
write.csv(daymet.info.plot,file="NEONdaymetPlotinfo.csv",row.names = FALSE)
df_batch <- download_daymet_batch(file_location = "NEONdaymetPlotinfo.csv", 
                                  start = 2013, end = 2021, simplify = TRUE)
head(df_batch)

unique(df_batch$measurement)

df_batch$date <- as.Date(df_batch$yday-1, origin=paste(df_batch$year,"01-01",sep="-"))
df_batch$month <- month(df_batch$date)

names(df_batch)[1] <- "plotID"

daymetPlotsummaries.long <- df_batch %>%
  group_by(plotID,measurement) %>%
  summarise(mean = mean(value))

daymetPlotmeans <- pivot_wider(daymetPlotsummaries.long, 
                               names_from = measurement,
                               values_from = mean)

write.csv(daymetPlotmeans, file="DaymetPlotEnvDataMeans.csv", row.names = F)

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

save(daymetPlotsummaries.long, daymetPlotmeans, file="DaymetPlotMeans.RData")


#### using the microbial plots
library(daymetr)
daymet.micro.plots <- soil.periodic.merge %>%
  dplyr::select(siteID, plotID, nlcdClass, decimalLatitude, decimalLongitude, elevation) %>%
  distinct()

daymet.batchplots <- daymet.micro.plots[,c(2,4,5)]


# https://daymet.ornl.gov/overview 
write.csv(daymet.batchplots,file="NEONdaymetMicroPlotinfo.csv",row.names = FALSE)
daymetsoilplots <- download_daymet_batch(file_location = "NEONdaymetMicroPlotinfo.csv", 
                                  start = 2018, end = 2022, simplify = TRUE)
head(daymetsoilplots)

unique(daymetsoilplots$measurement)

daymetsoilplots$date <- as.Date(daymetsoilplots$yday-1, origin=paste(daymetsoilplots$year,"01-01",sep="-"))
daymetsoilplots$month <- month(daymetsoilplots$date)

names(daymetsoilplots)[1] <- "plotID"

daymetsoilplots.sums <- daymetsoilplots %>%
  pivot_wider(names_from = measurement, values_from = value) %>%
  mutate(temp = (tmax..deg.c. + tmin..deg.c.)/2) %>%
  mutate(temp.14days = zoo::rollapply(temp, 14, mean, align="right", fill=NA),
         temp.30days = zoo::rollapply(temp, 30, mean, align="right", fill=NA),
         temp.365days = zoo::rollapply(temp, 365, mean, align="right", fill=NA),
         prcp.14days = zoo::rollapply(prcp..mm.day., 14, sum, align="right", fill=NA), # need to make sure there aren't any NAs in the dataset
         prcp.30days = zoo::rollapply(prcp..mm.day., 30, sum, align="right", fill=NA),
         prcp.365days = zoo::rollapply(prcp..mm.day., 365, sum, align="right", fill=NA))

save(daymetsoilplots.sums,file="daymetSoilSums.RData")
# daymetPlotsummaries.long <- daymetsoilplots %>%
#   group_by(plotID,measurement) %>%
#   summarise(mean = mean(value))
# 
# daymetPlotmeans <- pivot_wider(daymetPlotsummaries.long, 
#                                names_from = measurement,
#                                values_from = mean)
# 
# write.csv(daymetPlotmeans, file="DaymetPlotEnvDataMeans.csv", row.names = F)


###############################################################################
## Summarize plot / habitat info
###############################################################################

# which plots (within sites) overlap for the different data types?
plantplots <- sort(unique(plantsDD$plotID))
length(plantplots)
#1579
beetleplots <- sort(unique(beetlesDD$plotID))
length(beetleplots)
#514
birdplots <- sort(unique(birdsDD$plotID))
length(birdplots)
#648
smammalplots <- sort(unique(smammalsDD$plotID))
length(smammalplots)
#328
biomassplots <- sort(unique(biomass$plotID))
length(biomassplots)
#1628

# habitat types
sort(unique(plantsDD$nlcdClass))

plant.habitat.plots <- plantsDD %>%
  group_by(nlcdClass) %>%
  summarise(plant.plots = length(unique(plotID)))

habitat.plots <- full_join(plant.habitat.plots,
                           beetlesDD %>%
                             group_by(nlcdClass) %>%
                             summarise(beetle.plots = length(unique(plotID))))

habitat.plots <- full_join(habitat.plots,
                           birdsDD %>%
                             group_by(nlcdClass) %>%
                             summarise(bird.plots = length(unique(plotID))))

habitat.plots <- full_join(habitat.plots,
                           smammalsDD %>%
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
  birdsDD %>%
  group_by(siteID, nlcdClass) %>%
  summarise(bird.n = length(unique(plotID)))
)

site.habitat.plots <- full_join(site.habitat.plots,
                                beetlesDD %>%
                                  group_by(siteID, nlcdClass) %>%
                                  summarise(beetle.n = length(unique(plotID)))
)

site.habitat.plots <- full_join(site.habitat.plots,
                                smammalsDD %>%
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
  



  