#############################################################################
## Subset & Merge Datasets of plants, biomass, and soils
#############################################################################


# see dataorg.R for initial downloads of the data
# plant cover, beetles, birds, small mammals are from the neonDivData package
# plant productivity, soil initial characterization, and soil periodic data
# are directly from NEON via neonUtilities package

load("BasicNEONdata.RData")
library(tidyverse)





###############################################################################
### Clean / Subset Plant productivity data
###############################################################################

prod$hbp_massdata # biomass data by group per clip sample
prod$hbp_perbout # gives clipArea and whether exclosure, also subplot info (only plotID in mass data)
# match up the two datasets by sampleID

### need to deal with the issue of inconsistent clip area between samples
# most have clipArea == 0.2 with clipWidth=0.1 and clipLength=2
# but some don't
summary(factor(prod$hbp_perbout$nlcdClass[which(prod$hbp_perbout$clipArea>0.2)]))
# cultivatedCrops grasslandHerbaceous          pastureHay 
# 375                   7                  18 
# if we get rid of cultivated crops and pasture hay, then only 7 grasslandHerbaceous
# that have larger clipArea
# 2 of the grasslandHerbaceous say "Agricultural" under plotManagement

summary(prod$hbp_perbout$clipLength*prod$hbp_perbout$clipWidth==prod$hbp_perbout$clipArea)
#    Mode   FALSE    TRUE    NA's 
# logical     213   22023       3 
# 213 clipAreas don't each the product of width and length
# which ones are they?
# [NAs seem to be a few dates sampling didn't happen, so no sampleID - can remove]

summary(factor(prod$hbp_perbout$nlcdClass[which(prod$hbp_perbout$clipLength*prod$hbp_perbout$clipWidth!=prod$hbp_perbout$clipArea)]))
# cultivatedCrops grasslandHerbaceous          pastureHay 
#             192                   3                  18 
# again, removing cultivatedCrops and pastureHay will remove most of the problems. what are the 3 grassland? 
prod$hbp_perbout[which(prod$hbp_perbout$clipLength*prod$hbp_perbout$clipWidth!=prod$hbp_perbout$clipArea & 
                         prod$hbp_perbout$nlcdClass=="grasslandHerbaceous"),]
#   plotID subplotID       clipID  collectDate  exclosure clipLength clipWidth clipArea
# LAJA_001    31_400 LAJA_001_080   2017-04-28       <NA>       1.5       0.7       1
# LAJA_016    31_400 LAJA_016_202   2017-04-28       <NA>       1.5       0.7       1
# LAJA_009    31_400 LAJA_009_208   2017-04-28       <NA>       1.5       0.7       1
# There is no sample ID for these and they aren't in the mass data, so will be removed automatically. 


# first join the perbout info and the mass info by sampleID so just 1 dataframe
prod.mass <- as.tibble(prod$hbp_massdata)
prod.perbout <- as.tibble(prod$hbp_perbout)

# change some columns to factors
prod.mass <- prod.mass %>%
  mutate_at(vars(domainID, siteID, plotType, herbGroup), factor)
prod.perbout <- prod.perbout %>%
  mutate_at(vars(domainID, siteID, subplotID, nlcdClass, plotType, exclosure), factor)

# join dataframes
# remove before joining:  measuredBy, recordedBy, dataQF, samplingImpractical, remarks in the perbout data
# make datetime just day (sometimes have time of day, but don't always so exact datetimes may not line up)
productivity <- left_join(prod.mass %>%
                            mutate(collectDate=floor_date(collectDate, unit="days"), 
                                   setDate=floor_date(setDate, unit="days")), 
                          prod.perbout %>%
                            select(!c(measuredBy, recordedBy, dataQF, samplingImpractical, remarks)) %>%
                            mutate(collectDate=floor_date(collectDate, unit="days"), 
                                   setDate=floor_date(setDate, unit="days"))
                          , 
                          by=c("domainID", "siteID", "namedLocation", "plotID",
                                                        "plotType", "setDate", "collectDate",
                                                        "sampleID", "release", "publicationDate"), 
                          relationship = "many-to-one")

summary(productivity)


# Reduced Data ---------------------------------------------------------------#
# removing nlcdClass==cultivatedCrops; pastureHay, emergentHerbaceousWetlands,
# woodyWetlands
# only using 2019 onwards
# only using data from peak growing season of June, July, August -------------#

productivity <- productivity %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands",
         year(collectDate)>2018, 
         month(collectDate)==6 | month(collectDate)==7 | month(collectDate)==8 ) 
# now only a fourth of the original dataset


productivity[productivity$clipArea!=0.2,]
# leaves 1 KONZ sample that has clipArea==1. 
# Remove it.
productivity <- productivity %>%
  filter(clipArea == 0.2)

# now 19655 rows of data


###############################################################################
### Clean / Subset Plant Cover data
###############################################################################


# Reduced Data ---------------------------------------------------------------#
# removing nlcdClass==cultivatedCrops; pastureHay, emergentHerbaceousWetlands,
# woodyWetlands
# only using 2019 onwards
# only using data from peak growing season of June, July, August -------------#

plant.cover <- plantsDD %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands",
         year(observation_datetime)>2018, 
         month(observation_datetime)==6 | month(observation_datetime)==7 | 
           month(observation_datetime)==8 ) 
# now only a fourth of the original dataset

plant.cover <- plant.cover %>%
  mutate(observation_datetime=floor_date(observation_datetime, unit="days")) %>%
  mutate_at(vars(siteID, plotID, subplotID), factor) %>%
  mutate_at(vars(sample_area_m2, boutNumber), as.numeric)
  

# <---------------------------------------------------------------------------#
# See "ExploringSampling.Rmd"
# found that 7% of the plot-years don't have 6 1m^2 pots. Remove them.

plant.1mplotyears <- plant.cover[grep("_1_",plant.cover$subplotID),] %>%
  group_by(siteID, plotID, year=year(observation_datetime)) %>%
  summarise(num.1m2.plots=length(unique(subplotID)))
summary(factor(plant.1mplotyears$num.1m2.plots))
rm.plot.yr <- plant.1mplotyears[which(plant.1mplotyears$num.1m2.plots<6),]
rm.plot.yr$plot_year = paste(rm.plot.yr$plotID, rm.plot.yr$year, sep=" ")

reduced.plant.cover <- plant.cover %>%
  mutate(plot_year = paste(plotID, year(observation_datetime),sep=" ")) %>%
  filter( (plot_year %in% rm.plot.yr$plot_year) == FALSE)




# 
# ### How do plots line up between datasets ------------------------------------#
# cover.site.info <- plant.cover %>%
#   group_by(siteID) %>%
#   summarise(n.plots = length(unique(plotID)), mean.lat=mean(latitude), mean.long=mean(longitude), mean.elev=mean(elevation),
#             habitats = str_flatten(unique(nlcdClass),collapse=" | ", na.rm=TRUE))
# 
# cover.plot.info <- plant.cover %>%
#   group_by(siteID, plotID) %>%
#   summarise(n.subplots = length(unique(subplotID)), n.years = length(unique(year(observation_datetime))), 
#             habitat = unique(nlcdClass))
# 
# coverplots <- sort(unique(plant.cover$plotID)) #1044
# plantprodplots <- sort(unique(productivity$plotID)) #1431
# length(which(coverplots %in% plantprodplots)) 
# length(which(plantprodplots %in% coverplots)) #428 plots in common across cover and productivity datasets
# # so will want to just make sure plots and subplots are within the same
# # habitat type per site rather than matching up individual plots
# # ----------------------------------------------------------------------------#



###############################################################################
### Clean / Subset Periodic Soil data
###############################################################################
# match up data across spreadsheets with sampleID (and to microbial data)

soil_periodic$sls_soilCoreCollection
soil_periodic$sls_soilChemistry
soil_periodic$sls_soilMoisture
soil_periodic$sls_soilpH
# not all samples were tested for everything
# a lot fewer were tested for Chemistry
# and some samples have multiple rows per dataframe because were tested more 
# than once or using different methods.

# Join data frames ------------------------------------------------#
# First, only subsetting by year - 2019 onward
# not filtering yet by nlcdClass, growing season months, or sampleTiming
# will do that after merge

soil.cores <- soil_periodic$sls_soilCoreCollection %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  filter(
  #        nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
  #        nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands",
       year(collectDate)>2018) %>%
       # month(collectDate)==6 | month(collectDate)==7 | month(collectDate)==8 )  
       ### not subsetting by month because there were a lot more that sampled 
        ## outside this range - a lot in October - can do it later if want to
  select(c(3,4,6:16,19,20,23,25,29,32,34:37,41)) %>%
  mutate_at(vars(siteID, plotID, plotType, nlcdClass, subplotID,
                 horizon, boutType,sampleTiming), factor) 

soil.chem <-  soil_periodic$sls_soilChemistry %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  filter(year(collectDate)>2018) %>%
  # month(collectDate)==6 | month(collectDate)==7 | month(collectDate)==8 )  
  ### not subsetting yet
  select(c(5:7,9:10,16:21)) %>%
  mutate_at(vars(siteID, plotID, plotType, co2Trapped), factor)       
  
ls <- unlist(lapply(as.list(soil.chem$sampleID), function(x){length(which(soil.chem$sampleID==x))}))
summary(factor(ls))
# lots of repeat samples
# 4 rows with sample "BLAN_038-M-5.5-36.5-20200714". 2 different measurements of two types of data. 
# analyticalRepNumber goes from 1 to 4 if sample has 4 rows in data.

# I will take means before merging. 
soil.chem.ave <- soil.chem %>%
  group_by(siteID, plotID, plotType, collectDate, sampleID) %>%
  summarise(d15N = mean(d15N, na.rm=TRUE), 
            organicd13C = mean(organicd13C, na.rm=TRUE),
            nitrogenPercent = mean(nitrogenPercent, na.rm = TRUE),
            organicCPercent = mean(organicCPercent, na.rm = TRUE))
            
  

soil.mois <-  soil_periodic$sls_soilMoisture %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  filter(year(collectDate)>2018) %>%
  # month(collectDate)==6 | month(collectDate)==7 | month(collectDate)==8 )  
  ### not subsetting yet
  select(c(3,4,7,8,12,15:19)) %>%
  mutate_at(vars(siteID, plotID, horizon), factor)       
ls <- unlist(lapply(as.list(soil.mois$sampleID), function(x){length(which(soil.mois$sampleID==x))}))
summary(factor(ls))
# don't need to average. only 1 row per sample


soil.pH <-  soil_periodic$sls_soilpH %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  filter(year(collectDate)>2018) %>%
  # month(collectDate)==6 | month(collectDate)==7 | month(collectDate)==8 )  
  ### not subsetting yet
  select(c(3,4,7,9,13:21)) %>%
  mutate_at(vars(siteID, plotID, horizon), factor)       
ls <- unlist(lapply(as.list(soil.pH$sampleID), function(x){length(which(soil.pH$sampleID==x))}))
summary(factor(ls))
# don't need to average. only 1 row per sample


soil.periodic.merge <- left_join(soil.cores, soil.chem.ave) # 
soil.periodic.merge <- full_join(soil.periodic.merge, soil.mois)     
soil.periodic.merge <- full_join(soil.periodic.merge, soil.pH)     
summary(soil.periodic.merge)

# Still haven't subset to make easier to match up to microbes, etc ----------#
# Need to do at some point!
# Could subset to peakGreenness (instead of by month 6|7|8) -----------------#
# soil.periodic.merge <- soil.periodic.merge %>%
#   filter(sampleTiming=="peakGreenness")
# summary(month(soil.periodic.merge$collectDate))
# month still goes from 1 to 11 - not sure how that's peak greenness
# so may want to subset by month to be consistent?
# Sill haven't subset to nlcdClass ------------------------------------------#



dim(soil.periodic.merge) 
# 19117 rows (samples)
length(which(!is.na(soil.periodic.merge$geneticSampleID)))
# 7330 genetic samples



## See "ExploringSampling.Rmd" to look at how sampling align between soils,
# plant cover, productivity







 
###############################################################################
# Microbial Sampling
###############################################################################

# a lot of these are provisional and should probably be re-downloaded with 
# release versions

names(microbial)
# [1] "categoricalCodes_10108"           "citation_10108_PROVISIONAL"      
# [3] "citation_10108_RELEASE-2024"      "issueLog_10108"                  
# [5] "mmg_soilDnaExtraction"            "mmg_soilMarkerGeneSequencing_16S"
# [7] "mmg_soilMarkerGeneSequencing_ITS" "mmg_soilPcrAmplification_16S"    
# [9] "mmg_soilPcrAmplification_ITS"     "readme_10108"                    
# [11] "validation_10108"                 "variables_10108"                 

summary(microbial$mmg_soilMarkerGeneSequencing_ITS$collectDate)
# 2013-06-27 to 2022-11-07

# some rows in mmg_soilMarkerGeneSequencing_ITS are duplicated
dim(microbial$mmg_soilMarkerGeneSequencing_ITS)
#[1] 15796    41
# if remove uid (row identifier)
dim(distinct(microbial$mmg_soilMarkerGeneSequencing_ITS[,-1]))
#[1] 15711    40


microbe.ITS <- as.tibble(distinct(microbial$mmg_soilMarkerGeneSequencing_ITS[,-1])) %>% # remove duplicated rows before proceeding
  mutate(collectDate=floor_date(collectDate, unit="days"),
         geneticSampleID= unlist(lapply(as.list(dnaSampleID), 
                                 function(x){paste(str_split_1(x,pattern="-GEN")[1], 
                                                   "-GEN", sep ="")}))) %>%    #str_sub(dnaSampleID, 1, -6)) %>% 
  mutate_at(vars(domainID, siteID, plotID, qaqcStatus), factor) %>%
  filter(year(collectDate)>2018,
         grepl("COMP", dnaSampleID)==F) %>% # removed the one dnaSampleID with "COMP" in it because that means it was pooled and there were other IDs with that plot-date
  select(c(2:4,42,8:10,24)-1)

# check that my version of geneticSampleID looks ok
l.gen <- unlist(lapply(as.list(microbe.ITS$geneticSampleID), 
                       function(x){length(str_split_1(x,pattern="-"))}))
summary(factor(l.gen))
# 5    6 
# 1572 5213 
# string of length both 5 and 6 are normal

ls <- unlist(lapply(as.list(microbe.ITS$geneticSampleID), function(x){length(which(microbe.ITS$geneticSampleID==x))}))
summary(factor(ls))
# 1 
# 6785 
# only 1 row per geneticSampleID - that's good.

summary(microbe.ITS$qaqcStatus)
# Fail Not recorded         Pass 
#  150            0         6635 
# may want to only use the Pass ones? Worry about later
  
gs <- soil.periodic.merge$geneticSampleID[which(!is.na(soil.periodic.merge$geneticSampleID))]
summary(microbe.ITS$geneticSampleID %in% gs)
#    Mode    TRUE 
# logical    6785 
# all the microbial samples are in the soil.periodic.merge data - yay


microbe.ITS.metadata <- left_join(microbe.ITS, soil.cores)


length(unique(soil.periodic.merge$geneticSampleID))
length(unique(microbe.ITS$geneticSampleID))

# add to soil.periodic.merge whether there the geneticSampleID is in the 
# microbe.ITS and qaqc for it.
soil.periodic.merge <- left_join(soil.periodic.merge,
                                 microbe.ITS.metadata %>%
                                   select(geneticSampleID, qaqcStatus))
names(soil.periodic.merge)[which(names(soil.periodic.merge)=="qaqcStatus")] <- "ITS.qaqcStatus"
soil.periodic.merge$ITS.sequence = factor(ifelse(is.na(soil.periodic.merge$ITS.qaqcStatus), "N", "Y"))
soil.periodic.merge <- soil.periodic.merge[,c(1:41, 43, 42)]


### Remove crops, pasture, wetlands -----------------------------------------#
soil.periodic.merge <- soil.periodic.merge %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands")

microbe.ITS.metadata <- microbe.ITS.metadata %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands")

  
  
### Still haven't filtered Time of year -------------------------#


###############################################################################
# Initial Soil Characterization
###############################################################################

# soil is the initial characterization per site

soil.plots <- soil$spc_perplot %>%
  group_by(siteID, plotID, nlcdClass) %>%
  summarise(n.dates = length(unique(collectDate)))

# 15 plots in ABBY that were sampled. each plot was sampled 1 date and has
# different rows per soil horizon in the datasets.

# I won't cut these down at all by date

# need to figure out how to match up the plotID/horizon to microbial data

########### Merge datasets
# "spc_biogeochem"              "spc_bulkdensity"            
# "spc_particlesize"            "spc_perhorizon"              "spc_perplot"

# spc_perplot has one row per plot with plot info
# spc_perhorizon has 2-9 rows per plot, 1 for each horizon
# ditto for spc_biogeochem, particlesize
# But bulkdensity is missing some plots. Has 527/727


# join dataframes plot and horizon info
soil.initial <- left_join(as.tibble(soil$spc_perhorizon) %>%
                            mutate(collectDate=floor_date(collectDate, unit="days")) %>% 
                            mutate_at(vars(domainID, siteID, plotID, pitID, horizonID, horizonName), factor) %>%
                            select(c(3:7,10:13)),
                          
                          as.tibble(soil$spc_perplot) %>%
                            mutate(collectDate=floor_date(collectDate, unit="days")) %>%
                            mutate_at(vars(domainID, siteID, plotType, plotID, pitID, nlcdClass ), factor) %>%
                            select(c(2:4,6,7:18,20:22,27:32)),
                          
                          relationship = "many-to-one")

############### now add biochem
soil.initial.biogeochem <-as.tibble(soil$spc_biogeochem) %>%
                            mutate(collectDate=floor_date(collectDate, unit="days")) %>% 
                            mutate_at(vars(domainID, siteID, plotID, horizonID, horizonName), factor) %>%
                            select(c(3:5,7:9,13:15,31:39,43:45,49:50,54,58:73,79:83,87:92,
                                     96:99,103,107,111,115,119))
dim(soil.initial.biogeochem)
# [1] 3037   60
length(unique(soil.initial.biogeochem$horizonID ))
# 3037
# no duplicates
# FYI: horizonID is unique but, there can be several horizonIDs per plotID-horizonName (with diff top and bottom depths)

# merge
soil.initial <- left_join(soil.initial, soil.initial.biogeochem)

############## now add bulkdensity
soil.initial.bulkdensity <-as.tibble(soil$spc_bulkdensity) %>%
                            mutate(collectDate=floor_date(collectDate, unit="days")) %>% 
                            mutate_at(vars(domainID, siteID, plotID, horizonID, horizonName), factor) %>%
                            select(c(3:5,7:9,14:22))
dim(soil.initial.bulkdensity)
length(unique(soil.initial.bulkdensity$horizonID))
# a few horizonIDs repeated
ls <- unlist(lapply(as.list(soil.initial.bulkdensity$horizonID), 
                     function(x){length(which(soil.initial.bulkdensity$horizonID==x))}))
summary(factor(ls))
as.data.frame(soil.initial.bulkdensity[which(ls>1),])
# 3 are repeated because diff analyses with diff SampleType: clod vs compliant cavity
# take the mean (but it's there's only 1 entry per variable type)

soil.initial.bulkdensity <- soil.initial.bulkdensity %>%
  group_by(domainID, siteID,   plotID, collectDate, horizonID, horizonName) %>%
  summarise(across(2:9, \(x) mean(x, na.rm=TRUE), .names="{.col}"))
# returns NaN instead of NA for missing values now

# merge
soil.initial <- left_join(soil.initial,soil.initial.bulkdensity)




############ now add particlesize

soil.initial.particlesize <- as.tibble(soil$spc_particlesize) %>%
                              mutate(collectDate=floor_date(collectDate, unit="days")) %>% 
                              mutate_at(vars(domainID, siteID, plotID, horizonID, horizonName), factor) %>%
                              select(c(3:5,7:9,12:14,18:31))

dim(soil.initial.particlesize)
length(unique(soil.initial.particlesize$horizonID))
# no duplicates

# merge
soil.initial <- left_join(soil.initial, soil.initial.particlesize)


### Remove crops, pasture, wetlands -----------------------------------------#

soil.initial <- soil.initial %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands")


# microbe.ITS.metadata has sampleTopDepth and sampleBottomDepth
# need to think about how to match that up with initial sampling 
# "horizonTopDepth"     and     "horizonBottomDepth"
# because they aren't going to be exactly the same


###############################################################################
 
save(plant.cover, productivity, soil.periodic.merge, 
     soil.initial, microbe.ITS.metadata,
     file="ReducedMergedPlantSoilData.RData")

# All datasets have removed nlcdClass crops, pasture, wetlands.
# All but the soil.initial only include 2019 and onwards.
# The plant.cover an productivity datasets have been filtered to include
# only June, July, August
# the soil.periodic and microbe metadata haven't been filtered by time of year