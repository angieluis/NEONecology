#############################################################################
## Subset & Merge Datasets of plants, biomass, and soils
#############################################################################


# see dataorg.R for initial downloads of the data
# plant cover, beetles, birds, small mammals are from the neonDivData package
# plant productivity, soil initial characterization, and soil periodic data
# are directly from NEON via neonUtilities package

load("BasicNEONdata.RData")
library(tidyverse)


# We want to subset the data peak greenness. The periodic soil sampling
# has a variable called "sampleTiming" for which one category is "peakGrennness"
# none of the other datasets have this. Get the date range they called 
# "peakGreenness" for each site.

site.habitat.peakGreenness <- soil_periodic$sls_soilCoreCollection %>%
  mutate(MonthDay = format(as.Date(collectDate), "%m-%d")) %>%
  group_by(siteID, nlcdClass) %>%
  summarise(first.date = min(MonthDay), 
            last.date = max(MonthDay))

# This is not helpful. For some sites very wide. (mid-April through November)
write_csv(site.habitat.peakGreenness, file="SitePeakGreennessDates.csv")


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
# only using 2019 onwards  ---------------------------------------------------#

productivity <- productivity %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands",
         year(collectDate)>2018) 


productivity[productivity$clipArea!=0.2,]
# leaves 1 KONZ sample that has clipArea==1. 
# Remove it.
productivity <- productivity %>%
  filter(clipArea == 0.2)

# now 25653 rows of data





#### Need to subset to peak growing season -----------------------------------#
# Base peak growing season on when highest biomass was recorded per site-habitat.
# Use summed dryMass over different categories (e.g.,Cool Season Graminoids vs
# Bryophyte) per sampleID (after averaging over duplicates).
# Data is now per sampleID
# adding column for the date peak biomass was recorded and if the row's date
# is within x month(s) of the peak date

productivity.persample <-  productivity %>%
  # first group by subsampleID and take mean (because sometimes weighed twice for qa)
  group_by(siteID, nlcdClass, collectDate, year=year(collectDate), plotID, 
           sampleID, subsampleID, subplotID, plotType, plotSize, plotManagement, clipArea, 
           exclosure, elevation, slopeAspect, slopeGradient) %>%
  summarise(sub.dryMass = mean(dryMass)) %>%
  # then group by sampleID to get summed biomass for a sample
  group_by(siteID, nlcdClass, collectDate, year=year(collectDate), plotID, 
           sampleID, subplotID, plotType, plotSize, plotManagement, clipArea, 
           exclosure, elevation, slopeAspect, slopeGradient) %>%
  summarise(sum.dryMass = sum(sub.dryMass)) %>%
  # calculate the date that had the highest biomass for the site-habitat
  group_by(siteID, nlcdClass) %>%
  mutate(peak.date = collectDate[which(sum.dryMass==max(sum.dryMass))][1]) %>%
  ungroup() %>%
  mutate(peak.date = ymd(paste(year, 
                               format(as.Date(peak.date), "%m-%d"), 
                               sep="-")),
         peak.window.2mo = factor(if_else(abs(difftime(ymd(collectDate), peak.date, units="days"))<32 , "Y", "N")),
         # if within 31 days of peak time then "Y" for peak.window.2mo (2 month-window)
         peak.window.3mo = factor(if_else(abs(difftime(ymd(collectDate), peak.date, units="days"))<46 , "Y", "N")),
         # 3 month window
         peak.window.4mo = factor(if_else(abs(difftime(ymd(collectDate), peak.date, units="days"))<62 , "Y", "N")))
         # 4 month window

summary(productivity.persample$peak.window.2mo)
#    N    Y 
# 1437 3921  
# removes about a quarter of the data
summary(productivity.persample$peak.window.3mo)
#    N    Y 
# 1437 3921  
summary(productivity.persample$peak.window.4mo)
#   N    Y 
# 720 4638 






###############################################################################
### Clean / Subset Plant Cover data
###############################################################################


# Reduced Data ---------------------------------------------------------------#
# removing nlcdClass==cultivatedCrops; pastureHay, emergentHerbaceousWetlands,
# woodyWetlands
# only using 2019 onwards ----------------------------------------------------#

plant.cover <- plantsDD %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands",
         year(observation_datetime)>2018) 


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


# Subset to peak growing season based on productivity data ------------------#

# add peak dates

# pull out the peak.date for each site-habitat

site.habitat.peakDate <- productivity.persample %>%
  group_by(siteID, nlcdClass) %>%
  summarise(peakdate = unique(format(as.Date(peak.date), "%m-%d")))
# for a given site, the dates between habitats are similar. 
# So to fill in the missing ones take the mean of the other habitats at that site
site.peakDate <- site.habitat.peakDate %>%
  group_by(siteID) %>%
  summarise(site.peakdate = format(mean(mdy(paste(peakdate, 2019, sep="-") )), "%m-%d")) 





reduced.plant.cover <- left_join(left_join(reduced.plant.cover %>%
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



summary(reduced.plant.cover$peak.window.2mo)
#      N      Y 
# 164558 196218 

# How many additional site-habitats will be missing if we reduce to within the peak window?
sort(unique(paste(reduced.plant.cover$siteID[which(reduced.plant.cover$peak.window.2mo=="N")],
                  reduced.plant.cover$nlcdClass[which(reduced.plant.cover$peak.window.2mo=="N")])))
## Lots!!! 79 site-habitats will be removed!
# if increase the window to 3 months total, still 61 site-habitats removed.
# if 4 month window, 33 site-habitats removed


 
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
  filter(year(collectDate)>2018) %>%
  select(c(3,4,6:16,19,20,23,25,29,32,34:37,41, 45, 57:59)) %>%
  mutate_at(vars(siteID, plotID, plotType, nlcdClass, subplotID,
                 horizon, boutType,sampleTiming,biomassSampleCondition), factor) 

soil.chem <-  soil_periodic$sls_soilChemistry %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  filter(year(collectDate)>2018) %>%
  select(c(5:7,9:10,16:22,26)) %>%
  mutate_at(vars(siteID, plotID, plotType, co2Trapped, cnIsotopeQF), factor)       
  
ls <- unlist(lapply(as.list(soil.chem$sampleID), function(x){length(which(soil.chem$sampleID==x))}))
summary(factor(ls))
# lots of repeat samples
# 4 rows with sample "BLAN_038-M-5.5-36.5-20200714". 2 different measurements of two types of data. 
# analyticalRepNumber goes from 1 to 4 if sample has 4 rows in data.

# I will take means before merging. 
soil.chem.ave <- soil.chem %>%
  group_by(siteID, plotID, plotType, collectDate, sampleID) %>%
  summarise(d15N = mean(d15N[which(cnIsotopeQF=="OK")], na.rm=TRUE), # removes 44 samples where "sample d15N out of qaqc tolerance"
            organicd13C = mean(organicd13C, na.rm=TRUE),
            nitrogenPercent = mean(nitrogenPercent, na.rm = TRUE),
            organicCPercent = mean(organicCPercent, na.rm = TRUE),
            CNratio = mean(CNratio, na.rm = TRUE))
            #cnIsotopeQF = str_flatten(cnIsotopeQF, collapse = " | "))
            

soil.mois <-  soil_periodic$sls_soilMoisture %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  filter(year(collectDate)>2018) %>%
  select(c(3,4,7,8,12,15:19)) %>%
  mutate_at(vars(siteID, plotID, horizon), factor)       
ls <- unlist(lapply(as.list(soil.mois$sampleID), function(x){length(which(soil.mois$sampleID==x))}))
summary(factor(ls))
# don't need to average. only 1 row per sample


soil.pH <-  soil_periodic$sls_soilpH %>%
  mutate(collectDate = floor_date(collectDate, unit="days")) %>%
  filter(year(collectDate)>2018) %>%
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
# or nlcdClass or date. Do that below after merging -------------------------#



dim(soil.periodic.merge) 
# 19117 rows (samples)
length(which(!is.na(soil.periodic.merge$geneticSampleID)))
# 7330 genetic samples

length(which(!is.na(soil.periodic.merge$biomassID)))
#  6374 microbial biomass samples

## See "ExploringSampling.Rmd" to look at how sampling align between soils,
# plant cover, productivity







 
###############################################################################
# Microbial Sequencing Sampling
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


microbe.ITS <- as_tibble(distinct(microbial$mmg_soilMarkerGeneSequencing_ITS[,-1])) %>% # remove duplicated rows before proceeding
  mutate(collectDate=floor_date(collectDate, unit="days"),
         geneticSampleID= unlist(lapply(as.list(dnaSampleID), 
                                 function(x){paste(str_split_1(x,pattern="-GEN")[1], 
                                                   "-GEN", sep ="")}))) %>%    #str_sub(dnaSampleID, 1, -6)) %>% 
  mutate_at(vars(domainID, siteID, plotID, qaqcStatus), factor) %>%
  filter(year(collectDate)>2018,
         qaqcStatus == "Pass", # only keep those that pass qaqc (this removed 150 that failed)
         grepl("COMP", dnaSampleID)==F ) %>% # removed the one dnaSampleID with "COMP" 
  # in it because that means it was pooled and there were other IDs with that plot-date
  select(c("siteID", "plotID", "geneticSampleID", "collectDate", "processedDate", "dnaSampleID",
         "linkerPrimerSequence",	"reverseLinkerPrimerSequence", "illuminaIndex1",	
         "illuminaIndex2", "sampleTotalReadNumber",	"sampleFilteredReadNumber",	
         "minFilteredReadLength",	"maxFilteredReadLength", "sequencerRunID",  
         "qaqcStatus", "publicationDate",	"release"))
# Ylva wanted FilteredReadLength, but they are all NA for our date range

# check that my version of geneticSampleID looks ok
l.gen <- unlist(lapply(as.list(microbe.ITS$geneticSampleID), 
                       function(x){length(str_split_1(x,pattern="-"))}))
summary(factor(l.gen))
#    5    6 
# 1560 5075  
# string of length both 5 and 6 are normal

ls <- unlist(lapply(as.list(microbe.ITS$geneticSampleID), function(x){length(which(microbe.ITS$geneticSampleID==x))}))
summary(factor(ls))
# 1 
# 6635 
# only 1 row per geneticSampleID - that's good.


gs <- soil.periodic.merge$geneticSampleID[which(!is.na(soil.periodic.merge$geneticSampleID))]
summary(microbe.ITS$geneticSampleID %in% gs)
#    Mode    TRUE 
# logical    6635 
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


### Remove crops, pasture, wetlands -----------------------------------------#
soil.periodic.merge <- soil.periodic.merge %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands")

microbe.ITS.metadata <- microbe.ITS.metadata %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands")

  
### Indicate peak biomass date windows for later subsetting -----------------#


soil.periodic.merge <- left_join(left_join(soil.periodic.merge %>%
                                   mutate(year = year(collectDate))
                                 , site.habitat.peakDate),
                                  site.peakDate) %>%
  mutate(peak.date = mdy(paste(if_else(is.na(peakdate), site.peakdate, peakdate), year, sep="-")),
         peak.window.2mo = factor(if_else(abs(difftime(ymd(collectDate), 
                                                   peak.date, units="days"))<32 , "Y", "N")),
         peak.window.3mo = factor(if_else(abs(difftime(ymd(collectDate), 
                                                       peak.date, units="days"))<46 , "Y", "N")),
         peak.window.4mo = factor(if_else(abs(difftime(ymd(collectDate), 
                                                       peak.date, units="days"))<62 , "Y", "N")))
summary(soil.periodic.merge$peak.window.2mo)
# N     Y 
# 11923  4097
summary(soil.periodic.merge$peak.window.3mo)
# N     Y 
# 10337  5683 
summary(soil.periodic.merge$peak.window.4mo)
# N    Y 
# 8803 7217 

# How many site-habitats will be missing if we reduce to within the peak window?
sort(unique(paste(soil.periodic.merge$siteID[which(soil.periodic.merge$peak.window.2mo=="N")],
                  soil.periodic.merge$nlcdClass[which(soil.periodic.merge$peak.window.2mo=="N")])))
##  87 site-habitats will be removed!
# if increase the window to 3 months total, still 79 site-habitats removed.
# 78 with 4 month window




microbe.ITS.metadata <- left_join(microbe.ITS.metadata %>%
                                   mutate(year = year(collectDate))
                                 , site.habitat.peakDate) %>%
  mutate(peak.date = mdy(paste(peakdate, year, sep="-")),
         peak.window.2mo = factor(if_else(abs(difftime(ymd(collectDate), 
                                                       peak.date, units="days"))<32 , "Y", "N")),
         peak.window.3mo = factor(if_else(abs(difftime(ymd(collectDate), 
                                                       peak.date, units="days"))<46 , "Y", "N")),
         peak.window.4mo = factor(if_else(abs(difftime(ymd(collectDate), 
                                                       peak.date, units="days"))<62 , "Y", "N")))
summary(microbe.ITS.metadata$peak.window.2mo)
# removes 3/4 of data

# How many additional site-habitats will be missing if we reduce to within the peak window?
sort(unique(paste(microbe.ITS.metadata$siteID[which(microbe.ITS.metadata$peak.window.2mo=="N")],
                  microbe.ITS.metadata$nlcdClass[which(microbe.ITS.metadata$peak.window.2mo=="N")])))
# removes 66 site-habitats
# 4-month window removes 49 site-habitats


# --------------------------------------------------------------------------- #
# Within the 4-month window, we only want 1 ITS core per site-habitat per year 
# (so don't have to worry about averaging etc)
# Which site-habitats have more than 1 per year? For those, choose the one that is closest to 
# the plant sampling, prioritizing those in the same plot, then those with the 
# most recent cover data, then most recent productivity data.

ITS.site.hab.year.samples <- microbe.ITS.metadata %>%
  filter(peak.window.4mo == "Y") %>%
  group_by(siteID, nlcdClass, year) %>%
  summarise(n.dates = length(unique(collectDate)),
            n.plots = length(unique(plotID)))



ITS.site.hab.year.samples$use.plot <- character(length=dim(ITS.site.hab.year.samples)[1])
ITS.site.hab.year.samples$use.date <- Date(length=dim(ITS.site.hab.year.samples)[1])
for(i in 1:dim(ITS.site.hab.year.samples)[1]){
  
  dat <- microbe.ITS.metadata %>%
    filter(peak.window.4mo == "Y", siteID == ITS.site.hab.year.samples$siteID[i], 
           nlcdClass == ITS.site.hab.year.samples$nlcdClass[i], 
           year == ITS.site.hab.year.samples$year[i])
  
  plot.dates <- distinct(dat[,c("plotID","collectDate")])
  
  if(dim(plot.dates)[1] ==1 ){ # if there is only 1 plot date, use it
    ITS.site.hab.year.samples$use.plot[i] <- as.character(plot.dates$plotID[1])
    ITS.site.hab.year.samples$use.date[i] <- plot.dates$collectDate[1]
  } else { # if there is more than 1, then figure out which one to use:
  cover.plotdates <- reduced.plant.cover %>%
    filter(peak.window.4mo == "Y", siteID == ITS.site.hab.year.samples$siteID[i], 
           nlcdClass == ITS.site.hab.year.samples$nlcdClass[i], 
           year == ITS.site.hab.year.samples$year[i]) %>%
    select(c("plotID","observation_datetime")) %>%
    distinct()
  
  biomass.plotdates <- productivity.persample %>%
    filter(peak.window.4mo == "Y", siteID == ITS.site.hab.year.samples$siteID[i], 
           nlcdClass == ITS.site.hab.year.samples$nlcdClass[i], 
           year == ITS.site.hab.year.samples$year[i]) %>%
    select(c("plotID","collectDate")) %>%
    distinct()
  
  # for each plot, was it also sampled for cover or productivity?
  plot.dates$plotin.cover <- plot.dates$plotID %in% unique(cover.plotdates$plotID)
  plot.dates$plotin.prod <- plot.dates$plotID %in% unique(biomass.plotdates$plotID)
  
  both.plot.dates <- plot.dates%>%
    filter(plotin.cover == T, plotin.prod ==T)
  
  if(dim(both.plot.dates)[1]==1){ # if there is only 1 plot that is in both datasets, use that:
    ITS.site.hab.year.samples$use.plot[i] <- as.character(both.plot.dates$plotID[1])
    ITS.site.hab.year.samples$use.date[i] <- both.plot.dates$collectDate[1]
  }
  if(dim(both.plot.dates)[1] > 1){ # if there is more than 1 plot in both datasets, use the one that is closest in sampling date to plant cover data 
    both.plot.dates$since.cover <- rep(NA, dim(both.plot.dates)[1])
    for(j in 1: dim(both.plot.dates)[1]){
      cd <- cover.plotdates[which(as.character(cover.plotdates$plotID)==as.character(both.plot.dates$plotID[j])),] 
      both.plot.dates$since.cover[j] <-  min(abs(diff.Date(c(cd$observation_datetime, both.plot.dates$collectDate[j]))))
    }
      
    ITS.site.hab.year.samples$use.plot[i] <- as.character(both.plot.dates$plotID[which.min(both.plot.dates$since.cover)])
    ITS.site.hab.year.samples$use.date[i] <- both.plot.dates$collectDate[which.min(both.plot.dates$since.cover)]
  }
  if(dim(both.plot.dates)[1] == 0){ #if there is not a plot in both datasets, but in 1 of the datasets... 
    if(length(which(plot.dates$plotin.cover==TRUE))>0){
      d <- plot.dates[which(plot.dates$plotin.cover==T),]
      d$since.cover <- rep(NA, dim(d)[1])
      for(k in 1: dim(d)[1]){
        cd <- cover.plotdates[which(as.character(cover.plotdates$plotID)==as.character(d$plotID[k])),] 
        d$since.cover[k] <-  min(abs(diff.Date(c(cd$observation_datetime, d$collectDate[k]))))
      }
      ITS.site.hab.year.samples$use.plot[i] <- as.character(d$plotID[which.min(d$since.cover)])
      ITS.site.hab.year.samples$use.date[i] <- d$collectDate[which.min(d$since.cover)]
    } else{ 
      if(length(which(plot.dates$plotin.prod==TRUE))>0){
        d <- plot.dates[which(plot.dates$plotin.prod==T),]
        d$since.prod <- rep(NA, dim(d)[1])
        for(k in 1: dim(d)[1]){
          bd <- biomass.plotdates[which(as.character(biomass.plotdates$plotID)==as.character(d$plotID[k])),] 
          d$since.prod[k] <-  min(abs(diff.Date(c(bd$collectDate, d$collectDate[k]))))
        }
        ITS.site.hab.year.samples$use.plot[i] <- as.character(d$plotID[which.min(d$since.prod)])
        ITS.site.hab.year.samples$use.date[i] <- d$collectDate[which.min(d$since.prod)]
      } else{ # no overlap in plots, so just choose the date that is closest to a cover survey (regardless of plot) 
      if(dim(cover.plotdates)[1]>0){
          d <- plot.dates
          d$since.cover <- rep(NA, dim(d)[1])
          for(k in 1: dim(d)[1]){
            dates <- cover.plotdates$observation_datetime
            d$since.cover[k] <-  min(abs(dates - d$collectDate[k])) # this doesn't work if more than 1 date so need to fix!
          } 
        ITS.site.hab.year.samples$use.plot[i] <- as.character(d$plotID[which.min(d$since.cover)])
        ITS.site.hab.year.samples$use.date[i] <- d$collectDate[which.min(d$since.cover)]
        
      } else{# if there is no cover data use biomass date
        if(dim(biomass.plotdates)[1]>0){
          d <- plot.dates
          d$since.prod <- rep(NA, dim(d)[1])
          for(k in 1: dim(d)[1]){
            dates <- biomass.plotdates$collectDate
            d$since.prod[k] <-  min(abs(dates - d$collectDate[k])) # this doesn't work if more than 1 date so need to fix!
          } 
          ITS.site.hab.year.samples$use.plot[i] <- as.character(d$plotID[which.min(d$since.prod)])
          ITS.site.hab.year.samples$use.date[i] <- d$collectDate[which.min(d$since.prod)]
         } else {
        ITS.site.hab.year.samples$use.plot[i] <- NA
        ITS.site.hab.year.samples$use.date[i] <- NA
        }
        
      }
    } 
    }
    
  }
  }
}

# still 3 NAs - a site-habitat-year that had ITS data but no plant data


# Pull out ITS sampleIDs for these plot-dates



# --------------------------------------------------------------------------- #

names(ITS.site.hab.year.samples)[6:7] <- c("plotID", "collectDate")
x <- full_join(drop_na(ITS.site.hab.year.samples),microbe.ITS.metadata)
length(which(!is.na(x$n.dates)))
reduced.microbe.ITS.metadata <- x[which(!is.na(x$n.dates)),-(4:5)]
# this includes all mineral and organic horizons - not sure if we want that.

write_csv(reduced.microbe.ITS.metadata, "ReducedMicrobeITSmetadata.csv")

###############################################################################
# Microbial Biomass Sampling
###############################################################################

load("MicrobialBiomassRaw.RData")

names(microbial.biomass.raw)
names(microbial.biomass.raw$sme_microbialBiomass)
names(microbial.biomass.raw$sme_scaledMicrobialBiomass)

unscaledIDs <- unique(microbial.biomass.raw$sme_microbialBiomass$biomassID)
scaledIDs <- unique(microbial.biomass.raw$sme_scaledMicrobialBiomass$biomassID)
#biomassID is just sampleID with "-BM" on end

summary(unscaledIDs %in% scaledIDs)
summary(scaledIDs %in% unscaledIDs )
# quite a few don't overlap

summary(scaledIDs %in% soil.periodic.merge$biomassID)
#    Mode   FALSE    TRUE 
# logical     715    6347   
summary(unscaledIDs %in% soil.periodic.merge$biomassID)
#    Mode   FALSE    TRUE 
# logical    1828    4021 


dim(microbial.biomass.raw$sme_scaledMicrobialBiomass)
length(scaledIDs)
# 5 samples were run twice - need to take mean or first value (they're very similar)
ls <- unlist(lapply(as.list(microbial.biomass.raw$sme_scaledMicrobialBiomass$biomassID), 
                    function(x){length(which(microbial.biomass.raw$sme_scaledMicrobialBiomass$biomassID==x))}))
summary(factor(ls))
microbial.biomass.raw$sme_scaledMicrobialBiomass$biomassID[which(ls==2)]

cols <- names(microbial.biomass.raw$sme_scaledMicrobialBiomass)[c(13,24:84,86)]
microbial.biomass.scaled <- microbial.biomass.raw$sme_scaledMicrobialBiomass %>%
  mutate_at(vars(domainID, siteID, plotID), factor) %>%
  filter(year(collectDate)>2018) %>% 
  # nlcdClass isn't in this dataset - it will remove the extraneous ones when I merge with periodic soil data
  #filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
  #      nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands")
  group_by(siteID, plotID, sampleID, biomassID) %>%
  summarise(across(all_of(cols), ~ mean(.x, na.rm=TRUE)))
summary(microbial.biomass.scaled)
summary(unique(microbial.biomass.scaled$biomassID) %in% soil.periodic.merge$biomassID)
# totalLipidScaledConcentration

soil.periodic.merge <- left_join(soil.periodic.merge, microbial.biomass.scaled)
  



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
                            select(c("siteID",	"plotID","collectDate",	"horizonID",	"horizonName",
                                     "biogeoTopDepth",	"biogeoBottomDepth",	"biogeoCenterDepth",
                                     "gypsumConc",	"caco3MethodPub",	"caco3Method",	"caco3Conc",
                                     "caNh4d",	"kNh4d",	"mgNh4d",	"naNh4d",	"cecdNh4",	"alSatCecd33",
                                     "baseSumCecd10",	"bsesatCecd10", "alKcl",	"feKcl",	"mnKcl",
                                     "phCacl2",	"phH2o", "carbonTot",	"nitrogenTot",	"ctonRatio",	
                                     "estimatedOC",	"sulfurTot", "alOxalate",	"feOxalate",	"mnOxalate",
                                     "pOxalate",	"siOxalate", "acidity", "MehlichIIITotP","Bray1PExtractable", 
                                     "OlsenPExtractable"))
                              

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
                            select(c("siteID", "plotID", "collectDate", "horizonID", "horizonName",
                                     "bulkDensTopDepth",	"bulkDensBottomDepth","bulkDensCenterDepth",	
                                     "bulkDensOvenDry",	"bulkDensFieldMoist"))
dim(soil.initial.bulkdensity)
length(unique(soil.initial.bulkdensity$horizonID))
# a few horizonIDs repeated
ls <- unlist(lapply(as.list(soil.initial.bulkdensity$horizonID), 
                     function(x){length(which(soil.initial.bulkdensity$horizonID==x))}))
summary(factor(ls))
as.data.frame(soil.initial.bulkdensity[which(ls>1),])
# 3 are repeated because diff analyses with diff SampleType: clod vs compliant cavity
# take the mean (but there's only 1 entry per variable type)

soil.initial.bulkdensity <- soil.initial.bulkdensity %>%
  group_by(siteID,   plotID, collectDate, horizonID, horizonName) %>%
  summarise(across(1:5, \(x) mean(x, na.rm=TRUE), .names="{.col}"))
# returns NaN instead of NA for missing values now

# merge
soil.initial <- left_join(soil.initial,soil.initial.bulkdensity)




############ now add particlesize

soil.initial.particlesize <- as.tibble(soil$spc_particlesize) %>%
                              mutate(collectDate=floor_date(collectDate, unit="days")) %>% 
                              mutate_at(vars(domainID, siteID, plotID, horizonID, horizonName), factor) %>%
                              select(c("siteID",	"plotID", "collectDate",  "horizonID",	"horizonName",
                                       "biogeoTopDepth",	"biogeoBottomDepth","biogeoCenterDepth",
                                       "sandTotal",	"siltTotal",	"clayTotal"))

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

# and not same plots either

#584 plots with initial soil samples
#343 plots with microbe ITS data
summary(unique(microbe.ITS.metadata$plotID) %in% unique(soil.initial$plotID))
# 208 of which are the same plots in the initial soil sample 

# e.g., ABBY_070 has 8 microbial samples


summary(microbe.ITS.metadata[,c("plotID","collectDate","qaqcStatus","sampleTiming","sampleID","horizon","litterDepth","sampleTopDepth","sampleBottomDepth")])

soil.horizon.lineup <- data.frame(siteID = microbe.ITS.metadata$siteID,
                                  nlcdClass = microbe.ITS.metadata$nlcdClass,
                                  microbe.sampleID = microbe.ITS.metadata$sampleID, 
                                  microbe.horizon = microbe.ITS.metadata$horizon,
                                  n.initial.horizons.any = rep(NA, dim(microbe.ITS.metadata)[1]),
                                  n.initial.horizons.total = rep(NA, dim(microbe.ITS.metadata)[1]))
for(i in 1:dim(microbe.ITS.metadata)[1]){
  mt <- microbe.ITS.metadata$sampleTopDepth[i]
  mb <- microbe.ITS.metadata$sampleBottomDepth[i]
  
  mt1 <- as.Date(mt*10, origin = "2020-01-01")
  mb1 <- as.Date(mb*10, origin = "2020-01-01")
  int.m <- interval(start = mt1, end = mb1 )
  
  # which soil.initial horizon IDs are there for that site-habitat combo?
  sdat <- soil.initial %>%
    filter(siteID == microbe.ITS.metadata$siteID[i],
           nlcdClass == microbe.ITS.metadata$nlcdClass[i]) %>%
    select(horizonID, horizonTopDepth, horizonBottomDepth)
  
  int.s <- interval(start = as.Date(sdat$horizonTopDepth*10, origin = "2020-01-01"),
                    end = as.Date(sdat$horizonBottomDepth*10, origin = "2020-01-01"))
  
  
  # int_overlaps(int.m, int.s) # any overlap
  # int.m %within% int.s # total overlap 
  hor.any <- sdat$horizonID[which(int_overlaps(int.m, int.s))]
  soil.horizon.lineup$n.initial.horizons.any[i] <- length(hor.any)
  hor.total <- sdat$horizonID[which(int.m %within% int.s)]
  soil.horizon.lineup$n.initial.horizons.total[i] <- length(hor.total)
}

# if want total overlap - lots of 0s
summary(soil.horizon.lineup$n.initial.horizons.total)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   0.000   1.244   1.000  18.000 

# if want any overlap:
summary(soil.horizon.lineup$n.initial.horizons.any)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   11.00   18.00   20.47   29.00   49.00 
# 45 out of 6085 don't overlap at all.
# For each periodic soil sample could probably average 
# initial soil variables across the ones that overlap any




# need to create a list of horizonIDs for each row of microbe.ITS.metadata

soilhorizon.periodictoinitial <- list()
for(i in 1:dim(microbe.ITS.metadata)[1]){
  mt <- microbe.ITS.metadata$sampleTopDepth[i]
  mb <- microbe.ITS.metadata$sampleBottomDepth[i]
  
  mt1 <- as.Date(mt*10, origin = "2020-01-01")
  mb1 <- as.Date(mb*10, origin = "2020-01-01")
  int.m <- interval(start = mt1, end = mb1 )
  
  # which soil.initial horizon IDs are there for that site-habitat combo?
  sdat <- soil.initial %>%
    filter(siteID == microbe.ITS.metadata$siteID[i],
           nlcdClass == microbe.ITS.metadata$nlcdClass[i]) %>%
    select(horizonID, horizonTopDepth, horizonBottomDepth)
  
  int.s <- interval(start = as.Date(sdat$horizonTopDepth*10, origin = "2020-01-01"),
                    end = as.Date(sdat$horizonBottomDepth*10, origin = "2020-01-01"))
  
  
  # int_overlaps(int.m, int.s) # which horizons overlap at all (any overlap not total overlap)
  soilhorizon.periodictoinitial[[i]] <- sdat$horizonID[which(int_overlaps(int.m, int.s))]

}



####### average initial soil variables for any horizon in that site-habitat
# that overlaps with microbial sampling --------------------------------------#
scols <- c(33,36:71) # columns in the soil.initial data to average
microbe.ITS.metadata.initial.soils <- data.frame(microbe.ITS.metadata[,c(1:4,36,37)], # columns in microbe.ITS.metadata to match by (top and bottom depth)
                  soil.initial[dim(microbe.ITS.metadata)[1],scols]) # add the column names and NA
mcols <- 7:43 # which columns to put the new averages in
for(r in 1:dim(dat)[1]){
  hdat <- as.data.frame(soil.initial[match(soilhorizon.periodictoinitial[[r]], soil.initial$horizonID),])
  for(c in 1:length(scols)){
    microbe.ITS.metadata.initial.soils[r,mcols[c]] <- mean(hdat[,scols[c]], na.rm=T)
  }
}
summary(microbe.ITS.metadata.initial.soils)




###############################################################################
# Root biomass
###############################################################################

# depthIncrementID (from perDepth and perRoot) seems to match the first part of cnSampleID (CN)
# both start with pitProfileID (from perdepth and perpit)
# and cnSampleID from CN seems to match sampleID from perRootSample
# namedLocation from CN seems to match pitNamedLocation from all other files - change name

root.CN <- root.biomass.raw$mpr_carbonNitrogen %>%
  select("namedLocation","siteID", "collectDate", "cnSampleID", "sampleType", 
  "d15N",	"d13C",	"nitrogenPercent",	"carbonPercent",	"CNratio", "analyticalRepNumber") %>%
  mutate(pitProfileID = unlist(lapply(str_split(cnSampleID, "[.]",3), 
                                      function(y){paste(y[1],y[2],sep=".")}))) %>%
  rename(pitNamedLocation = namedLocation,
         sampleID =cnSampleID)
 

root.depth <- root.biomass.raw$mpr_perdepthincrement %>%
  select("pitNamedLocation","pitProfileID",	"depthIncrementID",	"topDepth",	
         "bottomDepth",	"depthIncrementVolume")  


root.pit <- root.biomass.raw$mpr_perpitprofile %>%
  select("pitNamedLocation", "siteID",	"decimalLatitude",	"decimalLongitude", 
         "elevation", "nlcdClass",   "pitID",	"pitProfileID",	"rootStatus",	
         "sizeCategory",	"maxProfileDepth",	"totalRootBiomass",	"depth100RootBiomass") %>%
  mutate(across(c(nlcdClass, rootStatus, sizeCategory), factor))

root.perroot <- root.biomass.raw$mpr_perrootsample %>%
  select("pitNamedLocation", "depthIncrementID",	"sampleID", "rootStatus",	
         "sizeCategory",	"rootDryMass", "incrementRootBiomass",	"incrementRootDensity") %>%
  mutate(pitProfileID = unlist(lapply(str_split(sampleID, "[.]",3), 
                                      function(y){paste(y[1],y[2],sep=".")}))) %>%
  mutate(across(c(rootStatus, sizeCategory), factor))

root_merge <- full_join(root.depth, root.perroot)
root_merge <- full_join(root_merge, root.CN)
root_merge <- full_join(root_merge, root.pit)

summary(root_merge)

root.site.habitat.summary <- root_merge %>%
  group_by(siteID, nlcdClass) %>%
  summarise(n_pits = length(unique(pitNamedLocation)),
            n_CN_samples = length(unique(sampleID[!is.na(CNratio)])),
            n_rootDryMass_samples = length(unique(sampleID[!is.na(rootDryMass)])))

# looks like one pit per site (not site-habitat)
# but there are some NAs. need to look at them... 

root_merge[which(is.na(root_merge$nlcdClass)),]

write_csv(root.site.habitat.summary, file="RootSampling.csv")


###############################################################################
# pull in the neonDivData on beetles and birds and smammals 
# and add peak productivity windows
###############################################################################


#### Beetles -----------------------------------------------------------------#
beetlesDD <- neonDivData::data_beetle

beetles <- beetlesDD %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands",
         year(observation_datetime)>2018) 

beetles <- left_join(left_join(beetles %>%
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



#### Small mammals -----------------------------------------------------------#

smammalsDD <- neonDivData::data_small_mammal

smammals <- smammalsDD %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands",
         year(observation_datetime)>2018) 

smammals <- left_join(left_join(smammals %>%
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



#### Birds  -----------------------------------------------------------------#
birdsDD <- neonDivData::data_bird

birds <- birdsDD %>%
  filter(nlcdClass != "cultivatedCrops" & nlcdClass !="pastureHay" & 
           nlcdClass !="emergentHerbaceousWetlands" & nlcdClass !="woodyWetlands",
         year(observation_datetime)>2018) 

birds <- left_join(left_join(birds %>%
                                 mutate(year = year(observation_datetime))
                               , site.habitat.peakDate),
                     site.peakDate) %>%
  mutate(peak.date = mdy(paste(if_else(is.na(peakdate), site.peakdate, peakdate), year, sep="-")), #if it's present, use the site-habitat peak date, otherwise use the site peak date
         # is this within 1 month on either side of peak date (so 2 month window?)
         peak.window.2mo = factor(if_else(abs(difftime(observation_datetime, 
                                                       peak.date, units="days"))<32 , "Y", "N")),
         # is this within 1.5 month on either side of peak date (so 3 month window?)
         peak.window.3mo = factor(if_else(abs(difftime(observation_datetime, 
                                                       peak.date, units="days"))<46 , "Y", "N")),
         peak.window.4mo = factor(if_else(abs(difftime(observation_datetime, 
                                                       peak.date, units="days"))<62 , "Y", "N")))




###############################################################################
 
save(plant.cover, reduced.plant.cover, productivity, productivity.persample,
     soil.periodic.merge, soil.initial, microbe.ITS.metadata, beetles, smammals, birds,
     file="ReducedMergedData.RData")

# All datasets have removed nlcdClass crops, pasture, wetlands.
# All but the soil.initial only include 2019 and onwards.

# not filtered by time of year but has column for peak window (2,3,4 months)
# based on plant biomass data
