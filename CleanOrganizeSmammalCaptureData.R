###############################################################################
## Cleaning and Organizing Raw SMammal Data for abundance estimation
###############################################################################

library(tidyverse)

source("FunctionsForNEONdata.R")

load("rawSmammalCaptureData.RData")
names(smammal.capturedata)

#### raw data with all captures from all sites and dates that weren't provisional as of 7/31/24
#### includes lines for empty traps 

NEONsmammalcaptures <- as_tibble(smammal.capturedata$mam_pertrapnight)
names(NEONsmammalcaptures)
head(NEONsmammalcaptures) 

trapping.info <-  trap.status.function(captures = NEONsmammalcaptures)
# traps_available include all traps that were set (even if disturbed, sign, etc)

hist(trapping.info$traps_available, breaks=20) 
hist(trapping.info$traps_available[which(trapping.info$traps_available<40)], breaks=20)



## Clean and organize capture data ------------------------------------------##

#collectDate and endDate are the the same

smammal.captures <- NEONsmammalcaptures %>%
  filter(!is.na(taxonID))
# remove lines with no taxon ID (many of the lines are empty traps)

# keep pertinent columns for captures
smammal.captures <- smammal.captures %>%
  select(5:7,18,20,22:25,29:48,50,53,56,58) %>%
  mutate_at(vars(siteID, plotID, tagID, taxonID, nativeStatusCode,
                 recapture, fate, replacedTag, lifeStage, testes, nipples,
                 pregnancyStatus, vagina, sex, larvalTicksAttached,
                 nymphalTicksAttached, adultTicksAttached, tickNumber), factor)

summary(smammal.captures)

length(which(!is.na(smammal.captures$bloodSampleID)))
# 37231 out of 146858 captures had blood sample IDs
# many were tested for pathogen (SNV through 2019 and after that tick borne pathogens)
# eventually I want to merge those data

smammal.species.info <- smammal.captures %>%
  group_by(scientificName, taxonRank, taxonID) %>%
  summarise(n=n())
print(smammal.species.info,n=200)
# lots IDed to just genus or to a group (several species)
# a few subspecies
# summary(factor(smammal.captures$taxonRank))
# class       family        genus        order      species speciesGroup   subspecies 
#   357           25         1912           46       141266         3248            4 

write_csv(smammal.species.info, file="NEONsmammalSpecies.csv")




## Summaries of Plots including trapnights --------------------------------- ##

smammal.plot.data <- NEONsmammalcaptures %>%
  group_by(siteID, plotID, plotType, nlcdClass) %>% # there area few with slightly diff elevation, long, lat so average
  summarise(latitude = mean(decimalLatitude),
            longitude = mean(decimalLongitude),
            elevation = mean(elevation),
            n.dates = length(unique(floor_date(collectDate, unit="days"))),
            first.date = min(collectDate),
            last.date = max(collectDate),
            total.richness = length(na.omit(unique(taxonID))))
trapstat <- trapping.info %>%
  group_by(siteID, plotID) %>%
  summarise(trapnights = sum(traps_available))
smammal.plot.data <- left_join(smammal.plot.data,
                               trapstat)

dat2 <- NEONsmammalcaptures %>%
  group_by(siteID,plotID) %>%
  summarise(uniqueIDs = length(na.omit(unique(tagID))), # a few plots had 0 captures, so need to make those 0 instead of NA
            uniquePM = length(unique(tagID[which(taxonID=="PEMA")])))
smammal.plot.data <- left_join(smammal.plot.data, dat2)
smammal.plot.data <- smammal.plot.data %>%
  mutate(IDsper100trapnights = uniqueIDs/(trapnights/100),
         PMper100trapnights = uniquePM/(trapnights/100))
# this is total unique individuals (not abundance) - they could have been present for more than 1 night

hist(smammal.plot.data$total.richness, breaks=20)

par(mfrow=c(2,1))
hist(smammal.plot.data$n.dates, breaks=30, main="number of dates trapped per plot")
hist(smammal.plot.data$trapnights, breaks=30, main="number of trap-nights per plot")

par(mfrow=c(1,1))


ggplot(smammal.plot.data,                               
       aes(x = plotID,
           y = PMper100trapnights,
           col = siteID)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position="top") +
  geom_point() +
  guides(col=guide_legend(ncol=12)) # makes the legend have more columns

hist(smammal.plot.data$PMper100trapnights, breaks = 50)
length(which(smammal.plot.data$PMper100trapnights==0))/ dim(smammal.plot.data)[1]
# a third have 0 deer mice

gsum <- smammal.plot.data %>% # first get sample sizes so can add
  group_by(nlcdClass) %>%
  summarise(n.plots = length(unique(plotID)))
ggplot(smammal.plot.data, aes(x = nlcdClass, y = IDsper100trapnights, col=nlcdClass)) +
  geom_violin(aes(fill = nlcdClass), show.legend = F) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_text(data =gsum, aes(nlcdClass, Inf, label=n.plots), vjust = 1, show.legend = FALSE) + # add sample sizes at top
  ylab("IDs per 100 trap nights") +
  ggtitle("Individuals per habitat")  


gsum2 <- smammal.plot.data %>% # first get sample sizes so can add
  group_by(siteID, nlcdClass) %>%
  summarise(n.plots = length(unique(plotID)))
ggplot(smammal.plot.data, aes(x = paste(siteID,nlcdClass), y = IDsper100trapnights, col=siteID)) +
  geom_violin(aes(fill = nlcdClass), show.legend = F) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, size=11)) +
  geom_text(data =gsum2, aes(paste(siteID,nlcdClass), Inf, label=n.plots), vjust = 1, show.legend = FALSE) + # add sample sizes at top
  ylab("IDs per 100 trap nights") +
  xlab("") + 
  ggtitle("Individuals per site-habitat")  

ggplot(smammal.plot.data, aes(x = paste(siteID,nlcdClass), y = PMper100trapnights, col=siteID)) +
  geom_violin(aes(fill = nlcdClass), show.legend = F) +
  theme_classic(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, size=11)) +
  #geom_text(data =gsum, aes(nlcdClass, Inf, label=n.sites), vjust = 1, show.legend = FALSE) + # add sample sizes at top
  ylab("unique deer mice per 100 trap nights") +
  xlab("") + 
  ggtitle("Deer mice per site-habitat")  



smammal.sitehabitat.data <- smammal.plot.data %>%
  group_by(siteID, nlcdClass) %>%
  summarise(n.plots = length(unique(plotID)),
            total.trapnights = sum(trapnights),
            sum.uniqueIDs = sum(uniqueIDs),
            mean.IDsper100trapnights = mean(IDsper100trapnights))


# save(smammal.captures, smammal.plot.data, 
#      trapping.info ,file="smammalcaptures.RData")


###############################################################################
# Clean data for abundance estimation
###############################################################################

# remove those without tag numbers
smammal.captures.clean <- smammal.captures %>%
  filter(!is.na(tagID))

# add genus and species columns to help with grouping for estimation
smammal.captures.clean$genus <- str_split_i(smammal.captures.clean$scientificName," ",1)
smammal.captures.clean$species <- paste(smammal.captures.clean$genus , 
                                        str_split_i(smammal.captures.clean$scientificName," ",2))

# remove some 'species' that aren't helpful
sp.rm <- c("CRSP", # crecitidae
           "DIVI", # opossum 6 captures
           "LEAM", # rabbit 5 captures
           "HESP", # Heteromyidae family
           "MRSP", # Muridae family
           "OTHE", # Mammalia
           "ROSP", # Rodentia
           "SRSP" # Soricidae
)
smammal.captures.clean <- smammal.captures.clean %>%
  filter( !taxonID %in% sp.rm ) 

# collapse a couple subspecies to species
smammal.captures.clean$species[which(smammal.captures.clean$taxonID=="SIHE")] <-"Sigmodon hispidus"
smammal.captures.clean$taxonID[which(smammal.captures.clean$taxonID=="SIHE")] <-"SIHI"
smammal.captures.clean$species[which(smammal.captures.clean$taxonID=="ICTM")] <-"Ictidomys tridecemlineatus"
smammal.captures.clean$taxonID[which(smammal.captures.clean$taxonID=="ICTM")] <-"ICTR"





#### Quality Control ---------------------------------------------------------#

# Are there duplicates of same tag/plot on same date? ------------------- ##
# Yes, lots:
cols <- names(smammal.captures.clean)
cols <- cols[-match(c("collectDate", "tagID", "plotID"),cols)]
dupes.date <- smammal.captures.clean %>%
  group_by(collectDate, tagID, plotID) %>%
  summarise(n = n(),
            across(all_of(cols), 
                   ~ length(unique(.x)), .names = "n_{.col}")) %>% 
  filter(n > 1)
dim(dupes.date)
# 78  36
# lots have multiple columns with n=2. trapCoordinate is most often, but
# taxonID is the biggest problem 
as.data.frame(dupes.date[which(dupes.date$n_taxonID==2),])
#28 rows where same tag number on same date in same plot is 2 different species

as.data.frame(
  smammal.captures.clean %>%
  filter(collectDate==ymd("2014-07-25") & tagID == "NEON.MAM.D01.R0433" &  plotID=="HARV_001")
)
# same tag number: male PEMA and female PELE,  same plot and date. Ugh.



# Is the same tag IDed as different species or sexes? ------------------- ##
# Yes, lots:
dupes.taxsex <- smammal.captures.clean %>%
  group_by(tagID, plotID) %>%
  summarise(across(c(taxonID, sex), 
                   ~ length(unique(.x)), .names = "n_{.col}")) %>%
  filter(n_taxonID>1 | n_sex>1)

dim(dupes.taxsex)
# dim(dupes.taxsex)
# [1] 5686    4
# Yikes. 5686 have multiple species or sexes
length(which(dupes.taxsex$n_taxonID>1))
# 3206 tagIDs have multiple species assigned (within the same plot)
# I guess I'll make unique ID be plot_taxon_tag 
# (it's possible some individuals move between plots within the same site,
# but here I will consider them separate.)
# right now I don't care about sex, so will leave those to deal with later.

smammal.captures.clean <- smammal.captures.clean %>%
  mutate(plot_taxon_tag = paste(plotID, taxonID, tagID, sep="_"))



##create a list of dataframes of sessions for each plot
trapping.session.info <- NEON.session.function(trapping.info, 
                                               int.break=10) #if interval is greater than 10 it's a different session

dat <- trapping.session.info %>%
  group_by(plotID, prim.session) %>%
  summarise(trapnights = sum(traps_available))
summary(dat)
hist(dat$trapnights, breaks=30)
hist(dat$trapnights[which(dat$trapnights<100)], breaks=30)


### Add hantavirus data ----------------------------------------------------- ##

load("HantavirusTestResults.RData")

summary(hanta.bloodtesting$bloodSampleID %in% smammal.captures.clean$bloodSampleID)
#    Mode   FALSE    TRUE 
# logical       3   15855 


length(unique(hanta.bloodtesting$bloodSampleID[which(!is.na(hanta.bloodtesting$testResult))]))
#[1] 14002
dim(hanta.bloodtesting[which(!is.na(hanta.bloodtesting$testResult)),])
# 14044    24
# 42 samples have multiple results (that aren't NA). Are they always the same? pos/neg?
dat <- hanta.bloodtesting %>%
  filter(!is.na(testResult)) %>%
  group_by(bloodSampleID) %>%
  summarise(n = n(),
            n_dates = length(unique(collectDate)),
            n_plots = length(unique(plotID)),
            n_sites = length(unique(siteID)),
            n_results = length(unique(testResult)),
            test = factor(str_flatten(sort(testResult), collapse = " "))) %>%
  filter(n>1)
summary(dat)
dat$bloodSampleID[which(dat$n>1)]
# all of them are negative both times.
# remove duplicates then merge with capture data
dat <- hanta.bloodtesting %>%
  filter(!is.na(testResult)) %>%
  mutate(hantaResult = testResult) %>%
  select(siteID, plotID, collectDate, bloodSampleID, hantaResult) 
dat <- distinct(dat)
smammal.captures.clean <- left_join(smammal.captures.clean, dat)

hanta.species.summary <- smammal.captures.clean %>%
  filter(!is.na(hantaResult)) %>%
  group_by(species) %>%
  summarise(n_pos = length(which(hantaResult=="Positive")),
            n_test = length(which(!is.na(hantaResult))),
            n_prev = n_pos/n_test)


### Add tick-borne pathogen data ------------------------------------------- ##

load("TickBornePathogenTestResults.RData")

### lots of blood samples don't line up to smammal data 

summary(TBD.bloodtesting$sampleID %in% smammal.captures.clean$bloodSampleID)
#    Mode   FALSE    TRUE 
# logical   47697   37582

summary(TBD.bloodtesting$sampleID %in% NEONsmammalcaptures$bloodSampleID)
#    Mode   FALSE    TRUE 
# logical   47697   37582 
# the problem is not that I've removed individuals in the 'clean data'

summary(TBD.bloodtesting$sampleCode %in% NEONsmammalcaptures$bloodSampleBarcode)
#    Mode   FALSE    TRUE 
# logical   47653   37626
# doesn't matter whether I use sampleID or sampleCode, sample code catches a few more

# the tick data contains 1 year more than the smammal data I'm working with
# but still there are lot not identified.
summary(TBD.bloodtesting$collectDate[which(TBD.bloodtesting$sampleID %in% smammal.captures.clean$bloodSampleID
  ==FALSE)])
max(smammal.captures.clean$collectDate)
# "2022-12-15 GMT"

ds <- TBD.bloodtesting$collectDate[which(TBD.bloodtesting$sampleID %in% smammal.captures.clean$bloodSampleID
                                         ==FALSE)]
length(which(ds > max(smammal.captures.clean$collectDate))) # how many after "2022-12-15 GMT"
# 24342 of the pathogen samples were after my smammal data ends
# late leaves 23355 that should be in the dataset.

# ".E" is supposed to be earSampleID and ".B" is bloodSampleID but in the 
# smammal data, bloodSampleID sometimes has .E on the end

# make another column that ignores the end. Indicates that animal that day:
TBD.bloodtesting$animalSampleID <- substr(TBD.bloodtesting$sampleID,1,nchar(TBD.bloodtesting$sampleID)-2)

# for smammal.captures, use blood and then if missing use ear
smammal.captures.clean$animalSampleID <- substr(smammal.captures.clean$bloodSampleID,
                                                1,nchar(smammal.captures.clean$bloodSampleID)-2)
inds <- which(is.na(smammal.captures.clean$animalSampleID))
smammal.captures.clean$animalSampleID[inds] <- substr(smammal.captures.clean$earSampleID[inds] ,
                                                1,nchar(smammal.captures.clean$earSampleID[inds] )-2)

summary(TBD.bloodtesting$animalSampleID %in% smammal.captures.clean$animalSampleID)
#    Mode   FALSE    TRUE 
# logical   24414   60865 
ds <- TBD.bloodtesting$collectDate[which(TBD.bloodtesting$animalSampleID %in% 
                                           smammal.captures.clean$animalSampleID
                                         ==FALSE)]
length(which(ds < max(smammal.captures.clean$collectDate))) # how many are missing before "2022-12-15 GMT"
# much better, still about 72 missing animals compared to TBD test results

### merge the tick borne pathogen data with the smammal capture data

# first want each animal capture to be a row and results of various tests as columns
TBD.bySample <- TBD.bloodtesting %>%
  mutate(testResult2 = ifelse(testResult=="Negative", 0, ifelse(testResult=="Positive", 1, NA))) %>%
  pivot_wider(id_cols = c(siteID, plotID, animalSampleID, collectDate),
              names_from = testPathogenName,
              values_from = testResult2,
              values_fn = ~ sum(.x, na.rm=TRUE)) # sometimes both blood and ear were tested so take sum

summary(TBD.bySample)
# a bunch of the pathogens weren't found at all in these data, so remove those columns
TBD.bySample <- TBD.bySample[ ,c(1:4, which(apply(TBD.bySample[,-(1:4)],2,sum, na.rm=T)!=0)+4)]
# 0 means all samples were negative, 1 means 1 sample (either blood or ear was positive), 2 means both positive
# replace all 2 with 1, so 1 will mean any sample was positive
TBD.bySample[,5:dim(TBD.bySample)[2]] <- replace(TBD.bySample[,5:dim(TBD.bySample)[2]], TBD.bySample[,5:dim(TBD.bySample)[2]]==2, 1)


smammal.captures.clean <- left_join(smammal.captures.clean, TBD.bySample)



################################################################################
### For capture rate estimation - Remove sessions that only have 1 night ---- ## 
### (won't help with estimation)
################################################################################

reduced.trapping.session.info <- left_join(
  trapping.session.info, 
  trapping.session.info %>%
    group_by(siteID, plotID, prim.session) %>%
    summarise(keep = ifelse(any(Sec>1), 1, 0)) 
  ) %>%
  filter(keep==1) %>%
  select(!keep)

# also remove those from the captures
reduced.smammal.captures <- smammal.captures.clean %>%
  filter((paste(plotID, collectDate) %in% paste(reduced.trapping.session.info$plotID, 
                                                reduced.trapping.session.info$collectDate)))

reduced.plots <- sort(unique(reduced.smammal.captures$plotID))
summary(reduced.plots %in% sort(unique(reduced.trapping.session.info$plotID)))
#    Mode    TRUE 
# logical     231 
summary(sort(unique(reduced.trapping.session.info$plotID)) %in% reduced.plots)
#    Mode   FALSE    TRUE 
# logical       7     231  # this is ok if no animals were caught. seems to be the case.

# function below needs them to be the same length, so remove plots from 
# reduced.trapping.session.info not in reduced.plots
reduced.trapping.session.info <- reduced.trapping.session.info %>%
 filter(plotID %in% reduced.plots)


save(smammal.plot.data, smammal.captures.clean, trapping.session.info,
     reduced.smammal.captures, reduced.trapping.session.info, reduced.plots,
     file="NEONsmammalcapturesClean.RData")


############################################################################
## Format NEON data from all sites all species so I can run Bayesian
## Closed Robust design models to estimate a species-specific capture
## rate, p. I will use this to estimate abundance for each species
## for each month for each site.
##
## see "FunctionsForAllSpeciesEstimations.R" for data functions
## See "ModelRunAllSpeciesEstimation.R" for model and run code
## Actually using "ModelRunSpeciesTrapMMCaptureRateEstimations.R" which accounts 
## variation in traps set per night
############################################################################


taxon.cap.sum <- reduced.smammal.captures %>%
  group_by(taxonID) %>%
  summarise(genus = unique(genus), 
            species = unique(species),
            n=n()) 
taxon.cap.sum$species.num <- 1:dim(taxon.cap.sum)[1]
genera <- sort(unique(taxon.cap.sum$genus))
taxon.cap.sum$genus.num <- match(taxon.cap.sum$genus,genera)
print(taxon.cap.sum, n=200)
# probably should group all the Peromyscus speciesGroups together with Peromyscus sp.

write_excel_csv(taxon.cap.sum, file="SmammalTaxonCaps.csv")

# Use functions to format data for all sites all years -----------------------#

NEON.CHlong.all <- NEON.pRDcapture.history.long.fun(
    cleaned.data = reduced.smammal.captures, 
    ID.col = "plot_taxon_tag", # column to identify individuals
    session.info = reduced.trapping.session.info, 
    cols = c("genus", "species"), # which column names to keep as individual covariates
    taxon = "all")
# takes < a minute

save(NEON.CHlong.all, taxon.cap.sum, file="NEONchForPmodels.RData")

