###############################################################################
## Cleaning and Organizing Raw SMammal Data for abundance estimation
###############################################################################

library(tidyverse)
load("rawSmammalCaptureData.RData")
names(smammal.capturedata)

#### raw data with all captures from all sites and dates that weren't provisional as of 7/31/24
#### includes lines for empty traps 

NEONsmammalcaptures <- as_tibble(smammal.capturedata$mam_pertrapnight)
names(NEONsmammalcaptures)
head(NEONsmammalcaptures) 

## Dates & trapping schedule & trap-nights ----------------------------------##

site.names <- sort(unique(NEONsmammalcaptures$siteID))
smammal.sitedates <- list()
trapStatuslist <- list()
for(i in 1:length(site.names)){
  dat <- NEONsmammalcaptures[which(NEONsmammalcaptures$siteID==site.names[i]),]
  plots <- sort(unique(dat$plotID))
  plotdates <- list()
  trapstat <- list()
  for(p in 1:length(plots)){
    inds <- which(dat$plotID==plots[p])
    dates <- dat$collectDate[inds]
    plotdates[[p]] <- sort(unique(dates))
    trapstat[[p]] <- data.frame(date=plotdates[[p]],trap_not_set=rep(NA,length(plotdates[[p]])),
                                trap_disturbed=rep(NA,length(plotdates[[p]])), 
                                trap_open_but_sign=rep(NA,length(plotdates[[p]])),
                                multiple_caps=rep(NA,length(plotdates[[p]])),
                                capture=rep(NA,length(plotdates[[p]])),
                                set_and_empty=rep(NA,length(plotdates[[p]])),
                                traps_available=rep(NA,length(plotdates[[p]]))) 
    
    # I counted traps available as total number minus those that weren't set (still counted disturbed/closed as available) 
    # so trap nights per session will be the sum of this over the session
    
    datx <- dat[inds,]  
    for(d in 1:length(plotdates[[p]])){
      datxd <- datx[which(dates==plotdates[[p]][d]),]
      trapstat[[p]]$trap_not_set[d] <-  length(which(datxd$trapStatus=="1 - trap not set"))
      trapstat[[p]]$trap_disturbed[d] <-  length(which(datxd$trapStatus=="2 - trap disturbed/door closed but empty"))
      trapstat[[p]]$trap_open_but_sign[d] <-  length(which(datxd$trapStatus=="3 - trap door open or closed w/ spoor left"))
      trapstat[[p]]$multiple_caps[d] <-  length(which(datxd$trapStatus=="4 - more than 1 capture in one trap"))
      trapstat[[p]]$capture[d] <-  length(which(datxd$trapStatus=="5 - capture"))
      trapstat[[p]]$set_and_empty[d] <-  length(which(datxd$trapStatus=="6 - trap set and empty"))
      trapstat[[p]]$traps_available[d] <-  sum(trapstat[[p]][d, 3:7]) 
      
    }
    ## sometimes no traps were set the whole day or even session. I don't see the point of keeping that in the dataset.
    ## Remove those dates
    rmind <- which(trapstat[[p]]$traps_available==0)
    if(length(rmind)>0){
      trapstat[[p]] <- trapstat[[p]][-rmind,]
      plotdates[[p]] <- plotdates[[p]][-rmind]
    }
  }
  names(plotdates) <- plots
  names(trapstat) <- plots
  smammal.sitedates[[i]] <- plotdates
  trapStatuslist[[i]] <- trapstat
}
names(smammal.sitedates) <- site.names
names(trapStatuslist) <- site.names





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
# eventually put those data in here (SNV through 2019 and after that tick borne pathogens)

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




## Plot Summaries including trapnights ------------------------------------- ##

smammal.plot.data <- NEONsmammalcaptures %>%
  group_by(siteID, plotID, plotType, nlcdClass) %>% # there area few with slightly diff elevation, long, lat so average
  summarise(latitude = mean(decimalLatitude),
            longitude = mean(decimalLongitude),
            elevation = mean(elevation),
            n.dates = length(unique(floor_date(collectDate, unit="days"))),
            first.date = min(collectDate),
            last.date = max(collectDate),
            total.richness = length(na.omit(unique(taxonID))))

sites <- sort(unique(smammal.captures$siteID))
plots <- sort(unique(smammal.captures$plotID))

trapnights <- numeric()
for(i in 1:dim(smammal.plot.data)[1]){
  site <- smammal.plot.data$siteID[i]
  plot <- smammal.plot.data$plotID[i]
  statl <- trapStatuslist[[site]]
  trapnights[i] <- sum(statl[[plot]]$traps_available)
}
smammal.plot.data$trapnights <- trapnights

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
            mean.IDsper100trapnight = mean(IDsper100trapnights))


save(smammal.captures, smammal.plot.data, smammal.sitedates, trapStatuslist ,file="smammalcaptures.RData")


###############################################################################
# Clean data for abundance estimation
###############################################################################

# remove those without tag numbers
reduced.smammal.captures <- smammal.captures %>%
  filter(!is.na(tagID))


# remove some 'species' that aren't helpful
sp.rm <- c("DIVI", # opossum 6 captures
           "LEAM", # rabbit 5 captures
           "HESP", # Heteromyidae family
           "OTHE", # Mammalia
           "ROSP", # Rodentia
           "SRSP" # Soricidae
)
reduced.smammal.captures <- reduced.smammal.captures %>%
  filter( !taxonID %in% sp.rm ) 

# collapse a couple subspecies to species
reduced.smammal.captures$taxonID[which(reduced.smammal.captures$taxonID=="SIHE")] <-"SIHI"
reduced.smammal.captures$taxonID[which(reduced.smammal.captures$taxonID=="ICTM")] <-"ICTR"


# there are 1982 "PELAPEMA" P. leu and P. man
# a lot more Sorex were just ID to genus that to species, so only estimate p for genus


# add genus and species columns to help with grouping for estimation
reduced.smammal.captures$genus <- str_split_i(reduced.smammal.captures$scientificName," ",1)
reduced.smammal.captures$species <- paste(reduced.smammal.captures$genus , 
                                          str_split_i(reduced.smammal.captures$scientificName," ",2))



#### Quality Control ---------------------------------------------------------#

# Are there duplicates of same tag/plot on same date? ------------------- ##
# Yes, lots:
cols <- names(reduced.smammal.captures)
cols <- cols[-match(c("collectDate", "tagID", "plotID"),cols)]
dupes.date <- reduced.smammal.captures %>%
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
  reduced.smammal.captures %>%
  filter(collectDate==ymd("2014-07-25") & tagID == "NEON.MAM.D01.R0433" &  plotID=="HARV_001")
)
# same tag number: male PEMA and female PELE,  same plot and date. Ugh.



# Is the same tag IDed as different species or sexes? ------------------- ##
# Yes, lots:
dupes.taxsex <- reduced.smammal.captures %>%
  group_by(tagID, plotID) %>%
  summarise(across(c(taxonID, sex), 
                   ~ length(unique(.x)), .names = "n_{.col}")) %>%
  filter(n_taxonID>1 | n_sex>1)

dim(dupes.taxsex)
# dim(dupes.taxsex)
# [1] 5686    4
# Yikes. 5686 have multiple species or sexes
length(which(dupes.taxsex$n_taxonID>1))
# 3207 tagIDs have multiple species assigned (within the same plot)
# I guess I'll make unique ID be plot_taxon_tag 
# (it's possible some individuals move between plots within the same site,
# but here I will consider them separate.)
# right now I don't care about sex, so will leave those to deal with later.

reduced.smammal.captures <- reduced.smammal.captures %>%
  mutate(plot_taxon_tag = paste(plotID, taxonID, tagID, sep="_"))




source("FunctionsForNEONdata.R")


##create a list of dataframes of sessions for each plot
prim.session.list <- NEON.primary.session.function(smammal.sitedates, int.break=10)


# remove sessions that only have 1 night (won't help with estimation)
reduced.prim.session.list <- lapply(prim.session.list, function(x){
  dat <- x %>%
    group_by(prim.session) %>%
    summarise(n.days = length(unique(date)))
  if(length(which(dat$n.days==1))>0){
    rm.session <- dat$prim.session[which(dat$n.days==1)]
    x[-match(rm.session, x$prim.session),]
  } else {
    x
  }
})
which.empty <- which(unlist(lapply(reduced.prim.session.list, function(x){
  dim(x)[1]}))==0)
reduced.prim.session.list <- reduced.prim.session.list[-which.empty]

removed.sessions <- data.frame(plotID = character(), date = Date(), prim.session = numeric() )
for(i in 1:length(prim.session.list)){
  x <- prim.session.list[[i]]
  dat <- x %>%
    group_by(prim.session) %>%
    summarise(n.days = length(unique(date)))
  rm.session <- dat$prim.session[which(dat$n.days==1)]
  if(length(rm.session)>0){
    removed.sessions <- rbind(removed.sessions,
                              data.frame(plotID = names(prim.session.list)[i], x[match(rm.session, x$prim.session),1:2]))
  }
}  
# remove those from the captures
reduced.smammal.captures <- reduced.smammal.captures %>%
  filter(!(paste(plotID, collectDate) %in% paste(removed.sessions$plotID, removed.sessions$date)))

reduced.plots <- sort(unique(reduced.smammal.captures$plotID))
summary(reduced.plots %in% names(reduced.prim.session.list))
#    Mode    TRUE 
# logical     231 
summary(names(reduced.prim.session.list) %in% reduced.plots)
#    Mode   FALSE    TRUE 
# logical       7     231  #<---------- this is ok if no animals were caught. seems to be the case.

# function below needs them to be the same length, so remove plots from 
# reduced.prim.session.list not in reduced.plots
rm <- which(!names(reduced.prim.session.list) %in% reduced.plots)
reduced.prim.session.list <- reduced.prim.session.list[-rm]


save(smammal.plot.data, smammal.captures, smammal.sitedates,
     trapStatuslist , prim.session.list,
     reduced.smammal.captures, reduced.prim.session.list, reduced.plots,
     file="NEONsmammalcapturesClean.RData")


############################################################################
## Format NEON data from all sites all species so I can run Bayesian
## Closed Robust design models to estimate a species-specific capture
## rate, p. I will use this to estimate abundance for each species
## for each month for each site.
##
## see "FunctionsForAllSpeciesEstimations.R" for data functions
## See "ModelRunAllSpeciesEstimation.R" for model and run code
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
    prim.session.list = reduced.prim.session.list, # list with plots, dates, and prim.session
    cols = c("genus", "species"), # which column names to keep as individual covariates
    taxon = "all")
# takes < a minute

save(NEON.CHlong.all, taxon.cap.sum, file="NEONchForPmodels.RData")

