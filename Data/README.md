## Description of data files... 

"SmammalTraitData.csv" contains small mammal diet and activity data from the EltonTraits database (Wilman et al. 2014), and life history data from the Amniote database (Myhrvold et al. 2015). Scientific2 may have a synonym the data were listed under.

For small mammal abundances, I estimated abundance per plot per primary capture occasion (usually over 3 days), accounting for species specific capture rates, using closed Bayesian mark-recpature models. Since number of traps set varied, I assumed daily probabiltiy of detection increased with number of traps set according to a Michaelis-Menten equation with asymptote of 1 (if enough traps are set, detection probability will be 1). All species have an intercept of 0 (if 0 traps are set, probabiltiy of detection is 0). I estimated a half saturation constant for each species. I used these estimates along with number of traps set and days trapped to get a per-primary occasion capture probaility, and divided number of individuals detected per primary occasion by the probability of detection to get a abundance estimate per species per primary occasion. I have these data in longdata form, but did not upload them here because large.

"SmammalSiteHabitatYear.csv" has the average of those abundance estimates per site-habitat per year but only during the 4-month window of peak growing season. (Where each row is a site-habitat-year and species as columns). Since I didn't estimate productivity for ag and wetland sites, those are not here.