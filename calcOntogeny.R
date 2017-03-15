## DIVERSITY DYNAMICS OF THE HAWAIIAN ISLANDS
# Authors: Jun Ying Lim & Charles Marshall
# Last modified: 27 May 2016

## DIRECTORIES =============================================
stem.dir <- "~/Dropbox/Hawaii Diversity Dynamics/"
data.dir <- file.path(stem.dir, "raw_data")
main.dir <- file.path(stem.dir, "analyses")
output.dir <- file.path(main.dir, "results")
figure.dir <- file.path(main.dir, "figures")

## PACKAGES =============================================
library(reshape2) # Utility functions for preparing data for ggplot2
library(plyr) # Rescaling of valbues
source(file.path(main.dir, "divDynamics.R")) # functions to calculate diversity and custom graphical functions

## CURRENT ISLAND AREAS =============================================
# (All numbers from data book; no longer from maps)
hawaii_currA <- 10433.55 # hawaii
maui_currA <- 1999.45 # maui
kaho_currA <- 115.5 #kahoolawe
lanai_currA <- 365.36 # lanai
molo_currA <- 674.58 # molokai
oahu_currA <- 1547.88 # oahu
kauai_currA <- 1430.59 # kauai
niihau_currA <- 175.09 # ni'ihau

## LAST GLACIAL MAXIMUM AREAS =============================================
# (All number from Weigelt et al 2016)
hawaii_lgmA <- 10^2.848 - 1 + 10^4.019
maui_lgmA <- 10^3.578 - 1 + 10^3.277 +
             10^2.085 - 1 + 10^2.067 # kahoolawe was still disconnected
#10^3.725 - 1 + 10^2.564 #lanai
#10^3.699 - 1 + 10^2.833 # molokai
oahu_lgmA <- 10^2.925 - 1 + 10^3.193 # includes mokoli'i islet
kauai_lgmA <- 10^2.711 - 1 + 10^3.159
complex_lgmA <- c(hawaii_lgmA, maui_lgmA, oahu_lgmA, kauai_lgmA)

## CURRENT ISLAND HEIGHTS  =============================================
## (All numbers from data book, except for west molokai)
hawaii_currH <- 4205 # hawaii (mauna kea)
emau_currH <- 3055 # east maui (haleakala)
wmau_currH <- 1764 # west maui (puukukui)
kaho_currH <- 452 # kahoolawe (puu moaulanui)
lanai_currH <- 1026 # lanai (lanaihale)
emolo_currH <- 1512 # east molokai (kamakou)
wmolo_currH <- 421 # west molokai (puu nana); taken from Clague & Sherrod (2014)
kool_currH <-  960 # east oahu (ko'olau range)
waia_currH <- 1220 # west oahu (ka'ala, wai'anae range)
kauai_currH <- 1598 # kauai (kawaikini)
niihau_currH <- coneheight(theta = 5, A = niihau_currA * 10^6) # niihau is a collapsed cone, so instead of its current height (i.e., 381 m), we estimate its "true" height from its area

## MAXIMUM ISLAND HEIGHTS  =============================================
# All numbers from Carson & Clague (1995) In: Wagner & Funk (1995))
# Estimated by adding deepest break-in slope with current height
hawaii_maxH <- 4205 # hawaii (mauna kea)  ## isn't it 4600?
emau_maxH <- 5000 # east maui (haleakala)  
wmau_maxH <- 3400 # west maui (puukukui)
kaho_maxH <- 2100 # kahoolawe (puu moaulanui)
lanai_maxH <- 2200 # lanai (lanaihale)
emolo_maxH <- 3300 # east molokai (kamakou)
wmolo_maxH <- 1600 # west molokai (puu nana); taken from Clague & Sherrod (2014)
kool_maxH <-  1900 # east oahu (ko'olau range)
waia_maxH <- 2200 # west oahu (ka'ala, wai'anae range)
kauai_maxH <- 2600 # kauai (kawaikini)
niihau_maxH <- 1400 # niihau (paniau)

# Calculate maximum areas of island complexes
hawaii_maxA <- onecone(maxH = hawaii_maxH, currH = hawaii_currH, currA = hawaii_currA)
maui_maxA <- 13500 # Price & Elliot-Fisk (2004) gives a high stand value of 14000 km2. We deducted 500 for Penguin Bank as we don't believe it was ever high enough
# NOTE: we also calculated them as individual cones and summed them up - you get virtually identical results

oahu_maxA <- twocone(maxHa = kool_maxH, maxHb = waia_maxH,
                     currHa = kool_currH, currHb = waia_currH,
                     currA = oahu_currA) #+ 
              #conearea(theta = 5, H = 1000) # Ka'ena (excl. because was never high enough to catch rainfall)

# Alternate method that calculates area from maximum height alone (and assuming slope of 5.5 deg)
# Assumes all volcanic islands are essentailly cones and does not account for shield overlap
oahu_maxA_small <- pi * (kool_maxH / tan(5.5*2*pi/360))^2 / 1E6 + 
                   pi * (waia_maxH / tan(5.5*2*pi/360))^2 / 1E6

# Alternate method that assumes a steeper height of repose for oahu and kauai
repose_ratio <- tan(7 * 2 * pi / 360) / tan(5.5*2*pi/360)
oahu_maxA_large <- twocone(maxHa = kool_maxH * repose_ratio,
                           maxHb = waia_maxH * repose_ratio,
                           currHa = kool_currH,
                           currHb = waia_currH,
                           currA = oahu_currA)

kauai_maxA <- onecone(maxH = kauai_maxH,
                      currH = kauai_currH,
                      currA = kauai_currA)
              #onecone(maxH = niihau_maxH,
              #        currH = niihau_currH,
              #        currA = niihau_currA) # Ni'ihau (excl. because Niihau was never connected to Kauai subaerially)

# Alternate method that calculates area from maximum height alone (and assuming slope of 5.5 deg)
kauai_maxA_small <- pi * (kauai_maxH / tan(5.5*2*pi/360))^2 / 1E6 

# Alternate method that assumes different angle of repose
repose_ratio <- tan(7 * 2 * pi / 360) / tan(5.5*2*pi/360)
kauai_maxA_large <- onecone(maxH = kauai_maxH * repose_ratio,
                            currH = kauai_currH,
                            currA = kauai_currA)

# Adjust maximum area with decay rate 
# Wai'anae will have lost area when Ko'olau is at max area
waia_currA <- oahu_currA * waia_currH^2 / (waia_currH^2 + kool_currH^2) # calculate current area of wai'anae range

waia_maxA <- onecone(currA = waia_currA, maxH = waia_maxH, currH = waia_currH) # calculate maximum area of wai'ane using height ratio method (original method)
waia_maxA_large <- onecone(currA = waia_currA, maxH = waia_maxH * repose_ratio, currH = waia_currH) # calculate maximum area of wai'ane using height ratio method (original method)
waia_maxA_small <- pi * (waia_maxH / tan(5.5*2*pi/360))^2 / 1E6

# Uncertainty in times
# Assuming maximum dates
adjOahuArea <- function(waia_maxA, oahu_maxA){
  waia_ArealDecayRate <- (waia_maxA - waia_currA) / 3.55 
  waia_lag <- 3.55 - 1.8 # time between max area achieved by wai'anae vs. ko'olau
  oahu_adjA_old <- waia_ArealDecayRate * waia_lag
  oahu_maxA_old <- oahu_maxA - oahu_adjA_old
  
  waia_ArealDecayRate <- (waia_maxA - waia_currA) / 3.06
  waia_lag <- 3.06 - 0.8 # time between max area achieved by wai'anae vs. ko'olau
  oahu_adjA_young <- waia_ArealDecayRate * waia_lag
  oahu_maxA_young <- oahu_maxA - oahu_adjA_young
  
  # Assuming mean dates
  waia_ArealDecayRate <- (waia_maxA - waia_currA) / 3.305
  waia_lag <- 3.305 - 1.3 # time between max area achieved by wai'anae vs. ko'olau
  oahu_adjA_mean <- waia_ArealDecayRate * waia_lag
  oahu_maxA_mean <- oahu_maxA - oahu_adjA_mean
  
  return(c(oahu_maxA_mean, oahu_maxA_young, oahu_maxA_old))
}

oahu_maxA_adj <- adjOahuArea(waia_maxA = waia_maxA, oahu_maxA = oahu_maxA)
oahu_maxA_large_adj <- adjOahuArea(waia_maxA = waia_maxA_large, oahu_maxA = oahu_maxA_large)
oahu_maxA_small_adj <- adjOahuArea(waia_maxA = waia_maxA_small, oahu_maxA = oahu_maxA_small)

# Combine data
complex_maxA <- c(hawaii_maxA, maui_maxA, oahu_maxA_adj[1], kauai_maxA)
complex_maxA_yngOa <- c(hawaii_maxA, maui_maxA, oahu_maxA_adj[2], kauai_maxA)
complex_maxA_oldOa <- c(hawaii_maxA, maui_maxA, oahu_maxA_adj[3], kauai_maxA)

complex_maxA_large <- c(hawaii_maxA, maui_maxA, oahu_maxA_large_adj[1], kauai_maxA_large)
complex_maxA_large_yngOa <- c(hawaii_maxA, maui_maxA, oahu_maxA_large_adj[2], kauai_maxA_large)
complex_maxA_large_oldOa <- c(hawaii_maxA, maui_maxA, oahu_maxA_large_adj[3], kauai_maxA_large)

complex_maxA_small <- c(hawaii_maxA, maui_maxA, oahu_maxA_small_adj[1], kauai_maxA_small)
complex_maxA_small_yngOa <- c(hawaii_maxA, maui_maxA, oahu_maxA_small_adj[2], kauai_maxA_small)
complex_maxA_small_oldOa <- c(hawaii_maxA, maui_maxA, oahu_maxA_small_adj[3], kauai_maxA_small)


complex_currA <- c(hawaii_currA,
                   sum(c(maui_currA, kaho_currA, lanai_currA, molo_currA)),
                   oahu_currA,
                   kauai_currA)

islTimes <- data.frame(complex_id = c("Hawaii", "Maui Nui", "Oahu", "Kauai"),
                       currA = complex_currA,
                       maxA = complex_maxA,
                       lgmA = complex_lgmA,
                       maxA_yngOa = complex_maxA_yngOa,
                       maxA_oldOa = complex_maxA_oldOa)


## ISLAND GROWTH AND DECAY TIMES =============================================
# numbers from Clague & Sherrod (2014)
#islTimes$t1_min = c(1.32, 1.98, 2.92, 2.02) # end of shield building
#islTimes$t1_max = c(1.32, 2.12, 2.92, 2.37) # until end of late stage of shield building (which is mostly alkalic volcanics)
#islTimes$t2 = c(1.32, 3.08, 4.72, 6.02) # not quite age of island, time for which it has been habitable

## New times
islTimes$t_habitable_min <- c(1.1, 2.1, 3.93, 6)
islTimes$t_habitable_max <- c(1.3, 3, 4.34, 6.3)
islTimes$t_maxA_min <- c(0, 0.96, 0.8, 3.65)
islTimes$t_maxA_max <- c(0, 1.1, 1.8, 4.0)

islTimes$t_habitable_mean <- c(1.2, 2.55, 4.135, 6.15)
islTimes$t_maxA_mean <- c(0, 1.03, 1.3, 3.825)

# The longest t1s are bracketed by the youngest t_max and the oldest t_habitable
# The shortest t1s by the oldest t_max and the youngest t_habitable
islTimes$t1_sG_oldIsl <- islTimes$t_habitable_max - islTimes$t_maxA_max
islTimes$t2_sG_oldIsl <- islTimes$t_maxA_max

islTimes$t1_lG_oldIsl <- islTimes$t_habitable_max - islTimes$t_maxA_min
islTimes$t2_lG_oldIsl <- islTimes$t_maxA_min

islTimes$t1_sG_yngIsl <- islTimes$t_habitable_min - islTimes$t_maxA_max
islTimes$t2_sG_yngIsl <- islTimes$t_maxA_max

islTimes$t1_lG_yngIsl <- islTimes$t_habitable_min - islTimes$t_maxA_min
islTimes$t2_lG_yngIsl <- islTimes$t_maxA_min

islTimes$t1_mean <- islTimes$t_habitable_mean - islTimes$t_maxA_mean
islTimes$t2_mean <- islTimes$t_maxA_mean

# Add island areas

## ALTERNATIVE AREAS
islTimes_large <- islTimes 
islTimes_large$maxA <- complex_maxA_large
islTimes_large$maxA_yngOa <- complex_maxA_large_yngOa
islTimes_large$maxA_oldOa <- complex_maxA_large_oldOa

islTimes_small <- islTimes 
islTimes_small$maxA <- complex_maxA_small
islTimes_small$maxA_yngOa <- complex_maxA_small_yngOa
islTimes_small$maxA_oldOa <- complex_maxA_small_oldOa

# Export data
write.csv(islTimes, file.path(output.dir, "islTimes.csv"), row.names = FALSE)
write.csv(islTimes_large, file.path(output.dir, "islTimes_largeIsl.csv"), row.names = FALSE)
write.csv(islTimes_small, file.path(output.dir, "islTimes_smallIsl.csv"), row.names = FALSE)