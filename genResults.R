## DIVERSITY DYNAMICS OF THE HAWAIIAN ISLANDS
# Authors: Jun Ying Lim & Charles Marahall
# Last modified: 3 Oct 2016 (Version 2)
# Version update notes:
# Removed all the SAR analyses; randomization tests of island times
# Added in LGM area analyses

## DIRECTORIES =============================================
stem.dir <- "~/Dropbox/Hawaii Diversity Dynamics/"
data.dir <- file.path(stem.dir, "raw_data")
main.dir <- file.path(stem.dir, "analyses")
output.dir <- file.path(main.dir, "results")
figure.dir <- file.path(main.dir, "figures")

## PACKAGES =============================================
library(reshape2) # Utility functions for preparing data for ggplot2
library(plyr) # Rescaling of valbues
library(minpack.lm) # uses a different algorithm to get convergence
library(MuMIn) # AICc function
library(pracma) #error functions
library(Rmpfr) #error functions
library(ggplot2); library(RColorBrewer)
source(file.path(main.dir, "divDynamics.R")) # functions to calculate diversity and custom graphical functions

## INPUT SPECIES DATA =============================================
# FOr different island area datasets, comment out the right lines
islClades <- read.csv(file.path(output.dir, "hawaii_taxon_data.csv"), stringsAsFactors = FALSE)
islTimes <- read.csv(file.path(output.dir, "islTimes.csv"), stringsAsFactors = FALSE)
islTimes_small <- read.csv(file.path(output.dir, "islTimes_smallIsl.csv"), stringsAsFactors = FALSE) # assuming kauai and oahu were smaller
islTimes_large <- read.csv(file.path(output.dir, "islTimes_largeIsl.csv"), stringsAsFactors = FALSE) # assuming kauai and oahu were larger

#islTimes <- read.csv(file.path(output.dir, "islTimes_areaSensitivity.csv"), stringsAsFactors = FALSE)
#islTimes <- read.csv(file.path(output.dir, "islTimes_smallIsl.csv"), stringsAsFactors = FALSE)
#islTimes <- read.csv(file.path(output.dir, "islTimes_smallKa.csv"), stringsAsFactors = FALSE)

# Comment out for default results
#islTimes <- islTimes_small
#islTimes <- islTimes_large


## GENERATE ONTOGENY PARAMETERS =============================================
# Exclude equivocal groups
excludeClades <- c("Nesosydne", "Viola")
islClades <- subset(islClades, ! Taxon %in% excludeClades)

# Split into taxonomic groups
arthroClades <- subset(islClades, Taxon %in% c("PW Drosophila", "Tetragnatha", "Orsonwelles", "Ariamnes", "Laupala", "Platynini", "Megalagrion"))
plantClades <- subset(islClades, Group == "Plant")
vertClades <- subset(islClades, Taxon == "Honeycreepers")
cladesList <- list(vertClades, arthroClades, plantClades)

# ANALYSES (CURRENT AREAS) =============================================
# For arthropods, z = 0.471 (incl. Ha), z = 0.491 (excl. Ha)
# For birds, z = 0.271 (incl Ha), z = 0.271 (excl. Ha)

z_list <- list(0.271, 0.491, 0.38)
tol <- 0.05; currA_name <- "currA"

# Young island, short growth / long decline
t1_name <- "t1_sG_yngIsl"; t2_name <- "t2_sG_yngIsl"; maxA_name <- "maxA_oldOa"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_sG_yngIsl <- do.call("rbind", resultsList)
#write.csv(resultsList_sG_yngIsl, file.path(output.dir, "resultsList_sG_yngIsl_smallIsl.csv"), row.names = FALSE)
# write.csv(resultsList_sG_yngIsl, file.path(output.dir, "resultsList_sG_yngIsl_largeIsl.csv"), row.names = FALSE)
#write.csv(resultsList_sG_yngIsl, file.path(output.dir, "resultsList_sG_yngIsl_smallKa.csv"), row.names = FALSE)
write.csv(resultsList_sG_yngIsl, file.path(output.dir, "resultsList_sG_yngIsl.csv"), row.names = FALSE)


# Old island, short growth / long decline
t1_name <- "t1_sG_oldIsl"; t2_name <- "t2_sG_oldIsl"; maxA_name <- "maxA_oldOa"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_sG_oldIsl <- do.call("rbind", resultsList)
# write.csv(resultsList_sG_oldIsl, file.path(output.dir, "resultsList_sG_oldIsl_smallIsl.csv"), row.names = FALSE)
# write.csv(resultsList_sG_oldIsl, file.path(output.dir, "resultsList_sG_oldIsl_largeIsl.csv"), row.names = FALSE)
#write.csv(resultsList_sG_oldIsl, file.path(output.dir, "resultsList_sG_oldIsl_smallKa.csv"), row.names = FALSE)
write.csv(resultsList_sG_oldIsl, file.path(output.dir, "resultsList_sG_oldIsl.csv"), row.names = FALSE)

# Young island, long growth / short decline
t1_name <- "t1_lG_yngIsl"; t2_name <- "t2_lG_yngIsl"; maxA_name <- "maxA_yngOa"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_lG_yngIsl <- do.call("rbind", resultsList)
#write.csv(resultsList_lG_yngIsl, file.path(output.dir, "resultsList_lG_yngIsl_smallIsl.csv"), row.names = FALSE)
# write.csv(resultsList_lG_yngIsl, file.path(output.dir, "resultsList_lG_yngIsl_largeIsl.csv"), row.names = FALSE)
#write.csv(resultsList_lG_yngIsl, file.path(output.dir, "resultsList_lG_yngIsl_smallKa.csv"), row.names = FALSE)
write.csv(resultsList_lG_yngIsl, file.path(output.dir, "resultsList_lG_yngIsl.csv"), row.names = FALSE)

# Old island, long growth / short decline
t1_name <- "t1_lG_oldIsl"; t2_name <- "t2_lG_oldIsl"; maxA_name <- "maxA_yngOa"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_lG_oldIsl <- do.call("rbind", resultsList)
#write.csv(resultsList_lG_oldIsl, file.path(output.dir, "resultsList_lG_oldIsl_smallIsl.csv"), row.names = FALSE)
# write.csv(resultsList_lG_oldIsl, file.path(output.dir, "resultsList_lG_oldIsl_largeIsl.csv"), row.names = FALSE)
#write.csv(resultsList_lG_oldIsl, file.path(output.dir, "resultsList_lG_oldIsl_smallKa.csv"), row.names = FALSE)
write.csv(resultsList_lG_oldIsl, file.path(output.dir, "resultsList_lG_oldIsl.csv"), row.names = FALSE)

# Based on mean t_maxA and t_habitable
t1_name <- "t1_mean"; t2_name <- "t2_mean"; maxA_name <- "maxA"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_mean <- do.call("rbind", resultsList)
#write.csv(resultsList_mean, file.path(output.dir, "resultsList_mean_smallIsl.csv"), row.names = FALSE)
# write.csv(resultsList_mean, file.path(output.dir, "resultsList_mean_largeIsl.csv"), row.names = FALSE)
#write.csv(resultsList_mean, file.path(output.dir, "resultsList_mean_smallKa.csv"), row.names = FALSE)
write.csv(resultsList_mean, file.path(output.dir, "resultsList_mean.csv"), row.names = FALSE)

# Just for comparison with previous results
# t1_name <- "t1_mean"; t2_name <- "t2_mean"; maxA_name <- "maxA_oldOa"
# ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
#                         maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
# resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
# resultsList_mean_2 <- do.call("rbind", resultsList)


## ANALYSES (LGM AREAS) =============================================
z_list <- list(0.271, 0.491, 0.38)
islTimes_lgm <- islTimes
islTimes_lgm[islTimes$complex_id == "Hawaii", names(islTimes_lgm) %in% c("maxA", "maxA_oldOa", "maxA_yngOa")] <- islTimes_lgm$lgmA[islTimes$complex_id == "Hawaii"] # assigning maximum area as lgmA
tol <- 0.05; currA_name <- "lgmA"
# Before, the code only fixed the maxA column but not the others...

# Young island, short growth / long decline
t1_name <- "t1_sG_yngIsl"; t2_name <- "t2_sG_yngIsl"; maxA_name <- "maxA_oldOa"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_sG_yngIsl_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_sG_yngIsl_lgm, file.path(output.dir, "resultsList_sG_yngIsl_lgm.csv"), row.names = FALSE)

# Old island, short growth / long decline
t1_name <- "t1_sG_oldIsl"; t2_name <- "t2_sG_oldIsl"; maxA_name <- "maxA_oldOa"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_sG_oldIsl_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_sG_oldIsl_lgm, file.path(output.dir, "resultsList_sG_oldIsl_lgm.csv"), row.names = FALSE)

# Young island, long growth / short decline
t1_name <- "t1_lG_yngIsl"; t2_name <- "t2_lG_yngIsl"; maxA_name <- "maxA_yngOa"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_lG_yngIsl_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_lG_yngIsl_lgm, file.path(output.dir, "resultsList_lG_yngIsl_lgm.csv"), row.names = FALSE)

# Old island, long growth / short decline
t1_name <- "t1_lG_oldIsl"; t2_name <- "t2_lG_oldIsl"; maxA_name <- "maxA_yngOa"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_lG_oldIsl_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_lG_oldIsl_lgm, file.path(output.dir, "resultsList_lG_oldIsl_lgm.csv"), row.names = FALSE)

# Based on mean t_maxA and t_habitable
t1_name <- "t1_mean"; t2_name <- "t2_mean"; maxA_name <- "maxA"
ontParameters <- lapply(X = z_list, FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_mean_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_mean_lgm, file.path(output.dir, "resultsList_mean_lgm.csv"), row.names = FALSE)
