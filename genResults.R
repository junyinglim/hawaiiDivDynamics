## DIVERSITY DYNAMICS OF THE HAWAIIAN ISLANDS
# Authors: Jun Ying Lim & Charles Marahall
# Last modified: 22 April 2016

## DIRECTORIES =============================================
main.dir <- "~/Dropbox/Hawaii Diversity Dynamics/analyses"
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
library(parallel) # to parallelize the lapply
source(file.path(main.dir, "divDynamics.R")) # functions to calculate diversity and custom graphical functions

## INPUT SPECIES DATA =============================================
islClades <- read.csv(file.path(output.dir, "hawaii_taxon_data.csv"), stringsAsFactors = FALSE)
islTimes <- read.csv(file.path(output.dir, "islTimes.csv"), stringsAsFactors = FALSE)

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
tol <- 0.05; currA_name <- "currA"
# Young island, short growth / long decline
t1_name <- "t1_sG_yngIsl"; t2_name <- "t2_sG_yngIsl"; maxA_name <- "maxA_oldOa"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_sG_yngIsl <- do.call("rbind", resultsList)
write.csv(resultsList_sG_yngIsl, file.path(output.dir, "resultsList_sG_yngIsl.csv"), row.names = FALSE)

# Old island, short growth / long decline
t1_name <- "t1_sG_oldIsl"; t2_name <- "t2_sG_oldIsl"; maxA_name <- "maxA_oldOa"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_sG_oldIsl <- do.call("rbind", resultsList)
write.csv(resultsList_sG_oldIsl, file.path(output.dir, "resultsList_sG_oldIsl.csv"), row.names = FALSE)

# Young island, long growth / short decline
t1_name <- "t1_lG_yngIsl"; t2_name <- "t2_lG_yngIsl"; maxA_name <- "maxA_yngOa"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_lG_yngIsl <- do.call("rbind", resultsList)
write.csv(resultsList_lG_yngIsl, file.path(output.dir, "resultsList_lG_yngIsl.csv"), row.names = FALSE)

# Old island, long growth / short decline
t1_name <- "t1_lG_oldIsl"; t2_name <- "t2_lG_oldIsl"; maxA_name <- "maxA_yngOa"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_lG_oldIsl <- do.call("rbind", resultsList)
write.csv(resultsList_lG_oldIsl, file.path(output.dir, "resultsList_lG_oldIsl.csv"), row.names = FALSE)

# Based on mean t_maxA and t_habitable
t1_name <- "t1_mean"; t2_name <- "t2_mean"; maxA_name <- "maxA"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        maxA_name = maxA_name, t1_name = t1_name, t2_name = t2_name,currA_name = currA_name, df = islTimes)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_mean <- do.call("rbind", resultsList)
write.csv(resultsList_mean, file.path(output.dir, "resultsList_mean.csv"), row.names = FALSE)


## ANALYSES (LGM AREAS) =============================================
islTimes_lgm <- islTimes
islTimes_lgm$maxA[islTimes$complex_id == "Hawaii"] <- islTimes_lgm$lgmA[islTimes$complex_id == "Hawaii"]
tol <- 0.05; currA_name <- "lgmA"

# Young island, short growth / long decline
t1_name <- "t1_sG_yngIsl"; t2_name <- "t2_sG_yngIsl"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        t1_name = t1_name, t2_name = t2_name, currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_sG_yngIsl_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_sG_yngIsl_lgm, file.path(output.dir, "resultsList_sG_yngIsl_lgm.csv"), row.names = FALSE)

# Old island, short growth / long decline
t1_name <- "t1_sG_oldIsl"; t2_name <- "t2_sG_oldIsl"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        t1_name = t1_name, t2_name = t2_name, currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_sG_oldIsl_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_sG_oldIsl_lgm, file.path(output.dir, "resultsList_sG_oldIsl_lgm.csv"), row.names = FALSE)

# Young island, long growth / short decline
t1_name <- "t1_lG_yngIsl"; t2_name <- "t2_lG_yngIsl"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        t1_name = t1_name, t2_name = t2_name, currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_lG_yngIsl_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_lG_yngIsl_lgm, file.path(output.dir, "resultsList_lG_yngIsl_lgm.csv"), row.names = FALSE)

# Old island, long growth / short decline
t1_name <- "t1_lG_oldIsl"; t2_name <- "t2_lG_oldIsl"
ontParameters <- lapply(X = list(0.271, 0.47, 0.38), FUN = calcOntogeny,
                        t1_name = t1_name, t2_name = t2_name, currA_name = currA_name, df = islTimes_lgm)
resultsList <- mapply(testModelList, .data = cladesList, df = ontParameters, MoreArgs = list(ftol = tol), SIMPLIFY = FALSE)
resultsList_lG_oldIsl_lgm <- do.call("rbind", resultsList)
write.csv(resultsList_lG_oldIsl_lgm, file.path(output.dir, "resultsList_lG_oldIsl_lgm.csv"), row.names = FALSE)