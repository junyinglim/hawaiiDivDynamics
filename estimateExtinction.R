## ESTIMATING TOTAL EXTINCTION USING THE ISLAND ONTOGENY MODEL
# Author: Charles Marshall and Jun Ying Lim
# Simulates species richness for time varying logistic growth with a carrying capacity (K_max) where both speciation and extinction are diversity dependent with equal but opposite strength (alpha), where the maximum equilibrium values of p (origination) and q (extinction) = 0.1*r_0_max, and the initial rate of extinction is 0.

## DIRECTORIES =============================================
stem.dir <- "~/Dropbox/Hawaii Diversity Dynamics/"
data.dir <- file.path(stem.dir, "raw_data")
main.dir <- file.path(stem.dir, "analyses")
output.dir <- file.path(main.dir, "results")
figure.dir <- file.path(main.dir, "figures")

# Import data
islClades <- read.csv(file.path(output.dir, "hawaii_taxon_data.csv"), stringsAsFactors = FALSE)
islTimes <- read.csv(file.path(output.dir, "islTimes.csv"), stringsAsFactors = FALSE)
islTimes_small <- read.csv(file.path(output.dir, "islTimes_smallIsl.csv"), stringsAsFactors = FALSE) # assuming kauai and oahu were smaller
islTimes_large <- read.csv(file.path(output.dir, "islTimes_largeIsl.csv"), stringsAsFactors = FALSE) # assuming kauai and oahu were larger


#  For Maui Nui for LOBELIADS

simulateDiversity <- function(nSim, t_g, t_d, p, K_max, r_0_max, dt =0.01){
  # Simulate diversity using the island ontogeny model for one island
  # Args:
  #     nSim = No of simulations
  #     t_g = Duraton of growth (Myr)
  #     t_d = Duration of decay (Myr)
  #     p = proportion decay, from eq 3. and 8
  #     r_0_max = intrinsic maximum species accumulation rate (immigration + speciation - extinction)
  #     K_max = clade level carrying capacity
  res <- list()
  
  for(i in 1:nSim){
    
    TaxonID <- vector(); TaxonID[1] <- 1
    Ancestor <- vector(); Ancestor[1] <- 0
    T_origin <- vector(); T_origin[1] <- 0
    T_extinct <- vector(); T_extinct[1] <- NA
    
    for (t1 in seq(0.01, t_g, by = dt)){	# This starts the simulation at time 0, not dt.
      S_alive <- sum(is.na(T_extinct))
      p = (r_0_max*(t1/t_g) - 0.9*(r_0_max/K_max)*S_alive)*dt # speciation rate
      q = 0.1*(r_0_max/K_max)*S_alive*dt # extinction rate
      for (j in 1:length(TaxonID)){			# For every taxon, extant or extinct.
        if(is.na(T_extinct[j])){  # if taxa is not extinct, continue..
          if (q > runif(1, 0, 1)){ 	 	# If true, taxon become extinct. 
            T_extinct[j] <- t1		# Assign its extinction time.
            #S_alive <- S_alive-1		# Decrease count of extant taxa.
          }		
          
          if (p > runif(1, 0, 1)){ 		# If true taxon speciates.
            S_total <- length(TaxonID) + 1
            TaxonID[S_total] <- S_total # New taxon given unique number as name.
            Ancestor[S_total] <- j	# The name of its ancestor.
            T_origin[S_total] <- t1	# Its time of origin
            T_extinct[S_total] <- NA
          }
        }
      }
    }
    
    # Speciation and extinction during the decay phase
    if(t_d > 0){
      for (t2 in seq(0.01, t_d, by = 0.01)){ # start at 0.01
        S_alive <- sum(is.na(T_extinct))
        p = (r_0_max*(1 - p * (t2/t_d)) - 0.9*(r_0_max/K_max)*S_alive)*dt 
        q = 0.1*(r_0_max/K_max)*S_alive*dt
        
        for (m in 1:length(TaxonID)){ # For every taxon, extant or extinct.
          if(is.na(T_extinct[m])){  # Thus, extinct taxa are ignored.
            if (q > runif(1, 0, 1)){    # If true, taxon become extinct. 
              T_extinct[m] <- t_g + t2    # Assign its extinction time.
            }       
            
            if (p > runif(1, 0, 1)){    # If true taxon speciates.
              S_total <- length(TaxonID)+1
              TaxonID[S_total] <- S_total      # New taxon given unique number as name.
              Ancestor[S_total] <- m         # The name of its ancestor.
              T_origin[S_total] <- t_g + t2    # Its time of origin
              T_extinct[S_total] <- NA                # Increase count of extant taxa.
            }
          }
        }
      }
    }
    res[[i]] <- data.frame(TaxonID, Ancestor, T_origin, T_extinct)
  }
  return(res)
}

simStats <- function(df){
  # Takes the species tables generated and gives basic summary statistics for each stochastic run
  S_total <- nrow(df)
  S_extant <- sum(is.na(df$T_extinct))
  S_extinct <- S_total - S_extant
  propExtinct <- S_extinct / S_total
  temp <- subset(df, !is.na(T_extinct)) # subset extinct species
  avgLifeSpan <- mean(temp$T_extinct - temp$T_origin)
  temp2 <- subset(df, is.na(T_extinct)) # subset extant species
  crownAge <- min(temp2$T_origin, na.rm = TRUE) # WRONG need to fix
  data.frame(S_total, S_extant, S_extinct, propExtinct, avgLifeSpan, crownAge)
}

summarizeSummary <- function(df){
  # Takes the summary statistics and generates the mean species diversity etc
  meanS_total <- mean(df$S_total)
  meanS_extinct <- mean(df$S_extinct)
  meanProp_extinct <- mean(df$propExtinct)
}

simulateDiversityAllIslands <- function(islTimes, r_0_max, K_max, p, nSim = 100){
  res <- list()
  for(i in 1:4){
    temp <- simulateDiversity(t_g = islTimes$t1_mean[i],
                              t_d = islTimes$t2_mean[i],
                              r_0_max = r_0_max,
                              K_max = K_max,
                              p = p,
                              nSim)
    res[i] <- do.call("rbind", lapply(temp, simStats))
  }
  resTotal <- do.call("rbind", lapply(res, summarizeSummary))
  return(resTotal)
}




test <- simulateDiversity(nSim = 10, t_g = 1.52, t_d = 1.03, p = 0.424, K_max= 81, r_0_max = 7.36, dt = 0.0001)

summarize <- 
res <- do.call("rbind", lapply(test, simStats))
res

test2 <- simulateDiversityAllIslands(islTimes = islTimes, nSim = 1, p = 0.424, K_max = 81, r_0_max = 7.36)
test2
