require("reshape2") # Utility functions for preparing data for ggplot2
require("plyr") # Rescaling of values
require("pracma") # error functions
require("Rmpfr") # for more precise floating points
require("ggplot2")
require("gridExtra")

## FUNCTIONS CALCULATING THE ONTOGENY PARAMETERS FROM GEOLOGIC INFORMATION ========================================================
twocone <- function(currA, maxHa, maxHb, currHa, currHb){
	# Calculate the maximum area of a single cone given
	# the ratio of the maximum and current height squared, and its current area
	#
	# Args:
	# 	currA = current area
	# 	maxHa = maximum height of cone A
	# 	currHa = current height of cone A
	# 	maxHb = maximum height of cone B
	# 	currHb = current height of cone B
	#
	# Returns:
	# 	The maximum area of a two-cone island complex given the proportional ratio of their maximum and current heights

	maxA <- currA * (maxHa^2 + maxHb^2) / (currHa ^2 + currHb^2)
	return(maxA)
}

onecone <- function(currA, maxH, currH){
	# Calculate the maximum area of a single cone given
	# the ratio of the maximum and current height squared, and its current area
	#
	# Args:
	# 	currA = current area
	# 	maxH = maximum height
	# 	currH = current height
	#
	# Returns:
	# 	The maximum area of a cone given the ratio of its maximum and current heights

	maxA <- currA * maxH ^2 / currH ^2
	return(maxA)
}

conearea <- function(theta, H){
	# Calculate surface area of a cone (excluding base) given angle of cone and height of cone
	#
	# Args:
	# 	theta: angle of cone
	# 	H: height of cone
	#
	# Returns:
	# 	The surface area of a cone (excluding base) of angle theta and height H.

	pi * H / tan(theta*pi/180) * sqrt(H^2 + (H / tan(theta*pi/180))^2) /10^6
}

coneheight <- function(theta, A){
	# Calculate height of cone given surface area (excluding base) and angle of cone
	#
	# Args:
	# 	theta: angle of cone
	#	  A: surface area of cone (excluding base)
	#
	# Returns:
	# 	The height of a cone of angle theta and surface area A.
	sqrt((A * tan(theta*pi/180)) / (pi * sqrt(1 + 1/(tan(theta*pi/180)^2))))

}
	
randIslTimes <- function(n, df, sd = 0.1, method){
  ## Generates randomly drawn island times from normal and exponential probability ditributions
  #
  # Args:
  #   n: number of randomizations
  #   df: original island times
  #   sd: standard deviation of normal distribution / multiplier of exponential distribution
  #   method: "exp" for exponential distribution, "normal" for normal distribution
  #
  # Returns:
  #   List of data frames containing island times (t1_rand and t2_rand)
  complex_id = c("Hawaii", "Maui Nui", "Oahu", "Kauai")
  res <- lapply(1:n, FUN = function(x){
    z <- lapply(1:4, FUN = function(x){
      if(method == "normal"){
        t1_rand <- rnorm(n = 1, sd = sd, mean = df[x,]$t1_min)
        t2_rand <- rnorm(n = 1, sd = sd, mean = df[x,]$t2)  
      } else if(method == "exp"){
        t1_rand <- df[x,]$t1_min + (rexp(n = 1) * sd)
        t2_rand <- df[x,]$t2 + (rexp(n = 1) * sd)
      }
      complex_id <- complex_id[x]
      data.frame(t1_rand, t2_rand, complex_id, stringsAsFactors = FALSE)
    })
    islTimes <- do.call("rbind", z)
    islTimes$t2_rand[islTimes$complex_id == "Hawaii"] <- islTimes$t1_rand[islTimes$complex_id == "Hawaii"] # Hawaii t1_rand = t2_rand
    islTimes <- merge(islTimes, df, by = "complex_id", sort = FALSE) # prevents sorting of islands; correct sequence of islands CRUCIAL for model fitting
    return(islTimes)
  })
  return(res)
}


calcOntogeny <- function(df, z, t1_name, t2_name, currA_name, maxA_name){
	# Calculates relative maximum K, relative current K and proportion of K lost, from a dataframe where each row represents a different island in the archipelago
	#
	# Args:
	#	df: data.frame containing current area, maximum area, 
	#	z: z-score for which areal parameters should be scaled
	# 	t1_name: name of column for t1
	# 	t2_name: name of column for t2
  #   currA_name: name of column for current area
  #   maxA_name: name of column for maximum area
	#
	# Returns:
	# 	data.frame containing all the ontogeny parameters, including input variables

	rel_currA <- df[[currA_name]] / max(df[[currA_name]])
	#p1 <- (1 - (df$currA / df$maxA))^z 		# proportion of carrying capacity lost
	p1 <- 1 - (df[[currA_name]] / df[[maxA_name]])^z 		# proportion of carrying capacity lost
	maxK <- df[[maxA_name]]^z 						# maximum carrying capacity
	currK <- df[[currA_name]]^z
	rel_currK <- currK / max(currK)		# relative current carrying capacity
	rel_maxK <- maxK / max(maxK)		# relative maximum carrying capacity
  t_total <- df[[t1_name]] + df[[t2_name]]
	 
	return(data.frame(currA = df[[currA_name]],
					  maxA = df[[maxA_name]],
					  p1 = p1,
					  rel_maxK = rel_maxK,
					  rel_currK = rel_currK,
					  rel_currA = rel_currA,
					  z = z,
					  t1 = df[[t1_name]],
					  t2 = df[[t2_name]],
				  	t_total = t_total))
}	


## DEFINE DIVERSIFICATION MODELS ========================================================
CR <- function(d, t){
	# Computes the expected diversity given diversification rate and time
	#
	# Args:
	# 	t: vector of time
	# 	d: diverisifcation rate
	# 
	# Returns: 
	#	The expected diversity for each element of t.
	exp(d * t)
}


DD <- function(d, K, t){
	# Computes the expected diversity given density-dependent diversification rate and time
	#
	# Args:
	# 	t: vector of time
	# 	d: diverisifcation rate
	# 	K: carrying capacity (i.e., strength of density dependence)
	# 
	# Returns: 
	#	The expected diversity for each element of t.
	
	res <- K * mpfr(exp(d * t), precBits = 500) / (K + mpfr(exp(d * t), precBits = 500) - 1)
	return(as.numeric(res))

}

DD_var <- function(d, K, t, c){
	# Computes the expected diversity given density-dependent diversification rate, time with varying carrying capacities
	#
	# Args:
	# 	t: vector of time
	# 	d: diverisifcation rate
	# 	K: carrying capacity (i.e., strength of density dependence)
	#	c: standardization for varying carrying capacity
	# 
	# Returns: 
	#	The expected diversity for each element of t.
	
	res <- c * K * mpfr(exp(d * t), precBits = 500) / (c * K + mpfr(exp(d * t), precBits = 500) - 1)
	return(as.numeric(res))
}


growth_phase <- function(a, b, t1){
	# Computes the expected diveristy given a diversity-dependent ontogeny model during the growth phase (see DD_ont)
	
	A <- sqrt(a) * t1 / sqrt(2)
	S <- (2*a^0.5*exp(a*t1^2/2))/(sqrt(2*pi) * b * erfi(A) + 2*sqrt(a))
	return(as.numeric(S))
}

decay_phase <- function(S_g, b, r_max, c, t2){
	# Computes the expected diveristy given a diversity-dependent ontogeny model during the decay phase (see DD_ont)
	
	err1 <- (c * t2 - r_max) / sqrt(2*c)
	err2 <- r_max / sqrt(2*c)
	S_numer <- (2*sqrt(c) * S_g * exp(r_max * t2 - (c*(t2^2)/2)))

	S_denom_1 <- sqrt(2*Const("pi")) * b * S_g * exp(mpfr((r_max^2 / (2*c)), precBits = 500))
	S_denom_2 <- erf(mpfr(err1, precBits = 500)) + erf(mpfr(err2, precBits = 500))
	S_denom_3 <- 2*sqrt(c)
	
	S_denom <- (S_denom_1 * S_denom_2) + S_denom_3
	S <- S_numer / S_denom	
	
	S[which(is.na(S) == TRUE)] <- S_g[which(is.na(S) == TRUE)] # if you are not in your decay phase yet, then you get your growth phase value
	
	return(as.numeric(S))
}

	
DD_ont <- function(r_max_0, K_max, t1, t2, p1, rel_maxK, t = NULL){
	# Computes the expected diversity given density-dependent diversification rate and time
	#
	# Args:
	# 	t1: vector of time at assumed maximum carrying capacity of island
	# 	t2: vector of time since formation
	# 	rel_maxK: relative maximum carrying capacity
	# 	p1: current remaining carrying capacity
	#	K_max: maximum carrying capacity
	# 	r_max_0: maximum intitial diverisifcation rate
	# 	t: vector of time
	# 
	# Returns: 
	#	The expected diversity for each element of t.

	rel_r_max <- r_max_0 * rel_maxK 	# Adjusts the maximum diversification rate to be relative to the maximum carrying capacity of the archipelago
	rel_K_max <- K_max * rel_maxK 		# Adjusts the maximum carrying capacity to be a function of the maximum carrying capacity of the archipelago
	a <- rel_r_max / t1
	b <- rel_r_max / rel_K_max
	c <- p1 * rel_r_max / t2
	
	if(is.null(t)){
		S_g <- growth_phase(a = a, b = b, t1 = t1)
		S_d <- decay_phase(S_g = S_g, b = b, r_max = rel_r_max, c = c, t2 = t2)
		return(S_d)
	} else {
		t <- round(t, 2) # R has a very limited floating point representation using seq(), so will just make sure here
		t1_plot <- t[t <= t1]
		t2_plot <- t[t > t1]
		t2_plot <- t2_plot - t1 # calculate time since max carrying capacity
		S_g <- growth_phase(a = a, b = b, t1 = t1_plot)
		S_d <- decay_phase(S_g = S_g[length(S_g)], b = b, r_max = rel_r_max, c = c, t2 = t2_plot)
		return(c(S_g, S_d))
	}
}


## FUNCTIONS CALCULATING THE REALIZED DIVERSIFICATION RATE OVER TIME ========================================================
CR_r <- function(d, t){
	# Returns realized diversification rates under a model of exponential diversification from time = 0 to time t
	#
	# Args:
	# 	d = maximum / intrinsic diversification rate
	# 	t = time
	#
	# Returns:
	# 	A vector of length(t)

	S <- CR(d = d, t = t) # calculate diversity
	return(log(S) /t)
}

DD_r <- function(d, K, t){
	# Returns realized diversification rates under a model of single carrying capacity from time = 0 to time t
	#
	# Args:
	# 	d = maximum / intrinsic diversification rate
	# 	K = maximum carrying capacity
	# 	t = time
	#
	# Returns:
	# 	 A vector of length(t)

	res <- d * (1 - ((K * exp(d*t)) / (K * (K + exp(d * t) - 1))))
	return(as.numeric(res))
}

DD_var_r <- function(K, d, t, c){
	# Returns realized diversification rates under a model of single carrying capacity from time = 0 to time t
	#
	# Args:
	# 	d = maximum / intrinsic diversification rate
	# 	K = maximum carrying capacity
	# 	c = relative carrying capacity
	# 	t = time
	#
	# Returns:
	# 	 A vector of length(t)

	res <- d * (1 - ((c * K * exp(d*t)) / (c * K * (c * K + exp(d * t) - 1))))
	return(as.numeric(res))
}

DD_ont_r <- function(r_max_0, K_max, t1, t2, p1, rel_maxK, t = NULL){
	# Calculates instantaneous diversification rate under an island ontogeny diversity-dependent model
	#
	# Args:
	# 	t1: vector of time at assumed maximum carrying capacity of island
	# 	t2: vector of time since formation
	# 	rel_maxK: relative maximum carrying capacity
	# 	p1: current remaining carrying capacity
	#	K_max: maximum carrying capacity
	# 	r_max_0: maximum intitial diverisifcation rate
	# 
	# Returns: 
	#	The instantaneous diversification rate for each element of t.
  
	rel_r_max <- r_max_0 * rel_maxK 	# Adjusts the maximum diversification rate to be relative to the maximum carrying capacity of the archipelago
	rel_K_max <- K_max * rel_maxK 		# Adjusts the maximum carrying capacity to be a function of the maximum carrying capacity of the archipelago
	a <- rel_r_max / t1[length(t1)]
	b <- rel_r_max / rel_K_max
	c <- p1 * rel_r_max / t2[length(t2)] # in this case, t2 is the time since, see DD_ont to see difference

	if(is.null(t)){
    	t1_plot <- t1
    	t2_plot <- t2 - t1
    	S_g <- growth_phase(a = a, b = b, t1 = t1_plot)
		  S_d <- decay_phase(S_g = S_g, b = b, r_max = rel_r_max, c = c, t2 = t2_plot) # need to take the last value
  	} else {
    	t <- round(t, 2)
    	t1_plot <- t[t <= t1]
    	t2_plot <- t[t > t1]
    	t2_plot <- t2_plot - t1 # in this case, t2 is the time since, see DD_ont to see difference
    	S_g <- growth_phase(a = a, b = b, t1 = t1_plot)
		  S_d <- decay_phase(S_g = S_g[length(S_g)], b = b, r_max = rel_r_max, c = c, t2 = t2_plot) # need to take the last value 
  	}
  
	
	d_g <- a*t1_plot - b*S_g # growth phase diversification rate
	d_d <- rel_r_max - c*t2_plot - b*S_d # decay phase diversification rate
	
 	if(is.null(t)){
    	return(d_d)
  	} else {
    	return(c(d_g, d_d))
  	}
}


## FUNCTIONS PERFORMING MODEL FITTING ========================================================
## NLS REGRESSION (syntax tutorial: http://robinlovelace.net/2013/10/23/nls-demonstation.html)
#(http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Nonlinear-Regression.pdf)

testModels <- function(df, t, t1, t2, c, p1, rel_maxK, z, ...){
	# Uses non-linear regression to fit 3 alternative diversification models
	#
	# Args:
	# 	df: 
	#
	# Returns:
	# 	dataframe with best-fit parameter estimates and model weights
	
	# Define non linear regression algorithm parameters
	control = list(maxiter = 1000, ...)

	print(paste0("Fitting models to ", df$Taxon))
	obs <- data.frame(y = unlist(df[,3:6]))

	# Informed starting parameters
	mean_div <- mean(log(obs$y) / t)  # remember to log N
	mean_obs <- mean(obs$y) # mean species diversity to initialize K

	# Total sum of squares
	TSS <- sum((obs$y - mean(obs$y))^2)

	# Fitting models
	CR_mod <- minpack.lm::nlsLM(y ~ CR(d, t = t), data = obs,
					start = list(d = mean_div),
					lower = c(d = 0.01),
					upper = c(d = 30),
					control = control, trace = F)
	#CR_params <- CR_mod$m$getPars()
	CR_params <- summary(CR_mod)$coefficients
	CR_AIC <- AIC(CR_mod)
	CR_R2 <- 1 - (sum((obs$y - fitted(CR_mod))^2) / TSS)
	

	DD_mod <- minpack.lm::nlsLM(y ~ DD(d, K, t = t), data = obs,
					start = list(d = mean_div, K = mean_obs),
					lower = c(d = 0.01, K = 0.01),
					upper = c(d = 30, K = 300),
					control = control, trace = F)
	#DD_AICc <- AICc(DD_mod)
	DD_AIC <- AIC(DD_mod)
	#DD_params <- DD_mod$m$getPars()
	DD_R2 <- 1 - (sum((obs$y - fitted(DD_mod))^2) / TSS)
	DD_params <- summary(DD_mod)$coefficients
		

	DD_var_mod <- minpack.lm::nlsLM(y ~ DD_var(d, K, t = t, c = c), data = obs,
						start = list(d = mean_div, K = mean_obs),
						lower = c(d = 0.01, K = 0.01),
						upper = c(d = 30, K = 300),
						control = control, trace = F)
	#DD_var_AICc <- AICc(DD_var_mod)
	DD_var_AIC <- AIC(DD_var_mod)
	#DD_var_params <- DD_var_mod$m$getPars()
	DD_var_R2 <- 1 - (sum((obs$y - fitted(DD_var_mod))^2) / TSS)
	DD_var_params <- summary(DD_var_mod)$coefficients
	
	DD_ont_mod <- minpack.lm::nlsLM(y ~ DD_ont(r_max_0, K_max, t1 = t1, t2 = t2,
						p1 = p1, rel_maxK = rel_maxK), data = obs, 
						start = list(r_max_0 = mean_div, K_max = mean_obs),
						lower = c(r_max_0 = 0.01, K_max = 0.01),
						upper = c(r_max_0 = 30, K_max = 300),
						control = control, trace = F)
	#DD_ont_AICc <- AICc(DD_ont_mod)
	DD_ont_AIC 	<- AIC(DD_ont_mod)
	DD_ont_R2 <- 1 - (sum((obs$y - fitted(DD_ont_mod))^2) / TSS)
	DD_ont_params <- summary(DD_ont_mod)$coefficients
	
	# Calculate model akaike weights
	## NOTE: Excluded exponential model as it has one fewer parameter and is hence not comparable
	min_AIC <- min(c(DD_AIC, DD_var_AIC, DD_ont_AIC))
	mod_delta <- c(DD_AIC, DD_var_AIC, DD_ont_AIC) - min_AIC
	mod_weight <- exp(-mod_delta / 2) / sum(exp(-mod_delta/2) )		

	DD_w <- signif(mod_weight[1])
	DD_var_w <- signif(mod_weight[2])
	DD_ont_w <- signif(mod_weight[3])
	
  	AIC <- c(DD_AIC, DD_var_AIC, DD_ont_AIC, CR_AIC)
  	r <- c(DD_params[1,1], DD_var_params[1,1], DD_ont_params[1,1], CR_params[1,1])
  	K <- c(DD_params[2,1], DD_var_params[2,1], DD_ont_params[2,1], NA)
  	r_SE <- c(DD_params[1,2], DD_var_params[1,2], DD_ont_params[1,2], CR_params[1,2])
  	K_SE <- c(DD_params[2,2], DD_var_params[2,2], DD_ont_params[2,2], NA)

  	R2 <- c(DD_R2, DD_var_R2, DD_ont_R2, CR_R2)

  	model <- c("Single K", "Varying K", "Varying K + Ontogeny", "Exponential")
  	w <- c(DD_w, DD_var_w, DD_ont_w, NA)
  
  	res <- data.frame(model, r, r_SE, K, K_SE, AIC, R2, w, z,
  					c_1  = c[1], c_2 = c[2], c_3 = c[3], c_4 = c[4],
  					t1_1 = t1[1], t1_2 = t1[2], t1_3 = t1[3], t1_4 = t1[4],
 				    t2_1 = t2[1], t2_2 = t2[2], t2_3 = t2[3], t2_4 = t2[4],
 				    p1_1 = p1[1], p1_2 = p1[2], p1_3 = p1[3], p1_4 = p1[4],
 				    rel_maxK_1 = rel_maxK[1], rel_maxK_2 = rel_maxK[2],
                    rel_maxK_3 = rel_maxK[3], rel_maxK_4 = rel_maxK[4],
                    S_1 = obs$y[1], S_2 = obs$y[2], S_3 = obs$y[3], S_4 = obs$y[4]) # species diversities
	return(res)
}


testModelList <- function(df, .data = islClades, ...){
  # Convenient wrapper function for fitting diversification models
  # Runs testModels on each row of species numbers in .data
  #
  # Args:
  # 	.data = dataframe of species numbers
  # 	df = dataframe with island ontogeny parameters (see calcOntogeny function)
  #		... = arguments pased onto model-fitting options
  #
  # Returns:
  # 	dataframe of model weights and parameter estimates

  res <- plyr::ddply(.data,
        .variables = .(Taxon), .fun = testModels, 
        t = df$t_total, c = df$rel_currK,
        t1 = df$t1, t2 = df$t2,
        rel_maxK = df$rel_maxK,
        p1 = df$p1, z = df$z[1], ...)
  
  return(res)
}

fitSAR <- function(x, df){
  temp <- subset(df, Taxon == x)
  mod <- lm(logS ~ logcurrA, data = temp)
  return(data.frame("Taxon" = x, "intercept" = coefficients(mod)[1], "slope" = coefficients(mod)[2]))
}

## GRAPHICAL FUNCTIONS ========================================================
CR_DivList <- function(df, t, isl){
  Div <- CR(d = df$r, t = t)
	data.frame(Taxon = df$Taxon, Div, t, isl)
}
DD_DivList <- function(df, t, isl){
	Div <- DD(d = df$r, K = df$K, t = t)
	data.frame(Taxon = df$Taxon, Div, t, isl)
}
DDvar_DivList <- function(df, c, t, isl){
  Div <- DD_var(c = df[, c], d = df$r, K = df$K, t = t)
  data.frame(Taxon = df$Taxon, Div, t, isl)
}
DDont_DivList <- function(df, t1, t2, p1, rel_maxK, isl, t){
  Div <- DD_ont(r_max_0 = df$r, K_max = df$K, t1 = df[,t1], t2 = df[,t2], p1 = df[, p1], rel_maxK = df[, rel_maxK], t = t)
  data.frame(Taxon = df$Taxon, Div, t, isl)
}

predictDivTrajectory <- function(df, model, haw_time = NULL, ma_time = NULL, oa_time = NULL, ka_time = NULL){
  # Produces a plot of species diversities
  # Takes dataframe produced from testModels or testModelList functions
  #
  # Args:
  # 	x: dataframe
  #
  # Returns:
  # 	ggplot2 object plotting relative model akaike weights	
	if(is.null(haw_time)){
 		haw_time <- seq(0.01, df[1, "t2_1"] + df[1, "t1_1"], 0.005)
 	}
 	if(is.null(ma_time)){
 		ma_time <- seq(0.01, df[1, "t2_2"] + df[1, "t1_2"], 0.005)
 	}
 	if(is.null(oa_time)){
 		oa_time <- seq(0.01, df[1, "t2_3"] + df[1, "t1_3"], 0.005)  
 	}
 	if(is.null(ka_time)){
 		ka_time <- seq(0.01, df[1, "t2_4"] + df[1, "t1_4"], 0.005)
 	}
  	if(model == "Exponential"){
  		HA <- ddply(.data = df[df$model == model, ],
  				.variable = .(Taxon), .fun = CR_DivList, t = haw_time, isl = "Hawaii")
  		MA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = CR_DivList, t = ma_time, isl = "Maui Nui")
  	  	OA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = CR_DivList, t = oa_time, isl = "Oahu")
  	  	KA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = CR_DivList, t = ka_time, isl = "Kauai")
  	}
  	if(model == "Single K"){
    	all <- ddply(.data = df[df$model == model, ],
        	        .variable = .(Taxon), .fun = DD_DivList, t = ka_time, isl = "All")
  	}
 	if(model == "Varying K"){
    	HA <- ddply(.data = df[df$model == model, ],
        	        .variable = .(Taxon), .fun = DDvar_DivList, t = haw_time, c = "c_1", isl = "Hawaii")
    	MA <- ddply(.data = df[df$model == model, ],
        	        .variable = .(Taxon), .fun = DDvar_DivList, t = ma_time, c = "c_2", isl = "Maui Nui")
    	OA <- ddply(.data = df[df$model == model, ],
        	        .variable = .(Taxon), .fun = DDvar_DivList, t = oa_time, c = "c_3", isl = "Oahu")
    	KA <- ddply(.data = df[df$model == model, ],
        	        .variable = .(Taxon), .fun = DDvar_DivList, t = ka_time, c = "c_4", isl = "Kauai")
  	}
  	if(model == "Varying K + Ontogeny"){
    	HA <- ddply(.data = df[df$model == model, ],
        	        .variable = .(Taxon), .fun = DDont_DivList, t = haw_time, 
            	    t1 = "t1_1", t2 = "t2_1", p1 = "p1_1", rel_maxK = "rel_maxK_1", isl = "Hawaii")
    	MA <- ddply(.data = df[df$model == model, ],
        	        .variable = .(Taxon), .fun = DDont_DivList, t = ma_time, 
            	    t1 = "t1_2", t2 = "t2_2", p1 = "p1_2", rel_maxK = "rel_maxK_2", isl = "Maui Nui")
    	OA <- ddply(.data = df[df$model == model, ],
            	    .variable = .(Taxon), .fun = DDont_DivList, t = oa_time, 
        	        t1 = "t1_3", t2 = "t2_3", p1 = "p1_3", rel_maxK = "rel_maxK_3", isl = "Oahu")
    	KA <- ddply(.data = df[df$model == model, ],
        	        .variable = .(Taxon), .fun = DDont_DivList, t = ka_time, 
            	    t1 = "t1_4", t2 = "t2_4", p1 = "p1_4", rel_maxK = "rel_maxK_4", isl = "Kauai")
  	}
	if(model %in% c("Varying K", "Varying K + Ontogeny", "Exponential")){
		Div <- rbind(HA, MA, OA, KA)
		Div$isl <- factor(Div$isl, levels = c("Hawaii", "Maui Nui", "Oahu", "Kauai"))	
	} else {
		Div <- all
	}
	
	return(Div)
}
Naive_DivRateList <- function(df, t, isl){
	# Calculates naive diversification rate, essentially number of species divided by age of island
	#
	if(isl == "Ha"){
		S <- df$S_1
	}
	if(isl == "Ma"){
		S <- df$S_2
	}
	if(isl == "Oa"){
		S <- df$S_3
	}
	if(isl == "Ka"){
		S <- df$S_4
	}
	DivRate <- rep(log(S) / t[length(t)], length(t))
	data.frame(Taxon = df$Taxon, DivRate, t, isl)
}

CR_DivRateList <- function(df, t, isl){
	DivRate <- CR_r(d = df$r, t = t)
	data.frame(Taxon = df$Taxon, DivRate, t, isl)
}
DD_DivRateList <- function(df, t, isl){
	DivRate <- DD_r(d = df$r, K = df$K, t = t)
	data.frame(Taxon = df$Taxon, DivRate, t, isl)
}

DDvar_DivRateList <- function(df, c, t, isl){
  DivRate <- DD_var_r(c = df[, c], d = df$r, K = df$K, t = t)
  data.frame(Taxon = df$Taxon, DivRate, t, isl)
}
DDont_DivRateList <- function(df, t, t1, t2, p1, rel_maxK, isl){
  DivRate <- DD_ont_r(r_max_0 = df$r, K_max = df$K, t1 = df[,t1], t2 = df[,t2], p1 = df[, p1], rel_maxK = df[, rel_maxK], t = t)
  data.frame(Taxon = df$Taxon, DivRate, t, isl)
}
predictDivRateTrajectory <- function(df, model, haw_time = NULL, ma_time = NULL, oa_time = NULL, ka_time = NULL){
  #	Generate island-specific realized diversification rate trajectories for different taxa under a specified diversification-rate model
  # Uses parameter estimates generated from TestModelList()
  #
  # Args:
  #		df: results data.frame generate TestModelList()
  # 	model: diversification rate model, i.e., "Single K", "Varying K" or "Varying K + Ontogeny"
  #
  # Returns:
  # 	data.frame containing variables diversification rate ("DivRate"), taxon ("Taxon"), time ("t") and island ("isl")
  if(is.null(haw_time)){
    haw_time <- seq(0.01, df[1, "t2_1"] + df[1, "t1_1"], 0.005)
  }
  if(is.null(ma_time)){
    ma_time <- seq(0.01, df[1, "t2_2"] + df[1, "t1_2"], 0.005)
  }
  if(is.null(oa_time)){
    oa_time <- seq(0.01, df[1, "t2_3"] + df[1, "t1_3"], 0.005)  
  }
  if(is.null(ka_time)){
    ka_time <- seq(0.01, df[1, "t2_4"] + df[1, "t1_4"], 0.005)
  }
  
  if(model == "Naive"){
   	HA <- ddply(.data = df[df$model == "Single K", ], # anything model will do since method only needs observed species numbers
    	        .variable = .(Taxon), .fun = Naive_DivRateList, t = haw_time, isl = "Ha")
    MA <- ddply(.data = df[df$model == "Single K", ],
        	    .variable = .(Taxon), .fun = Naive_DivRateList, t = ma_time, isl = "Ma")
    OA <- ddply(.data = df[df$model == "Single K", ],
        	    .variable = .(Taxon), .fun = Naive_DivRateList, t = oa_time, isl = "Oa")
    KA <- ddply(.data = df[df$model == "Single K", ],
        	    .variable = .(Taxon), .fun = Naive_DivRateList, t = ka_time, isl = "Ka")
  }
  if(model == "Exponential"){
  	HA <- ddply(.data = df[df$model == model, ],
  				.variable = .(Taxon), .fun = CR_DivRateList, t = haw_time, isl = "Ha")
  	MA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = CR_DivRateList, t = ma_time, isl = "Ma")
    OA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = CR_DivRateList, t = oa_time, isl = "Oa")
    KA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = CR_DivRateList, t = ka_time, isl = "Ka")
  }
  if(model == "Single K"){
  	HA <- ddply(.data = df[df$model == model, ],
  				.variable = .(Taxon), .fun = DD_DivRateList, t = haw_time, isl = "Ha")
  	MA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DD_DivRateList, t = ma_time, isl = "Ma")
    OA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DD_DivRateList, t = oa_time, isl = "Oa")
    KA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DD_DivRateList, t = ka_time, isl = "Ka")
  }
  if(model == "Varying K"){
    HA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DDvar_DivRateList, t = haw_time, c = "c_1", isl = "Ha")
    MA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DDvar_DivRateList, t = ma_time, c = "c_2", isl = "Ma")
    OA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DDvar_DivRateList, t = oa_time, c = "c_3", isl = "Oa")
    KA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DDvar_DivRateList, t = ka_time, c = "c_4", isl = "Ka")
  }
  if(model == "Varying K + Ontogeny"){
    HA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DDont_DivRateList, t = haw_time, 
                t1 = "t1_1", t2 = "t2_1", p1 = "p1_1", rel_maxK = "rel_maxK_1", isl = "Ha")
    MA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DDont_DivRateList, t = ma_time, 
                t1 = "t1_2", t2 = "t2_2", p1 = "p1_2", rel_maxK = "rel_maxK_2", isl = "Ma")
    OA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DDont_DivRateList, t = oa_time, 
                t1 = "t1_3", t2 = "t2_3", p1 = "p1_3", rel_maxK = "rel_maxK_3", isl = "Oa")
    KA <- ddply(.data = df[df$model == model, ],
                .variable = .(Taxon), .fun = DDont_DivRateList, t = ka_time, 
                t1 = "t1_4", t2 = "t2_4", p1 = "p1_4", rel_maxK = "rel_maxK_4", isl = "Ka")
  }

  DivRate <- rbind(HA, MA, OA, KA)
  DivRate$isl <- factor(DivRate$isl, levels = c("Ha", "Ma", "Oa", "Ka"))
  return(DivRate)
}


g_legend <- function(gg){
	# Extracts a legend from a ggplot2 object
	#
	# Args:
	# 	gg: ggplot2 object with a legend
	# 
	# Returns:
	# 	ggplot2 object
    g <- ggplotGrob(gg)$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    return(legend)
}


theme_hawaii <- function(){
	# Theme for diversification rate plot

    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.spacing = unit(2, "lines"), # margin between facet_wrap panels
          strip.background = element_blank(), # facet_wrap background
          strip.text.x = element_blank(), # facet_wrap text
          axis.title = element_text(colour = "grey50", size = 15),
          axis.text = element_text(colour = "grey50", size = 12),
          axis.line = element_line(colour = "grey50"),
          axis.line.x = element_blank())
}


theme_divplot <- function(base_size = 12, base_family = ""){
    theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.spacing = unit(2, "lines"), # margin between facet_wrap panels
        strip.background = element_blank(), # facet_wrap background
	    strip.text.x = element_blank(), # facet_wrap text
        axis.title = element_text(colour = "grey50", size = 15),
        axis.text = element_text(colour = "grey50", size = 12),
        axis.line = element_line(colour = "grey50"))
}