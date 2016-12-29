sink("jags_min.txt")
cat("
  model {
    
    # PRIORS for fixed detection parameters
    beta.a0 ~ dnorm(0,0.01)
    sigma.0 ~ dunif(0.01, 200)
    beta0 ~ dnorm(0, 0.01)
    
    # DETECTION PROBABILITY FUNCTIONS  
    for(k in 1:nsurveys){ 
      
      # add covariates to scale parameter  DISTANCE (perceptibility)
      log(sigma[k]) <- log(sigma.0)
      
      # add covariates for availability here TIME-REMOVAL (availability)
      p.a.mu[k] <- beta.a0
      
      p.a[k] <- exp(p.a.mu[k]) / (1 + exp(p.a.mu[k])) # manual logit above to avoid BUGS issues with logit function
      
      ######## Distance sampling detection probability estimation
      # Using summation technique - Pr(p of x)=exp(-x^2/2*sigma^2)*f(x)
      for(b in 1:nbreaks) {
        log(g[b,k]) <- -mdpts[b] * mdpts[b] / (2 * sigma[k] * sigma[k])  # half-normal detection function - first half of eq., 
        f[b,k] <- (2 * mdpts[b] * delta ) / (maxd * maxd) # this is f(x), the scaled radial density function
        ##ADD [b] TO DELTA IF INTERVALS ARE NOT ALL THE SAME AMONG BREAKS
        
        pi.pd[b,k] <- g[b,k] * f[b,k]  # this is the product Pr(detect)*Pr(distribution)
        pi.pd.c[b,k] <- pi.pd[b,k] / pdet[k]  # standardizing based on overall capture probability - conditional formulation
      }
      
      pdet[k] <- sum(pi.pd[,k])  # probability of detection is the sum of all rectangular areas
      
      ######## Time-removal detection probability estimation
      
      # Removal detection prob for each of 3 time periods - currently hard coded for 3,2, 5 min intervals
      pi.pa[1, k] <- 1 - pow((1 - p.a[k]), 3)
      pi.pa[2, k] <- (1 - pow((1 - p.a[k]), 2)) * pow((1 - p.a[k]), 3)
      pi.pa[3, k] <- (1 - pow((1 - p.a[k]), 5)) * pow((1 - p.a[k]), 5)
      
      for (j in 1:J) {
        pi.pa.c[j,k] <- pi.pa[j,k] / pavail[k] # standardizing based on overall availability - conditional formulation
      }
      pavail[k] <- sum(pi.pa[,k]) # probability of capture is the sum of all time intervals
      # 
    } # End detection loop
    
    ######## Observation-level model  
    for(i in 1:nobs) {  
      # single binomial trial with categorical distribution linking distance class and time interval to survey point
      dclass[i] ~ dcat(pi.pd.c[ , surveyid[i]]) 
      tinterval[i] ~ dcat(pi.pa.c[ , surveyid[i]])
    }
    
    ######## Abundance estimation    
    for(k in 1:nsurveys) { 
      # binomial model for # of captured individuals
      y[k] ~ dbin(pdet[k], navail[k]) # counts related to probability of detection, given availability
      navail[k] ~ dbin(pavail[k], N[k]) # probability of being available, given total population size available for sampling
      
      ## abundance model
      N[k] ~ dpois(lambda[k]) # predicted abundance per survey/site/point
      
      # Add site-level covariates to lambda
      log(lambda[k]) <- beta0
    }
    
  }
  
  ", fill=TRUE)
sink()