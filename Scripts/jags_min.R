sink("Scripts/Jags/jags_min.txt")
cat("
  model {
    # Prior distributions for basic parameters
    # Intercepts
    beta.a0 ~ dnorm(0,0.01)    # intercept for availability
    alpha0 ~ dnorm(0, 0.01)    # intercept for sigma
    # alpha1 ~ dnorm(0,0.01)     # slope on sigma covariate
    
    # Coefficients
    # beta.a1 ~ dnorm(0,0.01)  # slope for availability covariate
    # beta.a2 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0,0.01)      # intercept for lambda
    # beta1 ~ dnorm(0,0.01)       # slope for lambda covariate
    # beta2 ~ dnorm(0, 0.01)
    
    for(s in 1:nsites){
      # Add covariates to scale parameter DISTANCE (perceptibility)
      log(sigma[s]) <- alpha0
      
      # Add covariates for availability here TIME-REMOVAL (availability)
      p.a[s] <- exp(beta.a0) / (1+exp(beta.a0)) 
      # Optional covariates on availability
      # p.a[s] <- exp(beta.a0 + eps.time[point[s]] + beta.a1*day[s])/(1+exp(beta.a0 + eps.time[point[s]] + beta.a1*day[s]))
      
      # Distance sampling detection probability model
      for(b in 1:nD){
        log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  # Half-normal  
        f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) # Radial density function
        pi.pd[b,s] <- g[b,s]*f[b,s]  #  Product Pr(detect)*Pr(distribution)
        pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]  # Conditional probabilities
      }
      pdet[s] <- sum(pi.pd[,s])  # Probability of detection at all 
      
      # Time-removal probabilities
      pi.pa[1, s] <- 1 - pow((1 - p.a[s]), 3)
      pi.pa[2, s] <- (1 - pow((1 - p.a[s]), 2)) * pow((1 - p.a[s]), 3)
      pi.pa[3, s] <- (1 - pow((1 - p.a[s]), 5)) * pow((1 - p.a[s]), 5)
      
      for (k in 1:K){
        # pi.pa[k,s] <- p.a[s] * pow(1-p.a[s], (k-1)) # assumes equal length time bins  
        pi.pa.c[k,s] <- pi.pa[k,s]/phi[s] # Conditional probabilities of availability
      }
      phi[s] <- sum(pi.pa[,s]) # Probability of ever available
    }
    # Conditional observation model for categorical covariates
    for(i in 1:nobs){  
      dclass[i] ~ dcat(pi.pd.c[,site[i]]) 
      tint[i] ~ dcat(pi.pa.c[,site[i]])
    }
    # Abundance model
    for(s in 1:nsites){ 
      # Binomial model for # of captured individuals
      # n[s] ~ dbin(pmarg[s], M[s]) # Formulation b, see text
      # pmarg[s] <- pdet[s]*phi[s]
      n[s] ~ dbin(pdet[s], N[s])    # Formulation a, see text
      N[s] ~ dbin(phi[s], M[s])      # Number of available individuals
      M[s] ~ dpois(lambda[s])       # Abundance per survey/site/point
      # Add site-level covariates to lambda
      log(lambda[s]) <- beta0 #+ eps.lam[point[s]] + beta1*habitat[s] 
    }
    # Derived quantities
    # Mtot <- sum(M[])  # Total population size
    # Ntot <- sum(N[])  # Total available population size
    # PDETmean <- mean(pdet[]) # Mean perceptibility across sites
    # PHImean <- mean(phi[]) # Mean availability across sites
  }
  ", fill=TRUE)
sink()