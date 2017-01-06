sink("Scripts/Jags/tr-ds.txt")
cat("
  model {
    # Prior distributions for basic parameters
    # Intercepts
    beta.a0 ~ dnorm(0,0.01)    # intercept for availability
    alpha0 ~ dunif(0, 20) # dnorm(0, 0.0001)    # intercept for sigma
    beta0 ~ dnorm(0,0.01)       # intercept for lambda
    
    # Coefficients
    beta.a1 ~ dnorm(0,0.01)     # slope for availability day of the year
    beta.a2 ~ dnorm(0,0.01)     # slope for availability time of day
    beta.a3 ~ dnorm(0,0.01)     # slope for availability time of day^2
    #  beta.a4 ~ dnorm(0, 0.01)    # effect of year 2 on availability
    beta1 ~ dnorm(0,0.01)       # slope for lambda covariate: veg
    beta2 ~ dnorm(0,0.01)       # slope for lambda covariate: shrub
    beta3 ~ dnorm(0,0.01)       # slope for lambda covariate: tree
    beta4 ~ dnorm(0, 0.01)      # effect of year 2 on abundance
    alpha1 ~ dunif(-5, 5)     # slope on sigma covariate: tree
    # alpha2 ~ dunif(-5, 5)     # slope on sigma covariate: visit 2
    # alpha3 ~ dunif(-5, 5)     # slope on sigma covariate: visit 3
    # alpha4 ~ dunif(-5, 5)     # slope on sigma covariate: visit 4
    # alpha2 ~ dnorm(0, 0.001)    # effect of year 2 on detection (distance) - maybe not necessary
    
    # random point abundance  
    for(p in 1:npoints) {
    eps.lam[p] ~ dnorm(0, tau.lam)
    }
    sigma.lam ~ dunif(0, 2)
    tau.lam <- 1 / (sigma.lam * sigma.lam)
    
    # random overdispersed detection
    for(s in 1:nsites) {
    eps.dist[s] ~ dnorm(0, tau.dist)
    }
    sigma.dist ~ dunif(0, 10)
    tau.dist <- 1 / (sigma.dist * sigma.dist)
    
    # random point availability
    for(p in 1:npoints) {
    eps.time[p] ~ dnorm(0, tau.time)
    }
    sigma.time ~ dunif(0, 2)
    tau.time <- 1 / (sigma.time * sigma.time)
    
    for(s in 1:nsites){
    # Add covariates to scale parameter DISTANCE (perceptibility)
    log(sigma[s]) <- alpha0 + alpha1*tree[s] # + eps.dist[siteVisit[s]] # + alpha2*visit2[s] + alpha3*visit3[s] + alpha4*visit4[s] # + eps.dist[point[s]]
    # sigma[s] ~ dunif(0, 14.9) # must be constrained to <15 if distance not standardized
    
    # Add covariates for availability here TIME-REMOVAL (availability)
    # p.a[s] <- exp(beta.a0) / (1+exp(beta.a0)) 
    # Optional covariates on availability
    mu.avail[s] <- beta.a0 + beta.a1*day[s] + beta.a2*time[s] + beta.a3*time[s]*time[s] # + beta.a4*year[s] + eps.time[point[s]]
    p.a[s] <- exp(mu.avail[s])/(1+exp(mu.avail[s]))
    
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
    pi.pa.c[k,s] <- pi.pa[k,s]/pavail[s] # Conditional probabilities of availability
    }
    pavail[s] <- sum(pi.pa[,s]) # Probability of ever available
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
    # pmarg[s] <- pdet[s]*pavail[s]
    n[s] ~ dbin(pdet[s], N[s])    # Formulation a, see text
    N[s] ~ dbin(pavail[s],M[s])      # Number of available individuals
    M[s] ~ dpois(lambda[s])       # Abundance per survey/site/point
    # Add site-level covariates to lambda
    log(lambda[s]) <- beta0 + eps.lam[point[s]] + beta1*veg[s] + beta2*shrub[s] + beta4*year[s] # + beta3*tree[s] 
    }
    
    ######## Goodness of fit tests
    for(k in 1:nsites){
    navail.fit[k] ~ dbin(pavail[k], N[k]) # create new realization of model
    y.fit[k] ~ dbin(pdet[k], N[k]) # create new realization of model
    
    e.pd[k]<- pdet[k]*N[k] # original model prediction
    E.pd[k]<- pow(( n[k]- e.pd[k]),2)/(e.pd[k]+0.5)
    E.New.pd[k]<- pow((y.fit[k]-e.pd[k]),2)/(e.pd[k]+0.5)
    
    e.pa[k]<- pavail[k]*N[k] # original model prediction
    E.pa[k]<- pow(( N[k]- e.pa[k]),2)/(e.pa[k]+0.5)
    E.New.pa[k]<- pow((navail.fit[k]-e.pa[k]),2)/(e.pa[k]+0.5)
    }
    fit.pd<- sum(E.pd[])
    fit.new.pd<- sum(E.New.pd[])
    
    fit.pa<- sum(E.pa[])
    fit.new.pa<- sum(E.New.pa[])
    
    ######## Summary stats
    # meanpavail<-mean(pavail[]) # mean probability of availability
    # meanpdet<-mean(pdet[]) # mean probability of perceptibility
    bayesp.pd<-step(fit.new.pd-fit.pd) # Bayesian p-value for perceptibility model
    bayesp.pa<-step(fit.new.pa-fit.pa) # Bayesian p-value for availability model
    meanN<-mean(M[]) # mean site-level abundance
    # totN<-sum(N[])  # population size of total area surveyed
    meansig<-mean(sigma[]) # mean scale parameter across sites
    dens<-meanN/(maxd*maxd*3.14159) # density of birds per ha
    
    # Derived quantities
    Mtot <- sum(M[])  # Total population size
    Ntot <- sum(N[])  # Total available population size
    PDETmean <- mean(pdet[]) # Mean perceptibility across sites
    PAVAILmean <- mean(pavail[]) # Mean availability across sites
    }
    ", fill=TRUE)
sink()