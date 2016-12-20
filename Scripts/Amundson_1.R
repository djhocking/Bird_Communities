## Supplementary file 2 of Amundson et al. AUK
# C. Amundson camundson@usgs.gov 10_17_13

# BUGS/JAGS code for combined hierarchical distance-sampling and time-removal N-mixture model for 
# point-transect counts.

# Note: Model as written works only in JAGS (R2jags)

###### Variables
# y = count of birds per point
# surveyid = survey point/site ID for each individual observed
# dclass = distance class per observation
# tinterval = time interval per observation
# nsurveys = number of points surveyed
# nobs = total count of individuals observed (sum(y))
# delta = distance width for each bin (set to be equal intervals in this example)
# nbreaks = number of distance bins
# mdpts = midpoints of distance bins
# maxd = maximum truncated distance (300 m in this example)
# J = number of time periods
# tree, grass, ag , wet = standardized covariates representing proportion of each habitat within a 300-m radius buffer around each point
# ntrans = total number of transects
# tran = transect ID for each point
# date = Julian date of survey
# Nst = initial value for estimate of N - must be close to N or model will not run

###### Index
# k in 1:nsurveys # surveys
# b in 1:nbreaks # distance bins
# j in 1:J # time intervals
# i in 1:nobs # observations/individuals
# t in 1:ntrans # transects

###### Parameters to estimate
# meansig = mean scale parameter across sites (half normal shape in this example)
# meanpdet = mean probability of perceptibility
# meanpavail = mean probability of availability
# beta.a0 = intercept for availability
# beta.a1 = coefficient for 'date' in availability model
# sigma.0 = intercept for perceptibility
# beta.p1 = coefficient for 'tree' in perceptibility model
# beta1 = coefficient for 'tree' in abundance model
# beta2 = coefficient for 'ag' in abundance model 
# beta3 = coefficient for 'wet' in abundance model 
# beta4 = coefficient for 'grass' in abundance model 
# beta5 = coefficient for 'wet^2' in abundance model 
# meanN = mean site-level abundance
# mu.tran = mean abundance intercept across transects
# sd.tran = SD of random transect effect
# totN = population size of total area surveyed
# bayesp.pd = Bayesian p-value for pd model
# bayesp.pa = Bayesian p-value for pa model
# dens = density of birds per hectare = totN/area surveyed

sink("sim_model.txt")
cat("
    model {
    
    #PRIORS for fixed detection parameters
    # intercepts
    beta.a0~dnorm(0,0.01)
    sigma.0~dunif(0,500) # to upper bound of sigma (pd ~ 1)  
    # coefficients
    beta.a1~dnorm(0,0.01)
    beta.p1~dnorm(0,0.01)
    beta1~dnorm(0,0.01)
    beta2~dnorm(0,0.01)
    beta3~dnorm(0,0.01)
    beta4~dnorm(0,0.01)
    beta5~dnorm(0,0.01)
    
    # random transect (spatial autocorrelation) prior on abundance intercept - nested point within transect
    for(t in 1:ntrans){
    beta0[t] ~ dnorm(mu.tran, tau.tran) 
    beta.tran[t]<-beta0[t]-mu.tran
    }
    mu.tran~ dnorm(0,0.01)
    tau.tran<-pow(sd.tran,-2)
    sd.tran~dunif(0,2)
    
    ##DETECTION PROBABILITY FUNCTIONS  
    for(k in 1:nsurveys){ 
    # add covariates to scale parameter  DISTANCE (perceptibility)
    log(sigma[k])  <- log(sigma.0) + beta.p1*tree[k] 
    # add covariates for availability here TIME-REMOVAL (availability)
    p.a[k]<-exp(beta.a0+beta.a1*date[k])/(1+exp(beta.a0+beta.a1*date[k])) 
    # manual logit above to avoid BUGS issues with logit function
    
    ######## Distance sampling detection probability estimation
    # Using summation technique - Pr(p of x)=exp(-x^2/2*sigma^2)*f(x)
    for(b in 1:nbreaks){
    log(g[b,k])<- -mdpts[b]*mdpts[b]/(2*sigma[k]*sigma[k])  # half-normal detection function - first half of eq., 
    f[b,k]<-  ( 2*mdpts[b]*delta )/(maxd*maxd) # this is f(x), the scaled radial density function
    ##ADD [b] TO DELTA IF INTERVALS ARE NOT ALL THE SAME AMONG BREAKS
    
    pi.pd[b,k]<- g[b,k]*f[b,k]  #this is the product Pr(detect)*Pr(distribution)
    pi.pd.c[b,k]<- pi.pd[b,k]/pdet[k]  # standardizing based on overall capture probability - conditional formulation
    }
    
    pdet[k]<-sum(pi.pd[,k])  # probability of detection is the sum of all rectangular areas
    
    ######## Time-removal detection probability estimation
    
    for (j in 1:J){
    pi.pa[j,k] <- p.a[k] * pow(1-p.a[k], (j-1)) # see salamander example, Royle and Dorazio (2008)
    pi.pa.c[j,k]<- pi.pa[j,k]/pavail[k] # standardizing based on overall availability - conditional formulation
    }
    pavail[k]<-sum(pi.pa[,k]) # probability of capture is the sum of all time intervals
    
    }
    
    ######## Observation-level model  
    for(i in 1:nobs){  
    #single binomial trial with categorical distribution linking distance class and time interval to survey point
    dclass[i] ~ dcat(pi.pd.c[,surveyid[i]]) 
    tinterval[i] ~ dcat(pi.pa.c[,surveyid[i]])
    }
    
    ######## Abundance estimation    
    for(k in 1:nsurveys){ 
    
    # binomial model for # of captured individuals
    y[k]~ dbin(pdet[k],navail[k]) # counts related to probability of detection, given availability
    navail[k]~dbin(pavail[k],N[k]) # probability of being available, given total population size available for sampling
    
    ## abundance model
    N[k]~dpois(lambda[k])# predicted abundance per survey/site/point
    
    # Add site-level covariates to lambda
    log(lambda[k])<- beta0[tran[k]] + beta1*tree[k] + beta2*ag[k] + beta3*grass[k] + beta4*wet[k] +beta5*wet[k]*wet[k]
    
    }
    
    ######## Goodness of fit tests
    for(k in 1:nsurveys){
    
    navail.fit[k] ~ dbin(pavail[k],N[k]) # create new realization of model
    y.fit[k] ~ dbin(pdet[k],navail[k]) # create new realization of model
    
    e.pd[k]<- pdet[k]*navail[k] # original model prediction
    E.pd[k]<- pow(( y[k]- e.pd[k]),2)/(e.pd[k]+0.5)
    E.New.pd[k]<- pow((y.fit[k]-e.pd[k]),2)/(e.pd[k]+0.5)
    
    e.pa[k]<- pavail[k]*N[k] # original model prediction
    E.pa[k]<- pow(( navail[k]- e.pa[k]),2)/(e.pa[k]+0.5)
    E.New.pa[k]<- pow((navail.fit[k]-e.pa[k]),2)/(e.pa[k]+0.5)
    }
    fit.pd<- sum(E.pd[])
    fit.new.pd<- sum(E.New.pd[])
    
    fit.pa<- sum(E.pa[])
    fit.new.pa<- sum(E.New.pa[])
    
    ######## Summary stats
    meanpavail<-mean(pavail[]) # mean probability of availability
    meanpdet<-mean(pdet[]) # mean probability of perceptibility
    bayesp.pd<-step(fit.new.pd-fit.pd) # Bayesian p-value for perceptibility model
    bayesp.pa<-step(fit.new.pa-fit.pa) # Bayesian p-value for availability model
    meanN<-mean(N[]) # mean site-level abundance
    totN<-sum(N[])  # population size of total area surveyed
    meansig<-mean(sigma[]) # mean scale parameter across sites
    dens<-meanN/(maxd*maxd*3.14159/10000) # density of birds per ha
    
    }
    
    
    ", fill=TRUE)
sink()

# End code.
