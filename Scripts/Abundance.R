###############################################################################
# Estimate Abundance of birds from point counts with dist and removal time
# Daniel J. Hocking
# Collaborators: Meta Griffin, Frank Ammer
# October 2016
###############################################################################

# Load Packages
library(dplyr)
library(lubridate)
library(tidyr)
library(readr)
library(jagsUI)

# Load Data
df_counts <- read_csv("Data/R_2015_pointcount.csv")
str(df_counts)

sink("jags_model.txt")
cat("
  model {
    
    # PRIORS for fixed detection parameters
    beta.a0 ~ dnorm(0,0.01)
    
    # DETECTION PROBABILITY FUNCTIONS  
    for(k in 1:nsurveys){ 
      # add covariates to scale parameter  DISTANCE (perceptibility)
      # log(sigma[k]) <- log(sigma.0) + beta.p1*tree[k] 
      sigma[k] ~ dunif(0, 100)
      # add covariates for availability here TIME-REMOVAL (availability)
      p.a[k] <- exp(beta.a0) / (1 + exp(beta.a0)) 
      # manual logit above to avoid BUGS issues with logit function
      
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
      
      for (j in 1:J) {
        pi.pa[j,k] <- p.a[k] * pow(1-p.a[k], (j-1)) # see salamander example, Royle and Dorazio (2008)
        pi.pa.c[j,k] <- pi.pa[j,k] / pavail[k] # standardizing based on overall availability - conditional formulation
      }
      pavail[k] <- sum(pi.pa[,k]) # probability of capture is the sum of all time intervals
      
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
      # log(lambda[k]) <- beta0[tran[k]] + beta1*tree[k] + beta2*ag[k] + beta3*grass[k] + beta4*wet[k] + beta5*wet[k]*wet[k]
      lambda[k] ~ dunif(0, 1000)
    }
    
    ######## Goodness of fit tests
    # for(k in 1:nsurveys){
    #   navail.fit[k] ~ dbin(pavail[k],N[k]) # create new realization of model
    #   y.fit[k] ~ dbin(pdet[k],navail[k]) # create new realization of model
    #   
    #   e.pd[k]<- pdet[k]*navail[k] # original model prediction
    #   E.pd[k]<- pow(( y[k]- e.pd[k]),2)/(e.pd[k]+0.5)
    #   E.New.pd[k]<- pow((y.fit[k]-e.pd[k]),2)/(e.pd[k]+0.5)
    #   
    #   e.pa[k]<- pavail[k]*N[k] # original model prediction
    #   E.pa[k]<- pow(( navail[k]- e.pa[k]),2)/(e.pa[k]+0.5)
    #   E.New.pa[k]<- pow((navail.fit[k]-e.pa[k]),2)/(e.pa[k]+0.5)
    # }
    # fit.pd<- sum(E.pd[])
    # fit.new.pd<- sum(E.New.pd[])
    # 
    # fit.pa<- sum(E.pa[])
    # fit.new.pa<- sum(E.New.pa[])
    # 
    # ######## Summary stats
    # meanpavail<-mean(pavail[]) # mean probability of availability
    # meanpdet<-mean(pdet[]) # mean probability of perceptibility
    # bayesp.pd<-step(fit.new.pd-fit.pd) # Bayesian p-value for perceptibility model
    # bayesp.pa<-step(fit.new.pa-fit.pa) # Bayesian p-value for availability model
    # meanN<-mean(N[]) # mean site-level abundance
    # totN<-sum(N[])  # population size of total area surveyed
    # meansig<-mean(sigma[]) # mean scale parameter across sites
    # dens<-meanN/(maxd*maxd*3.14159/10000) # density of birds per ha
    
  }
    
    ", fill=TRUE)
sink()


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


# Bundle data

# tinterval = 1,2, or 3
# dclass = 1, 2, or 3

# example may assume equal time intervals

jags_data<-list(y=y,
               surveyid=as.numeric(surveyid),
               dclass=as.numeric(dclass),
               nsurveys=nsites,
               nobs=sum(y, na.rm = TRUE),
               delta=c(50, 50, 100),
               nbreaks=3,
               mdpts=c(25, 75, 150),
               maxd=200,
               J=max(unique(tinterval)),
               tinterval=as.numeric(tinterval),
               )

# create initial values for N and navail that are very close to actual values or model will not run!
Nst <- jags_data$y + 1

# Inits function
inits <- function(){list(N=Nst,navail=Nst,sigma.0=runif(1,120,150), beta0=runif(ntrans,-1,1),beta.a1=runif(1,-1,1),
                         beta1=runif(1,-1,1),beta2=runif(1,-1,1),mu.tran=runif(1,-1,1),sd.tran=runif(1,0,2),
                         beta3=runif(1,-1,1),beta4=runif(1,-1,1),beta5=runif(1,-1,1),beta.p1=runif(1,-1,1),
                         beta.a0=runif(1,-1,1))} #

# parameters to estimate
# careful printing pavail, pdet and N - will have nsites (e.g., in this example 100) values
# params<-c("meansig","meanpdet","meanpavail","beta.a0","beta.a1","sigma.0","beta.p1","beta1","beta2","beta3","beta4","beta5","meanN","mu.tran","sd.tran","totN","bayesp.pa","bayesp.pd","beta0","beta.tran","N")
params<-c("beta.a0","sigma.0","N")

# MCMC settings
# pavail can be subject to poor mixing in field data - keep thin high, burn-in long, and conduct sufficient number of iterations
nc<-3
ni<-100
nb<-1
nt<-1

## ONLY WORKS IN JAGS
# A bug fix for JAGS - model may produce error without this fix
set.factory("bugs::Conjugate", FALSE, type="sampler")

sim.fit<-jags(data=win.data, inits=inits, parameters.to.save=params, model.file="sim_model.txt",
              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)



