sink("jags_full_2015.txt")
  cat("
      model {
      
      # PRIORS for fixed detection parameters
      beta.a0 ~ dnorm(0,0.01)
        beta.a1 ~ dnorm(0,0.01)
        beta.a2 ~ dnorm(0,0.01)
        beta.a3 ~ dnorm(0,0.01)
        beta.a4 ~ dnorm(0,0.01)
       ln.beta0 ~ dnorm(0,0.001)
 # beta0 ~ dunif(log(1), log(10000))
beta0 <- exp(ln.beta0)
        beta1 ~ dnorm(0,0.01)
        beta2 ~ dnorm(0,0.01)
        beta3 ~ dnorm(0,0.01)
        beta4 ~ dnorm(0,0.01)
        beta.p1 ~ dnorm(0,0.01)
 sigma.0~dunif(0, 200)

# random effect of location (Point)
for(l in 1:n_points) {
  eps.n[l] ~ dnorm(0, tau.eps.n)
}
tau.eps.n <- 1 / (sigma.eps.n * sigma.eps.n)
sigma.eps.n ~ dunif(0, 10)
      
      # DETECTION PROBABILITY FUNCTIONS  
      for(k in 1:nsurveys){ 

      # add covariates to scale parameter  DISTANCE (perceptibility)
       log(sigma[k]) <- log(sigma.0) + beta.p1*VegHgt[point[k]] # + beta.p1*tree[k] 
      # sigma[k] ~ dunif(0, 100)

      # add covariates for availability here TIME-REMOVAL (availability)
        p.a.mu[k] <- beta.a0 + beta.a1*day[k] + beta.a3*time[k] # day of year, time of day, weather, wind (hard because categorical) beta.a2*day[k]*day[k] + + beta.a4*time[k]*time[k] 

      p.a[k] <- exp(p.a.mu[k]) / (1 + exp(p.a.mu[k])) 
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
      
      # Removal detection prob for each of 3 time periods - currently hard coded for 3,2, 5 min intervals
pi.pa[1, k] <- 1 - pow((1 - p.a[k]), 3)
pi.pa[2, k] <- (1 - pow((1 - p.a[k]), 2)) * pow((1 - p.a[k]), 3)
pi.pa[3, k] <- (1 - pow((1 - p.a[k]), 5)) * pow((1 - p.a[k]), 5)

        for (j in 1:J) {
      # pi.pa[j,k] <- p.a[k] * pow(1-p.a[k], (j-1)) # if intervals are all the same length of time
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
      log(lambda[k]) <- exp(ln.beta0) + beta3*Shrub_stm_total[k] + beta4*Tree_stm_total[k] # + eps.n[k] + beta1*VegHgt[k] + beta2*VegHgt[k]*VegHgt[k]
       # lambda[k] ~ dunif(0, 10000) 
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
      # # bayesp.pd<-step(fit.new.pd-fit.pd) # Bayesian p-value for perceptibility model
      # # bayesp.pa<-step(fit.new.pa-fit.pa) # Bayesian p-value for availability model
      # meanN<-mean(N[]) # mean site-level abundance
      # totN<-sum(N[])  # population size of total area surveyed
      # # meansig<-mean(sigma[]) # mean scale parameter across sites
      # dens<-meanN/(maxd*maxd*3.14159/10000) # density of birds per ha
      
      }
      
      ", fill=TRUE)
  sink()