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
library(rjags)
library(jagsUI)

# Load Data
df_counts <- read_csv("Data/point_counts_2015.csv")
str(df_counts)
summary(df_counts)

df_counts[6:ncol(df_counts)] %>%
  dplyr::summarise_each(funs(sum)) %>%
  t(.)

df_counts[6:ncol(df_counts)] %>%
  dplyr::summarise_each(funs(sd)) %>%
  t(.)

# GRCA one of most abundant

# Get detection covariates
df_detect <- read.csv("Data/covs_detect.csv", stringsAsFactors = FALSE, header = TRUE)

df_detect$Date <- mdy(df_detect$Date)
df_detect$Time <- strptime(df_detect$Time, "%H:%M")
df_detect$timeofday <-as.numeric(df_detect$Time)
df_detect$day <- as.integer(strftime(df_detect$Date, format = "%j"))

str(df_detect)

# Standardize continuous covariates
# beta.a0 + beta.a1*day[k] + beta.a2*day[k]*day[k] + beta.a3*time[k] + beta.a4*time[k]*time[k]
std_covs <- function(x) {
  std <- (x - mean(x)) / sd(x)
  return(std)
}
# std_covs(df_detect$day) # test
df_detect <- df_detect %>%
  dplyr::select(-Time) %>%
  dplyr::filter(Year == "2015") %>%
  dplyr::mutate(day_std = std_covs(day),
                time_std = std_covs(timeofday))

summary(df_detect)

# Get Abundance covariates
df_abund <- read.csv("Data/covs_abund.csv", stringsAsFactors = FALSE, header = TRUE)

str(df_abund)

# Standardize abundance covariates
#beta0 + beta1*GMU[k] + beta2*VegHgt[k] + beta3*VegHgt_SD[k] 
df_abund2 <- df_abund[1:8]

df_abund3 <- df_abund[9:ncol(df_abund)] %>%
  dplyr::mutate()

df_abund_std <- df_abund %>%
  dplyr::mutate_each(., funs(std_covs), VegHgt_Avg, VegHgt_Stdv)

# check for excessive correlation (>0.7 once standardized)
### Scatterplot Matrix ###
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
}


panel.diff <- function(x,y, digits=2, prefix="", cex.cor){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Diff <- max(abs(x - y))
  txt <- format(c(Diff, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  #if(missing(cex.cor)) cex <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1)
}

Pairs <- df_abund_std[ , c("GMU", "VegHgt_Avg", "VegHgt_Stdv", "LD_Avg", "Avg_Canopy", "GCVR_Forb", "Shrub_stm_total", "Tree_stm_total")]
pairs(Pairs, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)

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

# example may assume equal time intervals **************************

# GRCA one of most abundant
y <- df_counts %>%
  dplyr::filter(Visit == 1) %>%
  dplyr::select(GRCA) %>%
  .[["GRCA"]]

surveyid <- df_counts %>%
  dplyr::filter(Visit == 1) %>%
  .[["ID"]]

dclass <- df_counts %>%
  dplyr::filter(Visit == 1) %>%
  .[["Dist_bin"]]

nsites <- length(unique(surveyid))

tinterval <- df_counts %>%
  dplyr::filter(Visit == 1) %>%
  .[["Time_bin"]]
  
jags_data<-list(y=y,
               surveyid=as.numeric(surveyid),
               dclass=as.numeric(dclass),
               nsurveys=nsites,
               nobs=sum(y, na.rm = TRUE),
               delta=50, # c(50, 50, 100),
               nbreaks=3,
               mdpts=c(25, 75, 125),
               maxd=150,
               J=max(unique(tinterval)),
               tinterval=as.numeric(tinterval)
               )

# create initial values for N and navail that are very close to actual values or model will not run!
Nst <- jags_data$y + 1

# Inits function
inits <- function() {
  list(N=Nst,
       navail=Nst,
       sigma.0=runif(1,120,150), 
       beta.a1=runif(1,-1,1),
       beta1=runif(1,-1,1),
       beta2=runif(1,-1,1),
       mu.tran=runif(1,-1,1),
       sd.tran=runif(1,0,2),
       beta3=runif(1,-1,1),
       beta4=runif(1,-1,1),
       beta5=runif(1,-1,1),
       beta.p1=runif(1,-1,1),
       beta.a0=runif(1,-1,1)
       )
  } #

# parameters to estimate
# careful printing pavail, pdet and N - will have nsites (e.g., in this example 100) values
# params<-c("meansig","meanpdet","meanpavail","beta.a0","beta.a1","sigma.0","beta.p1","beta1","beta2","beta3","beta4","beta5","meanN","mu.tran","sd.tran","totN","bayesp.pa","bayesp.pd","beta0","beta.tran","N")
params<-c("beta.a0", "meanpavail", "meanpdet", "meanN", "totN", "dens", "N")

# MCMC settings
# pavail can be subject to poor mixing in field data - keep thin high, burn-in long, and conduct sufficient number of iterations
nc<-3
ni<-6000
nb<-3000
nt<-3

## ONLY WORKS IN JAGS
# A bug fix for JAGS - model may produce error without this fix
set.factory("bugs::Conjugate", FALSE, type="sampler")

sim_fit<-jags(data=jags_data,parameters.to.save=params, model.file="jags_unequal_time.txt",
              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
              parallel = TRUE,
              n.cores = 3) #  inits=inits, 

names(sim_fit)
str(sim_fit)
summary(sim_fit)
jagsUI::traceplot(sim_fit, parameters = c("beta.a0", "N[1]", "dens"))



