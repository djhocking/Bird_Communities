###############################################################################
# Estimate Abundance of birds from point counts with dist and removal time
# Daniel J. Hocking
# Collaborators: Meta Griffin, Frank Ammer
# October 2016
###############################################################################

# set testing
testing <- TRUE

# Load Packages
library(dplyr)
library(lubridate)
library(tidyr)
library(readr)
library(rjags)
library(jagsUI)

# Load Data
df_2015 <- read_csv("Data/point_counts_2015.csv")
str(df_2015)
summary(df_2015)
length(unique(df_2015$Point))

df_2016 <- read_csv("Data/point_counts_2016.csv")
str(df_2016)
length(unique(df_2016$Point))

# a couple sites were not surveyed each year
unique(df_2015$Point) %in% unique(df_2016$Point) # all 2015 sites visited again in 2016
unique(df_2016$Point) %in% unique(df_2015$Point) # 2 new 2016 sites
(new_points <- unique(df_2016$Point)[(unique(df_2016$Point) %in% unique(df_2015$Point)) == FALSE])
unique(df_2015$Point)

# combine 2 years of point counts with 2 visits each year
df_counts <- bind_rows(df_2015, df_2016, .id = "id")

# replace NA for all species not found in both years with 0
df_counts[is.na(df_counts)] <- 0

# return NA for all species at 2 sites not surveyed in 2015
# df_counts2 <- df_counts %>%
#   dplyr::select(-id, -ID, -Visit, -Time_bin, Dist_bin) %>%
#   dplyr::mutate_each(., funs(ifelse()))

#-----skip for now and just exclude these 2 sites-------
length(unique(df_counts$Point))
df_counts <- df_counts %>%
  dplyr::filter(!(Point %in% new_points))
length(unique(df_counts$Point))

dim(df_2015)
dim(df_2016)
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

length(unique(df_detect$Point))
df_detect <- df_detect %>%
  dplyr::filter(!(Point %in% new_points))
length(unique(df_detect$Point))

str(df_detect)
summary(df_detect) # can't have any NA

# Standardize continuous covariates
# beta.a0 + beta.a1*day[k] + beta.a2*day[k]*day[k] + beta.a3*time[k] + beta.a4*time[k]*time[k]
std_covs <- function(x) {
  std <- (x - mean(x)) / sd(x)
  return(std)
}
# std_covs(df_detect$day) # test
df_detect <- df_detect %>%
  dplyr::select(-Time) %>%
  # dplyr::filter(Year == "2015") %>%
  dplyr::mutate(day_std = std_covs(day),
                time_std = std_covs(timeofday))

summary(df_detect)

# Get Abundance covariates
df_abund <- read.csv("Data/covs_abund.csv", stringsAsFactors = FALSE, header = TRUE)

length(unique(df_abund$Point))
df_abund <- df_abund %>%
  dplyr::filter(!(Point %in% new_points))
length(unique(df_abund$Point))

str(df_abund)
summary(df_abund) # can't have any NA if used in the model

# Standardize abundance covariates
#beta0 + beta1*GMU[k] + beta2*VegHgt[k] + beta3*VegHgt_SD[k] 
df_abund2 <- df_abund[1:8]

df_abund3 <- df_abund[9:ncol(df_abund)] %>%
  dplyr::mutate()

df_abund_std <- df_abund %>%
  dplyr::mutate_each(., funs(std_covs), VegHgt_Avg, VegHgt_Stdv, LD_Avg, LD_Stdv, Avg_Canopy, GCVR_Forb, GCVR_ShrubTree, Shrub_stm_total, Tree_stm_total)

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

Pairs <- df_abund_std[ , c("GMU", "VegHgt_Avg", "VegHgt_Stdv", "LD_Avg", "LD_Stdv", "Avg_Canopy", "GCVR_Forb", "Shrub_stm_total", "Tree_stm_total")]
pairs(Pairs, upper.panel=panel.smooth, lower.panel=panel.cor, diag.panel=panel.hist)

# Use (GMU - maybe), VegHgt_Avg, Shrub_stm_total, Tree_stm_total

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
  # dplyr::filter(Visit == 1) %>%
  dplyr::select(YEWA) %>%
  .[["YEWA"]]

# make variable for point-year
df_counts$Pyear <- paste0(df_counts$ID, "_", df_counts$Year)
df_detect$Pyear <- paste0(df_detect$ID, "_", df_detect$Year)
df_abund_std$Pyear <- paste0(df_abund_std$ID, "_", df_abund_std$Year)

length(unique(df_counts$Pyear))
length(unique(df_detect$Pyear))
length(unique(df_abund_std$Pyear))

missing_abund_data <- unique(df_abund_std$Pyear)[!(unique(df_counts$Pyear) %in% unique(df_abund_std$Pyear))]

# missing abundance data for one point-year combination - use previous year for that one location or exclude
df_counts <- df_counts %>%
  dplyr::filter(!(Pyear %in% missing_abund_data))
df_detect <- df_detect %>%
  dplyr::filter(!(Pyear %in% missing_abund_data))

# make variable for point-visit
df_counts$Survey <- paste0(df_counts$ID, "_", df_counts$Year, "_", df_counts$Visit)
df_detect$Survey <- paste0(df_detect$ID, "_", df_detect$Year, "_", df_detect$Visit)
length(unique(df_detect$Survey))

surveyid <- df_counts %>%
  # dplyr::filter(Visit == 1) %>%
  .[["Survey"]]

detect_surveys <- unique(df_detect$Survey)

length(unique(surveyid))
length(detect_surveys)
(new_surveys <- unique(detect_surveys)[unique(detect_surveys) %in% (unique(surveyid)) == FALSE])

df_detect <- df_detect %>%
  dplyr::filter(!(Survey %in% new_surveys))

dclass <- df_counts %>%
  # dplyr::filter(Visit == 1) %>%
  .[["Dist_bin"]]

nsurveys <- length(unique(surveyid))


tinterval <- df_counts %>%
  # dplyr::filter(Visit == 1) %>%
  .[["Time_bin"]]

# df_abund_std <- df_abund_std %>%
#   dplyr::filter(Year == 2015)

# point <- as.integer(df_counts[which(df_counts$Visit == 1), ]$ID)
point <- as.integer(df_counts$ID)
n_points <- length(unique(point))
  
jags_data<-list(y=y,
               surveyid=as.numeric(as.factor(surveyid)),
               dclass=as.numeric(dclass),
               nsurveys=nsurveys,
               pyear = as.numeric(as.factor(df_detect$Pyear)),
               nobs=sum(y, na.rm = TRUE),
               n_points = n_points,
               point = point,
               delta=50, # c(50, 50, 100),
               nbreaks=3,
               mdpts=c(25, 75, 125),
               maxd=150,
               J=max(unique(tinterval)),
               tinterval=as.numeric(tinterval),
               day = df_detect$day_std,
               time = df_detect$time_std,
               VegHgt = df_abund_std$VegHgt_Avg,
               Shrub_stm_total = df_abund_std$Shrub_stm_total,
               Tree_stm_total = df_abund_std$Tree_stm_total,
               Year = df_detect$Year - 2015
               )

# create initial values for N and navail that are very close to actual values or model will not run!
Nst <- jags_data$y + 2
navail <- jags_data$y + 1

# Inits function
inits <- function() {
  list(N=Nst,
       # navail=Nst,
       # sigma.0=runif(1,120,150), 
       # beta.a1=runif(1,-1,1),
       # beta1=runif(1,-1,1),
       # beta2=runif(1,-1,1),
       # mu.tran=runif(1,-1,1),
       # sd.tran=runif(1,0,2),
       # beta3=runif(1,-1,1),
       # beta4=runif(1,-1,1),
       # beta5=runif(1,-1,1),
       # beta.p1=runif(1,-1,1),
        beta.a0=runif(1,0,9)
       )
  } #

# parameters to estimate
# careful printing pavail, pdet and N - will have nsites (e.g., in this example 100) values
# params<-c("meansig","meanpdet","meanpavail","beta.a0","beta.a1","sigma.0","beta.p1","beta1","beta2","beta3","beta4","beta5","meanN","mu.tran","sd.tran","totN","bayesp.pa","bayesp.pd","beta0","beta.tran","N")
params<-c("beta.a0", "beta.a1", "beta.a2", "beta.a3", "beta.a4", "beta0", "ln.beta0", "beta1", "beta2", "beta3", "beta4", "meanpavail", "meanpdet", "meanN", "totN", "dens", "N", "pavail", "pdet", "sigma.0", "beta.p1", "sigma") # "sigma.eps.n", 

# MCMC settings
# pavail can be subject to poor mixing in field data - keep thin high, burn-in long, and conduct sufficient number of iterations
nc<-3
ni<-18000
nb<-9000
nt<-3

if(testing == TRUE) {
  ni<-60
  nb<-30
}

## ONLY WORKS IN JAGS
# A bug fix for JAGS - model may produce error without this fix
set.factory("bugs::Conjugate", FALSE, type="sampler")

sim_fit<-jags(data=jags_data,parameters.to.save=params, model.file="jags_full.txt",
              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
              parallel = TRUE,
              n.cores = 3) #  inits=inits, 

# names(sim_fit)
# str(sim_fit)
# summary(sim_fit)
jagsUI::traceplot(sim_fit, parameters = c("beta.a0", "beta.a1", "beta.a2", "beta.a3", "beta.a4", "beta0", "beta1", "beta2", "beta3", "beta4", "sigma.eps.n", "N[1]", "dens"))

jagsUI::traceplot(sim_fit, parameters = c("beta.a0", "beta.a1", "beta.a2", "beta.a3", "beta.a4", "beta0", "N[1]", "sigma.0", "beta.p1"))

# minimum number of birds (total observed)
N_min <- df_counts %>%
  # dplyr::filter(Visit == 1) %>%
  dplyr::group_by(id, Point, ID, Year, Visit, Pyear, Survey) %>%
  dplyr::summarise_each(., funs(sum)) %>%
  dplyr::select(Point, AMGO) %>%
  .[["AMGO"]]

# sim_fit$summary[ , "50%"]
N_est <- sim_fit$summary %>%
  as.data.frame(.) %>%
  dplyr::mutate(parameter = rownames(sim_fit$summary)) %>%
  dplyr::filter(., grepl("N[", parameter, fixed = TRUE))

N_est <- N_est[ , c("parameter", "50%", "2.5%", "97.5%")]

cbind(N_est, N_min)
