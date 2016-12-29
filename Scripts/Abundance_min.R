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

df_counts <- df_counts %>%
  dplyr::filter(Visit == 1)

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

(total_counts <- df_counts[6:ncol(df_counts)] %>%
    dplyr::summarise_each(funs(sum)) %>%
    t(.))

variability <- df_counts[6:ncol(df_counts)] %>%
  dplyr::summarise_each(funs(sd)) %>%
  t(.)

# Get detection covariates
df_detect <- read.csv("Data/covs_detect.csv", stringsAsFactors = FALSE, header = TRUE) 

df_detect <- df_detect %>%
  dplyr::filter(Visit == 1)

df_detect$Date <- mdy(df_detect$Date)
df_detect$Time <- strptime(df_detect$Time, "%H:%M")
df_detect$timeofday <-as.numeric(df_detect$Time)
df_detect$day <- as.integer(strftime(df_detect$Date, format = "%j"))

length(unique(df_detect$Point))
df_detect <- df_detect %>%
  dplyr::select(-Time) %>%
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
  # dplyr::filter(Year == "2015") %>%
  dplyr::mutate(day_std = std_covs(day),
                time_std = std_covs(timeofday))

summary(df_detect)

# Get Abundance covariates
df_abund <- read.csv("Data/covs_abund.csv", stringsAsFactors = FALSE, header = TRUE) 

df_abund <- df_abund %>%
  dplyr::filter(Year == 2015)

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

# Use (GMU - maybe), VegHgt_Avg, Shrub_stm_total, Tree_stm_total

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

detect_surveys <- unique(df_detect$Survey)
surveyid <- unique(df_counts$Survey)

length(unique(surveyid))
length(detect_surveys)
(new_surveys <- unique(detect_surveys)[unique(detect_surveys) %in% (unique(surveyid)) == FALSE])

df_detect <- df_detect %>%
  dplyr::filter(!(Survey %in% new_surveys))


#-------------------- Create observation level data ---------------

## Now create observation level dataset that can be sampled
# NEED TO RESTRICT COUNT TO INDIVIDUALS IN POP THEN LATER ADD ZEROS

df_z_spp <- df_counts %>%
  dplyr::group_by(Point, ID, Year, Visit, Pyear, Survey) %>%
  dplyr::select(-Time_bin, -Dist_bin, -id) %>%
  dplyr::summarise_each(., funs(sum))

df_bins <- df_counts %>%
  dplyr::select(Point, ID, Year, Visit, Pyear, Survey, Time_bin, Dist_bin) %>%
  dplyr::distinct()

df_z_spp$SurveyID <- 1:nrow(df_z_spp)

# make surveyid df
df_surveyid <- df_z_spp %>%
  dplyr::select(Point, ID, Year, Visit, Pyear, Survey, SurveyID) 

# loop
i <- 2 #sp on col 8-81
spp <- names(df_counts[ , 8:81])
sp <- spp[i]
vars <- c("Point", "ID", "Year", "Visit", "Pyear", "Survey", "SurveyID", sp)
vars_full <- c("Point", "ID", "Year", "Visit", "Pyear", "Survey", "SurveyID", "Time_bin", "Dist_bin", sp)

df_obs <- df_counts %>%
  dplyr::left_join(df_surveyid) %>%
  dplyr::select(one_of(vars_full)) %>%
  dplyr::rename_(count = sp) %>%
  dplyr::filter(count > 0)

df_z <- df_z_spp %>%
  # dplyr::filter(Visit == 1) %>%
  dplyr::select(one_of(vars)) %>%
  dplyr::rename_(count = sp)
# 
# df_zz <- subset(df_z, count > 0) # only individuals in population-sites with zero counts are removed temporarily
# str(df_zz)

# df_obsonly <- data.frame()
# for(i in 1:length(df_zz$count)){ # basically make 1s for every individual at a site, then link to site id
#   tmp=df_zz[i,]
#   ind=rep(1,tmp$count)
#   SurveyID=rep(tmp$SurveyID, tmp$count)
#   df_obsonly=rbind(df_obsonly, data.frame(ind, SurveyID))
# }

df_ind <- data.frame()
for(i in 1:length(df_obs$count)){ # basically make 1s for every individual at a site, then link to site id
  tmp=df_obs[i,]
  ind=rep(1,tmp$count)
  SurveyID=rep(tmp$SurveyID, tmp$count)
  Time_bin = rep(tmp$Time_bin, tmp$count)
  Dist_bin = rep(tmp$Dist_bin, tmp$count)
  df_ind=rbind(df_ind, data.frame(ind, SurveyID, Time_bin, Dist_bin))
}

df_ind <- left_join(df_ind, df_surveyid)

# Z_vec <- df_obsonly$ind  # create vector of only individuals in population (no zero counts)

# Create distances for each observation - translate them to pdet (perceptibility)
maxd <- 150
# n_cut <- sum(Z_vec) # total number of observations
# delta <- maxd / (n_cut + 1)
# dclass <- seq(delta, maxd, , n_cut)

# data <- NULL # create null dataset to insert values

## Create data frame with all our info
# data <- data.frame(y,u,v,d=as.numeric(d3),tint,sigma=obsonly$sigma,siteid=obsonly$siteid)

# Now add in the zero counts from sites where no birds observed - see creation of 'z' above for reference
# fulldata<-merge(df_z, data, all=TRUE) # observations plus sites with zero counts that are NA
# x<-subset(data,y==1) # break it down to observations again, but with additional info from other columns

# df_z <- df_z_spp %>%
#   # dplyr::filter(Visit == 1) %>%
#   dplyr::select(one_of(vars)) %>%
#   dplyr::rename_(count = sp)

############################ FORMAT DATA FOR JAGS/BUGS ENTRY ##########################################################

# Subset to individuals with greater than zero count again for JAGS data entry
# fulldata[is.na(fulldata)]=0 # change NAs to zeros; count doesn't make sense in this data frame - need to summarize by site into new.count
# # get observed counts per site
# install.packages("plyr")
# library(plyr)
# new.count<-ddply(fulldata,~siteid,summarize,y=sum(y),count=median(count),mtint=mean(tint),siteid=max(siteid))
# str(new.count) # make sure it has R rows and sum of count = pop
# sum(new.count$count) # must equal A - pop size
# sum(new.count$y) # total number of individuals observed
# 
# y<-new.count$y # this is now the number of birds seen per site/point
# y # make sure it looks appropriate

y <- df_z$count
n_breaks <- 3

tinterval <- df_ind$Time_bin # time interval for each observation ONLY (not sites with zero counts obviously)
nobs <- length(tinterval) # total count of individuals observed
if(nobs != sum(df_z$count)) stop # ensure that correct number of obs
nsurveys <- length(y) # number of surveys (site visits)
surveyid <- df_ind$SurveyID # survey point/site ID only for individuals observed
delta <- maxd / n_breaks # bin size or width
mdpts <- seq(delta / 2, maxd , delta) # generates midpoint distance of bins up to max distance
dclass <- df_ind$Dist_bin # distance class for each observation
point <- as.integer(df_z$ID)
n_points <- length(unique(point))

#------------

jags_data<-list(y=y,
                surveyid=as.numeric(surveyid),
                dclass=as.numeric(dclass),
                nsurveys=as.integer(nsurveys),
                pyear = as.numeric(as.factor(df_z$Pyear)),
                nobs=nobs,
                n_points = as.integer(n_points),
                point = as.integer(point),
                delta=delta, # c(50, 50, 100),
                nbreaks=n_breaks,
                mdpts=mdpts,
                maxd=maxd,
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
Nst <- jags_data$y + 1
navail <- jags_data$y + 1

# Inits function
inits <- function() {
  list(N=Nst,
        navail=Nst#,
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
       #beta.a0=runif(1,0,9)
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

sim_fit<-jagsUI::jags(data=jags_data,parameters.to.save=params, model.file="jags_min.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
                      parallel = TRUE,
                      n.cores = 3) #  inits=inits, 

# names(sim_fit)
# str(sim_fit)
# summary(sim_fit)

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

df_year <- df_counts %>%
  select(Point, ID, Year, Visit, Survey) %>%
  distinct()
df_year

df_survey <- dplyr::left_join(data.frame(Survey = unique(surveyid), SurveyID = unique(as.numeric(as.factor((surveyid)))), stringsAsFactors = FALSE), df_year)
data_frame(df_survey)
cbind(df_survey, N_est, N_min)

# when not testing - change these to save rather than print
if(testing) {
  jagsUI::traceplot(sim_fit, parameters = c("beta.a0", "beta.a1", "beta.a2", "beta.a3", "beta.a4", "beta0", "beta1", "beta2", "beta3", "beta4", "sigma.eps.n", "N[1]", "dens"))
  
  jagsUI::traceplot(sim_fit, parameters = c("beta.a0", "beta.a1", "beta.a2", "beta.a3", "beta.a4", "beta0", "N[1]", "sigma.0", "beta.p1"))
}




