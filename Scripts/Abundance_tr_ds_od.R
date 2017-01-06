# install.packages("AHMbook", repos = NULL, type="source")
# install.packages("plotrix")
# library(AHMbook)

# set testing
testing <- FALSE

# Load Packages
library(dplyr)
library(lubridate)
library(tidyr)
library(readr)
library(rjags)
library(jagsUI)
set.factory("bugs::Congugate", FALSE, type = "sampler")

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

(total_counts <- df_counts[6:ncol(df_counts)] %>%
    dplyr::summarise_each(funs(sum)) %>%
    t(.))

variability <- df_counts[6:ncol(df_counts)] %>%
  dplyr::summarise_each(funs(sd)) %>%
  t(.)

# plot counts over time and distance to see effects
library(ggplot2)
df_time <- df_counts %>%
  dplyr::group_by(Time_bin, Visit, Year) %>%
  select(Year, Visit, Time_bin, AMGO) %>%
  dplyr::mutate(count_min = ifelse(Time_bin == 1, AMGO/3, ifelse(Time_bin == 2, AMGO/2, AMGO/5))) %>%
  summarise(count_min = sum(count_min))

ggplot(data = df_time, aes(Time_bin, count_min)) + geom_line() + geom_point() + facet_wrap(~Year + Visit)

area1 <- pi * 50 ^ 2
area2 <- (pi * 100 ^ 2) - area1
area3 <- (pi * 150 ^ 2) - area1 - area2

df_dist <- df_counts %>%
  dplyr::group_by(Dist_bin, Visit, Year) %>%
  select(Year, Visit, Dist_bin, AMGO) %>%
  dplyr::mutate(count_min = ifelse(Dist_bin == 1, AMGO/area1, ifelse(Dist_bin == 2, AMGO/area2, AMGO/area3))) %>%
  summarise(count_min = sum(count_min))

ggplot(data = df_dist, aes(Dist_bin, count_min)) + geom_line() + geom_point() + facet_wrap(~Year + Visit)

# Get detection covariates
df_detect <- read.csv("Data/covs_detect.csv", stringsAsFactors = FALSE, header = TRUE)

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

str(df_abund_std)

# check for excessive correlation (>0.7 once standardized)
# Use VegHgt_Avg, Shrub_stm_total, Tree_stm_total

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
length(unique(df_detect$Survey))

# summarise counts to survey (rather than survey-time-dist)
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

# expand covariates so one per survey
df_abund_std <- df_surveyid %>%
  dplyr::left_join(df_abund_std)

df_abund_std[is.na(df_abund_std)] <- 0

df_detect <- df_surveyid %>%
  dplyr::left_join(df_detect) %>%
  dplyr::select(-Date)

df_detect[is.na(df_detect)] <- 0

df_bins <- df_counts %>%
  dplyr::select(Point, ID, Year, Visit, Pyear, Survey, Time_bin, Dist_bin) %>%
  dplyr::distinct()

#-------------------- Create observation level data ---------------

## Now create observation level dataset that can be sampled
# NEED TO RESTRICT COUNT TO INDIVIDUALS IN POP THEN LATER ADD ZEROS
converge <- data.frame(species = NA, converge = NA)

# loop

# i <- 2 #sp on col 8-81
spp <- names(df_counts[ , 8:81])

for(i in 1:length(unique(spp))) {

sp <- spp[i]
vars <- c("Point", "ID", "Year", "Visit", "Pyear", "Survey", "SurveyID", sp)
vars_full <- c("Point", "ID", "Year", "Visit", "Pyear", "Survey", "SurveyID", "Time_bin", "Dist_bin", sp)

df_obs <- df_counts %>%
  dplyr::left_join(df_surveyid) %>%
  dplyr::select(one_of(vars_full)) %>%
  dplyr::rename_(count = sp) %>%
  dplyr::filter(count > 0)

df_z <- df_z_spp %>%
  dplyr::select(one_of(vars)) %>%
  dplyr::rename_(count = sp) %>%
  dplyr::mutate(visit2 = ifelse(Visit == 2 & Year == 2015, 1, 0),
                visit3 = ifelse(Visit == 1 & Year == 2016, 1, 0),
                visit4 = ifelse(Visit == 2 & Year == 2016, 1, 0))

df_ind <- data.frame()
for(j in 1:length(df_obs$count)){ # basically make 1s for every individual at a site, then link to site id
  tmp=df_obs[j,]
  ind=rep(1,tmp$count)
  SurveyID=rep(tmp$SurveyID, tmp$count)
  Time_bin = rep(tmp$Time_bin, tmp$count)
  Dist_bin = rep(tmp$Dist_bin, tmp$count)
  df_ind=rbind(df_ind, data.frame(ind, SurveyID, Time_bin, Dist_bin))
}

df_ind <- left_join(df_ind, df_surveyid)


K <- max(df_counts$Time_bin) # number of removal time periods
nsurveys <- length(unique(df_surveyid$SurveyID))


# Create the observed encounter frequencies per site (include the zeros! )
n <- df_z$count               # The full site vector
names(n) <- df_surveyid$SurveyID
survey <- df_ind$SurveyID
nobs <- length(survey) 
npoints <- length(unique(df_z$ID))
df_points <- data.frame(ID = unique(df_z$ID), stringsAsFactors = FALSE)
df_points$point <- 1:npoints
df_z <- dplyr::left_join(df_z, df_points)
point <- df_z$point

y <- df_z$count
n_breaks <- 3

# Create the distance class data
maxd <- 150
maxd <- maxd / 100
nD <- 3             # Number of distance classes 
delta <- maxd / nD        # bin size or width
mdpts <- seq(delta / 2, maxd , delta) # midpoint distance of bins up to max distance
dclass <- df_ind$Dist_bin # distance class for each observation
tint <- df_ind$Time_bin

# Bundle data and summarize Use VegHgt_Avg, Shrub_stm_total, Tree_stm_total
str(jags.data <- list(n=n, 
                      site=survey, 
                      dclass=as.numeric(dclass),
                      nsites=nsurveys, 
                      nobs=nobs, 
                      npoints = npoints,
                      point = point,
                      delta=delta,
                      maxd = maxd,
                      nD=nD,
                      mdpts=mdpts,
                      B=3, 
                      K=K, 
                      tint=tint, 
                      year = df_abund_std$Year - min(df_abund_std$Year),
                      time = as.numeric(df_detect$time_std),
                      day = as.numeric(df_detect$day_std),
                      visit2 = as.numeric(df_z$visit2),
                      visit3 = as.numeric(df_z$visit3),
                      visit4 = as.numeric(df_z$visit4),
                      siteVisit = as.numeric(df_z$SurveyID),
                      veg=as.numeric(df_abund_std$VegHgt_Avg),
                      shrub=as.numeric(df_abund_std$Shrub_stm_total),
                      tree=as.numeric(df_abund_std$Tree_stm_total)))

# Create initial values (including for M and N) and list parameters to save
Nst <- n + 1
Mst <- n + 1 
inits <- function(){
  list(M=Mst, 
       alpha0=runif(1, 1, 5),
       beta0=runif(1,-1,1),
       beta.a1=runif(1,-1,1),
       beta1=runif(1,-1,1),
       alpha1=runif(1,-1,1),
       beta.a0=runif(1,-1,1),
       N=Nst
  )
}
params <- c("beta.a0", "beta.a1", "beta.a2", "beta.a3", "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "beta4", "sigma.lam", "sigma.dist", "sigma.time", "N", "M", "PDETmean", "PAVAILmean", "Mtot", "Ntot", "meansig", "dens", "bayesp.pd", "bayesp.pa", "pavail", "pdet") #, "sigma") # "beta.a4", "alpha2", "alpha1", "alpha2", "alpha3", "alpha4",

# MCMC settings
ni <- 150000
nb <- 50000  # 10000
nt <- 10     # 18
nc <- 3

if(testing) {
  ni = 1000
  nb = 500
  nt = 1
  nc = 3
}

# Run JAGS in parallel (ART 7.3 min), check convergence and summarize posteriors
start_t <- proc.time()
sim_fit <- try(jags(data=jags.data, inits=inits, parameters=params, 
              model.file ="Scripts/Jags/tr-ds-od.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, 
              parallel = TRUE))
proc.time() - start_t # 9 hours

if(testing) {
  dev.off()
jagsUI::traceplot(sim_fit, parameters = c("beta.a0", "beta.a1", "beta.a2", "beta.a3", "alpha0", "beta0", "beta1", "beta2", "beta3", "beta4", "sigma.lam", "sigma.dist", "sigma.time", "M[1]", "N[1]", "Mtot", "Ntot", "meansig", "dens", "bayesp.pd", "bayesp.pa", "PDETmean", "PAVAILmean")) #, "sigma[1]")) # "beta.a4", "alpha2","alpha1", "alpha2", "alpha3", "alpha4", 
  print(sim_fit, 3)
} else {
  if(class(sim_fit) != "try-error") {
  saveRDS(sim_fit, file = paste0("Output/MCMC/", sp, ".Rds"))
  library(ggmcmc)
  ggmcmc(ggs(sim_fit$samples, family = "^alpha|^beta|^sigma|^Mtot|^dens"), file = paste0("Output/traceplots_od_", sp, ".pdf"), plot=c("traceplot"))

# ci.median <- ci(ggs(sim_fit$samples, family="^N")) %>%
#   select(Parameter, median)
# 
# L.radon <- data.frame(
#   Parameter=c(
#     paste("N[", df_surveyid$SurveyID, "]", sep="")),
#   Label=rep(df_surveyid$Point, 1),
#   Year=rep(as.character(df_surveyid$Year, 1)),
#   Visit=rep(as.character(df_surveyid$Visit, 1)),
#   Coefficient=gl(1, length(df_surveyid$SurveyID), 
#                  labels=c("Abundance")))
# 
# ggs_caterpillar(ggs(sim_fit$samples, par_labels=L.radon, family="^N")) +
#   facet_wrap(~ Year + Visit, scales="free") #+ aes(color=Uranium)

# Record convergence
# converge[i, "converge"] <- max(unlist(sim_fit$Rhat)) <= 1.1
# converge[i, "species"] <- sp

# minimum number of birds (total observed)
N_min <- df_counts %>%
  # dplyr::filter(Visit == 1) %>%
  dplyr::group_by(id, Point, ID, Year, Visit, Pyear, Survey) %>%
  dplyr::summarise_each(., funs(sum)) %>%
  dplyr::select_("Point", sp) %>%
  .[[sp]]

# sim_fit$summary[ , "50%"]
bayesp <- sim_fit$summary %>%
  as.data.frame(.) %>%
  dplyr::mutate(parameter = rownames(sim_fit$summary)) %>%
  dplyr::filter(., grepl("bayes", parameter, fixed = TRUE))
bayesp

sim_fit$summary %>%
  as.data.frame(.) %>%
  dplyr::mutate(parameter = rownames(sim_fit$summary)) %>%
  dplyr::filter(., grepl("sigma", parameter, fixed = TRUE))

N_est <- sim_fit$summary %>%
  as.data.frame(.) %>%
  dplyr::mutate(parameter = rownames(sim_fit$summary)) %>%
  dplyr::filter(., grepl("M[", parameter, fixed = TRUE))

N_est <- N_est[ , c("parameter", "50%", "2.5%", "97.5%")]

df_year <- df_counts %>%
  select(Point, ID, Year, Visit, Survey) %>%
  distinct()
df_year

df_survey <- dplyr::left_join(df_surveyid, df_year)
df_survey <- data.frame(df_survey, stringsAsFactors = FALSE)
N_est <- cbind(df_survey, N_est, N_min) %>%
  dplyr::select(-Pyear, -Survey, -parameter) %>%
  dplyr::mutate(Species = sp)
names(N_est) <- c("Point", "ID", "Year", "Visit", "SurveyID", "Median", "LCRI", "UCRI", "Nmin", "Species")
if(i == 1) {
  N_table <- N_est
} else {
  N_table <- bind_rows(N_table, N_est)
}
write.csv(N_table, file = paste0("Output/N_table_od.csv"), row.names = FALSE)

library(ggplot2)
ggplot(N_est, aes(N_min, Median)) + geom_point() + theme_bw() + ggtitle(sp)
ggsave(file = paste0("Output/Figures/Obs_Pred_od_", sp, ".pdf"))

if(i == 1) {
  summary_table <- data.frame(matrix(NA, nrow = 1, ncol = 9))
  names(summary_table) <- c("Species", "PDETmean", "PAVAILmean", "meansig", "dens", "bayesp.pd", "bayesp.pa", "Rhat_max", "Converged")
  
  summary_table[i, "Species"] <- sp
  summary_table[i, "PDETmean"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("PDETmean", parameter, fixed = TRUE)) %>% dplyr::select(5))
  summary_table[i, "PAVAILmean"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("PAVAILmean", parameter, fixed = TRUE)) %>% dplyr::select(5))
  summary_table[i, "meansig"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("meansig", parameter, fixed = TRUE)) %>% dplyr::select(5))
  summary_table[i, "dens"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("dens", parameter, fixed = TRUE)) %>% dplyr::select(5))
  summary_table[i, "bayesp.pd"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("bayesp.pd", parameter, fixed = TRUE)) %>% dplyr::select(mean))
  summary_table[i, "bayesp.pa"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("bayesp.pa", parameter, fixed = TRUE)) %>% dplyr::select(mean))
  summary_table[i, "Rhat_max"] <- max(unlist(sim_fit$Rhat), na.rm = TRUE)
  summary_table[i, "Converged"] <- ifelse(max(unlist(sim_fit$Rhat), na.rm = TRUE) <= 1.1, 1, 0)
} else {
  summary_table[i, "Species"] <- sp
  summary_table[i, "PDETmean"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("PDETmean", parameter, fixed = TRUE)) %>% dplyr::select(5))
  summary_table[i, "PAVAILmean"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("PAVAILmean", parameter, fixed = TRUE)) %>% dplyr::select(5))
  summary_table[i, "meansig"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("meansig", parameter, fixed = TRUE)) %>% dplyr::select(5))
  summary_table[i, "dens"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("dens", parameter, fixed = TRUE)) %>% dplyr::select(5))
  summary_table[i, "bayesp.pd"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("bayesp.pd", parameter, fixed = TRUE)) %>% dplyr::select(mean))
  summary_table[i, "bayesp.pa"] <- as.numeric(as.data.frame(sim_fit$summary) %>% dplyr::mutate(parameter = rownames(sim_fit$summary)) %>% dplyr::filter(., grepl("bayesp.pa", parameter, fixed = TRUE)) %>% dplyr::select(mean))
  summary_table[i, "Rhat_max"] <- max(unlist(sim_fit$Rhat), na.rm = TRUE)
  summary_table[i, "Converged"] <- ifelse(max(unlist(sim_fit$Rhat), na.rm = TRUE) <= 1.1, 1, 0)
}
  } else {
    summary_table[i, "Species"] <- sp
    summary_table[i, "PDETmean"] <- NA
    summary_table[i, "PAVAILmean"] <- NA
    summary_table[i, "meansig"] <- NA
    summary_table[i, "dens"] <- NA
    summary_table[i, "bayesp.pd"] <- NA
    summary_table[i, "bayesp.pa"] <- NA
    summary_table[i, "Rhat_max"] <- NA
    summary_table[i, "Converged"] <- NA
  } # end try
  write.csv(summary_table, file = paste0("Output/summary_table_od.csv"), row.names = FALSE)
} # end testing else statement


} # end species for loop