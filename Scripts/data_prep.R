# install.packages("AHMbook", repos = NULL, type="source")
# install.packages("plotrix")
# library(AHMbook)

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

