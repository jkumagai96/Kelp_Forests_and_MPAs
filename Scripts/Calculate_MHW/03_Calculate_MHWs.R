# Date: November 22nd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Calculate MHWs
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# Load Packages 
library(dplyr) # For basic data manipulation
library(ggplot2) # For visualising data
library(heatwaveR) # For detecting MHWs

# Load data
SST_data_full <- readRDS("Processed_data/SST/SST_1984_2021.rds")

##### Format Data ##############################################################
# Filter
SST_data <- SST_data_full %>% 
  filter(lat <= 32.6) %>% 
  filter(lon >= -117.3) %>% # Filter to test the code 
  filter(t >= "1982-01-01")

##### Check about the time series continuity 
SST_data %>% 
  mutate(year = format(t, "%Y")) %>% 
  group_by(year) %>% 
  summarize(length(unique(t)))
# Not all days have data every year 

##### Declare Functions ########################################################
annual_intensity <- function(df){
  # First calculate the climatologies
  clim <- heatwaveR::ts2clm(data = df, climatologyPeriod = c("1985-01-01", "2014-12-31"))
  
  # Then the events
  mhw <- detect_event(data = clim)
  
  # Return only the annual intensity, mhw events, and mhw days
  mhw$climatology %>% 
    dplyr::select(-c(doy, threshCriterion, durationCriterion, event)) %>% # Q: seas vs. thresh? 
    # Define cumulative intensity as only experienced during mhw or just above seasonality?
    drop_na() %>% 
    mutate(year = format(t, "%Y")) %>% 
    mutate(intensity = temp - thresh) %>% 
    group_by(year) %>%
    summarise(
      mhw_events = length(unique(event_no)),
      mhw_days = length(t),
      mhw_int_cumulative = sum(intensity),
      .groups = "drop"
    )
}

##### Calculate Heatwaves ######################################################
MHW_dplyr <- SST_data %>% 
  # Then we group the data by the 'lon' and 'lat' columns
  group_by(lon, lat) %>% 
  # Then we run our MHW detecting function on each group
  group_modify(~annual_intensity(.x))


SST_data %>% 
  mutate(year = format(t, "%Y")) %>% 
  group_by(year) %>% 
  summarize(length(unique(t)))


##### Example ##################################################################
head(heatwaveR::sst_WA)

# Detect the events in a time series
ts <- ts2clm(sst_WA, climatologyPeriod = c("1982-01-01", "2011-12-31"))
mhw <- detect_event(ts)

mhw$climatology %>% 
  dplyr::select(-c(doy, threshCriterion, durationCriterion, event)) %>% # Q: seas vs. thresh? 
  # Define cumulative intensity as only experienced during mhw or just above seasonality?
  drop_na() %>% 
  mutate(year = format(t, "%Y")) %>% 
  mutate(intensity = temp - thresh) %>% 
  group_by(year) %>%
  summarise(
    mhw_events = length(unique(event_no)),
    mhw_days = length(t),
    mhw_int_cumulative = sum(intensity),
    .groups = "drop"
  )




# # Calculate climatolgoies
# annual_intensity <- function(df){
#   # First calculate the climatologies
#   clim <- heatwaveR::ts2clm(data = df, climatologyPeriod = c("1985-01-01", "2014-12-31"))
# 
#   # Then the events
#   mhw <- detect_event(data = clim)
# 
#   # Return only the annual intensity, mhw events, and mhw days
#   mhw$climatology %>%
#     dplyr::select(-c(doy, threshCriterion, durationCriterion, event)) %>% # Q: seas vs. thresh?
#     # Define cumulative intensity as only experienced during mhw or just above seasonality?
#     drop_na() %>%
#     mutate(year = format(t, "%Y")) %>%
#     mutate(intensity = temp - thresh) %>%
#     group_by(year) %>%
#     summarise(
#       mhw_events = length(unique(event_no)),
#       mhw_days = length(t),
#       mhw_int_cumulative = sum(intensity),
#       .groups = "drop"
#     )
# }
# 
# 
# ##### Processing ###############################################################
# 
# # Save each unqiue strip of longitude
# for(i in 1:length(unique(SST_data$lon))){
#   SST_data_sub <- SST_data %>%
#     filter(lon == unique(lon)[i])
#   saveRDS(object = SST_data_sub, file = paste0("Processed_data/SST/longitude_strip/SST_data_lon_",i,".Rds"))
# }
# 
# # The 'dplyr' wrapper function to pass to 'plyr'
# dplyr_wraper <- function(file_name){
#   MHW_dplyr <- readRDS(file_name) %>%
#     group_by(lon, lat) %>%
#     group_modify(~annual_intensity(.x))
# }
# 
# # Create a vector of the files we want to use
# SST_data_files <- dir("Processed_data/SST/longitude_strip", pattern = "SST_data_lon_*", full.names = T)
# 
# # Use 'plyr' technique to run 'dplyr' technique with multiple cores
# system.time(
#   MHW_result <- plyr::ldply(SST_data_files, .fun = dplyr_wraper, .parallel = T)
# )
# 
# # Save for later use as desired
# saveRDS(MHW_result, "Processed_data/SST/MHW_result_test.Rds")
# 
# ##### Attempt 2 ################################################################
# # Define processing funciton ---------------------------------------------------
# anual_intensity <- function(data) {
#   dplyr::select(-c(doy, threshCriterion, durationCriterion, event)) %>% # Q: seas vs. thresh?
#     # Define cumulative intensity as only experienced during mhw or just above seasonality?
#     drop_na() %>%
#     mutate(year = format(t, "%Y")) %>%
#     mutate(intensity = temp - thresh) %>%
#     group_by(year) %>%
#     summarise(
#       mhw_events = length(unique(event_no)),
#       mhw_days = length(t),
#       mhw_int_cumulative = sum(intensity),
#       .groups = "drop"
#     )
# }
# 
# # Get HWs ----------------------------------------------------------------------
# mhw <- SST_dat %>%
#   nest(data = c(t, temp)) %>%
#   mutate(
#     ts = map(data, ts2clm, climatologyPeriod = c("1984-01-01", "2013-12-31")),
#     mhw = map(ts, detect_event),
#     summary = map(mhw, anual_intensity)
#   )
# 
# # Save for later use as desired
# saveRDS(mhw, "Processed_data/SST/MHW_result_test.Rds")
