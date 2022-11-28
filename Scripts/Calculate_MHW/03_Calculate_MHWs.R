# Date: November 22nd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Calculate MHWs
# BIO 202: Ecological Statistics

##### Set up ###################################################################
# Load Packages 
library(dplyr) # For basic data manipulation
library(ggplot2) # For visualizing data
library(heatwaveR) # For detecting MHWs

# Load data
SST_data_full <- readRDS("Processed_data/SST/SST_1984_2021.rds")

##### Format Data ##############################################################
# Filter
# SST_data <- SST_data_full %>%
#   filter(lat <= 33.2) %>%
#   filter(lon >= -119)  # Filter to test the code

SST_data <- SST_data_full 

##### Check about the time series continuity 
SST_data %>% 
  mutate(year = format(t, "%Y")) %>% 
  group_by(year) %>% 
  summarize(length(unique(t))) %>% view()
# Not all days have data every year 

##### Declare Functions ########################################################
annual_intensity_hw <- function(df){
  # First calculate the climatologies
  clim <- heatwaveR::ts2clm(data = df, climatologyPeriod = c("1985-01-01", "2014-12-31"))
  # Then the events
  mhw <- detect_event(data = clim)
  
  # Return only the annual intensity, mhw events, and mhw days
  mhw$climatology %>% 
    dplyr::select(-c(doy, threshCriterion, durationCriterion, event)) %>% 
    # Define cumulative intensity as only experienced during mhw 
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


annual_intensity_cold_spells <- function(df){
  # First calculate the climatologies
  clim <- heatwaveR::ts2clm(data = df, climatologyPeriod = c("1985-01-01", "2014-12-31"), pctile = 10)
  # Then the events
  cs <- detect_event(data = clim, coldSpells = TRUE)
  
  # Return only the annual intensity, cold spell events, and cold spell days
  cs$climatology %>% 
    dplyr::select(-c(doy, threshCriterion, durationCriterion, event)) %>% 
    # Define cumulative intensity as only experienced during cold spell 
    drop_na() %>% 
    mutate(year = format(t, "%Y")) %>% 
    mutate(intensity = temp - thresh) %>% 
    group_by(year) %>%
    summarise(
      cs_events = length(unique(event_no)),
      cs_days = length(t),
      cs_int_cumulative = sum(intensity),
      .groups = "drop"
    )
}

##### Calculate Heatwaves and Cold Spells ######################################
### Heatwaves
MHW_df <- SST_data %>% 
  # Then we group the data by the 'lon' and 'lat' columns
  group_by(lon, lat) %>% 
  # Then we run our MHW detecting function on each group
  group_modify(~annual_intensity_hw(.x)) %>% 
  ungroup()

### Cold spells
CS_df <- SST_data %>% 
  # Then we group the data by the 'lon' and 'lat' columns
  group_by(lon, lat) %>% 
  # Then we run our MHW detecting function on each group
  group_modify(~annual_intensity_cold_spells(.x)) %>% 
  ungroup()


##### Export ###################################################################
saveRDS(object = MHW_df,
        file = "Processed_data/SST/MHW_1984_2021.rds")

saveRDS(object = CS_df,
        file = "Processed_data/SST/CS_1984_2021.rds")

