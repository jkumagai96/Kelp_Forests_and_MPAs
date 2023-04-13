# Date: November 22nd 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Calculate MHWs
# BIO 202: Ecological Statistics

print('begining of script')
##### Set up ###################################################################
# Load Packages 
library(tidyverse) # For basic data manipulation

library(heatwaveR) # For detecting MHWs

print('packages loaded')

# Load data
SST_data <- readRDS("Processed_data/SST/SST_1983_2021.rds")

print('data loaded')
##### Declare Functions ########################################################
# Marine cold spell detect event and summarize per year function
mhw_event_only <- function(df){
  # First calculate the climatologies
  clim <- ts2clm(data = df, climatologyPeriod = c("1983-01-01", "2012-12-31"))
  # Then the events
  event <- detect_event(data = clim)
  #Then the annual
  annual <- block_average(event, x = t, y = temp, report = "full")
  # Return only the event metric dataframe of results
  return(annual)
}

# Marine cold spell detect event and summarize per year functino 
mcs_event_only <- function(df){
  # First calculate the climatologies
  clim <- ts2clm(data = df, climatologyPeriod = c("1983-01-01", "2012-12-31"))
  # Then the cold spell events
  event <- detect_event(data = clim, coldSpells = TRUE)
  #Then the annual
  annual <- block_average(event, x = t, y = temp, report = "full")
  # Return only the event metric dataframe of results
  return(annual)
}

print('function loaded')
##### Calculate Heatwaves and Cold Spells ######################################
### Heatwaves
MHW_df <- SST_data %>% 
  # Then we group the data by the 'lon' and 'lat' columns
  group_by(lon, lat) %>% 
  # Then we run our MHW detecting function on each group
  group_modify(~mhw_event_only(.x)) %>% 
  # Select columns we are interested in 
  dplyr::select(lat, lon, year, count, total_days, total_icum) %>%
 # Change NAs to 0's when there are no detected events
  dplyr::mutate(total_days = replace_na(total_days, 0)) %>% 
  dplyr::mutate(total_icum = replace_na(total_icum, 0))

print('marine heat waves calculated')

### Cold spells
CS_df <- SST_data %>% 
  # Then we group the data by the 'lon' and 'lat' columns
  group_by(lon, lat) %>% 
  # Then we run our marine cold spell detecting function on each group
  group_modify(~mcs_event_only(.x)) %>% 
  # Select columns we are interested in 
  dplyr::select(lat, lon, year, count, total_days, total_icum) %>% 
  # Change NAs to 0's when there are no detected events
  dplyr::mutate(total_days = replace_na(total_days, 0)) %>% 
  dplyr::mutate(total_icum = replace_na(total_icum, 0))

print('cold spells calculated')

##### Export ###################################################################
saveRDS(object = MHW_df,
        file = "Processed_data/SST/MHW_1983_2021.rds")

saveRDS(object = CS_df,
        file = "Processed_data/SST/CS_1983_2021.rds")

print('script finished')
