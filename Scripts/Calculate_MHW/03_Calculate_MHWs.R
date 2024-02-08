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
# Years
year <- 1983:2021
yrs <- tibble(year)
# Marine cold spell detect event and summarize per year function
mhw_event_only <- function(df){
  
  # First calculate the climatologies
  clim <- ts2clm(data = df, climatologyPeriod = c("1983-01-01", "2012-12-31"))
  
  # Then the events
  event <- detect_event(data = clim)
  
  #Then the cummulative intensity and days per year 
  out <- event$climatology %>%
    mutate(intensity = temp - seas, 
           year = year(t)) %>% 
    filter(event == TRUE) %>% 
    group_by(year, .drop = FALSE) %>%
    summarise(total_icum = sum(intensity, na.rm = TRUE),
              total_days = n()) 
  
  annual <- left_join(yrs, out, by = "year") %>% 
    mutate(total_icum = ifelse(is.na(total_icum), 0, total_icum),
           total_days = ifelse(is.na(total_days), 0, total_days))
  
  # Return only the event metric dataframe of results
  return(annual)
}

# Marine cold spell detect event and summarize per year function 
mcs_event_only <- function(df){
  # First calculate the climatologies
  clim <- ts2clm(data = df, climatologyPeriod = c("1983-01-01", "2012-12-31"), pctile = 10)
  
  # Then the cold spell events
  event <- detect_event(data = clim, coldSpells = TRUE)
  
  #Then the cummulative intensity and days per year 
  out <- event$climatology %>%
    mutate(intensity = seas - temp, 
           year = year(t)) %>% 
    filter(event == TRUE) %>% 
    group_by(year, .drop = FALSE) %>%
    summarise(total_icum = sum(intensity, na.rm = TRUE),
              total_days = n()) 
  
  annual <- left_join(yrs, out, by = "year") %>% 
      mutate(total_icum = ifelse(is.na(total_icum), 0, total_icum),
             total_days = ifelse(is.na(total_days), 0, total_days))
    
  # Return only the event metric dataframe of results
  return(annual)
}


##### Calculate Heatwaves and Cold Spells ######################################
### Heatwaves
MHW_df <- SST_data %>% 
  # Then we group the data by the 'lon' and 'lat' columns
  group_by(lon, lat) %>% 
  # Then we run our MHW detecting function on each group
  group_modify(~mhw_event_only(.x)) 

print('marine heat waves calculated')

### Cold spells
CS_df <- SST_data %>% 
  # Then we group the data by the 'lon' and 'lat' columns
  group_by(lon, lat) %>% 
  # Then we run our marine cold spell detecting function on each group
  group_modify(~mcs_event_only(.x)) 

print('cold spells calculated')

##### Export ###################################################################
saveRDS(object = MHW_df,
        file = "Processed_data/SST/MHW_1983_2021.rds")

saveRDS(object = CS_df,
        file = "Processed_data/SST/CS_1983_2021.rds")

print('script finished')

##### Checking Values ##########################################################
MHW_df %>% 
  #filter(lat <= 34.4486) %>% # southern california
  group_by(year) %>% 
  summarize(values = mean(total_days)) %>% 
  ggplot(aes(x = year, y = values)) +
  geom_point()
  