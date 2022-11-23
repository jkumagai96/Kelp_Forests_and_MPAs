# Date: November 16th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Calculate marine heat wave
# BIO 202: Ecological Statistics

##### Set up ###################################################################
library(dplyr)
library(ggplot2)
library(heatwaveR)
library(lubridate)

##### Example first ############################################################
head(heatwaveR::sst_WA)

# Detect the events in a time series
ts <- ts2clm(sst_WA, climatologyPeriod = c("1982-01-01", "2011-12-31")) # Exactly thirty years 
mhw <- detect_event(ts)

# View just a few metrics
mhw$event %>% 
  dplyr::ungroup() %>%
  dplyr::select(event_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) %>% 
  dplyr::arrange(-intensity_max) %>% 
  head(5)

# Calculate number of mhw events per year 
mhw$event$date_peak %>% 
  format(format = "%Y") %>% 
  table() %>% 
  as.data.frame() %>% 
  rename(year = ".", n_mhw = Freq)


# Calculate cummulative intensity (degrees x days) per year 
mhw$event %>% 
  mutate(start_yr <- format(date_start))
