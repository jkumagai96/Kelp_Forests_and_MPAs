# Date: February 6th 2023
# Author: Joy Kumagai (kumagaij@stanford.edu) 
# Purpose: Bootstrap - randomization approach 
# BIO 202: Ecological Statistics

###### Research Question #######################################################
# Do marine protected areas increase the resilience of kelp forests (increase in kelp area)
# following a marine heat wave?

# Is there more kelp area in marine protected areas? Partial vs. fully protected?

###### Set up ##################################################################
# Load Packages
library(tidyverse)
library(sf)

# Load Data
kelp_data_all <- read.csv("Processed_data/data_tables/kelp_data_all_variables_and_mpa_status_per_year.csv")
points_in_mpas <- read.csv("Processed_data/data_tables/Spatial_intersect_mpas_and_station_points.csv") 

##### Structure Data ###########################################################
kelp_data_south <- kelp_data_all %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full"))) %>% 
  select(PixelID, year, area, mpa_status, region) %>% 
  filter(region == "South_Coast") %>% 
  drop_na() 
  
kelp_data_central <- kelp_data_all %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full"))) %>% 
  select(PixelID, year, area, mpa_status, region) %>% 
  filter(region == "Central_Coast") %>% 
  drop_na() 

kelp_data_north_central <- kelp_data_all %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full"))) %>% 
  select(PixelID, year, area, mpa_status, region) %>% 
  filter(region == "North_Central_Coast") %>% 
  drop_na() 

kelp_data_north <- kelp_data_all %>% 
  mutate(mpa_status = factor(mpa_status, levels = c("None", "Partial", "Full"))) %>% 
  select(PixelID, year, area, mpa_status, region) %>% 
  filter(region == "North_Coast") %>% 
  drop_na() 

a <- kelp_data_all %>% select(region, PixelID) %>% distinct()

points_in_mpas_south <- points_in_mpas %>% 
  left_join(a, by = "PixelID") %>% 
  filter(region == "South_Coast")

points_in_mpas_central <- points_in_mpas %>% 
  left_join(a, by = "PixelID") %>% 
  filter(region == "Central_Coast")

points_in_mpas_north_central <- points_in_mpas %>% 
  left_join(a, by = "PixelID") %>% 
  filter(region == "North_Central_Coast")

points_in_mpas_north <- points_in_mpas %>% 
  left_join(a, by = "PixelID") %>% 
  filter(region == "North_Coast")

##### South Coast ############################################################
## Calculate True Values
# Change global option to not print out message about group summaries 
options(dplyr.summarise.inform = FALSE)

# Calculate average kelp area per category per year 
true_values <- kelp_data_south %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_area = mean(area)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = avg_area) %>% 
  mutate(P_N = None - Partial,
         F_N = None - Full,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

### Bootstrapping 
## Set up:
set.seed(20) # So the results are repeatable

# Set up variables outside of the loops
bootstrap_list <- list()
years <- 1984:2021

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data_south %>% select(-mpa_status)

all_pixels <- kelp_data_south$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_south), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_south$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas_south, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  not_protected <- kelp_w_mpas_r %>% 
    filter(mpa_status == "None") %>% 
    arrange(PixelID)
  
  subset <- kelp_w_mpas_r %>% 
    filter(mpa_status != "None")
  
  for (i in 1:nrow(subset)) {
    status <- subset$mpa_status[i]
    establishment_yr <- subset$Estab_Yr_1[i]
    
    if (subset$year[i] >= establishment_yr) {
      subset$mpa_status[i] <- status
    } else {
      subset$mpa_status[i] <- "None"
      subset$AreaMar_12[i] <- 0
    }
  }
  
  subset <- subset %>% arrange(PixelID)
  new_data <- rbind(not_protected, subset)
  
  # Calculate table 
  values <- new_data %>% 
    group_by(year, mpa_status) %>% 
    summarise(avg_area = mean(area)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = avg_area) %>% 
    mutate(P_N = Partial - None,
           F_N = Full - None,
           F_P = Full - Partial) %>% 
    select(-c(None, Partial, Full)) %>% 
    ungroup()
  
  bootstrap_list[[j]] <- values
}
end <- Sys.time()

start - end

### Processing Results 
bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results_south <- data.frame()

for (i in 1:length(years)) {
  
  # Filter bootstrap values by year
  df <- bootstrap_df %>% 
    filter(year == years[i]) 
  
  tv <- true_values %>% 
    filter(year == years[i])
  
  if (length(unique(df$F_N)) == 1) {
    sig_F_N <- NA } else {
      sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
    }
  
  if (length(unique(df$P_N)) == 1) {
    sig_P_N <- NA } else {
      sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
    }
  
  if (length(unique(df$F_P)) == 1) {
    sig_F_P <- NA } else {
      sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)
    }
  
  t <- data.frame(years[i],
                  sig_F_N,
                  sig_P_N,
                  sig_F_P)
  
  colnames(t) <- c("year", "F_N", "P_N", "F_P")
  
  final_results_south <- rbind(final_results_south, t)
  
}


### Export 
final_results_south

write.csv(final_results_south, 
          "Processed_data/data_tables/bootstrap_south_coast.csv", 
          row.names = F)

rm(kelp_data_south)

##### Central Coast ############################################################
# Calculate average kelp area per category per year 
true_values <- kelp_data_central %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_area = mean(area)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = avg_area) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

### Bootstrapping 
## Set up:
set.seed(20) # So the results are repeatable

# Set up variables outside of the loops
bootstrap_list <- list()
years <- 1984:2021

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data_central %>% select(-mpa_status)

all_pixels <- kelp_data_central$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_central), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_central$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas_central, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  not_protected <- kelp_w_mpas_r %>% 
    filter(mpa_status == "None") %>% 
    arrange(PixelID)
  
  subset <- kelp_w_mpas_r %>% 
    filter(mpa_status != "None")
  
  for (i in 1:nrow(subset)) {
    status <- subset$mpa_status[i]
    establishment_yr <- subset$Estab_Yr_1[i]
    
    if (subset$year[i] >= establishment_yr) {
      subset$mpa_status[i] <- status
    } else {
      subset$mpa_status[i] <- "None"
      subset$AreaMar_12[i] <- 0
    }
  }
  
  subset <- subset %>% arrange(PixelID)
  new_data <- rbind(not_protected, subset)
  
  # Calculate table 
  values <- new_data %>% 
    group_by(year, mpa_status) %>% 
    summarise(avg_area = mean(area)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = avg_area) %>% 
    mutate(P_N = Partial - None,
           F_N = Full - None,
           F_P = Full - Partial) %>% 
    select(-c(None, Partial, Full)) %>% 
    ungroup()
  
  bootstrap_list[[j]] <- values
}
end <- Sys.time()

start - end

### Processing Results 
bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results_central <- data.frame()

for (i in 1:length(years)) {
  
  # Filter bootstrap values by year
  df <- bootstrap_df %>% 
    filter(year == years[i]) 
  
  tv <- true_values %>% 
    filter(year == years[i])
  
  if (length(unique(df$F_N)) == 1) {
    sig_F_N <- NA } else {
      sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
    }
  
  if (length(unique(df$P_N)) == 1) {
    sig_P_N <- NA } else {
      sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
    }
  
  if (length(unique(df$F_P)) == 1) {
    sig_F_P <- NA } else {
      sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)
    }
  
  t <- data.frame(years[i],
                  sig_F_N,
                  sig_P_N,
                  sig_F_P)
  
  colnames(t) <- c("year", "F_N", "P_N", "F_P")
  
  final_results_central <- rbind(final_results_central, t)
  
}


### Export 
final_results_central
write.csv(final_results_central, 
          "Processed_data/data_tables/bootstrap_central_coast.csv", 
          row.names = F)

##### North Central Coast ##############################################################
# Calculate average kelp area per category per year 
true_values <- kelp_data_north_central %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_area = mean(area)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = avg_area) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

### Bootstrapping 
## Set up:
set.seed(20) # So the results are repeatable

# Set up variables outside of the loops
bootstrap_list <- list()
years <- 1984:2021

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data_north_central %>% select(-mpa_status)

all_pixels <- kelp_data_north_central$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_north_central), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_north_central$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas_north_central, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  not_protected <- kelp_w_mpas_r %>% 
    filter(mpa_status == "None") %>% 
    arrange(PixelID)
  
  subset <- kelp_w_mpas_r %>% 
    filter(mpa_status != "None")
  
  for (i in 1:nrow(subset)) {
    status <- subset$mpa_status[i]
    establishment_yr <- subset$Estab_Yr_1[i]
    
    if (subset$year[i] >= establishment_yr) {
      subset$mpa_status[i] <- status
    } else {
      subset$mpa_status[i] <- "None"
      subset$AreaMar_12[i] <- 0
    }
  }
  
  subset <- subset %>% arrange(PixelID)
  new_data <- rbind(not_protected, subset)
  
  # Calculate table 
  values <- new_data %>% 
    group_by(year, mpa_status) %>% 
    summarise(avg_area = mean(area)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = avg_area) %>% 
    mutate(P_N = Partial - None,
           F_N = Full - None,
           F_P = Full - Partial) %>% 
    select(-c(None, Partial, Full)) %>% 
    ungroup()
  
  bootstrap_list[[j]] <- values
}
end <- Sys.time()

start - end

### Processing Results 
bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results_north_central <- data.frame()

for (i in 1:length(years)) {
  
  # Filter bootstrap values by year
  df <- bootstrap_df %>% 
    filter(year == years[i]) 
  
  tv <- true_values %>% 
    filter(year == years[i])
  
  if (length(unique(df$F_N)) == 1) {
    sig_F_N <- NA } else {
      sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
    }
  
  if (length(unique(df$P_N)) == 1) {
    sig_P_N <- NA } else {
      sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
    }
  
  if (length(unique(df$F_P)) == 1) {
    sig_F_P <- NA } else {
      sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)
    }
  
  t <- data.frame(years[i],
                  sig_F_N,
                  sig_P_N,
                  sig_F_P)
  
  colnames(t) <- c("year", "F_N", "P_N", "F_P")
  
  final_results_north_central <- rbind(final_results_north_central, t)
  
}


### Export 
final_results_north_central
write.csv(final_results_north_central, 
          "Processed_data/data_tables/bootstrap_north_central_coast.csv", 
          row.names = F)

##### North Coast ##############################################################
# Calculate average kelp area per category per year 
true_values <- kelp_data_north %>% 
  group_by(year, mpa_status) %>% 
  summarise(avg_area = mean(area)) %>% 
  pivot_wider(names_from = mpa_status,
              values_from = avg_area) %>% 
  mutate(P_N = Partial - None,
         F_N = Full - None,
         F_P = Full - Partial) %>% 
  select(-c(None, Partial, Full))

### Bootstrapping 
## Set up:
set.seed(20) # So the results are repeatable

# Set up variables outside of the loops
bootstrap_list <- list()
years <- 1984:2021

# Remove original mpa_status from kelp_data
kelp_data_r <- kelp_data_north %>% select(-mpa_status)

all_pixels <- kelp_data_north$PixelID %>% unique()

# within the for loop 
start <- Sys.time()
for (j in 1:10000) {
  
  # Sample random pixels that will overlap w/ MPAs 
  r_pixels <- sample(all_pixels, size = nrow(points_in_mpas_north), replace = F) # number of pixels originally overlapping with MPAs 
  
  # Assign the random pixels to the points in the mpas instead of the original values
  points_in_mpas_north$PixelID <- r_pixels
  
  # Join randomized mpa data with kelp data 
  kelp_w_mpas_r <- left_join(kelp_data_r, points_in_mpas_north, by = "PixelID") 
  
  # Ensure that the mpa status accounts for the year of establishment
  kelp_w_mpas_r$mpa_status[is.na(kelp_w_mpas_r$mpa_status)] <- "None"
  
  not_protected <- kelp_w_mpas_r %>% 
    filter(mpa_status == "None") %>% 
    arrange(PixelID)
  
  subset <- kelp_w_mpas_r %>% 
    filter(mpa_status != "None")
  
  for (i in 1:nrow(subset)) {
    status <- subset$mpa_status[i]
    establishment_yr <- subset$Estab_Yr_1[i]
    
    if (subset$year[i] >= establishment_yr) {
      subset$mpa_status[i] <- status
    } else {
      subset$mpa_status[i] <- "None"
      subset$AreaMar_12[i] <- 0
    }
  }
  
  subset <- subset %>% arrange(PixelID)
  new_data <- rbind(not_protected, subset)
  
  # Calculate table 
  values <- new_data %>% 
    group_by(year, mpa_status) %>% 
    summarise(avg_area = mean(area)) %>% 
    pivot_wider(names_from = mpa_status,
                values_from = avg_area) %>% 
    mutate(P_N = Partial - None,
           F_N = Full - None,
           F_P = Full - Partial) %>% 
    select(-c(None, Partial, Full)) %>% 
    ungroup()
  
  bootstrap_list[[j]] <- values
}
end <- Sys.time()

start - end

### Processing Results 
bootstrap_df <- do.call(rbind.data.frame, bootstrap_list)

# Each loop will be a row in the new dataframe
final_results_north <- data.frame()

for (i in 1:length(years)) {
  
  # Filter bootstrap values by year
  df <- bootstrap_df %>% 
    filter(year == years[i]) 
  
  tv <- true_values %>% 
    filter(year == years[i])
  
  if (length(unique(df$F_N)) == 1) {
    sig_F_N <- NA } else {
      sig_F_N <- 1 - ecdf(df$F_N)(tv$F_N)
    }
  
  if (length(unique(df$P_N)) == 1) {
    sig_P_N <- NA } else {
      sig_P_N <- 1 - ecdf(df$P_N)(tv$P_N)
    }
  
  if (length(unique(df$F_P)) == 1) {
    sig_F_P <- NA } else {
      sig_F_P <- 1 - ecdf(df$F_P)(tv$F_P)
    }
  
  t <- data.frame(years[i],
                  sig_F_N,
                  sig_P_N,
                  sig_F_P)
  
  colnames(t) <- c("year", "F_N", "P_N", "F_P")
  
  final_results_north <- rbind(final_results_north, t)
  
}


### Export 
final_results_north
write.csv(final_results_north, 
          "Processed_data/data_tables/bootstrap_north_coast.csv", 
          row.names = F)