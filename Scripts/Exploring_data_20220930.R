# Date: September 30th 2022
# Author: Joy Kumagai (kumagaij@stanford.edu)
# Purpose: Exploring Data for Final Project
# BIO 202: Ecological Statistics

###### Load in Data and packages #####
library(tidyverse)
library(sf)
library(ncdf4)
library(raster)

mpas_all <- st_read("Data/Marine Protected Areas/NOAA_MPAI_2020_IUCN_gdb/NOAA_MPAI_v2020.gdb")
kelp_cover_file <- "Data/Kelp/knb-lter-sbc.74.17/LandsatKelpBiomass_2022_Q2_withmetadata.nc"
kelp_cover <- nc_open(kelp_cover_file)
print(kelp_cover)
dname <- "Year"

##### Exploring Data #####
colnames(mpas_all)
mpas <- mpas_all %>% 
  filter(State == "CA" & Estab_Yr <= 2014) # 2014 is the year the marine heat wave started

glimpse(mpas)
unique(mpas$Prot_Lvl)
plot(mpas[,1])

# Kelp cover
print(kelp_cover)
dname <- "Year"

lon <- ncvar_get(kelp_cover,"longitude") # X variable
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(kelp_cover,"latitude") # Y variable 
nlat <- dim(lat)
head(lat)


year <- ncvar_get(kelp_cover,"year") # T variable
quarter <- ncvar_get(kelp_cover, "quarter") # Also a T variable 
area <- ncvar_get(kelp_cover, "area") # matrix
area_se <- ncvar_get(kelp_cover, "area_se") # matrix 
biomass <- ncvar_get(kelp_cover, "biomass")
biomass_se <- ncvar_get(kelp_cover, "biomass_se")

# create a simple dataframe 
# area of kelp in Q1 of 1994 
area_kelp_Q1_1994 <- data.frame(lon, lat, area[,1]) 
colnames(area_kelp_Q1_1994) <- c("longitude", "latitude", "Area_m2_1994_Q1")
sum(area_kelp_Q1_1994$Area_m2_1994_Q1, na.rm = T)


# Load data from tiff 
kelp <- raster("Data/Kelp/LandsatKelp_Quarterly_2022_2.tif")


