################################################################################
# screate standard grid
################################################################################
#
# Juan Carlos Villaseñor-Derbez
# juancvd@stanford.edu
# Sept 29, 2022
#
# DCreates a standard grid at the 0.01 resolution, using EPSG 4326
#
################################################################################

## SET UP ######################################################################

# Load packages ----------------------------------------------------------------
library(raster)


## PROCESSING ##################################################################

# Build the raster -------------------------------------------------------------
standard_grid <- raster(
  nrows = 2200, ncol = 26000,
  xmn = -130, xmx = -104,
  ymn = 21, ymx = 43,
  crs = "EPSG:4326",
  vals = 0L,
  resolution = 0.01
)

## EXPORT ######################################################################

# Export the raster ------------------------------------------------------------

writeRaster(
  x = standard_grid,
  filename = here::here("Data", "standard_grid.tif"),
  overwrite = T
)
#*** An alternative here, is to simply save the raster as an R data object:
saveRDS(standard_grid, "Data/standard_grid.Rda") # Including the path is in my opinion a better option than using here
#*** NOTE that you can save just about ANY R object in this way, which is quite neat
