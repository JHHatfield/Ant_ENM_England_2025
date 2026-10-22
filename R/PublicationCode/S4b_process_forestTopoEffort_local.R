# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Process non-climate data: Topography, forest data, sampling effort

library(tidyverse)
library(tidyterra)
library(terra)
source('source/misc_functions.R')

ROI <- vect('spatial_other/ROI_kmproj.shp') 
ROI_27700 <- vect('spatial_other/ROI_outline_27700.shp') 

# Directory where the large raster files are located
covar_dir <- '/Users/joriswiethase/Library/CloudStorage/GoogleDrive-j.wiethase@gmail.com/My Drive/Ant modelling/Ant_ENM_England_old/covariates'

# Forest ----------------------------------------------------
## Forest characteristics ----------------------------------------------------
## Aggregated from 1m VOM Lidar data layer, Environment Agency
cover_VOM_30m <- rast(file.path(covar_dir, 'raw/cover_VOM_30m_27700.tif')) %>%    
      terra::project(crs(km_proj)) %>%
      mask(ROI)

perc09_height_VOM_30m <- rast(file.path(covar_dir, 'raw/perc09_VOM_30m_27700.tif')) %>% 
      terra::project(crs(km_proj)) %>%
      mask(ROI)

sd_height_VOM_30m <- rast(file.path(covar_dir, 'raw/sd_VOM_30m_27700.tif')) %>%    
      terra::project(crs(km_proj)) %>%
      mask(ROI)

mean_height_VOM_30m <- rast(file.path(covar_dir, 'raw/mean_VOM_30m_27700.tif')) %>%    
      terra::project(crs(km_proj)) %>%
      mask(ROI)

## Distance to forest patches ----------------------------------------------------
forest_mask <- cover_VOM_30m
values(forest_mask) <- ifelse(values(forest_mask) < 0.3, NA, 1)

## Create a distance to forest patch layer
distance_forest <- terra::distance(forest_mask, unit = "m") %>%
      mask(ROI)

## Forest mask layer ----------------------------------------------------
## How far away were most ant records found?
sporadic <- vect('species_data/processed_shp/sporadic_combined.shp')
vals <- terra::extract(distance_forest, sporadic, ID = F)
distance_upper <- quantile(vals$cover_VOM, probs = 0.95)
print(distance_upper)

## 95% of nests are within 67.1 meters of the forest edge. Round up and use as buffer
forest_mask_buff <- buffer(forest_mask, width = 0.07, background = 0)
values(forest_mask_buff) <- ifelse(values(forest_mask_buff) == TRUE, 1, 0)
forest_mask_buff <- mask(forest_mask_buff, ROI)

# How many nests are outside this buffer?
nests_captured <- terra::extract(forest_mask_buff, sporadic, ID = F)
table(nests_captured) # 23 fall outside the buffered zone, 525 are inside

## Final forest stack ----------------------------------------------------
forest_stack_30m <- c(cover_VOM_30m, perc09_height_VOM_30m, sd_height_VOM_30m, mean_height_VOM_30m, 
                      forest_mask_buff)
names(forest_stack_30m) <- c('cover_VOM', 'perc09_height_VOM', 'sd_height_VOM', 'mean_height_VOM',
                             'forest_mask_buff')

# Topography ----------------------------------------------------
# Derived from NASA DEM layer
NASA_dem_30m <- rast(file.path(covar_dir, 'raw/nasa_dem_30m.tif')) %>% 
   crop(ROI_27700) %>% # Use original projection, as elevation units have to be same as map units
   mask(ROI_27700) 

aspect_30m <- terrain(NASA_dem_30m, v = "aspect")
northness_30m <- cos(aspect_30m * pi / 180)
eastness_30m <- sin(aspect_30m * pi / 180)

slope_30m <- terrain(NASA_dem_30m, v = "slope")
hillshade_30m <- shade(slope_30m, aspect_30m, angle=45)

topo_stack_30m <- c(northness_30m, eastness_30m, hillshade_30m, slope_30m) %>%    
   terra::project(crs(km_proj)) 

names(topo_stack_30m) <- c('northness', 'eastness', 'hillshade', 'slope')

# Effort raster ----------------------------------------------------
effort_lgcp_vector <- vect("spatial_other/effort_lgcp_10km.shp") 

rast_OS_grid <- rast(ext = ext(effort_lgcp_vector), crs = crs(effort_lgcp_vector), res = 10, vals = 1)

effort_rast_lgcp_10km <- rasterize(effort_lgcp_vector, rast_OS_grid, field = "days_sampl") %>% 
   crop(ROI) %>% 
   mask(ROI) %>% 
   +1 %>% 
   log() 

# Export ----------------------------------------------------
writeRaster(forest_stack_30m, "covariates/processed/forest_stack_30m.tif", overwrite = TRUE)
writeRaster(distance_forest, "covariates/processed/distance_forest_30m.tif", overwrite = TRUE)
writeRaster(topo_stack_30m, "covariates/processed/topo_stack_30m.tif", overwrite=TRUE)
writeRaster(effort_rast_lgcp_10km, "covariates/processed/effort_rast_lgcp_10km.tif", overwrite=TRUE)


