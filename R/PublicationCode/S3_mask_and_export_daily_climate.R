# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Pre-process previously downloaded climate rasters, by cropping, masking, and re-projecting

library(tidyverse)
library(tidyterra)
library(terra)
source('source/misc_functions.R')

ROI <- vect('spatial_other/ROI_kmproj.shp')
list <- list.files("covariates/raw/rcp85", full.names = T, recursive = T)

crop_and_export <- function(filepath){
      img <- rast(filepath) %>% 
            terra::project(crs(km_proj)) %>% 
            crop(ROI) %>% 
            mask(ROI) 
      
      writeRaster(img,
                  paste0("covariates/processed/rcp85_masked/",
                         basename(filepath),
                         ".tif"),
                  overwrite = T)
      
}

lapply(list, crop_and_export)
