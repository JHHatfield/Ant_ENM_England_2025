# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
#
# Edited by Jack Hatfield
# 
# Script Description:  
# Prepare the covariate layers for Maxent. Check covariate colinearity, and select a final set for use in the Maxent models

library(terra)
library(raster)
library(tidyverse)
library(tidyterra)
library(data.table)
library(tidyr)
source('source/misc_functions.R')

## Covariates ------------------------------------------
forest_covariates <- rast("../Data/covariates/processed/forest_PCA_30m.tif")
 
#forest_mask_buff <- rast("../Data/covariates/processed/forest_mask_buff_30m.tif") %>%
#       terra::resample(forest_covariates) %>%
#       subst(0, NA)
#writeRaster(forest_mask_buff, '../Data/covariates/processed/forest_mask_buff.tif', overwrite = T)

forest_mask_buff <- rast("../Data/covariates/processed/forest_mask_buff.tif")
 
 
climate_current <- rast("../Data/covariates/processed/current_clim_stack.tif") %>%
       resample(forest_covariates)
 climate_rcp45 <- rast("../Data/covariates/processed/cp45_clim_stack.tif") %>%
       resample(forest_covariates)
 climate_rcp85 <- rast("../Data/covariates/processed/rcp85_clim_stack.tif") %>%
       resample(forest_covariates)

# writeRaster(climate_current, '../Data/covariates/processed/climate_current.tif', overwrite = T)
# writeRaster(climate_rcp45, '../Data/covariates/processed/climate_rcp45.tif', overwrite = T)
# writeRaster(climate_rcp85, '../Data/covariates/processed/climate_rcp85.tif', overwrite = T)

print('climate done')

 topo_covariates <- rast("../Data/covariates/processed/topo_stack_30m.tif") %>%
       resample(forest_covariates)
 topo_covariates$slope <- log(topo_covariates$slope+1)
# writeRaster(topo_covariates, '../Data/covariates/processed/topo_covariates.tif', overwrite = T)
# topo_covariates <- rast("Ant_ENM/data/topo_covariates.tif")
 
 distance_forest <- rast("../Data/covariates/processed/distance_forest_30m.tif") %>%
       terra::resample(forest_covariates) %>%
       +1 %>%
       log()
names(distance_forest) <- "distance_forest"
print('distance_forest done')
 
covariates_stack <- c(forest_covariates, climate_current, climate_rcp45, climate_rcp85, topo_covariates, distance_forest) %>%
       terra::mask(forest_mask_buff)
print('covariates_stack done')

# check and log
range(covariates_stack$current_perc09_total_rain_coldest)
covariates_stack$current_perc09_total_rain_coldest <- log(covariates_stack$current_perc09_total_rain_coldest+1)
range(covariates_stack$current_perc09_total_rain_coldest)

range(covariates_stack$current_perc09_total_rain_hottest)
covariates_stack$current_perc09_total_rain_hottest <- log(covariates_stack$current_perc09_total_rain_hottest+1)
range(covariates_stack$current_perc09_total_rain_hottest)

range(covariates_stack$slope)
covariates_stack$slope <- log(covariates_stack$slope+1)
range(covariates_stack$slope)

range(covariates_stack$distance_forest)
covariates_stack$distance_forest <- log(covariates_stack$distance_forest+1)
range(covariates_stack$distance_forest)

# writeRaster(covariates_stack, '../Data/covariates/processed/full_covariates_stack_future.tif', overwrite = T)

print("Predictors done")

keep <- which(values(forest_mask_buff) == 1)

set.seed(1)
cells <- sample(keep, 100000) # very big so just take a large sample
sample_df  <- terra::extract(covariates_stack, cells)  %>%
  select(current_perc09_total_rain_coldest, current_perc09_total_rain_hottest,
         current_perc01_temp_coldest, current_perc09_temp_hottest,
         current_dry_duration_09perc, current_annual_mean_temp,
         current_annual_sd_temp, current_temp_seasonality, northness, eastness,
         hillshade, slope, forest_PC1, forest_PC2,
         distance_forest) %>%
  as.data.table() %>%
  drop_na()

cat("Number of complete cases:", nrow(sample_df), "\n")

usdm::vif(combined_df)
print('VIF threshold: 5')
usdm::vifstep(combined_df, th = 5)

print('VIF threshold: 7')
usdm::vifstep(combined_df, th = 7)

print('VIF threshold: 10')
usdm::vifstep(combined_df, th = 10)

cor_matrix <- cor(combined_df, use = "complete.obs")

pdf('../Figures/covariate_correlation.pdf', width = 15, height = 15)
corrplot::corrplot(cor_matrix,
                   method = "circle",
                   type = "upper",
                   order = "original",
                   addCoef.col = "black",
                   number.cex = 0.65,
                   tl.cex = 1,
                   tl.offset = 3,
                   tl.col = "black", tl.srt = 45,
                   diag = FALSE,
                   col = colorRampPalette(c("red", "white", "blue"))(200),
                   na.label = NULL,
                   addgrid.col = "grey")
dev.off()

# extract the predictor combination with acceptable VIF - do for each time period
covariates_stack_pVIF_current <- covariates_stack %>%
  select(current_perc09_total_rain_coldest, current_perc09_total_rain_hottest,
         current_perc01_temp_coldest, 
         current_dry_duration_09perc, 
         current_temp_seasonality, northness, eastness,
         hillshade, slope, forest_PC1, forest_PC2,
         distance_forest)

print(names(covariates_stack_pVIF_current))

covariates_stack_pVIF_rcp45 <- covariates_stack %>%
  select(rcp45_perc09_total_rain_coldest, rcp45_perc09_total_rain_hottest,
         rcp45_perc01_temp_coldest, 
         rcp45_dry_duration_09perc, 
         rcp45_temp_seasonality, northness, eastness,
         hillshade, slope, forest_PC1, forest_PC2,
         distance_forest)

print(names(covariates_stack_pVIF_rcp45))

covariates_stack_pVIF_rcp85 <- covariates_stack %>%
  select(rcp85_perc09_total_rain_coldest, rcp85_perc09_total_rain_hottest,
         rcp85_perc01_temp_coldest, 
         rcp85_dry_duration_09perc, 
         rcp85_temp_seasonality, northness, eastness,
         hillshade, slope, forest_PC1, forest_PC2,
         distance_forest)

print(names(covariates_stack_pVIF_rcp85))

writeRaster(covariates_stack_pVIF_current, '../Data/covariates/processed/full_covariates_stack_future_pro_VIF_current.tif', overwrite = T)
writeRaster(covariates_stack_pVIF_rcp45, '../Data/covariates/processed/full_covariates_stack_future_pro_VIF_rcp45.tif', overwrite = T)
writeRaster(covariates_stack_pVIF_rcp85, '../Data/covariates/processed/full_covariates_stack_future_pro_VIF_rcp85.tif', overwrite = T)
