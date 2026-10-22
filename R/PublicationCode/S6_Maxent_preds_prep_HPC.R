library(terra)
library(raster)
library(tidyverse)
library(tidyterra)
library(data.table)
library(tidyr)
source('Ant_ENM/source/misc_functions.R')

## Covariates ------------------------------------------
# forest_covariates <- rast("Ant_ENM/data/forest_PCA_30m.tif")
# 
# forest_mask_buff <- rast("Ant_ENM/data/forest_mask_buff_30m.tif") %>%
#       terra::resample(forest_covariates) %>%
#       subst(0, NA)
# writeRaster(forest_mask_buff, 'Ant_ENM/data/forest_mask_buff.tif', overwrite = T)
# forest_mask_buff <- rast("Ant_ENM/data/forest_mask_buff.tif")
# 
# 
# climate_current <- rast("Ant_ENM/data/current_clim_stack.tif") %>%
#       resample(forest_covariates)
# climate_rcp45 <- rast("Ant_ENM/data/rcp45_clim_stack.tif") %>%
#       resample(forest_covariates)
# climate_rcp85 <- rast("Ant_ENM/data/rcp85_clim_stack.tif") %>%
#       resample(forest_covariates)
# writeRaster(climate_current, 'Ant_ENM/data/climate_current.tif', overwrite = T)
# writeRaster(climate_rcp45, 'Ant_ENM/data/climate_rcp45.tif', overwrite = T)
# writeRaster(climate_rcp85, 'Ant_ENM/data/climate_rcp85.tif', overwrite = T)
# 
# print('climate done')
# topo_covariates <- rast("Ant_ENM/data/topo_stack_30m.tif") %>%
#       resample(forest_covariates)
# topo_covariates$slope <- log(topo_covariates$slope+1)
# writeRaster(topo_covariates, 'Ant_ENM/data/topo_covariates.tif', overwrite = T)
# topo_covariates <- rast("Ant_ENM/data/topo_covariates.tif")
# 
# distance_forest <- rast("Ant_ENM/data/distance_forest_30m.tif") %>%
#       terra::resample(forest_covariates) %>%
#       +1 %>%
#       log()
# names(distance_forest) <- "distance_forest"
# print('distance_forest done')
# 
# covariates_stack <- c(forest_covariates, climate_current, climate_rcp45, climate_rcp85, topo_covariates, distance_forest) %>%
#       terra::mask(forest_mask_buff)
# print('covariates_stack done')
# writeRaster(covariates_stack, 'Ant_ENM/data/full_covariates_stack.tif', overwrite = T)
# print("Predictors done")

combined_df <- rast('Ant_ENM/data/full_covariates_stack.tif') %>%
   select(current_perc09_total_rain_coldest, current_perc09_total_rain_hottest,
          current_perc01_temp_coldest, current_perc09_temp_hottest,
          current_dry_duration_09perc, current_annual_mean_temp,
          current_annual_sd_temp, current_temp_seasonality, northness, eastness,
          hillshade, slope, forest_PC1, forest_PC2,
          distance_forest) %>%
   as.data.table() %>%
   drop_na()

usdm::vif(combined_df)
print('VIF threshold: 5')
usdm::vifstep(combined_df, th = 5)

print('VIF threshold: 7')
usdm::vifstep(combined_df, th = 7)

print('VIF threshold: 10')
usdm::vifstep(combined_df, th = 10)

cor_matrix <- cor(combined_df, use = "complete.obs")

pdf('Ant_ENM/figures/covariate_correlation.pdf', width = 15, height = 15)
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



