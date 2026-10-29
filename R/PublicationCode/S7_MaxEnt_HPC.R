options(java.parameters = "-Xmx64g") 
library(rJava)
library(ENMeval)
library(terra)
library(raster)
library(tidyverse)
library(tidyterra)
library(data.table)
library(dismo)
library(tidyr)
library(parallel)
library(ecospat)

n_cores = detectCores()-1
source('source/misc_functions.R')

## Covariates ------------------------------------------

covariates_stack <- rast('../Data/covariates/processed/full_covariates_stack_future_pro_VIF_current.tif')

# The resulting dataframe (df) contains only the variables with VIF < 10
keep <- c("northness", "eastness", "hillshade", "slope",
          "distance_forest", "forest_PC1", "forest_PC2",
          "current_perc09_total_rain_coldest", 
          "current_perc09_total_rain_hottest", 
          "current_perc01_temp_coldest", 
          "current_dry_duration_09perc",      
          "current_temp_seasonality")

predictors <- subset(covariates_stack, 
                     subset = keep)

range(predictors$slope) # check these variables are already logged
range(predictors$current_perc09_total_rain_coldest)
range(predictors$current_perc09_total_rain_hottest)
range(predictors$distance_forest)

# log those that are not
# predictors$slope <- log(predictors$slope+1)
predictors$current_perc09_total_rain_coldest <- log(predictors$current_perc09_total_rain_coldest+1)
predictors$current_perc09_total_rain_hottest <- log(predictors$current_perc09_total_rain_hottest+1)
# predictors$distance_forest <- log(predictors$distance_forest+1)

print("Predictors done")

for (job in 1:36){

n_background_choices = c(10000, 20000)
thin_choices = c(0, 100, 1000)
seed_choices = c(42, 43, 44)
species_choices = c("Formica rufa", "Formica lugubris")

mult_combs <- crossing(species_choices, n_background_choices, thin_choices, seed_choices) 
comb_values <- mult_combs[job, ]
print(comb_values)

species_choice = comb_values$species_choices
thin_dist = comb_values$thin_choices 
seed = comb_values$seed_choices 
n_background = comb_values$n_background_choices 

# LOAD DATA FILES ------------------------------------------

## Observations --------------------------------------------
sporadic <- read.csv('../Data/species_data/processed_csv/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)

occs <- sporadic %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>% 
      thin_spatial(., dist_meters = thin_dist, seed = seed) %>% 
      as.data.frame(geom = "XY")

unseen_records_vect <- vect('../Data/species_data/processed_shp/exhaustive_combined.shp') %>%
      filter(species == species_choice) 

unseen_records <- unseen_records_vect %>%
      as.data.frame(geom = "XY") %>%
      dplyr::select(x, y) 
gc()

# Background points
bg_bias <- read.csv(paste0('../Data/bg_bias/bg_bias_', gsub(" ", "_", species_choice), "_seed", seed, '_n', n_background, '_thin', thin_dist, '.csv'), row.names = 1)

occs_pred <- terra::extract(predictors, occs, ID = F)
occs <- cbind(as.data.frame(occs, geom = 'xy'), occs_pred)
bg_bias_pred <- terra::extract(predictors, bg_bias, ID = F)
bg_bias <- cbind(bg_bias, bg_bias_pred) 

gc()

set.seed(seed)         

print(paste("Number of occurrence points:", nrow(occs)))
print(paste("Number of background points:", nrow(bg_bias)))         


cb2 <- get.checkerboard2(occs[,1:2], bg=bg_bias[,1:2], env=raster::raster(predictors), aggregation.factor = c(10,10))

# perform model with or without hinge features depending on if there are enough occurrence points

if(nrow(occs)>70){
e.mx <- ENMevaluate(occs = occs, bg = bg_bias,
                    algorithm = 'maxent.jar', partitions = "user", user.grp = cb2,
                    tune.args = list(fc = c("Q", "LQ", "LQP", "LQT", "QH", "LQH"), rm = seq(1, 5, 0.5)),
                    parallel = T, numCores = n_cores,
                    other.settings = list(path = '../tmp_files',
                                          other.args = list(pictures = F,
                                                            outputdirectory = '../tmp_files',
                                                            writeclampgrid = F,
                                                            writemess = F,
                                                            plots = F,
                                                            outputgrids = F,
                                                            visible = F)))
}else{
e.mx <- ENMevaluate(occs = occs, bg = bg_bias,
                    algorithm = 'maxent.jar', partitions = "user", user.grp = cb2,
                    tune.args = list(fc = c("Q", "LQ", "LQP", "LQT"), rm = seq(1, 5, 0.5)),
                    parallel = T, numCores = n_cores, 
                    other.settings = list(path = '../tmp_files',
                                          other.args = list(pictures = F,
                                                            outputdirectory = '../tmp_files',
                                                            writeclampgrid = F,
                                                            writemess = F,
                                                            plots = F,
                                                            outputgrids = F,
                                                            visible = F)))
}
print("Maxent model done")

# Model selection
res <- eval.results(e.mx)

# Using or.100p, CBI and AICc stepwise
opt.mult.step <- res %>% 
   dplyr::filter(AICc == min(AICc)) %>%  # Akaike Information Criterion corrected for small sample size
   dplyr::filter(or.10p.avg == min(or.10p.avg)) %>% 
   dplyr::filter(auc.val.avg == max(auc.val.avg)) 
  
print(opt.mult.step)

if(NROW(opt.mult.step) > 1){
   opt.mult.step <- opt.mult.step[1, ]
}
best_model <- eval.models(e.mx)[[opt.mult.step$tune.args]]

saveRDS(best_model, paste0("../model_out/", gsub(" ", "_", species_choice), '_thin', thin_dist, '_All_30m_seed_', seed, '_', n_background, ".RDS"))
saveRDS(e.mx, paste0("../model_out/ENMeval_", gsub(" ", "_", species_choice), '_thin', thin_dist, '_All_30m_seed_', seed, '_', n_background, ".RDS"))

var_importance <- data.frame(best_model@results[grep("permutation.importance", rownames(best_model@results)),]) 
names(var_importance) <- "perm_imp"
print(arrange(var_importance, desc(perm_imp)))

pdf(paste0("../figures/effects_", gsub(" ", "_", species_choice), '_thin', thin_dist, '_All_30m_seed_', seed, '_', n_background, ".pdf"), width = 14, height = 14)
par(mfrow = c(1, 1))
plot_maxent_effects(best_model)
dev.off()

pdf(paste0("../figures/envSim_", gsub(" ", "_", species_choice), '_thin', thin_dist, '_All_30m_seed_', seed, '_', n_background, ".pdf"), width = 14, height = 14)
par(mfrow = c(1, 1))
evalplot.envSim.hist(sim.type = "mess", ref.data = "occs", occs.z = occs, 
                     bg.z = bg_bias, occs.grp = cb2$occs.grp, bg.grp = cb2$bg.grp)
dev.off()

}
