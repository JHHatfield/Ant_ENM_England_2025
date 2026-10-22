# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Apply PCA to forest structure data, to isolate distinct forest types

library(tidyverse)
library(tidyterra)
library(terra)

# Load data ----------------------------------------------------
forest_stack <- rast('Ant_ENM/data/forest_stack_30m.tif')

forest_mask_buff <- forest_stack$forest_mask_buff %>%
   subst(0, NA)

forest_stack$perc09_height_VOM_sqrt <- sqrt(forest_stack$perc09_height_VOM)
forest_stack$sd_height_VOM_sqrt <- sqrt(forest_stack$sd_height_VOM)

## Forest PCA ----------------------------------------------------
forest_stack_PCA <- forest_stack %>%
      tidyterra::select(cover_VOM, perc09_height_VOM_sqrt, sd_height_VOM_sqrt) %>% 
      mask(forest_mask_buff)

rpc_forest <- terra::prcomp(forest_stack_PCA, center = T, scale. = T)
summary(rpc_forest)

# pdf("Ant_ENM/figures/forest_PCA_scree.pdf", width = 7, height = 7)
# fviz_eig(rpc_forest, addlabels = TRUE)
# dev.off()
# 
# pdf("Ant_ENM/figures/forest_PCA_variables.pdf", width = 7, height = 7)
# fviz_pca_var(rpc_forest, col.var = "black")
# dev.off()

# forest_cumulative <- rpc_forest$sdev^2 / sum(rpc_forest$sdev^2)
# forest_cumulative_df <- data.frame(cumsum(forest_cumulative)) %>%
#       mutate(PCA_comp = paste0("PC", rownames(.)),
#              prop_var = forest_cumulative) %>%
#       rename(cum_var = cumsum.forest_cumulative.)

# forest_scores_df <- data.frame(rpc_forest$x)

forest_loadings_df <- data.frame(rpc_forest$rotation) %>%
      mutate(PCA_comp = rownames(.)) %>%
      dplyr::select(PCA_comp, everything())
print(forest_loadings_df)

forest_PC_stack <- predict(forest_stack, rpc_forest) %>%
      subset(1:2) %>% 
      c(., forest_mask_buff)

names(forest_PC_stack) <- gsub("PC", "forest_PC", names(forest_PC_stack))


# Export ----------------------------------------------------
# write.csv(forest_cumulative_df, "covariates/processed/HPC_forest_cumulative_df.csv")
# fwrite(forest_scores_df, "covariates/processed/HPC_forest_scores_df.csv")
write.csv(forest_loadings_df, "Ant_ENM/data/HPC_forest_loadings_df.csv")
writeRaster(forest_PC_stack, paste0("Ant_ENM/data/forest_PCA_30m.tif"), overwrite=TRUE)


