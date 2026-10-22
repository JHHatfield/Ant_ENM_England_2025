# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Apply further processing to the pre-processed climate layers

library(tidyverse)
library(tidyterra)
library(terra)
library(lubridate)
source('source/misc_functions.R')


# https://data.ceda.ac.uk/badc/deposited2021/chess-scape
# https://catalogue.ceda.ac.uk/uuid/8194b416cbee482b89e0dfbe17c5786c/
# Climate downloaded from: https://nex-gddp-cmip6.s3.us-west-2.amazonaws.com/index.html#NEX-GDDP-CMIP6/EC-Earth3/
# Chose EC-Earth3 due to good performance in Europe in Palmer et al. 2023

for(period in c('current', 'rcp45', 'rcp85')){
      # Import daily climate and set parameters ----------------------------------------------------
      if(period == 'current'){
            start_year = 2010
            end_year = 2013
            temp_files <- list.files(path = paste0("covariates/processed/", period, "_masked"), pattern = "1km_mon", full.names = TRUE)
            rain_files <- list.files(path = paste0("covariates/processed/", period, "_masked"), pattern = "1km_day", full.names = TRUE)
            
            temp_stack <- rast(temp_files)
            names(temp_stack) <- sprintf("%d%02d", rep(start_year:end_year, each = 12), 1:12)
            
      } else {
            start_year = 2070
            end_year = 2073
            temp_files <- list.files(path = paste0("covariates/processed/", period, "_masked"), pattern = "tas_uk_1km_monthly", full.names = TRUE)
            rain_files <- list.files(path = paste0("covariates/processed/", period, "_masked"), pattern = "pr_uk_1km_daily", full.names = TRUE)
            
            temp_stack <- rast(temp_files) 
            temp_stack <- temp_stack %>% subset(which(year(time(.)) %in% seq(start_year, end_year, 1)))
            names(temp_stack) <- sprintf("%d%02d", rep(start_year:end_year, each = 12), 1:12)
            values(temp_stack) <- values(temp_stack)-273.15
      }
      
      rainfall_threshold <- 1 # Rain threshold in mm, used to calculate dry spell duration
      n_quarters = 4  # The number of individual coldest and hottest quarters to look for, 4 since we have data for 4 years
      
      # Identify the hottest and coldest quarters ----------------------------------------------------
      quarters <- ceiling(1:nlyr(temp_stack) / 3)
      temp_stack_qtr <- tapp(temp_stack, quarters, fun = median)
      avg_qtr_temp <- global(temp_stack_qtr, stat = 'median', na.rm = T)
      
      ## Coldest quarters ----------------------------------------------------
      coldest_qtrs_indices <- order(avg_qtr_temp$mean)[1:n_quarters]
      coldest_qtrs_dates <- lapply(coldest_qtrs_indices, function(q) {
            year <- start_year + ((q - 1) %/% 4)
            qtr <- (q - 1) %% 4 + 1
            return(paste(year, qtr, sep = "-Q"))
      })

      ## Hottest quarters ----------------------------------------------------
      hottest_qtrs_indices <- order(avg_qtr_temp$mean, decreasing = TRUE)[1:n_quarters]
      hottest_qtrs_dates <- lapply(hottest_qtrs_indices, function(q) {
            year <- start_year + ((q - 1) %/% 4)
            qtr <- (q - 1) %% 4 + 1
            return(paste(year, qtr, sep = "-Q"))
      })
      
      # Get total rainfall of coldest and hottest quarters ----------------------------------------------------
      get_quarterly_rainfall <- function(qtrs_dates){
            quarterly_rainfall_rasters <- list()
            for (i in seq_along(qtrs_dates)) {
                  q = qtrs_dates[[i]]
                  year_qtr <- strsplit(q, "-Q")[[1]]
                  year <- as.numeric(year_qtr[1])
                  qtr <- as.numeric(year_qtr[2])
                  
                  q_start <- (qtr - 1) * 3 + 1
                  q_end <- q_start + 2
                  
                  q_rainfall_sum <- NULL
                  
                  for (m in q_start:q_end) {
                        month_pattern <- sprintf("%04d%02d", year, m)
                        month_rain_files <- grep(month_pattern, rain_files, value = TRUE)
                        if (length(month_rain_files) > 0) {
                              month_rain_stack <- rast(month_rain_files)
                              if(period != 'current'){
                                    values(month_rain_stack) <- values(month_rain_stack)*86400
                              } 
                              month_rain_sum <- sum(month_rain_stack, na.rm = TRUE)
                              if (is.null(q_rainfall_sum)) {
                                    q_rainfall_sum <- month_rain_sum
                              } else {
                                    q_rainfall_sum <- q_rainfall_sum + month_rain_sum
                              }
                        }
                  }
                  quarterly_rainfall_rasters[[i]] <- q_rainfall_sum
            }
            
            total_rainfall_09perc <- app(rast(quarterly_rainfall_rasters), fun = calc_quantiles)
            return(total_rainfall_09perc)
      }
      
      get_quarterly_temperature <- function(qtrs_dates, quantile = 0.9){
            quarterly_temp_rasters <- list()
            for (i in seq_along(qtrs_dates)) {
                  q = qtrs_dates[[i]]
                  year_qtr <- strsplit(q, "-Q")[[1]]
                  year <- as.numeric(year_qtr[1])
                  qtr <- as.numeric(year_qtr[2])
                  
                  q_start <- (qtr - 1) * 3 + 1
                  q_end <- q_start + 2
                  
                  month_temp_rasters <- list()
                  
                  for (k in 1:length(q_start:q_end)) {
                        month <- (q_start:q_end)[k]
                        month_pattern <- sprintf("%04d%02d", year, month)
                        month_temp_sub <- temp_stack %>% subset(month_pattern)
                        month_temp_rasters[[k]] <- month_temp_sub
                  }
                  quarterly_temp_rasters[[i]] <- rast(month_temp_rasters) %>% median()
            }
            
            temp_qtr_09perc <- app(rast(quarterly_temp_rasters), fun = calc_quantiles, probs = quantile) 
            return(temp_qtr_09perc)
      }
      
      perc09_total_rain_coldest <- get_quarterly_rainfall(coldest_qtrs_dates)
      perc09_total_rain_hottest <- get_quarterly_rainfall(hottest_qtrs_dates)
      
      perc01_temp_coldest <- get_quarterly_temperature(coldest_qtrs_dates, quantile = 0.1)
      perc09_temp_hottest <- get_quarterly_temperature(hottest_qtrs_dates)
      
      # Get dryspell duration ----------------------------------------------------
      year_files <- list()
      dryspell_rasters <- list()
      
      # Organize files by year
      for (file in rain_files) {
            if(period == 'current'){
                  year <- substr(basename(file),  31, 34) 
            } else {
                  year <- substr(basename(file),  53, 56) 
            }
            
            if (!year %in% names(year_files)) {
                  year_files[[year]] <- c(file)
            } else {
                  year_files[[year]] <- c(year_files[[year]], file)
            }
      }
      
      calculate_longest_dryspell <- function(year_stack, threshold) {
            dry_day_count <- terra::init(year_stack %>% subset(1), 0) %>% 
                  mask(year_stack %>% subset(1))
            
            for (i in 1:nlyr(year_stack)) {
                  daily_rain <- terra::subset(year_stack, i)
                  is_dry <- daily_rain < threshold
                  # Increment dry day count where it's dry, reset to 0 where it's not
                  dry_day_count <- ifel(is_dry, dry_day_count + 1, 0)
                  
                  if (i == 1) {
                        max_dry_spell <- dry_day_count
                  } else {
                        max_dry_spell <- max(max_dry_spell, dry_day_count)
                  }
            }
            return(max_dry_spell)
      }
      
      # Process each year
      for (year in names(year_files)) {
            year_stack <- rast(year_files[[year]])
            if(period != 'current'){
                  values(year_stack) <- values(year_stack)*86400
            } 
            names(year_stack) <- paste0("day", 1:nlyr(year_stack))
            hottest_qtr <- hottest_qtrs_dates[grep(year, hottest_qtrs_dates)]
            qtr <- as.numeric(substr(hottest_qtr, 7, 7))
            start_month <- (qtr - 1) * 3 + 1
            end_month <- start_month + 2
            start_date <- as.Date(paste(year, start_month, "01", sep="-"))
            end_date <- as.Date(paste(year, end_month + 1, "01", sep="-")) - 1
            days_in_hottest_qtr <- yday(seq.Date(start_date, end_date, by="day"))
            year_stack_hottest <- year_stack %>% subset(days_in_hottest_qtr)
            longest_dryspell <- calculate_longest_dryspell(year_stack_hottest, rainfall_threshold)
            dryspell_rasters[[year]] <- longest_dryspell
      }
      
      dry_duration_09perc <- app(rast(dryspell_rasters), fun = calc_quantiles)
      
      # Annual temperature mean, sd (seasonality) and range ----------------------------------------------------
      yearly_temp <- tapp(temp_stack, index =  unique(substr(names(temp_stack), 1, 4)), fun = mean, na.rm = TRUE)
      
      annual_mean_temp <- mean(yearly_temp, na.rm = T)
      year_sd_temp_list <- list()
      year_temp_range_list <- list()
      
      for (year in unique(substr(names(temp_stack), 1, 4))) {
            year_sd <- temp_stack[[substr(names(temp_stack), 1, 4) == year]] %>% 
                  app(., fun = sd, na.rm = TRUE)
            year_sd_temp_list[[as.character(year)]] <- year_sd
            temp_max <- temp_stack[[substr(names(temp_stack), 1, 4) == year]] %>% 
                  max()
            temp_min <- temp_stack[[substr(names(temp_stack), 1, 4) == year]] %>% 
                  min()
            year_range <- temp_max - temp_min
            year_temp_range_list[[as.character(year)]] <- year_range
      }
      
      mean_year_sd_temp <- mean(rast(year_sd_temp_list))
      mean_year_range_temp <- mean(rast(year_temp_range_list))
      
      clim_stack <- c(perc09_total_rain_coldest, perc09_total_rain_hottest,
                      perc01_temp_coldest, perc09_temp_hottest,
                      dry_duration_09perc, annual_mean_temp, mean_year_sd_temp, mean_year_range_temp)
      
      names(clim_stack) <- paste0(period, '_', c('perc09_total_rain_coldest', 'perc09_total_rain_hottest',
                                                 'perc01_temp_coldest', 'perc09_temp_hottest',
                                                 'dry_duration_09perc', 'annual_mean_temp', 'annual_sd_temp', 'temp_seasonality'))
      
      writeRaster(clim_stack, paste0("covariates/processed/", period, "_clim_stack.tif"), overwrite = TRUE)
}

