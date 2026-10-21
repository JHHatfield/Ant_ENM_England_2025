# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Download climate data for any time period. 

library(httr)
library(rvest)
library(lubridate)

# Set year range
start_year <- 2010
end_year <- 2013

# Base URL for the directory, un-comment as needed
base_url <- "https://dap.ceda.ac.uk/badc/deposited2021/chess-scape/data/rcp45_bias-corrected/01/daily/pr/"
download_dir <- "'/Users/joriswiethase/Library/CloudStorage/GoogleDrive-j.wiethase@gmail.com/My Drive/Ant modelling/Ant_ENM_new/covariates/raw/rcp45'"

# base_url <- "https://dap.ceda.ac.uk/badc/deposited2021/chess-scape/data/rcp85_bias-corrected/01/daily/pr/"
# download_dir <- "'/Users/joriswiethase/Library/CloudStorage/GoogleDrive-j.wiethase@gmail.com/My Drive/Ant modelling/Ant_ENM_new/covariates/raw/rcp85'"


# Retrieve directory listing
response <- GET(base_url)
if (http_status(response)$category != "Success") {
      stop("Failed to access directory. Check URL and permissions.")
}

# Parse HTML content
page_content <- read_html(content(response, as = "text"))

# Extract all .nc file links
file_links <- page_content %>% 
      html_nodes("a") %>% 
      html_attr("href") %>% 
      grep("\\.nc$", ., value = TRUE)

files_to_download <- file_links[grepl(paste0(paste0('_', seq(start_year, end_year, 1)), collapse = '|'), file_links)]

# Download files using wget
for (file in files_to_download) {
      download_url <- paste0(base_url, file)
      system2("wget", args = c(
            "-e", "robots=off",
            "--no-parent",
            "--mirror",
            "-q", 
            "-P", download_dir,
            shQuote(download_url)
      ))
      message("Downloaded: ", file)
}

