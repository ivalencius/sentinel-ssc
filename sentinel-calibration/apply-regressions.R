### LIBRARY IMPORTS ###
library(dataRetrieval)
library(nhdplusTools)
# library(tidyhydat)

library(readr)
library(readxl)

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(ggpubr)
library(gstat)
library(ggspatial)
library(svglite)
library(plotly)
library(gridExtra)
library(grid)
library(scales)
library(clue)

library(data.table)
library(gdata) # for cbindX
library(plyr) # for round any
library(gtools) # smartbind

library(dplyr)
library(tidyverse)
library(tidyquant)
library(tidyr)
library(broom)
library(modelr)

library(scales)
library(kdensity)
library(NbClust)
library(zoo)
library(segmented)
library(lubridate)
library(reshape2)
library(matrixStats)
library(smoother)

library(glmnet)
library(boot)
library(kernelboot)
library(np)

library(automap)
library(sp)
library(USAboundaries)
library(sf)
library(rgeos)
library(raster)
library(rgdal)
library(maptools)
library(PBSmapping)
library(dict)
library(clue)

### SET PARAMETERS ###
# file2Open <- "C:/Users/ilanv/Desktop/sentinel-ssc/imports/sentinel_2_SR__transect_2017_2022_20scale.csv"
file2Open <- "C:/Users/ilanv/Desktop/sentinel-ssc/flint_river/sentinel_flint_2017_2024.csv"
clusters2Open <- "C:/Users/ilanv/Desktop/sentinel-ssc/exports/cluster_regress/"
saveName <- "FLINT_20m_SSC.csv"

# Need to alter band regressors to match regression used
load(file = paste0(clusters2Open, "cluster_regressors.RData"))
load(file = paste0(clusters2Open, "cluster_bands.RData"))

### SET DIRECTORIES ###

print('SETTING UP IMPORTS')
# Set root directory
# wd_root <- "D:/valencig/Thesis/sentinel-ssc/sentinel-calibration"
wd_root <- "C:/Users/ilanv/Desktop/sentinel-ssc"
setwd(wd_root)

# Imports folder (store all import files here)
wd_imports <- paste0(wd_root, "/imports/")

# Exports folder (save all figures, tables here)
wd_exports <- paste0(wd_root, "/exports/")

# Store all transect data here
wd_transect <- paste0(wd_exports, "transects/")

### HElPER FUNCTION ###

# From https://github.com/evandethier/satellite-ssc/blob/master/landsat-calibration/landsat-57-calibration.R
fancy_scientific <- function(l) { 
  # turn in to character string in scientific notation 
  l <- log10(l)
  # return(parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep="")))
  return(parse(text = paste("10^",as.character(l),sep = "")))
}

# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_exports, wd_transect)
for(i in seq_along(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}
### IMPORT SENTINEL-2 DATA ###

# helper function -> pads a character string with leading 0s (for USGS codes)
# Columns of Sentinel-2 data export we want to remove
# cols2Go <- c(
#   "system.index",
  # ".geo",
#   "lat",
#   "lon",
#   "snow_ice",
#   "solar_az",
#   "solar_zen",
#   "thin_cirrus_percentage",
#   "water_median",
#   "sensor"
  # )
# Import seperate Sentinel-2 .csv files

# for (file in gee_files){
sentinel <- read.csv(file2Open)
# Remove system.index and .geo columns
# gee_clean <- gee_raw[ , !(names(gee_raw) %in% cols2Go)]
# Convert string to date object
sentinel[,'date'] <- as.Date(sentinel$date)

setnames(sentinel,
          old = c('B1_median',
                  'B2_median',
                  'B3_median',
                  'B4_median',
                  'B5_median',
                  'B6_median',
                  'B7_median',
                  'B8_median',
                  'B8A_median',
                  'B9_median',
                  'B11_median',
                  'B12_median'), new = c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12'))
# sentinel <- rbind(sentinel, setDT(gee_clean))
# }

# Site numbers stored as numbers in gee so need to replace leading zeros
# sentinel$site_no <- unlist(lapply(sentinel$site_no, pad0))
# sentinel$site_no <- as.character(sentinel$site_no)

### APPLY CLUSTERING ###

# Load Clusters
load(file = paste0(clusters2Open, 'clusters_calculated.RData'))

c_bands<- c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')

# For every 2km segment, cluster based on median reflectance
cluster_medians <- data.table(matrix(nrow = 0, ncol = length(c_bands) + 1))
names(cluster_medians) <- c('distance_km', c_bands)
unique_dist <- unique(sentinel$distance_km)

pb <- txtProgressBar(0, length(unique_dist), style = 3)
for (i in seq_along(unique_dist)){
  setTxtProgressBar(pb, i)
  dist <- unique_dist[i]
  # Extract data
  dist_data <- setDT(sentinel[sentinel$distance_km == dist,])[, ..c_bands] %>% as.matrix()
  # Reduce by median
  dist_median <- colMedians(dist_data, na.rm=T)
  # Create table to hold site data
  dist_df <- data.table(t(data.table(dist_median)))
  names(dist_df) <- c_bands
  dist_df$distance_km <- dist
  # Generate regression
  cluster_medians <- rbind(cluster_medians, dist_df)
}

# Assign cluster to each distance
cluster_medians$cluster <- cl_predict(clusters_calculated, cluster_medians)

# Extract cluster for each dist
cluster_dist <- cluster_medians[, .(distance_km, cluster)]

getMode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Quick plot to ensure majority cluster is clear
hist(cluster_dist$cluster)

# Get majority cluster
majCluster <- getMode(cluster_dist$cluster)

# Apply majority cluster to the whole dataset
# sentinel$cluster <- majCluster
sentinel <- left_join(sentinel, cluster_dist, by = ("distance_km" = "distance_km"))

### Apply Regression ###

# Add all band variables
sentinelBands <- setDT(copy(sentinel))[, ":="(
  # Add squared columns
  B1.2 = B1^2,
  B2.2 = B2^2,
  B3.2 = B3^2,
  B4.2 = B4^2,
  B5.2 = B5^2,
  B6.2 = B6^2,
  B7.2 = B7^2,
  B8.2 = B8^2,
  B8A.2 = B8A^2,
  B9.2 = B9^2,
  B11.2 = B11^2,
  B12.2 = B12^2,
  # Add square root columns
  B1.0.5 = B1^0.5,
  B2.0.5 = B2^0.5,
  B3.0.5 = B3^0.5,
  B4.0.5 = B4^0.5,
  B5.0.5 = B5^0.5,
  B6.0.5 = B6^0.5,
  B7.0.5 = B7^0.5,
  B8.0.5 = B8^0.5,
  B8A.0.5 = B8A^0.5,
  B9.0.5 = B9^0.5,
  B11.0.5 = B11^0.5,
  B12.0.5 = B12^0.5,
  # Add band ratios
  B2.B1=B2/B1,
  B3.B1=B3/B1,
  B4.B1=B4/B1,
  B5.B1=B5/B1,
  B6.B1=B6/B1,
  B7.B1=B7/B1,
  B8.B1=B8/B1,
  B8A.B1=B8A/B1,
  B9.B1=B9/B1,
  B11.B1=B11/B1,
  B12.B1=B12/B1,
  
  B3.B2=B3/B2,
  B4.B2=B4/B2,
  B5.B2=B5/B2,
  B6.B2=B6/B2,
  B7.B2=B7/B2,
  B8.B2=B8/B2,
  B8A.B2=B8A/B2,
  B9.B2=B9/B2,
  B11.B2=B11/B2,
  B12.B2=B12/B2,
  
  B4.B3=B4/B3,
  B5.B3=B5/B3,
  B6.B3=B6/B3,
  B7.B3=B7/B3,
  B8.B3=B8/B3,
  B8A.B3=B8A/B3,
  B9.B3=B9/B3,
  B11.B3=B11/B3,
  B12.B3=B12/B3,
  
  B5.B4=B5/B4,
  B6.B4=B6/B4,
  B7.B4=B7/B4,
  B8.B4=B8/B4,
  B8A.B4=B8A/B4,
  B9.B4=B9/B4,
  B11.B4=B11/B4,
  B12.B4=B12/B4,
  
  B6.B5=B6/B5,
  B7.B5=B7/B5,
  B8.B5=B8/B5,
  B8A.B5=B8A/B5,
  B9.B5=B9/B5,
  B11.B5=B11/B5,
  B12.B5=B12/B5,
  
  B7.B6=B7/B6,
  B8.B6=B8/B6,
  B8A.B6=B8A/B6,
  B9.B6=B9/B6,
  B11.B6=B11/B6,
  B12.B6=B12/B6,
  
  B8.B7=B8/B7,
  B8A.B7=B8A/B7,
  B9.B7=B9/B7,
  B11.B7=B11/B7,
  B12.B7=B12/B7,
  
  B8A.B8=B8A/B8,
  B9.B8=B9/B8,
  B11.B8=B11/B8,
  B12.B8=B12/B8,
  
  B9.B8A=B9/B8A,
  B11.B8A=B11/B8A,
  B12.B8A=B12/B8A,
  
  B11.B9=B11/B9,
  B12.B9=B12/B9,
  
  B12.B11=B12/B11
)]

bands2Reg <- cluster_bands[[sentinelBands$cluster[1]]]

# Should only be 1 cluster per river
cluster <- sentinelBands$cluster[1]

ssc_pred <- predict(
    cluster_regressors[[cluster]],
      newx = as.matrix(sentinelBands %>% dplyr::select(all_of(bands2Reg))),
      s = "lambda.min"
    )

sentinel$log10_SSC_mgL <- ssc_pred
sentinel$SSC_mgL <- 10^ssc_pred

# Quick sanity check that everything looks okay
sanityCheck <- ggplot(sentinel, aes(x = distance_km, y = SSC_mgL)) +
  stat_summary(fun = mean, geom = "line")

ggsave(sanityCheck, filename = paste0(wd_transect, "transect_sanity_check.pdf"),
       width = 12, height = 8)

# Save Data
write.csv(sentinel, file = paste0(wd_transect, saveName), row.names = FALSE)
