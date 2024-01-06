### LIBRARY IMPORTS ###

library(ggplot2)
library(data.table)
library(matrixStats)
library(glmnet)
library(hydroGOF)
library(clue)
library(dplyr)

### SET PARAMETERS ###
file2Open <- "/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/imports/chattahoochee_2017_2024.csv"
clusters2Open <- "exports/cluster_regress/lag_4/lag_4_RATING_clusters_calculated.RData"
saveName <- "CHAT_SSC.csv"


### SET DIRECTORIES ###

print('SETTING UP IMPORTS')
# Set root directory
wd_root <- "/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/"
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

sentinel <- read.csv(file2Open)
sentinel$date <- as.Date(sentinel$date)

setnames(sentinel, old = c('B2_median',
                  'B3_median',
                  'B4_median',
                  'B5_median',
                  'B6_median',
                  'B7_median',
                  'B8_median',
                  'B8A_median',
                  'B11_median',
                  'B12_median'), 
                  new = c('B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'))

setnames(sentinel, old = colnames(sentinel), new=toupper(colnames(sentinel)))

# Site numbers stored as numbers in gee so need to replace leading zeros
# sentinel$site_no <- unlist(lapply(sentinel$site_no, pad0))
# sentinel$site_no <- as.character(sentinel$site_no)

### APPLY CLUSTERING ###

# Load Clusters -> "clusters_calculated"
load(file = clusters2Open)

c_bands<- c('B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B11','B12')

# For every 2km segment, cluster based on median reflectance
cluster_medians <- data.table(matrix(nrow = 0, ncol = length(c_bands) + 1))
names(cluster_medians) <- c('DISTANCE_KM', c_bands)
unique_dist <- unique(sentinel$DISTANCE_KM)

pb <- txtProgressBar(0, length(unique_dist), style = 3)
for (i in seq_along(unique_dist)){
  setTxtProgressBar(pb, i)
  dist <- unique_dist[i]
  # Extract data
  dist_data <- setDT(sentinel[sentinel$DISTANCE_KM == dist,])[, ..c_bands] %>% as.matrix()
  # Reduce by median
  dist_median <- colMedians(dist_data, na.rm=T)
  # Create table to hold site data
  dist_df <- data.table(t(data.table(dist_median)))
  names(dist_df) <- c_bands
  dist_df$DISTANCE_KM <- dist
  # Generate regression
  cluster_medians <- rbind(cluster_medians, dist_df)
}

# Assign cluster to each distance
cluster_medians$cluster <- cl_predict(clusters_calculated, cluster_medians)

# Extract cluster for each dist
cluster_dist <- cluster_medians[, .(DISTANCE_KM, cluster)]
setnames(cluster_dist, old=colnames(cluster_dist), new=toupper(colnames(cluster_dist)))

# Quick plot to ensure majority cluster is clear
hist(cluster_dist$CLUSTER)

# Plot cluster over distance
ggplot(cluster_dist, aes(x=DISTANCE_KM, y=CLUSTER)) +
geom_line()

# Apply majority cluster to the whole dataset
# sentinel$cluster <- majCluster
sentinel <- left_join(sentinel, cluster_dist, by = "DISTANCE_KM")

### Apply Regression ###

# If already clustered read file in
# sentinel <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/exports/transects/CHATTAHOOCHEE_2017_2024.csv")
# Add all band variables
sentinelBands <- setDT(copy(sentinel))[, ":="(
  # Add squared columns
  # B1.2 = B1^2,
  B2.2 = B2^2,
  B3.2 = B3^2,
  B4.2 = B4^2,
  B5.2 = B5^2,
  B6.2 = B6^2,
  B7.2 = B7^2,
  B8.2 = B8^2,
  B8A.2 = B8A^2,
  # B9.2 = B9^2,
  B11.2 = B11^2,
  B12.2 = B12^2,
  # Add square root columns
  # B1.0.5 = B1^0.5,
  B2.0.5 = B2^0.5,
  B3.0.5 = B3^0.5,
  B4.0.5 = B4^0.5,
  B5.0.5 = B5^0.5,
  B6.0.5 = B6^0.5,
  B7.0.5 = B7^0.5,
  B8.0.5 = B8^0.5,
  B8A.0.5 = B8A^0.5,
  # B9.0.5 = B9^0.5,
  B11.0.5 = B11^0.5,
  B12.0.5 = B12^0.5,
  # Add band ratios
  # B2.B1=B2/B1,
  # B3.B1=B3/B1,
  # B4.B1=B4/B1,
  # B5.B1=B5/B1,
  # B6.B1=B6/B1,
  # B7.B1=B7/B1,
  # B8.B1=B8/B1,
  # B8A.B1=B8A/B1,
  # B9.B1=B9/B1,
  # B11.B1=B11/B1,
  # B12.B1=B12/B1,
  
  B3.B2=B3/B2,
  B4.B2=B4/B2,
  B5.B2=B5/B2,
  B6.B2=B6/B2,
  B7.B2=B7/B2,
  B8.B2=B8/B2,
  B8A.B2=B8A/B2,
  # B9.B2=B9/B2,
  B11.B2=B11/B2,
  B12.B2=B12/B2,
  
  B4.B3=B4/B3,
  B5.B3=B5/B3,
  B6.B3=B6/B3,
  B7.B3=B7/B3,
  B8.B3=B8/B3,
  B8A.B3=B8A/B3,
  # B9.B3=B9/B3,
  B11.B3=B11/B3,
  B12.B3=B12/B3,
  
  B5.B4=B5/B4,
  B6.B4=B6/B4,
  B7.B4=B7/B4,
  B8.B4=B8/B4,
  B8A.B4=B8A/B4,
  # B9.B4=B9/B4,
  B11.B4=B11/B4,
  B12.B4=B12/B4,
  
  B6.B5=B6/B5,
  B7.B5=B7/B5,
  B8.B5=B8/B5,
  B8A.B5=B8A/B5,
  # B9.B5=B9/B5,
  B11.B5=B11/B5,
  B12.B5=B12/B5,
  
  B7.B6=B7/B6,
  B8.B6=B8/B6,
  B8A.B6=B8A/B6,
  # B9.B6=B9/B6,
  B11.B6=B11/B6,
  B12.B6=B12/B6,
  
  B8.B7=B8/B7,
  B8A.B7=B8A/B7,
  # B9.B7=B9/B7,
  B11.B7=B11/B7,
  B12.B7=B12/B7,
  
  B8A.B8=B8A/B8,
  # B9.B8=B9/B8,
  B11.B8=B11/B8,
  B12.B8=B12/B8,
  
  # B9.B8A=B9/B8A,
  B11.B8A=B11/B8A,
  B12.B8A=B12/B8A,
  
  # B11.B9=B11/B9,
  # B12.B9=B12/B9,
  
  B12.B11=B12/B11
)]


# For just 1 cluster per river
# Load regression
# load(file = "exports/cluster_regress/lag_4/RATING_regression/ssc_lm_4.RData")
# sentinel_pred <- copy(sentinel)
# # Should only be 1 cluster per river
# cluster <- 4
# sentinel_pred$CLUSTER_SEL <- cluster
# ssc <- predict(
#     ssc_lm,
#     newx = as.matrix(sentinelBands %>% select(starts_with("B") & -B2_COUNT)),
#     s = "lambda.min"
#     )
# sentinel_pred$LOG10_SSC_MGL <- ssc
# sentinel_pred$SSC_MGL <- 10^ssc

# For multiple clusters per river
sentinel_pred <- copy(sentinel)
sentinel_pred$PRED_LOG10_SSC_MGL <- NA
# For one cluster uncomment
cluster <- 4
sentinelBands$CLUSTER <- cluster
sentinel_pred$CLUSTER <- cluster

path <- "exports/cluster_regress/lag_4/RATING_regression/ssc_lm_"
for (clust in unique(sentinelBands$CLUSTER)) {
  # Load in "ssc_lm"
  load(file = paste0(path, clust, ".RData"))
  ssc <- predict(
    ssc_lm,
    newx = as.matrix(sentinelBands %>% filter(CLUSTER == clust) %>% select(starts_with("B") & -B2_COUNT)),
    s = "lambda.min"
    )
  sentinel_pred[sentinel_pred$CLUSTER == clust, ]$PRED_LOG10_SSC_MGL = ssc
}
sentinel_pred$PRED_SSC_MGL <- 10^sentinel_pred$PRED_LOG10_SSC_MGL
sentinel_pred <- sentinel_pred %>% filter(PRED_LOG10_SSC_MGL <= 4)

# Quick sanity check that everything looks okay
ggplot(sentinel_pred, aes(x = DISTANCE_KM, y = PRED_SSC_MGL)) +
  stat_summary(fun = mean, geom = "line")

#sentinel_pred %>% select(DISTANCE_KM, CLUSTER) %>% distinct() %>% filter(CLUSTER != 4)

# ggsave(filename = paste0(wd_transect, "transect_sanity_check.pdf"),
#        width = 12, height = 8)

# Save Data
write.csv(sentinel_pred, file = "/Users/ilanvalencius/Documents/River-Sed-Manuscript/figure-data/CHAT4.csv", row.names = FALSE)
