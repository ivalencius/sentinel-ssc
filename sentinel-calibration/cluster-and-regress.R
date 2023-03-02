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

wd_figures <- paste0(wd_exports, "figures/")
wd_cluster <- paste0(wd_exports, "cluster_regress/")
# Sub folders
# wd_extent<- paste0(wd_exports, 'station_transects/')
wd_gee <- paste0(wd_exports, "GEE_raw/")

### HElPER FUNCTION ###

# From https://github.com/evandethier/satellite-ssc/blob/master/landsat-calibration/landsat-57-calibration.R
fancy_scientific <- function(l) { 
  # turn in to character string in scientific notation 
  l <- log10(l)
  # return(parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep="")))
  return(parse(text = paste("10^",as.character(l),sep = "")))
} 

# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_exports, wd_figures, wd_gee, wd_cluster)
                         # , wd_exports_gc,wd_station_standalone, 
                         # wd_standalone_models, wd_standalone_figures, wd_autocorrelation)
for(i in seq_along(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}

### IMPORT IN SITU DATA ###

# Load in-situ data
ratingSSC <- read.csv(file = paste(wd_exports, "ratingSSCQ.csv")) # rating curve derived
inSituSSC <- read.csv(file = paste(wd_exports, "ONEDAYSAMPLE_SSC.csv", sep="")) # in-situ only

# Extract only log10SSC, site_no, and sample_dt
sscCols <- c("site_no", "log10_SSC_mgL", "sample_dt")
rSSC <- data.frame(
  site_no = as.character(ratingSSC[, sscCols]$site_no),
  log10_SSC_mgL = as.numeric(ratingSSC[, sscCols]$log10_SSC_mgL),
  sample_dt = as.Date(ratingSSC[, sscCols]$sample_dt)
)

iSSC <- data.frame(
  site_no = as.character(inSituSSC[, sscCols]$site_no),
  log10_SSC_mgL = as.numeric(inSituSSC[, sscCols]$log10_SSC_mgL),
  sample_dt = as.Date(inSituSSC[, sscCols]$sample_dt)
)

rm(ratingSSC)
rm(inSituSSC)

### IMPORT SENTINEL-2 DATA ###

# helper function -> pads a character string with leading 0s (for USGS codes)
# pad0 <- function(x) {
#   while (nchar(x) < 8) {
#     x <- paste0("0", x)
#   }
#   return(x)
# }
lag_days <- 4
gee_files <- list.files(path = wd_gee, pattern = 'sentinel*', full.names = TRUE)
# Columns of Sentinel-2 data export we want to remove
cols2Go <- c(
  "system.index",
  ".geo",
  "lat",
  "lon",
  "snow_ice",
  "solar_az",
  "solar_zen",
  "thin_cirrus_percentage",
  "water_median",
  "sensor",
  "B2_count"
  )
# Import seperate Sentinel-2 .csv files
sentinel <- data.table()
for (file in gee_files){
  gee_raw <- read.csv(file)
  # Remove system.index and .geo columns
  gee_clean <- gee_raw[ , !(names(gee_raw) %in% cols2Go)]
  # Convert string to date object
  gee_clean[,'date'] <- as.Date(gee_clean$date)
  setnames(gee_clean, 
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
  sentinel <- rbind(sentinel, setDT(gee_clean))
}
rm(gee_raw)
rm(gee_clean)

# Site numbers stored as numbers in gee so need to replace leading zeros
# sentinel$site_no <- unlist(lapply(sentinel$site_no, pad0))
sentinel$site_no <- as.character(sentinel$site_no)

### Harmonize USGS AND SENTINEL-2 DATA ###

# Harmonize Rating Curve Derived Samples
rHarmonized <- setDT(sentinel)[
  # Create lead lag times and add to DF
  , ":="(
    match_dt_start = date - lag_days,
    match_dt_end = date + lag_days
    )][
    # Match by dates inside lead-lag range
    setDT(rSSC)[, ":="(match_dt = sample_dt)],
    on = .(
      site_no == site_no,
      match_dt_start <= match_dt,
      match_dt_end >= match_dt
      )
    ][!is.na(B1)][,
    lag_days := as.numeric(difftime(sample_dt, date), "days")
    ][, ":="(
    # Remove duplicated columns
    match_dt_start = NULL, match_dt_end = NULL
    )]
          
# Harmonize InSitu Samples
iHarmonized <- setDT(sentinel)[
  # Create lead lag times and add to DF
  , ":="(
    match_dt_start = date - lag_days,
    match_dt_end = date + lag_days
    )][
    # Match by dates inside lead-lag range
    setDT(iSSC)[, ":="(match_dt = sample_dt)],
    on = .(
      site_no == site_no,
      match_dt_start <= match_dt,
      match_dt_end >= match_dt
      )
    ][!is.na(B1)][,
    lag_days := as.numeric(difftime(sample_dt, date), "days")
    ][, ":="(
    # Remove duplicated columns
    match_dt_start = NULL, match_dt_end = NULL
    )]

# ISSUE: in-situ samples not stored for sites with no linear regression
missingSites <- unique(iHarmonized$site_no)[!(unique(iHarmonized$site_no) %in% unique(rHarmonized$site_no))]
missingData <- iHarmonized %>%
  filter(site_no %in% missingSites)
rHMerge <- rbind(rHarmonized, missingData)

# Reduce harmonized samples (rating curve derived SSC)
rHClean <- rHMerge %>%
  # For each image --> isolate one with minimum lag day
  group_by(site_no, product_id) %>%
  slice(which.min(abs(lag_days))) %>%
  ungroup() %>%
  # Now there can be multiple images per day --> reduce to one
  group_by(site_no, sample_dt) %>%
  summarise(across(where(is.numeric), mean)) %>%
  ungroup()

# Reduce harmonized samples (in-situ derived SSC)
iHClean <- iHarmonized %>% 
  # For each image --> isolate one with minimum lag day
  group_by(site_no, product_id) %>%
  slice(which.min(abs(lag_days))) %>% 
  ungroup() %>%
  # Now there can be multiple images per day --> reduce to one
  group_by(site_no, sample_dt) %>%
  summarise(across(where(is.numeric), mean)) %>%
  ungroup()

# Save data
write.csv(rHClean, file = paste0(wd_cluster, "RATING_TRAINING_SAMPLES.csv"), row.names = FALSE)
write.csv(iHClean, file = paste0(wd_cluster, "INSITU_TRAINING_SAMPLES.csv"), row.names = FALSE)

### FROM NOW ON THE PIPELINE IS ONLY WORKING WITH RATING CURVE DERIVED DATA ###

### CLUSTER STATIONS ###

# Bands used for regression
cluster_regressors <- c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')

# Extract median spectral information at each site
cluster_medians <- data.table(matrix(nrow=0, ncol=length(cluster_regressors)+1))
names(cluster_medians) <- c('site_no', cluster_regressors)
unique_sites <- unique(as.character(rHClean$site_no))
pb <- txtProgressBar(0, length(unique_sites), style = 3)
for (i in seq_along(unique_sites)){
  setTxtProgressBar(pb, i)
  site <- unique_sites[i]
  # Extract data
  site_data <- as.data.table(
    rHClean[rHClean$site_no == site,]
    )[, ..cluster_regressors] %>% as.matrix()
  # Reduce by median
  site_median <- colMedians(site_data, na.rm = TRUE) # avoid NA in any column
  # Create table to hold site data
  site_df <- data.table(t(data.table(site_median)))
  names(site_df) <- cluster_regressors
  site_df$site_no <- site
  # Generate regression
  cluster_medians <- rbind(cluster_medians, site_df)
}

# 6 clusters --> following Dethier et. al (2020)
clusters_calculated <- kmeans(cluster_medians[, ..cluster_regressors], centers = 6,
                              nstart = 100, iter.max = 50)

# Compute cluster centers
cluster_centers <- clusters_calculated$centers
# Assign cluster to each site
cluster_medians$cluster <- clusters_calculated$cluster
# Extract cluster for each station
cluster_stations <- cluster_medians[, .(site_no, cluster)]
rHClustered <- left_join(rHClean, cluster_stations, by = c("site_no" = "site_no"))

# Visualize Clusters
ssc_categories <- c(0,50,100,250,500,750,1e6)
# ssc_categories <- c(0,50,100,200,500,1e6)
# ssc_categories <- c(0,10,25,50,75,100,150,200,250,300,350, 400, 450, 500,600, 700, 800,900,1000,1100,1500, 1e6)

# Generate SSC labels as 'low value' - 'high value'
ssc_category_labels <- paste0(ssc_categories[-length(ssc_categories)],'-',c(ssc_categories[-1]))
# Make highest SSC category "> highest value"
ssc_category_labels[length(ssc_category_labels)] <- paste0('> ', ssc_categories[length(ssc_category_labels)])

viz_df <- as.data.table(rHClustered)[
  ,':='(cluster_sel = cluster,
        # # Categorize SSC value as one of selected categories
        ssc_category = cut(10^log10_SSC_mgL,
                           breaks = ssc_categories,
                           labels = ssc_category_labels))][]

# Generate median B,G,R for each SSC category and each cluster or site
ssc_category_color <- viz_df %>% group_by(cluster_sel, ssc_category) %>% 
  summarise_at(vars("B4", "B3", "B2"), mean)

raster_color_types <- geom_raster(aes(fill = rgb(B4/2200,B3/2200,B2/2200))) # true color

cluster_viz <- ggplot(ssc_category_color, aes(x = cluster_sel, y = ssc_category)) +
  raster_color_types +
  scale_fill_identity() +
    theme_classic() +
  # scale_x_continuous(expand_scale(add = c(0,0))) + 
  # scale_y_discrete(expand_scale(mult = c(0,0))) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(
    y = 'SSC range (mg/L)',
    x = 'River grouping'
    )
ggsave(cluster_viz, filename = paste0(wd_cluster, "cluster_sanity_check.pdf"),
       width = 12, height = 8)

### Save variables ###
save(clusters_calculated, file = paste0(wd_cluster, 'clusters_calculated.RData'))
write.csv(cluster_stations, paste0(wd_cluster, "STATION_CLUSTER_ID.csv"), row.names = FALSE)
write.csv(rHClustered, paste0(wd_cluster, "RATING_TRAINING_CLUSTERED.csv"), row.names = FALSE)


### DEVELOP REGRESSION -> REMEMBER: ONLY INSITU DATA ###

# Load variables
rHClustered <- read.csv(paste0(wd_cluster, "RATING_TRAINING_CLUSTERED.csv"))
rHClustered$site_no <- as.character(rHClustered$site_no)

# Add all band variables
rToReg <- setDT(rHClustered)[, ":="(
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

# Create 25% holdout set on a per-cluster basis
set.seed(23)
rToReg <- subset(rToReg, select=-c(B2_count))
holdout25 <- rToReg %>% group_by(cluster) %>% sample_frac(0.25) %>% ungroup() %>% setDT()
train75 <- anti_join(rToReg, holdout25, by = c("site_no", "sample_dt", "log10_SSC_mgL"))
rm(rToReg)

regBands <- dict()
regBands[["raw_bands"]] <- c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')
regBands[["raw_ratio"]] <- c(
    c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12'),
    colnames(train75)[unlist(lapply(colnames(train75), grepl, pattern=".B"))]
  )
regBands[["raw_square"]] <- c(
    c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12'),
    colnames(train75)[endsWith(colnames(train75), ".2")]
  )
regBands[["raw_sqrt"]] <- c(
    c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12'),
    colnames(train75)[endsWith(colnames(train75), "0.05")]
  )
regBands[["full_band"]] <- colnames(train75)[startsWith(colnames(train75), "B")]


relative_error <- function(true, pred){
  return(10^(
      median(
        abs(
          log10(
            pred/true
        )
        )
      )
    ) - 1)
}

# Create dataframe to store error metrics
relativeErrors <- data.frame(
  raw_bands = c(0, 0, 0, 0, 0, 0, 0),
  raw_ratio = c(0, 0, 0, 0, 0, 0, 0),
  raw_square = c(0, 0, 0, 0, 0, 0, 0),
  raw_sqrt = c(0, 0, 0, 0, 0, 0, 0),
  full_band = c(0, 0, 0, 0, 0, 0, 0)
)
rownames(relativeErrors) <- c(
  "cluster_1",
  "cluster_2",
  "cluster_3",
  "cluster_4",
  "cluster_5",
  "cluster_6",
  "net"
)

# Modified from https://github.com/evandethier/satellite-ssc/blob/master/landsat-calibration/landsat-57-calibration.R
# Try band combos
for(bands in c("raw_bands", "raw_ratio", "raw_square", "raw_sqrt", "full_band")){ 
# for(bands in c("raw_bands")) {

  # Directory to save model outputs
  bandDir <- paste0(wd_cluster, bands, "/")
  if(!dir.exists(bandDir)) {
    dir.create(bandDir)
  }

  # Create variable to save predicted values over each cluster
  errorDf <- data.frame(matrix(ncol = 2))
  colnames(errorDf) <- c("log10_SSC_mgL", "predlog10_SSC")

  # Relative Errors
  clusterErrors <- c()

  # Dataframe to store regression equations
  cluster_funs <- list()

  # Regress for each cluster
  for (i in c(1, 2, 3, 4, 5, 6)){
  # for (i in c(2)) {
    # regressors_sel <- regressors[-which(regressors == 'site_no')]
    
    lm_data <- train75[cluster == i] # only chooses sites within cluster
    # lm_data_hold <- lm_data_lm[-which(lm_data_lm$site_no %in% holdout_sts),] # for eliminating certain sites
    
    # Extract training data
    glm_y <- as.matrix(lm_data$log10_SSC_mgL)
    glm_x <- as.matrix(lm_data %>% dplyr::select(regBands[[bands]]))
    
    ssc_lm <- cv.glmnet(x = glm_x, y = glm_y, family = 'gaussian', type.measure = "mse", nfolds = 10)

    # Save regression equation for the cluster
    cluster_funs[[i]] <- ssc_lm
    
    # Select model within one standard error with min coefficients
    cv.opt <- coef(ssc_lm, s = "lambda.min")
    coef_ex <- cbind(rownames(cv.opt), as.numeric(cv.opt))
    colnames(coef_ex) <- c('variable', 'value')
    
    write.table(coef_ex, sep = ",", file = paste0(bandDir, "cluster_", i, "_lasso_fit_coeff.csv"), row.names = FALSE)
    
    glm_x <- NA
    glm_y <- NA

    # NEED TO STORE REGRESSOR
    # cluster_funs[[i]] <- ssc_lm

    forError <- holdout25[cluster == i]

    glm_pred <- predict(ssc_lm,
      newx = as.matrix(forError %>% dplyr::select(regBands[[bands]])),
      # s = "lambda.1se"
      s = "lambda.min"
    )

    holdoutPlot <- data.frame(
      log10_SSC_mgL = forError$log10_SSC_mgL,
      lambda.min = glm_pred
    )

    colnames(holdoutPlot) <- c("log10_SSC_mgL", "predlog10_SSC")

    # Remove NA values that can be introduced --> WHY
    holdoutPlot <- holdoutPlot[is.finite(rowSums(holdoutPlot)), ]

    # Combine data with all clusters
    errorDf <- rbind(errorDf, holdoutPlot)

    # Determine relative error
    RE <- relative_error(10^holdoutPlot$log10_SSC_mgL, 10^holdoutPlot$predlog10_SSC)
    clusterErrors <- append(clusterErrors, RE)

    regPlot <- ggplot(holdoutPlot, aes(x = 10^log10_SSC_mgL, y = 10^predlog10_SSC)) +
      geom_point(na.rm = TRUE) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      scale_x_log10(
        limits = c(1, 10000), 
        labels = fancy_scientific,
        breaks = c(10, 100, 1000, 10000)) +
      scale_y_log10(
        limits = c(1, 10000),
        labels = fancy_scientific,
        breaks = c(10, 100, 1000, 10000)) +
      # scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
      # season_facet +
      theme(legend.position = 'right') +
      theme_bw() +
      annotation_logticks() +
      labs(
        title = paste0("Prediction for Cluster ", i),
        caption = paste0("Band Combination: '", bands, "'"),
        subtitle = paste0("Relative Error: ", round(RE, 3)),
        x = 'Actual SSC (mg/L)',
        y = 'Satellite Estimated SSC (mg/L)'
      )
    
    # Save the regression plot
    ggsave(regPlot, filename = paste0(bandDir, "error-cluster-", i, ".pdf"),
       width = 8, height = 8)
  }
  
  # Remove first row of dataframe (caused by creating a matrix)
  errorDf <- errorDf[is.finite(rowSums(errorDf)), ]

  # Convert dataframe to numeric
  errorDf[1, ] <- as.numeric(errorDf[1, ])
  errorDf[2, ] <- as.numeric(errorDf[2, ])

  # Evaluate relative error over all clusters
  netErr <- relative_error(10^errorDf$log10_SSC_mgL, 10^errorDf$predlog10_SSC)
  clusterErrors <- append(clusterErrors, netErr)
  relativeErrors[[bands]] <- clusterErrors

  # Plot all cluster data
  totalPlot <- ggplot(errorDf, aes(x = 10^log10_SSC_mgL, y = 10^predlog10_SSC)) +
      geom_point(na.rm = TRUE) +
      stat_density_2d(
        aes(fill = ..level.., alpha = ..level..),
        bins = 10,
        geom = "polygon",
        colour = "black"
        ) +
      guides(alpha = "none", fill = "none") +
      scale_fill_gradient(low = "black", high = "#00ffd971") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      scale_x_log10(
        limits = c(1, 10000), 
        labels = fancy_scientific,
        breaks = c(10, 100, 1000, 10000)) +
      scale_y_log10(
        limits = c(1, 10000),
        labels = fancy_scientific,
        breaks = c(10, 100, 1000, 10000)) +
      # scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
      # season_facet +
      theme(legend.position = 'right') +
      theme_bw() +
      annotation_logticks() +
      labs(
        title = "Prediction for all clusters",
        caption = paste0("Band Combination: '", bands, "'"),
        subtitle = paste0("Relative Error: ", round(netErr, 3)),
        x = 'Actual SSC (mg/L)',
        y = 'Satellite Estimated SSC (mg/L)'
      )

  # Save the regression plot
  ggsave(totalPlot, filename = paste0(bandDir, "total-error.pdf"),
      width = 8, height = 8)

  # Save regression Equations
  save(cluster_funs, file = paste0(bandDir, "regressors.RData"))
}

# Save Relative Errors
write.table(relativeErrors, sep = ",", file = paste0(wd_cluster, "RELATIVE_ERRORS.csv"))

# Investigate relative errors and for each cluster, save the best model
# Cluster 1 --> full bands
# Cluster 2 --> raw_square
# Cluster 3 --> raw_ratio
# Cluster 4 --> raw_square
# Cluster 5 --> raw_ratio
# Cluster 6 --> full_band
cluster_regressors <- list()
load(paste0(wd_cluster, "full_band/regressors.RData"))
cluster_regressors[[1]] <- cluster_funs[[1]]
cluster_regressors[[6]] <- cluster_funs[[6]]
load(paste0(wd_cluster, "raw_square/regressors.RData"))
cluster_regressors[[2]] <- cluster_funs[[2]]
cluster_regressors[[4]] <- cluster_funs[[4]]
load(paste0(wd_cluster, "raw_ratio/regressors.RData"))
cluster_regressors[[3]] <- cluster_funs[[3]]
cluster_regressors[[5]] <- cluster_funs[[5]]
save(cluster_regressors, file=paste0(wd_cluster, "cluster_regressors.RData"))

# Save cluster bands along
cluster_bands <- list()
cluster_bands[[1]] <- regBands[["full_band"]]
cluster_bands[[6]] <- regBands[["full_band"]]
cluster_bands[[2]] <- regBands[["raw_square"]]
cluster_bands[[4]] <- regBands[["raw_square"]]
cluster_bands[[3]] <- regBands[["raw_ratio"]]
cluster_bands[[5]] <- regBands[["raw_ratio"]]
save(cluster_bands, file=paste0(wd_cluster, "cluster_bands.RData"))
# Create total error plot

# Create variable to save predicted values over each cluster
errorDf <- data.frame(matrix(ncol = 2))
colnames(errorDf) <- c("log10_SSC_mgL", "predlog10_SSC")
# Regress for each cluster
for (i in c(1, 2, 3, 4, 5, 6)){
  forError <- holdout25[cluster == i]

  glm_pred <- predict(cluster_regressors[[i]],
    newx = as.matrix(forError %>% dplyr::select(cluster_bands[[i]])),
    # s = "lambda.1se"
    s = "lambda.min"
  )

  holdoutPlot <- data.frame(
    log10_SSC_mgL = forError$log10_SSC_mgL,
    lambda.min = glm_pred
  )

  colnames(holdoutPlot) <- c("log10_SSC_mgL", "predlog10_SSC")

  # Remove NA values that can be introduced
  holdoutPlot <- holdoutPlot[is.finite(rowSums(holdoutPlot)), ]

  # Combine data with all clusters
  errorDf <- rbind(errorDf, holdoutPlot)
}

# Remove first row of dataframe (caused by creating a matrix)
errorDf <- errorDf[is.finite(rowSums(errorDf)), ]

# Convert dataframe to numeric
errorDf[1, ] <- as.numeric(errorDf[1, ])
errorDf[2, ] <- as.numeric(errorDf[2, ])

# Evaluate relative error over all clusters
netErr <- relative_error(10^errorDf$log10_SSC_mgL, 10^errorDf$predlog10_SSC)

# Plot all cluster data
totalPlot <- ggplot(errorDf, aes(x = 10^log10_SSC_mgL, y = 10^predlog10_SSC)) +
    geom_point(na.rm = TRUE, alpha=0.5) +
    stat_density_2d(
      aes(fill = ..level.., alpha = ..level..),
      bins = 10,
      geom = "polygon",
      colour = "black",
      alpha = 0.25
      ) +
    guides(alpha = "none", fill = "none") +
    scale_fill_gradient(low = "black", high = "#00ffd971") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_x_log10(
      limits = c(1, 10000), 
      labels = fancy_scientific,
      breaks = c(10, 100, 1000, 10000)) +
    scale_y_log10(
      limits = c(1, 10000),
      labels = fancy_scientific,
      breaks = c(10, 100, 1000, 10000)) +
    # scale_fill_brewer(palette = 'PuOr') + scale_color_brewer(palette = 'PuOr') +
    # season_facet +
    theme(legend.position = 'right') +
    theme_bw() +
    annotation_logticks() +
    labs(
      title = "Prediction for all clusters",
      subtitle = paste0("Relative Error: ", round(netErr, 3)),
      x = 'In-Situ SSC (mg/L)',
      y = 'Satellite Estimated SSC (mg/L)'
    )

# Save the regression plot
ggsave(totalPlot, filename = paste0(wd_cluster, "net-error.pdf"),
    width = 8, height = 8)
