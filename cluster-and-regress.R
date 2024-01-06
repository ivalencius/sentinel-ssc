### LIBRARY IMPORTS ###

# Trying alpha = 0.75

library(ggplot2)
library(data.table)
library(matrixStats)
library(glmnet)
library(hydroGOF)
library(dplyr)

### SET DIRECTORIES ###

print("SETTING UP IMPORTS")
# Set root directory
wd_root <- "/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc"
setwd(wd_root)

# Imports folder (store all import files here)
wd_imports <- paste0(wd_root, "/imports/")
# Exports folder (save all figures, tables here)
wd_exports <- paste0(wd_root, "/exports/")

wd_figures <- paste0(wd_exports, "figures/")
wd_cluster <- paste0(wd_exports, "cluster_regress/")
wd_gee <- paste0(wd_exports, "GEE_raw/")

export_folder_paths <- c(wd_figures, wd_cluster, wd_gee)
# , wd_exports_gc,wd_station_standalone,
# wd_standalone_models, wd_standalone_figures, wd_autocorrelation)
for (i in seq_along(export_folder_paths)) {
  path_sel <- export_folder_paths[i]
  if (!dir.exists(path_sel)) {
    dir.create(path_sel)
  }
}

### HElPER FUNCTION ###
# From https://github.com/evandethier/satellite-ssc/blob/master/landsat-calibration/landsat-57-calibration.R
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- log10(l)
  # return(parse(text=paste("'Discharge [m'", "^3* s", "^-1 ", "*']'", sep="")))
  return(parse(text = paste("10^", as.character(l), sep = "")))
}

### IMPORT IN SITU DATA ###
print("LOADING DATA")
# Load in-situ data
ratingSSC <- read.csv(file = paste0(wd_exports, "RATINGCURVE_SSC.csv")) # rating curve derived
inSituSSC <- read.csv(file = paste0(wd_exports, "ONEDAYSAMPLE_SSC.csv", sep = "")) # in-situ only
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
rSSC <- ratingSSC %>% mutate(
  SITE_NO = as.character(SITE_NO),
  SAMPLE_DT = as.Date(SAMPLE_DT))
iSSC <- inSituSSC %>% mutate(
  RATINGCURVE = FALSE,
  SITE_NO = as.character(SITE_NO),
  SAMPLE_DT = as.Date(SAMPLE_DT))
# Get full site numbers (with leading zeros)
rSSC$SITE_NO <- unlist(lapply(rSSC$SITE_NO, pad0))
iSSC$SITE_NO <- unlist(lapply(iSSC$SITE_NO, pad0))

# Change source for rating curve derived samples
rSSC$SOURCE <- ifelse(rSSC$R2 == 0, "USGS", "USGS RATING CURVE")
iSSC$SOURCE <- "USGS"

### IMPORT SENTINEL-2 DATA ###
# gee_files <- list.files(path = wd_gee, pattern = "sentinel*", full.names = TRUE)
# Columns of Sentinel-2 data export we want to remove
cols2Go <- c(
  "system.index",
  ".geo",
  "lat",
  "lon",
  # "snow_ice",
  "solar_az",
  "solar_zen"
  # "thin_cirrus_percentage",
  # "water_median",
  # "sensor",
  # "B2_count"
)

gee_raw <- read.csv(paste0(wd_gee, "REFLECTANCE-SSC-2017_2024_60m_drainage-2000km2.csv"))
stations <- read.csv(paste0(wd_gee, "VALID_SITES_width-60m_drainage-2000km2.csv")) %>%
select("MERIT_WIDTH_M", "GRWL_WIDTH_M","DRAINAGE_KM2", "STATION_NM", "HUC", 
"SITE_TYPE", "SITE_NO")
stations$SITE_NO <- as.character(stations$SITE_NO)
stations$SITE_NO <- unlist(lapply(stations$SITE_NO, pad0))
# gee_raw <- left_join(gee_raw, stations)
sentinel <- gee_raw[, !(names(gee_raw) %in% cols2Go)]
# No longer using band 1 or 9
setnames(sentinel,
    old = c(
      # "B1_median",
      "B2_median",
      "B3_median",
      "B4_median",
      "B5_median",
      "B6_median",
      "B7_median",
      "B8_median",
      "B8A_median",
      # "B9_median",
      "B11_median",
      "B12_median",
      "date",
      "product_id",
      "snow_ice",
      "thin_cirrus_percentage",
      "B2_count",
      "sensor",
      "water_median"
    ), new = c( "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12",
    "DATE","PRODUCT_ID", "SNOW_ICE", "THIN_CIRRUS_PERCENTAGE", "B2_COUNT",
    "SENSOR", "WATER_MEDIAN")
  )
sentinel$DATE <- as.Date(sentinel$DATE)
sentinel$SITE_NO <- as.character(sentinel$SITE_NO)
sentinel$SITE_NO <- unlist(lapply(sentinel$SITE_NO, pad0))
sentinel <- left_join(sentinel, stations, by = "SITE_NO")

# Load Landsat Data
landsat_ssc <- read.csv(paste0(wd_gee, "river_landsat_pred.csv"))
landsat_ssc <- landsat_ssc[, !(names(landsat_ssc) %in% c("month", "decade", "NUM_PIX"))]
landsat_ssc$SOURCE <- "LANDSAT"
landsat_ssc$SAMPLE_DT <- as.Date(landsat_ssc$SAMPLE_DT)
landsat_ssc$SITE_NO <- unlist(lapply(landsat_ssc$SITE_NO, pad0))
setnames(landsat_ssc, old = c("SSC_mgL","cluster"),
new = c("SSC_MGL", "LANDSAT_CLUSTER"))
landsat_ssc$LOG10_SSC_MGL <- log10(landsat_ssc$SSC_MGL)
landsat_ssc <- left_join(landsat_ssc, stations, by = "SITE_NO")
landsat_ssc2017 <- landsat_ssc %>% filter(SAMPLE_DT >= as.Date("2017-01-01"))

# Merge landsat and rating curve data
rSSC <- bind_rows(rSSC, landsat_ssc2017)

# Remove outlier stations --> too much ice or strange artificats
badSites <- c(
  "WHITE R NEAR OACOMA,SD",
  "BRAZOS RV NR HEMPSTEAD, TX",
  "BRAZOS RV AT RICHMOND, TX",
  "COLORADO RV AT WHARTON, TX",
  "ALSEK R AT DRY BAY NR YAKUTAT AK",
  "TAKU R NR JUNEAU AK"
)
sentinel <- sentinel[!(sentinel$STATION_NM %in% badSites),]
rSSC <- rSSC[!(rSSC$STATION_NM %in% badSites),]
iSSC <- iSSC[!(iSSC$STATION_NM %in% badSites),]

### Harmonize USGS AND SENTINEL-2 DATA ###
print("HARMONIZING DATA")
lag_days <- c(4, 3)

for (lag_day in lag_days) {
  # Create folder for specific lag times
  lag_dir <- paste0(wd_cluster, "lag_", lag_day, "/")
  if (!dir.exists(lag_dir)) {
    dir.create(lag_dir)
  }

  # Harmonize Rating Curve Derived Samples
  # Create lead lag times and add to DF
  sentinel_copy <- setDT(copy(sentinel))[
    , ":="(
      MATCH_DT_START = DATE - lag_day,
      MATCH_DT_END = DATE + lag_day)
  ]
  # Match by dates inside lead-lag range
  rHarmonized <- sentinel_copy[setDT(copy(rSSC))[, ":="(MATCH_DT = SAMPLE_DT)],
    on = .(
      SITE_NO == SITE_NO,
      MATCH_DT_START <= MATCH_DT,
      MATCH_DT_END >= MATCH_DT
    )
  ][
    !is.na(B2)
  ][
    , LAG_DAYS := as.numeric(difftime(SAMPLE_DT, DATE), "days")
  ][
    , ":="(MATCH_DT_START = NULL, MATCH_DT_END = NULL)
  ]

  # Add landsat dates

  iHarmonized <- sentinel_copy[setDT(copy(iSSC))[, ":="(MATCH_DT = SAMPLE_DT)],
    on = .(
      SITE_NO == SITE_NO,
      MATCH_DT_START <= MATCH_DT,
      MATCH_DT_END >= MATCH_DT
    )
  ][
    !is.na(B2)
  ][
    , LAG_DAYS := as.numeric(difftime(SAMPLE_DT, DATE), "days")
  ][
    , ":="(MATCH_DT_START = NULL, MATCH_DT_END = NULL)
  ]

  # ISSUE: in-situ samples not stored for sites with no linear regression
  missingSites <- unique(iHarmonized$SITE_NO)[!(unique(iHarmonized$SITE_NO) %in% unique(rHarmonized$SITE_NO))]
  missingData <- iHarmonized %>%
    filter(SITE_NO %in% missingSites)
  rHMerge <- rbind(rHarmonized, missingData, fill=TRUE)

  # Reduce harmonized samples (rating curve derived SSC)
  rHClean <- rHMerge %>%
    # For each image --> isolate one with minimum lag day
    group_by(SITE_NO, PRODUCT_ID) %>%
    slice(which.min(abs(LAG_DAYS))) %>%
    ungroup() %>%
    # Now there can be multiple images per day --> reduce to one
    group_by(SITE_NO, SAMPLE_DT) %>%
    summarise(across(where(is.numeric), mean)) %>%
    ungroup()

  # Reduce harmonized samples (in-situ derived SSC)
  iHClean <- iHarmonized %>%
    # For each image --> isolate one with minimum lag day
    group_by(SITE_NO, PRODUCT_ID) %>%
    slice(which.min(abs(LAG_DAYS))) %>%
    ungroup() %>%
    # Now there can be multiple images per day --> reduce to one
    group_by(SITE_NO, SAMPLE_DT) %>%
    summarise(across(where(is.numeric), mean)) %>%
    ungroup()

  # Averaging across samples you lose:
  # "DATE"        "PRODUCT_ID"  "SENSOR" "STATION_NM"  "SITE_TYPE"  
  # "AGENCY_CD"   "RATINGCURVE"
  rHClean <- left_join(rHClean, rHMerge %>% 
  select(SITE_NO, STATION_NM, SAMPLE_DT, PRODUCT_ID, SITE_TYPE, RATINGCURVE, SOURCE) 
  %>% distinct(), by = c("SITE_NO", "SAMPLE_DT"))
  iHClean <- left_join(iHClean, iHarmonized %>% 
  select(SITE_NO, STATION_NM, SAMPLE_DT, PRODUCT_ID, SITE_TYPE, RATINGCURVE, SOURCE) 
  %>% distinct(), by = c("SITE_NO", "SAMPLE_DT"))

  # Clean up fractional lag days
  rHClean$LAG_DAYS <- floor(rHClean$LAG_DAYS)
  iHClean$LAG_DAYS <- floor(iHClean$LAG_DAYS)

  # Remove excess snow
  rHClean <- rHClean %>% filter(SNOW_ICE <= 5)
  iHClean <- iHClean %>% filter(SNOW_ICE <= 5)

  # Save data
  write.csv(rHClean,
    file = paste0(lag_dir, "lag_", lag_day, "_RATING_TRAINING_SAMPLES.csv"),
    row.names = FALSE
  )
  write.csv(iHClean,
    file = paste0(lag_dir, "lag_", lag_day, "_INSITU_TRAINING_SAMPLES.csv"),
    row.names = FALSE
  )
}

# Now just cluster based on landsat clusters
# site_clusters <- rSSC %>% filter(SOURCE == "LANDSAT") %>% 
# select(SITE_NO, LANDSAT_CLUSTER) %>% group_by(SITE_NO) %>% 
# summarise(CLUSTER = unique(LANDSAT_CLUSTER)) %>% ungroup()
# switch <- c("RATING", "INSITU")
# for (lag_day in lag_days) {
#   lag_dir <- paste0(wd_cluster, "lag_", lag_day, "/")
#   for (type in switch) {
#     data <- read.csv(paste0(lag_dir, "lag_", lag_day, "_", type, "_TRAINING_SAMPLES.csv"))
#     data$SITE_NO <- as.character(data$SITE_NO)
#     dataClustered <- left_join(data, site_clusters, by = "SITE_NO")
#     # Create 25% holdout set on a per-cluster basis
#     set.seed(23)
#     # rToReg <- subset(data, select = -c(B2_count))
#     holdout25 <- dataClustered %>%
#       group_by(CLUSTER) %>%
#       sample_frac(0.25) %>%
#       ungroup() %>%
#       setDT()
#     train75 <- anti_join(dataClustered, holdout25, by = c("SITE_NO", "SAMPLE_DT", "LOG10_SSC_MGL"))
#     holdout25$TYPE <- "HOLDOUT"
#     train75$TYPE <- "TRAINING"
#     export <- rbind(holdout25, train75)
#     # Save variables ### --> SUBSET TRAINING SET
#     write.csv(export, paste0(lag_dir, "lag_", lag_day, "_", type, "_CLUSTERED.csv"), row.names = FALSE)
#   }
# }

extract_cluster_data <- function(cluster_regressors, data) {
  cluster_medians <- data.table(matrix(nrow = 0, ncol = length(cluster_regressors) + 1))
  names(cluster_medians) <- c("SITE_NO", cluster_regressors)
  unique_sites <- unique(as.character(data$SITE_NO))
  pb <- txtProgressBar(0, length(unique_sites), style = 3)
  for (i in seq_along(unique_sites)) {
    setTxtProgressBar(pb, i)
    site <- unique_sites[i]
    # Extract data
    site_data <- as.data.table(
      data[data$SITE_NO == site, ]
    )[, ..cluster_regressors] %>% as.matrix()
    # Reduce by median
    site_median <- colMedians(site_data, na.rm = TRUE) # avoid NA in any column
    # Create table to hold site data
    site_df <- data.table(t(data.table(site_median)))
    names(site_df) <- cluster_regressors
    site_df$SITE_NO <- site
    # Generate regression
    cluster_medians <- rbind(cluster_medians, site_df)
  }
  return(cluster_medians)
}

# Cluster data
print("Clustering data")
switch <- c("RATING", "INSITU")
for (lag_day in lag_days) {
  lag_dir <- paste0(wd_cluster, "lag_", lag_day, "/")
  # Bands used for regression --> no more B1 or B9
  cluster_regressors <- c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12")

  for (type in switch) {
    data <- read.csv(paste0(lag_dir, "lag_", lag_day, "_", type, "_TRAINING_SAMPLES.csv"))
    data$SITE_NO <- as.character(data$SITE_NO)
    # Extract median spectral information at each site
    data_median <- extract_cluster_data(cluster_regressors, data)
    # 6 clusters --> following Dethier et. al (2020)
    clusters_calculated <- kmeans(data_median[, ..cluster_regressors],
      centers = 6,
      nstart = 100, iter.max = 50
    )
    # Compute cluster centers
    cluster_centers <- clusters_calculated$centers
    # Assign cluster to each site
    data_median$CLUSTER <- clusters_calculated$cluster
    # Extract cluster for each station
    cluster_stations <- data_median[, .(SITE_NO, CLUSTER)]
    dataClustered <- left_join(data, cluster_stations,
      by = c("SITE_NO" = "SITE_NO")
    )
    # Create 25% holdout set on a per-cluster basis
    set.seed(23)
    # rToReg <- subset(data, select = -c(B2_count))
    holdout25 <- dataClustered %>%
      group_by(CLUSTER) %>%
      sample_frac(0.25) %>%
      ungroup() %>%
      setDT()
    train75 <- anti_join(dataClustered, holdout25)
    holdout25$TYPE <- "HOLDOUT"
    train75$TYPE <- "TRAINING"
    export <- rbind(holdout25, train75)
    ### Save variables ### --> SUBSET TRAINING SET
    save(clusters_calculated, file = paste0(lag_dir, "lag_", lag_day, "_", type, "_clusters_calculated.RData"))
    # write.csv(cluster_stations, paste0(lag_dir, "STATION_CLUSTER_ID.csv"), row.names = FALSE)
    write.csv(export, paste0(lag_dir, "lag_", lag_day, "_", type, "_CLUSTERED.csv"), row.names = FALSE)
  }
}

# Quick cluster sanity check
clust_sanity_check <- function(lag_oi, type_oi) {
  test_df <- read.csv(paste0(wd_cluster, "lag_", lag_oi, "/lag_", lag_oi, "_", type_oi, "_CLUSTERED.csv"))
  # Visualize Clusters
  # ssc_categories <- c(0, 50, 100, 250, 500, 750, 1e6)
  ssc_categories <- c(0,50,100,200,500,1e6)
  # ssc_categories <- c(0,10,25,50,75,100,150,200,250,300,350, 400, 450, 500,600, 700, 800,900,1000,1100,1500, 1e6)

  # Generate SSC labels as 'low value' - 'high value'
  ssc_category_labels <- paste0(ssc_categories[-length(ssc_categories)], "-", c(ssc_categories[-1]))
  # Make highest SSC category "> highest value"
  ssc_category_labels[length(ssc_category_labels)] <- paste0("> ", ssc_categories[length(ssc_category_labels)])

  viz_df <- as.data.table(test_df)[
    , ":="(cluster_sel = CLUSTER,
      # # Categorize SSC value as one of selected categories
      ssc_category = cut(10^LOG10_SSC_MGL,
        breaks = ssc_categories,
        labels = ssc_category_labels
    ))
  ][]

  # Generate median B,G,R for each SSC category and each cluster or site
  ssc_category_color <- viz_df %>%
    group_by(cluster_sel, ssc_category) %>%
    summarise_at(vars("B4", "B3", "B2"), mean)

  raster_color_types <- geom_raster(aes(fill = rgb(B4 / 2600, B3 / 2600, B2 / 2600))) # false color?

  ggplot(ssc_category_color, aes(x = cluster_sel, y = ssc_category)) +
    raster_color_types +
    scale_fill_identity() +
    theme_classic() +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    # theme(axis.text.x = element_text(angle = 90)) +
    labs(
      y = "SSC range (mg/L)",
      x = "River grouping"
    ) +
    theme(text = element_text(family = "JetBrains Mono NL"))
}

# Generate cluster sanity check figures
for (lag_day in lag_days) {
  for (type in switch) {
    fname <- paste0(wd_cluster, "lag_", lag_day, "/lag_", lag_day, "_", type, "_cluster_sanity_check.png")
    fig <- clust_sanity_check(lag_day, type)
    ggsave(fig, height = 5, width = 5, dpi = 300, filename = fname)
  }
}

# Investigate snow and ice across clusters
# f <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/exports/cluster_regress/lag_4/lag_4_INSITU_CLUSTERED.csv")
# ggplot(f, aes(x=SNOW_ICE)) +
# geom_histogram() + 
# facet_wrap(~CLUSTER)

# Function to calculate relative error
relative_error <- function(true, pred) {
  return(10^(median(abs(log10(pred / true)), na.rm=TRUE)) - 1)
}

# Function to generate lasso regression
gen_reg <- function(lag_day, type) {

  lag_dir <- paste0(wd_cluster, "lag_", lag_day, "/")
  data <- read.csv(paste0(lag_dir, "lag_", lag_day, "_", type, "_CLUSTERED.csv"))

  data$SITE_NO <- as.character(data$SITE_NO)

  # Add all band variables
  dataPrepped <- setDT(data)[, ":="(
  # Add squared columns --> No more B1 or B9
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
  # B2.B1 = B2 / B1,
  # B3.B1 = B3 / B1,
  # B4.B1 = B4 / B1,
  # B5.B1 = B5 / B1,
  # B6.B1 = B6 / B1,
  # B7.B1 = B7 / B1,
  # B8.B1 = B8 / B1,
  # B8A.B1 = B8A / B1,
  # B9.B1 = B9 / B1,
  # B11.B1 = B11 / B1,
  # B12.B1 = B12 / B1,

  B3.B2 = B3 / B2,
  B4.B2 = B4 / B2,
  B5.B2 = B5 / B2,
  B6.B2 = B6 / B2,
  B7.B2 = B7 / B2,
  B8.B2 = B8 / B2,
  B8A.B2 = B8A / B2,
  # B9.B2 = B9 / B2,
  B11.B2 = B11 / B2,
  B12.B2 = B12 / B2,

  B4.B3 = B4 / B3,
  B5.B3 = B5 / B3,
  B6.B3 = B6 / B3,
  B7.B3 = B7 / B3,
  B8.B3 = B8 / B3,
  B8A.B3 = B8A / B3,
  # B9.B3 = B9 / B3,
  B11.B3 = B11 / B3,
  B12.B3 = B12 / B3,

  B5.B4 = B5 / B4,
  B6.B4 = B6 / B4,
  B7.B4 = B7 / B4,
  B8.B4 = B8 / B4,
  B8A.B4 = B8A / B4,
  # B9.B4 = B9 / B4,
  B11.B4 = B11 / B4,
  B12.B4 = B12 / B4,

  B6.B5 = B6 / B5,
  B7.B5 = B7 / B5,
  B8.B5 = B8 / B5,
  B8A.B5 = B8A / B5,
  # B9.B5 = B9 / B5,
  B11.B5 = B11 / B5,
  B12.B5 = B12 / B5,

  B7.B6 = B7 / B6,
  B8.B6 = B8 / B6,
  B8A.B6 = B8A / B6,
  # B9.B6 = B9 / B6,
  B11.B6 = B11 / B6,
  B12.B6 = B12 / B6,

  B8.B7 = B8 / B7,
  B8A.B7 = B8A / B7,
  # B9.B7 = B9 / B7,
  B11.B7 = B11 / B7,
  B12.B7 = B12 / B7,

  B8A.B8 = B8A / B8,
  # B9.B8 = B9 / B8,
  B11.B8 = B11 / B8,
  B12.B8 = B12 / B8,

  # B9.B8A = B9 / B8A,
  B11.B8A = B11 / B8A,
  B12.B8A = B12 / B8A,

  # B11.B9 = B11 / B9,
  # B12.B9 = B12 / B9,

  B12.B11 = B12 / B11
  )]

  # If any band has a 0 value, it will create infinite values in the data, need to remove
  # Only check complete cases over the bands and SSC
  dataFinite <- dataPrepped[complete.cases(dataPrepped %>% 
  select(starts_with("B") | LOG10_SSC_MGL))]
  

  # Extract holdout and training data
  train75 <- dataFinite %>% filter(TYPE == "TRAINING")
  holdout25 <- dataFinite %>% filter(TYPE == "HOLDOUT")

  # Add columns to store holdout predictions
  holdout25$PRED_LOG10_SSC_MGL <- 0

  reg_dir <- paste0(lag_dir, type, "_regression/")
  if (!dir.exists(reg_dir)) {
    dir.create(reg_dir)
  }

  errors <- data.frame(
    CLUSTER = c(1:6, "NET"),
    RMSE = rep(0, 7),
    MAE = rep(0, 7),
    R2 = rep(0, 7),
    RE = rep(0, 7)
  )

  regVars <- data.frame(
    Variable = c("Intercept", colnames(train75 %>% select(starts_with("B") 
    & -starts_with("B2_COUNT")))))

  # Generate regression for all clusters
  for (clust in 1:6) {
    print(paste0("Cluster: ",clust))
    lm_data <- train75 %>% filter(CLUSTER == clust)
    holdout_data <- holdout25 %>% filter(CLUSTER == clust)
    # Extract training data
    glm_y <- as.matrix(lm_data$LOG10_SSC_MGL)
    glm_x <- as.matrix(lm_data %>% select(starts_with("B") & -starts_with("B2_COUNT")))

    ssc_lm <- cv.glmnet(x = glm_x, y = glm_y, family = "gaussian", alpha = 0.5,
    type.measure = "mse", nfolds = 10, nlambda = 100)
    
    save(ssc_lm, file = paste0(lag_dir, type, "_regression/ssc_lm_", clust, ".RData"))

    # Save plot of lambda
    png(filename=paste0(reg_dir, "cluster_", clust, "_lambda_plot.png"),
    width=6, height=6, units="in", res=300)
    plot(ssc_lm)
    dev.off()

    # Select model within one standard error with min coefficients --> No longer do this
    cv.opt <- coef(ssc_lm, s = "lambda.min")
    # cv.opt <- coef(ssc_lm, s = "lambda.1se")
    regVars[, paste0("CLUSTER_", clust)] <- as.numeric(cv.opt)

    # Predict on holdout data
    glm_pred <- predict(ssc_lm,
    newx = as.matrix(holdout_data %>% select(starts_with("B") & -starts_with("B2_COUNT"))),
    # s = "lambda.1se"
    s = "lambda.min"
    )

    holdout_data$PRED_LOG10_SSC_MGL <- as.numeric(glm_pred)
    # Set Inf values to NA
    holdout_data$PRED_LOG10_SSC_MGL[holdout_data$PRED_LOG10_SSC_MGL == Inf] <- NA

    # Calculate errors
    RMSE <- rmse(10^holdout_data$PRED_LOG10_SSC_MGL, 10^holdout_data$LOG10_SSC_MGL, na.rm = TRUE)
    MAE <- mae(10^holdout_data$PRED_LOG10_SSC_MGL, 10^holdout_data$LOG10_SSC_MGL, na.rm = TRUE)
    R2 <- cor(10^holdout_data$LOG10_SSC_MGL, 10^holdout_data$PRED_LOG10_SSC_MGL, use="complete.obs")^2
    RE <- relative_error(10^holdout_data$LOG10_SSC_MGL,
    10^holdout_data$PRED_LOG10_SSC_MGL)

    # Add to errors data frame
    errors[errors$CLUSTER == clust, "RMSE"] <- RMSE
    errors[errors$CLUSTER == clust, "MAE"] <- MAE
    errors[errors$CLUSTER == clust, "R2"] <- R2
    errors[errors$CLUSTER == clust, "RE"] <- RE
    
    tryCatch({
    ggplot(holdout_data, aes(x = 10^LOG10_SSC_MGL, y = 10^PRED_LOG10_SSC_MGL)) +
      geom_point(aes(shape = SOURCE), na.rm = TRUE, alpha=0.7) +
      stat_density_2d(
        aes(fill = ..level.., alpha = ..level..),
        bins = 10,
        geom = "polygon",
        colour = "black",
        alpha = 0.25,
        na.rm = TRUE
        ) +
      guides(alpha = "none", fill = "none") +
      scale_fill_gradient(low = "black", high = "#00ffd971") +
      scale_shape_manual(values=c(15, 16, 17)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      scale_x_log10(
        limits = c(1, 10000), 
        labels = fancy_scientific,
        breaks = c(10, 100, 1000, 10000)) +
      scale_y_log10(
        limits = c(1, 10000),
        labels = fancy_scientific,
        breaks = c(10, 100, 1000, 10000)) +
      theme_bw() +
      theme(legend.position = c(0.8, 0.2)) +
      annotation_logticks() +
      labs(
        shape = "Data Source",
        title = paste0("Maximum Lead/Lag Time: ±", lag_day, " days | Source: ", type, " | Cluster ", clust),
        subtitle = paste0("Relative Error: ", round(RE, 4)),
        x = 'In-Situ SSC (mg/L)',
        y = 'Satellite Estimated SSC (mg/L)') +
    theme(text = element_text(family = "JetBrains Mono NL"))
    ggsave(
      filename = paste0(reg_dir, "cluster_", clust, "_regression_plot.png"),
      width = 10, height = 10, units = "in", dpi = 300)
      }, error = function(e) {
      print(paste0("Error in cluster: ", clust, " -> most likely intercept only variable"))},
      warning = function(w) {
        ggplot(holdout_data, aes(x = 10^LOG10_SSC_MGL, y = 10^PRED_LOG10_SSC_MGL)) +
        geom_point(aes(shape = SOURCE), na.rm = TRUE, alpha=0.7) +
        stat_density_2d(
          aes(fill = ..level.., alpha = ..level..),
          bins = 10,
          geom = "polygon",
          colour = "black",
          alpha = 0.25,
          na.rm = TRUE
          ) +
        guides(alpha = "none", fill = "none") +
        scale_fill_gradient(low = "black", high = "#00ffd971") +
        scale_shape_manual(values=c(15, 16, 17)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        scale_x_log10(
          limits = c(1, 10000), 
          labels = fancy_scientific,
          breaks = c(10, 100, 1000, 10000)) +
        scale_y_log10(
          limits = c(1, 10000),
          labels = fancy_scientific,
          breaks = c(10, 100, 1000, 10000)) +
        theme_bw() +
        theme(legend.position = c(0.8, 0.2)) +
        annotation_logticks() +
        labs(
          shape = "Data Source",
          title = paste0("Maximum Lead/Lag Time: ±", lag_day, " days | Source: ", type, " | Cluster ", clust),
          subtitle = paste0("Relative Error: ", round(RE, 4)),
          x = 'In-Situ SSC (mg/L)',
          y = 'Satellite Estimated SSC (mg/L)') +
      theme(text = element_text(family = "JetBrains Mono NL"))
      ggsave(
        filename = paste0(reg_dir, "cluster_", clust, "_regression_plot.png"),
        width = 10, height = 10, units = "in", dpi = 300)
    })

    # Add predicted data to holdoutData
    holdout25[CLUSTER == clust]$PRED_LOG10_SSC_MGL <- holdout_data$PRED_LOG10_SSC_MGL
  }

  write.table(regVars, file = paste0(reg_dir, "regression_variables.csv"), sep = ",", row.names = FALSE)

  # Need to remove any infinite values (maybe 1)
  holdout25 <- holdout25[complete.cases(holdout25 %>% select(starts_with("B") & -starts_with("B2_COUNT"))),]

  # Calculate errors
  RMSE <- rmse(10^holdout25$PRED_LOG10_SSC_MGL, 10^holdout25$LOG10_SSC_MGL, na.rm = TRUE)
  MAE <- mae(10^holdout25$PRED_LOG10_SSC_MGL, 10^holdout25$LOG10_SSC_MGL, na.rm = TRUE)
  R2 <- cor(10^holdout25$LOG10_SSC_MGL, 10^holdout25$PRED_LOG10_SSC_MGL, use="complete.obs")^2
  RE <- relative_error(10^holdout25$LOG10_SSC_MGL,
  10^holdout25$PRED_LOG10_SSC_MGL)

  # Add to errors data frame
  errors[errors$CLUSTER == "NET", "RMSE"] <- RMSE
  errors[errors$CLUSTER == "NET", "MAE"] <- MAE
  errors[errors$CLUSTER == "NET", "R2"] <- R2
  errors[errors$CLUSTER == "NET", "RE"] <- RE

  write.table(errors, file = paste0(reg_dir, "regression_errors.csv"), sep = ",", row.names = FALSE)
 
  write.csv(holdout25, file = paste0(reg_dir, "holdout_data.csv"), row.names = FALSE)

  ggplot(holdout25, aes(x = 10^LOG10_SSC_MGL, y = 10^PRED_LOG10_SSC_MGL)) +
    geom_point(aes(shape = SOURCE), na.rm = TRUE, alpha=0.5) +
    stat_density_2d(
      aes(fill = ..level.., alpha = ..level..),
      bins = 10,
      geom = "polygon",
      colour = "black",
      alpha = 0.25,
      na.rm = TRUE
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
    scale_shape_manual(values=c(15, 16, 17)) +
    theme_bw() +
    theme(legend.position = c(0.8, 0.2)) +
    annotation_logticks() +
    labs(
      shape = "Data Source",
      title = paste0("Maximum Lead/Lag Time: ±", lag_day, " days | Source: ", type),
      subtitle = paste0("Relative Error: ", round(RE, 4)),
      x = 'In-Situ SSC (mg/L)',
      y = 'Satellite Estimated SSC (mg/L)') +
  theme(text = element_text(family = "JetBrains Mono NL"))

  ggsave(filename = paste0(reg_dir, "All_cluster_regression_plot.png"),
    width = 10, height = 10, units = "in", dpi = 300)

  # Save individual station calibration regression plots
  station_dir <- paste0(reg_dir, "station_calibration_plots/")
  if (!dir.exists(station_dir)) {
    dir.create(station_dir)
  }

  for (site in unique(holdout25$SITE_NO)) {
    data <- holdout25 %>% filter(SITE_NO == site)
    data$LAG_DAYS <- abs(data$LAG_DAYS)

    RE <- relative_error(10^data$LOG10_SSC_MGL,
    10^data$PRED_LOG10_SSC_MGL)
  
    ggplot(data, aes(x = 10^LOG10_SSC_MGL, y = 10^PRED_LOG10_SSC_MGL)) +
    geom_point(aes(color = factor(LAG_DAYS), shape = SOURCE), na.rm = TRUE, size = 3, alpha = 0.8) +
    scale_color_brewer(limits = factor(c(0,1,2,3,4)), palette = "Purples", direction = -1) +
    scale_shape_manual(values=c(15, 16, 17)) +
    geom_abline(slope = 1, intercept = 0, color = "black") +
    scale_x_log10(
      limits = c(1, 10000), 
      labels = fancy_scientific,
      breaks = c(10, 100, 1000, 10000)) +
    scale_y_log10(
      limits = c(1, 10000),
      labels = fancy_scientific,
      breaks = c(10, 100, 1000, 10000)) +
    theme_bw() +
    annotation_logticks() +
    theme(legend.position = c(0.85, 0.3)) +
    labs(
      shape = "Data Source",
      color = "|Lead/Lag Days|",
      title = paste0("Maximum Lead/Lag Time: ±", lag_day, " days | Source: ", type),
      subtitle = paste0("Relative Error: ", round(RE, 4)),
      x = 'In-Situ SSC (mg/L)',
      y = 'Satellite Estimated SSC (mg/L)') +
    theme(text = element_text(family = "JetBrains Mono NL"))
  
    ggsave(filename = paste0(station_dir, "indiv_calib_", site, ".png"),
    width = 10, height = 10, units = "in", dpi = 300)
  }
}

# Generate regressions
print("GENERATING REGRESSIONS")
# Not enough samples in each cluster
for (lag_day in lag_days) {
  # for (lag_day in c(3,4)) {
  for (type in switch) {
    print(paste0("Lag_day: ", lag_day, " | Type: ", type))
    gen_reg(lag_day, type)
  }
}

# Examine all net errors
switch <- c("RATING", "INSITU")
lag_days <- c(3, 4)
allErrors <- data.frame(
  type = rep(switch, each = length(lag_days)),
  lag_day = rep(lag_days, length(switch)),
  net_Error = rep(0, length(lag_days) * length(switch))
) %>% arrange(lag_day)
for (lag_day in lag_days) {
  for (type in switch) {
    errorDf <- read.csv(paste0(wd_cluster, "lag_", lag_day, "/", type, "_regression/regression_errors.csv"))
    allErrors[allErrors$type == type & allErrors$lag_day == lag_day, "net_Error"] <- errorDf[errorDf$CLUSTER == "NET", "RE"]
  }
}
allErrors <- allErrors %>% arrange(net_Error)
View(allErrors)
write.csv(file=paste0(wd_cluster, "Errors.csv"), allErrors, row.names = FALSE)