################################ TO DO #########################################
# - Histogram of SSC samples
# - Histogram of SSC samples *after* rating curve
# - Some fluxes are negative --> take absolute value?

#### LIBRARY IMPORTS ####
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

#### SET DIRECTORIES ####

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
# Sub folders
# wd_extent<- paste0(wd_exports, 'station_transects/')
wd_gee <- paste0(wd_exports, "GEE_raw/")
wd_rating <- paste0(wd_exports, "rating_curves/")

# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_exports, wd_figures, wd_gee, wd_rating)
                         # , wd_exports_gc,wd_station_standalone, 
                         # wd_standalone_models, wd_standalone_figures, wd_autocorrelation)
for(i in seq_along(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}

projection <- CRS("+proj=longlat +datum=WGS84 +no_defs")

#### INITIALIZE MAP DATA FOR N.AMERICA ####

cat("\t", "-> Initializing USA map","\n")
setwd(wd_imports)
# get states shapefile for clipping/display
us_states_spatial <-  us_states(map_date = NULL, resolution = c("low"), states = NULL) %>% as('Spatial')
# us_states_spatial <- us_states_spatial[us_states_spatial$state_abbr != "AK" & us_states_spatial$state_abbr != "HI" & us_states_spatial$state_abbr != "PR",]
# us_states_spatial <- us_states_spatial[us_states_spatial$state_abbr != "AK",]

#plot(us_states_spatial)

# import canadian provinces
# cat('\t','-> Initializing CA map','\n')
# canada_prov<-getData('GADM', country="CAN", level=1) # provinces
# canada_prov <- spTransform(canada_prov, projection)
# proj4string(canada_prov) <- proj4string(us_states_spatial)

# canada_prov <- canada_prov[,c('NAME_1')] # only select column with province name
# names(canada_prov) <- c('name') # rename columns labeling canadian provinces
us_states_spatial <- us_states_spatial[,c('state_abbr')] # only select column with province name
# us_states_spatial <- rbind(us_states_spatial,canada_prov) # combine canadian and us shapefiles

# fortify state shapefile
us_ca <- fortify(us_states_spatial)

#### IMPORT IN SITU SSC DATA ####

# From USGS
print("IMPORTING IN SITU DATA")
startDate <- "2017-01-01"
endDate <- "2023-02-07"
cat("\t","-> Start Date:", startDate, "\n")
cat("\t","-> End Date:", endDate, "\n")

#ssc_codes <- c('00530','80154')
ssc_codes <- c("80154")
cat("\t","-> Downloading USGS SSC Data", "\n")
wqp_insitu_raw <- data.table(readWQPdata(parameterCd = ssc_codes,
                                          startDate=startDate,
                                          endDate=endDate))[,":=" (
                                            agency_cd = "USGS",
                                            site_no = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
                                            sample_dt = ActivityStartDate,
                                            SSC_mgL = ResultMeasureValue,
                                            sample_depth_m = ActivityDepthHeightMeasure.MeasureValue * 0.3048,
                                            sample_method = SampleCollectionMethod.MethodName
                                            )][,.(agency_cd, site_no, sample_dt, SSC_mgL, sample_depth_m, sample_method)]

# Need to use WQP data as well
nwis_insitu_raw <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("agency_cd", "site_no", "sample_dt", "SSC_mgL"))))
for (state in us_states_spatial$state_abbr){
  nwis_data <- data.table(readNWISdata(stateCd=state,
                                     parameterCd = ssc_codes, 
                                            startDt=startDate,
                                            endDt=endDate))
  if(nrow(nwis_data)==0){} # Do nothing if no data returned
  else{
    nwis_data_clean <- nwis_data[, ":=" (agency_cd = agency_cd,
    site_no = site_no,
    sample_dt = dateTime,
    SSC_mgL = X_80154_00003
    )][,.(agency_cd, site_no, sample_dt, SSC_mgL)]
    nwis_insitu_raw <- rbind(nwis_insitu_raw, nwis_data_clean)
  }
}
nwis_insitu_raw <- nwis_insitu_raw %>% filter(SSC_mgL >= 0)

# Merge nwis and usgs
usgs_insitu_raw <- bind_rows(data.frame(wqp_insitu_raw), data.frame(nwis_insitu_raw))
# Filter out nonsense data
usgs_insitu_raw <- usgs_insitu_raw %>% filter(SSC_mgL >=0)
# Convert to all Caps
usgs_insitu_raw <- mutate_all(usgs_insitu_raw, toupper)

# For WQP need USGS- prefix
usgs_stations <- paste0("USGS-", unique(usgs_insitu_raw$site_no))

# Get properties of usgs WQP stations
cat("\t","-> Downloading USGS Station Metadata", "\n")
station_info <- data.table(whatWQPsites(siteid = usgs_stations))[,":="(
  station_nm = MonitoringLocationName,
  drainage_area_km2 = DrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
  contributing_drainage_area_km2 = ContributingDrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
  lat = LatitudeMeasure,
  lon = LongitudeMeasure,
  elevation_m = VerticalMeasure.MeasureValue *0.3048, # ft -> m # nolint
  site_no = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
  site_type = MonitoringLocationTypeName,
  HUC = HUCEightDigitCode
)][, .(site_no,
      station_nm,
      lat,
      lon,
      drainage_area_km2,
      contributing_drainage_area_km2,
      elevation_m,
      site_type,
      HUC
      )]
# HUC8 code for use with https://developers.google.com/earth-engine/datasets/catalog/USGS_WBD_2017_HUC08#table-schema

#is_w_station <- in-situ data harmonized with station data
is_w_stations <- left_join(usgs_insitu_raw, station_info, by=("site_no"="site_no"))

# Remove station data with less than 10 occurances --> for rating curve
# usgs_insitu_raw <- usgs_insitu_raw[usgs_insitu_raw$site_no %in% names(which(table(usgs_insitu_raw$site_no) >= 10)), ]

# All data stored as characters, convert some columns to numbers'
cols <- c("lat", "lon", "drainage_area_km2", "SSC_mgL", "sample_depth_m")
# cols2 <- c()
# station_info[, cols] <- lapply(cols, function(x) as.numeric(station_info[[x]]))
is_w_stations[, cols] <- lapply(cols, function(x) as.numeric(is_w_stations[[x]]))


# Save and load variables for quick execution
# save(usgs_insitu_raw, file=paste0(wd_root, '/tmp_vars/usgs_insitu_raw.RData'))
# load(paste0(wd_root, '/tmp_vars/usgs_insitu_raw.RData'))

# Reduce in-situ SSC measurements to one sample per day
# samplesPerDay <- usgs_insitu_raw %>% 
#   group_by(site_no, sample_dt) %>%
#   summarise(samples = n()) %>%
#   filter(samples > 1)

# Reduce in-situ SSC measurements to one sample per station per day
oneDaySSC <- usgs_insitu_raw %>% 
  group_by(site_no, sample_dt) %>%
  summarise(SSC_mgL = mean(as.numeric(SSC_mgL))) %>%
  left_join(station_info, by = ("site_no"="site_no"))

# Apply log to SSC measurements
is_w_stations$log10_SSC_mgL <- log10(is_w_stations$SSC_mgL)
is_w_stations$log10_SSC_mgL[is_w_stations$log10_SSC_mgL <= 0] <- 0.01 # Roughly 1 mg/L
hist(is_w_stations$log10_SSC_mgL)

oneDaySSC$log10_SSC_mgL <- log10(oneDaySSC$SSC_mgL)
oneDaySSC$log10_SSC_mgL[oneDaySSC$log10_SSC_mgL <= 0] <- 0.01 # Roughly 1 mg/L
hist(oneDaySSC$log10_SSC_mgL)

# Save raw data
write.csv(usgs_insitu_raw, file = paste(wd_exports, "NWIS_WQP_SSC.csv", sep=""), row.names=F)
write.csv(is_w_stations, file = paste(wd_exports, "NWIS_WQP_SSC_STATIONINFO.csv", sep=""), row.names=F)
write.csv(oneDaySSC, file = paste(wd_exports, "ONEDAYSAMPLE_SSC.csv", sep=""), row.names=F)
write.csv(station_info, file = paste(wd_exports, "STATIONINFO.csv", sep=""), row.names=F)

# Also save station info as a shape file
coordinates(station_info) <- ~lon+lat
proj4string(station_info) <- projection
shapefile(station_info, paste0(wd_exports, "STATION_INFO.shp"), overwrite=T)

# Free up some space
rm(usgs_insitu_raw)
rm(is_w_stations)
rm(wqp_insitu_raw)

### Create SSC vs. Q Rating Curve ###

### IMPORT IN SITU DISCHARGE DATA ###
# Import data and make rating curves for all stations
# --> Takes a LONG time
discharge_data <- data.table(
  readNWISdata(
    sites = unique(oneDaySSC$site_no),
    service = "dv",
    parameterCd = "00060",
    statCd="00003",
    startDt = startDate,
    endDt = endDate
  ))[,":=" (
    sample_dt = dateTime,
    Q_m3s = X_00060_00003 * 0.02832 # ft3/s to m3/s
  )][,.(
    site_no,
    sample_dt,
    Q_m3s
    )]
# REMEMBER: flux can be negative

# Convert to all Caps
discharge_data <- mutate_all(discharge_data, toupper)
discharge_data <- discharge_data %>% filter(!is.na(Q_m3s))

# Save and load variables for quick execution
# save(discharge_data, file=paste0(wd_root, '/tmp_vars/discharge_data.RData'))
# load(paste0(wd_root, '/tmp_vars/discharge_data.RData'))
write.csv(discharge_data, file = paste(wd_exports, "NWIS_Q.csv"), row.names = FALSE)

# For loading discharge data
discharge_data <- read.csv(file = paste(wd_exports, "NWIS_Q.csv"))

# Coerce data to numeric after loading CSV
cols <- c("site_no")
discharge_data[, cols] <- lapply(cols, function(x) as.character(discharge_data[[x]]))
cols <- c("Q_m3s")
discharge_data[, cols] <- lapply(cols, function(x) as.numeric(discharge_data[[x]]))

# Coerce 0 discharge to 0.001 to prevent issues with log
discharge_data[discharge_data$Q_m3s == 0, ] <- 0.001

# For negative discharge --> take absolute value
discharge_data[discharge_data$Q_m3s < 0, ] <- abs(discharge_data[discharge_data$Q_m3s < 0, ]$Q_m3s)

# ISSUES MERGING DATA --> Need to tack on days without site data
SSCQ <- left_join(discharge_data, oneDaySSC, by=c("site_no"="site_no", "sample_dt"="sample_dt"))

# Merge station info and discharge data
# discharge_ssc <- left_join(discharge_data, station_info, by=c('site_no'))
ratingSSCQ <- setDT(copy(SSCQ))[
  ,':='(
    agency_cd = 'USGS',
    log10_Q_m3s = log10(Q_m3s),
    log10_SSC_flux_MTyr = log10(Q_m3s * SSC_mgL * 3.10585 * 10**-5),
    SSC_flux_MTyr = Q_m3s * SSC_mgL * 3.10585 * 10**-5,
    R2 = 0.0,
    ratingCurve = FALSE
    )]

# Get site numbers for stations with Q data and >=10 SSC observations
station_nums <- ratingSSCQ %>%
  filter(is.finite(log10_Q_m3s) & is.finite(SSC_mgL)) %>%
  group_by(site_no) %>%
  summarise(samples = n()) %>%
  filter(samples >= 10) %>% # UNCOMMENT THIS
  ungroup() %>%
  pull("site_no")

lm_eqn <- function(m){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
}

pb <- txtProgressBar(0, length(station_nums), style = 3)

# TESTING <- NEED TO EXTRACT DATA FROM DISCHARGE DATA
for (i in seq_along(station_nums)){
  setTxtProgressBar(pb, i)

  # Extract data for the site
  station_data <- ratingSSCQ %>% 
    filter(site_no == station_nums[i])
  
  # Extract data for the site with flux measurements
  station_flux <- station_data %>%
    filter(is.finite(log10_SSC_flux_MTyr) & is.finite(log10_Q_m3s))

  # Remove bad values which mess up calibration
  # station_data = station_data[(is.finite(log10_SSC_flux_MTyr) & is.finite(log10_SSC_flux_MTyr)),]
  # Generate regression
  reg <- lm(log10_SSC_flux_MTyr ~ log10_Q_m3s, data = station_flux)

  # Only keep discharge data within the range of regression --> NOT IMPLEMENTED

  # ??? HOW TO APPLY REGRESSION TO UNSEEN DATA ???
  R2 <- summary(reg)$r.squared
  if (R2 >= 0.80) {
    # Get indexes of samples where there is no SSC (i.e. no in-situ sample)
    goodIDs <- ratingSSCQ$site_no == station_nums[i] & is.na(ratingSSCQ$SSC_mgL)
    # Apply regression
    ssc_flux <- predict.lm(reg, newdata = ratingSSCQ[goodIDs])
    ratingSSCQ[goodIDs]$log10_SSC_flux_MTyr <- ssc_flux
    # Get secondary metrics
    ratingSSCQ[goodIDs]$SSC_flux_MTyr <- 10**ssc_flux
    ratingSSCQ[goodIDs]$SSC_mgL <- (10**ssc_flux) / (ratingSSCQ[goodIDs]$Q_m3s * 3.10585 * 10**-5)
    ratingSSCQ[goodIDs]$log10_SSC_mgL <- log10(ratingSSCQ[goodIDs]$SSC_mgL)
    ratingSSCQ[goodIDs]$R2 <- R2
    ratingSSCQ[goodIDs]$ratingCurve <- TRUE
  }

  # Export regression
  reg_plot <- ggplot(station_flux, aes(x = log10_Q_m3s, y = log10_SSC_flux_MTyr)) +
    geom_point(color = 'green') +
    geom_smooth(method='lm', formula= y~x, color='black',aes(fill = 'standard error')) +
    annotate(geom='text',label = lm_eqn(reg), parse = TRUE, x = -Inf, y = Inf, hjust = -0.2, vjust = 2) +
    theme_bw() +
    scale_fill_manual(values = c('gray'), name = "1 standard Error")  +
    labs(title=paste0('Flux vs Q rating Curve for Station USGS-',station_nums[i]),
         caption=paste0(startDate, " to ", endDate))
  ggsave(reg_plot, filename = paste0(wd_rating, paste0('USGS-',station_nums[i],'.png')),
         width = 10, height = 8)
}

# NUMBER OF TRAINING SAMPLES: 162084
# Check number of modified samples -->
print(nrow(ratingSSCQ %>% filter(is.finite(SSC_mgL))))

# Extract rating curve stations with bad R2
ratingSSCQFiltered <- ratingSSCQ[is.finite(ratingSSCQ$SSC_mgL),]

# Remove stations with excessive SSC
ratingSSCQFiltered <- ratingSSCQFiltered[ratingSSCQFiltered$SSC_mgL <= 5000,]

print("Size of raw in-situ data")
print(nrow(oneDaySSC))

print("Size of rating-curve in-situ data")
print(nrow(ratingSSCQFiltered))
# NOTE: Rating Curve Derived stations DO NOT have station data -> just site_no and values
write.csv(ratingSSCQFiltered, file = paste(wd_exports, "ratingSSCQ.csv"), row.names = FALSE)