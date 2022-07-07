#### LIBRARY IMPORTS ####
library(dataRetrieval)
library(nhdplusTools)
library(tidyhydat)

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

library(data.table)
library(gdata) # for cbindX
library(plyr) # for round any

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

# TO DO/ISSUES
# - Log transform data
# - Extract discharge data
# - Create rating curve

### FOR QUICK USE WITH SENTINEL-CALIBRATION.py ###
# Standardizing
# d1 <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/regression/reg_linear.csv")
# d2 <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/regression/reg_lasso.csv")
# d3 <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/regression/reg_ridge.csv")
# d4 <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/regression/reg_elasticNet.csv")
# No standardizing
# d1 <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/no_standard/regression/reg_linear.csv")
# d2 <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/no_standard/regression/reg_lasso.csv")
# d3 <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/no_standard/regression/reg_ridge.csv")
# d4 <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/no_standard/regression/reg_elasticNet.csv")
# 
# d5 <- read.csv("D:/valencig/Thesis/chattahoochee-dams/chattahoochee-dams-exports/SSC_pred.csv")
# d5 <- d5[as.Date(d5$landsat_dt) > '2017-01-01',]
# d5 <- data.frame(d5$distance_10km, d5$SSC_mgL)
# names(d5) <- c('landsat_distance_10km', 'landsat')
# data = data.table(distance_10km = round_any(d1$distance_km, 10),
#                   linear = d1$pred_SSC_mgL,
#                   lasso = d2$pred_SSC_mgL,
#                   ridge = d3$pred_SSC_mgL,
#                   elasticNet = d4$pred_SSC_mgL)
# 
# # Remove columns without enough appearances
# data <- data[data$distance_10km %in% names(which(table(data$distance_10km) >=200)), ]

# Test filterings
# data$linear <- replace(data$linear, data$linear >= 1000, NA)
# data$linear <- replace(data$linear, data$linear <= -200, NA)
# data$lasso <- replace(data$lasso, data$lasso >= 1000, NA)
# data$lasso <- replace(data$lasso, data$lasso <= -200, NA)
# data$ridge <- replace(data$ridge, data$ridge >= 1000, NA)
# data$ridge <- replace(data$ridge, data$ridge <= -200, NA)
# data$elasticNet <- replace(data$elasticNet, data$elasticNet >= 1000, NA)
# data$elasticNet <- replace(data$elasticNet, data$elasticNet <= -200, NA)

# data <- data[data$linear < 1000 & data$lasso < 1000 & data$ridge < 1000 & data$elasticNet < 1000 & data$linear > -1000 & data$lasso > -1000 & data$ridge > -1000 & data$elasticNet > -1000,]
#data <- join(data, d5, by=('distance_km'='distance_km'))
# data <- cbindX(d5, data)
# 
# ggplot(data, aes(x = distance_10km, y = linear)) +
#   stat_summary(geom = 'line', fun = 'median', aes(color = 'linear')) +
#   stat_summary(geom = 'line', fun = 'median', aes(x = distance_10km, y = lasso, color = 'lasso')) +
#   stat_summary(geom = 'line', fun = 'median', aes(x = distance_10km, y = ridge, color = 'ridge')) +
#   stat_summary(geom = 'line', fun = 'median', aes(x = distance_10km, y = elasticNet, color = 'elasticNet')) +
#   stat_summary(geom = 'line', fun = 'median', aes(x = landsat_distance_10km, y = landsat, color = 'Landsat')) +
#   #facet_wrap(.~paste0('Decade: ', decade)) + # comment out to display mean across all time periods
#   theme_bw()+
#   #geom_hline(yintercept=0)+
#   labs(title='Regressions',
#     x = 'distance downstream',
#        y = 'SSC (mg/L)')

#### SET DIRECTORIES ####

print('SETTING UP IMPORTS')
# Set root directory
wd_root <- "D:/valencig/Thesis/sentinel-ssc/sentinel-calibration"
setwd(wd_root)

# Imports folder (store all import files here)
wd_imports <- paste0(wd_root,"/imports/")
# Exports folder (save all figures, tables here)
wd_exports <- paste0(wd_root,"/exports/")

wd_figures <- paste0(wd_exports, "figures/")
# Sub folders
wd_extent<- paste0(wd_exports, 'station_transects/')
wd_gee <- paste0(wd_exports,'GEE_raw/')
wd_rating <- paste0(wd_exports, 'rating_curves/')

# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_exports, wd_figures, wd_extent, wd_gee, wd_rating)
                         # , wd_exports_gc,wd_station_standalone, 
                         # wd_standalone_models, wd_standalone_figures, wd_autocorrelation)
for(i in 1:length(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}

projection <- CRS("+proj=longlat +datum=WGS84 +no_defs")

#### INITIALIZE MAP DATA FOR N.AMERICA ####

cat('\t','-> Initializing USA map','\n')
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
us_states_spatial <- us_states_spatial[,c('name')] # only select column with province name
# us_states_spatial <- rbind(us_states_spatial,canada_prov) # combine canadian and us shapefiles

# fortify state shapefile
us_ca <- fortify(us_states_spatial)

#### IMPORT IN SITU SSC DATA ####

# From USGS
print('IMPORTING IN SITU DATA')
startDate <- '2017-01-01'
endDate <- '2022-07-07'
cat('\t','-> Start Date:', startDate,'\n')
cat('\t','-> End Date:', endDate,'\n')

#ssc_codes <- c('00530','80154')
ssc_codes <- c('80154')
cat('\t','-> Downloading USGS SSC Data','\n')
usgs_insitu_raw <- data.table(readWQPdata(parameterCd = ssc_codes, 
                                          startDate=startDate,
                                          endDate=endDate))[,":=" (
                                            agency_cd = 'USGS',
                                            site_no = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
                                            sample_dt = ActivityStartDate,
                                            SSC_mgL = ResultMeasureValue,
                                            sample_depth_m = ActivityDepthHeightMeasure.MeasureValue * 0.3048,
                                            sample_method = SampleCollectionMethod.MethodName
                                            )][,.(agency_cd, site_no, sample_dt, SSC_mgL, sample_depth_m, sample_method)]

# Need lat, lon, channel width
usgs_stations <- paste0("USGS-", unique(usgs_insitu_raw$site_no))

# Get properties of usgs WQP stations
cat('\t','-> Downloading USGS Station Metadata','\n')
wqp_info <- data.table(whatWQPsites(siteid=usgs_stations))[,":="(
  site_no = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
  station_nm = MonitoringLocationName,
  drainage_area_km2 = DrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
  lat = LatitudeMeasure,
  lon = LongitudeMeasure,
  elevation_m = VerticalMeasure.MeasureValue *0.3048 # ft -> m
)][,.(site_no, station_nm, lat, lon, drainage_area_km2, elevation_m)]
# Get stream width data '00004'
usgs_width <- data.table(readWQPdata(siteNumbers=usgs_stations,
                                     parameterCd = '00004'))[,":=" (
                                     site_no = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
                                     width_m = ResultMeasureValue * 0.3048
                                      )][,.(site_no, width_m)]
usgs_width <- usgs_width %>% group_by(site_no) %>% summarize(width_m = mean(width_m))
# Match stream width to station data
wqp_info <- left_join(wqp_info, usgs_width, by=("site_no"="site_no"))
# Get rows with no basin area and try to extract it from NLDI extent (may only work for a few sites)
# no_drainage <- wqp_info %>% filter(is.na(wqp_info$drainage_area_km2))
# for (row in 1:nrow(no_drainage)){
#   tryCatch({
#     station_name <- no_drainage$station_nm[i]
#     lat <- no_drainage$lat[i]
#     lon <- no_drainage$lon[i]
#     basin <- findNLDI(location = c(lon, lat), find = 'basin')
#     drainage_m2 <- basin$basin$geometry %>% st_area()
#     drainage_km2 <- drainage_m2 / 1000000
#     no_drainage[i]$drainage_area_km2 <- drainage_km2
#   },
#   error = function(e) {
#     print('no basin found')
#   })
# }
# Join station data and SSC data
usgs_insitu_raw <- left_join(usgs_insitu_raw, wqp_info, by=("site_no"="site_no"))
# Remove na data
#usgs_insitu_raw <- usgs_insitu_raw[!is.na(usgs_insitu_raw$width_m)]
# Remove station data with less than 10 occurances
usgs_insitu_raw <- usgs_insitu_raw[usgs_insitu_raw$site_no %in% names(which(table(usgs_insitu_raw$site_no) >= 10)), ]
# Remove stations less than 70 m wide and drainage area < 3000 km2
usgs_insitu_raw <- usgs_insitu_raw %>% filter(usgs_insitu_raw$width_m >= 70 & usgs_insitu_raw$drainage_area_km2 >= 3000)
# All data stored as characters, convert some columns to numbers
usgs_insitu_raw$lat <- as.numeric(usgs_insitu_raw$lat)
usgs_insitu_raw$lon <- as.numeric(usgs_insitu_raw$lon)
usgs_insitu_raw$drainage_area_km2 <- as.numeric(usgs_insitu_raw$drainage_area_km2)
usgs_insitu_raw$elevation_m <- as.numeric(usgs_insitu_raw$elevation_m)
usgs_insitu_raw$width_m <- as.numeric(usgs_insitu_raw$width_m)
usgs_insitu_raw$SSC_mgL <- as.numeric(usgs_insitu_raw$SSC_mgL)
# Apply log to SSC measurements
usgs_insitu_raw$log10_SSC_mgL <- data.table(log10(usgs_insitu_raw$SSC_mgL))
usgs_insitu_raw$log10_SSC_mgL[usgs_insitu_raw$log10_SSC_mgL == 0] <- 0.01 # Roughly 1 mg/L

### IMPORT IN SITU DISCHARGE DATA ###
discharge_data <- data.table(readNWISdata(
  sites = usgs_insitu_raw$site_no,
  parameterCd = '00060',
  startDt = startDate))[,":=" (
    sample_dt = dateTime,
    discharge_m3s = X_00060_00003 * 0.02832 # ft3/s to m3/s
    )][,.(site_no, sample_dt, discharge_m3s)]
# Filter negative and NA discharge
discharge_data <- discharge_data[discharge_data$discharge_m3s > 0, ]

usgs_insitu_raw <- left_join(usgs_insitu_raw, discharge_data, by=c('site_no','sample_dt'))

# Remove negative discharge and NA discharge
usgs_insitu_raw <- usgs_insitu_raw[usgs_insitu_raw$discharge_m3s >= 0,]

# Save and load variables for quick execution
save(usgs_insitu_raw, file=paste0(wd_root, '/tmp_vars/usgs_insitu_raw.RData'))
load(paste0(wd_root, '/tmp_vars/usgs_insitu_raw.RData'))

### Create Rating Curve ###

# Merge station info and discharge data
discharge_ssc <- left_join(discharge_data, wqp_info, by=c('site_no'))
discharge_ssc <- setDT(discharge_ssc)[
  ,':='(
    agency_cd = 'USGS',
    log10_discharge_m3s = log10(discharge_m3s),
    rating_log10_SSC_flux_MTyr = 0.0,
    rating_SSC_flux_MTyr = 0.0,
    rating_SSC_mgL = 0.0,
    rating_log10_SSC_mgL = 0.0,
    rating_R2 = 0.0
    )]

lm_eqn <- function(m){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

station_nums <- unique(usgs_insitu_raw$site_no)
pb <- txtProgressBar(0, length(station_nums), style = 3)
for (i in 1:length(station_nums)){
  setTxtProgressBar(pb, i)
  # Extract data
  station_data = setDT(usgs_insitu_raw[usgs_insitu_raw$site_no == station_nums[i],])[,":="(
    SSC_flux_MTyr = discharge_m3s * SSC_mgL * 3.10585 * 10**-5,
    log10_SSC_flux_MTyr = log10(discharge_m3s * SSC_mgL * 3.10585 * 10**-5),
    log10_discharge_m3s = log10(discharge_m3s)
  )]
  # Generate regression
  reg = lm(log10_SSC_flux_MTyr~log10_discharge_m3s, data=station_data)
  # Plot regression
  reg_plot = ggplot(station_data, aes(x = log10_discharge_m3s, y=log10_SSC_flux_MTyr)) +
    geom_point(color='green') +
    geom_smooth(method='lm', formula= y~x, color='black',aes(fill = 'standard error')) +
    annotate(geom='text',label = lm_eqn(reg), parse = TRUE, x = -Inf, y = Inf, hjust = -0.2, vjust = 2) +
    theme_bw() +
    scale_fill_manual(values = c('gray'), name = "Metrics")  +
    labs(title=paste0('Rating Curve for Station USGS-',station_nums[i]),
         caption=paste0(startDate,' to present'))
  ggsave(reg_plot, filename = paste0(wd_rating, paste0('USGS-',station_nums[i],'.pdf')),
         width = 10, height = 8)
  # Apply regression to discharge data
  ssc_flux = predict.lm(reg, newdata=discharge_ssc[discharge_ssc$site_no == station_nums[i],])
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$rating_log10_SSC_flux_MTyr <- ssc_flux
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$rating_SSC_flux_MTyr <- 10**ssc_flux
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$rating_SSC_mgL <- (10**ssc_flux) / (discharge_ssc[discharge_ssc$site_no == station_nums[i],]$discharge_m3s * 3.10585 * 10**-5)
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$rating_log10_SSC_mgL <- log10((10**ssc_flux) / (discharge_ssc[discharge_ssc$site_no == station_nums[i],]$discharge_m3s * 3.10585 * 10**-5))
  # Tack on R2 value
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$rating_R2 <- summary(reg)$r.squared
}

# Export data
save(usgs_insitu_raw, file=paste0(wd_root, '/tmp_vars/discharge_ssc.RData'))
load(paste0(wd_root, '/tmp_vars/discharge_ssc.RData'))

### PLOT USGS SITE LOCATIONS###
bare_station <- data.frame(usgs_insitu_raw$lon, usgs_insitu_raw$lat)
names(bare_station) <- c("lon","lat")
bare_station <- bare_station %>% group_by(lon, lat) %>% summarize(num_samples = n())
bare_station_sf <- st_as_sf(bare_station, coords = c("lon", "lat")) %>% st_set_crs(projection)

usgs_no_satellite_plot <- ggplot() +
  geom_map(data = us_ca, map = us_ca, aes(map_id = id), 
           color = "grey30", 
           fill = "grey95", 
           size = 0.25) +
  xlim(-170, -70) + # Cut off Alaska islands for plotting
  geom_sf(data = bare_station_sf,
          aes(size = num_samples),
          color = '#5ab4ac',
          shape = 1,
          fill = NA) +
  scale_size_continuous(breaks= pretty_breaks()) +
  labs(title = 'Viable USGS Stations',
       caption = paste0(startDate, ' to ', endDate, ': total of ', length(unique(usgs_insitu_raw$site_no)), ' sites'))

ggsave(usgs_no_satellite_plot, filename = paste0(wd_figures, 'usgs_stations_no_satellite.pdf'),
       width = 10, height = 8)

### GET NLDI EXTENTS FROM ALL STATIONS ###
# This code only needs to be run once, the shape files it creates are used to pull Sentinel Data from GEE
print('IMPORTING NLDI EXTENTS FOR LANDSAT ACQUISITION')
if(!dir.exists(wd_extent)){
  dir.create(wd_extent)}
cat('\t','-> Saving them to:', wd_extent,'\n')

# Save station points as shape file
station_points <- data.frame(station=wqp_info$site_no, lat=as.numeric(wqp_info$lat), lon=as.numeric(wqp_info$lon))
station_points <- SpatialPointsDataFrame(coords = station_points[,c('lon','lat')], 
                                         data=station_points,
                                         proj4string = projection)
st_write(st_as_sf(station_points), paste0(wd_extent,'station_points.shp'), quiet=TRUE, append=FALSE)

# To test querying sentinel data at varying distances from stations
# station_distances <- c(2, 5, 8, 10, 15, 20)
# station_distances <- c(5)
# # Save first station extent from each distance to make comparison plot
# distance_transects <- data.frame(row.names = c('station','geometry'))
# # Loop over all distances from station
# for (dist in station_distances) {
#   cat('\t','-> Getting transects', dist,'km up/downstream from stations\n')
#   # To save station name and assosiated NLDI extent
#   transects_df <- data.frame(row.names = c('station','geometry'))
#   points_df <- data.frame(row.names = c('station','geometry'))
#   # Progress bar
#   pb <- txtProgressBar(0, nrow(bare_station), style = 3)
#     # Loop over all unique stations
#     for (row in 1:nrow(bare_station)){
#       # print(row)
#       setTxtProgressBar(pb, row)
#       # Extract lat and lon from dataframe of unique station data
#       lon1 <- as.numeric(bare_station[row, "lon"][[1]])
#       lat1 <- as.numeric(bare_station[row, "lat"][[1]])
#       station_num <- usgs_insitu_raw[which(usgs_insitu_raw$lon==lon1 & usgs_insitu_raw$lat==lat1),]$site_no[1]
#       # Get NLDI stream extent up and down main channel
#       tryCatch({
        # extent <- findNLDI(wqp = station_num,
        #                    nav = c('UM', 'DM'),
        #                    find = c("flowlines"),
        #                    distance_km = dist)
#         # Combine upper and lower main channel
#         merged_extent <- st_join(extent$UM_flowlines, extent$DM_flowlines)
#         transects_df <<- rbind(transects_df, data.frame(station = station_num, st_combine(merged_extent$geometry)))
#         },
#       error = function(e) {
#         # Some sites don't have NLDI extents so just store location of station
#         station_point <- st_sfc(st_point(c(lon1, lat1)))
#         st_crs(station_point) <- projection
#         points_df <<- rbind(points_df, data.frame(station = station_num, station_point))
#       })
#       # Combine line extents into one, apply label of station
#       # if (row == nrow(bare_station)){
#       #   distance_transects <- rbind(distance_transects, reduce_extent)
#       # }
#     }
#   close(pb)
#   # Create sf object from dataframe holding station nums and geometry
#   transects_sf <- st_as_sf(transects_df)
#   points_df <- st_as_sf(points_df)
#   # Save transects as one shapefile --> NEED TO FIX
#   st_write(transects_sf, paste0(wd_extent,'transects_',dist,'km.shp'), quiet=TRUE, append=FALSE)
#   st_write(points_sf, paste0(wd_extent,'station_points',dist,'km.shp'), quiet=TRUE, append=FALSE)
# }
# Make comparison plot 
# theme to remove tick marks and axes labels
# tick_theme <- theme_bw() + theme(axis.text.x=element_blank(),
#                       axis.ticks.x=element_blank(),
#                       axis.text.y=element_blank(),
#                       axis.ticks.y=element_blank())
# p1 <- ggplot() + 
#   geom_sf(data = distance_transects[1,],
#           aes(geometry = geometry),
#           color = 'red') +
#   tick_theme +
#   labs(title = '2 km') 
# p2 <- ggplot() + 
#   geom_sf(data = distance_transects[2,],
#           aes(geometry = geometry),
#           color = 'red') +
#   tick_theme +
#   labs(title = '5 km')
# p3 <- ggplot() + 
#   geom_sf(data = distance_transects[3,],
#           aes(geometry = geometry),
#           color = 'red') +
#   tick_theme +
#   labs(title = '8 km')
# p4 <- ggplot() + 
#   geom_sf(data = distance_transects[4,],
#           aes(geometry = geometry),
#           color = 'red') +
#   tick_theme +
#   labs(title = '10 km')
# p5 <- ggplot() + 
#   geom_sf(data = distance_transects[5,],
#           aes(geometry = geometry),
#           color = 'red') +
#   tick_theme +
#   labs(title = '15 km')
# p6 <- ggplot() + 
#   geom_sf(data = distance_transects[6,],
#           aes(geometry = geometry),
#           color = 'red') +
#   tick_theme +
#   labs(title = '20 km')
# 
# distance_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, 
#                               top=textGrob("River Transect Lengths", gp=gpar(fontsize=20,font=3)),
#                               ncol=3)
# ggsave(distance_plot, filename = paste0(wd_figures, 'transect_lengths.pdf'),
#        width = 10, height = 8)

### Clustering ###
# print('RUNNING CLUSTERING ALGORITHMS')
# wd_clusters<- paste0(wd_exports, 'clusters/')
# if(!dir.exists(wd_clusters)){
#   dir.create(wd_clusters)}
# cat('\t','-> Saving cluster data to:', wd_clusters,'\n')
# 
# # Need to not use conda env for final code running
# use_condaenv('sentinel-ssc', required = TRUE)
# system('python -c print(hello world)')

### GET NLDI EXTENTS FOR ALL RIVERS ###

### IMPORT AND HARMONIZE SENTINEL DATA ###
lag_days <- 4
lag_days_examine <- c(1,2,3,4,5,6,7,8,9,10)
gee_files <- list.files(path = wd_gee, pattern = 'sentinel*', full.names = TRUE)
gee_data <- data.table()
for (file in gee_files){
  gee_raw <- read.csv(file)
  gee_clean <- gee_raw[ , !(names(gee_raw) %in% c('system.index','.geo'))]
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
  gee_data <- rbind(gee_data, setDT(gee_clean))
}

num_samples <- data.frame(lag_days_examine, rep(NA, length(lag_days_examine)))
names(num_samples) <- c('Lag_days','Num_samples')

# Raw usgs_data
for (row in 1:nrow(num_samples)){
  day <- num_samples[row,]$Lag_days
  cleaned <- setDT(gee_data)[
    # Create lead lag times and add to DF
    ,':='(match_dt_start = date - day,
          match_dt_end = date + day,
          site_no = unlist(strsplit(station, split='-'))[c(FALSE, TRUE)])][
            # Match by dates inside lead-lag range
            usgs_insitu_raw[,':='(match_dt = as.Date(usgs_insitu_raw$sample_dt))],
            on = .(site_no == site_no, match_dt_start <= match_dt, match_dt_end >= match_dt)][
              ,lag_days := as.numeric(difftime(sample_dt, date),'days')][# Remove rows without images
                !is.na(B1)
              ]
  num_samples[row, 'Num_samples'] <- nrow(cleaned)
}

lag_day_plot <- ggplot(num_samples, aes(x = Lag_days, y = Num_samples, fill = Lag_days)) + 
  theme_bw() +
  geom_bar(stat="identity") +
  labs(title = 'Number of Viable Samples vs. Lag days')

ggsave(lag_day_plot, filename = paste0(wd_figures, 'lag_day_plot.pdf'),
       width = 8, height = 10)

usgs_sentinel_harmonzied <- setDT(gee_data)[
  # Create lead lag times and add to DF
  ,':='(match_dt_start = date - lag_days,
        match_dt_end = date + lag_days,
        site_no = unlist(strsplit(station, split='-'))[c(FALSE, TRUE)])][
          # Match by dates inside lead-lag range
          usgs_insitu_raw[,':='(match_dt = as.Date(usgs_insitu_raw$sample_dt))],
          on = .(site_no == site_no, match_dt_start <= match_dt, match_dt_end >= match_dt)][
            ,lag_days := as.numeric(difftime(sample_dt, date),'days')][
              # Remove duplicated columns
              ,':='(station = NULL, match_dt_start = NULL, match_dt_end = NULL, i.lat = NULL, i.lon = NULL)
              ][ # Add squared columns
                ,':='(B1.2 = B1^2,
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
                      )][# Remove rows without images
                        !is.na(B1)
                        ]
write.csv(usgs_sentinel_harmonzied, paste0(wd_gee, 'ssc_harmonized.csv'), row.names = FALSE)

usgs_sentinel_harmonized <- read.csv(paste0(wd_gee, 'ssc_harmonized.csv'))

# Discharge SSC data
num_samples2 <- data.frame(lag_days_examine, rep(NA, length(lag_days_examine)))
names(num_samples2) <- c('Lag_days','Num_samples')

for (row in 1:nrow(num_samples)){
  day <- num_samples2[row,]$Lag_days
  cleaned2 = setDT(gee_data)[
    # Create lead lag times and add to DF
    ,':='(match_dt_start = date - day,
          match_dt_end = date + day,
          site_no = unlist(strsplit(station, split='-'))[c(FALSE, TRUE)])][
            # Match by dates inside lead-lag range
            discharge_ssc[,':='(match_dt = as.Date(discharge_ssc$sample_dt))],
            on = .(site_no == site_no, match_dt_start <= match_dt, match_dt_end >= match_dt)][
              ,lag_days := as.numeric(difftime(sample_dt, date),'days')][# Remove rows without images
                !is.na(B1)
              ]

  num_samples2[row, 'Num_samples'] <- nrow(cleaned2)
}

lag_day_plot2 <- ggplot(num_samples2, aes(x = Lag_days, y = Num_samples, fill = Lag_days)) + 
  theme_bw() +
  geom_bar(stat="identity") +
  labs(title = 'Number of Viable Samples vs. Lag days')

ggsave(lag_day_plot2, filename = paste0(wd_figures, 'rating_lag_day_plot.pdf'),
       width = 8, height = 10)

usgs_sentinel_harmonzied2 <- setDT(gee_data)[
  # Create lead lag times and add to DF
  ,':='(match_dt_start = date - lag_days,
        match_dt_end = date + lag_days),][
          # Match by dates inside lead-lag range
          discharge_ssc[,':='(match_dt = as.Date(discharge_ssc$sample_dt))],
          on = .(site_no == site_no, match_dt_start <= match_dt, match_dt_end >= match_dt)][
            ,lag_days := as.numeric(difftime(sample_dt, date),'days')][
              # Remove duplicated columns
              ,':='(station = NULL, match_dt_start = NULL, match_dt_end = NULL, i.lat = NULL, i.lon = NULL)
            ][ # Add squared columns
              ,':='(B1.2 = B1^2,
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
              )][# Remove rows without images
                !is.na(B1)
              ]
write.csv(usgs_sentinel_harmonzied2, paste0(wd_gee, 'rating_ssc_harmonized.csv'), row.names = FALSE)

usgs_sentinel_harmonized2 <- read.csv(paste0(wd_gee, 'rating_ssc_harmonized.csv'))

# Make histogram of sample dates
date_vec <- as.Date(usgs_sentinel_harmonized$sample_dt)
date_frame <- data.frame(Lag_days = usgs_sentinel_harmonized$lag_days, sample_dt = date_vec, year = format(date_vec, format='%Y'))
date_plot <- ggplot(date_frame, aes(x=sample_dt)) +
  theme_bw() +
  facet_wrap(~year, scales = "free") +
  geom_bar(stat="count") +
  scale_x_date(breaks="4 month", labels=date_format("%b")) +
  labs(title='Dates of Viable Sentinel-2 Imagery Data',
       caption='From in situ USGS SSC data.')

ggsave(date_plot, filename = paste0(wd_figures, 'acquisition_day.pdf'),
       width = 12, height = 8)

# Make histogram of sample dates for discharge ssc data
date_vec <- as.Date(usgs_sentinel_harmonized2$sample_dt)
date_frame <- data.frame(Lag_days = usgs_sentinel_harmonized2$lag_days, sample_dt = date_vec, year = format(date_vec, format='%Y'))
date_plot <- ggplot(date_frame, aes(x=sample_dt)) +
  theme_bw() +
  facet_wrap(~year, scales = "free") +
  geom_bar(stat="count") +
  scale_x_date(breaks="4 month", labels=date_format("%b")) +
  labs(title='Dates of Viable Sentinel-2 Imagery Data',
       caption='From rating curve derived SSC data.')

ggsave(date_plot, filename = paste0(wd_figures, 'rating_acquisition_day.pdf'),
       width = 12, height = 8)

### CLEAN TRANSECT DATA ###
transect_file <- "D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/exports/GEE_raw/transect/sentinel_2_SR__transect_2017_2022.csv"
transect_raw <- read.csv(transect_file)
transect_clean <- transect_raw[ , !(names(transect_raw) %in% c('system.index','.geo'))]
transect_clean[,'date'] <- as.Date(transect_clean$date)
setnames(transect_clean, 
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

transect_sentinel_harmonzied <- setDT(transect_clean)[ # Add squared columns
              ,':='(B1.2 = B1^2,
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
write.csv(transect_sentinel_harmonzied, paste0(wd_gee, '/transect/transect_harmonized.csv'), row.names = FALSE)
