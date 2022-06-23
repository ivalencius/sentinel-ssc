#### LIBRARY IMPORTS ####
library(dataRetrieval)
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
# - Saving extents only saves first row
# - Some stations don't have width/other information
# - Get more sources of data, only have 85k right now
# - Fix data acquisition dates for Sentinel
# - Implement various clutering methods in python once Sentinel and SSC harmonized

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


# Create folders within root directory to organize outputs if those folders do not exist
export_folder_paths <- c(wd_exports, wd_figures)
                         # , wd_exports_gc,wd_station_standalone, 
                         # wd_standalone_models, wd_standalone_figures, wd_autocorrelation)
for(i in 1:length(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if(!dir.exists(path_sel)){
    dir.create(path_sel)}
}

#### INITIALIZE MAP DATA FOR N.AMERICA ####

cat('\t','-> Initializing USA map','\n')
setwd(wd_imports)
# get states shapefile for clipping/display
us_states_spatial <-  us_states(map_date = NULL, resolution = c("low"), states = NULL) %>% as('Spatial')
# us_states_spatial <- us_states_spatial[us_states_spatial$state_abbr != "AK" & us_states_spatial$state_abbr != "HI" & us_states_spatial$state_abbr != "PR",]
# us_states_spatial <- us_states_spatial[us_states_spatial$state_abbr != "AK",]

#plot(us_states_spatial)

# set projection for states shapefile
projection <- CRS("+proj=longlat +datum=WGS84 +no_defs")

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

#### --- ####
#### IMPORT AND CLEAN -- LANDSAT DATA ####
set.seed(1)

#### IMPORT IN SITU DATA ####

# From USGS
print('IMPORTING IN SITU DATA')
startDate <- '2014-01-01'
endDate <- '2022-06-17'
cat('\t','-> Start Date:', startDate,'\n')
cat('\t','-> End Date:', endDate,'\n')

#ssc_codes <- c('00530','80154')
ssc_codes <- c('80154')
cat('\t','-> Downloading USGS SSC Data','\n')
usgs_insitu_raw <- data.table(readWQPdata(parameterCd = ssc_codes, 
                                          startDate=startDate, # sentinel 1 launched in 2014
                                          endDate=endDate))[,":=" (
                                            agency_cd = 'USGS',
                                            site_no = MonitoringLocationIdentifier,
                                            sample_dt = ActivityStartDate,
                                            SSC_mgL = ResultMeasureValue,
                                            sample_depth_m = ActivityDepthHeightMeasure.MeasureValue * 0.3048,
                                            sample_method = SampleCollectionMethod.MethodName
                                            )][,.(agency_cd, site_no, sample_dt, SSC_mgL, sample_depth_m, sample_method)]

# Need lat, lon, channel width
usgs_stations <- unique(usgs_insitu_raw$site_no)
# Get properties of usgs WQP stations
cat('\t','-> Downloading USGS Station Metadata','\n')
wqp_info <- data.table(whatWQPsites(siteid=usgs_stations))[,":="(
  site_no = MonitoringLocationIdentifier,
  station_nm = MonitoringLocationName,
  drainage_area_km2 = DrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
  lat = LatitudeMeasure,
  lon = LongitudeMeasure,
  elevation_m = VerticalMeasure.MeasureValue *0.3048 # ft -> m
)][,.(site_no, station_nm, lat, lon, drainage_area_km2, elevation_m)]
# Get stream width data '00004'
usgs_width <- data.table(readWQPdata(siteNumbers=usgs_stations,
                                     parameterCd = '00004'))[,":=" (
                                     site_no = MonitoringLocationIdentifier,
                                     width_m = ResultMeasureValue * 0.3048
                                      )][,.(site_no, width_m)]
usgs_width <- usgs_width %>% group_by(site_no) %>% summarize(width_m = mean(width_m))
# Match stream width to station data
wqp_info <- left_join(wqp_info, usgs_width, by=("site_no"="site_no"))
# Join station data and SSC data
usgs_insitu_raw <- left_join(usgs_insitu_raw, wqp_info, by=("site_no"="site_no"))
# Remove na data
usgs_insitu_raw <- usgs_insitu_raw[!is.na(usgs_insitu_raw$width_m)]
# Remove stations less than 30 m wide
usgs_insitu_raw <- usgs_insitu_raw[usgs_insitu_raw$width_m >= 30,]
# All data stored as characters, convert some columns to numbers
usgs_insitu_raw$lat <- as.numeric(usgs_insitu_raw$lat)
usgs_insitu_raw$lon <- as.numeric(usgs_insitu_raw$lon)
usgs_insitu_raw$drainage_area_km2 <- as.numeric(usgs_insitu_raw$drainage_area_km2)
usgs_insitu_raw$elevation_m <- as.numeric(usgs_insitu_raw$elevation_m)
usgs_insitu_raw$width_m <- as.numeric(usgs_insitu_raw$width_m)
usgs_insitu_raw$SSC_mgL <- as.numeric(usgs_insitu_raw$SSC_mgL)

# Save and load variables for quick execution
save(usgs_insitu_raw, file=paste0(wd_root, '/tmp_vars/usgs_insitu_raw.RData'))
load(paste0(wd_root, '/tmp_vars/usgs_insitu_raw.RData'))

### PLOT USGS SITE LOCATIONS: No Landsat check###
bare_station <- data.frame(usgs_insitu_raw$lon, usgs_insitu_raw$lat)
names(bare_station) <- c("lon","lat")
bare_station <- bare_station %>% group_by(lon, lat) %>% summarize(num_samples = n())
bare_station_sf <- st_as_sf(bare_station, coords = c("lon", "lat")) %>% st_set_crs(projection)

usgs_no_landsat_plot <- ggplot() +
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
  labs(title = 'Viable USGS Stations (>= 30m stream width)',
       caption = paste0(startDate, ' to ', endDate, ': total of ', length(unique(usgs_insitu_raw$site_no)), ' sites'))

ggsave(usgs_no_landsat_plot, filename = paste0(wd_figures, 'usgs_stations_no_landsat.pdf'),
       width = 10, height = 8)

### GET NLDI EXTENTS FROM ALL STATIONS ###
print('IMPORTING NLDI EXTENTS FOR LANDSAT ACQUISITION')

wd_extent<- paste0(wd_exports, 'station_transects/')
if(!dir.exists(wd_extent)){
  dir.create(wd_extent)}
cat('\t','-> Saving them to:', wd_extent,'\n')

# To test querying sentinel data at varying distances from stations
station_distances <- c(2, 5, 8, 10, 15, 20)
# station_distances <- c(5)
# Save first station extent from each distance to make comparison plot
distance_transects <- data.frame(row.names = c('station','geometry'))
# Loop over all distances from station
for (dist in station_distances) {
  cat('\t','-> Getting transects', dist,'km up/downstream from stations\n')
  # To save station name and assosiated NLDI extent
  transects_df <- data.frame(row.names = c('station','geometry'))
  # Progress bar
  pb <- txtProgressBar(0, nrow(bare_station), style = 3)
    # Loop over all unique stations
    for (row in 1:nrow(bare_station)){
      setTxtProgressBar(pb, row)
      # Extract lat and lon from dataframe of unique station data
      lon1 <- as.numeric(bare_station[row, "lon"][[1]])
      lat1 <- as.numeric(bare_station[row, "lat"][[1]])
      station_num <- usgs_insitu_raw[which(usgs_insitu_raw$lon==lon1 & usgs_insitu_raw$lat==lat1),]$site_no[1]
      # Get NLDI stream extent up and down main channel
      tryCatch({
        extent <- findNLDI(wqp = station_num,
                           nav = c('UM', 'DM'),
                           distance_km = dist)
        # Combine upper and lower main channel
        merged_extent <- st_join(extent$UM_flowlines, extent$DM_flowlines)
      },
      error = function(e) {
        # Some sites don't have NLDI extents so just store location of station
        merged_extent <- findNLDI(wqp = station_num)
      })
      # Combine line extents into one, apply label of station
      reduce_extent <- data.frame(station = station_num, st_combine(merged_extent$geometry))
      transects_df <- rbind(transects_df, reduce_extent)
      if (row == 1){
        distance_transects <- rbind(distance_transects, reduce_extent)
      }
    }
  close(pb)
  # Create sf object from dataframe holding station nums and geometry
  transects_sf <- st_as_sf(transects_df)
  # Save transects as one shapefile --> NEED TO FIX
  st_write(transects_sf, paste0(wd_extent,'transects_',dist,'km.shp'), quiet=TRUE, append=FALSE)
}
#ggplot() + geom_sf(data=transects_sf)
# Make comparison plot 
# theme to remove tick marks and axes labels
tick_theme <- theme_bw() + theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank())
p1 <- ggplot() + 
  geom_sf(data = distance_transects[1,],
          aes(geometry = geometry),
          color = 'red') +
  tick_theme +
  labs(title = '2 km') 
p2 <- ggplot() + 
  geom_sf(data = distance_transects[2,],
          aes(geometry = geometry),
          color = 'red') +
  tick_theme +
  labs(title = '5 km')
p3 <- ggplot() + 
  geom_sf(data = distance_transects[3,],
          aes(geometry = geometry),
          color = 'red') +
  tick_theme +
  labs(title = '8 km')
p4 <- ggplot() + 
  geom_sf(data = distance_transects[4,],
          aes(geometry = geometry),
          color = 'red') +
  tick_theme +
  labs(title = '10 km')
p5 <- ggplot() + 
  geom_sf(data = distance_transects[5,],
          aes(geometry = geometry),
          color = 'red') +
  tick_theme +
  labs(title = '15 km')
p6 <- ggplot() + 
  geom_sf(data = distance_transects[6,],
          aes(geometry = geometry),
          color = 'red') +
  tick_theme +
  labs(title = '20 km')

distance_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, 
                              top=textGrob("River Transect Lengths", gp=gpar(fontsize=20,font=3)),
                              ncol=3)
ggsave(distance_plot, filename = paste0(wd_figures, 'transect_lengths.pdf'),
       width = 10, height = 8)