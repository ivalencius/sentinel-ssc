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

#### SET DIRECTORIES ####

print('SETTING UP IMPORTS')
# Set root directory
wd_root <- "D:/Thesis/sentinel-ssc/sentinel-calibration"

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
# us_states <- us_states[us_states$state_abbr != "AK" & us_states$state_abbr != "HI" & us_states$state_abbr != "PR",]

plot(us_states_spatial)
# plot(us_states_merge)

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
