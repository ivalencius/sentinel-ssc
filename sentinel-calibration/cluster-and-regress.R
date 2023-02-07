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
us_states_spatial <- us_states_spatial[,c('state_abbr')] # only select column with province name
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
wqp_insitu_raw <- data.table(readWQPdata(parameterCd = ssc_codes,
                                          startDate=startDate,
                                          endDate=endDate))[,":=" (
                                            agency_cd = 'USGS',
                                            site_no = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
                                            sample_dt = ActivityStartDate,
                                            SSC_mgL = ResultMeasureValue,
                                            sample_depth_m = ActivityDepthHeightMeasure.MeasureValue * 0.3048,
                                            sample_method = SampleCollectionMethod.MethodName
                                            )][,.(agency_cd, site_no, sample_dt, SSC_mgL, sample_depth_m, sample_method)]

# Need to use WQP data as well
nwis_insitu_raw <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c('agency_cd', 'site_no', 'sample_dt', 'SSC_mgL'))))
for (state in us_states_spatial$state_abbr){
  nwis_data <- data.table(readNWISdata(stateCd=state,
                                     parameterCd = ssc_codes, 
                                            startDt=startDate,
                                            endDt=endDate))
  if(nrow(nwis_data)==0){} # Do nothing if no data returned
  else{
    nwis_data_clean <- nwis_data[,":=" (agency_cd = agency_cd, 
                                              site_no = site_no,
                                              sample_dt = dateTime,
                                              SSC_mgL = X_80154_00003
                                            )][,.(agency_cd, site_no, sample_dt, SSC_mgL)]
    nwis_insitu_raw <- rbind(nwis_insitu_raw, nwis_data_clean)
  }
}
nwis_insitu_raw <- nwis_insitu_raw %>% filter(SSC_mgL >= 0)

# Remove discrete wqp data recorded on same day/location as nwis data (to prevent over counting one day)
# pb <- txtProgressBar(0, nrow(wqp_insitu_raw), style = 3)
# for (i in 1:nrow(wqp_insitu_raw)){
#   setTxtProgressBar(pb, i)
#   nwis_sample <- nwis_insitu_raw %>% filter(sample_dt==wqp_insitu_raw$sample_dt[i] & site_no==wqp_insitu_raw$site_no[i])
#   if (nrow(nwis_sample==0)){
#     wqp_insitu_raw$site_no[i] <- -1 # Set -1 flag and filter after
#   }
# }
# wqp_insitu_raw <- wqp_insitu_raw %>% filter(site_no > 0)

# Merge nwis and usgs
usgs_insitu_raw <- bind_rows(data.frame(wqp_insitu_raw), data.frame(nwis_insitu_raw))
# Filter out nonsense data
usgs_insitu_raw <- usgs_insitu_raw %>% filter(SSC_mgL >=0)
# Convert to all Caps
usgs_insitu_raw <- mutate_all(usgs_insitu_raw, toupper)
# # Save raw data
# write.csv(usgs_insitu_raw, file = paste(wd_exports,"NWIS_WQP_SSC.csv"), row.names=F)
# Aggregate data taken on the same day
# usgs_insitu_unique <- usgs_insitu_raw %>% group_by(site_no, )
# usgs_insitu_raw <- nwis_insitu_raw

# For WQP need USGS- prefix
usgs_stations <- paste0("USGS-", unique(usgs_insitu_raw$site_no))

# Get properties of usgs WQP stations
cat('\t','-> Downloading USGS Station Metadata','\n')
station_info <- data.table(whatWQPsites(siteid=usgs_stations))[,":="(
  station_nm = MonitoringLocationName,
  drainage_area_km2 = DrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
  contributing_drainage_area_km2 = ContributingDrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
  lat = LatitudeMeasure,
  lon = LongitudeMeasure,
  elevation_m = VerticalMeasure.MeasureValue *0.3048, # ft -> m
  site_no = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
  site_type = MonitoringLocationTypeName,
  HUC = HUCEightDigitCode
)][,.(site_no,
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
is_w_stations[, cols] <- lapply(cols, function(x) as.numeric(is_w_stations[[x]]))

# Apply log to SSC measurements
is_w_stations$log10_SSC_mgL <- log10(is_w_stations$SSC_mgL)
is_w_stations$log10_SSC_mgL[is_w_stations$log10_SSC_mgL <= 0] <- 0.01 # Roughly 1 mg/L
hist(is_w_stations$log10_SSC_mgL)

# Save and load variables for quick execution
# save(usgs_insitu_raw, file=paste0(wd_root, '/tmp_vars/usgs_insitu_raw.RData'))
# load(paste0(wd_root, '/tmp_vars/usgs_insitu_raw.RData'))

# Save raw data
write.csv(usgs_insitu_raw, file = paste(wd_exports,"NWIS_WQP_SSC.csv"), row.names=F)
write.csv(is_w_stations, file = paste(wd_exports,"NWIS_WQP_SSC_STATIONINFO.csv"), row.names=F)

### Create SSC vs. Q Rating Curve ###

### Investigate SSC Samples ###
samplesPerDay <- group_by(usgs_insitu_raw, site_no, sample_dt)

### IMPORT IN SITU DISCHARGE DATA ###
# Import data and make rating curves for all stations, remove unused ones later (to save time if multiple gee code running)
discharge_data <- data.table(readNWISdata( # TAKES A LONG TIME
  sites = unique(is_w_stations$site_no),
  parameterCd = '00060',
  startDt = startDate))[,":=" (
    sample_dt = dateTime,
    Q_m3s = X_00060_00003 * 0.02832 # ft3/s to m3/s
  )][,.(site_no, 
        sample_dt, 
        Q_m3s
        )]

# Convert to all Caps
discharge_data <- mutate_all(discharge_data, toupper)

# Save and load variables for quick execution
# save(discharge_data, file=paste0(wd_root, '/tmp_vars/discharge_data.RData'))
# load(paste0(wd_root, '/tmp_vars/discharge_data.RData'))
write.csv(discharge_data, file = paste(wd_exports,"NWIS_Q.csv"), row.names=F)

# Filter negative and NA discharge
# discharge_data <- discharge_data %>% filter(Q_m3s>=0)

SSCQ <- left_join(is_w_stations, discharge_data, by=c('site_no','sample_dt'))

# Remove stations without discharge data or na SSC columns
# usgs_insitu_raw <- usgs_insitu_raw[!is.na(usgs_insitu_raw$discharge_m3s),]
# usgs_insitu_raw <- usgs_insitu_raw[!is.na(usgs_insitu_raw$SSC_mgL),]

# Merge station info and discharge data
# discharge_ssc <- left_join(discharge_data, station_info, by=c('site_no'))
ratingSSCQ <- setDT(SSCQ)[
  ,':='(
    agency_cd = 'USGS',
    log10_discharge_m3s = log10(discharge_m3s), # already filtered for Q>=0
    log10_SSC_flux_MTyr = 0.0,
    SSC_flux_MTyr = 0.0,
    SSC_mgL = 0.0,
    log10_SSC_mgL = 0.0,
    R2 = 0.0
    )]

lm_eqn <- function(m){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# Only perform regression on stations with >10 samples per day
station_nums <- na.omit(unique(ratingSSCQ$site_no))

pb <- txtProgressBar(0, length(station_nums), style = 3)
for (i in 1:length(station_nums)){
  setTxtProgressBar(pb, i)
  # Extract data
  station_data = as.data.table(usgs_insitu_raw[usgs_insitu_raw$site_no == station_nums[i],])[,":="(
    SSC_flux_MTyr = discharge_m3s * SSC_mgL * 3.10585 * 10**-5,
    log10_SSC_flux_MTyr = log10(discharge_m3s * SSC_mgL * 3.10585 * 10**-5),
    log10_discharge_m3s = log10(discharge_m3s)
  )]
  # Remove bad values which mess up calibration
  station_data = station_data[(is.finite(log10_SSC_flux_MTyr) & is.finite(log10_SSC_flux_MTyr)),]
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
  newdata = discharge_ssc[discharge_ssc$site_no == station_nums[i],]
  # Only keep discharge data within the range of regression
  newdata[!(newdata$log10_discharge_m3s <= max(station_data$log10_discharge_m3s) & newdata$log10_discharge_m3s >= min(station_data$log10_discharge_m3s)),] = NA
  ssc_flux = predict.lm(reg, newdata=newdata)
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$log10_SSC_flux_MTyr <- ssc_flux
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$SSC_flux_MTyr <- 10**ssc_flux
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$SSC_mgL <- (10**ssc_flux) / (discharge_ssc[discharge_ssc$site_no == station_nums[i],]$discharge_m3s * 3.10585 * 10**-5)
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$log10_SSC_mgL <- log10((10**ssc_flux) / (discharge_ssc[discharge_ssc$site_no == station_nums[i],]$discharge_m3s * 3.10585 * 10**-5))
  # Tack on R2 value
  discharge_ssc[discharge_ssc$site_no == station_nums[i],]$R2 <- summary(reg)$r.squared
}

# Remove stations with bad fit
discharge_ssc <- discharge_ssc[discharge_ssc$R2>=0.8,]

# Remove stations with excessive ssc
discharge_ssc <- discharge_ssc[discharge_ssc$SSC_mgL<50000,]

# Export data
# save(discharge_ssc, file=paste0(wd_root, '/tmp_vars/discharge_ssc.RData'))
# load(paste0(wd_root, '/tmp_vars/discharge_ssc.RData'))

### GET NLDI EXTENTS FROM ALL STATIONS ###

### IMPORT AND HARMONIZE SENTINEL DATA ###
# helper function -> pads a character string with leading 0s (for USGS codes)
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
lag_days <- 4
lag_days_examine <- c(0,1,2,3,4,5,6,7,8)
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
# Site numbers stored as numbers in gee so lose leading zero
gee_data$site_no <- unlist(lapply(gee_data$site_no, pad0))

### Cluster each site ###
cluster_data<-setDT(gee_data)[
  # Create lead lag times and add to DF
  ,':='(site_no = as.character(site_no)),]
# All band ratios 

### PLOT USGS SITE LOCATIONS###
# bare_station <- data.frame(gee_data$lon, gee_data$lat)
# names(bare_station) <- c("lon","lat")
# bare_station <- bare_station %>% group_by(lon, lat) %>% summarize(num_samples = n())
# bare_station_sf <- st_as_sf(bare_station, coords = c("lon", "lat")) %>% st_set_crs(projection)
# 
# usgs_no_satellite_plot <- ggplot() +
#   geom_map(data = us_ca, map = us_ca, aes(map_id = id),
#            color = "grey30",
#            fill = "grey95",
#            size = 0.25) +
#   xlim(-170, -70) + # Cut off Alaska islands for plotting
#   geom_sf(data = bare_station_sf,
#           aes(size = num_samples, color=mean_width),
#           color = '#5ab4ac',
#           shape = 1,
#           fill = NA) +
#   scale_size_continuous(breaks = pretty_breaks()) +
#   labs(title = 'Viable USGS Stations',
#        caption = paste0(startDate, ' to ', endDate, ': total of ', length(unique(gee_data$site_no)), ' sites'))

# ggsave(usgs_no_satellite_plot, filename = paste0(wd_figures, 'usgs_stations_no_satellite.pdf'),
#        width = 10, height = 8)
# 
# num_samples <- data.frame(lag_days_examine, rep(NA, length(lag_days_examine)))
# names(num_samples) <- c('Lag_days','Num_samples')

# Raw usgs_data
# for (row in 1:nrow(num_samples)){
#   day <- num_samples[row,]$Lag_days
#   cleaned <- gee_data[
#     # Create lead lag times and add to DF
#     ,':='(match_dt_start = date - day,
#           match_dt_end = date + day,
#           site_no = site_no)][
#             # Match by dates inside lead-lag range
#             as.data.table(usgs_insitu_raw)[,':='(match_dt = as.Date(usgs_insitu_raw$sample_dt))],
#             on = .(site_no == site_no, match_dt_start <= match_dt, match_dt_end >= match_dt)][
#               ,lag_days := as.numeric(difftime(sample_dt, date),'days')][# Remove rows without images
#                 !is.na(B1)
#               ]
#   num_samples[row, 'Num_samples'] <- nrow(cleaned)
# }
# 
# lag_day_plot <- ggplot(num_samples, aes(x = Lag_days, y = Num_samples, fill = Lag_days)) + 
#   theme_bw() +
#   geom_bar(stat="identity") +
#   labs(title = 'Number of Viable Samples vs. Lag days')
# 
# ggsave(lag_day_plot, filename = paste0(wd_figures, 'lag_day_plot.pdf'),
#        width = 8, height = 10)

usgs_sentinel_harmonized <- setDT(gee_data)[
  # Create lead lag times and add to DF
  ,':='(match_dt_start = date - lag_days,
      match_dt_end = date + lag_days,
      site_no = site_no)][
        # Match by dates inside lead-lag range
        as.data.table(usgs_insitu_raw)[,':='(match_dt = as.Date(usgs_insitu_raw$sample_dt))],
        on = .(site_no == site_no, match_dt_start <= match_dt, match_dt_end >= match_dt)][
          ,lag_days := as.numeric(difftime(sample_dt, date),'days')][
            # Remove duplicated columns
            ,':='(match_dt_start = NULL, match_dt_end = NULL, i.lat = NULL, i.lon = NULL)
          ][ # Add squared columns
                ,':='(
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
                      )][# Remove rows without images
                        !is.na(B1)
                        ]
# Assign clusters
# usgs_sentinel_harmonized <- left_join(usgs_sentinel_harmonized, cluster_stations, by=('site_no'='site_no'))
# Remove NA data
usgs_sentinel_harmonized <- usgs_sentinel_harmonized[usgs_sentinel_harmonized$log10_SSC_mgL <=4.5,]
usgs_sentinel_harmonized[usgs_sentinel_harmonized$SSC_mgL==0,]$log10_SSC_mgL <- 0.01

# Reduce to one in-situ station per product_id
usgs_sentinel_harmonized <- usgs_sentinel_harmonized %>% group_by(site_no, product_id) %>% 
    slice(which.min(abs(lag_days))) # just take first sample if there are multiple on one

# Reduce to one image per site --> NEED TO FIX
usgs_sentinel_harmonized <- usgs_sentinel_harmonized %>% group_by(site_no, sample_dt) %>%
  summarise(across(everything(), mean),
            .groups = 'drop')  %>%
  as.data.frame()


# Discharge SSC data
# num_samples2 <- data.frame(lag_days_examine, rep(NA, length(lag_days_examine)))
# names(num_samples2) <- c('Lag_days','Num_samples')
# 
# for (row in 1:nrow(num_samples2)){
#   day <- num_samples2[row,]$Lag_days
#   cleaned2 = setDT(gee_data)[
#     # Create lead lag times and add to DF
#     ,':='(match_dt_start = date - day,
#           match_dt_end = date + day,
#           site_no = site_no)][
#             # Match by dates inside lead-lag range
#             discharge_ssc[,':='(match_dt = as.Date(discharge_ssc$sample_dt))],
#             on = .(site_no == site_no, match_dt_start <= match_dt, match_dt_end >= match_dt)][
#               ,lag_days := as.numeric(difftime(sample_dt, date),'days')][# Remove rows without images
#                 !is.na(B1)
#               ]
# 
#   num_samples2[row, 'Num_samples'] <- nrow(cleaned2)
# }
# 
# lag_day_plot2 <- ggplot(num_samples2, aes(x = Lag_days, y = Num_samples, fill = Lag_days)) + 
#   theme_bw() +
#   geom_bar(stat="identity") +
#   labs(title = 'Number of Viable Samples vs. Lag days')
# 
# ggsave(lag_day_plot2, filename = paste0(wd_figures, 'rating_lag_day_plot.pdf'),
#        width = 8, height = 10)

rating_ssc_harmonized <- setDT(gee_data)[
  # Create lead lag times and add to DF
  ,':='(match_dt_start = date - lag_days,
        match_dt_end = date + lag_days),][
          # Match by dates inside lead-lag range
          discharge_ssc[,':='(match_dt = as.Date(discharge_ssc$sample_dt))],
          on = .(site_no == site_no, match_dt_start <= match_dt, match_dt_end >= match_dt)][
            ,lag_days := as.numeric(difftime(sample_dt, date),'days')][
              # Remove duplicated columns
              ,':='(match_dt_start = NULL, match_dt_end = NULL, i.lat = NULL, i.lon = NULL)
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
              )][# Remove rows without images
                !is.na(B1)
              ]
# Assign clusters
# rating_ssc_harmonized <- left_join(rating_ssc_harmonized, cluster_stations, by=('site_no'='site_no'))
# Take minimum lag day
rating_ssc_harmonized <- rating_ssc_harmonized %>% group_by(site_no, product_id) %>% 
  slice(which.min(abs(lag_days)))

# Delete data already present in in-situ data
# Cast site_no to double and pad
usgs_sentinel_harmonized$site_no <- unlist(lapply(usgs_sentinel_harmonized$site_no, pad0))
rating_ssc_harmonized$sample_dt <- as.character(rating_ssc_harmonized$sample_dt)
rating_ssc_harmonized$date <- as.character(rating_ssc_harmonized$date)
usgs_sentinel_harmonized$date <- as.character(usgs_sentinel_harmonized$date)
usgs_sentinel_harmonized$sample_dt <- as.character(usgs_sentinel_harmonized$sample_dt)
rating_ssc_harmonized <- anti_join(rating_ssc_harmonized, usgs_sentinel_harmonized, 
                                   by=c("site_no"="site_no","sample_dt"="sample_dt","product_id"="product_id"))

# Merge with in-situ data
rating_ssc_harmonized <- rbind(usgs_sentinel_harmonized, rating_ssc_harmonized)

# Remove data > 4.5 log 10
rating_ssc_harmonized <- rating_ssc_harmonized %>% filter(log10_SSC_mgL <=4.5)
rating_ssc_harmonized[rating_ssc_harmonized$SSC_mgL==0,]$log10_SSC_mgL <- 0.01

# Reduce to one sample per day

### CLUSTER --> CURRENTLY USING USGS DATA ###
# regressors <- c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11', 
#                 'B12','B2.B1','B3.B1','B4.B1',
#                 'B2.B1','B3.B1','B4.B1','B5.B1','B6.B1','B7.B1','B8.B1',
#                 'B8A.B1','B9.B1','B11.B1','B12.B1','B3.B2','B4.B2','B5.B2',
#                 'B6.B2','B7.B2','B8.B2','B8A.B2','B9.B2','B11.B2','B12.B2',
#                 'B4.B3','B5.B3','B6.B3','B7.B3','B8.B3','B8A.B3','B9.B3',
#                 'B11.B3','B12.B3','B5.B4','B6.B4','B7.B4','B8.B4','B8A.B4',
#                 'B9.B4','B11.B4','B12.B4','B6.B5','B7.B5','B8.B5','B8A.B5',
#                 'B9.B5','B11.B5','B12.B5','B7.B6','B8.B6','B8A.B6','B9.B6',
#                 'B11.B6','B12.B6','B8.B7','B8A.B7','B9.B7','B11.B7','B12.B7',	
#                 'B8A.B8','B9.B8','B11.B8','B12.B8','B9.B8A','B11.B8A',
#                 'B12.B8A','B11.B9','B12.B9','B12.B11')
regressors <- c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')
cluster_medians <- data.table(matrix(nrow=0, ncol=length(regressors)+1))
names(cluster_medians) <- c('site_no', regressors)
unique_sites <- unique(as.character(usgs_sentinel_harmonized$site_no))
pb <- txtProgressBar(0, length(unique_sites), style = 3)
for (i in 1:length(unique_sites)){
  setTxtProgressBar(pb, i)
  site = unique_sites[i]
  # Extract data
  site_data <- as.data.table(usgs_sentinel_harmonized[usgs_sentinel_harmonized$site_no == site,])[, ..regressors] %>% as.matrix()
  # Reduce by median
  site_median <- colMedians(site_data, na.rm=T)
  # Create table to hold site data
  site_df <- data.table(t(data.table(site_median)))
  names(site_df) <- regressors
  site_df$site_no <- site
  # Generate regression
  cluster_medians <- rbind(cluster_medians, site_df)
}
# Replace NA values with 0
# 6 clusters to harmonize with Evans data
clusters_calculated <- kmeans(cluster_medians[, ..regressors], centers = 6,
                              nstart = 100, iter.max = 1000)

# Compute cluster centers
cluster_centers <- clusters_calculated$centers
# Assign cluster to each site
cluster_medians$cluster <- clusters_calculated$cluster
# Extract cluster for each station
cluster_stations <- cluster_medians[,.(site_no, cluster)]
usgs_sentinel_harmonized$site_no <- as.character(usgs_sentinel_harmonized$site_no)
usgs_sentinel_harmonized <- left_join(usgs_sentinel_harmonized, cluster_stations, by=c("site_no"="site_no"))

# Visualize Data'
ssc_categories <- c(0,50,100,250,500,750,1e6)
# ssc_categories <- c(0,50,100,200,500,1e6)
# ssc_categories <- c(0,10,25,50,75,100,150,200,250,300,350, 400, 450, 500,600, 700, 800,900,1000,1100,1500, 1e6)

# Generate SSC labels as 'low value' - 'high value'
ssc_category_labels <- paste0(ssc_categories[-length(ssc_categories)],'-',c(ssc_categories[-1]))
# Make highest SSC category "> highest value"
ssc_category_labels[length(ssc_category_labels)] <- paste0('> ', ssc_categories[length(ssc_category_labels)])

viz_df <- as.data.table(usgs_sentinel_harmonized)[
  ,':='(cluster_sel = cluster,
        # # Categorize SSC value as one of selected categories
        ssc_category = cut(10^log10_SSC_mgL, 
                           breaks = ssc_categories,
                           labels = ssc_category_labels))][]

# Generate median B,G,R for each SSC category and each cluster or site
ssc_category_color <- viz_df %>% group_by(cluster_sel, ssc_category) %>% 
  summarise_at(vars("B4","B3","B2"), mean)

raster_color_types <- geom_raster(aes(fill = rgb(B4/2200,B3/2200,B2/2200))) # true color

ggplot(ssc_category_color, aes(x = cluster_sel, y = ssc_category)) +
raster_color_types +
scale_fill_identity() +
  theme_classic()+
# scale_x_continuous(expand_scale(add = c(0,0))) + 
# scale_y_discrete(expand_scale(mult = c(0,0))) +
theme(axis.text.x = element_text(angle = 90)) + 
labs(
  y = 'SSC range (mg/L)',
  x = 'River grouping')

### Cluster based on rating curve derived data ###
cluster_medians2 <- data.table(matrix(nrow=0, ncol=length(regressors)+1))
names(cluster_medians2) <- c('site_no', regressors)
unique_sites2 <- unique(rating_ssc_harmonized$site_no)
pb <- txtProgressBar(0, length(unique_sites2), style = 3)
for (i in 1:length(unique_sites2)){
  setTxtProgressBar(pb, i)
  site = unique_sites2[i]
  # Extract data
  site_data <- as.data.table(rating_ssc_harmonized[rating_ssc_harmonized$site_no == site,])[, ..regressors] %>% as.matrix()
  # Reduce by median
  site_median <- colMedians(site_data, na.rm=T)
  # Create table to hold site data
  site_df <- data.table(t(data.table(site_median)))
  names(site_df) <- regressors
  site_df$site_no <- site
  # Generate regression
  cluster_medians2 <- rbind(cluster_medians2, site_df)
}
# Replace NA values with 0
# 6 clusters to harmonize with Evans data
clusters_calculated2 <- kmeans(cluster_medians2[, ..regressors], centers = 6,
                              nstart = 100, iter.max = 1000)

# Compute cluster centers
cluster_centers2 <- clusters_calculated2$centers
# Assign cluster to each site
cluster_medians2$cluster <- clusters_calculated2$cluster
# Extract cluster for each station
cluster_stations2 <- cluster_medians2[,.(site_no, cluster)]
rating_ssc_harmonized <- left_join(rating_ssc_harmonized, cluster_stations2, by=c("site_no"="site_no"))
# Visualize Data'
ssc_categories <- c(0,50,100,250,500,750,1e6)
# ssc_categories <- c(0,50,100,200,500,1e6)
# ssc_categories <- c(0,10,25,50,75,100,150,200,250,300,350, 400, 450, 500,600, 700, 800,900,1000,1100,1500, 1e6)

# Generate SSC labels as 'low value' - 'high value'
ssc_category_labels <- paste0(ssc_categories[-length(ssc_categories)],'-',c(ssc_categories[-1]))
# Make highest SSC category "> highest value"
ssc_category_labels[length(ssc_category_labels)] <- paste0('> ', ssc_categories[length(ssc_category_labels)])

viz_df <- as.data.table(rating_ssc_harmonized)[
  ,':='(cluster_sel = cluster,
        # # Categorize SSC value as one of selected categories
        ssc_category = cut(10^log10_SSC_mgL,
                           breaks = ssc_categories,
                           labels = ssc_category_labels))][]

# Generate median B,G,R for each SSC category and each cluster or site
ssc_category_color <- viz_df%>%group_by(cluster_sel, ssc_category) %>%
  summarise_at(vars("B4","B3","B2"), mean)

raster_color_types <- geom_raster(aes(fill = rgb(B4/2200,B3/2200,B2/2200))) # true color

ggplot(ssc_category_color, aes(x = cluster_sel, y = ssc_category)) +
  raster_color_types +
  scale_fill_identity() +
  theme_classic()+
  # scale_x_continuous(expand_scale(add = c(0,0))) +
  # scale_y_discrete(expand_scale(mult = c(0,0))) +
  # theme(axis.text.x = element_text(angle = 90)) +
  labs(
    y = 'SSC range (mg/L)',
    x = 'River grouping')

### Save variables ###
save(clusters_calculated, file=paste0(wd_root, '/tmp_vars/clusters_calculated.RData'))
save(clusters_calculated2, file=paste0(wd_root, '/tmp_vars/clusters_calculated2.RData'))
write.csv(rating_ssc_harmonized, paste0(wd_gee, 'rating_ssc_harmonized.csv'), row.names = FALSE)
write.csv(usgs_sentinel_harmonized, paste0(wd_gee, 'ssc_harmonized.csv'), row.names = FALSE)

### Load variables ###
load(paste0(wd_root, '/tmp_vars/clusters_calculated.RData'))
load(paste0(wd_root, '/tmp_vars/clusters_calculated2.RData'))
rating_ssc_harmonized <- read.csv(paste0(wd_gee, 'rating_ssc_harmonized.csv'))
usgs_sentinel_harmonized <- read.csv(paste0(wd_gee, 'ssc_harmonized.csv'))

# Make histogram of sample dates
# date_vec <- as.Date(usgs_sentinel_harmonized$sample_dt)
# date_frame <- data.frame(Lag_days = usgs_sentinel_harmonized$lag_days, sample_dt = date_vec, year = format(date_vec, format='%Y'))
# date_plot <- ggplot(date_frame, aes(x=sample_dt)) +
#   theme_bw() +
#   facet_wrap(~year, scales = "free") +
#   geom_bar(stat="count") +
#   scale_x_date(breaks="4 month", labels=date_format("%b")) +
#   labs(title='Dates of Viable Sentinel-2 Imagery Data',
#        caption='From in situ USGS SSC data.')
# 
# ggsave(date_plot, filename = paste0(wd_figures, 'acquisition_day.pdf'),
#        width = 12, height = 8)

# Make histogram of sample dates for discharge ssc data
# date_vec <- as.Date(rating_ssc_harmonized$sample_dt)
# date_frame <- data.frame(Lag_days = rating_ssc_harmonized$lag_days, sample_dt = date_vec, year = format(date_vec, format='%Y'))
# date_plot <- ggplot(date_frame, aes(x=sample_dt)) +
#   theme_bw() +
#   facet_wrap(~year, scales = "free") +
#   geom_bar(stat="count") +
#   scale_x_date(breaks="4 month", labels=date_format("%b")) +
#   labs(title='Dates of Viable Sentinel-2 Imagery Data',
#        caption='From rating curve derived SSC data.')
# 
# ggsave(date_plot, filename = paste0(wd_figures, 'rating_acquisition_day.pdf'),
#        width = 12, height = 8)

# Make histogram of sample dates for discharge ssc data (no multiples)
# date_vec <- as.Date(rating_no_multiples$sample_dt)
# date_frame <- data.frame(Lag_days = rating_no_multiples$lag_days, sample_dt = date_vec, year = format(date_vec, format='%Y'))
# date_plot <- ggplot(date_frame, aes(x=sample_dt)) +
#   theme_bw() +
#   facet_wrap(~year, scales = "free") +
#   geom_bar(stat="count") +
#   scale_x_date(breaks="4 month", labels=date_format("%b")) +
#   labs(title='Dates of Viable Sentinel-2 Imagery Data',
#        caption='From rating curve derived SSC data (no multiples).')
# 
# ggsave(date_plot, filename = paste0(wd_figures, 'rating_acquisition_day_no_multiples.pdf'),
#        width = 12, height = 8)

### RUN BELOW TO IMPORT AND CLEAN NEW DATA ###
transect_file <- "D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/exports/GEE_raw/transect/sentinel_2_SR__transect_2017_2022_20scale.csv"
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

transect_clean <- as.data.table(transect_clean)[,':='(
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

cluster_medians <- data.table(matrix(nrow=0, ncol=length(regressors)+1))
names(cluster_medians) <- c('distance_km', regressors)
unique_dist <- unique(transect_clean$distance_km)
pb <- txtProgressBar(0, length(unique_dist), style = 3)
for (i in 1:length(unique_dist)){
  setTxtProgressBar(pb, i)
  dist = unique_dist[i]
  # Extract data
  dist_data <- transect_clean[transect_clean$distance_km == dist,][, ..regressors] %>% as.matrix()
  # Reduce by median
  dist_median <- colMedians(dist_data, na.rm=T)
  # Create table to hold site data
  dist_df <- data.table(t(data.table(dist_median)))
  names(dist_df) <- regressors
  dist_df$distance_km <- dist
  # Generate regression
  cluster_medians <- rbind(cluster_medians, dist_df)
}

# Load Clusters
load(paste0(wd_root, '/tmp_vars/clusters_calculated.RData'))
load(paste0(wd_root, '/tmp_vars/clusters_calculated2.RData'))

# DO YOU WANT TO USE RATING CURVE DERIVED CLUSTERS
clusters_calculated <- clusters_calculated2

# Assign cluster to each site
cluster_medians$cluster <- cl_predict(clusters_calculated, cluster_medians)
# Extract cluster for each dist
cluster_dist <- cluster_medians[,.(distance_km, cluster)]

transect_sentinel_harmonized <- left_join(transect_clean, cluster_dist, by=('distance_km'='distance_km'))

# Generate median B,G,R for each SSC category and each cluster or site
ssc_category_color <- transect_sentinel_harmonized  %>%group_by(cluster) %>% 
  summarise_at(vars("B4","B3","B2"), mean)
raster_color_types <- geom_raster(aes(fill = rgb(B4/2200,B3/2200,B2/2200))) # true color

ggplot(ssc_category_color, aes(x = cluster, y=cluster)) +
  raster_color_types +
  scale_fill_identity() +
  theme_classic()+
  # scale_x_continuous(expand_scale(add = c(0,0))) + 
  # scale_y_discrete(expand_scale(mult = c(0,0))) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs()
# width_drainage <- read.csv("D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/exports/GEE_raw/transect/chattahoochee_centerline_buffered_width.csv")
# transect_sentinel_harmonized <- left_join(transect_sentinel_harmonized, width_drainage, by=('distance_km'='distance_km'))
# Cloud and ice filter
# transect_sentinel_harmonized <- transect_sentinel_harmonized[
#   transect_sentinel_harmonized$thin_cirrus_percentage < 5 &
#   transect_sentinel_harmonized$snow_ice < 5,
# ]

# Bar plot of cluster
cluster_freq <- as.data.frame(table(transect_sentinel_harmonized$cluster))
colnames(cluster_freq) <- c('Cluster', 'Freq')
ggplot(cluster_freq, aes(Cluster, Freq)) +
  geom_bar(stat='identity')

# Set cluster to majority cluster (should be 1 clear winner)
transect_sentinel_harmonized$cluster <- 3

write.csv(transect_sentinel_harmonized, paste0(wd_gee, '/transect/chattahoochee_transect_20scale_ratingclust3.csv'), row.names = FALSE)
transect_sentinel_harmonized <- read.csv(paste0(wd_gee, '/transect/transect_harmonized_20scale.csv'))
