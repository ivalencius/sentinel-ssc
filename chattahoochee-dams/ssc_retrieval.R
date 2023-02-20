---
title: "SSC_Retrieval"
output: github_document
author: "Ilan Valencius"
date: "2022-06-16"
---

# Use: "For retrieving USGS suspended sediment concentration and discharge data along a river transect and integrating with landsat derived SSC measurements"

## Import Packages

```{r packages}
library(dataRetrieval)
library(sp)
library(sf)
library(maptools)
library(magrittr)
library(maps)
library(rgdal)
library(ggplot2)
library(usmap)
library(stringr)
library(dataRetrieval)
library(dplyr)
library(patchwork)
library(dams)
library(data.table)
library(nhdplusTools)
library(riverdist)
library(GISTools)
library(rgeos)
library(tidyr)
library(qpcR)
library(scales)
library(vioplot)
library(raster)
library(units)
library(extrafont)
loadfonts(device = "win")
library(latex2exp) # use TeX() to render latex
library(paletteer) # color palettes, 
```

## Set up imports

-   00060: discharge, cubic feet per second
-   80154: Suspended sediment concentration, milligrams per liter
-   Dates should be in in YYYY-MM-DD format

```{r imports}
wd_root <- 'D:/valencig/Thesis/USGS_dataquery'
centerline_kml <- "D:/valencig/Thesis/Data/chattahoochee_shape/chattahoochee_centerline.kml"

# Location of dams saved as kml files in dam_read_dir
# dam_names <- paste0(c("Buford Dam",
#                       "Morgan Falls",
#                       "Langdale Dam",
#                       "Eagle and Phoenix Dam",
#                       "West Point Dam",
#                       "Riverview Dam",
#                       "Bartletts Ferry Dam",
#                       "Goat Rock Dam",
#                       "Oliver Dam",
#                       "North Highlands Dam",
#                       "Walter F. George Dam",
#                       "George W Andrew Lock and Dam",
#                       "Jim Woodruff Dam"), ".kml")
# dam_read_dir <- "D:/valencig/Thesis/Data/chattahoochee_shape/dams/"

# When you want station info
startDate <- "1985-01-01"
endDate <- "2023-01-01"

## Query parameters for station data ##

# Daily value parameters
dv_parameters <- c("00060")
# Water quality parameters
wq_parameters <- c("80154") #(SSC, mgL)
# for WQP portal, characteristic name is hardcoded in
# Stat Codes
statCd <- c("00003") # stat code: mean

## If using station numbers from file ##

# 8 digit station numbers of interest, column name must be 'site_no'
station_file <- "D:/valencig/Thesis/USGS_dataquery/imports/station_numbers.csv"
usgs_prefix <- FALSE # do station numbers have 'USGS-' in title

## If getting station data from transect ##
upstream_distance_km <- 686

## Dam info ## 
use_dam_file <- TRUE # if set to false will autoquery data from NID website
dam_file <- 'D:/valencig/Thesis/USGS_dataquery/imports/nation.csv'

## Predicted SSC from landsat ##
landsat_csv <- 'D:/valencig/Thesis/chattahoochee-dams/chattahoochee-dams-exports/SSC_pred.csv'

## Don't touch this code ##

# helper function -> pads a character string with leading 0s (for USGS codes)
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}

centerline <- as.data.frame(maptools::getKMLcoordinates(centerline_kml, 
                                                        ignoreAltitude = TRUE))
names(centerline) <- c('lon','lat')
# Convert centerline to sf for plotting
centerline_sf <- st_as_sf(centerline, coords=c("lon","lat")) %>% st_set_crs(4326)

wd_root <- paste0(wd_root, '/', startDate, '_', endDate)
if(!dir.exists(wd_root)){
    dir.create(wd_root)
}
save_dir <- paste0(wd_root, '/variables/')
```

## Plot centerline transect

```{r transect}
# plot_dam <- function(dam_name){
#   dam <- as.data.frame(t(as.data.frame(maptools::getKMLcoordinates(paste0(dam_read_dir, dam_name),
#                                                         ignoreAltitude = TRUE))))
#   names(dam) <- c('lon','lat')
#   dam_transform <- usmap_transform(dam)
#   geom_point(data = dam_transform,
#              aes(x = x, y = y),
#              color = "green",
#              shape=2)
# }
# 
# centerline_transform <- usmap_transform(centerline)
# 
# # Create Dam Map
# dam_map = plot_usmap(fill = "black", alpha = 0.5, color = "white", size = 1,
#            include= c("GA","FL", "AL"), labels = TRUE) +
#   labs(title = "Chattahooche Centerline", subtitle = "Buford Dam to Brothers River") +
#   geom_path(data = centerline_transform,
#             aes(x = x, y = y),
#             color = "red",
#             size = 1) +
#   lapply(dam_names, plot_dam)
# 
# dam_map
# 
# # Save Map
# ggsave(dam_map,
#        filename = 'exports/rivertransect.png',
#        width = 5, height = 8)

```

## Get station info and plot station locations (from file)

```{r station plot}
# Function to remove empty columns of query
# rm_empty <- function(data){
#   empty_columns <- colSums(is.na(data) | data == "") == nrow(data)
#   return(data[, !empty_columns])
# }
# stations <- read.csv(station_file,colClasses=c("character")) # read as character to keep leading zeros
# if (usgs_prefix) {
#   # remove USGS-prefix
#   station_nums <- unlist(strsplit(stations$site_no, split='-'))[c(FALSE, TRUE)]
# } else {
#   station_nums <- stations$site_no
# }
# station_nums <- unlist(lapply(station_nums, pad0))
# station_data <- readNWISsite(station_nums)
# station_data <- rm_empty(station_data)

# To export information on stations to csv
# write.csv(station_data, "exports/valencius_chattahoochee_usgs_stations.csv", row.names=FALSE)
# 
# station_transform = usmap_transform(station_data, 
#                                       input_names = c("dec_long_va", "dec_lat_va"))
# 
# station_map <- plot_usmap(fill = "black", alpha = 0.5, color = "white", size = 1, 
#            include= c("GA","FL", "AL"), labels = TRUE) +
#   labs(title = "USGS Station Locations") +
#   geom_path(data = centerline_transform,
#             aes(x = x, y = y),
#             color = "red",
#             size = 1) +
#   geom_point(data = station_transform,
#               aes(x = x, y = y),
#               color = "yellow",
#              shape=4)
# 
# station_map

# Save Map
# ggsave(station_map,
#        filename = 'exports/valencius_stationlocs.png',
#        width = 5, height = 8)
```

## Query each station (from file)

```{r station info}
## Get daily discharge data ##
# dv <- readNWISdv(station_nums, dv_parameters, startDate, endDate, statCd = statCd)
# dv_filtered <- data.frame(dv["site_no"], 
#                           dv["Date"],
#                           dv["X_00060_00003"]) # X_00060_00003 will need to change depending on your stat codes
## Get water quality data ##
# wq <- readWQPqw(paste0("USGS-", station_nums), wq_parameters, startDate, endDate)

#stations <- read.table(station_file, header = TRUE, colClasses = "character")

## Get daily discharge data ##
# names(dv_filtered) <- c("site_no", "date", "discharge")

# Strip 'USGS-' identifier from station numbers
# wq["MonitoringLocationIdentifier"] <- lapply(list(wq[["MonitoringLocationIdentifier"]]), substring,
#                                                first = 6)
# wq_filtered <- data.frame(wq["MonitoringLocationIdentifier"],
#                             wq["ActivityStartDate"],
#                             wq["ResultMeasureValue"])
# names(wq_filtered) <- c("site_no","date", "sscmgL")
## To export station data to csv ##
# write.csv(dv_filtered, "exports/stations_from_file/valencius_discharge.csv", row.names=FALSE)
# write.csv(wq_filtered, "exports/stations_from_file/valencius_ssc.csv", row.names=FALSE)
```

## Automatically extract SSC and Q from stations

-   **Note:** this may take a while

```{r load data}
## Run if you have previously saved data from {r get basin} ##
load(file=paste0(save_dir, 'riverdata.Rdata'))
load(file=paste0(save_dir, 'nwis_stations.Rdata'))
load(file=paste0(save_dir, 'wqp_stations.Rdata'))
load(file=paste0(save_dir, 'SSC_full.Rdata'))
load(file=paste0(save_dir, 'Q_full.Rdata'))
```

```{r fix line2network fxn}
line2network <- function (sp = NA, path = ".", layer = NA, tolerance = 100, 
  reproject = NULL, supplyprojection = NULL) 
{
  if (suppressWarnings(is.na(sp))) {
    sp <- suppressWarnings(rgdal::readOGR(dsn = path, layer = layer, 
      verbose = F))
  }
  if (class(sp) != "SpatialLinesDataFrame") {
    stop("Specified shapefile is not a linear feature.")
  }
  if (is.na(sp@proj4string@projargs) & !is.null(supplyprojection)) {
    sp@proj4string@projargs <- supplyprojection
  }
  if (is.na(sp@proj4string@projargs)) {
    stop("Shapefile projection information is missing.  Use supplyprojection= to specify a Proj.4 projection to use.  If the input shapefile is in WGS84 geographic (long-lat) coordinates, this will be +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 (in double-quotes).  If so, it must also be reprojected using reproject=.")
  }
  proj4 <- strsplit(sp@proj4string@projargs, split = " ")
  projected <- sp::is.projected(sp)
  if (is.null(reproject) & !projected) 
    stop("Distances can only be computed from a projected coordinate system.  Use reproject= to specify a Proj.4 projection to use.")
  if (!is.null(reproject)) {
    sp <- sp::spTransform(sp, sp::CRS(reproject))
    proj4 <- strsplit(sp@proj4string@projargs, split = " ")
  }
  units <- "unknown"
  for (i in 1:length(proj4[[1]])) {
    if (proj4[[1]][i] != "") {
      proj4arg <- strsplit(proj4[[1]][i], split = "=")
      if (proj4arg[[1]][1] == "+units") {
        units <- proj4arg[[1]][2]
        cat("\n", "Units:", proj4arg[[1]][2], "\n")
      }
    }
  }
  if (length(sp@lines) > 1) {
    sp_line <- NA
    sp_seg <- NA
    lines <- list()
    j <- 1
    for (i in 1:length(sp@lines)) {
      for (k in 1:length(sp@lines[i][[1]]@Lines)) {
        lines[[j]] <- sp@lines[i][[1]]@Lines[[k]]@coords
        sp_line[j] <- i
        sp_seg[j] <- k
        j <- j + 1
      }
    }
  }
  if (length(sp@lines) == 1) {
    lines <- sp@lines[1][[1]]@Lines
    length <- length(lines)
    lines.new <- list()
    for (i in 1:length) {
      lines.new[[i]] <- lines[[i]]@coords
    }
    lines <- lines.new
    sp_line <- rep(1, length)
    sp_seg <- 1:length
  }
  length <- length(lines)
  rivID <- 1:length
  lineID <- data.frame(rivID, sp_line, sp_seg)
  connections <- calculateconnections(lines = lines, tolerance = tolerance)
  if (any(connections %in% 5:6)) 
    braided <- TRUE
  lengths <- rep(NA, length)
  for (i in 1:length) {
    lengths[i] <- pdisttot(lines[[i]])
  }
  names <- rep(NA, length)
  mouth.seg <- NA
  mouth.vert <- NA
  mouth <- list(mouth.seg, mouth.vert)
  names(mouth) <- c("mouth.seg", "mouth.vert")
  sequenced <- FALSE
  braided <- NA
  cumuldist <- list()
  for (i in 1:length) {
    xy <- lines[[i]]
    n <- dim(xy)[1]
    cumuldist[[i]] <- c(0, cumsum(sqrt(((xy[1:(n - 1), 1] - 
      xy[2:n, 1])^2) + ((xy[1:(n - 1), 2] - xy[2:n, 2])^2))))
  }
  out.names <- c("sp", "lineID", "lines", "connections", "lengths", 
    "names", "mouth", "sequenced", "tolerance", "units", 
    "braided", "cumuldist")
  out <- list(sp, lineID, lines, connections, lengths, names, 
    mouth, sequenced, tolerance, units, braided, cumuldist)
  names(out) <- out.names
  class(out) <- "rivernetwork"
  length1 <- length(out$lengths)
  length2 <- length(out$lengths)
  if (length2 < length1) 
    cat("\n", "Removed", length1 - length2, "duplicate segments.", 
      "\n")
  length3 <- length(out$lengths)
  if (length3 < length2) 
    cat("\n", "Removed", length2 - length3, "segments with lengths shorter than the connectivity tolerance.", 
      "\n")
  return(out)
}
```

```{r station extract}
save_dir <- paste0(wd_root,'/variables/')
if(!dir.exists(save_dir)){
    dir.create(save_dir)
}
if(!dir.exists(paste0(wd_root,'/exports/'))){
    dir.create(paste0(wd_root,'/exports/'))
}
## Determine Watershed Geometries ##
summarize.nldi = function(input){
  data.frame(name = names(input), 
             class = sapply(input, class)[1], 
             row.names = NULL) %>% 
    mutate(feature_count = ifelse(class == "sf", sapply(input, nrow), 
                                  sapply(input, length)))
}
# Get last lat, long pair from centerline
start_loc <- tail(centerline, 1)
start_loc <- c(start_loc[[1]], start_loc[[2]]) 
riverdata <- findNLDI(location = start_loc,
              nav = c("UT", "UM"),
              find = c("nwis", "WQP","basin", "flowlines"),
              distance_km = upstream_distance_km)

summarize.nldi(riverdata)

## Merge Site Names ##
# get all WQP sites
mc_names <- paste0("USGS-",lapply(riverdata$UM_WQP$comid, pad0))
ut_names <- paste0("USGS-",lapply(riverdata$UT_WQP$comid, pad0))
# merge with nwis sites
mc_names <- unique(append(mc_names, riverdata$UM_nwissite$identifier))
ut_names <- unique(append(ut_names, riverdata$UT_nwissite$identifier))

# Remove mislabelled sites outside of the watershed (remove from ut and mc) 
mc_names <- mc_names[- match("USGS-02297310",mc_names)]
mc_names <- mc_names[- match("USGS-02297272",mc_names)]
mc_names <- mc_names[- match("USGS-03298550",mc_names)]
ut_names <- ut_names[- match("USGS-02297310",ut_names)]
ut_names <- ut_names[- match("USGS-02297272",ut_names)]
ut_names <- ut_names[- match("USGS-03298550",ut_names)]
ut_names <- ut_names[- match("USGS-02296750",ut_names)]
ut_names <- ut_names[- match("USGS-02047783",ut_names)]
ut_names <- ut_names[- match("USGS-03435970",ut_names)]
ut_names <- ut_names[- match("USGS-03291310",ut_names)]
ut_names <- ut_names[- match("USGS-03437900",ut_names)]
ut_names <- ut_names[- match("USGS-02130840",ut_names)]
ut_names <- ut_names[- match("USGS-02308865",ut_names)]
ut_names <- ut_names[- match("USGS-02349900",ut_names)]
ut_names <- ut_names[- match("USGS-02350080",ut_names)]

## Get Dicharge Data ##riverdata$UM_nwissite$identifier
Q_mc <- data.table(readNWISdv(gsub("USGS-", "", mc_names), 
                parameterCd = dv_parameters, 
                startDate = startDate,
                endDate = endDate) %>% renameNWISColumns())[,":=" (
                  site_no = paste0("USGS-", site_no),
                  date = Date,
                  flow = Flow,
                  portal = 'NWIS',
                  channel = 'MAIN STEM'
                )][,.(site_no, date, flow, portal, channel)]
Q_ut <- data.table(readNWISdv(gsub("USGS-", "", ut_names), 
                parameterCd = dv_parameters, 
                startDate = startDate,
                endDate = endDate) %>% renameNWISColumns())[,":=" (
                  site_no = paste0("USGS-", site_no),
                  date = Date,
                  flow = Flow,
                  portal = 'NWIS',
                  channel = 'TRIBUTARY'
                )][,.(site_no, date, flow, portal, channel)]

# Remove main channel sites from ut sites
Q_ut <- anti_join(Q_ut, Q_mc, by="site_no")

## Get SSC data ##
SSC_wqp_mc <- data.table(readWQPdata(siteNumbers=mc_names, 
                parameterCd = wq_parameters,
                startDate=startDate,
                endDate=endDate
                ))[,":=" (
                  site_no = MonitoringLocationIdentifier,
                  date = ActivityStartDate,
                  sscmg_L = ResultMeasureValue,
                  portal = 'WQP',
                  channel = 'MAIN STEM'
                )][,.(site_no, date, sscmg_L, portal, channel)]

# May need to save as a few different variables (ex foo1, foo2, foo3) and rbind() 
# because WQP can't handle large requests (limit to 1000 site numbers at a time)
SSC_wqp_ut <- data.table(readWQPdata(siteNumbers=ut_names, 
                parameterCd = wq_parameters,
                startDate = startDate,
                endDate = endDate
                ))[,":=" (
                  site_no = MonitoringLocationIdentifier,
                  date = ActivityStartDate,
                  sscmg_L = ResultMeasureValue,
                  portal = 'WQP',
                  channel = 'TRIBUTARY'
                )][,.(site_no, date, sscmg_L, portal, channel)]

# Remove main channel sites from ut sites
SSC_wqp_ut <- anti_join(SSC_wqp_ut, SSC_wqp_mc, by="site_no")

## Combine upper tributaries and main channel ##
Q_full <- rbind(Q_mc, Q_ut)
SSC_full <- rbind(SSC_wqp_mc, SSC_wqp_ut)

## Identify Stations of Interest ##
# NWIS sites
mc_nwis_names <- unique(Q_mc$site_no)
ut_nwis_names <- unique(Q_ut$site_no)

mc_nwis <- data.table(whatNWISsites(site=gsub("USGS-", "",mc_nwis_names)))[,":="(
    portal = "NWIS",
    site_no = paste0("USGS-", site_no),
    lat = dec_lat_va,
    lon = dec_long_va,
    channel="MAIN STEM"
  )][,.(portal, site_no, station_nm,site_tp_cd, lat, lon, channel)]
ut_nwis <- data.table(whatNWISsites(site=gsub("USGS-", "",ut_nwis_names)))[,":="(
    portal = "NWIS",
    site_no = paste0("USGS-", site_no),
    lat = dec_lat_va,
    lon = dec_long_va,
    channel="TRIBUTARY"
  )][,.(portal, site_no, station_nm,site_tp_cd, lat, lon, channel)]

# WQP sites
mc_wqp_names <- unique(SSC_wqp_mc$site_no)
ut_wqp_names <- unique(SSC_wqp_ut$site_no)

mc_wqp <- data.table(whatWQPsites(siteid=mc_wqp_names))[,":="(
    portal = "WQP",
    site_no = MonitoringLocationIdentifier,
    station_nm = MonitoringLocationName,
    site_code = MonitoringLocationTypeName,
    drainage_area_km2 = DrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
    lat = LatitudeMeasure,
    lon = LongitudeMeasure,
    elevation_m = VerticalMeasure.MeasureValue *0.3048, # ft -> m
    channel="MAIN STEM"
  )][,.(portal, site_no, station_nm, site_code, lat, lon, drainage_area_km2, elevation_m, channel)]

ut_wqp <- data.table(whatWQPsites(siteid=ut_wqp_names))[,":="(
    portal = "WQP",
    site_no = MonitoringLocationIdentifier,
    station_nm = MonitoringLocationName,
    site_code = MonitoringLocationTypeName,
    drainage_area_km2 = DrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
    lat = LatitudeMeasure,
    lon = LongitudeMeasure,
    elevation_m = VerticalMeasure.MeasureValue *0.3048, # ft -> m
    channel="TRIBUTARY"
  )][,.(portal, site_no, station_nm, site_code, lat, lon, drainage_area_km2, elevation_m, channel)]  

## Map gauge location to river distance ##
# trace(line2network, edit=T) --> comment out following two lines
# suppressMessages(out <- removemicrosegs(out))
# suppressMessages(out <- removeduplicates(out))
S1 <- as(riverdata$UM_flowlines, 'Spatial')
nldi_proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
AKalbers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154
    +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
river_center <- line2network(S1, reproject=AKalbers)

# Apply projection to mc nwis sites then convert into meter projection
mc_nwis_spatial <- mc_nwis
mc_nwis_spatial$lat <- as.numeric(mc_nwis_spatial$lat)
mc_nwis_spatial$lon <- as.numeric(mc_nwis_spatial$lon)
coordinates(mc_nwis_spatial) <- ~lon + lat
proj4string(mc_nwis_spatial) <- '+init=epsg:4326'
mc_nwis_spatial <- spTransform(mc_nwis_spatial, CRS(AKalbers))
mc_nwis_spatial <- as.data.frame(mc_nwis_spatial)
nwis_riv <- xy2segvert(x=mc_nwis_spatial$lon, y=mc_nwis_spatial$lat, rivers=river_center)

# Do the same for mc wqp sites
mc_wqp_spatial <- mc_wqp
mc_wqp_spatial$lat <- as.numeric(mc_wqp_spatial$lat)
mc_wqp_spatial$lon <- as.numeric(mc_wqp_spatial$lon)
coordinates(mc_wqp_spatial) <- ~lon + lat
proj4string(mc_wqp_spatial) <- '+init=epsg:4326'
mc_wqp_spatial <- spTransform(mc_wqp_spatial, CRS(AKalbers))
mc_wqp_spatial <- as.data.frame(mc_wqp_spatial)
wqp_riv <- xy2segvert(x=mc_wqp_spatial$lon, y=mc_wqp_spatial$lat, rivers=river_center)

# Determine distance from start of channel for nwis sites
distances_nwis <- list()
for (i in 1:length(nwis_riv$vert)){
  # Distance to mouth of river
  dist_m <- riverdistance(startseg = nwis_riv$seg[i], startvert = nwis_riv$vert[i],
                            endvert = 1, endseg = 1, rivers=river_center)
  # Convert from degrees to km
  distances_nwis <- append(distances_nwis, upstream_distance_km-dist_m/1000)
}
nwis_riv$dist_downstream_km <- distances_nwis

# Determine distance from start of channel for wqp sites
distances_wqp <- list()
for (i in 1:length(wqp_riv$vert)){
  # Distance to mouth of river
  dist_m <- riverdistance(startseg = wqp_riv$seg[i], startvert = wqp_riv$vert[i],
                            endvert = 1, endseg = 1, rivers=river_center)
  # Convert from degrees to km
  distances_wqp <- append(distances_wqp, upstream_distance_km-dist_m/1000)
}
wqp_riv$dist_downstream_km <- distances_wqp

## Merge station info ##
mc_nwis$dist2river_m <- nwis_riv$snapdist
mc_nwis$dist_downstream_km <- nwis_riv$dist_downstream_km
mc_wqp$dist2river_m <- wqp_riv$snapdist
mc_wqp$dist_downstream_km <- wqp_riv$dist_downstream_km

nwis_stations <- rbind(mc_nwis, ut_nwis, fill=TRUE)
wqp_stations <- rbind(mc_wqp, ut_wqp, fill=TRUE)

## Map station distances to Q and SSC measurements ##
Q_full <- left_join(Q_full, nwis_stations)
SSC_full <- left_join(SSC_full, wqp_stations)

## Convert variable data tables to upper case ##
nwis_stations <- mutate_all(nwis_stations, toupper)
wqp_stations <- mutate_all(wqp_stations, toupper)
Q_full <- mutate_all(Q_full, toupper)
SSC_full <- mutate_all(SSC_full, toupper)

## Save variables to save time ##
save(riverdata, file=paste0(save_dir, 'riverdata.Rdata'))
save(nwis_stations, file=paste0(save_dir, 'nwis_stations.Rdata'))
save(wqp_stations, file=paste0(save_dir, 'wqp_stations.Rdata'))
save(SSC_full, file=paste0(save_dir, 'SSC_full.Rdata'))
save(Q_full, file=paste0(save_dir, 'Q_full.Rdata'))

## Export data to csv ##
write.csv(Q_full, paste0(wd_root,"/exports/Q_full.csv"), row.names=FALSE)
write.csv(SSC_full, paste0(wd_root,"/exports/SSC_full.csv"), row.names=FALSE)
write.csv(nwis_stations, paste0(wd_root, '/exports/Q_stations.csv'), row.names=FALSE)
write.csv(wqp_stations, paste0(wd_root, '/exports/SSC_stations.csv'), row.names=FALSE)
```

## Determine locations of all dams in the watershed

-   Dam data cane be queried automatically or downloaded in csv form from the NID Website [<https://nid.usace.army.mil/#/downloads>]
-   **Note:** If you download the nation.csv which contains all NID data, you must remove the first row which contains the last date the file was modified (R uses the first row as column names)

```{r dam query}
## Issue --> not all dams stored in NID database
## If importing from file ##
if (use_dam_file) {
  nid_dams <- read.csv(dam_file)
} else {
  ## CURRENTLY NOT WORKING
  nid_dams <- get_nid()
}

# Remove entries which are missing lat, lon
nid_dams <- na.omit(nid_dams, cols=c('Longitude','Latitude'))
# Convert to sf object
nid_dams_sf <- st_as_sf(nid_dams, coords = c("Longitude", "Latitude")) %>% st_set_crs(4326)

# Get dams in the basin
basin_dams <- st_intersection(nid_dams_sf,riverdata$basin$geometry)
```

# Start Making Plots

## Compare Centerline KML and NLDI

```{r NLDI sanity check}
centerline_compare <-
  ggplot() + 
    geom_sf(data = riverdata$basin) + 
    geom_sf(data = centerline_sf, col = "red", alpha = .5, size = 1)+
    theme_void() +
    labs(caption = "Input extent") +
  ggplot() +
    geom_sf(data = riverdata$basin, col = NA) + 
    geom_sf(data = riverdata$UM_flowlines, col = "blue",size = 1) +
    theme_void() +
    labs(title = "Centerlines",
         caption = "NLDI extent") + 
    theme(legend.position = "none",
           plot.title = element_text(face = "bold", hjust = .5))

centerline_compare

# Save Map
ggsave(centerline_compare, 
       filename = paste0(wd_root, '/exports/centerline_compare.png'),
       width = 5, height = 8)
```

## Plot all NWIS stations (time independent)

```{r NWIS stations}
## Plot Station Data ##
nwis_station_plot <- ggplot() + 
  ## Basin and Channels ##
  geom_sf(data = riverdata$basin) + 
  geom_sf(data = riverdata$UM_flowlines, col = "blue") + 
  geom_sf(data = riverdata$UT_flowlines, col = "blue", alpha = .5, lwd=0.25) + 
  theme_bw() +
  coord_sf(label_graticule = 'NESW') +
  
  geom_sf(data = riverdata$UM_nwissite, shape = 21, 
          fill = NA, stroke = 0.75, size = 1.5, color="red") + # Shape aes: USGS gage, set with scale_shape_manual
  geom_sf(data = riverdata$UT_nwissite, shape = 21, 
          fill = NA, stroke = 0.25, size = 0.5, color="red") + # Shape aes: USGS gage, set with scale_shape_manual 
  scale_size_continuous(range = c(1,7)) +
  scale_color_brewer(palette = "Set1") +
  #theme_void() +
  labs(title = "USGS Station Locations") + 
  theme(legend.position = "none",
         plot.title = element_text(face = "bold", hjust = .5))

nwis_station_plot

# Save Map
ggsave(nwis_station_plot,
       filename = paste0(wd_root, '/exports/usgs_stations.png'),
       width = 5, height = 8)
```

## Plot locations of USGS stations which measure discharge (in time period)

```{r discharge stations}
## Plot Station Data ##
nwis_plot <- nwis_stations%>%st_as_sf(coords=c('lon','lat')) %>% st_set_crs(4326)

discharge_locs <- ggplot() + 
  ## Basin and Channels ##
  geom_sf(data = riverdata$basin, show.legend = F) + 
  geom_sf(data = riverdata$UM_flowlines, col = "blue", show.legend = F) + 
  geom_sf(data = riverdata$UT_flowlines, col = "blue", alpha = .5, lwd=0.25, show.legend = F) + 
  theme_bw() +
  scale_size_continuous(range = c(1,7)) +
  
  ## Upper main channel#
  geom_sf(data = nwis_plot, aes(color=channel, shape=channel, size=channel),  fill = NA, stroke = 0.75) + 
  #scale_size_continuous(range = c(1,7)) +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = c('green', 'red')) +
  scale_shape_manual(values = c(21, 1)) +
  scale_size_manual(values = c(3, 1)) +
  #theme_void() +
  coord_sf(label_graticule = 'NESW') +
  labs(title = "USGS Station Locations Measuring Discharge",
       caption = paste0(startDate, ' to ', endDate))
  #theme(legend.position = "none",
  #       plot.title = element_text(face = "bold", hjust = .5))

discharge_locs

# Save Map
ggsave(discharge_locs, 
       filename = paste0(wd_root, '/exports/discharge_stations.png'),
       width = 5, height = 8)

```

## Plot locations of SSC sampling (in time period)

```{r ssc stations}
## Plot Station Data ##
wqp_plot<-wqp_stations%>%st_as_sf(coords=c('lon','lat')) %>% st_set_crs(4326)

ssc_locs <- ggplot() + 
  ## Basin and Channels ##
  geom_sf(data = riverdata$basin, show.legend = F) + 
  geom_sf(data = riverdata$UM_flowlines, col = "blue", show.legend = F) + 
  geom_sf(data = riverdata$UT_flowlines, col = "blue", alpha = .5, lwd=0.25, show.legend = F) + 
  theme_bw() +
  scale_size_continuous(range = c(1,7)) +
  
  ## Upper main channel#
  geom_sf(data = wqp_plot, aes(color=channel, shape=channel, size=channel),  fill = NA, stroke = 0.75) + 
  #scale_size_continuous(range = c(1,7)) +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = c('green', 'red')) +
  scale_shape_manual(values = c(21, 1)) +
  scale_size_manual(values = c(3, 1)) +
  #theme_void() +
  coord_sf(label_graticule = 'NESW') +
  labs(title = "USGS Station Locations Measuring SSC",
       caption = paste0(startDate, ' to ', endDate))
  #theme(legend.position = "none",
  #       plot.title = element_text(face = "bold", hjust = .5))
ssc_locs

# Save Map
ggsave(ssc_locs, 
       filename = paste0(wd_root, '/exports/ssc_sites.png'),
       width = 5, height = 8)

```

## Plot dam locations and flowlines

```{r dam locations}
dam_flowline_map <- ggplot() + 
  ## Basin and Channels ##
  geom_sf(data = riverdata$basin) + 
  geom_sf(data = riverdata$UM_flowlines, col = "blue") + 
  geom_sf(data = riverdata$UT_flowlines, col = "blue", alpha = .5, lwd=0.25) + 
  theme_bw() +
  coord_sf(label_graticule = 'NESW') +
  geom_sf(data = basin_dams, 
          # Set size by storage capacity
          aes(size = Volume..Cubic.Yards.),
          color = 'red',
          shape = 2,
          fill = NA, stroke = 1) +# No fill, stroke (same as lwd for lines) is equal to 1 
  scale_size_continuous(range = c(0.5,3)) +
  labs(title = "NID Dam Locations")

dam_flowline_map

# Save Map
ggsave(dam_flowline_map, 
       filename = paste0(wd_root, '/exports/dam_flowline_map.png'),
       width = 5, height = 8)
```

## Import landsata/Sentinel-2 data and harmonize with USGS data

```{r landsat import}
landsat_data <- data.table(read.csv(landsat_csv))[,':='(
  lat=Latitude,
  lon=Longitude,
  date=sample_dt
)][,.(SSC_mgL, date, year, distance_km, lat, lon)]

landsat_data$year_chunk <- cut(as.numeric(landsat_data$year), seq(from = 1985, to = 2026, by = 5), right=FALSE)
landsat_data <- landsat_data %>% mutate(SSC_mgL=replace(SSC_mgL, SSC_mgL>=10000, NA))

## Filter data to time of interest stamp ##
# Remember dates are in YYYY-MM-DD format
landsat_data_filter <- landsat_data %>% 
  filter(as.Date(date)>=as.Date(startDate)) %>%
  filter(as.Date(date)<=as.Date(endDate))
  # # Filter by year
  # filter(year<=as.numeric(substr(endDate, 1, 4))) %>%
  # filter(year>=as.numeric(substr(startDate, 1, 4))) %>%
  # # Filter by month for end years
  # filter(!(year==as.numeric(substr(startDate, 1, 4)) & month<as.numeric(substr(startDate, 6, 7)))) %>%
  # filter(!(year==as.numeric(substr(endDate, 1, 4)) & month>as.numeric(substr(endDate, 6, 7))))

### THEMES- From Evan Dethier ###
theme_facet <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    # panel.grid = element_blank(),
    # legend.position = 'none',
    panel.border = element_rect(size = 0.5),
    strip.background = element_rect(fill = 'white'),
    text = element_text(size=12),
    axis.text = element_text(size = 12), 
    plot.title = element_text(size = 13)
  )

# season_facet <- theme_facet + theme(
#   #legend.position = 'none', 
#   strip.background = element_blank(),
#   strip.text = element_text(hjust = 0, margin = margin(0,0,0,0, unit = 'pt'))
# )

lm_eqn <- function(m){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 4),
                        b = format(unname(coef(m)[2]), digits = 4),
                        r2 = format(summary(m)$r.squared, digits = 4)))
  as.character(as.expression(eq));
}
```

# DEALING WITH SENTINEL DATA

```{r sentinel import}
# Import data trained on USGS In-Situ data
# sentinel_data <- data.table(
#   read.csv("D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_scaling\\train_USGS_clust_USGS\\regression\\lasso_full_ratio_clust6.csv"))[,':='(
#   SSC_mgL = pred_SSC_mgL,
#   year = format(as.Date(date), "%Y")
# )][,.(SSC_mgL, date, year, distance_km, lat, lon, cluster)]

# Import data trained on USGS In-Situ Data + Rating curve derived SSC
sentinel_dataRating <- data.table(
  read.csv("D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_scaling\\train_ratingCurve_clust_ratingCurve\\regression\\lasso_full_bands_ratingclust3.csv"))[,':='(
  SSC_mgL_R = pred_SSC_mgL, # SSC (mg/L) Rating curve --> data trained and clusters using rating curve derived SSC
  year = format(as.Date(date), "%Y")
)][,.(SSC_mgL_R, date, year, distance_km, lat, lon, cluster)]

# Join In-Situ and rating curve derived
# sentinel_data <- left_join(sentinel_data, sentinel_dataRating)
# Because columns area already aligned (from same csv)
sentinel_data$SSC_mgL_R <- sentinel_dataRating$SSC_mgL_R

# Filter out SSC measurements above 10*4.5
sentinel_data <- sentinel_data %>% filter(SSC_mgL <31622, SSC_mgL_R<31662)
# sentinel_data <- sentinel_data %>% filter(SSC_mgL < max(landsat_data_filter$SSC_mgL, na.rm=T))
# Add year chunk
sentinel_data$year_chunk <- cut(as.numeric(sentinel_data$year), seq(from = 1985, to = 2026, by = 5), right=FALSE)

# ggplot(landsat_data_filter, aes(x=SSC_mgL)) + geom_density() +labs(title="Landsat")
# ggplot(sentinel_data, aes(x=SSC_mgL)) + geom_density() +labs(title="Sentinel")
# 
# # Try splitting into high and low data
# # s_high <- sentinel_data %>% filter(SSC_mgL > 50)
# # s_low <- sentinel_data %>% filter(SSC_mgL <= 50)
# # foo <- sentinel_data
# # foo <- s_low %>% group_by(distance_km) %>% tally()
# # foo2 <- s_low %>% group_by(distance_km) %>% summarize_at("SSC_mgL",mean)
# # foo2$SSC_mgL <- as.integer(foo2$SSC_mgL)
# # foo$mov_avg <- rollmean(foo$SSC_mgL, 15, na.rm=T,na.pad = TRUE, align = "right")
# # Need to go through by date and eliminate outliers for every day
# # Example look at 44 km and 2020-12-20, two wonky observations at that day
# remove_outliers <- function(x) {
#   qnt <- quantile(x$SSC_mgL, probs=c(.25, .75), na.rm = T)
#   # # H <- 1.5 * IQR(x$SSC_mgL, na.rm = T)
#   y <- x
#   # y[x$SSC_mgL < (qnt[1]),]$SSC_mgL <- NA
#   y[x$SSC_mgL > (qnt[2]),]$SSC_mgL <- NA
#   # y <- x
#   # y$SSC_mgL <- mean(x$SSC_mgL)
#   return(y)
# }
# 
# sentinel_data_filter <- sentinel_data %>% group_by(date) %>%  
#   group_modify(~remove_outliers(.x)) %>%
#   filter(!is.na(SSC_mgL))
#   # mutate(new = case_when(action == 0 ~ remove_outliers(stuff), TRUE ~ stuff))
# 
# s_high <- sentinel_data_filter %>% group_by(date) %>% summarize_at("SSC_mgL", mean) %>% filter(SSC_mgL > 50)
# s_low <- sentinel_data_filter  %>% group_by(date) %>% summarize_at("SSC_mgL", mean) %>% filter(SSC_mgL <= 50)

ggplot(sentinel_data, aes(x = distance_km, y = SSC_mgL, color="in-situ")) +
    # geom_point() +
    stat_summary(geom = 'line', fun = 'mean', na.rm=T) +
    stat_summary(aes(y=SSC_mgL_R, color='using rating curve'),geom = 'line', fun = 'mean', na.rm=T) +
    # coord_cartesian(ylim=c(0,60)) +
    # geom_line(aes(color='Red', y=frollmean(SSC_mgL, 7, na.rm=T))) +
    # season_facet +
    theme_minimal()+
    labs(title = 'SSC Concentration - Chattahoochee River',
         caption = paste0(startDate, ' to ', endDate),
         x = 'Distance Downstream (km)', 
         y = 'SSC (mg/L)')

# Plot a random day
# dates <- sample(sentinel_data$date, 5)
# for (date in dates){
#   ggplot(sentinel_data %>% filter(date == date), aes(x = distance_km, y = SSC_mgL)) +
#     # geom_point() +
#     stat_summary(geom = 'point', fun = 'mean', na.rm=T) +
#     # coord_cartesian(ylim=c(0,60)) +
#     # geom_line(aes(color='Red', y=frollmean(SSC_mgL, 7, na.rm=T))) +
#     # season_facet +
#     theme_minimal()+
#     labs(title = 'SSC Concentration - Chattahoochee River',
#          caption = unique(date),
#          x = 'Distance Downstream (km)', 
#          y = 'SSC (mg/L)')
# }
```

# Start making plots

```{r}
SSC_full$dist_downstream_km <- as.numeric(SSC_full$dist_downstream_km)
SSC_full$sscmg_L <- as.numeric(SSC_full$sscmg_L)

## Create Plots ##
ssc_landsat_plot <- ggplot(landsat_data_filter, aes(x = distance_km, y = SSC_mgL)) +
    stat_summary(geom = 'line', fun = 'mean') +
    season_facet +
    labs(title = 'SSC Concentration - Chattahoochee River',
         caption = paste0(startDate, ' to ', endDate),
         x = 'Distance Downstream (km)', 
         y = 'SSC (mg/L)')
ssc_landsat_plot

ssc_sentinel_plot <- ggplot(sentinel_data, aes(x = distance_km, y = SSC_mgL)) +
    stat_summary(geom = 'line', fun = 'mean') +
    season_facet +
    labs(title = 'SSC Concentration - Chattahoochee River',
         caption = paste0(startDate, ' to ', endDate),
         x = 'Distance Downstream (km)', 
         y = 'SSC (mg/L)')
ssc_sentinel_plot

ssc_usgs_plot <- ggplot(data=drop_na(SSC_full), aes(x = dist_downstream_km, y = sscmg_L)) +
    geom_point() +
    stat_summary(
      geom = 'point',
      fun= 'mean',
      col = 'red',
      size = 3,
      shape = 24,
      fill = 'red',
      na.rm=T) +
    season_facet +
    labs(title = 'SSC Concentration - Chattahoochee River',
         caption = paste0(startDate, ' to ', endDate),
         x = 'Distance Downstream (km)', 
         y = 'SSC (mg/L)')
ssc_usgs_plot

# Combine data into one dataframe for one plot
# plot_data <- landsat_data_filter
# usgs_dist <- drop_na(SSC_full)$dist_downstream_km
# usgs_sscmg_L <- drop_na(SSC_full)$sscmg_L
# plot_data <- qpcR:::cbind.na(plot_data, usgs_dist)
# plot_data <- qpcR:::cbind.na(plot_data, usgs_sscmg_L)
# plot_data$year_chunk <- cut(as.numeric(plot_data$year), seq(from = 1985, to = 2026, by = 5), right=FALSE)

SSC_main <- SSC_full %>% filter(channel == 'MAIN STEM')
SSC_main$year_chunk <- cut(as.numeric(format(as.Date(SSC_main$date), "%Y")), seq(from = 1985, to = 2026, by = 5), right=FALSE)
# Get mean usgs ssc for each year
# Need to get mean for each distance then across the whole transect
ssc_mean <- data.table(year_chunk = unique(SSC_main$year_chunk), MN = as.numeric()) %>% arrange(year_chunk)
ssc_mean$false_int <- c(1,2,3,4,5,6,7,8)

for (year_subset in unique(landsat_data_filter$year_chunk)){
  year_data <- landsat_data_filter %>% filter(year_chunk == year_subset)
  x_means <-  year_data %>% group_by(distance_km) %>% summarize(x_ssc = mean(SSC_mgL, na.rm=T))
  ssc_mean[ssc_mean$year_chunk==year_subset,]$MN <- mean(x_means$x_ssc)
}

# plot_data
ssc_combined <- 
    ggplot(landsat_data_filter, aes(x = distance_km, y = SSC_mgL)) +
      stat_summary(
        aes(x = distance_km, y = SSC_mgL, color='Landsat'),
        geom = 'line', 
        fun = 'mean',
        inherit.aes = FALSE,
        show.legend=TRUE,
        na.rm=T) +
    stat_summary(
        data = sentinel_data,
        aes(x = distance_km, y = SSC_mgL_R, color='Sentinel'),
        geom = 'line', 
        fun = 'mean',
        inherit.aes = FALSE,
        show.legend=TRUE,
        na.rm=T) +
      #geom_point(data = drop_na(plot_data), aes(x = usgs_dist, y = usgs_sscmg_L)) +
        stat_summary(
          data = SSC_main,
          aes(x = dist_downstream_km, y = sscmg_L, color='1 standard Error'),
          geom = 'errorbar',
          fun.data= mean_se,
          fun.args = list(mult = 1), # mult = # of standard errors
          size = 0.5,
          width =  30,
          #color = 'green',
          inherit.aes = FALSE,
          show.legend=TRUE,
          na.rm=T) +
        stat_summary(
          data = SSC_main,
          aes(x = dist_downstream_km, y = sscmg_L, color='USGS'),
          geom = 'point',
          fun= 'mean',
          #size = 3,
          #shape=24,
          inherit.aes = FALSE,
          show.legend=TRUE,
          fill = 'purple',
          na.rm=T) +
        stat_summary(
          data = SSC_main,
          aes(x = dist_downstream_km, y = sscmg_L, color='USGS'),
          geom = 'line',
          fun= 'mean',
          #size = 3,
          #shape=24,
          inherit.aes = FALSE,
          show.legend=TRUE,
          na.rm=T) +
      # Mean value of USGS ssc
      # geom_hline(data=ssc_mean, aes(yintercept = MN, color='Mean Landsat'), linetype="dotted") +
      theme_minimal() +
      # facet_wrap(~year_chunk, scales = "free_y") +
      facet_wrap(~year_chunk) +
      coord_cartesian(ylim=c(0,150)) +
      theme(
        legend.position = "bottom"
      ) +
      guides(color=guide_legend("Key")) +
      scale_color_manual(values = c('green', '#999999','#E69F00','#9CDB79')) +
      labs(title = 'SSC Concentration - Chattahoochee River',
           caption = paste0(startDate, ' to ', endDate),
           x = 'Distance Downstream (km)', 
           y = 'SSC (mg/L)')

ssc_combined

# Plot mean SSC over time
lm_ssc <- function(df){
    m <- lm(MN ~ false_int, df);
    eq <- substitute(italic('Mean SSC') == b %.% italic('(Years since 1985)/5')+a*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

mean_ssc_plot <- ggplot(data = ssc_mean, aes(x=year_chunk, y=MN)) +
  geom_line(group=1) +
  geom_smooth(mapping=aes(x=false_int, y=MN), method = "lm", se=FALSE, color="red", formula = y ~ x) +
  geom_label(x = 6, y = 55, label = lm_ssc(ssc_mean), parse = T) +
  season_facet +
  labs(title='Landsat Derived Mean SSC Along Chattahoochee Transect',
       x='Year Grouping',
       y='Average SSC (mg/L)')

mean_ssc_plot

## Save combined plot ##
ggsave(ssc_combined, 
       filename = paste0(wd_root, '/exports/ssc_comparisons.png'),
       width = 12, height = 6)

## Save mean_ssc plot ##
ggsave(mean_ssc_plot, 
       filename = paste0(wd_root, '/exports/landsat_mean_ssc.png'),
       width = 12, height = 6)

```

## Switch to using Rating curve derived Sentinel Data
Delete USGS column and keep rating curve column

```{r}
sentinel_data <- sentinel_data %>% select(-c("SSC_mgL"))
colnames(sentinel_data)[colnames(sentinel_data) == 'SSC_mgL_R'] <- 'SSC_mgL'
```
## Investigate SSC flux at the mouth of the river

```{r}
# Prepare Q data
Q_full_mc <- Q_full %>% filter(channel=='MAIN STEM')
Q_full_mc$year_chunk <- cut(as.numeric(format(as.Date(Q_full_mc$date), "%Y")), seq(from = 1985, to = 2026, by = 5), right=FALSE)
Q_full_mc$flow <- as.numeric(Q_full_mc$flow)
Q_full_mc$dist_downstream_km <- as.integer(as.numeric(Q_full_mc$dist_downstream_km))

# Prepare SSC data
SSC_full_mc <- SSC_full %>% filter(channel=='MAIN STEM')
SSC_full_mc$year_chunk <- cut(as.numeric(format(as.Date(SSC_full_mc$date), "%Y")), seq(from = 1985, to = 2026, by = 5), right=FALSE)
SSC_full_mc$sscmg_L <- as.numeric(SSC_full_mc$sscmg_L)
SSC_full_mc$dist_downstream_km <- as.integer(as.numeric(SSC_full_mc$dist_downstream_km))

# List unique ssc distances/q distances to create a rating curve
q_dists <- unique(Q_full_mc$dist_downstream_km)
ssc_dists <- unique(SSC_full_mc$dist_downstream_km)
similar_dists <- sort(intersect(q_dists, ssc_dists))
last_dist <- similar_dists[length(similar_dists)]

# Get station number at that distance
mouth_site <- SSC_full_mc[dist_downstream_km==last_dist,][1]$site_no
mouth_ssc <- SSC_full_mc %>% filter(site_no==mouth_site)
mouth_q <- Q_full_mc %>% filter(site_no==mouth_site)

# Need to convert  flow ft^3/2 to m^3/2
mouth_q$flow_m3s <- mouth_q$flow * 0.02832

# Create training data
ssc_q <- left_join(mouth_ssc, dplyr::select(mouth_q, c(date, flow_m3s)), by= "date")
ssc_q$SSC_flux_MTyr <- ssc_q$flow_m3s * ssc_q$sscmg_L * 3.10585 * 10**-5
ssc_q$log10_SSC_flux_MTyr <- log10(ssc_q$SSC_flux_MTyr)
ssc_q$log10_discharge_m3s  <- log10(ssc_q$flow_m3s)

# To get accurate means reduce to get average per day then average over the whole year
ssc_day_avg <- ssc_q %>% group_by(date) %>% 
  summarise_at(c("SSC_flux_MTyr","log10_SSC_flux_MTyr","log10_discharge_m3s","flow_m3s","sscmg_L"), mean) %>%
  left_join(subset(ssc_q, 
                   select=-c(SSC_flux_MTyr, log10_SSC_flux_MTyr, log10_discharge_m3s, flow_m3s, sscmg_L)), by="date") %>%
  distinct()

# Save training data
write.csv(ssc_q, paste0(wd_root,"/exports/", mouth_site,".csv"), row.names=FALSE)
write.csv(ssc_day_avg, paste0(wd_root,"/exports/", mouth_site,"_reduced_by_day.csv"), row.names=FALSE)

# Create flux plot not on rating curve
ssc_flux <- ggplot(ssc_day_avg, aes(x=format(as.Date(ssc_day_avg$date), "%Y"), y=SSC_flux_MTyr)) +
  # geom_point() +
  stat_summary(geom = 'line', fun = 'mean', group=1, aes(color="Average SSC Flux (MT/yr)")) +
  stat_summary(geom = 'errorbar', aes(color='1 standard Error'), fun.data= mean_se,
               fun.args = list(mult = 1), # mult = # of standard errors
               size = 0.5,
               width =  0.2,
               show.legend=TRUE) +
  # coord_cartesian(ylim=c(0,3)) +
  scale_color_manual(values = c("green", "red")) +
  season_facet +
  labs(title="SSC Flux at the Mouth of \nthe Chattahoochee (MT/yr)",
       caption=paste0('USGS Site ',mouth_site))

ssc_flux

ggsave(ssc_flux, 
       filename = paste0(wd_root, paste0('/exports/USGS-',mouth_site,'_flux.png')),
       width = 12, height = 8)

# Create rating curve
reg = lm(log10_SSC_flux_MTyr~log10_discharge_m3s, data=ssc_q)

# Plot regression

reg_plot = ggplot(ssc_q, aes(x = log10_discharge_m3s, y=log10_SSC_flux_MTyr)) +
  geom_point(color='green') +
  geom_smooth(method='lm', formula= y~x, color='black') +
  annotate(geom='text',label = lm_eqn(reg), parse = TRUE, x = -Inf, y = Inf, hjust = -0.2, vjust = 2) +
  theme_bw() +
  scale_fill_manual(values = c('gray'), name = "Metrics")  +
  labs(title=paste0('Rating Curve for Station USGS-',mouth_site),
       caption=paste0(startDate,' to ', endDate))

ggsave(reg_plot, 
       filename = paste0(wd_root, paste0('/exports/USGS-',mouth_site,'.png')),
       width = 10, height = 8)

# Apply regression to discharge data
mouth_q$log10_discharge_m3s <- log10(mouth_q$flow_m3s)

ssc_flux = predict.lm(reg, newdata=mouth_q)
mouth_q$log10_SSC_flux_MTyr <- ssc_flux
mouth_q$SSC_flux_MTyr <- 10**ssc_flux
mouth_q$SSC_mgL <- (10**ssc_flux) / (mouth_q$flow_m3s * 3.10585 * 10**-5)
mouth_q$log10_SSC_mgL <- log10((10**ssc_flux) / (mouth_q$flow_m3s * 3.10585 * 10**-5))

pred <- mouth_q
# Save predicted data
write.csv(pred, paste0(wd_root,"/exports/", mouth_site,"_ssc_pred.csv"), row.names=FALSE)

# Plot SSC flux over time
ssc_flux_plot <- ggplot(pred, aes(x=year_chunk, y=SSC_flux_MTyr)) +
  # geom_point() +
  stat_summary(geom = 'line', fun = 'mean', group=1, aes(color="Average SSC Flux (MT/yr)")) +
  stat_summary(geom = 'errorbar', aes(color='1 standard Error'), fun.data= mean_se,
               fun.args = list(mult = 1), # mult = # of standard errors
               size = 0.5,
               width =  0.2,
               show.legend=TRUE) +
  coord_cartesian(ylim=c(0,1.2)) +
  scale_color_manual(values = c("green", "red")) +
  season_facet +
  labs(title="Rating Curve Derived SSC\nFlux at the Mouth of \nthe Chattahoochee (MT/yr)",
       caption=paste0('USGS Site ',mouth_site))

ssc_flux_plot

ggsave(ssc_flux_plot, 
       filename = paste0(wd_root, paste0('/exports/USGS-',mouth_site,'_rating_flux.png')),
       width = 12, height = 8)

# Plot discharge over time
# Got annual flow from https://waterdata.usgs.gov/nwis/annual?referred_module=sw&amp;site_no=02359170&amp;por_02359170_26942=2396790,00060,26942,1977,2022&amp;start_dt=1984&amp;end_dt=2022&amp;year_type=W&amp;format=html_table&amp;date_format=YYYY-MM-DD&amp;rdb_compression=file&amp;submitted_form=parameter_selection_list
annual_q <- data.table(year=seq(1984,2021), flow_ft3s=c(
  31460,
  16070,
  19050,
  25600,
  16380,
  19270,
  29040,
  26200,
  21590,
  28180,
  30900,
  26820,
  26300,
  24250,
  36860,
  18020,
  10990,
  16890,
  10560,
  28960,
  18170,
  35230,
  16330,
  12420,
  14010,
  23590,
  32160,
  12100,
  9715,
  22340,
  25360,
  18680,
  28020,
  19050,
  30210,
  32090,
  29210
))
annual_q$flow_m3s <- annual_q$flow_ft3s * 0.02832

discharge_plot<- ggplot(annual_q, aes(x=year, y=flow_m3s)) +
  stat_summary(geom = 'line', fun = 'mean', group=1) +
  season_facet +
  coord_cartesian(ylim=c(0,1100)) +
  labs(title="Average Discharge at the Mouth of the Chattahoochee (m3/s)",
       x="year") +
  stat_summary(aes(label=round(..y..,0)), fun='mean', geom="text", size=2,
             vjust = -1) 

ggsave(discharge_plot, 
       filename = paste0(wd_root, paste0('/exports/USGS-',mouth_site,'_discharge.png')),
       width = 20, height = 8)

```

## Investigate SSC FLux along the whole river transect

```{r}
# Import transect topography
transect_topo <- read.csv("D:\\valencig\\Thesis\\Data\\chattahoochee_landsat\\chattahoochee_transect_topo.csv")
# Cant convert to km2 because cell size varies with latitude
transect_topo <- data.frame(
  distance_km = as.numeric(transect_topo$distance_km),
  elev_min = as.numeric(transect_topo$elev_min),
  flow_acc_max= as.numeric(transect_topo$flow_acc_max) 
  
)

ggplot(transect_topo, aes(x=distance_km, y=flow_acc_max))+
  geom_line()

# Need to convert  flow ft^3/2 to m^3/2
Q_full_mc$flow_m3s <- Q_full_mc$flow * 0.02832

# Convert distance downstream to nearest 2 km to match with landsat data
Q_full_mc$dist_downstream_2km <- round(Q_full_mc$dist_downstream_km/2)*2

# ggplot(Q_full_mc, aes(x=dist_downstream_km, y=flow_m3s)) +
#   stat_summary(geom='line', fun='mean', na.rm=T) +
#   facet_wrap(~year_chunk)

# Merge flow data with landsat data
# landsat_flux <- landsat_data_filter %>% dplyr::select(-c(lat, lon))
landsat_flux <- left_join(landsat_data_filter, 
                          subset(Q_full_mc, select=c(date, dist_downstream_2km, flow_m3s)), 
                          by=c('date', 'distance_km'='dist_downstream_2km'))
# Merge topographic data with landat data
landsat_flux <- left_join(landsat_flux, transect_topo, by="distance_km")
# Set identifier for distances with flux
landsat_flux$in_situ_flow <- F
landsat_flux[!is.na(flow_m3s),]$in_situ_flow <- T
landsat_flux$flow_m3s_log10 <- log10(landsat_flux$flow_m3s)
landsat_flux$flow_acc_log10 <- log10(landsat_flux$flow_acc_max)

# Merge flow data with sentinel data
sentinel_flux <- left_join(sentinel_data, 
                          subset(Q_full_mc, select=c(date, dist_downstream_2km, flow_m3s)), 
                          by=c('date', 'distance_km'='dist_downstream_2km'))
# Merge topographic data with landat data
sentinel_flux <- left_join(sentinel_flux, transect_topo, by="distance_km")
# Set identifier for distances with flux
sentinel_flux$in_situ_flow <- F
sentinel_flux[!is.na(flow_m3s),]$in_situ_flow <- T
sentinel_flux$flow_m3s_log10 <- log10(sentinel_flux$flow_m3s)
sentinel_flux$flow_acc_log10 <- log10(sentinel_flux$flow_acc_max)

# Create folder to hold regression plots
wd_landsat_Q <- paste0(wd_root, '/exports/landsat_Q_reg/')
if(!dir.exists(wd_landsat_Q)){
    dir.create(wd_landsat_Q)
}
wd_sentinel_Q <- paste0(wd_root, '/exports/sentinel_Q_reg/')
if(!dir.exists(wd_sentinel_Q)){
    dir.create(wd_sentinel_Q)
}

# Fill in Q_full_mc data
Q_full_mc$flow_m3s_log10 <- log10(Q_full_mc$flow_m3s)
names(Q_full_mc)[names(Q_full_mc) == "dist_downstream_2km"] <- "distance_km"
Q_full_mc <- left_join(Q_full_mc, transect_topo, by=("distance_km"))
Q_full_mc$flow_acc_log10 <- log10(as.numeric(Q_full_mc$flow_acc_max))
# Remove NaN values created from log10
Q_full_mc <- Q_full_mc %>% filter(is.finite(flow_m3s_log10))

# Baseline Attempt --> create discharge vs. drainage curve every day
naive_reg <- function(.x, dir){
  # in_situ_flow <- .x %>% filter(in_situ_flow == T)
  # Case 1: no in-situ flow measurements at that day
  # Issue: around 700 of 1650 dates of landsat sample do not have in-situ flow information
  # if (nrow(in_situ_flow) == 0) {
  #   .x$added_flow = F
  #   return(.x)
  # }
  in_situ_flow <- Q_full_mc %>% filter(date==.x$date2[1])
  # Case 2: There are flow measurements, perform linear regression
  reg <- lm(flow_m3s_log10~flow_acc_log10, data = in_situ_flow)
  eq <- substitute(italic('Flow_m3s') == 10^a %.%italic('Hydroshed_pixels')^b, 
         list(a = format(unname(coef(reg)[1]), digits = 3),
              b = format(unname(coef(reg)[2]), digits = 3)))
    
  reg_plot = ggplot(in_situ_flow, aes(x=flow_acc_log10, y=flow_m3s_log10)) +
    geom_point() +
    geom_smooth(method = "lm", se=FALSE, color="green", formula = y ~ x) +
    annotate(geom='text', label = as.character(as.expression(eq)), 
             parse = T, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    theme_minimal() +
    labs(title=paste0('Date', .x$date2[1]),
        x='Drainage Area (log10 Hydroshed pixels)',
        y='Discharge (log10 m3s)')
  
  # Save regression plot
  ggsave(reg_plot, 
       filename = paste(dir, .x$date2[1], '.pdf'),
       width = 8, height = 8)
  .x[.x$in_situ_flow == F, ]$flow_m3s_log10 <- predict(reg, .x[.x$in_situ_flow == F,])
  .x[.x$in_situ_flow == F, ]$flow_m3s <- 10**.x[.x$in_situ_flow == F,]$flow_m3s_log10
  return(.x)
}

# Add dummy date to access date inside of group modify
landsat_flux$date2 <- landsat_flux$date

landsat_flux <- landsat_flux %>%
  group_by(date) %>%
  group_modify(~naive_reg(., wd_landsat_Q)) %>% 
  subset(select=-c(date2))
# %>% 
#   filter(flow_m3s_predict>=0)

# Add dummy date to access date inside of group modify
sentinel_flux$date2 <- sentinel_flux$date

sentinel_flux <- sentinel_flux %>%
  group_by(date) %>%
  group_modify(~naive_reg(., wd_sentinel_Q)) %>% 
  subset(select=-c(date2))
# %>% 
#   filter(flow_m3s_predict>=0)

sentinel_flux <- sentinel_flux %>% filter(SSC_mgL >0)

ggplot(landsat_flux, aes(x=distance_km, y=flow_m3s)) +
  # geom_point() + 
  stat_summary(fun='mean', geom="line", aes(color="Landsat")) +
  stat_summary(data = sentinel_flux, fun='mean', geom="line", aes(color="Sentinel")) +
  stat_summary(data = Q_full_mc, fun='mean', geom="line", aes(color="in-situ")) +
  theme_minimal() +
  # coord_cartesian(ylim=c(0,750)) +
  labs(title="Comparison of in-situ and regressed discharge")

```

# Continue investigating flux at every point

```{r}
# Determine flux at every point
landsat_flux$SSC_flux_MTyr <- landsat_flux$flow_m3s * landsat_flux$SSC_mgL * 3.10585 * 10**-5
sentinel_flux$SSC_flux_MTyr <- sentinel_flux$flow_m3s * sentinel_flux$SSC_mgL * 3.10585 * 10**-5

# Reduce and get average for one sample per day
landsat_day_avg <- landsat_flux %>% group_by(distance_km, date) %>% 
  summarise_at(c("SSC_flux_MTyr","flow_m3s","SSC_mgL"), mean) %>%
  left_join(subset(landsat_flux, 
                   select=-c(SSC_flux_MTyr, flow_m3s, SSC_mgL)), by=c("date", "distance_km")) %>%
  distinct()

# Reduce and get average for one sample per day
sentinel_day_avg <- sentinel_flux %>% group_by(distance_km, date) %>% 
  summarise_at(c("SSC_flux_MTyr","flow_m3s","SSC_mgL"), mean) %>%
  left_join(subset(sentinel_flux, 
                   select=-c(SSC_flux_MTyr, flow_m3s, SSC_mgL)), by=c("date", "distance_km")) %>%
  distinct()
# Add year chunk to day avg
sentinel_day_avg$year_chunk <- cut(as.numeric(sentinel_day_avg$year), seq(from = 1985, to = 2026, by = 5), right=FALSE)

transect_flux_plot <- ggplot(landsat_day_avg, 
                             aes(x=distance_km, y=SSC_flux_MTyr, 
                                 color="Landsat Derived\nSSC Flux (MT/yr)")) +
  # geom_point() +
  stat_summary(geom = 'line', fun = 'mean', group=1, 
               aes(color="Average Landsat SSC Flux (MT/yr)"), na.rm=T) +
  stat_summary(data = sentinel_day_avg, geom = 'line', fun = 'mean', group=1, 
               aes(color="Average Sentinel SSC Flux (MT/yr)"), na.rm=T) +
  # stat_summary(geom = 'errorbar', aes(color='1 standard Error'), fun.data= mean_se,
  #              fun.args = list(mult = 1), # mult = # of standard errors
  #              size = 0.5,
  #              width =  0.2,
  #              show.legend=TRUE) +
  # coord_cartesian(ylim=c(0,2)) +
  scale_color_manual(values = c('#999999','#E69F00')) +
  theme_minimal() +
  facet_wrap(~year_chunk) +
  theme(legend.position = "bottom") +
  labs(title="SSC Flux along the Chattahoochee Transect",
       caption="Reduced to one average sample per day per transect section.")

transect_flux_plot

transect_flow_plot <- ggplot(landsat_day_avg, 
                             aes(x=distance_km, y=flow_m3s, color=" Average Landsat Q")) +
  # geom_point() +
  stat_summary(geom = 'line', fun = 'mean', group=1, na.rm=T) +
  stat_summary(data=sentinel_day_avg, aes(color="Average Sentinel Q",), geom = 'line', fun = 'mean', group=1, na.rm=T) +
  scale_color_manual(values = c('#999999','#E69F00')) +
  theme_minimal() +
  facet_wrap(~year_chunk) +
  theme(legend.position = "bottom") +
  labs(title="Discharge along the Chattahoochee Transect (m3/s)",
       caption="Reduced to one average sample per day per transect section.")

transect_flow_plot

# Take last 10 km of river
landsat_mouth <- landsat_day_avg %>% filter(distance_km >= (max(landsat_day_avg$distance_km)-10))
# %>%filter(distance_km<=max(landsat_day_avg$distance_km)-8 & distance_km>=max(landsat_day_avg$distance_km)-12)
sentinel_mouth <- sentinel_day_avg %>% filter(distance_km >= (max(sentinel_day_avg$distance_km)-10))

flux_mouth_plot <- ggplot(landsat_mouth, aes(x=year_chunk, y=SSC_flux_MTyr, color="Landsat Derived\nSSC Flux (MT/yr)")) +
  # geom_point() +
  stat_summary(geom = 'line', fun = 'mean', group=1, 
               aes(color="Landsat Derived\nSSC Flux (MT/yr)")) +
  stat_summary(data=sentinel_mouth, geom = 'line', fun = 'mean', group=1, 
               aes(color="Sentinel Derived\nSSC Flux (MT/yr)")) +
  # stat_summary(geom = 'errorbar', aes(color='1 standard Error'), fun.data= mean_se,
  #              fun.args = list(mult = 1), # mult = # of standard errors
  #              size = 0.5,
  #              width =  0.2,
  #              show.legend=TRUE) +
  # coord_cartesian(ylim=c(0,2)) +
  scale_color_manual(values = c('#999999','#E69F00')) +
  theme_minimal()+
  theme(legend.position = "bottom") +
  labs(title="SSC Flux at the mouth of the chattahoochee",
       caption="Reduced to one average sample per day.")

flux_mouth_plot

# Combine with no rating curve in-situ measurements
mouth_flux_compare <- ggplot(landsat_mouth, aes(x=year_chunk, y=SSC_flux_MTyr)) +
  # geom_point() +
  stat_summary(geom = 'line', fun = 'mean', group=1, 
               aes(color="Landsat Derived\nSSC Flux (MT/yr)")) +
  stat_summary(data = sentinel_mouth, geom = 'line', fun = 'mean', group=1, 
               aes(color="Sentinel Derived\nSSC Flux (MT/yr)")) +
  # stat_summary(geom = 'errorbar', 
  #              aes(color='1 standard deviation'), 
  #              fun.data= mean_sdl,
  #              fun.args = list(mult = 1), # mult = # of standard errors
  #              size = 0.5,
  #              width =  0.2,
  #              show.legend=TRUE) +
  # geom_point() +
  stat_summary(data=ssc_day_avg, geom = 'line', fun = 'mean', group=1, aes(x=year_chunk, y=SSC_flux_MTyr, color="In-Situ SSC Flux (MT/yr)")) +
  # stat_summary(data=ssc_day_avg, 
  #              geom = 'errorbar', 
  #              aes(x=year_chunk, y=SSC_flux_MTyr, color='1 standard deviation'), 
  #              fun.data= mean_sdl,
  #              fun.args = list(mult = 1), # mult = # of standard errors
  #              size = 0.5,
  #              width =  0.2,
  #              show.legend=TRUE) +
  theme_minimal()+
  theme(legend.position = "bottom") +
  coord_cartesian(ylim=c(0, 5)) +
  scale_color_manual(values = c('#9CDB79', '#999999','#E69F00')) +
  labs(title="SSC Flux at the mouth of the chattahoochee",
       caption="Reduced to one average sample per day.",
       x="year")

mouth_flux_compare 

# Save plots
ggsave(transect_flux_plot, 
       filename = paste0(wd_root, '/exports/transect_flux.png'),
       width=12, height=8)

ggsave(transect_flow_plot, 
       filename = paste0(wd_root, '/exports/transect_flow.png'),
       width = 12, height = 8)

ggsave(flux_mouth_plot, 
       filename = paste0(wd_root, '/exports/lmouth_flux.png'),
       width = 12, height = 8)

ggsave(mouth_flux_compare , 
       filename = paste0(wd_root, '/exports/mouth_flux_compare_year_chunk.png'),
       width = 20, height = 8)
```

# Save variables

```{r}
# save(Q_full_mc, file=paste0(save_dir, "Q_full_mc.RData"))
# save(landsat_data_filter, file=paste0(save_dir, "landsat_data_filter.RData"))
# save(landsat_flux, file=paste0(save_dir, "landsat_flux.RData"))
# save(landsat_day_avg, file=paste0(save_dir, "landsat_day_avg.RData"))
# save(landsat_mouth, file=paste0(save_dir, "landsat_mouth.RData"))
# save(sentinel_data, file=paste0(save_dir, "sentinel_data.RData"))
# save(sentinel_flux, file=paste0(save_dir, "sentinel_flux.RData"))
# save(sentinel_day_avg, file=paste0(save_dir, "sentinel_day_avg.RData"))
# save(sentinel_mouth, file=paste0(save_dir, "sentinel_mouth.RData"))
# save(ssc_day_avg, file=paste0(save_dir, "ssc_day_avg.RData"))
# save(ssc_q, file=paste0(save_dir, "ssc_q.RData"))

load(file=paste0(save_dir, "Q_full_mc.RData"))
load(file=paste0(save_dir, "landsat_data_filter.RData"))
load(file=paste0(save_dir, "landsat_flux.RData"))
load(file=paste0(save_dir, "landsat_day_avg.RData"))
load(file=paste0(save_dir, "landsat_mouth.RData"))
load(file=paste0(save_dir, "sentinel_data.RData"))
load(file=paste0(save_dir, "sentinel_flux.RData"))
load(file=paste0(save_dir, "sentinel_day_avg.RData"))
load(file=paste0(save_dir, "sentinel_mouth.RData"))
load(file=paste0(save_dir, "ssc_day_avg.RData"))
load(file=paste0(save_dir, "ssc_q.RData"))
load(file=paste0(save_dir, "mouth_q.RData"))

# Load SSC data
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "SSC_full.RData"))
SSC_full$year_chunk <- cut(as.numeric(format(as.Date(SSC_full$date), "%Y")), 
                           seq(from = 1985, to = 2026, by = 5), right=FALSE)
SSC_full$dist_downstream_km <- as.numeric(SSC_full$dist_downstream_km)
SSC_full$distance_km <- 
SSC_full$sscmg_L <- as.numeric(SSC_full$sscmg_L)
```

# Investigate some stats using continuous discharge

Use 'landsat_mouth' and 'ssc_day_avg' to avoid sampling bias

```{r}
t1 <- levels(landsat_flux$year_chunk)[1]
t2 <- levels(landsat_flux$year_chunk)[8]

ylims_landsat <- landsat_mouth %>%
  group_by(year_chunk) %>%
  summarise(Q1 = quantile(SSC_flux_MTyr, 1/4, na.rm=T), Q3 = quantile(SSC_flux_MTyr, 3/4, na.rm=T)) %>%
  ungroup() %>%
  #get lowest Q1 and highest Q3
  summarise(lowQ1 = min(Q1), highQ3 = max(Q3))

ylims_ssc <- ssc_day_avg %>%
  group_by(year_chunk) %>%
  summarise(Q1 = quantile(SSC_flux_MTyr, 1/4, na.rm=T), Q3 = quantile(SSC_flux_MTyr, 3/4, na.rm=T)) %>%
  ungroup() %>%
  #get lowest Q1 and highest Q3
  summarise(lowQ1 = min(Q1), highQ3 = max(Q3))

# Violin Plots 

ggplot(landsat_mouth, aes(x=year_chunk, y=SSC_flux_MTyr)) +
  # coord_cartesian(ylim=c(0,3)) +
  coord_cartesian(ylim = as.numeric(ylims_landsat)*1.05) +
  geom_violin(scale="count", aes(fill = year_chunk)) +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  scale_fill_brewer(palette="Blues") + theme_classic()

ggplot(ssc_day_avg, aes(x=year_chunk, y=SSC_flux_MTyr)) +
  # coord_cartesian(ylim=c(0,3)) +
  coord_cartesian(ylim = as.numeric(ylims_ssc)*1.05) +
  geom_violin(scale="count", aes(fill = year_chunk)) +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  scale_fill_brewer(palette="Blues") + theme_classic()

# Side by side box plot
mouth_boxplots <- ggplot(ssc_day_avg, aes(x=year_chunk, y=SSC_flux_MTyr, color="In-Situ")) +
  coord_cartesian(ylim=c(0,5)) +
  # coord_cartesian(ylim = c(0,as.numeric(ylims_landsat)[2]*1.10)) +
  # geom_violin(scale="count", aes(fill = year_chunk)) +
  geom_boxplot(width=0.2, outlier.shape=NA, position= position_nudge(x=-.13)) +
  geom_boxplot(data=landsat_mouth, aes(color="Landsat"), width=0.2, outlier.shape=NA, 
               position= position_nudge(x=+.13)) +
  geom_boxplot(data=sentinel_mouth, aes(color="Sentinel"), width=0.2, outlier.shape=NA, 
               position= position_nudge(x=+.39)) +
  # Add labels for amount of observations
  # geom_text(data = landsat_flux, aes(label=count), 
  #           position=position_dodge(width=1.0)) +
  scale_color_manual(values = c('#9CDB79', '#999999','#E69F00')) +
  labs(title="Comparison of Flux at the Mouth of the Chattahoochee",
       caption="USGS samples are taken at site USGS-02359170, APALACHICOLA RIVER NR SUMATRA,FLA.\nFor Landsat/Sentinel, the mouth of the river is defined as the last 10km of river.")

mouth_boxplots

# ggsave(mouth_boxplots, 
#        filename = paste0(wd_root, '/exports/mouth_flux_boxplots.png'),
#        width = 20, height = 8)

```

# Investigate Creating a Rating Curve (at mouth) --> REMAKE PLOTS ABOVE
--> REDO WITH 

```{r}
# log10 flux vs log10 discharge

### Using landsat data <- landsat_mouth ###
landsat_mouth$SSC_flux_MTyr_log10 <- log10(landsat_mouth$SSC_flux_MTyr)
fluxFlow_landsat <- lm(SSC_flux_MTyr_log10~flow_m3s_log10, data=landsat_mouth)

fluxFlow_landsat_eq <- substitute(italic('SSC Flux (Mt/yr)') == 
                                    10^a %.%italic('Discharge (m^3/s)')^b~italic('R'^2)~':'~r2, 
       list(a = format(unname(coef(fluxFlow_landsat)[1]), digits = 3),
            b = format(unname(coef(fluxFlow_landsat)[2]), digits = 3),
            r2 = format(summary(fluxFlow_landsat)$r.squared, digits = 4)))

### Using Sentinel data <- sentinel_mouth ###
sentinel_mouth$SSC_flux_MTyr_log10 <- log10(sentinel_mouth$SSC_flux_MTyr)
fluxFlow_sentinel<- lm(SSC_flux_MTyr_log10~flow_m3s_log10, data=sentinel_mouth)

fluxFlow_sentinel_eq <- substitute(italic('SSC Flux (Mt/yr)') == 
                                    10^a %.%italic('Discharge (m^3/s)')^b~italic('R'^2)~':'~r2, 
       list(a = format(unname(coef(fluxFlow_sentinel)[1]), digits = 3),
            b = format(unname(coef(fluxFlow_sentinel)[2]), digits = 3),
            r2 = format(summary(fluxFlow_sentinel)$r.squared, digits = 4)))

### Using USGS data <- ssc_day_avg ###
# Add duplicate columns for ease of use
ssc_day_avg$SSC_flux_MTyr_log10 <- ssc_day_avg$log10_SSC_flux_MTyr
ssc_day_avg$flow_m3s_log10 <- ssc_day_avg$log10_discharge_m3s
fluxFlow_USGS<- lm(SSC_flux_MTyr_log10~flow_m3s_log10, data=ssc_day_avg)

fluxFlow_USGS_eq <- substitute(italic('SSC Flux (Mt/yr)') == 
                                    10^a %.%italic('Discharge (m^3/s)')^b~italic('R'^2)~':'~r2, 
       list(a = format(unname(coef(fluxFlow_USGS)[1]), digits = 3),
            b = format(unname(coef(fluxFlow_USGS)[2]), digits = 3),
            r2 = format(summary(fluxFlow_USGS)$r.squared, digits = 4)))

### Display Rating Curve Plots ###
ggplot(landsat_mouth, aes(x=flow_m3s_log10, y=SSC_flux_MTyr_log10)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, color="green", formula = y ~ x) +
  annotate(geom='text', label = as.character(as.expression(fluxFlow_landsat_eq)), 
           parse = T, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
  theme_minimal()+
  labs(title=paste0('Landsat Derived Flux-Discharge Rating Curve'),
      x='Sediment Flux (log10 Mt/yr)',
      y='Discharge (log10 m3s)')

ggplot(sentinel_mouth, aes(x=flow_m3s_log10, y=SSC_flux_MTyr_log10)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, color="green", formula = y ~ x) +
  annotate(geom='text', label = as.character(as.expression(fluxFlow_sentinel_eq)), 
           parse = T, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
  theme_minimal()+
  labs(title=paste0('Sentinel Derived Flux-Discharge Rating Curve'),
      x='Sediment Flux (log10 Mt/yr)',
      y='Discharge (log10 m3s)')

ggplot(ssc_day_avg, aes(x=flow_m3s_log10, y=SSC_flux_MTyr_log10)) +
  geom_point() +
  geom_smooth(method = "lm", se=FALSE, color="green", formula = y ~ x) +
  annotate(geom='text', label = as.character(as.expression(fluxFlow_USGS_eq)), 
           parse = T, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
  theme_minimal()+
  labs(title=paste0('USGS Derived Flux-Discharge Rating Curve'),
      x='Sediment Flux (log10 Mt/yr)',
      y='Discharge (log10 m3s)')

### Apply rating curves to discharge data at mouth of river -> mouth_q ###
mouth_q$rating_flux_landsat <- 10**predict(fluxFlow_landsat, mouth_q)
mouth_q$rating_flux_sentinel <- 10**predict(fluxFlow_sentinel, mouth_q)
#sentinel clip <- only to valid time chunk where sentinel data exists
# mouth_q$rating_flux_sentinelClip <- -100000 # need column to be numeric
# mouth_q[mouth_q$year_chunk == unique(sentinel_mouth$year_chunk),]$rating_flux_sentinelClip <- 10**predict(fluxFlow_sentinel, mouth_q[mouth_q$year_chunk == unique(sentinel_mouth$year_chunk),])
# mouth_q[mouth_q$rating_flux_sentinelClip == -100000, ]$rating_flux_sentinelClip <- NA # filter our Null data
mouth_q$rating_flux_USGS<- 10**predict(fluxFlow_USGS, mouth_q)

# Boxplot of using reverse rating curve

ylims_flux <- mouth_q %>%
  group_by(year_chunk) %>%
  summarise(Q1 = quantile(rating_flux_landsat, 1/4, na.rm=T), Q3 = quantile(rating_flux_landsat, 3/4, na.rm=T)) %>%
  ungroup() %>%
  #get lowest Q1 and highest Q3
  summarise(lowQ1 = min(Q1), highQ3 = max(Q3))

ggplot(mouth_q, aes(x=year_chunk, y=rating_flux_USGS, color="USGS")) +
  # coord_cartesian(ylim=c(0,3)) +
  coord_cartesian(ylim = c(0, as.numeric(ylims_flux[2])*1.10)) +
  # geom_violin(scale="count", aes(fill = year_chunk)) +
  geom_boxplot(width=0.2, outlier.shape=NA, position= position_nudge(x=-.13)) +
  geom_boxplot(aes(y=rating_flux_landsat, color="Landsat"), width=0.2, outlier.shape=NA, 
               position= position_nudge(x=+.13)) +
  geom_boxplot(aes(y=rating_flux_sentinel, color="Sentinel"), width=0.2, outlier.shape=NA, 
               position= position_nudge(x=+.39)) +
  # Add labels for amount of observations
  # geom_text(data = landsat_flux, aes(label=count), 
  #           position=position_dodge(width=1.0)) +
  scale_color_manual(values = c('#9CDB79', '#999999','#E69F00')) +
  theme_minimal()+
  labs(title="FILL",
       caption="USGS samples are taken at site USGS-02359170, APALACHICOLA RIVER NR SUMATRA,FLA.\nFor Landsat/Sentinel, the mouth of the river is defined as the last 10km of river.",
       y="Sediment FLux (Mt/yr)",
       x="Year Grouping")

### Compare to No rating Curve Data ###
ggplot(landsat_mouth, aes(x=year_chunk, y=SSC_flux_MTyr, color="Landsat Derived\nSSC Flux (MT/yr)")) +
  # geom_point() +
  coord_cartesian(ylim=c(0,3)) +
  stat_summary(geom = 'line', fun = 'mean', group=1, 
               aes(color="Landsat")) +
  # Add number of observations from landsat
  geom_label(aes(label=..count.., size=..count..), y=2.75,stat='count', color='red') +
  # stat_summary(fun= 'count', geom='text') +

  stat_summary(data=sentinel_mouth, geom = 'line', fun = 'mean', group=1, 
               aes(color="Sentinel-2")) +
    
  stat_summary(data = mouth_q, geom = 'line', fun = 'mean', group=1, 
               aes(y=rating_flux_landsat, color="Landsat\n(Rating Curve Derived)")) +
    
  stat_summary(data=mouth_q[mouth_q$year_chunk == unique(sentinel_mouth$year_chunk),], geom = 'line', fun = 'mean', group=1, 
               aes(y=rating_flux_sentinel, color="Sentinel-2\n(Rating Curve Derived)")) +
  
  stat_summary(data=mouth_q, geom = 'line', fun = 'mean', group=1, 
               aes(y=rating_flux_USGS, color="USGS\n(Rating Curve Derived)")) +
  
  scale_color_manual(values = c('#999999','#E69F00','blue','green','purple')) +
  geom_label(data=mouth_q[mouth_q$year_chunk == unique(sentinel_mouth$year_chunk),],
            aes(label=..count.., size=..count..), y=2.75-.3, stat='count', color='black') +
  theme_minimal()+
  guides(size = 'none') +
  theme(legend.position ="bottom") +
  labs(title="SSC Flux at the mouth of the chattahoochee",
       caption="Reduced to one average sample per day.",
       y="Sediment Flux (Mt/yr)",
       x="Year Grouping")

```

# Investigate using rating curve across the whole transect
CAN REAPPLY Q_date_regs ABOVE TO MAKE LIFE EASIER

```{r}
# Stores regression for every date with USGS data (checked consistency with landsat_Q_plots)
Q_date_regs <- Q_full_mc %>% group_by(date) %>% do(model=lm(flow_m3s_log10~flow_acc_log10, .))
Q_date_regs$year_chunk <- cut(as.numeric(format(as.Date(Q_date_regs$date), "%Y")), 
                           seq(from = 1985, to = 2026, by = 5), right=FALSE)
drainage_dist <- as.data.frame(landsat_day_avg) %>% dplyr::select(c("flow_acc_log10", "distance_km")) %>% distinct()

# Create folder to hold regression plots
wd_landsat_T <- paste0(wd_root, '/exports/landsat_transect_reg/')
if(!dir.exists(wd_landsat_T)){
    dir.create(wd_landsat_T)
}

wd_sentinel_T <- paste0(wd_root, '/exports/sentinel_transect_reg/')
if(!dir.exists(wd_sentinel_T)){
    dir.create(wd_sentinel_T)
}

# wd_USGS_T <- paste0(wd_root, '/exports/usgs_transect_reg/')
# if(!dir.exists(wd_USGS_T)){
#     dir.create(wd_USGS_T)
# }

# Baseline Attempt --> create discharge vs. drainage curve every day
naive_ssc_regL <- function(.x, dir){
  # print(paste0("Year Chunk: ",.x$year_chunk2[1]))
  # Get regression between log10 Q and log10 ssc_flux
  QF_reg <- lm(SSC_flux_MTyr_log10~flow_m3s_log10, data = .x)
  
  eq <- substitute(italic('Sediment Flux') == 10^a %.%italic('Discharge (m3/s)')^b, 
         list(a = format(unname(coef(QF_reg)[1]), digits = 3),
              b = format(unname(coef(QF_reg)[2]), digits = 3)))
    
  reg_plot = ggplot(.x, aes(x=flow_m3s_log10, y=SSC_flux_MTyr_log10)) +
    geom_density2d_filled(aes(fill = ..level..),bins = 10, show.legend = F) +
    geom_smooth(method = "lm", se=FALSE, color="green", formula = y ~ x) +
    annotate(geom='text', label = as.character(as.expression(eq)), 
             parse = T, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    labs(title=paste0('Year Chunk ', .x$year_chunk2[1]),
        x='Discharge (log10 m3/s)',
        y='Sediment Flux (log10 Mt/yr)',
        caption=paste0("R^2: ", format(summary(QF_reg)$r.squared, digits=4))) +
    # Themes
    theme_minimal() +
    # theme(text=element_text(family="LM Roman 10")) +
    # Remove background grid
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # scale_fill_paletteer_d("cartography::taupe.pal")
    scale_fill_manual(values=paletteer_dynamic("cartography::taupe.pal",10))
  
  ggsave(reg_plot,
       filename = paste0(dir, .x$year_chunk2[1], '.pdf'),
       width = 12, height = 8)
  
  # Get all Q vs Drainage Regs for the time period
  QD_regs <- Q_date_regs %>% filter(year_chunk == .x$year_chunk2[1])
  # Get dates of USGS rating curves
  dates <- QD_regs$date
  # Get unique dates (not in landsat data)
  unique_dates <- intersect(dates, unique(.x$date))
  
  # Create new dataframe to hold data
  new_df <- data.frame(date=rep(dates, each=nrow(drainage_dist)),
                       distance_km = rep(drainage_dist$distance_km, length(dates)),
                       flow_acc_log10 = rep(drainage_dist$flow_acc_log10, length(dates)))
  # Takes a while to fill in the dataset
  # print("making new df")
  # pb = txtProgressBar(min = 0, max = length(dates), initial = 0)
  # for (i in 1:length(dates)) {
  #   date <- dates[i]
  #   setTxtProgressBar(pb, i)
  #   # Add new rows to satellite dataframe
  #   new_df[new_df$date==date, ]$distance_km <- drainage_dist$distance_km
  #   new_df[new_df$date==date, ]$flow_acc_log10<- drainage_dist$flow_acc_log10
  # }
  # close(pb)
  # For all dates in satellite dataset get flow vs drainage
  new_df$flow_m3s_log10 <- 0
  new_df$SSC_flux_MTyr_log10 <- 0
  # Also takes a while
  print(paste0("Applying Regressions: ",.x$year_chunk2[1]))
  pb = txtProgressBar(min = 0, max = length(dates), initial = 0)
  offset <- nrow(drainage_dist)
  for (i in 1:length(dates)){
    setTxtProgressBar(pb, i)
    start <- (i-1)*offset+1
    end <- i*offset
    row <- new_df[i,]
    date <- dates[i]
    reg_eq <- QD_regs[QD_regs$date == row$date, ]$model[[1]]
    
    new_df[start:end,]$flow_m3s_log10 <- predict(reg_eq, new_df[start:end,])
    new_df[start:end,]$SSC_flux_MTyr_log10 <- predict(QF_reg, new_df[start:end,])
  }
  close(pb)
  # NEED TO CHANGE THIS FOR SENTINEL
  landsat_regressed_flux <<- rbind(landsat_regressed_flux, new_df)
  return(.x)
}

naive_ssc_regS <- function(.x, dir){
  # print(paste0("Year Chunk: ",.x$year_chunk2[1]))
  # Get regression between log10 Q and log10 ssc_flux
  QF_reg <- lm(SSC_flux_MTyr_log10~flow_m3s_log10, data = .x)
  
  eq <- substitute(italic('Sediment Flux') == 10^a %.%italic('Discharge (m3/s)')^b, 
         list(a = format(unname(coef(QF_reg)[1]), digits = 3),
              b = format(unname(coef(QF_reg)[2]), digits = 3)))
    
  reg_plot = ggplot(.x, aes(x=flow_m3s_log10, y=SSC_flux_MTyr_log10)) +
    geom_density2d_filled(aes(fill = ..level..),bins = 10, show.legend = F) +
    geom_smooth(method = "lm", se=FALSE, color="green", formula = y ~ x) +
    annotate(geom='text', label = as.character(as.expression(eq)), 
             parse = T, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    labs(title=paste0('Year Chunk ', .x$year_chunk2[1]),
        x='Discharge (log10 m3/s)',
        y='Sediment Flux (log10 Mt/yr)',
        caption=paste0("R^2: ", format(summary(QF_reg)$r.squared, digits=4))) +
    # Themes
    theme_minimal() +
    # theme(text=element_text(family="LM Roman 10")) +
    # Remove background grid
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # scale_fill_paletteer_d("cartography::taupe.pal")
    scale_fill_manual(values=paletteer_dynamic("cartography::taupe.pal",10))
  
  ggsave(reg_plot,
       filename = paste0(dir, .x$year_chunk2[1], '.pdf'),
       width = 12, height = 8)
  
  # Get all Q vs Drainage Regs for the time period
  QD_regs <- Q_date_regs %>% filter(year_chunk == .x$year_chunk2[1])
  # Get dates of USGS rating curves
  dates <- QD_regs$date
  # Get unique dates (not in landsat data)
  unique_dates <- intersect(dates, unique(.x$date))
  
  # Create new dataframe to hold data
  new_df <- data.frame(date=rep(dates, each=nrow(drainage_dist)),
                       distance_km = rep(drainage_dist$distance_km, length(dates)),
                       flow_acc_log10 = rep(drainage_dist$flow_acc_log10, length(dates)))
  # Takes a while to fill in the dataset
  # print("making new df")
  # pb = txtProgressBar(min = 0, max = length(dates), initial = 0)
  # for (i in 1:length(dates)) {
  #   date <- dates[i]
  #   setTxtProgressBar(pb, i)
  #   # Add new rows to satellite dataframe
  #   new_df[new_df$date==date, ]$distance_km <- drainage_dist$distance_km
  #   new_df[new_df$date==date, ]$flow_acc_log10<- drainage_dist$flow_acc_log10
  # }
  # close(pb)
  # For all dates in satellite dataset get flow vs drainage
  new_df$flow_m3s_log10 <- 0
  new_df$SSC_flux_MTyr_log10 <- 0
  # Also takes a while
  print(paste0("Applying Regressions: ",.x$year_chunk2[1]))
  pb = txtProgressBar(min = 0, max = length(dates), initial = 0)
  offset <- nrow(drainage_dist)
  for (i in 1:length(dates)){
    setTxtProgressBar(pb, i)
    start <- (i-1)*offset+1
    end <- i*offset
    row <- new_df[i,]
    date <- dates[i]
    reg_eq <- QD_regs[QD_regs$date == row$date, ]$model[[1]]
    
    new_df[start:end,]$flow_m3s_log10 <- predict(reg_eq, new_df[start:end,])
    new_df[start:end,]$SSC_flux_MTyr_log10 <- predict(QF_reg, new_df[start:end,])
  }
  close(pb)
  # NEED TO CHANGE THIS FOR SENTINEL
  sentinel_regressed_flux <<- rbind(sentinel_regressed_flux, new_df)
  return(.x)
}

# Add dummy distance to access date inside of group modify
landsat_day_avg_fluxrating <- landsat_day_avg %>% copy()
# landsat_day_avg_fluxrating$dist2 <- landsat_day_avg_fluxrating$distance_km
landsat_day_avg_fluxrating$year_chunk2 <- landsat_day_avg_fluxrating$year_chunk
landsat_day_avg_fluxrating$SSC_flux_MTyr_log10 <- log10(landsat_day_avg_fluxrating$SSC_flux_MTyr)
# landsat_day_avg_fluxrating <- landsat_day_avg_fluxrating %>% dplyr::select(c("date", "distance_km","flow_acc_log10","flow_m3s_log10","SSC_flux_MTyr_log10"))

sentinel_day_avg_fluxrating <- sentinel_day_avg %>% copy()
# landsat_day_avg_fluxrating$dist2 <- landsat_day_avg_fluxrating$distance_km
sentinel_day_avg_fluxrating$year_chunk2 <- sentinel_day_avg_fluxrating$year_chunk
sentinel_day_avg_fluxrating$SSC_flux_MTyr_log10 <- log10(sentinel_day_avg_fluxrating$SSC_flux_MTyr)
# sentinel_day_avg_fluxrating <- sentinel_day_avg_fluxrating %>% dplyr::select(c("date", "distance_km","flow_acc_log10","flow_m3s_log10","SSC_flux_MTyr_log10","year_chunk2"))

# landsat_day_avg_fluxrating$dist2 <- landsat_day_avg_fluxrating$distance_km
# SSC_full_fluxrating <- SSC_full %>% copy()
# SSC_full_fluxrating <- left_join(SSC_full_fluxrating, Q_full_mc)
# SSC_full_fluxrating$year_chunk2 <- SSC_full_fluxrating$year_chunk
# SSC_full_fluxrating$SSC_flux_MTyr_log10 <- log10(as.numeric(SSC_full_fluxrating$SSC_flux_MTyr))

# Need to do it per year per 2km segement
landsat_regressed_flux <- data.frame(date="", distance_km=0, flow_acc_log10=0, flow_m3s_log10=0, SSC_flux_MTyr_log10=0)
sentinel_regressed_flux <- data.frame(date="", distance_km=0, flow_acc_log10=0, flow_m3s_log10=0, SSC_flux_MTyr_log10=0)

landsat_day_avg_fluxrating <- landsat_day_avg_fluxrating %>%
  group_by(year_chunk) %>%
  group_modify(~naive_ssc_regL(., wd_landsat_T)) 

sentinel_day_avg_fluxrating <- sentinel_day_avg_fluxrating %>%
  group_by(year_chunk) %>%
  group_modify(~naive_ssc_regS(., wd_sentinel_T)) 
# 
# SSC_full_fluxrating <- SSC_full_fluxrating %>%
#   group_by(year_chunk) %>%
#   group_modify(~naive_ssc_reg(., wd_USGS_T)) 

# %>% 
#   dplyr::select(-c(dist2))
landsat_regressed_flux$SSC_flux_MTyr <- 10**landsat_regressed_flux$SSC_flux_MTyr_log10
sentinel_regressed_flux$SSC_flux_MTyr <- 10**sentinel_regressed_flux$SSC_flux_MTyr_log10

# Save regressed flux (rating curve for flow and flux)
save(landsat_regressed_flux, file=paste0(save_dir, "landsat_regressed_flux.RData"))
save(sentinel_regressed_flux, file=paste0(save_dir, "sentinel_regressed_flux.RData"))

# Save flux only using actual ssc samples
# THIS IS THE SAME AS LANDSAT_DAY_AVG
save(landsat_day_avg_fluxrating, file=paste0(save_dir, "landsat_day_avg_fluxrating.RData"))
save(sentinel_day_avg_fluxrating, file=paste0(save_dir, "sentinel_day_avg_fluxrating.RData"))

# save(usgs_regressed_flux, file=paste0(save_dir, "landsat_regressed_flux.RData"))
```
# Initalize data for time series comparison --> TESTING USING RATING CURVE DERIVED DATA

```{r}
load(file=paste0(save_dir, 'riverdata.Rdata'))
basin_area <- st_area(riverdata$basin) %>% set_units(km^2)

### PRE-ANTHROPOGENIC ###

# From octopusdata.org 
# Octupus ID: S156WTS010 Sample ID: SAP66
Be10E <- 9.15 %>% set_units(mm/kyr) %>% set_units(m/yr)
Be10dE <- 2.15 %>% set_units(mm/kyr) %>% set_units(m/yr)

# Get total flux over one year (using density of quartz, 2.648 g/cm^3) --> switched to 1440 from Trimble 1977
Be10flux <- (Be10E * set_units(basin_area, m^2) * set_units(1440, kg/m^3)) %>% set_units(Mt/yr)
Be10dflux <- (Be10dE * set_units(basin_area, m^2) * set_units(1440, kg/m^3)) %>% set_units(Mt/yr)

### PRE-DAM ###

# From Trimble (1977), sourced originally from Dole and Stabler (1909)
# They used assumed bulk density of 1440 kg/m^3 
# Included 10% additional to account for Bed Load (B.L/)
Tflux <- ((set_units(0.057, mm/yr) %>% set_units(m/yr)) 
          * set_units(basin_area, m^2) 
          * set_units(1440, kg/m^3)) %>% set_units(Mt/yr)

### POST-DAM ###

#USGS in-situ data (at station closest to mouth - USGS-02359170)
# Increase by 10% to account for bed load
USGSflux <- mean(ssc_day_avg$SSC_flux_MTyr * 1.1) %>% set_units(Mt/yr)
USGSdflux <- sd(ssc_day_avg$SSC_flux_MTyr * 1.1) %>% set_units(Mt/yr)
USGSQ <- mean(ssc_day_avg$flow_m3s) %>% set_units(m^3/s)
USGSdQ <- sd(ssc_day_avg$flow_m3s) %>% set_units(m^3/s)
USGSDisYield <- (USGSflux / basin_area) %>% set_units(t/(km^2*yr))
USGSSSC <- mean(ssc_day_avg$sscmg_L) %>% set_units(mg/L)
USGSdSSC <- sd(ssc_day_avg$sscmg_L) %>% set_units(mg/L)
USGSsamp <- nrow(ssc_day_avg)

# Landsat Data (use last 10 km of the river)
# Increase by 10% to account for bed load
# Landsatflux <- mean(landsat_mouth$SSC_flux_MTyr * 1.1, na.rm=T) %>% set_units(Mt/yr)
# Landsatdflux <- sd(landsat_mouth$SSC_flux_MTyr * 1.1, na.rm=T) %>% set_units(Mt/yr)
# LandsatQ <- mean(landsat_mouth$flow_m3s, na.rm=T) %>% set_units(m^3/s)
# LandsatdQ <- sd(landsat_mouth$flow_m3s, na.rm=T) %>% set_units(m^3/s)
LandsatSSC <- mean(landsat_mouth$SSC_mgL, na.rm=T) %>% set_units(mg/L)
LandsatdSSC <- mean(landsat_mouth$SSC_mgL, na.rm=T) %>% set_units(mg/L)
Landsatsamp <- nrow(landsat_mouth)
# To use rating curve derived data uncomment below
Landsatflux <- mean(mouth_q$rating_flux_landsat * 1.1, na.rm=T) %>% set_units(Mt/yr)
# Landsatflux5 <- mean(mouth_q[mouth_q$year_chunk %in% unique(sentinel_mouth$year_chunk),]$rating_flux_landsat * 1.1, na.rm=T) %>% set_units(Mt/yr)
Landsatdflux <- sd(mouth_q$rating_flux_landsat * 1.1, na.rm=T) %>% set_units(Mt/yr)
LandsatQ <- mean(mouth_q$flow_m3s, na.rm=T) %>% set_units(m^3/s)
LandsatdQ <- sd(mouth_q$flow_m3s, na.rm=T) %>% set_units(m^3/s)


# Sentinel Data (use last 10 km of the river)
# Increase by 10% to account for bed load
# Sentinelflux <- mean(sentinel_mouth$SSC_flux_MTyr * 1.1, na.rm=T) %>% set_units(Mt/yr)
# Sentineldflux <- sd(sentinel_mouth$SSC_flux_MTyr * 1.1, na.rm=T) %>% set_units(Mt/yr)
# SentinelQ <- mean(sentinel_mouth$flow_m3s, na.rm=T) %>% set_units(m^3/s)
# SentineldQ <- sd(sentinel_mouth$flow_m3s, na.rm=T) %>% set_units(m^3/s)
SentinelSSC <- mean(sentinel_mouth$SSC_mgL, na.rm=T) %>% set_units(mg/L)
SentineldSSC <- mean(sentinel_mouth$SSC_mgL, na.rm=T) %>% set_units(mg/L)
Sentinelsamp <- nrow(sentinel_mouth)

# To use rating curve derived data uncomment below
Sentinelflux <- mean(mouth_q[mouth_q$year_chunk %in% unique(sentinel_mouth$year_chunk),]$rating_flux_sentinel * 1.1, na.rm=T) %>% set_units(Mt/yr)
Sentineldflux <- sd(mouth_q[mouth_q$year_chunk %in% unique(sentinel_mouth$year_chunk),]$rating_flux_sentinel * 1.1, na.rm=T) %>% set_units(Mt/yr)
SentinelQ <- mean(mouth_q[mouth_q$year_chunk %in% unique(sentinel_mouth$year_chunk),]$flow_m3s, na.rm=T) %>% set_units(m^3/s)
SentineldQ <- sd(mouth_q[mouth_q$year_chunk %in% unique(sentinel_mouth$year_chunk),]$flow_m3s, na.rm=T) %>% set_units(m^3/s)



# Data from Milliman et. al
# Increase by 10% to account for bed load
MillimanQ <- (set_units(22, km^3/yr)) %>% set_units(m^3/s) # Average discharge
MillimanTDS <- set_units(1.1 * 1.1, Mt/yr) # TDS load
MillimanDisYield <- set_units(0.021*1000, t/(km^2*yr)) # They use basin area of 52,000 km^2 but use 52 as divisor, need to multiply by 1,000
MillimanSSC <- set_units(50, mg/L) # Average SSC concentration

# Flux Dataframe
Timeflux <- data.frame(source=c("Reusser (2015)", "Trimble (1977)", "USGS*", "Milliman (2013)*", "Dethier (2020)*","Sentinel-2**"),
                       flux=c(Be10flux, Tflux, USGSflux, MillimanTDS, Landsatflux, Sentinelflux),
                       period=c("Pre-Antropogenic", "Pre-Dam", "Post-Dam", "Post-Dam", "Post-Dam","Post-Dam"))
Timeflux$period <- factor(Timeflux$period, levels = c("Pre-Antropogenic", "Pre-Dam", "Post-Dam"))
Timeflux$source <- factor(Timeflux$source, levels = Timeflux$source)

# Only past two years
Timeflux <- data.frame(source=c("Reusser (2015)", "Trimble (1977)", "USGS*", "Milliman (2013)*", "Dethier (2020)*","Sentinel-2**"),
                       flux=c(Be10flux, Tflux, USGSflux, MillimanTDS, Landsatflux, Sentinelflux),
                       period=c("Pre-Antropogenic", "Pre-Dam", "Post-Dam", "Post-Dam", "Post-Dam","Post-Dam"))
Timeflux$period <- factor(Timeflux$period, levels = c("Pre-Antropogenic", "Pre-Dam", "Post-Dam"))
Timeflux$source <- factor(Timeflux$source, levels = Timeflux$source)

# Sampling dataframe
Samples <- data.frame(source = c("USGS", "Landsat","Sentinel-2"),
                      samples = c(USGSsamp, Landsatsamp, Sentinelsamp))

# Q dataframe
Qs <- data.frame(source = c("USGS", "Landsat","Sentinel-2", "Milliman"),
                      q = c(USGSQ, LandsatQ, SentinelQ,MillimanQ))
```

# Create plots over time periods --> ADD LAST TWO YEAR CHUNKS FOR ACCURATE SENTINEL COMPARISON
```{r}
# Compare flux over three time periods
time_compare <- ggplot(Timeflux, aes(x=source, y=flux, fill=period)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label=sprintf("%0.2f", round(flux, digits = 2))), vjust=-0.2) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c('#999999','#E69F00','#9CDB79')) +
  labs(title="Sediment Flux in the Chattahoochee Basin",
       x = "",
       caption = "Using bulk density of 1440 kg/m^3 from Trimble (1977) and assuming bed load accounts for an extra 10% across SSC derived flux.\n*Defined as 1985-present.\n**Data only available from 2018-2022.")

time_compare 

# Compare number of samples
samp_compare <- ggplot(Samples, aes(x=source, y=samples, fill=source)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label=samples), vjust=-0.2) +
  theme_minimal() +
  scale_fill_manual(values=c('#999999','#E69F00','#9CDB79')) +
  labs(title="Number of Samples (1985-2023)",
       x = "",
       caption = "Samples are taken at site USGS-02359170, APALACHICOLA RIVER NR SUMATRA,FLA.\nFor Landsat/Sentinel, the mouth of the river is defined as the last 10km of river.")

samp_compare

Q_compare <- ggplot(Qs, aes(x=source, y=q, fill=source)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label=sprintf("%0.0f", round(q, digits = 0))), vjust=-0.2) +
  theme_minimal() +
  scale_fill_manual(values=c('#999999','#9CDB79','#E69F00',"blue")) +
  labs(title="Average Discharge at the Mouth of the Chattahoochee",
       x = "",
       caption = "Samples are taken at site USGS-02359170, APALACHICOLA RIVER NR SUMATRA,FLA.\nFor Landsat/Sentinel, the mouth of the river is defined as the last 10km of river.")

Q_compare

# Save plots
ggsave(time_compare, 
       filename = paste0(wd_root, '/exports/time_compare.png'),
       width = 20, height = 8)
ggsave(samp_compare, 
       filename = paste0(wd_root, '/exports/samp_compare.png'),
       width = 20, height = 8)
ggsave(Q_compare, 
       filename = paste0(wd_root, '/exports/Q_compare.png'),
       width = 20, height = 8)
```

# Investigate large storms (now using rating curve derived data -> mouth_q)

```{r}
### LANDSAT ###

landsatQ_top <- mouth_q %>%
  # group_by(year_chunk) %>%
  summarise(top = quantile(flow_m3s, 0.95, na.rm=T)) %>%
  # ungroup() %>%
  # highest Q3
  summarise(high = max(top))

landsat_mouth_qB <- mouth_q %>% filter(flow_m3s < landsatQ_top$high)
landsat_mouth_qT <- mouth_q %>% filter(flow_m3s >= landsatQ_top$high)

# Assume 10% bedload
landsat_qB_fluxraw <- sum(landsat_mouth_qB$rating_flux_landsat * 1.1, na.rm=T) %>% set_units(Mt/yr)
landsat_qT_fluxraw <- sum(landsat_mouth_qT$rating_flux_landsat * 1.1, na.rm=T) %>% set_units(Mt/yr)

# Normalize 
landsat_qB_flux <- landsat_qB_fluxraw/(landsat_qB_fluxraw+landsat_qT_fluxraw)
landsat_qT_flux <- landsat_qT_fluxraw/(landsat_qB_fluxraw+landsat_qT_fluxraw)

# Get number of samples
landsat_qB_days <- nrow(landsat_mouth_qB)
landsat_qT_days <- nrow(landsat_mouth_qT)

### SENTINEL ###

SQ_top <- mouth_q %>% filter(year_chunk %in% unique(sentinel_day_avg$year_chunk)) %>%
  # group_by(year_chunk) %>%
  summarise(top = quantile(flow_m3s, 0.95, na.rm=T)) %>%
  # ungroup() %>%
  # highest Q3
  summarise(high = max(top))

# Assume 10% bedload
S_mouth_qB <- mouth_q %>% filter(year_chunk %in% unique(sentinel_day_avg$year_chunk)) %>% filter(flow_m3s < SQ_top$high)
S_mouth_qT <- mouth_q %>% filter(year_chunk %in% unique(sentinel_day_avg$year_chunk)) %>% filter(flow_m3s >= SQ_top$high)

# Assume 10% bedload
S_qB_fluxraw <- sum(S_mouth_qB$rating_flux_sentinel * 1.1, na.rm=T) %>% set_units(Mt/yr)
S_qT_fluxraw <- sum(S_mouth_qT$rating_flux_sentinel * 1.1, na.rm=T) %>% set_units(Mt/yr)

# Normalize 
S_qB_flux <- S_qB_fluxraw/(S_qB_fluxraw+S_qT_fluxraw)
S_qT_flux <- S_qT_fluxraw/(S_qB_fluxraw+S_qT_fluxraw)

# Get number of samples
S_qB_days <- nrow(S_mouth_qB)
S_qT_days <- nrow(S_mouth_qT)

### USGS ###
USGSQ_top <- mouth_q %>%
  # group_by(year_chunk) %>%
  summarise(top = quantile(flow_m3s, 0.95, na.rm=T)) %>%
  # ungroup() %>%
  # highest Q3
  summarise(high = max(top))

USGS_mouth_qB <- mouth_q %>% filter(flow_m3s < USGSQ_top$high)
USGS_mouth_qT <- mouth_q %>% filter(flow_m3s >= USGSQ_top$high)

# Assume 10% bedload
USGS_qB_fluxraw <- sum(USGS_mouth_qB$rating_flux_USGS * 1.1, na.rm=T) %>% set_units(Mt/yr)
USGS_qT_fluxraw <- sum(USGS_mouth_qT$rating_flux_USGS * 1.1, na.rm=T) %>% set_units(Mt/yr)

# Normalize 
USGS_qB_flux <- USGS_qB_fluxraw/(USGS_qB_fluxraw+USGS_qT_fluxraw)
USGS_qT_flux <- USGS_qT_fluxraw/(USGS_qB_fluxraw+USGS_qT_fluxraw)

# Get number of samples
USGS_qB_days <- nrow(USGS_mouth_qB)
USGS_qT_days <- nrow(USGS_mouth_qT)

### MAKE PLOTS ###

# Quartile dataset
storm_quartiles <- data.frame(flow=c("Landsat Bottom 95%", "Landsat Top 5%", 
                                     "Sentinel Bottom 95%", "Sentinel Top 5%",
                                     "USGS Bottom 95%", "USGS Top 5%"),
                              # Flux is proportion of flux
                              flux=c(landsat_qB_flux, landsat_qT_flux, 
                                     S_qB_flux, S_qT_flux, 
                                     USGS_qB_flux, USGS_qT_flux),
                              samp=c(landsat_qB_days, landsat_qT_days,
                                     S_qB_days, S_qT_days,
                                     USGS_qB_days, USGS_qT_days))

q_storms <- ggplot(storm_quartiles, aes(x=flow, y=flux, fill=flow)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label=sprintf("%0.1f", round(flux*100, digits = 3))), vjust=-0.2,) +
  # geom_bar(aes(y=Uflux), stat="identity", color="black") +
  # geom_text(aes(label=paste0("Number of Days ",Usamp)), vjust=-0.2) +
  theme_minimal() +
  theme(legend.position="none") +
  scale_fill_manual(values=c('#999999','#E69F00','#999999','#E69F00', '#999999','#E69F00')) +
  labs(title="Proportion of sediment flux in large storms - determined via discharge",
       x = "Discharge by Day")
       #caption = "Landsat: The top 5% of discharge days account for 30.2% of total flux (5244 [Mt/yr]).\nSentinel: The top 5% of discharge days account for 19.8% of total flux (1156 [Mt/yr]).\nUSGS: The top 5% of discharge days account for 38.3% of total flux (578 [Mt/yr]).")

q_storms

# ggsave(q_storms, 
#        filename = paste0(wd_root, '/exports/q_storms.png'),
#        width = 20, height = 8)
```

# Investigate Trapping Efficiencies of Dams 
USE EVANS CODE TO GET TRAPPING EFFICIENCY: https://github.com/evandethier/us-dam-trapping/blob/main/usgs-dam-trapping-imports/NID_dams_by_USGS_station.csv

- Trapping efficiency of 1 -> no sediment gets past
- For dam at a distance x:
- Get avg flux at start of resevoir-10 km (pre-dam) to x+10km (post-dam)
- TE = (pre-post)/(post)
- Now using SSC not flux
```{r}
load("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\riverdata.Rdata")
load(file=paste0(save_dir, "erdc_dams.RData"))

# nid_dams <- read.csv("D:\\valencig\\Thesis\\Data\\NID_dataset\\NID_raw.csv")
# nid_dams <- na.omit(nid_dams, cols=c('Longitude','Latitude'))
# # Convert to sf object
# nid_dams_sf <- st_as_sf(nid_dams, coords = c("Longitude", "Latitude")) %>% st_set_crs(4326)
# 
# # Get dams in the basin
# basin_dams <- st_intersection(nid_dams_sf, riverdata$basin$geometry)
# # Get their rows from nid_dams (geom_point needs lat and lon coordinates)
# dams <- nid_dams %>% filter(NID.ID %in% basin_dams$NID.ID)
# # One small dam (<25 feet) doesnt have a height category
# # dams[is.na(dams$NID.Height.Category),]$NID.Height.Category <- "Less than 25 feet"
# # dams$NID.Height.Category <- factor(dams$NID.Height.Category, levels=
# #                  c("Less than 25 feet",
# #                    "51-100 feet",
# #                    "Greater than 100 feet"))
# erdc_dams <- dams %>% filter(Owner.Types == "Federal")
# erdc_dams$distance_km <- 0
# erdc_dams$res_start <- 0
# erdc_dams[erdc_dams$Dam.Name == "Buford Dam",]$distance_km = 0
# erdc_dams[erdc_dams$Dam.Name == "Buford Dam",]$res_start = 0
# erdc_dams[erdc_dams$Dam.Name == "West Point Dam",]$distance_km = 240
# erdc_dams[erdc_dams$Dam.Name == "West Point Dam",]$res_start = 192
# erdc_dams[erdc_dams$Dam.Name == "Walter F. George Lock and Dam",]$distance_km = 432	
# erdc_dams[erdc_dams$Dam.Name == "Walter F. George Lock and Dam",]$res_start = 374
# erdc_dams[erdc_dams$Dam.Name == "Jim Woodruff Lock and Dam",]$distance_km = 552
# erdc_dams[erdc_dams$Dam.Name == "Jim Woodruff Lock and Dam",]$res_start = 520
# erdc_dams[erdc_dams$Dam.Name == "George W. Andrews Lock and Dam",]$distance_km = 478
# erdc_dams[erdc_dams$Dam.Name == "George W. Andrews Lock and Dam",]$res_start = 476
# erdc_dams$V_m3 <- erdc_dams$Volume..Cubic.Yards. * 0.7646

trapping_eff <- function(dist, res_start, years,  df){
  if (dist == 0) { # Buford Dam
    last_dist <- 0
    next_dist <- 10
  } else {
    last_dist <- res_start-10
    next_dist <- dist+10
  }
  # Don't look at dist --> where dam is
  pre <<- df %>% filter(distance_km >= last_dist & distance_km <= res_start) %>% filter(year_chunk == years)
  post <<- df %>% filter(distance_km <= next_dist & distance_km > dist) %>% filter(year_chunk == years)
  pre_ssc <- mean(pre$SSC_mgL, na.rm=T)
  post_ssc <- mean(post$SSC_mgL, na.rm=T)
  dssc <- post_ssc-pre_ssc
  te <- -dssc / pre_ssc
  l <- list("pre_ssc"=pre_ssc, "post_ssc"=post_ssc, "te"=te, "dssc"=dssc)
  return(l)
}

# trapping_eff(erdc_dams[2,], landsat_day_avg)

# Trapping efficiency for the whole time period
te_whole <- data.frame(Dam.Name=erdc_dams$Dam.Name, 
                       distance_km=erdc_dams$distance_km,
                       res_start=erdc_dams$res_start,
                       Volume_m3 = erdc_dams$V_m3,
                       years = rep(unique(landsat_day_avg$year_chunk), each=5))
# Filter to west point for testing
# te_whole <- te_whole %>% filter(Dam.Name == "West Point Dam")

# 1440 kg/m^3 from Trimble et al.
Mt_per_m3 <- 1.440*10**(-6)

# Filter to 2015-present year chunk
# sentinel_day_avg <- sentinel_day_avg %>% filter(year_chunk %in% c("[2015,2020)","[2020,2025)"))
# landsat_day_avg <- landsat_day_avg %>% filter(year_chunk %in% c("[2015,2020)","[2020,2025)"))

te_landsat <- te_whole %>% rowwise() %>% 
  mutate(te = trapping_eff(distance_km, res_start, years, landsat_day_avg)$te) %>%
  mutate(pre_ssc_mgL = trapping_eff(distance_km, res_start, years, landsat_day_avg)$pre_ssc) %>%
  mutate(post_ssc_mgL = trapping_eff(distance_km, res_start, years, landsat_day_avg)$post_ssc) %>%
  mutate(dssc_mgL = trapping_eff(distance_km, res_start, years, landsat_day_avg)$dssc)
  # mutate(time2fill_yr = (1.440*10**(-6) * Volume_m3) / -dflux_Mt)

te_sentinel <- te_whole %>% rowwise() %>% 
  mutate(te = trapping_eff(distance_km, res_start, years, sentinel_day_avg)$te) %>%
  mutate(pre_ssc_mgL = trapping_eff(distance_km, res_start, years, sentinel_day_avg)$pre_ssc) %>%
  mutate(post_ssc_mgL = trapping_eff(distance_km, res_start, years, sentinel_day_avg)$post_ssc) %>%
  mutate(dssc_mgL = trapping_eff(distance_km, res_start, years, sentinel_day_avg)$dssc)
  # mutate(time2fill_yr = (1.440*10**(-6) * Volume_m3) / -dflux_Mt)

# Investigate data
# ggplot() + 
#   stat_summary(
#     data = pre,
#     aes(x = distance_km, y = SSC_mgL, color='Landsat'),
#     geom = 'line', 
#     fun = 'mean',
#     inherit.aes = FALSE,
#     # show.legend=TRUE,
#     na.rm=T) +
#   geom_vline(
#     data = te_landsat,
#     aes(xintercept = res_start, color='res_start')) +
#   geom_vline(
#     data = te_landsat,
#     aes(xintercept = distance_km, color='dam'))

# Plot Trapping efficiencies over time
ggplot() +
  geom_line(data = te_landsat %>% filter(Dam.Name %in% c("West Point Dam","Walter F. George Lock and Dam")), 
            aes(x=years,  y = te, color=Dam.Name, group=Dam.Name)) +
  geom_line(data = te_sentinel %>% filter(Dam.Name %in% c("West Point Dam","Walter F. George Lock and Dam")), 
            aes(x=years,  y = te, color=Dam.Name, group=Dam.Name), na.rm=T)

# Investigate Walter F. George and West Point
print(paste("Walter F. George (Landsat)",mean(te_landsat[te_landsat$Dam.Name == "Walter F. George Lock and Dam",]$te)))
print(paste("Walter F. George (Sentinel-2)",mean(te_sentinel[te_sentinel$Dam.Name == "Walter F. George Lock and Dam",]$te, na.rm=T)))
print(paste("West Point Dam (Landsat)",mean(te_landsat[te_landsat$Dam.Name == "West Point Dam",]$te)))
print(paste("West Point Dam (Sentinel-2)",mean(te_sentinel[te_sentinel$Dam.Name == "West Point Dam",]$te, na.rm=T)))
# Save variables
# save(erdc_dams, file=paste0(save_dir, "erdc_dams.RData"))
# load(file=paste0(save_dir, "erdc_dams.RData"))
```

# Extract one image above/below west point dam
```{r}
sentinel_dataRating <- data.table(
  read.csv("D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_scaling\\train_ratingCurve_clust_ratingCurve\\regression\\lasso_full_bands_ratingclust3.csv"))[,':='(
  SSC_mgL_R = pred_SSC_mgL, # SSC (mg/L) Rating curve --> data trained and clusters using rating curve derived SSC
  year = format(as.Date(date), "%Y")
)]
west_point_img <- sentinel_dataRating %>%
  filter(product_id =="S2A_MSIL2A_20190306T162131_N0211_R040_T16SFB_20190306T204355") %>%
  dplyr::select(date, distance_km, pred_SSC_mgL)
west_point <- erdc_dams %>% filter(Dam.Name == "West Point Dam")

pre_res <- west_point_img %>% filter(distance_km >= west_point$res_start-10 & distance_km <= west_point$res_start)
post_res <- west_point_img %>% filter(distance_km <= west_point$distance_km+10 & distance_km > west_point$distance_km)

pre_ssc <- mean(pre_res$pred_SSC_mgL)
post_ssc <- mean(post_res$pred_SSC_mgL)


```