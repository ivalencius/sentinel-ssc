######################### Load packages ########################################

library(dataRetrieval)
library(sf)
library(sp)
library(ggplot2)
library(usmap)
library(dplyr)
library(dams)
library(data.table)
library(riverdist)

####################### Set up paramters #######################################

# Start point of the river - WGS84 SRS ordered lng/lat (X,Y) (Decimal)
startPoint <-  c(-85.03114856911581, 29.84561518630128)
# startPoint <- c(-85, 30)
# Distance to look upstream (set to be very large for whole transect)
dist2Look <- 686 # must be an integer (no expressions)

# Working Directory
wd_root <- "C:/Users/ilanv/Desktop/sentinel-ssc/chattahoochee-dams"
setwd(wd_root)

# Exports folder
wd_exports <- paste0(wd_root, "/exports/")
# Figures folder
wd_figure <- paste0(wd_exports, "figures/")

# NID dataset
damFileName <- paste0(wd_exports, "NID.csv")
# For extracting dam information
damsOfInterest <- c(
    "Buford Dam",
    "West Point Dam",
    "Walter F. George Lock and Dam",
    "Jim Woodruff Lock and Dam",
    "George W. Andrews Lock and Dam"
)

# Create folders within root directory
export_folder_paths <- c(wd_exports, wd_figure)
for (i in seq_along(export_folder_paths)){
  path_sel <- export_folder_paths[i]
  if (!dir.exists(path_sel)) { dir.create(path_sel) }
}

####################### Extract Site Info ######################################

summarize.nldi = function(input){
  data.frame(name = names(input), 
             class = sapply(input, class)[1], 
             row.names = NULL) %>% 
    mutate(feature_count = ifelse(class == "sf", sapply(input, nrow), 
                                  sapply(input, length)))
}

riverData <- findNLDI(
    location = startPoint,
    # wqp = lastStation,
    nav = "UM",
    find = c("nwis", "WQP", "basin", "flowlines"),
    distance_km = dist2Look
    )

basin_area <- st_area(riverData$basin) # m^2

summarize.nldi(riverData)

# Quick sanity check of station location
ggplot() +
    labs(title = "USGS Station Locations") +
    geom_sf(data = riverData$basin) + 
    geom_sf(data = riverData$UM_flowlines, col = "blue") + 
    geom_point(data = riverData$UM_WQP,
        aes(x = X, y = Y, color = "WQP"),
        size = 1) +
    geom_point(data = riverData$UM_nwissite,
        aes(x = X, y = Y, color = "NWIS"),
        size = 1) +
    theme_bw() +
    labs(
        title = "USGS Station Locations",
        x = "Lon",
        y = "Lat"
    )

# Extract site names

# Pads a character string with leading 0s (for USGS codes)
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
# get all WQP sites
mc_names <- paste0("USGS-", lapply(riverData$UM_WQP$comid, pad0))
# merge with nwis sites
mc_names <- unique(append(mc_names, riverData$UM_nwissite$identifier))
mc_numbers <- gsub("USGS-", "", mc_names)

####################### Extract SSC and Q  #####################################

# Extract SSC data
SSC <- data.table(readWQPdata(
    siteNumbers = mc_names,
    parameterCd = "80154" # SSC (mg/L)
    ))[, ":=" (
        agency_cd = "USGS",
        site_no = unlist(
            strsplit(
                MonitoringLocationIdentifier, split='-')
                )[c(FALSE, TRUE)],
        sample_dt = ActivityStartDate,
        SSC_mgL = ResultMeasureValue,
        sample_depth_m = ActivityDepthHeightMeasure.MeasureValue * 0.3048,
        sample_method = SampleCollectionMethod.MethodName)][,.(
            agency_cd, site_no, sample_dt, SSC_mgL,
            sample_depth_m, sample_method
        )]

Q <- data.table(readNWISdata(
    sites = mc_numbers,
    service = "dv",
    parameterCd = "00060",
    statCd = "00003",
    startDate = as.Date("1800-01-01")
    ))[, ":=" (
    agency_cd = "USGS",
    sample_dt = dateTime,
    Q_m3s = X_00060_00003 * 0.02832 # ft3/s to m3/s
  )][, .(agency_cd, site_no, sample_dt, Q_m3s)]

# Get all site _names 
USGS_sites <- lapply(Q$site_no, pad0)
USGS_sites <- append(USGS_sites, lapply(SSC$site_no, pad0))
USGS_sites <- unique(USGS_sites)

# Extract all site information
USGS_siteinfo <- readNWISsite(USGS_sites)

S1 <- as(riverData$UM_flowlines, 'Spatial')
nldi_proj4 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
AKalbers <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154
    +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
# Fix line2network function
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
river_center <- line2network(S1, reproject=AKalbers)

# Apply projection to USGS sites then convert into meter projection
stations <- data.frame(
  site_no = USGS_siteinfo$site_no,
  lat = as.numeric(USGS_siteinfo$dec_lat_va),
  lon = as.numeric(USGS_siteinfo$dec_long_va)
)

coordinates(stations) <- ~lon + lat
proj4string(stations) <- '+init=epsg:4326'
stations <- spTransform(stations, CRS(AKalbers))
stations <- as.data.frame(stations)
riverVert <- xy2segvert(x=stations$lon, y=stations$lat, rivers=river_center)

# Determine distance from start of channel for USGS sites
distances_km<- list()
for (i in 1:length(riverVert$vert)){
  # Distance to mouth of river
  dist_m <- riverdistance(startseg = riverVert$seg[i], startvert = riverVert$vert[i],
                            endvert = 1, endseg = 1, rivers=river_center)
  # Convert from degrees to km
  distances_km <- append(distances_km, dist2Look-dist_m/1000)
}
riverVert$dist_downstream_km <- distances_km

# Merge station info
stations$dist2river_m <- riverVert$snapdist
stations$dist_downstream_km <- riverVert$dist_downstream_km

# Some stations are errors -> not in the basin
stations <- stations %>% filter(dist2river_m <= 1000)
stations <- apply(stations,2,as.character)

# Apply stats to SSC and Q data
SSC <- left_join(SSC, stations) %>% filter(is.finite(dist2river_m))
Q <- left_join(Q, stations) %>% filter(is.finite(dist2river_m))

######################### Extract DAMS  ########################################

# Download National Inventory of Dams
# tryCatch({
#     nid_dams <- read.csv(damFileName, skip = 1)
#     }, error = function(cond) {
#         download.file("https://nid.usace.army.mil/api/nation/csv", damFileName)
#     }
# )
nid_dams <- read.csv(damFileName)

# Get dams in the basin --> now can Just filter by name
# # Remove dams with NA lon/lat
# nid_dams <- na.omit(nid_dams, cols = c("Longitude", "Latitude"))
# # Convert to sf object
# nid_dams_sf <- st_as_sf(nid_dams, coords = c("Longitude", "Latitude")) %>%
#     st_set_crs(4326)
# basin_dams <- st_intersection(nid_dams_sf, riverData$basin$geometry)
# # Get their rows from nid_dams (geom_point needs lat and lon coordinates)
# dams <- nid_dams %>% filter(NID.ID %in% basin_dams$NID.ID)
# dams <- dams %>% filter(Owner.Types == "Federal")

dams <- nid_dams %>% filter(
    Dam.Name %chin% damsOfInterest)

dams$distance_km <- 0
dams$res_start <- 0

dams[dams$Dam.Name == "Buford Dam", ]$distance_km <- 0
dams[dams$Dam.Name == "Buford Dam", ]$res_start <- 0
dams[dams$Dam.Name == "West Point Dam", ]$distance_km <- 240
dams[dams$Dam.Name == "West Point Dam", ]$res_start <- 192
dams[dams$Dam.Name == "Walter F. George Lock and Dam", ]$distance_km <- 432
dams[dams$Dam.Name == "Walter F. George Lock and Dam", ]$res_start <- 374
dams[dams$Dam.Name == "Jim Woodruff Lock and Dam", ]$distance_km <- 552
dams[dams$Dam.Name == "Jim Woodruff Lock and Dam", ]$res_start <- 520
dams[dams$Dam.Name == "George W. Andrews Lock and Dam", ]$distance_km <- 478
dams[dams$Dam.Name == "George W. Andrews Lock and Dam", ]$res_start <- 476
dams$V_m3 <- dams$Volume..Cubic.Yards. * 0.7646


########################## Save data  ##########################################
# Save Data
write.csv(apply(SSC,2,as.character),
    file = paste0(wd_exports, "PROCESSED_CHATTAHOOCHEE_SSC.csv"), row.names = FALSE)
write.csv(apply(Q,2,as.character),
    file = paste0(wd_exports, "PROCESSED_CHATTAHOOCHEE_Q.csv"), row.names = FALSE)
write.csv(stations,
    file = paste0(wd_exports, "CHATTAHOOCHEE_STATIONINFO.csv"), row.names = FALSE)
write.csv(dams,
    file = paste0(wd_exports, "ERDC_CHATTAHOOCHEE_DAMS.csv"), row.names = FALSE)
