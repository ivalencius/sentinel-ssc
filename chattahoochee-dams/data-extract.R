######################### Load packages ########################################

library(dataRetrieval)
# library(sp)
library(sf)
# library(maptools)
# library(magrittr)
# library(maps)
# library(rgdal)
library(ggplot2)
library(usmap)
# library(stringr)
library(dplyr)
# library(patchwork)
library(dams)
library(data.table)
# library(nhdplusTools)
# library(riverdist)
# library(GISTools)
# library(rgeos)
# library(tidyr)
# library(qpcR)
# library(scales)
# library(vioplot)
# library(raster)
# library(units)
# library(extrafont)
# loadfonts(device = "win")
# library(latex2exp) # use TeX() to render latex
# library(paletteer) # color palettes

####################### Set up paramters #######################################

# Start point of the river - WGS84 SRS ordered lng/lat (X,Y) (Decimal)
startPoint <-  c(-85.03114856911581, 29.84561518630128)
# startPoint <- c(-85, 30)
# Distance to look upstream (set to be very large for whole transect)
dist2Look <- 686 # must be an integer (no expressions)
# For extracting dam information
damsOfInterest <- c(
    "Buford Dam",
    "West Point Dam",
    "Walter F. George Lock and Dam",
    "Jim Woodruff Lock and Dam",
    "George W. Andrews Lock and Dam"
)

# Working Directory
wd_root <- "C:/Users/ilanv/Desktop/sentinel-ssc/chattahoochee-dams"
setwd(wd_root)

# Exports folder
wd_exports <- paste0(wd_root, "/exports/")
# Figures folder
wd_figure <- paste0(wd_exports, "figures/")

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

######################### Extract DAMS  #######################################

damFileName <- paste0(wd_exports, "NID_DAMS.csv")

# Download National Inventory of Dams
tryCatch({
    nid_dams <- read.csv(damFileName, skip = 1)
    }, error = function(cond) {
        download.file("https://nid.usace.army.mil/api/nation/csv", damFileName)
    }
)
nid_dams <- read.csv(damFileName, skip = 1)

# Get dams in the basin --> now can Just filter by name
# # Remove dams with NA lon/lat
# nid_dams <- na.omit(nid_dams, cols = c("Longitude", "Latitude"))
# # Convert to sf object
# nid_dams_sf <- st_as_sf(nid_dams, coords = c("Longitude", "Latitude")) %>%
#     st_set_crs(4326)
# basin_dams <- st_intersection(nid_dams_sf, riverData$basin$geometry)
# # Get their rows from nid_dams (geom_point needs lat and lon coordinates)
# dams <- nid_dams %>% filter(NID.ID %in% basin_dams$NID.ID)
# erdc_dams <- dams %>% filter(Owner.Types == "Federal")

dams <- nid_dams %>% filter(
    Dam.Name %chin% damsOfInterest ||
    Other.Names %chin% damsOfInterest)

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


########################## Save data  ##########################################
write.csv(SSC,
    file = paste0(wd_exports, "CHATTAHOOCHEE_SSC.csv"), row.names = FALSE)
write.csv(Q,
    file = paste0(wd_exports, "CHATTAHOOCHEE_Q.csv"), row.names = FALSE)