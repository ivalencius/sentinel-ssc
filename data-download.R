#### LIBRARY IMPORTS ####
library(dataRetrieval)
library(ggplot2)
library(data.table)
library(USAboundaries)
library(sp)
library(raster)
library(dplyr)

#### SET DIRECTORIES ####
startDate <- "2017-03-28" # Start of Sentinel-2 Record in GEE
endDate <- "2024-01-01"

print('SETTING UP IMPORTS')
# Set root directory
wd_root <- "/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc"
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
us_states_spatial <- us_states_spatial[,c('state_abbr')] # only select column with province name

# fortify state shapefile
us_ca <- fortify(us_states_spatial)

#### IMPORT IN SITU SSC DATA ####

# From USGS
print("IMPORTING IN SITU DATA")
cat("\t","-> Start Date:", startDate, "\n")
cat("\t","-> End Date:", endDate, "\n")

#ssc_codes <- c('00530','80154')
ssc_codes <- c("80154")
cat("\t","-> Downloading USGS SSC Data", "\n")
wqp_insitu_raw <- data.table(readWQPdata(
  parameterCd = ssc_codes,
  startDate=startDate,
  endDate=endDate))[,":=" (
    # AGENCY_CD = "USGS",
    SITE_NO = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
    SAMPLE_DT = ActivityStartDate,
    SSC_MGL = ResultMeasureValue,
    SAMPLE_DEPTH_M = ActivityDepthHeightMeasure.MeasureValue * 0.3048,
    SAMPLE_METHOD = SampleCollectionMethod.MethodName
  )][,.(SITE_NO, SAMPLE_DT, SSC_MGL, SAMPLE_DEPTH_M, SAMPLE_METHOD)]

# Need to use WQP data as well
nwis_insitu_raw <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("AGENCY_CD", "SITE_NO", "SAMPLE_DT", "SSC_MGL"))))

for (state in us_states_spatial$state_abbr){
  nwis_data <- data.table(readNWISdata(stateCd=state,
  parameterCd = ssc_codes, 
  startDt=startDate,
  endDt=endDate))
  if(nrow(nwis_data)==0){} # Do nothing if no data returned
  else{
    nwis_data_clean <- nwis_data[, ":=" (
    # AGENCY_CD = agency_cd,
    SITE_NO = site_no,
    SAMPLE_DT = dateTime,
    SSC_MGL = X_80154_00003
    )][,.(SITE_NO, SAMPLE_DT, SSC_MGL)]
    nwis_insitu_raw <- rbind(nwis_insitu_raw, nwis_data_clean)
  }
}
# Filter out nonsense data
# nwis_insitu_clean<- nwis_insitu_raw %>% filter(is.finite(SSC_MGL))

# Merge nwis and usgs
usgs_insitu_raw <- bind_rows(data.frame(wqp_insitu_raw), data.frame(nwis_insitu_raw))

# Filter out nonsense data
usgs_insitu_clean <- usgs_insitu_raw %>% filter(is.finite(SSC_MGL))

# Convert to all Caps
usgs_insitu_upper <- mutate_all(usgs_insitu_clean, toupper)

# For WQP need USGS- prefix
usgs_stations <- paste0("USGS-", unique(usgs_insitu_upper$SITE_NO))

# Get properties of usgs WQP stations
cat("\t","-> Downloading USGS Station Metadata", "\n")
station_info <- data.table(whatWQPsites(siteid = usgs_stations))[,":="(
  STATION_NM = MonitoringLocationName,
  # Now getting drainage area in GEE
  # DRAINAGE_AREA_KM2 = DrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
  # CONTRIBUTING_DRAINAGE_AREA_KM2 = ContributingDrainageAreaMeasure.MeasureValue*(2.58999/1.00000073), # sq mi -> sq km
  LAT = LatitudeMeasure,
  LON = LongitudeMeasure,
  # ELEVATION_M = VerticalMeasure.MeasureValue *0.3048, # ft -> m # nolint
  SITE_NO = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
  SITE_TYPE = MonitoringLocationTypeName,
  HUC = HUCEightDigitCode
)][, .(SITE_NO, STATION_NM, LAT, LON, SITE_TYPE, HUC)]
# HUC8 code for use with https://developers.google.com/earth-engine/datasets/catalog/USGS_WBD_2017_HUC08#table-schema

#is_w_station <- in-situ data harmonized with station data
is_w_stations <- left_join(usgs_insitu_clean, mutate_all(station_info, toupper), by=("SITE_NO"="SITE_NO"))

# Keep only stations with relevant data
types <- c("STREAM", "LAKE, RESERVOIR, IMPOUNDMENT", "ESTUARY", "RIVER/STREAM", "STREAM: CANAL")
print(paste0("Excluded columns :", paste0(setdiff(unique(is_w_stations$SITE_TYPE), types), collapse = ", ")))
water_only <- is_w_stations %>% filter(SITE_TYPE %in% types)

################################ CHECK THAT ALL STATIONS HAVE DRAINAGE DATA

# Reduce in-situ SSC measurements to one sample per station per day
oneDaySSC <- water_only %>% 
  group_by(SITE_NO, SAMPLE_DT) %>%
  summarise(across(where(is.numeric), mean)) %>%
  # slice_head() %>%
  ungroup() %>%
  left_join(water_only %>% select(-c(SAMPLE_DEPTH_M, SSC_MGL)),
  by = c("SITE_NO"="SITE_NO", "SAMPLE_DT"="SAMPLE_DT")) %>% 
  distinct() # Sumarize won't remove duplicate rows

# View(oneDaySSC %>% filter(SITE_NO == "453250122494501")) # Check that it worked

# Lat Lon are currently characters --> convert to numeric
oneDaySSC$LAT <- as.numeric(oneDaySSC$LAT)
oneDaySSC$LON <- as.numeric(oneDaySSC$LON)

# Filter SSC data that is erroneous (some streams have zero SSC when they are dry)
oneDaySSC <- oneDaySSC %>% filter(SSC_MGL > 0)
oneDaySSC$LOG10_SSC_MGL <- log10(oneDaySSC$SSC_MGL)
hist(oneDaySSC$LOG10_SSC_MGL)

# Save raw data
# write.csv(usgs_insitu_raw, file = paste(wd_exports, "NWIS_WQP_SSC.csv", sep=""), row.names=F)
# write.csv(is_w_stations, file = paste(wd_exports, "NWIS_WQP_SSC_STATIONINFO.csv", sep=""), row.names=F)
write.csv(oneDaySSC, file = paste0(wd_exports, "ONEDAYSAMPLE_SSC.csv", sep=""), row.names=F)

# For loading SSC
oneDaySSC <- read.csv(file = paste0(wd_exports, "ONEDAYSAMPLE_SSC.csv"))
oneDaySSC$SITE_NO <- as.character(oneDaySSC$SITE_NO)
oneDaySSC$SAMPLE_DT <- as.Date(oneDaySSC$SAMPLE_DT)
oneDaySSC$SSC_MGL <- as.numeric(oneDaySSC$SSC_MGL)
oneDaySSC$LOG10_SSC_MGL <- as.numeric(oneDaySSC$LOG10_SSC_MGL)
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
oneDaySSC$SITE_NO <- unlist(lapply(oneDaySSC$SITE_NO, pad0))


# Save station data for GEE
gee_station_data <- oneDaySSC %>% select(SITE_NO, STATION_NM, HUC, SITE_TYPE, LAT, LON) %>% distinct()
# Also save station info as a shape file
coordinates(gee_station_data) <- ~LON+LAT
proj4string(gee_station_data) <- projection
shapefile(gee_station_data, paste0(wd_gee, "STATION_DATA.shp"), overwrite=T)

### Create SSC vs. Q Rating Curve ###

### IMPORT IN SITU DISCHARGE DATA ###
# Import data and make rating curves for all stations
# --> Takes a LONG time
cat("\t","-> Downloading USGS Discharge Data", "\n")
discharge_data <- data.table(
  readNWISdata(
    sites = unique(oneDaySSC$SITE_NO),
    service = "dv",
    parameterCd = "00060",
    statCd="00003",
    startDt = startDate,
    endDt = endDate
  ))[,":=" (
    SITE_NO = site_no,
    SAMPLE_DT = dateTime,
    Q_M3S = X_00060_00003 * 0.02832 # ft3/s to m3/s
  )][,.(SITE_NO, SAMPLE_DT, Q_M3S)]

# Convert all discharges to positive
discharge_pos <- discharge_data %>% 
mutate(Q_M3S = abs(as.numeric(Q_M3S))) %>% 
filter(Q_M3S > 0)

# Save and load variables for quick execution
# save(discharge_data, file=paste0(wd_root, '/tmp_vars/discharge_data.RData'))
# load(paste0(wd_root, '/tmp_vars/discharge_data.RData'))
write.csv(discharge_pos, file = paste0(wd_exports, "NWIS_POSITIVE_Q.csv"), row.names = FALSE)

# For loading discharge data
discharge_pos <- read.csv(file = paste0(wd_exports, "NWIS_POSITIVE_Q.csv"))
discharge_pos$SITE_NO <- as.character(discharge_pos$SITE_NO)
discharge_pos$Q_M3S <- as.numeric(discharge_pos$Q_M3S)
discharge_pos$SAMPLE_DT <- as.Date(discharge_pos$SAMPLE_DT)

SSC_with_Q <- left_join(oneDaySSC, discharge_pos,
by=c("SITE_NO"="SITE_NO", "SAMPLE_DT"="SAMPLE_DT"))
Q_no_SSC <- anti_join(discharge_pos, SSC_with_Q, 
by=c("SITE_NO"="SITE_NO", "SAMPLE_DT"="SAMPLE_DT", "Q_M3S"="Q_M3S"))
SSCQ <- bind_rows(SSC_with_Q, Q_no_SSC)

# NEED TO HAVE POSITIVE FLUXES THE WHOLE TIME
ratingSSCQ <- setDT(copy(SSCQ))[
  ,':='(
    LOG10_Q_M3S = log10(Q_M3S),
    LOG10_SSC_FLUX_MTYR = log10(Q_M3S * SSC_MGL * 3.10585 * 10**-5),
    SSC_FLUX_MTYR = Q_M3S * SSC_MGL * 3.10585 * 10**-5,
    R2 = 0.0,
    RATINGCURVE = FALSE
    )]

# Get site numbers for stations with Q data and >=10 SSC observations
station_nums <- ratingSSCQ %>%
  filter(is.finite(LOG10_SSC_FLUX_MTYR)) %>%
  group_by(SITE_NO) %>%
  summarise(samples = n()) %>%
  filter(samples >= 10) %>%
  ungroup() %>%
  pull("SITE_NO")

lm_eqn <- function(m){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,
  list(a = format(unname(coef(m)[1]), digits = 2),
  b = format(unname(coef(m)[2]), digits = 2),
  r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
}

print("GENERATING FLUX V. Q RATING CURVES")
pb <- txtProgressBar(0, length(station_nums), style = 3)
for (i in seq_along(station_nums)){
# for (i in 1) {
  setTxtProgressBar(pb, i)

  # Extract data for the site
  station_data <- ratingSSCQ %>% 
    filter(SITE_NO == station_nums[i])
  
  # Extract data for the site with flux measurements
  station_flux <- station_data %>%
    filter(is.finite(LOG10_SSC_FLUX_MTYR))

  # Generate regression
  reg <- lm(LOG10_SSC_FLUX_MTYR ~ LOG10_Q_M3S, data = station_flux)

  R2 <- summary(reg)$r.squared
  if (R2 >= 0.80) {
    # Get indexes of samples where there is no SSC (i.e. no in-situ sample)
    goodIDs <- (ratingSSCQ$SITE_NO == station_nums[i]) & is.na(ratingSSCQ$SSC_MGL)
    # print(length(goodIDs[goodIDs == TRUE]))
    # Apply regression
    ssc_flux <- predict.lm(reg, newdata = ratingSSCQ[goodIDs])
    ratingSSCQ[goodIDs]$LOG10_SSC_FLUX_MTYR <- ssc_flux
    # Get secondary metrics
    ratingSSCQ[goodIDs]$SSC_FLUX_MTYR <- 10^ssc_flux
    ratingSSCQ[goodIDs]$SSC_MGL <- (10^ssc_flux) / (ratingSSCQ[goodIDs]$Q_M3S * 3.10585 * 10**-5)
    ratingSSCQ[goodIDs]$LOG10_SSC_MGL <- log10(ratingSSCQ[goodIDs]$SSC_MGL)
    ratingSSCQ[goodIDs]$R2 <- R2
    ratingSSCQ[goodIDs]$RATINGCURVE <- TRUE

    # Export regression
    reg_plot <- ggplot(station_flux, aes(x = LOG10_Q_M3S, y = LOG10_SSC_FLUX_MTYR)) +
      geom_point(color = 'green', ) +
      geom_smooth(method='lm', formula= y~x, color='black',
      aes(fill = '1 standard error')) +
      annotate(geom='text',label = lm_eqn(reg), parse = TRUE, x = -Inf, y = Inf, 
      hjust = -0.2, vjust = 2, family="JetBrains Mono NL") +
      theme_bw() +
      scale_fill_manual(values = c('gray'), name = "1 standard Error")  +
      theme(text = element_text(family = "JetBrains Mono NL")) +
      theme(legend.position = c(0.8, 0.2)) +
      theme(legend.title = element_blank()) +
      labs(title=paste0('Flux vs Q rating Curve for Station USGS-',
      station_nums[i]), subtitle=paste0(startDate, " to ", endDate))
    ggsave(reg_plot, filename = paste0(wd_rating, paste0('USGS-',
    station_nums[i],'.png')), width = 6, height = 6)
  }
}

# NUMBER OF TRAINING SAMPLES: 162084
# Check number of modified samples -->
regressed <-ratingSSCQ %>% filter(is.finite(SSC_MGL))
# print(nrow(ratingSSCQ %>% filter(is.finite(SSC_MGL))))

print("Size of raw in-situ data")
print(nrow(oneDaySSC))

print("Size of rating-curve in-situ data")
print(nrow(regressed))

write.csv(regressed, file=paste0(wd_exports, "RATINGCURVE_SSC.csv"), row.names = FALSE)