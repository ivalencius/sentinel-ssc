############################### TO DO ##########################################
# Redo plots with adjusted clusters

################################################################################
# Replacement for ssc_retrieval.R
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
library(zoo) # moving averages
library(cowplot)
library(showtext)
library(reshape2)
library(Hmisc)
library(glmnet)
library(latex2exp)
font_add("Avenir", regular = "C:/Users/ilanv/AppData/Local/Microsoft/Windows/Fonts/AvenirNextLTPro-Regular.otf")
# font_add_google()
showtext_auto()


# NEED TO LOAD LANDASAT D:/valencig/Thesis/chattahoochee-dams/chattahoochee-dams-exports/SSC_pred.csv

wd_root <- "C:/Users/ilanv/Desktop/sentinel-ssc/chattahoochee-dams"
setwd(wd_root)
# Exports folder
wd_exports <- paste0(wd_root, "/exports/")
# Figures folder
wd_figure <- paste0(wd_exports, "figures/")

sentinel <- read.csv(file="C:\\Users\\ilanv\\Desktop\\sentinel-ssc\\exports\\transects\\SENINEL2_20m_SSC.csv")
landsat <- read.csv(file="dethier_ls57_chattahoochee.csv")
USGS_SSC <- read.csv(file=paste0(wd_exports, "PROCESSED_CHATTAHOOCHEE_SSC.csv"))
USGS_Q <- read.csv(file=paste0(wd_exports, "PROCESSED_CHATTAHOOCHEE_Q.csv"))
# Note -> distances in sites are off by 
USGS_SITES <- read.csv(file=paste0(wd_exports, "CHATTAHOOCHEE_STATIONINFO.csv"))
DAMS <- read.csv(file=paste0(wd_exports, "ERDC_CHATTAHOOCHEE_DAMS.csv"))

landsat_flux <- read.csv(file=paste0(wd_exports, "LANDSAT_FLUX.csv"))
sentinel_flux <- read.csv(file=paste0(wd_exports, "SENTINEL_FLUX.csv"))

############################ PROCESS DATA ######################################

# Commented out sections already done --> preprocessing to landsat_flux and sentinel_flux
# Reduce SSC samples to one per day
# sentinel <- sentinel %>%
#   group_by(distance_km, date) %>%
#   dplyr::summarize(SSC_mgL = mean(SSC_mgL))

# landsat$date <- landsat$sample_dt
# landsat <- landsat %>%
#   group_by(distance_km, date) %>%
#   summarize(SSC_mgL = mean(SSC_mgL))

USGS_SSC <- USGS_SSC %>%
  group_by(site_no, sample_dt) %>%
  dplyr::summarize(SSC_mgL = mean(SSC_mgL))

# Remove outliers -> justified because very few measurements > 300
# maxSSC <- quantile(sentinel$SSC_mgL, probs = 0.9998)
# sentinel_outliers <- sentinel %>% filter(SSC_mgL > maxSSC)
# sentinel <- sentinel %>% filter(SSC_mgL <= maxSSC)
# maxLandsatSSC <- quantile(landsat$SSC_mgL, probs = 0.9998, na.rm=T)
# landsat_outliers <- landsat %>% filter(SSC_mgL > maxLandsatSSC)
# landsat <- landsat %>% filter(SSC_mgL <= maxLandsatSSC)

USGS_SSC$year <- as.numeric(format(as.Date(USGS_SSC$sample_dt), "%Y"))
USGS_Q$year <- as.numeric(format(as.Date(USGS_Q$sample_dt), "%Y"))
# sentinel$year <- as.numeric(format(as.Date(sentinel$date), "%Y"))
# landsat$year <- as.numeric(format(as.Date(landsat$date), "%Y"))

# Split into five year increments
USGS_SSC$year_chunk <- cut(USGS_SSC$year, seq(from = 1895, to = 2026, by = 5), right=FALSE)
USGS_Q$year_chunk <- cut(USGS_Q$year, seq(from = 1895, to = 2026, by = 5), right=FALSE)
# sentinel$year_chunk <- cut(sentinel$year, seq(from = 1895, to = 2026, by = 5), right=FALSE)
# landsat$year_chunk <- cut(landsat$year, seq(from = 1895, to = 2026, by = 5), right=FALSE)

# Merge with distance data
USGS_SSC <- left_join(USGS_SSC, USGS_SITES)
USGS_Q <- left_join(USGS_Q, USGS_SITES)

########################### Q vs. Drainage #####################################

sat_dists <- read.csv(file="C:\\Users\\ilanv\\Desktop\\sentinel-ssc\\exports\\transects\\SENINEL2_20m_SSC.csv")
sat_dists <- sat_dists %>% dplyr::select(distance_km, lat, lon) %>% unique()
sat_dists$drainage_m2 <- NA

for (i in 1:nrow(sat_dists)) {
  startPoint <- c(sat_dists[i, ]$lon, sat_dists[i, ]$lat)
  riverData <- findNLDI(
    location = startPoint,
    nav = "UM",
    find = c("basin"),
    distance_km = 1000
    )
  basin_area <- st_area(riverData$basin) # m^2
  sat_dists[i,]$drainage_m2 <- basin_area
}

# Go through and manually fix spots
ggplot(sat_dists, aes(x=distance_km, y=drainage_m2))+geom_line()
bad_dists <- sat_dists %>% filter(drainage_m2 < 10^8) %>% dplyr::select(distance_km)
# If there are two bad spots as neighbors, manually set first to be drainage before it
sat_dists[sat_dists$distance_km == 252, ]$drainage_m2 <- 9414080727.7062
sat_dists[sat_dists$distance_km == 418, ]$drainage_m2 <- 18054187914.0187
sat_dists[sat_dists$distance_km == 552, ]$drainage_m2 <- 21886962240.2116
sat_dists[sat_dists$distance_km == 400, ]$drainage_m2 <- (17425594951.7694 + 17338039822.1749)/2
sat_dists[sat_dists$distance_km == 438, ]$drainage_m2 <- (19581177883.1112 + 19629707695.549)/2
sat_dists[sat_dists$distance_km == 150, ]$drainage_m2 <- (6394035318.2442 + 6266584278.402)/2
sat_dists[sat_dists$distance_km == 456, ]$drainage_m2 <- (20056978385.8478 + 20010834637.6835)/2
sat_dists[sat_dists$distance_km == 420, ]$drainage_m2 <- (19142774242.5028 + 18054187914.0187)/2
sat_dists[sat_dists$distance_km == 424, ]$drainage_m2 <- (19199048948.6904 + 19142774242.5028)/2
sat_dists[sat_dists$distance_km == 256, ]$drainage_m2 <- (9414080727.7062 + 10051554840.489)/2
sat_dists[sat_dists$distance_km == 258, ]$drainage_m2 <- (9414080727.7062 + 10051554840.489)/2
sat_dists[sat_dists$distance_km == 550, ]$drainage_m2 <- 22548106021.7803
sat_dists[sat_dists$distance_km == 552, ]$drainage_m2 <- 22548106021.7803
sat_dists[sat_dists$distance_km == 554, ]$drainage_m2 <- 22548106021.7803
sat_dists[sat_dists$distance_km == 556, ]$drainage_m2 <- 22548106021.7803
sat_dists[sat_dists$distance_km == 558, ]$drainage_m2 <- 22548106021.7803
sat_dists[sat_dists$distance_km == 560, ]$drainage_m2 <- 22548106021.7803
sat_dists[sat_dists$distance_km == 562, ]$drainage_m2 <- 22548106021.7803
ggplot(sat_dists, aes(x=distance_km, y=drainage_m2))+geom_line()

sentinel <- left_join(sentinel, sat_dists)
landsat <- left_join(landsat, sat_dists)

# Now extract drainage for in-situ stations
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
USGS_siteinfo <- readNWISsite(lapply(unique(USGS_Q$site_no), pad0))

station_dists <- data.frame(
  site_no = USGS_siteinfo$site_no,
  lat = as.numeric(USGS_siteinfo$dec_lat_va),
  lon = as.numeric(USGS_siteinfo$dec_long_va)
)
station_dists$drainage_m2 <- NA

for (i in 1:nrow(station_dists)) {
  startPoint <- c(station_dists[i, ]$lon, station_dists[i, ]$lat)
  riverData <- findNLDI(
    location = startPoint,
    nav = "UM",
    find = c("basin"),
    distance_km = 1000
    )
  basin_area <- st_area(riverData$basin) # m^2
  station_dists[i,]$drainage_m2 <- basin_area
}

station_dists$site_no <- as.numeric(station_dists$site_no)
USGS_Q <- left_join(USGS_Q, station_dists, by=("site_no"="site_no"))
# Fix site USGS-02335500
USGS_Q[USGS_Q$site_no == 2335500, ]$drainage_m2 <- 45606853279

naive_reg <- function(.x, dir){
  in_situ_flow <- USGS_Q %>% filter(sample_dt ==.x$date2[1])
  # Case 2: There are flow measurements, perform linear regression
  reg <- lm(log10(Q_m3s)~log10(drainage_m2), data = in_situ_flow)
  eq <- substitute(italic("Q (m^3/s) ") == 10^a %.%italic("Drainage (m^2)")^b, 
         list(a = format(unname(coef(reg)[1]), digits = 3),
              b = format(unname(coef(reg)[2]), digits = 3)))
    
  reg_plot = ggplot(in_situ_flow, aes(x=log10(drainage_m2), y=log10(Q_m3s))) +
    geom_point() +
    geom_smooth(method = "lm", se=FALSE, color="green", formula = y ~ x) +
    annotate(geom='text', label = as.character(as.expression(eq)), 
             parse = T, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    theme_bw() +
    labs(title=paste0('Date: ', .x$date2[1]),
        x=TeX('Drainage Area ($m^2$)'),
        y=TeX('Discharge ($m^3$/s)'))
  
  # Save regression plot
  ggsave(reg_plot, 
       filename = paste(dir, .x$date2[1], ".pdf"),
       width = 8, height = 8)
  .x[.x$in_situ_flow == FALSE, ]$Q_m3s <-10^ predict(reg, .x[.x$in_situ_flow == FALSE,])
  return(.x)
}

# Join landsat and flux
landsat_flux <- left_join(landsat, 
  USGS_Q %>% dplyr::select(Q_m3s, drainage_m2, sample_dt),
  by=c("drainage_m2"="drainage_m2", "sample_dt"="sample_dt"))
landsat_flux$in_situ_flow <- FALSE
landsat_flux[is.finite(landsat_flux$Q_m3s), ]$in_situ_flow <- TRUE

# Add dummy date to access date inside of group modify
landsat_flux$date2 <- landsat_flux$sample_dt

# Create directory and apply regression
wd_landsat_Q <- paste0(wd_exports, "landsat_Q_reg/")
landsat_flux <- landsat_flux %>%
  group_by(sample_dt) %>%
  group_modify(~naive_reg(., wd_landsat_Q)) %>% 
  subset(select=-c(date2))

# Join sentinel and flux
sentinel$sample_dt <- sentinel$date
sentinel_flux <- left_join(sentinel, 
  USGS_Q %>% dplyr::select(Q_m3s, drainage_m2, sample_dt),
  by=c("drainage_m2"="drainage_m2", "sample_dt"="sample_dt"))
sentinel_flux$in_situ_flow <- FALSE
sentinel_flux[is.finite(sentinel_flux$Q_m3s), ]$in_situ_flow <- TRUE

# Add dummy date to access date inside of group modify
sentinel_flux$date2 <- sentinel_flux$sample_dt

# Create directory and apply regression
wd_sentinel_Q <- paste0(wd_exports, "sentinel_Q_reg/")
sentinel_flux <- sentinel_flux %>%
  group_by(sample_dt) %>%
  group_modify(~naive_reg(., wd_sentinel_Q)) %>% 
  subset(select=-c(date2))

# Calculate fluxes
landsat_flux$flux_MTyr <- landsat_flux$Q_m3s * landsat_flux$SSC_mgL * 3.10585 * 10^-5
sentinel_flux$flux_MTyr <- sentinel_flux$Q_m3s * sentinel_flux$SSC_mgL * 3.10585 * 10^-5

# Save flux variables so they can be used later
write.csv(landsat_flux,
  file = paste0(wd_exports, "LANDSAT_FLUX.csv"), row.names = FALSE)
write.csv(sentinel_flux,
  file = paste0(wd_exports, "SENTINEL_FLUX.csv"), row.names = FALSE)

########################### TOTAL SSC/FLUX #####################################
# Investigate pre 1990 site
# One site is farther down but stopped operating around 1990 (29523008500460X)
mouthSites <- c(2359170, 295230085004601, 295230085004602, 295230085004603)
# maxDist <- 676

# 02359170 and 29523008500460X
m_INSITU_SSC <- USGS_SSC %>% filter(site_no %in% mouthSites)
m_INSITU_Q <- USGS_Q %>% filter(site_no %in% mouthSites)

# For satellite data, get the last 15km of river
m_SENTINEL <- sentinel_flux %>% filter(distance_km >= 695)
m_LANDSAT <- landsat_flux %>% filter(distance_km >= 695)

# Secondary figure --> reanalyze using USGS-02359170 cluster from regression
# clusters2Open <- "C:/Users/ilanv/Desktop/sentinel-ssc/exports/cluster_regress/"
# load(file = paste0(clusters2Open, "cluster_regressors.RData"))
# load(file = paste0(clusters2Open, "cluster_bands.RData"))

# sentinelBands <- setDT(copy(m_SENTINEL))[, ":="(
#   # Add squared columns
#   B1.2 = B1^2,
#   B2.2 = B2^2,
#   B3.2 = B3^2,
#   B4.2 = B4^2,
#   B5.2 = B5^2,
#   B6.2 = B6^2,
#   B7.2 = B7^2,
#   B8.2 = B8^2,
#   B8A.2 = B8A^2,
#   B9.2 = B9^2,
#   B11.2 = B11^2,
#   B12.2 = B12^2,
#   # Add square root columns
#   B1.0.5 = B1^0.5,
#   B2.0.5 = B2^0.5,
#   B3.0.5 = B3^0.5,
#   B4.0.5 = B4^0.5,
#   B5.0.5 = B5^0.5,
#   B6.0.5 = B6^0.5,
#   B7.0.5 = B7^0.5,
#   B8.0.5 = B8^0.5,
#   B8A.0.5 = B8A^0.5,
#   B9.0.5 = B9^0.5,
#   B11.0.5 = B11^0.5,
#   B12.0.5 = B12^0.5,
#   # Add band ratios
#   B2.B1=B2/B1,
#   B3.B1=B3/B1,
#   B4.B1=B4/B1,
#   B5.B1=B5/B1,
#   B6.B1=B6/B1,
#   B7.B1=B7/B1,
#   B8.B1=B8/B1,
#   B8A.B1=B8A/B1,
#   B9.B1=B9/B1,
#   B11.B1=B11/B1,
#   B12.B1=B12/B1,
  
#   B3.B2=B3/B2,
#   B4.B2=B4/B2,
#   B5.B2=B5/B2,
#   B6.B2=B6/B2,
#   B7.B2=B7/B2,
#   B8.B2=B8/B2,
#   B8A.B2=B8A/B2,
#   B9.B2=B9/B2,
#   B11.B2=B11/B2,
#   B12.B2=B12/B2,
  
#   B4.B3=B4/B3,
#   B5.B3=B5/B3,
#   B6.B3=B6/B3,
#   B7.B3=B7/B3,
#   B8.B3=B8/B3,
#   B8A.B3=B8A/B3,
#   B9.B3=B9/B3,
#   B11.B3=B11/B3,
#   B12.B3=B12/B3,
  
#   B5.B4=B5/B4,
#   B6.B4=B6/B4,
#   B7.B4=B7/B4,
#   B8.B4=B8/B4,
#   B8A.B4=B8A/B4,
#   B9.B4=B9/B4,
#   B11.B4=B11/B4,
#   B12.B4=B12/B4,
  
#   B6.B5=B6/B5,
#   B7.B5=B7/B5,
#   B8.B5=B8/B5,
#   B8A.B5=B8A/B5,
#   B9.B5=B9/B5,
#   B11.B5=B11/B5,
#   B12.B5=B12/B5,
  
#   B7.B6=B7/B6,
#   B8.B6=B8/B6,
#   B8A.B6=B8A/B6,
#   B9.B6=B9/B6,
#   B11.B6=B11/B6,
#   B12.B6=B12/B6,
  
#   B8.B7=B8/B7,
#   B8A.B7=B8A/B7,
#   B9.B7=B9/B7,
#   B11.B7=B11/B7,
#   B12.B7=B12/B7,
  
#   B8A.B8=B8A/B8,
#   B9.B8=B9/B8,
#   B11.B8=B11/B8,
#   B12.B8=B12/B8,
  
#   B9.B8A=B9/B8A,
#   B11.B8A=B11/B8A,
#   B12.B8A=B12/B8A,
  
#   B11.B9=B11/B9,
#   B12.B9=B12/B9,
  
#   B12.B11=B12/B11
# )]
# ssc_pred <- predict(
#     cluster_regressors[[2]],
#       newx = as.matrix(sentinelBands %>% dplyr::select(all_of(cluster_bands[[2]]))),
#       s = "lambda.min"
#     )

# m_SENTINEL$SSC2 <- 10^ssc_pred

# SSC at mouth
# ggplot(m_INSITU_SSC %>% filter(year >= 2018), aes(x = 1, y = SSC_mgL)) +
#   geom_boxplot(width=0.2, outlier.shape=NA, position= position_nudge(x=-.13),
#     aes(color="USGS")) +
#   geom_boxplot(data=m_LANDSAT %>% filter(year >= 2018), aes(color="Landsat"), width=0.2, outlier.shape=NA, 
#                position= position_nudge(x=+.13)) +
#   geom_boxplot(data=m_SENTINEL, aes(color="Sentinel (cluster-1)"), width=0.2, outlier.shape=NA, 
#                position= position_nudge(x=+.39)) +
#   geom_boxplot(data=m_SENTINEL, aes(color="Sentinel (cluster-2)"), width=0.2, outlier.shape=NA, 
#                position= position_nudge(x=+.65)) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   coord_cartesian(ylim=c(0, 100)) +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) +
#   # scale_color_manual(values = c('#9CDB79','#E69F00', '#999999')) +
#   theme(
#     legend.position = c(.97, .97),
#     legend.justification = c("right", "top"),
#     legend.box.background = element_rect(color = "black", size=0.55),
#     legend.box.margin = margin(6, 6, 6, 6)) +
#   labs(title = "SSC at the mouth of the Apalachicola",
#       #  caption="Reduced to one average sample per day.",
#        x = "",
#        y = "SSC (mg/L)",
#        color = "Data Source")

# ggsave(
#   filename = paste0(wd_figure, "/MOUTH_SSCADJUSTED.png"),
#   width = 6,
#   height = 6,
#   units = "in",
#   dpi = 200
# )

ggplot(m_INSITU_SSC, aes(x = year_chunk, y = SSC_mgL)) +
  geom_boxplot(width=0.2, outlier.shape=NA, position= position_nudge(x=-.13),
    aes(color="USGS")) +
  geom_boxplot(data=m_LANDSAT, aes(color="Landsat"), width=0.2, outlier.shape=NA, 
               position= position_nudge(x=+.13)) +
  geom_boxplot(data=m_SENTINEL, aes(color="Sentinel"), width=0.2, outlier.shape=NA, 
               position= position_nudge(x=+.39)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim=c(0, 150)) +
  scale_color_manual(values = c('#9CDB79','#E69F00', '#999999')) +
  theme(
    legend.position = c(.97, .97),
    legend.justification = c("right", "top"),
    legend.box.background = element_rect(color = "black", size=0.55),
    legend.box.margin = margin(6, 6, 6, 6)) +
  labs(title = "SSC at the mouth of the Apalachicola",
      #  caption="Reduced to one average sample per day.",
       x = "Year",
       y = "SSC (mg/L)",
       color = "Data Source")

ggsave(
  filename = paste0(wd_figure, "/MOUTH_SSCTIME.png"),
  width = 6,
  height = 6,
  units = "in",
  dpi = 200
)

# SSC at site 2338000 (lon, lat) -> (84.901194, 33.476528)
# Site that is more upstream
# R derived distance_downstream_km -> 134.9
# Switch to cluster 4 for this station

# sentinelBands <- setDT(copy(sentinel_flux))[, ":="(
#   # Add squared columns
#   B1.2 = B1^2,
#   B2.2 = B2^2,
#   B3.2 = B3^2,
#   B4.2 = B4^2,
#   B5.2 = B5^2,
#   B6.2 = B6^2,
#   B7.2 = B7^2,
#   B8.2 = B8^2,
#   B8A.2 = B8A^2,
#   B9.2 = B9^2,
#   B11.2 = B11^2,
#   B12.2 = B12^2,
#   # Add square root columns
#   B1.0.5 = B1^0.5,
#   B2.0.5 = B2^0.5,
#   B3.0.5 = B3^0.5,
#   B4.0.5 = B4^0.5,
#   B5.0.5 = B5^0.5,
#   B6.0.5 = B6^0.5,
#   B7.0.5 = B7^0.5,
#   B8.0.5 = B8^0.5,
#   B8A.0.5 = B8A^0.5,
#   B9.0.5 = B9^0.5,
#   B11.0.5 = B11^0.5,
#   B12.0.5 = B12^0.5,
#   # Add band ratios
#   B2.B1=B2/B1,
#   B3.B1=B3/B1,
#   B4.B1=B4/B1,
#   B5.B1=B5/B1,
#   B6.B1=B6/B1,
#   B7.B1=B7/B1,
#   B8.B1=B8/B1,
#   B8A.B1=B8A/B1,
#   B9.B1=B9/B1,
#   B11.B1=B11/B1,
#   B12.B1=B12/B1,
  
#   B3.B2=B3/B2,
#   B4.B2=B4/B2,
#   B5.B2=B5/B2,
#   B6.B2=B6/B2,
#   B7.B2=B7/B2,
#   B8.B2=B8/B2,
#   B8A.B2=B8A/B2,
#   B9.B2=B9/B2,
#   B11.B2=B11/B2,
#   B12.B2=B12/B2,
  
#   B4.B3=B4/B3,
#   B5.B3=B5/B3,
#   B6.B3=B6/B3,
#   B7.B3=B7/B3,
#   B8.B3=B8/B3,
#   B8A.B3=B8A/B3,
#   B9.B3=B9/B3,
#   B11.B3=B11/B3,
#   B12.B3=B12/B3,
  
#   B5.B4=B5/B4,
#   B6.B4=B6/B4,
#   B7.B4=B7/B4,
#   B8.B4=B8/B4,
#   B8A.B4=B8A/B4,
#   B9.B4=B9/B4,
#   B11.B4=B11/B4,
#   B12.B4=B12/B4,
  
#   B6.B5=B6/B5,
#   B7.B5=B7/B5,
#   B8.B5=B8/B5,
#   B8A.B5=B8A/B5,
#   B9.B5=B9/B5,
#   B11.B5=B11/B5,
#   B12.B5=B12/B5,
  
#   B7.B6=B7/B6,
#   B8.B6=B8/B6,
#   B8A.B6=B8A/B6,
#   B9.B6=B9/B6,
#   B11.B6=B11/B6,
#   B12.B6=B12/B6,
  
#   B8.B7=B8/B7,
#   B8A.B7=B8A/B7,
#   B9.B7=B9/B7,
#   B11.B7=B11/B7,
#   B12.B7=B12/B7,
  
#   B8A.B8=B8A/B8,
#   B9.B8=B9/B8,
#   B11.B8=B11/B8,
#   B12.B8=B12/B8,
  
#   B9.B8A=B9/B8A,
#   B11.B8A=B11/B8A,
#   B12.B8A=B12/B8A,
  
#   B11.B9=B11/B9,
#   B12.B9=B12/B9,
  
#   B12.B11=B12/B11
# )]

# ssc_pred <- predict(
#     cluster_regressors[[4]],
#       newx = as.matrix(sentinelBands %>% dplyr::select(all_of(cluster_bands[[4]]))),
#       s = "lambda.min"
#     )

# sentinel_flux$SSC4 <- 10^ssc_pred
# sentinel_flux$flux_MTyr4 <- sentinel_flux$SSC4 * sentinel_flux$Q_m3s * 3.10585 * 10^-5

# GEE derived distance_km -> 144 (use for satellite data)
ggplot(USGS_SSC %>% filter(site_no == 2338000 & year >= 1980), aes(x = year_chunk, y = SSC_mgL)) +
  geom_boxplot(width=0.2, outlier.shape=NA, position= position_nudge(x=-.13),
    aes(color="USGS")) +
  geom_boxplot(data=landsat_flux %>% filter(distance_km >= 140 & distance_km <= 148),
    aes(color="Landsat"), width=0.2, outlier.shape=NA, 
    position= position_nudge(x=+.13)) +
  geom_boxplot(data=sentinel_flux %>% filter(distance_km >= 140 & distance_km <= 148),
    aes(color="Sentinel"), width=0.2, outlier.shape=NA, 
    position= position_nudge(x=+.39)) +
  # theme_minimal() +
  theme_bw() +
  theme(
    legend.position = c(.97, .97),
    legend.justification = c("right", "top"),
    legend.box.background = element_rect(color = "black", size=0.55),
    legend.box.margin = margin(6, 6, 6, 6)) +
  coord_cartesian(ylim=c(0, 500)) +
  scale_color_manual(values = c('#9CDB79','#E69F00', '#999999')) +
  labs(title = "SSC at USGS-02338000, Chattahoochee River near Whitesburg, GA",
      #  caption="Reduced to one average sample per day.",
       x = "year",
       y = "SSC (mg/L)",
       color = "Data Source")

ggsave(
  filename = paste0(wd_figure, "/USGS-2338000.png"),
  width = 6,
  height = 6,
  units = "in",
  dpi = 200
)

### THIS SECTION NEEDS CHANGING, USE LANDSAT_FLUX AND SENTINEL_FLUX ###

# Now get flux
m_INSITU_FLUX <- left_join(m_INSITU_SSC, 
  m_INSITU_Q %>% dplyr::select(c("Q_m3s", "sample_dt")), 
  by=c("sample_dt"="sample_dt"))
m_INSITU_FLUX$flux_MTyr <- m_INSITU_FLUX$Q_m3s * m_INSITU_FLUX$SSC_mgL * 3.10585 * 10^-5

# # Assume discharge at site USGS-02359270 doesn't change till end of river
# m_SENTINEL_FLUX <- left_join(m_SENTINEL, m_INSITU_Q 
#   %>% filter(site_no == 2359170) 
#   %>% dplyr::select("Q_m3s", "sample_dt", "site_no"),
#   by=c("date"="sample_dt"))
# m_SENTINEL_FLUX$flux_MTyr <- m_SENTINEL_FLUX$Q_m3s * m_SENTINEL_FLUX$SSC_mgL * 3.10585 * 10^-5

# m_LANDSAT_FLUX <- left_join(m_LANDSAT, m_INSITU_Q
#   %>% filter(site_no == 2359170) 
#   %>% dplyr::select("Q_m3s", "sample_dt", "site_no"),
#   by = c("date" = "sample_dt"))
# m_LANDSAT_FLUX$flux_MTyr <- m_LANDSAT_FLUX$Q_m3s * m_LANDSAT_FLUX$SSC_mgL * 3.10585 * 10^-5

# NOW USING CLUSTER 2 FOR SENTINEL DATA
# m_SENTINEL$flux_MTyr <- m_SENTINEL$Q_m3s * m_SENTINEL$SSC2 * 3.10585 * 10^-5

# Plot flux at mouth of the Chattahoochee
ggplot(m_INSITU_FLUX, aes(x = year_chunk, y = flux_MTyr)) +
  geom_boxplot(width=0.2, outlier.shape=NA, position= position_nudge(x=-.13),
    aes(color="USGS")) +
  geom_boxplot(data=m_LANDSAT, aes(color="Landsat"), width=0.2, outlier.shape=NA, 
               position= position_nudge(x=+.13)) +
  geom_boxplot(data=m_SENTINEL, aes(color="Sentinel"), width=0.2, outlier.shape=NA, 
               position= position_nudge(x=+.39)) +
  # theme_minimal() +
  theme_bw() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim=c(0, 5)) +
  theme(
    legend.position = c(.97, .97),
    legend.justification = c("right", "top"),
    legend.box.background = element_rect(color = "black", size=0.55),
    legend.box.margin = margin(6, 6, 6, 6)) +
  scale_color_manual(values = c('#9CDB79','#E69F00', '#999999')) +
  labs(title = "Sediment flux at the mouth of the Apalachicola",
      #  caption="Reduced to one average sample per day.",
       x = "year",
       y = "Suspended Sediment Flux (Mt/yr)",
       color = "Data Source")

ggsave(
  filename = paste0(wd_figure, "/MOUTH_FLUX.png"),
  width = 6,
  height = 6,
  units = "in",
  dpi = 200
)

# Plot of number of samples at mouth of the Chattahoochee
# 2019 - 2021
sampdf <- data.frame(
  names = c("Sentinel", "Landsat", "USGS"),
  counts = c(
    nrow(m_SENTINEL %>% filter(year >= 2019 & year <= 2021)),
    nrow(m_LANDSAT %>% filter(year >= 2019 & year <= 2021)),
    nrow(m_INSITU_SSC %>% filter(year >= 2019 & year <= 2021)))
)

ggplot(sampdf, aes(x = names, y = counts, fill = names)) +
  geom_histogram(stat = "identity") +
  # theme_classic() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Number of Samples at the mouth of the Apalachicola River [2019-2021]",
    x = "",
    y = "Number of Samples"
  )

ggsave(
  filename = paste0(wd_figure, "/MOUTH_SAMP.png"),
  width = 4,
  height = 4,
  units = "in",
  dpi = 300
)

basin_area <- 49866869646  %>% set_units(m^2)#m^2
# Pre-anthropogenic estimates from octopusdata.org 
# Octupus ID: S156WTS010 Sample ID: SAP66
Be10E <- 9.15 %>% set_units(mm/kyr) %>% set_units(m/yr)
Be10dE <- 2.15 %>% set_units(mm/kyr) %>% set_units(m/yr)

# Pre-Dam From Trimble (1977), sourced originally from Dole and Stabler (1909)
# They used assumed bulk density of 1440 kg/m^3 
# Included 10% additional to account for Bed Load (B.L/)
Tflux <- ((set_units(0.057, mm/yr) %>% set_units(m/yr)) 
          * set_units(basin_area, m^2) 
          * set_units(1440, kg/m^3)) %>% set_units(Mt/yr)

# Get total flux over one year (using density of quartz, 2.648 g/cm^3) --> switched to 1440 from Trimble 1977
Be10flux <- (Be10E * set_units(basin_area, m^2) * set_units(1440, kg/m^3)) %>% set_units(Mt/yr)
Be10dflux <- (Be10dE * set_units(basin_area, m^2) * set_units(1440, kg/m^3)) %>% set_units(Mt/yr)
# Compare to Milliman and Reusser
# Increase by 10% to account for bed load
MillimanQ <- (set_units(22, km^3/yr)) %>% set_units(m^3/s) # Average discharge
MillimanTDS <- set_units(1.1 * 1.1, Mt/yr) # TDS load
MillimanDisYield <- set_units(0.021*1000, t/(km^2*yr)) # They use basin area of 52,000 km^2 but use 52 as divisor, need to multiply by 1,000
MillimanSSC <- set_units(50, mg/L) # Average SSC concentration

usgsavg <- mean(m_INSITU_FLUX$flux_MTyr, na.rm = TRUE) %>% set_units(Mt/yr)
landsatavg <- mean(m_LANDSAT$flux_MTyr, na.rm = TRUE) %>% set_units(Mt/yr)
sentinelavg <- mean(m_SENTINEL$flux_MTyr, na.rm = TRUE) %>% set_units(Mt/yr)
# Flux Dataframe --> increase satellite by 10% to account for bedload
Timeflux <- data.frame(
  source=c("Reusser (2015)", "Trimble (1977)", "USGS*", "Milliman (2013)", "Dethier (2020)**","Sentinel-2***"),
  flux=c(Be10flux, Tflux, usgsavg, MillimanTDS, landsatavg*1.1, sentinelavg*1.1),
  period=c("Pre-Colonial", "Pre-Dam", "Post-Dam", "Post-Dam", "Post-Dam","Post-Dam"),
  cl.boot = c(0, 0, smean.cl.boot(m_INSITU_FLUX$flux_MTyr)[["Upper"]]-smean.cl.boot(m_INSITU_FLUX$flux_MTyr)[["Mean"]], 0, smean.cl.boot(m_LANDSAT$flux_MTyr)[["Upper"]]-smean.cl.boot(m_LANDSAT$flux_MTyr)[["Mean"]], smean.cl.boot(m_SENTINEL$flux_MTyr)[["Upper"]]-smean.cl.boot(m_SENTINEL$flux_MTyr)[["Mean"]]),
  nsamp = c(1, 1, nrow(m_INSITU_FLUX), 1, nrow(m_LANDSAT), nrow(m_SENTINEL))
)
# Timeflux$se <- Timeflux$sd / sqrt(Timeflux$nsamp)
Timeflux$period <- factor(Timeflux$period, levels = c("Pre-Colonial", "Pre-Dam", "Post-Dam"))
Timeflux$source <- factor(Timeflux$source, levels = Timeflux$source)
Timeflux$flux <- drop_units(Timeflux$flux)

write.csv(Timeflux,
  file = paste0(wd_exports, "BASIN_FLUX.csv"), row.names = FALSE
)

ggplot(Timeflux, aes(x=source, y=flux, fill=period)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label=sprintf("%0.2f", round(flux, digits = 2))), vjust=-0.2) +
  stat_summary(
    data = m_INSITU_FLUX,
    aes(x = "USGS*", y = flux_MTyr, fill = "Post-Dam"),
    fun.data = mean_se, geom = "errorbar", width = 0.5) +
  stat_summary(
    data = m_LANDSAT,
    aes(x = "Dethier (2020)**", y = flux_MTyr * 1.1, fill = "Post-Dam"),
    fun.data = mean_se, geom = "errorbar", width = 0.5) +
  stat_summary(
    data = m_SENTINEL,
    aes(x = "Sentinel-2***", y = flux_MTyr * 1.1, fill = "Post-Dam"),
    fun.data = mean_se, geom = "errorbar", width = 0.5) +
  # theme_classic() +
  theme_bw() +
  # theme(legend.position = "bottom") +
   theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.background = element_rect(color = "black", size=0.75),
    legend.box.margin = margin(6, 6, 6, 6)) +
  scale_fill_manual(values=c('#999999','#E69F00','#9CDB79')) +
  labs(title="Sediment Flux out of the Apalachicola-Chattahoochee Basin",
       x = "",
       y = "Sediment Flux (Mt/yr)",
       fill = "Time Period",
       caption = "Using bulk density of 1440 kg/m^3 from Trimble (1977) and assuming bed load accounts for an extra 10% across SSC derived flux.\n*1981-present\n**1985-present\n***2018-present")

ggsave(
  filename = paste0(wd_figure, "FLUX_PERIODS.png"),
  width = 8,
  height = 8,
  units = "in",
  dpi = 150
)

################ INVESTIGATE SEASONAL VARIATION OF SSC/FLUX ####################

m_INSITU_FLUX$month <- format(as.Date(m_INSITU_FLUX$sample_dt), "%m")
m_LANDSAT$month <- format(as.Date(m_LANDSAT$sample_dt), "%m")
m_SENTINEL$month <- format(as.Date(m_SENTINEL$date), "%m")

# Now using the median
usgsMonFlux <- m_INSITU_FLUX %>% 
  filter(year >= 2000) %>% 
  group_by(month) %>%
  dplyr::summarize(meanFlux = median(flux_MTyr), samp = n(), meanQ = median(Q_m3s), meanSSC = median(SSC_mgL))
# Only present day data (2018 - present)
landsatMonFlux <- m_LANDSAT %>% 
  filter(year >= 2018) %>% 
  group_by(month) %>%
  dplyr::summarize(meanFlux = median(flux_MTyr, na.rm=T),samp = n(), meanQ = median(Q_m3s, na.rm=T), meanSSC = median(SSC_mgL,na.rm=T))
sentinelMonFlux <- m_SENTINEL %>% group_by(month) %>%
  dplyr::summarize(meanFlux = median(flux_MTyr, na.rm=T),samp = n(), meanQ = median(Q_m3s, na.rm=T), meanSSC = median(SSC_mgL, na.rm=T))

monDf <- data.frame(
  month =c(usgsMonFlux$month, landsatMonFlux$month, sentinelMonFlux$month), 
  flux_Mtyr = c(usgsMonFlux$meanFlux, landsatMonFlux$meanFlux, sentinelMonFlux$meanFlux),
  Q_m3s = c(usgsMonFlux$meanQ, landsatMonFlux$meanQ, sentinelMonFlux$meanQ),
  SSC_mgL = c(usgsMonFlux$meanSSC, landsatMonFlux$meanSSC, sentinelMonFlux$meanSSC),
  samples = c(usgsMonFlux$samp, landsatMonFlux$samp, sentinelMonFlux$samp),
  source = c(rep("USGS", 12), rep("Landsat", 12), rep("Sentinel", 12))
)

write.csv(monDf, 
  file = paste0(wd_exports, "MEDIAN_MONTHLY_BREAKDOWN.csv"), row.names = FALSE
)

monPlot <- melt(monDf, id = c("source", "month"))

ggplot(monPlot, aes(x = month, y = value, color = source, group = source, linetype = source)) +
  geom_line() +
  scale_y_continuous(limits=c(0, NA)) + 
  facet_wrap(variable ~ ., scales = "free_y", ncol=1, labeller = labeller(variable = 
    c("flux_Mtyr" = "Sediment Flux (Mt/yr)",
      "Q_m3s" = "Discharge (m^3/s)",
      "samples" = "Number of Samples",
      "SSC_mgL" = "Suspended Sediment Concentration (mg/L)")
     )) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(linetype = "none") +
  scale_color_manual(values=c('#9CDB79','#E69F00','#999999')) +
  labs(
    y = "", color = "Data Source", linetype = "",
    title = "Seasonal Changes in Median Flux Parameters at the Mouth of Apalachicola-Chattahoochee Basin",
    caption = "USGS Data from 2000 onwards due to sparse sampling"
  )

ggsave(
  filename = paste0(wd_figure, "SEASONAL_MOUTH.png"),
  width = 6,
  height = 10,
  units = "in",
  dpi = 150
)
# DONT USE PINWHEEL PLOTS 
# u <- ggplot(usgsMonFlux, aes(x=meanFlux/2, y=meanFlux, fill=month, width=meanFlux)) +
#   # scale_fill_manual(values = sentinel_colors) +
#  geom_bar(stat="identity", width = 1, color="white") +
#   coord_polar("y", start=0) +
#   scale_y_reverse() +
#   theme_void() +
#   theme(panel.border = element_blank(),panel.grid = element_blank()) +
#   labs(
#     caption = "July from Tropical Storm Alberto (7-1994):\nhttps://gwri.gatech.edu/sites/default/files/files/docs/1995/StameyT-95.pdf"
#   )

# s <- ggplot(sentinelMonFlux, aes(x=meanFlux/2, y=meanFlux, fill=month, width=meanFlux)) +
#   # scale_fill_manual(values = sentinel_colors) +
#  geom_bar(stat="identity", width = 1, color="white") +
#   coord_polar("y", start=0) +
#   scale_y_reverse() +
#   theme_void() +
#   theme(panel.border = element_blank(),panel.grid = element_blank())

# l <- ggplot(landsatMonFlux, aes(x=meanFlux/2, y=meanFlux, fill=month, width=meanFlux)) +
#   # scale_fill_manual(values = sentinel_colors) +
#  geom_bar(stat="identity", width = 1, color="white") +
#   coord_polar("y", start=0) +
#   scale_y_reverse() +
#   theme_void() +
#   theme(panel.border = element_blank(),panel.grid = element_blank()) +
#   labs(
#     caption = "Extremely wet spring-fall 2009 (4-2009):\nhttps://gwri.gatech.edu/sites/default/files/files/docs/2011/2.6.1_Knaak_15.pdf"
#   )

# p <- plot_grid(u, l, s, labels=c("USGS", "Landsat", "Sentinel"))

# save_plot(
#   p,
#   filename = paste0(wd_figure, "/MOUTH_SEASONAL_FLUX.png"),
#   base_width = 30,
#   base_height = 30,
#   units = "cm"
# )

# Another visualization
# ggplot() +
#   # Landsat data
#   stat_summary(
#     data = m_LANDSAT,
#     aes(x = month, y = flux_MTyr, color='Landsat', group=1),
#     geom = 'line', 
#     fun = 'median',
#     inherit.aes = FALSE,
#     na.rm=T) +
#   # Sentinel data
#   stat_summary(
#     data = m_SENTINEL,
#     aes(x = month, y = flux_MTyr, color='Sentinel', group=1),
#     geom = 'line', 
#     fun = 'median',
#     inherit.aes = FALSE,
#     na.rm=T) +
#   # USGS data
#   stat_summary(
#     data = m_INSITU_FLUX,
#     aes(x = month, y = flux_MTyr, color='USGS', group=1),
#     fun = 'median',
#     geom = "line",
#     inherit.aes = FALSE,
#     na.rm=T) +
#   theme_classic() +
#   scale_x_discrete(labels=c(
#     "01"="Jan",
#     "02"="Feb",
#     "03"="Mar",
#     "04"="Apr",
#     "05"="May", 
#     "06"="Jun",
#     "07"="Jul",
#     "08"="Aug",
#     "09"="Sep",
#     "10"="Oct",
#     "11"="Nov", 
#     "12"="Dec")) +
#   labs(
#     title = "Median Sediment Flux through the Chattahoochee River Watershed",
#     x = "Month",
#     y = "Sediment Flux (Mt/yr)",
#     color = "Data Source"
#   )

# ggsave(
#   filename = paste0(wd_figure, "MOUTH_SEASONAL_FLUXLINE.png"),
#   width = 20,
#   height = 14,
#   units = "cm"
# )

# ggplot() +
#   # Landsat data
#   stat_summary(
#     data = m_LANDSAT_FLUX,
#     aes(x = month, y = SSC_mgL, color='Landsat', group=1),
#     geom = 'line', 
#     fun = 'median',
#     inherit.aes = FALSE,
#     na.rm=T) +
#   # Sentinel data
#   stat_summary(
#     data = m_SENTINEL_FLUX,
#     aes(x = month, y = SSC_mgL, color='Sentinel', group=1),
#     geom = 'line', 
#     fun = 'median',
#     inherit.aes = FALSE,
#     na.rm=T) +
#   # USGS data
#   stat_summary(
#     data = m_INSITU_FLUX,
#     aes(x = month, y = SSC_mgL, color='USGS', group=1),
#     fun = 'median',
#     geom = "line",
#     inherit.aes = FALSE,
#     na.rm=T) +
#   theme_classic() +
#   scale_x_discrete(labels=c(
#     "01"="Jan",
#     "02"="Feb",
#     "03"="Mar",
#     "04"="Apr",
#     "05"="May", 
#     "06"="Jun",
#     "07"="Jul",
#     "08"="Aug",
#     "09"="Sep",
#     "10"="Oct",
#     "11"="Nov", 
#     "12"="Dec")) +
#   labs(
#     title = "Median SSC at the mouth of Chattahoochee River Watershed",
#     x = "Month",
#     y = "SSC (mg/L)",
#     color = "Data Source"
#   )

# ggsave(
#   filename = paste0(wd_figure, "MOUTH_SEASONAL_SSC.png"),
#   width = 20,
#   height = 14,
#   units = "cm"
# )

# # Got annual flow from https://waterdata.usgs.gov/nwis/annual?referred_module=sw&amp;site_no=02359170&amp;por_02359170_26942=2396790,00060,26942,1977,2022&amp;start_dt=1984&amp;end_dt=2022&amp;year_type=W&amp;format=html_table&amp;date_format=YYYY-MM-DD&amp;rdb_compression=file&amp;submitted_form=parameter_selection_list
# yrs <- seq(1984,2021)
# annual_q <- data.table(year=yrs[!yrs==2017], flow_ft3s=c(
#   31460,
#   16070,
#   19050,
#   25600,
#   16380,
#   19270,
#   29040,
#   26200,
#   21590,
#   28180,
#   30900,
#   26820,
#   26300,
#   24250,
#   36860,
#   18020,
#   10990,
#   16890,
#   10560,
#   28960,
#   18170,
#   35230,
#   16330,
#   12420,
#   14010,
#   23590,
#   32160,
#   12100,
#   9715,
#   22340,
#   25360,
#   18680,
#   28020,
#   19050,
#   30210,
#   32090,
#   29210
# ))
# annual_q$flow_m3s <- annual_q$flow_ft3s * 0.02832
# annual_q$year_chunk <- cut(annual_q$year, seq(from = 1895, to = 2026, by = 5), right=FALSE)
# annual_q <- annual_q %>% group_by(year_chunk) %>% summarize(flow_m3s = mean(flow_m3s))

# ggplot(annual_q, aes(x=year_chunk, y=flow_m3s)) +
#   geom_line(group=1) +
#   theme_classic() +
#   theme(legend.position = "bottom") +
#   labs(
#     x = "",
#     y = "Average Q (m^3/s)"
#   )


######################### HISTOGRAM OF SSC SAMPLES #############################
# Only 13 samples are > 99.98 percentile -> they are omitted initially
# c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
# c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
c1 <- "#8dc7a0"
c2 <- "#999999"
# s1 <- hist(sentinel$SSC_mgL, breaks = 100,plot = FALSE)
# s2 <- hist(landsat[landsat$SSC_mgL < 500, ]$SSC_mgL, breaks = 100, plot = FALSE)
# plot(s2, col=c2)
# plot(s1, col=c1, add = TRUE,
#   main = "Histogram of ")
ggplot() +
  geom_histogram(
    data = landsat %>% filter(SSC_mgL < 300),
    aes(x=SSC_mgL, fill="Landsat"),
    binwidth = 5
  ) +
  geom_histogram(
    data = sentinel,
    aes(x=SSC_mgL, fill="Sentinel"),
    binwidth = 5
  ) +
  theme_classic() + 
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.background = element_rect(color = "black", size=0.75),
    legend.box.margin = margin(6, 6, 6, 6)) +
  labs(
    title = "Comparison of Sentinel-2 and Landsat 5/7 Samples",
    x = "SSC (mg/L)",
    y = "Count",
    fill = "Data Source"
  )

ggsave(
  filename = paste0(wd_figure, "SAMPLE_HIST.png"),
  width = 20,
  height = 14,
  units = "cm"
)

####################### REFLECTANCE/Q PINWHEEL #################################
# Get sentinel working then can repeat with landsat, etc
monCol <- data.frame(
  month = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
  mean_Q = NA
)

mouthSites <- c(2359170, 295230085004601, 295230085004602, 295230085004603)
# 02359170 and 29523008500460X
m_INSITU_SSC <- USGS_SSC %>% filter(site_no %in% mouthSites)
m_INSITU_Q <- USGS_Q %>% filter(site_no %in% mouthSites)
m_INSITU_Q$month <- format(as.Date(m_INSITU_Q$sample_dt), "%m")
USGS_Q_bymonth <- m_INSITU_Q %>% group_by(month) %>% summarize(meanQ_m3s = mean(Q_m3s, na.rm=T))
monCol$mean_Q <- USGS_Q_bymonth$meanQ_m3s

# Get sentinel pinwheel
sentinel_raw <- read.csv(file="C:\\Users\\ilanv\\Desktop\\sentinel-ssc\\exports\\transects\\SENTINEL2_20m_SSC.csv")
sentinel_raw$month <- format(as.Date(sentinel_raw$date), "%m")
# Get reflectance at mouth of river (last 15 km)
sentinel_bymonth <- sentinel_raw %>% filter(distance_km >= 695) %>%
  group_by(month) %>%
  summarize(
    r = mean(B4),
    g = mean(B3),
    b = mean(B2)
  )

sentinel_colors <- rgb(sentinel_bymonth$r/2200, sentinel_bymonth$g/2200, sentinel_bymonth$b/2200)

ggplot(monCol, aes(x=mean_Q/2, y=mean_Q, fill=month, width=mean_Q)) +
  scale_fill_manual(values = sentinel_colors) +
 geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) +
  scale_y_reverse() +
  theme_void() +
  theme(legend.position = "none",panel.border = element_blank(),panel.grid = element_blank())

ggsave(
  filename = paste0(wd_figure, "SENTINEL_PINWHEEL.png"),
  width = 10,
  height = 10,
  units = "cm"
)

# Get landsat pinwheel
landsat_raw <- read.csv(file="dethier_ls57_chattahoochee.csv")
landsat_raw$month <- format(as.Date(landsat_raw$sample_dt), "%m")
# Get reflectance at mouth of river (last 15 km)
landsat_bymonth <- landsat_raw %>% filter(distance_km >= 695) %>%
  group_by(month) %>%
  summarize(
    r = mean(B3),
    g = mean(B2),
    b = mean(B1)
  )

# /3000 for true color
landsat_colors <- rgb(landsat_bymonth$r/3000, landsat_bymonth$g/3000, landsat_bymonth$b/3000)

ggplot(monCol, aes(x=mean_Q/2, y=mean_Q, fill=month, width=mean_Q, group=interaction(mean_Q, month))) +
  scale_fill_manual(values = landsat_colors) +
  geom_bar(stat="identity", width = 1, color="white") +
  coord_polar("y", start=0) +
  scale_y_reverse() +
  theme_void() +
  theme(legend.position = "none",panel.border = element_blank(),panel.grid = element_blank())

ggsave(
  filename = paste0(wd_figure, "LANDSAT_PINWHEEL.png"),
  width = 10,
  height = 10,
  units = "cm"
)

# Now do SSC Samples over span of 1 year
# sentinel_monthssc <- sentinel_raw %>% filter(distance_km >= 695) %>%
#   group_by(month) %>%
#   summarize(
#     monSSC = mean(SSC_mgL)
#   )
# landsat_monthssc <- landsat_raw %>% filter(distance_km >= 695 & year > 2017) %>%
#   group_by(month) %>%
#   summarize(
#     monSSC = mean(SSC_mgL)
#   )

######################## TRAPPING EFFICIENCIES #################################

# - Trapping efficiency of 1 -> no sediment gets past
# - For dam at a distance x:
# - Get avg flux at start of resevoir-10 km (pre-dam) to x+10km (post-dam)
# - TE = -(pre-post)/(pre)
# - Now using SSC not flux

# Function to compute 0.95 confidence interval for the difference in two means
# https://thomaselove.github.io/2018-431-book/CI-Indep-Samples.html
# g is grouping variable --> from https://www.rdocumentation.org/packages/Hmisc/versions/5.0-0/topics/smean.sd
bootdif <- function(y, g) {
 g <- as.factor(g)
 a <- attr(smean.cl.boot(y[g==levels(g)[1]], B=500, reps=TRUE),'reps')
 b <- attr(smean.cl.boot(y[g==levels(g)[2]], B=500, reps=TRUE),'reps')
 meandif <- diff(tapply(y, g, mean, na.rm=TRUE))
 a.b <- quantile(b-a, c(.025,.975))
 res <- c(meandif, a.b)
 names(res) <- c('Mean Difference','.025','.975')
 res
}

sentinel_dam <- copy(sentinel_flux)
sentinel_dam$dam <- "Dam"
sentinel_dam$damName <- NA
landsat_dam <- copy(landsat_flux)
landsat_dam$dam <- "Dam"
landsat_dam$damName <- NA
for (i in 1:nrow(DAMS)) {
  dam <- DAMS[i, ]
  name <- dam$Dam.Name
  # Don't look at dist --> where dam is
  damdist <- dam$distance_km
  resdist <- dam$res_start
  if (damdist == 0) { # Buford Dam case --> ignore
    # Only alter dam data
    # sentinel_dam <- sentinel_dam %>% mutate(distance_km =
    #   ifelse(distance_km <= damdist + 10 & distance_km > damdist,
    #     damdist,
    #     distance_km
    #   ),
    #   damName = name
    # )
  } else {
    # Alter reservoir data
    sentinel_dam <- sentinel_dam %>% mutate(
      distance_km = ifelse(distance_km >= resdist - 10 & distance_km <= resdist,
        resdist,
        distance_km),
      damName = ifelse(distance_km == resdist,
        name,
        damName)
    )
    # Alter dam data
    sentinel_dam <- sentinel_dam %>% mutate(
      distance_km = ifelse(distance_km <= damdist + 10 & distance_km > damdist,
        damdist,
        distance_km),
      damName = ifelse(distance_km == damdist,
        name,
        damName)
    )
  }
}

# Now remove data not at the dam or reservoir locations
# sentinel_dam <- sentinel_dam %>% filter(distance_km %in% c(DAMS$distance_km, DAMS$res_start))
sentinel_dam <- sentinel_dam %>% filter(!is.na(damName))
sentinel_dam <- sentinel_dam %>% mutate(
  dam = ifelse(distance_km %in% DAMS$res_start,
    "reservoir",
    "dam"
  )
)

for (i in 1:nrow(DAMS)) {
  dam <- DAMS[i, ]
  name <- dam$Dam.Name
  # Don't look at dist --> where dam is
  damdist <- dam$distance_km
  resdist <- dam$res_start
  if (damdist == 0) { # Buford Dam case --> ignore
    # Only alter dam data
    # sentinel_dam <- sentinel_dam %>% mutate(distance_km =
    #   ifelse(distance_km <= damdist + 10 & distance_km > damdist,
    #     damdist,
    #     distance_km
    #   ),
    #   damName = name
    # )
  } else {
    # Alter reservoir data
    landsat_dam <- landsat_dam %>% mutate(
      distance_km = ifelse(distance_km >= resdist - 10 & distance_km <= resdist,
        resdist,
        distance_km),
      damName = ifelse(distance_km == resdist,
        name,
        damName)
    )
    # Alter dam data
    landsat_dam <- landsat_dam %>% mutate(
      distance_km = ifelse(distance_km <= damdist + 10 & distance_km > damdist,
        damdist,
        distance_km),
      damName = ifelse(distance_km == damdist,
        name,
        damName)
    )
  }
}

# Now remove data not at the dam or reservoir locations
# sentinel_dam <- sentinel_dam %>% filter(distance_km %in% c(DAMS$distance_km, DAMS$res_start))
landsat_dam <- landsat_dam %>% filter(!is.na(damName))
landsat_dam <- landsat_dam %>% mutate(
  dam = ifelse(distance_km %in% DAMS$res_start,
    "reservoir",
    "dam"
  )
)

# Definition of whiskers
# f <- function(x) {
#   r <- quantile(x, probs = c(0.10, 0.25, 0.5, 0.75, 0.90))
#   names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
#   r
# }

# ggplot() +
#   stat_summary(
#     data = sentinel_flux,
#     aes(x = distance_km, y = SSC_mgL, color='Sentinel'),
#     geom = 'line', 
#     fun = 'mean',
#     inherit.aes = FALSE,
#     na.rm=T) +
#     stat_summary(data = sentinel_dam,
#       aes(x = distance_km, y = SSC_mgL, color = dam),
#       fun.data = "mean_cl_boot",
#       geom = "crossbar",
#       width = 20) +
#     # geom_vline(xintercept = 26) +
#     # stat_summary(data = sentinel_dam,
#     #   aes(x = distance_km, y = SSC_mgL, group=distance_km, color=dam),
#     #   fun.data = f,
#     #   geom="boxplot",
#     #   width = 20
#     # ) +
#     coord_cartesian(ylim=c(0, 40)) +
#     theme_bw() + 
#     labs(
#       title = "testing",
#       x = "Distance Downstream of Buford Dam",
#       y = "SSC (mg/L)",
#       color = ""
#     )

# ggsave(
#   filename = paste0(wd_figure, "/SENTINEL_TE_BOXES.png"),
#   width = 10,
#   height = 8,
#   units = "in",
#   dpi = 150
# )

# Trapping efficiencies for everything but Buford and Woodruff (just use SSC)

sentinel_damstat <- copy(sentinel_dam)
sentinel_damstat <- sentinel_damstat %>% group_by(damName, dam) %>%
  summarise(
    avgSSC_mgL = smean.cl.boot(SSC_mgL)[["Mean"]], 
    boot95 = abs(smean.cl.boot(SSC_mgL)[["Mean"]] - smean.cl.boot(SSC_mgL)[["Lower"]]), # bootstrap error 
    n = n()
  )
View(sentinel_damstat)

landsat_damstat <- copy(landsat_dam)
landsat_damstat <- landsat_damstat %>% group_by(damName, dam) %>%
  summarise(
    avgSSC_mgL = smean.cl.boot(SSC_mgL)[["Mean"]], 
    boot95 = abs(smean.cl.boot(SSC_mgL)[["Mean"]] - smean.cl.boot(SSC_mgL)[["Lower"]]), # bootstrap error 
    n = n()
  )

View(landsat_damstat)

landsat_damstat2018 <- copy(landsat_dam) %>% filter(year >= 2018)
landsat_damstat2018 <- landsat_damstat2018 %>% group_by(damName, dam) %>%
  summarise(
    avgSSC_mgL = smean.cl.boot(SSC_mgL)[["Mean"]], 
    boot95 = abs(smean.cl.boot(SSC_mgL)[["Mean"]] - smean.cl.boot(SSC_mgL)[["Lower"]]), # bootstrap error 
    n = n()
  )
  
View(landsat_damstat2018)

# sentinel_te <- sentinel_damstat %>% group_by(damName) %>%
#   do(mutate(., 
#     TE = -(.[dam == "dam",]$avgSSC_mgL - .[dam =="reservoir",]$avgSSC_mgL)/.[dam =="reservoir",]$avgSSC_mgL
#   ))

# View(sentinel_te)

# sentinel_te <- copy(sentinel_dam)
# sentinel_te <- sentinel_te %>% group_by(damName) %>%
#   summarise(
#     dssc = bootdif(SSC_mgL, dam)[["Mean Difference"]]
#   )

# trapping_eff <- function(dist, res_start, years,  df){
#   last_dist <- res_start - 10
#   next_dist <- dist + 10
#   # Don't look at dist --> where dam is
#   pre <<- df %>% filter(distance_km >= last_dist & distance_km <= res_start) %>% filter(year_chunk == years)
#   post <<- df %>% filter(distance_km <= next_dist & distance_km > dist) %>% filter(year_chunk == years)
#   pre_ssc <- mean(pre$SSC_mgL, na.rm = TRUE)
#   sdPre <- sd(pre$SSC_mgL, na.rm = TRUE)
#   post_ssc <- mean(post$SSC_mgL, na.rm = TRUE)
#   sdPost <- sd(post$SSC_mgL, na.rm = TRUE)
#   dssc <- post_ssc-pre_ssc
#   dsscSD <- sqrt((sdPre)^2 + (sdPost)^2)
#   te <- -dssc / pre_ssc
#   l <- list("pre_ssc"=pre_ssc, "preSD" = sdPre, "post_ssc"=post_ssc, "postSD" = sdPost, "te"=te, "dssc"=dssc, "dsscSD" = dsscSD)
#   return(l)
# }

# teDAMS <- DAMS %>% filter(distance_km != 0 & distance_km != 552)
# teDAMS <- data.frame(
#   Dam.Name = teDAMS$Dam.Name, 
#   distance_km = teDAMS$distance_km,
#   res_start = teDAMS$res_start,
#   Volume_m3 = teDAMS$V_m3,
#   years = rep(unique(landsat_flux$year_chunk), each=3))

# te_landsat <- teDAMS %>% rowwise() %>% 
#   mutate(te = trapping_eff(distance_km, res_start, years, landsat_flux)$te) %>%
#   mutate(pre_ssc_mgL = trapping_eff(distance_km, res_start, years, landsat_flux)$pre_ssc) %>%
#   mutate(preSD = trapping_eff(distance_km, res_start, years, landsat_flux)$preSD) %>%
#   mutate(post_ssc_mgL = trapping_eff(distance_km, res_start, years, landsat_flux)$post_ssc) %>%
#   mutate(postSD = trapping_eff(distance_km, res_start, years, landsat_flux)$postSD) %>%
#   mutate(dssc_mgL = trapping_eff(distance_km, res_start, years, landsat_flux)$dssc) %>%
#   mutate(dsscSD = trapping_eff(distance_km, res_start, years, landsat_flux)$dsscSD)
#   # mutate(time2fill_yr = (1.440*10**(-6) * Volume_m3) / -dflux_Mt)

# te_sentinel <- teDAMS %>% rowwise() %>% 
#   mutate(te = trapping_eff(distance_km, res_start, years, sentinel_flux)$te) %>%
#   mutate(pre_ssc_mgL = trapping_eff(distance_km, res_start, years, sentinel_flux)$pre_ssc) %>%
#   mutate(preSD = trapping_eff(distance_km, res_start, years, sentinel_flux)$preSD) %>%
#   mutate(post_ssc_mgL = trapping_eff(distance_km, res_start, years, sentinel_flux)$post_ssc) %>%
#   mutate(postSD = trapping_eff(distance_km, res_start, years, sentinel_flux)$postSD) %>%
#   mutate(dssc_mgL = trapping_eff(distance_km, res_start, years, sentinel_flux)$dssc) %>%
#   mutate(dsscSD = trapping_eff(distance_km, res_start, years, sentinel_flux)$dsscSD)

# write.csv(te_landsat,
#   file = paste0(wd_exports, "TE_LANDSAT_BYYEAR.csv"), row.names = FALSE
# )
# write.csv(te_sentinel,
#   file = paste0(wd_exports, "TE_SENTINEL_BYYEAR.csv"), row.names = FALSE
# )

# Plot Trapping efficiencies over time
# ggplot() +
#   geom_line(data = te_landsat, 
#     aes(x=years,  y = te, color=Dam.Name, group=Dam.Name)) +
#   geom_line(data = te_sentinel, 
#     aes(x=years,  y = te, color=Dam.Name, group=Dam.Name), na.rm=T)

# Lake seminole specific (dam 552km, resevoir 520km)
chattahoochee_res <- sentinel_dam %>% filter(distance_km == 520)
chattahoochee_dam <- sentinel_dam %>% filter(distance_km == 552)
flint_res <- read.csv("C:/Users/ilanv/Desktop/sentinel-ssc/exports/transects/FLINT_20m_SSC.csv")
flint_res <- flint_res %>% group_by(date) %>% summarize(SSC_mgL = mean(SSC_mgL))

# Get drainage relationship
flintStart <- c(-84.63772, 30.80044)
riverData <- findNLDI(
    location = flintStart,
    nav = "UM",
    find = c("basin"),
    distance_km = 1000
    )
flint_basin_area <- st_area(riverData$basin) %>% drop_units() # m^2 
flint_res$drainage_m2 <- flint_basin_area 

# Redo Q vs drainage to get flow on the Flint
# Now extract drainage for in-situ stations
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
USGS_siteinfo <- readNWISsite(lapply(unique(USGS_Q$site_no), pad0))

station_dists <- data.frame(
  site_no = USGS_siteinfo$site_no,
  lat = as.numeric(USGS_siteinfo$dec_lat_va),
  lon = as.numeric(USGS_siteinfo$dec_long_va)
)
station_dists$drainage_m2 <- NA

for (i in 1:nrow(station_dists)) {
  startPoint <- c(station_dists[i, ]$lon, station_dists[i, ]$lat)
  riverData <- findNLDI(
    location = startPoint,
    nav = "UM",
    find = c("basin"),
    distance_km = 1000
    )
  basin_area <- st_area(riverData$basin) # m^2
  station_dists[i,]$drainage_m2 <- basin_area
}

station_dists$site_no <- as.numeric(station_dists$site_no)
USGS_Q <- left_join(USGS_Q, station_dists, by=("site_no"="site_no"))
# Fix site USGS-02335500
USGS_Q[USGS_Q$site_no == 2335500, ]$drainage_m2 <- 45606853279

flint_res$Q_m3s <- NA
for (d in unique(flint_res$date)) {
  in_situ_flow <- USGS_Q %>% filter(sample_dt == d)
  if (nrow(in_situ_flow) > 0) {
    reg <- lm(log10(Q_m3s)~log10(drainage_m2), data = in_situ_flow)
    flint_res[flint_res$date == d, ]$Q_m3s <-10^predict(reg, flint_res[flint_res$date == d,])
  }
}

flint_res$flux_MTyr <- flint_res$SSC_mgL * flint_res$Q_m3s * 3.10585 * 10^-5
flint_res <- flint_res %>% filter(is.finite(flux_MTyr))
# Calculate basin fraction compared to Chattahoochee 
# chattahoochee_basin_area <- sentinel_dam[sentinel_dam$distance_km == 520, ]$drainage_m2[1]
# flintScale <-  flint_basin_area / chattahoochee_basin_area

# Combine all data above the resevoir (lake seminole)
seminole_res <- data.frame(
  date = unique(c(chattahoochee_res$date, flint_res$date)),
  SSC_flux = NA
)

# Combine SSC on days when there is a sample on each river
for (i in 1:nrow(seminole_res)) {
  match_dt <- seminole_res[i, ]$date
  seminole_res[i, ]$SSC_flux <- mean(chattahoochee_res[chattahoochee_res$date == match_dt, ]$flux_MTyr) + mean(flint_res[flint_res$date == match_dt, ]$flux_MTyr)
}

# Remove days without both samples (will be NaN)
seminole_res <- seminole_res %>% filter(is.finite(SSC_flux))

# For Jim Woodruff -> "avgSSC_mgL" is actually sediment flux
sentinel_damstat[sentinel_damstat$damName == "Jim Woodruff Lock and Dam" & sentinel_damstat$dam == "reservoir", ] <- sentinel_damstat %>%
    filter(damName == "Jim Woodruff Lock and Dam" & dam == "reservoir") %>%
    mutate(
      avgSSC_mgL = smean.cl.boot(seminole_res$SSC_flux)[["Mean"]],
      boot95 = abs(smean.cl.boot(seminole_res$SSC_flux)[["Mean"]] - smean.cl.boot(seminole_res$SSC_flux)[["Lower"]]),
      n = length(seminole_res$SSC_flux)
      )

sentinel_damstat[sentinel_damstat$damName == "Jim Woodruff Lock and Dam" & sentinel_damstat$dam == "dam", ] <- sentinel_damstat %>%
  filter(damName == "Jim Woodruff Lock and Dam" & dam == "dam") %>%
  mutate(
    avgSSC_mgL = smean.cl.boot(chattahoochee_dam$flux_MTyr)[["Mean"]],
    boot95 = abs(smean.cl.boot(chattahoochee_dam$flux_MTyr)[["Mean"]] - smean.cl.boot(chattahoochee_dam$flux_MTyr)[["Lower"]]),
    n = length(chattahoochee_dam$flux_MTyr)
    )

sentinel_te <- sentinel_damstat %>% group_by(damName) %>%
  do(mutate(., 
    TE = -(.[dam == "dam",]$avgSSC_mgL - .[dam =="reservoir",]$avgSSC_mgL)/.[dam =="reservoir",]$avgSSC_mgL
  ))

View(sentinel_te)

landsat_te <- landsat_damstat %>% group_by(damName) %>%
  do(mutate(., 
    TE = -(.[dam == "dam",]$avgSSC_mgL - .[dam =="reservoir",]$avgSSC_mgL)/.[dam =="reservoir",]$avgSSC_mgL
  ))

View(landsat_te)

landsat_te2018 <- landsat_damstat2018 %>% group_by(damName) %>%
  do(mutate(., 
    TE = -(.[dam == "dam",]$avgSSC_mgL - .[dam =="reservoir",]$avgSSC_mgL)/.[dam =="reservoir",]$avgSSC_mgL
  ))

View(landsat_te2018)

sentinel_te <- left_join(sentinel_te, DAMS[c("Dam.Name", "V_m3", "distance_km")], by=c("damName"="Dam.Name"))
landsat_te <- left_join(landsat_te, DAMS[c("Dam.Name", "V_m3", "distance_km")], by=c("damName"="Dam.Name"))
landsat_te2018 <- left_join(landsat_te2018, DAMS[c("Dam.Name", "V_m3", "distance_km")], by=c("damName"="Dam.Name"))

write.csv(sentinel_te,
  file = paste0(wd_exports, "BOOTSTRAP_TE_SENTINEL.csv"), row.names = FALSE
)
write.csv(landsat_te,
  file = paste0(wd_exports, "BOOTSTRAP_TE_LANDSAT.csv"), row.names = FALSE
)
write.csv(landsat_te2018,
  file = paste0(wd_exports, "BOOTSTRAP_TE_LANDSAT2018.csv"), row.names = FALSE
)


# 2018-now
ggplot() +
  # Landsat Error
  stat_summary(
    data = landsat_flux %>% filter(year >= 2018),
    aes(x = distance_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#d52442", alpha = 0.3) +
  # USGS data
  stat_summary(
    data = USGS_SSC %>% filter(year >= 2018),
    aes(x = dist_downstream_km, y = SSC_mgL, color='USGS'),
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Landsat Error
  stat_summary(
    data = USGS_SSC %>% filter(year >= 2018),
    aes(x = dist_downstream_km, y = SSC_mgL, group = 1), 
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "errorbar",
    width = 20) +
  # Sentinel data
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = SSC_mgL, color='Sentinel'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Sentinel Error bar
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ), 
    geom = "ribbon", 
    fill = "#49b07b", alpha = 0.3) +
  # Landsat data
  stat_summary(
    data = landsat_flux %>% filter(year >= 2018),
    aes(x = distance_km, y = SSC_mgL, color='Landsat'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Sentinel TE Boxplot
  stat_summary(data = sentinel_dam %>% filter(damName == "West Point Dam"),
      aes(x = distance_km, y = SSC_mgL),
      fun.data = "mean_cl_boot",
      geom = "crossbar",
      width = 20) +
  # Landsat TE Boxplot 
  # stat_summary(data = landsat_dam %>% filter(year >= 2018),
  #     aes(x = distance_km, y = SSC_mgL, color = damName),
  #     fun.data = "mean_cl_boot",
  #     geom = "crossbar",
  #     width = 20) +
  coord_cartesian(ylim = c(0, 200)) +
  # Add Dam locations
  geom_vline(xintercept = DAMS$distance_km, linetype = "dotted") +
  theme_bw() +
  labs(
    title = "Suspended Sediment Concentration from 2018 - Present",
    x = "Distance Downstream From Buford Dam (km)",
    y = "SSC (mg/L)",
    color = "Data Source"
  )

ggsave(
  filename = paste0(wd_figure, "/TE_2018_2019.png"),
  width = 10,
  height = 8,
  units = "in",
  dpi = 150
)

# Rename some columns
# seminole_res$section <- "Upstream of Dam"
# chattahoochee_dam$section <- "Downstream of Dam"

# # Display side by side boxplots investigating effects of the Jim Woodruff Dam
# woodruffInfo <- rbind(seminole_res[c("SSC_mgL", "section")], chattahoochee_dam[c("SSC_mgL", "section")])

# # Calculate trapping efficiency using the median 
# # woodruffTEM <- (median(seminole_res$SSC_mgL) - median(chattahoochee_dam$SSC_mgL)) /median(seminole_res$SSC_mgL)
# # Calculate trapping efficiency using the average
# woodruffTEA <- (mean(seminole_res$SSC_mgL) - mean(chattahoochee_dam$SSC_mgL)) /mean(seminole_res$SSC_mgL)

# # Make plot
# ggplot(woodruffInfo, aes(x = section, y = SSC_mgL, fill=section)) +
#   geom_boxplot() +
#   stat_summary(fun="mean", linewidth = 2, size = 2) +
#   coord_cartesian(ylim = c(0, 70)) +
#   theme_classic() +
#   theme(legend.position = "none") +
#   labs(
#     title = "Jim Woodruff Lock and Dam",
#     x = "",
#     y = "SSC (mg/L)",
#     caption = paste0("\nTrapping Efficiency  = ", round(woodruffTEA,3))
#   )

# ggsave(
#   filename = paste0(wd_figure, "/WOODRUFF_TE.png"),
#   width = 20,
#   height = 20,
#   units = "cm"
# )

# All dams

# for landsat and sentinel, reduce to find the mean
# te_landsat_mean <- te_landsat %>% group_by(Dam.Name) %>%
#   summarize(te_mean = mean(te))
# te_sentinel_mean <- te_sentinel %>% group_by(Dam.Name) %>%
#   summarize(te_mean = mean(te, na.rm = TRUE))
# # get mean trap efficiency for the 
# meanTE <- data.frame(
#  name = c(te_landsat_mean$Dam.Name, "Jim Woodruff Lock and Dam"),
#  teSentinel = c(te_sentinel_mean$te_mean, woodruffTEA),
#  teLandsat = c(te_landsat_mean$te_mean, NA)
# )

# meanTE <- left_join(meanTE, DAMS %>% dplyr::select(Dam.Name, V_m3, Hydraulic.Height..Ft.), by=c("name"="Dam.Name"))
# meanTE$hydraulic_height_m <- meanTE$Hydraulic.Height..Ft. * 0.3048
# write.csv(meanTE,
#   file = paste0(wd_exports, "TE_MEANS.csv"), row.names = FALSE
# )

# ggplot(meanTE, aes(x = hydraulic_height_m, y = teLandsat)) +
#   geom_point() +
#   geom_point(aes(y = teSentinel))

########################### WHOLE TRANSECT #####################################

# All decades
ggplot() +
  # Landsat data
  stat_summary(
    data = landsat_flux,
    aes(x = distance_km, y = SSC_mgL, color='Landsat'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Landsat Error
  stat_summary(
    data = landsat_flux,
    aes(x = distance_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se, 
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#d52442", alpha = 0.3) +
  # Sentinel data
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = SSC_mgL, color='Sentinel'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Sentinel Error bar
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se, 
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#49b07b", alpha = 0.3) +
  # USGS data
  stat_summary(
    data = USGS_SSC %>% filter(year >= 1980),
    aes(x = dist_downstream_km, y = SSC_mgL, color='USGS'),
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # USGS Error
  stat_summary(
    data = USGS_SSC %>% filter(year >= 1980),
    aes(x = dist_downstream_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se, 
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "errorbar") +
  # Facet wrap by decade
  facet_wrap(~year_chunk) +
  coord_cartesian(ylim=c(0, 230)) +
  theme_bw() +
  labs(
    x = "Distance Downstream From Buford Dam (km)",
    y = "SSC (mg/L)",
    color = "Data Source"
  )

ggsave(
  filename = paste0(wd_figure, "/SSC_BYDECADE.png"),
  width = 10,
  height = 8,
  units = "in",
  dpi = 150
)

# 2018-now
ggplot() +
  # Landsat Error
  stat_summary(
    data = landsat_flux %>% filter(year >= 2018),
    aes(x = distance_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#d52442", alpha = 0.3) +
  # USGS data
  stat_summary(
    data = USGS_SSC %>% filter(year >= 2018),
    aes(x = dist_downstream_km, y = SSC_mgL, color='USGS'),
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Landsat Error
  stat_summary(
    data = USGS_SSC %>% filter(year >= 2018),
    aes(x = dist_downstream_km, y = SSC_mgL, group = 1), 
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "errorbar",
    width = 20) +
  # Sentinel data
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = SSC_mgL, color='Sentinel'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Sentinel Error bar
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#49b07b", alpha = 0.3) +
  # Landsat data
  stat_summary(
    data = landsat_flux %>% filter(year >= 2018),
    aes(x = distance_km, y = SSC_mgL, color='Landsat'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  coord_cartesian(ylim = c(0, 200)) +
  # Add Dam locations
  geom_vline(xintercept = DAMS$distance_km, linetype = "dotted") +
  theme_bw() +
  theme(
    legend.position = c(.05, .97),
    legend.justification = c("left", "top"),
    legend.box.background = element_rect(color = "black", size=0.75),
    legend.box.margin = margin(6, 6, 6, 6)) +
  labs(
    title = "Suspended Sediment Concentration from 2018 - Present",
    x = "Distance Downstream From Buford Dam (km)",
    y = "SSC (mg/L)",
    color = "Data Source"
  )

ggsave(
  filename = paste0(wd_figure, "/SSC_2018_2019.png"),
  width = 8,
  height = 6,
  units = "in",
  dpi = 200
)

# All decades flux
ggplot() +
  # Landsat Error
  stat_summary(
    data = landsat_flux,
    aes(x = distance_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#d52442", alpha = 0.3) +
  # USGS data
  stat_summary(
    data = USGS_SSC %>% filter(year >= 1980),
    aes(x = dist_downstream_km, y = SSC_mgL, color='USGS'),
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Landsat Error
  stat_summary(
    data = USGS_SSC %>% filter(year >= 1980),
    aes(x = dist_downstream_km, y = SSC_mgL, group = 1), 
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "errorbar",
    width = 20) +
  # Sentinel data
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = SSC_mgL, color='Sentinel'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Sentinel Error bar
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = SSC_mgL, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#49b07b", alpha = 0.3) +
  # Landsat data
  stat_summary(
    data = landsat_flux,
    aes(x = distance_km, y = SSC_mgL, color='Landsat'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  coord_cartesian(ylim = c(0, 200)) +
  # Add Dam locations
  # geom_vline(xintercept = DAMS$distance_km, linetype = "dotted") +
  facet_wrap(~year_chunk) + 
  theme_bw() +
  labs(
    title = "",
    x = "Distance Downstream From Buford Dam (km)",
    y = "SSC (mg/L)",
    color = "Data Source"
  )

ggsave(
  filename = paste0(wd_figure, "/SSC_BYDECADE.png"),
  width = 30,
  height = 20,
  units = "cm"
)

USGS_flux <- left_join(USGS_SSC, USGS_Q)
USGS_flux$flux_MTyr <- USGS_flux$Q_m3s * USGS_flux$SSC_mgL * 3.10585 * 10^-5

# 2018-now --> Flux
ggplot() +
  # Landsat data
  stat_summary(
    data = landsat_flux %>% filter(year >= 2018),
    aes(x = distance_km, y = flux_MTyr, color='Landsat'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +  
  # Landsat Error
  stat_summary(
    data = landsat_flux %>% filter(year >= 2018),
    aes(x = distance_km, y = flux_MTyr, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#9CDB79", alpha = 0.3) +
  # USGS data
  stat_summary(
    data = USGS_flux %>% filter(year >= 2018),
    aes(x = dist_downstream_km, y = flux_MTyr, color='USGS'),
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Landsat Error
  stat_summary(
    data = USGS_flux %>% filter(year >= 2018),
    aes(x = dist_downstream_km, y = flux_MTyr, group = 1), 
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "errorbar",
    width = 20) +
  # Sentinel data
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = flux_MTyr, color='Sentinel'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Sentinel Error bar
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = flux_MTyr, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ),
    geom = "ribbon", 
    fill = "#E69F00", alpha = 0.3) +
  coord_cartesian(ylim = c(0, 2)) +
  # Add Dam locations
  geom_vline(xintercept = DAMS$distance_km, linetype = "dotted") +
  scale_color_manual(values = c('#9CDB79','#E69F00', '#999999')) +
  theme_bw() +
  theme(
    legend.position = c(.05, .97),
    legend.justification = c("left", "top"),
    legend.box.background = element_rect(color = "black", size=0.75),
    legend.box.margin = margin(6, 6, 6, 6)) +
  labs(
    title = "Sediment Flux from 2018 - Present",
    x = "Distance Downstream From Buford Dam (km)",
    y = "Sediment Flux (Mt/yr)",
    color = "Data Source"
  )

ggsave(
  filename = paste0(wd_figure, "/FLUX_2018_2019.png"),
  width = 8,
  height = 6,
  units = "in",
  dpi = 200
)

# Flux by Decade
ggplot() +
  # Landsat data
  stat_summary(
    data = landsat_flux,
    aes(x = distance_km, y = flux_MTyr, color='Landsat'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +  
  # Landsat Error
  stat_summary(
    data = landsat_flux,
    aes(x = distance_km, y = flux_MTyr, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ), 
    geom = "ribbon", 
    fill = "#9CDB79", alpha = 0.3) +
  # USGS data
  stat_summary(
    data = USGS_flux %>% filter(year >= 1980),
    aes(x = dist_downstream_km, y = flux_MTyr, color='USGS'),
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Landsat Error
  stat_summary(
    data = USGS_flux %>% filter(year >= 1980),
    aes(x = dist_downstream_km, y = flux_MTyr, group = 1), 
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ), 
    geom = "errorbar",
    width = 20) +
  # Sentinel data
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = flux_MTyr, color='Sentinel'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    na.rm=T) +
  # Sentinel Error bar
  stat_summary(
    data = sentinel_flux,
    aes(x = distance_km, y = flux_MTyr, group = 1), 
    # fun.data = mean_se,
    fun.data = "mean_cl_boot",
    # fun.args = list(
    #   mult = 1
    # ), 
    geom = "ribbon", 
    fill = "#E69F00", alpha = 0.3) +
  coord_cartesian(ylim = c(0, 3)) +
  # Add Dam locations
  # geom_vline(xintercept = DAMS$distance_km, linetype = "dotted") +
  scale_color_manual(values = c('#9CDB79','#E69F00', '#999999')) +
  facet_wrap(~year_chunk) + 
  theme_bw() +
  labs(
    title = "",
    x = "Distance Downstream From Buford Dam (km)",
    y = "Sediment Flux (Mt/yr)",
    color = "Data Source"
  )

ggsave(
  filename = paste0(wd_figure, "/FLUX_BYDECADE.png"),
  width = 10,
  height = 8,
  units = "in",
  dpi = 150
)
