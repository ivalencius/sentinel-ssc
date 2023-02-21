# Code to make presentable plots

# Import packages
library(data.table)
library(dplyr)
library(scales)
library(units)
library(ggmap)
library(ggplot2)
library(ggrepel)
library(rjson)
library(gridExtra)
library(factoextra)
library(matrixStats)
library(psych)
library(clue)
library(GGally)
library(tidyverse)
library(hydroMap)
library(leaflet)
library(rgdal)
library(sf)
library(gridExtra)
library(ggpubr)
library(htmlwidgets)
library(htmltools)
# library(mapview)
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/Pandoc")

# Various themes
library(paletteer) # color palettes, 
library(latex2exp) # use TeX() to render latex

# Latex fonts
# library(extrafont)
# loadfonts(device = "win")
# font_import()
library(showtext)
font_add("Avenir", regular = "C:/Users/ilanv/AppData/Local/Microsoft/Windows/Fonts/AvenirNextLTPro-Regular.otf")
# font_add_google()
showtext_auto()

# Working Directory
wd_root <- "C:/Users/ilanv/Desktop/sentinel-ssc/chattahoochee-dams"
setwd(wd_root)

# Exports folder
wd_exports <- paste0(wd_root, "/exports/")
# Figures folder
wd_figure <- paste0(wd_exports, "figures/")
# Regress export folder
train_files <- "C:/Users/ilanv/Desktop/sentinel-ssc/exports/"

############################   TO-DO  ##########################################
# - Sentinel-landsat-USGS flux over time
# - Remake all old plots
#######################   RIVER GROUPING   #####################################
df <- read.csv(paste0(train_files, "/cluster_regress/RATING_TRAINING_CLUSTERED.csv"))
ssc_categories <- c(0,50,100,250,500,750,1e6)
# ssc_categories <- c(0,50,100,200,500,1e6)
# ssc_categories <- c(0,10,25,50,75,100,150,200,250,300,350, 400, 450, 500,600, 700, 800,900,1000,1100,1500, 1e6)

# Generate SSC labels as 'low value' - 'high value'
ssc_category_labels <- paste0(ssc_categories[-length(ssc_categories)],'-',c(ssc_categories[-1]))
# Make highest SSC category "> highest value"
ssc_category_labels[length(ssc_category_labels)] <- paste0('> ', ssc_categories[length(ssc_category_labels)])

df <- as.data.table(df)[
  ,':='(cluster_sel = cluster,
        # # Categorize SSC value as one of selected categories
        ssc_category = cut(10^log10_SSC_mgL, 
                           breaks = ssc_categories,
                           labels = ssc_category_labels))][]

# Generate median B,G,R for each SSC category and each cluster or site
ssc_category_color <- df %>% group_by(cluster_sel, ssc_category) %>% 
  summarise_at(vars("B4","B3","B2"), mean)

raster_color_types <- geom_raster(aes(fill = rgb(B4/2200,B3/2200,B2/2200))) # true color

ggplot(ssc_category_color, aes(x = cluster_sel, y = ssc_category)) +
  raster_color_types +
  scale_fill_identity() +
  # theme_classic()+
  # theme(axis.text.x = element_text(angle = 0)) + 
  labs(
    y = 'SSC range (mg/L)',
    x = 'River grouping') +
  # THEMES
  theme_bw() +
  # theme(text=element_text(family="LM Roman 10")) +
  theme(text=element_text(family="Avenir")) +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(
  filename = paste0(wd_figure, "CLUSTER_GROUPINGS.png"),
  width = 8,
  height = 8,
  units = "cm"
)

########################   CLUSTER CENTROIDS  ##################################
df <- read.csv(paste0(train_files, "/cluster_regress/RATING_TRAINING_CLUSTERED.csv"))

# load clusters (variable = cluster_calculated)
load(paste0(train_files, "/cluster_regress/clusters_Calculated.RData"))

# Bands used for regression
cluster_regressors <- c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')
# one row has one NA value --> filter it
# df <- df[, cluster_regressors]
# df <- df[is.finite(rowSums(df)), ]
cluster_medians <- data.table(matrix(nrow=0, ncol=length(cluster_regressors)+1))
names(cluster_medians) <- c('site_no', cluster_regressors)
unique_sites <- unique(as.character(df$site_no))
pb <- txtProgressBar(0, length(unique_sites), style = 3)
for (i in seq_along(unique_sites)){
  setTxtProgressBar(pb, i)
  site <- unique_sites[i]
  # Extract data
  site_data <- as.data.table(
    df[df$site_no == site,]
    )[, ..cluster_regressors] %>% as.matrix()
  # Reduce by median
  site_median <- colMedians(site_data, na.rm = TRUE) # avoid NA in any band column
  # Create table to hold site data
  site_df <- data.table(t(data.table(site_median)))
  names(site_df) <- cluster_regressors
  site_df$site_no <- site
  # Generate regression
  cluster_medians <- rbind(cluster_medians, site_df)
}

cluster_medians$cluster <- cl_predict(clusters_calculated, cluster_medians[, ..cluster_regressors])

cluster_medians$cluster <- as.character(cluster_medians$cluster)

# dataLong <- cluster_medians %>%
#   select(-site_no) %>%
#   pivot_longer(!cluster, names_to = "type", values_to = "count") 

# reduced <- cluster_medians %>%
#   group_by(cluster) %>%
#   summarise(meanB1 = mean(B1))

# ggplot(reduced, aes(x = cluster, y = meanB1)) +
#   # stat_summary(geom = "line") + 
#   # facet_wrap(cluster) +
#   # geom_density(alpha=0.4) +
#   geom_path() + 
#   theme_bw()
# upper.panel<-function(x, y){
#   points(x,y, pch=19, col=cluster_medians$site_no)
#   r <- round(cor(x, y), digits=2)
#   txt <- paste0("R = ", r)
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   text(0.5, 0.9, txt)
# }

# pairs(
#   cluster_medians[, ..cluster_regressors],
#   upper.panel = upper.panel,
#   # col = cluster_medians$site_no,
#   lower.panel=NULL
#   )
# pairs.panels(cluster_medians[, ..cluster_regressors], 
#              method = "pearson", # correlation method
#              hist.col = "#00AFBB",
#              density = TRUE,  # show density plots
#              ellipses = TRUE # show correlation ellipses
#              )

# ggpairs(
#   cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
#   aes(color = cluster, alpha = 0.4), 
#   columns = cluster_regressors,
#   upper = "blank"
#   # legends = TRUE
#   ) + 
#   theme_bw() +
#   theme(text=element_text(family="Avenir")) +
#   labs(title="Reflectances of Clustered Stations")

# # Cluster over PCA
# fviz_cluster(
#   object = clusters_calculated, 
#   data = cluster_medians[, ..cluster_regressors],
#   palette = "jco",
#   geom = "point",
#   ellipse.type = "convex", 
#   ggtheme = theme_bw()
#   )

B1 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B1), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B2 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B2), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B3 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B3), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B4 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B4), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B5 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B5), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B6 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B6), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B7 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B7), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B8 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B8), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B8A <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B8A), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B9 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B9), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B11 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B11), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")
B12 <- ggally_densityDiag(
  cluster_medians[, c("cluster",'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12')], 
  aes(color = cluster, x = B12), 
  alpha = 0.4) +
  theme_bw() +
  theme(legend.position="none")

gridded <- ggarrange(B1, B2, B3, B4, B5, B6, B7, B8, B8A, B9, B11, B12, nrow = 3, ncol = 4, 
  common.legend = TRUE, 
  legend="right")

ggsave(
  filename = paste0(wd_figure, "BAND_DENSITIES.png"),
  plot = gridded,
  width = 20,
  height = 14,
  units = "cm"
)
##########################   SAMPLE TYPE  ######################################

# df <- read.csv(paste0(train_files, "/cluster_regress/RATING_TRAINING_CLUSTERED.csv"))
# df[is.na(df$sample_method),] <- 'Q Rating Curve'
# samp_types <- df %>% group_by(sample_method) %>% 
#   tally() %>% 
#   filter(n>5) %>%
#   arrange(n)
# samp_types$sample_method <- factor(samp_types$sample_method, levels=samp_types$sample_method)
# samp_types$prop <- samp_types$n/nrow(df)
# samp_types$ymax <- cumsum(samp_types$prop)
# samp_types$ymin <- c(0, head(samp_types$ymax, n=-1))

# ggplot(samp_types, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=sample_method)) +
#   geom_rect() +
#   coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
#   xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart +
#   labs(fill="SSC Source Methodology",
#        # title="Proportion of Samples",
#        caption="Across every training SSC sample collected for Sentinel-2 Regression.") +
#   # THEMES
#   theme_void() +
#   theme(text=element_text(family="LM Roman 10")) +
#   # Remove background grid
#   # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#   #       axis.text.y=element_blank())+
#   scale_fill_manual(values=paletteer_dynamic("cartography::multi.pal",11))

##########################   CLUSTER TYPE  #####################################

df <- read.csv(paste0(train_files, "/cluster_regress/RATING_TRAINING_CLUSTERED.csv"))
df$cluster <- as.character(df$cluster)

ggplot(df, aes(x = cluster, fill = cluster)) +
  geom_histogram(stat="count", color = "black") +
  # coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  # xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart +
  labs(fill="Cluster",
       # title="Proportion of Samples",
       caption="Across every training SSC sample collected for Sentinel-2 Regression.") +
  # THEMES
  theme_bw() +
  theme(text=element_text(family="Avenir"))
  # Remove background grid
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # axis.text.y=element_blank())
  # scale_fill_manual(values=paletteer_dynamic("cartography::multi.pal",6))

ggplot(df, aes(x = cluster, group=1)) +
  stat_summary(aes(y = B1, color = "B1"), fun="median", geom="path") +
  stat_summary(aes(y = B2, color = "B2"), fun="median", geom="path") +
  stat_summary(aes(y = B3, color = "B3"), fun="median", geom="path") +
  stat_summary(aes(y = B4, color = "B4"), fun="median", geom="path") +
  stat_summary(aes(y = B5, color = "B5"), fun="median", geom="path") +
  stat_summary(aes(y = B6, color = "B6"), fun="median", geom="path") +
  stat_summary(aes(y = B7, color = "B7"), fun="median", geom="path") +
  stat_summary(aes(y = B8, color = "B8"), fun="median", geom="path") +
  stat_summary(aes(y = B8A, color = "B8A"), fun="median", geom="path") +
  stat_summary(aes(y = B9, color = "B9"), fun="median", geom="path") +
  stat_summary(aes(y = B11, color = "B11"), fun="median", geom="path") +
  stat_summary(aes(y = B12, color = "B12"), fun="median", geom="path") +
  theme_bw() +
  theme(text=element_text(family="Avenir")) +
  theme(panel.grid.major.x = 
    element_line(size=.1, color="black", linetype="dashed")) +
  theme(legend.position = "right") +
  theme(legend.key.height = unit(1, 'cm')) +
  labs(x = "Cluster", y="Reflectance", color = "Sentinel-2 Band")

ggsave(
  filename = paste0(wd_figure, "CLUSTER_REFLECTANCES.png"),
  height = 15,
  width = 15,
  units = "cm"
)
########################### SAMPLE BASINS  #####################################
# USE GNIS_ID TO GET RIVER NAME

df <- read.csv(paste0(train_files, "/cluster_regress/RATING_TRAINING_CLUSTERED.csv"))
station_df <- read.csv(paste0(train_files, "STATIONINFO.csv"))

stations <- read.csv(paste0(train_files, "GEE_Raw/SSC_stations_gte90m.csv"))
# stations$site_no <- as.character(stations$site_no)
# Joins stations
station_sites <- left_join(df, station_df, by = ("site_no"="site_no"))
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
station_sites$huc8 <- as.character(lapply(as.character(station_sites$HUC), pad0))
station_sites$site_no <- as.character(lapply(as.character(station_sites$site_no), pad0))
stations$site_no <- as.character(lapply(as.character(stations$site_no), pad0))
# Summarize training samples by basin
basinSamples <- station_sites %>% 
  group_by(huc8) %>%
  summarize(
    sites = length(unique(site_no))
  )
# Summarize training samples by sites
siteSamples <- station_sites %>% 
  group_by(site_no) %>%
  summarize(
    samples = n()
  )

stations <- left_join(stations, siteSamples, by = c("site_no"="site_no"))
stations <- stations %>% filter(samples >= 1)
stations$station_nm <- paste0(
  "<strong>", stations$sttn_nm, "</strong><br><strong>Site ID:</strong> USGS-", 
  stations$site_no,
  "<br><strong>Training Samples:</strong> ", stations$samples)
# Import basins
basins <- readOGR( 
  dsn= paste0(train_files, "/GEE_RAW/SSC_stations_HUC8Basins.shp"), 
  layer="SSC_stations_HUC8Basins",
  verbose=FALSE
) %>% st_as_sf()
basins <- left_join(basins, basinSamples, by = ("huc8"="huc8"))
basins <- basins %>% filter(sites >= 1)
basins$basin_nm <- paste0(
  "<strong>", basins$name, "</strong><br><strong>HUC8:</strong> ", basins$huc8, 
  "<br><strong>Area:</strong> ", round(as.numeric(basins$areasqkm)), " km^2")
reproj <- st_set_crs(basins, CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Create visualization
quantile(siteSamples$samples)
pal <- colorQuantile(
  palette = "Reds",
  domain = stations$samples,
  n = 4
)


tag.map.title <- tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(-50%,20%);
    position: fixed !important;
    left: 50%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 20px;
  }
"))

title <- tags$div(
  tag.map.title, HTML("HUC8 Basins Covered by Sentinel-2 Training Samples")
) 

labels <- c("1-3 (0-25%)", "3-11 (25-50%)", "12-86 (50-75%)", "87-315 (75-100%)")
hucBasins <- leaflet() %>% addTiles(
  urlTemplate = "http://mt0.google.com/vt/lyrs=m&hl=en&x={x}&y={y}&z={z}&s=Ga", 
  attribution = 'Google',
  ) %>%
  addPolygons(
    data = reproj, 
    stroke = FALSE,
    smoothFactor = 0.2,
    fillOpacity = 0.75,
    color = "blue",
    popup = ~basin_nm
    # color = ~pal(samples)
  ) %>% 
  addCircleMarkers(
    lng = stations$lon,
    lat = stations$lat,
    stroke = FALSE, 
    # color = "Black",
    color = pal(stations$samples),
    radius = 5,
    fillOpacity = 1,
    popup = stations$station_nm
  ) %>% addLegend(
    pal = pal, 
    values = stations$samples,
    # group = "circles",
    position = "topright",
    title = "Number of Samples (Quantiles)",
    labFormat = function(type, cuts, p) {paste0(labels)}) %>%
  addControl(title, position = "topleft", className="map-title")

# mapshot(hucBasins, url = paste0(wd_figure, "TRAINING_STATIONS.html"))
saveWidget(hucBasins, file = paste0(wd_figure, "TRAINING_STATIONS.html"))
############################ SSC STATIONS  #####################################
df <- read.csv(paste0(train_files, "/cluster_regress/SSC_stations_gte90m.csv"))

# Joins stations
station_sites <- left_join(df, station_df, by = ("site_no"="site_no"))
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
# Create visualization
pal <- colorQuantile(
  palette = "Reds",
  domain = reproj$samples,
  n = 4
)

leaflet() %>% addTiles(
  urlTemplate = "http://mt0.google.com/vt/lyrs=m&hl=en&x={x}&y={y}&z={z}&s=Ga", 
  attribution = 'Google') %>%
  addPolygons(
    data=reproj, 
    stroke = FALSE,
    smoothFactor = 0.2,
    fillOpacity = 1,
    color = ~pal(samples)
  ) %>% addLegend(pal = pal, values = basins$samples, group = "circles", position = "topright")

#########################   SSC/Q Sites   ######################################

# load riverdata
load("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\riverdata.Rdata")
# load Q stations
load("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\nwis_stations.Rdata")
# nwis_stations <- nwis_stations %>%
#   st_as_sf(coords=c('lon','lat')) %>%
#   st_set_crs(4326)
# load SSC stations
load("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\wqp_stations.Rdata")
# wqp_stations <- wqp_stations %>% 
  # st_as_sf(coords=c('lon','lat')) %>% 
  # st_set_crs(4326)

nwis_stations <- nwis_stations %>% select(-c('site_tp_cd','portal'))
wqp_stations <- wqp_stations %>% select(colnames(nwis_stations))

wqp_stations$lat <- as.numeric(wqp_stations$lat)
wqp_stations$lon <- as.numeric(wqp_stations$lon)
nwis_stations$lat <- as.numeric(nwis_stations$lat)
nwis_stations$lon <- as.numeric(nwis_stations$lon)

# Get stations which monitor both
SSCQ_stations <- union(nwis_stations, wqp_stations)

# Remove stations so nwis=Q only and wqp=SSC only
nwis_stations <- anti_join(SSCQ_stations, nwis_stations)
wqp_stations <- anti_join(SSCQ_stations, wqp_stations)

wqp_main <- wqp_stations %>% filter(channel == "MAIN STEM")
nwis_main <- nwis_stations %>% filter(channel == "MAIN STEM")
SSCQ_main <- SSCQ_stations %>% filter(channel == "MAIN STEM")

# Extract basemap
ch_bbox <- riverdata$basin %>% st_bbox()
ch_basemap <- get_map(location=c(-85.58610, 29.84263, -83.40837, 34.82789),
                      color="color",force=T) # manually input bbox

# All stations
ggmap(ch_basemap) +
  geom_sf(data = riverdata$basin, fill=NA, inherit.aes = FALSE) +
  geom_sf(data = riverdata$UM_flowlines, color="#1067b3", inherit.aes = FALSE) +
  geom_point(data = SSCQ_stations, aes(x=lon, y=lat, color = "Q/SSC Monitoring"), inherit.aes = FALSE) +
  geom_point(data = nwis_stations, aes(x=lon, y=lat, color = "Q Monitoring"), inherit.aes = FALSE) +
  geom_point(data = wqp_stations, aes(x=lon, y=lat, color = "SSC Monitoring"), inherit.aes = FALSE) +
  labs(color="USGS Stations",
       caption="The outline delineates the Chattahoochee/Apalachicola river basin.\nThe upstream extent of the starts at Buford Dam, GA.") +
  # THEMES
  theme_bw() +
  theme(text=element_text(family="LM Roman 10")) +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  scale_color_paletteer_d("nord::aurora")

# Only main stem stations
ggmap(ch_basemap) +
  geom_sf(data = riverdata$basin, fill=NA, inherit.aes = FALSE) +
  geom_sf(data = riverdata$UM_flowlines, inherit.aes = FALSE) +
  geom_point(data = SSCQ_main, aes(x=lon, y=lat, color = "Q/SSC Monitoring"), inherit.aes = FALSE) +
  geom_point(data = nwis_main, aes(x=lon, y=lat, color = "Q Monitoring"), inherit.aes = FALSE) +
  geom_point(data = wqp_main, aes(x=lon, y=lat, color = "SSC Monitoring"), inherit.aes = FALSE) +
  labs(color="USGS Stations",
       caption="The outline delineates the Chattahoochee/Apalachicola river basin.\nThe upstream extent of the study transect starts at Buford Dam, GA.") +
  # THEMES
  theme_bw() +
  theme(text=element_text(family="LM Roman 10")) +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_fill_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  scale_color_paletteer_d("nord::aurora")
  
  
#########################   Sentinel Sites   ###################################
sentinel <- read.csv("D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\exports\\GEE_raw\\rating_ssc_harmonized.csv") %>%
  select(lat,lon,mean_width,station_nm) %>%
  count(lat,lon,mean_width,station_nm) # n = number of samples

# Extract basemap
sentinel_bbox <- sentinel %>% st_bbox() 
map <- get_map(location=c(-125, 25.8, -69.6, 50)) # manually input bbox

ggmap(map) + 
  # coord_sf(crs = 4326) +
  geom_point(data = sentinel, aes(x=lon, y=lat, size=n, color=mean_width), inherit.aes = FALSE) +
  labs(caption="Data availability is from 2018-12-10 to 2022-07-20.",
       # title="USGS Sites Used for Sentinel-2 Regression",
       size="Number of Samples",
       color="Mean Width of Channel") +
  # THEMES
  theme_bw() +
  # theme(text=element_text(family="LM Roman 10")) +
  theme(text=element_text(family="Helvetica")) +
  # for no-axis labels
  theme(
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank()
  ) +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_paletteer_c("pals::ocean.tempo",-1)

#######################   Regression Coefficients   ############################

j <- fromJSON(file="D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_scaling\\train_ratingCurve_clust_ratingCurve\\regression\\lasso_full_bands\\info.JSON", simplify=TRUE)

# Load bands
bands <- c('B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11', 
          'B12','B2.B1','B3.B1','B4.B1','B1^2','B2^2','B3^2','B4^2',
          'B5^2','B6^2','B7^2','B8^2','B8A^2','B9^2','B11^2','B12^2',
          'B1^0.5','B2^0.5','B3^0.5','B4^0.5',
          'B5^0.5','B6^0.5','B7^0.5','B8^0.5','B8A^0.5','B9^0.5','B11^0.5','B12^0.5',
          'B2/B1','B3/B1','B4/B1','B5/B1','B6/B1','B7/B1','B8/B1',
          'B8A/B1','B9/B1','B11/B1','B12/B1','B3/B2','B4/B2','B5/B2',
          'B6/B2','B7/B2','B8/B2','B8A/B2','B9/B2','B11/B2','B12/B2',
          'B4/B3','B5/B3','B6/B3','B7/B3','B8/B3','B8A/B3','B9/B3',
          'B11/B3','B12/B3','B5/B4','B6/B4','B7/B4','B8/B4','B8A/B4',
          'B9/B4','B11/B4','B12/B4','B6/B5','B7/B5','B8/B5','B8A/B5',
          'B9/B5','B11/B5','B12/B5','B7/B6','B8/B6','B8A/B6','B9/B6',
          'B11/B6','B12/B6','B8/B7','B8A/B7','B9/B7','B11/B7','B12/B7',	
          'B8A/B8','B9/B8','B11/B8','B12/B8','B9/B8A','B11/B8A',
          'B12/B8A','B11/B9','B12/B9','B12/B11')

band_tmp <- integer(105)

# Get indices where ALL bands are not == 0
for (i in 1:length(bands)) {
  band_tmp <- j[["1_coefficients"]] + 
    j[["2_coefficients"]] + 
    j[["3_coefficients"]] +
    j[["4_coefficients"]] +
    j[["5_coefficients"]] +
    j[["6_coefficients"]]
}
band_mask <- band_tmp > 0

# Extract bands
b <- data.frame("1" = integer(105),
                "2" = integer(105),
                "3" = integer(105),
                "4" = integer(105),
                "5" = integer(105),
                "6" = integer(105))
rownames(b) <- bands
for (i in 1:length(bands)) {
  b[[1]][i] <- j[["1_coefficients"]][i]
  b[[2]][i] <- j[["2_coefficients"]][i]
  b[[3]][i] <- j[["3_coefficients"]][i]
  b[[4]][i] <- j[["4_coefficients"]][i]
  b[[5]][i] <- j[["5_coefficients"]][i]
  b[[6]][i] <- j[["6_coefficients"]][i]
}
bm <- data.frame("Cluster.1" = b[[1]][band_mask],
                 "Cluster.2" = b[[2]][band_mask],
                 "Cluster.3" = b[[3]][band_mask],
                 "Cluster.4" = b[[4]][band_mask],
                 "Cluster.5" = b[[5]][band_mask],
                 "Cluster.6" = b[[6]][band_mask])
rownames(bm) <- bands[band_mask]
bm <- format(bm, scientific = FALSE)
# bm <- format(bm, format="e", digits=4)
# bm[bm == 0] <- 0
write.csv(bm, "D:\\valencig\\Thesis\\Make_figs\\band_table.csv")

#########################   Satellite+USGS SSC   ###############################

load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "sentinel_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "SSC_full.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "landsat_day_avg.RData"))

SSC_full <- SSC_full %>% filter(channel == "MAIN STEM")
SSC_full$year_chunk <- cut(as.numeric(format(as.Date(SSC_full$date), "%Y")), 
                                   seq(from = 1985, to = 2026, by = 5), right=FALSE)
SSC_full$dist_downstream_km <- as.numeric(SSC_full$dist_downstream_km)
SSC_full$sscmg_L <- as.numeric(SSC_full$sscmg_L)

ggplot() + 
  stat_summary(
    data = landsat_day_avg,
    aes(x = distance_km, y = SSC_mgL, color='Landsat'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    # show.legend=TRUE,
    na.rm=T) +
  stat_summary(
    data = sentinel_day_avg,
    aes(x = distance_km, y = SSC_mgL, color='Sentinel-2'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    # show.legend=TRUE,
    na.rm=T) +
  stat_summary(
    data = SSC_full,
    aes(x = dist_downstream_km, y = sscmg_L, color='1 standard deviation'),
    geom = 'errorbar',
    fun.data= mean_sdl,
    fun.args = list(mult = 1), # mult = # of standard errors
    # size = 0.5,
    # width =  30,
    inherit.aes = FALSE,
    show.legend=FALSE,
    na.rm=T) +
  stat_summary(
    data = SSC_full,
    aes(x = dist_downstream_km, y = sscmg_L, color='USGS'),
    geom = 'point',
    fun= 'mean',
    #size = 3,
    #shape=24,
    inherit.aes = FALSE,
    # show.legend=TRUE,
    # fill = 'purple',
    na.rm=T) +
  stat_summary(
    data = SSC_full,
    aes(x = dist_downstream_km, y = sscmg_L, color='USGS'),
    geom = 'point',
    fun= 'mean',
    #size = 3,
    #shape=24,
    inherit.aes = FALSE,
    # show.legend=TRUE,
    na.rm=T) +
  facet_wrap(~year_chunk) +
  labs(
    # title="Suspended Sediment Concentration Along the Chattahoochee River",
       y="SSC (mg/L)",
       x="Distance from Buford Dam (km)",
       color="Data Source",
       caption="Landsat and Sentinel-2 data averaged to 1 sample per day per 2km.") +
  coord_cartesian(ylim = c(0, 275)) +
  # THEMES
  theme_minimal() +
  # theme(text=element_text(family="LM Roman 10")) +
  theme(text=element_text(family="Helvetica")) + 
  # Remove background grid
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_paletteer_d("nord::algoma_forest")
  scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",4))

#####################   Satellite+USGS 2015-present   ##########################

load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "sentinel_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "SSC_full.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "landsat_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "erdc_dams.RData"))

SSC_full <- SSC_full %>% filter(channel == "MAIN STEM")
SSC_full$year_chunk <- cut(as.numeric(format(as.Date(SSC_full$date), "%Y")), 
                           seq(from = 1985, to = 2026, by = 5), right=FALSE)
SSC_full$dist_downstream_km <- as.numeric(SSC_full$dist_downstream_km)
SSC_full$sscmg_L <- as.numeric(SSC_full$sscmg_L)

# Filter to 2015-present year chunk
sentinel_day_avg <- sentinel_day_avg %>% filter(year_chunk %in% c("[2015,2020)","[2020,2025)"))
landsat_day_avg <- landsat_day_avg %>% filter(year_chunk %in% c("[2015,2020)","[2020,2025)"))
SSC_full <- SSC_full %>% filter(year_chunk %in% c("[2015,2020)","[2020,2025)"))

ggplot() + 
  stat_summary(
    data = landsat_day_avg,
    aes(x = distance_km, y = SSC_mgL, color='Landsat'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    # show.legend=TRUE,
    na.rm=T) +
  stat_summary(
    data = sentinel_day_avg,
    aes(x = distance_km, y = SSC_mgL, color='Sentinel-2'),
    geom = 'line', 
    fun = 'mean',
    inherit.aes = FALSE,
    # show.legend=TRUE,
    na.rm=T) +
  stat_summary(
    data = SSC_full,
    aes(x = dist_downstream_km, y = sscmg_L, color='1 standard deviation'),
    geom = 'linerange',
    fun.data= mean_sdl,
    fun.args = list(mult = 1), # mult = # of standard errors
    # size = 0.5,
    # width =  30,
    inherit.aes = FALSE,
    # show.legend=FALSE,
    na.rm=T) +
  stat_summary(
    data = SSC_full,
    aes(x = dist_downstream_km, y = sscmg_L, color='USGS'),
    geom = 'point',
    fun= 'mean',
    #size = 3,
    #shape=24,
    inherit.aes = FALSE,
    # show.legend=TRUE,
    # fill = 'purple',
    na.rm=T) +
  stat_summary(
    data = SSC_full,
    aes(x = dist_downstream_km, y = sscmg_L, color='USGS'),
    geom = 'point',
    fun= 'mean',
    #size = 3,
    #shape=24,
    inherit.aes = FALSE,
    # show.legend=TRUE,
    na.rm=T) +
  labs(title="Suspended Sediment Concentration Along the Chattahoochee River (2015-present)",
       y="SSC (mg/L)",
       x="Distance from Buford Dam (km)",
       color="Data Source",
       caption="Landsat and Sentinel-2 data averaged to 1 sample per day per 2km.") +
  # Label dams 
  geom_label_repel(data = erdc_dams, 
                   aes(x=distance_km, label=Dam.Name, y=c(35-12, 30-12, 18, 25-12, 25-7)),
                   family="LM Roman 10",
                   nudge_y = 50,
                   size=3) +
  coord_cartesian(ylim = c(0, 100)) +
  # THEMES
  theme_minimal() +
  theme(text=element_text(family="LM Roman 10")) +
  # Remove background grid
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_paletteer_d("nord::algoma_forest")
  scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",4))

######################   Satellite Flux (2015-Present)  ########################
# Rating curves for flow and flux
# load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
#                  "landsat_regressed_flux.RData"))
# load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
#                  "sentinel_regressed_flux.RData"))
# 
# landsat_regressed_flux$SSC_flux_MTyr <- 10**landsat_regressed_flux$SSC_flux_MTyr_log10
# sentinel_regressed_flux$SSC_flux_MTyr <- 10**sentinel_regressed_flux$SSC_flux_MTyr_log10
# 
# # First row is null data, delete
# landsat_regressed_flux <- landsat_regressed_flux[-1,]
# sentinel_regressed_flux <- sentinel_regressed_flux[-1,]

# Rating curves just for flow
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "landsat_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "sentinel_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "erdc_dams.RData"))
erdc_dams <- erdc_dams[order(erdc_dams$distance_km),]
landsat_day_avg$SSC_flux_MTyr <- landsat_day_avg$SSC_flux_MTyr*1.1
sentinel_day_avg$SSC_flux_MTyr <- sentinel_day_avg$SSC_flux_MTyr*1.1

# Filter to 2015-present year chunk
sentinel_day_avg <- sentinel_day_avg %>% filter(year_chunk %in% c("[2015,2020)","[2020,2025)"))
landsat_day_avg <- landsat_day_avg %>% filter(year_chunk %in% c("[2015,2020)","[2020,2025)"))

# Increase flux by 10%
landsat_day_avg$SSC_flux_MTyr <- landsat_day_avg$SSC_flux_MTyr*1.1
sentinel_day_avg$SSC_flux_MTyr <- sentinel_day_avg$SSC_flux_MTyr*1.1

ggplot() +
  stat_summary(data = landsat_day_avg,
               aes(x = distance_km, y = SSC_flux_MTyr, color='Landsat'),
               geom = 'line',
               fun = 'mean',
               inherit.aes = FALSE,
               # show.legend=TRUE,
               na.rm=T) +
  stat_summary(data = sentinel_day_avg,
               aes(x = distance_km, y = SSC_flux_MTyr, color='Sentinel-2'),
               geom = 'line',
               fun = 'mean',
               inherit.aes = FALSE,
               # show.legend=TRUE,
               na.rm=T) +
  # Label dams 
  geom_label_repel(data = erdc_dams, 
                   aes(x=distance_km, label=Dam.Name, y=c(0, 0.05, 0.1, 0.45, 0.8)),
                   family="LM Roman 10",
                   nudge_y = .6,
                   size=3) +
  labs(
    # title="Satellite Derived Sediment Flux Along the Chattahoochee River (2015-Present)",
                                y="Sediment Flux (Mt/yr)",
                                x="Distance from Buford Dam (km)",
                                color="Data Source",
                                caption="Landsat and Sentinel-2 data averaged to 1 sample per day per 2km.\nFlux increased by 10% to account for Bed Load (Trimble, 1977)") +
  # THEMES
  theme_minimal() +
  # theme(text=element_text(family="LM Roman 10")) +
  theme(text=element_text(family="Helvetica")) + 
  # Remove background grid
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_paletteer_d("nord::algoma_forest")
  scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",2))

############################   Satellite Flux   ################################

# Rating curves just for flow
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "landsat_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "sentinel_day_avg.RData"))

# Increase flux by 10%
landsat_day_avg$SSC_flux_MTyr <- landsat_day_avg$SSC_flux_MTyr*1.1
sentinel_day_avg$SSC_flux_MTyr <- sentinel_day_avg$SSC_flux_MTyr*1.1

ggplot() +
  stat_summary(data = landsat_day_avg,
               aes(x = distance_km, y = SSC_flux_MTyr, color='Landsat'),
               geom = 'line',
               fun = 'mean',
               inherit.aes = FALSE,
               # show.legend=TRUE,
               na.rm=T) +
  stat_summary(data = sentinel_day_avg,
               aes(x = distance_km, y = SSC_flux_MTyr, color='Sentinel-2'),
               geom = 'line',
               fun = 'mean',
               inherit.aes = FALSE,
               # show.legend=TRUE,
               na.rm=T) +
  # stat_summary(data = landsat_day_avg,
  #              aes(x = distance_km, y = SSC_flux_MTyr, color='C'),
  #              geom = 'line',
  #              fun = 'mean',
  #              inherit.aes = FALSE,
  #              # show.legend=TRUE,
  #              na.rm=T) +
  # stat_summary(data = sentinel_regressed_flux,
  #              aes(x = distance_km, y = SSC_flux_MTyr, color='D'),
  #              geom = 'line',
  #              fun = 'mean',
  #              inherit.aes = FALSE,
  #              # show.legend=TRUE,
  #              na.rm=T) +
  facet_wrap(~year_chunk) +
  labs(
    title="Satellite Derived Sediment Flux Along the Chattahoochee River",
       y="Sediment Flux (Mt/yr)",
       x="Distance from Buford Dam (km)",
       color="Data Source",
       caption="Landsat and Sentinel-2 data averaged to 1 sample per day per 2km.\nFlux increased by 10% to account for Bed Load (Trimble, 1977)") +
  # THEMES
  theme_minimal() +
  theme(text=element_text(family="LM Roman 10")) +
  # Remove background grid
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_paletteer_d("nord::algoma_forest")
  scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",2))

############################## Dam Locations ###################################
# LOAD ERDC DAMS

# LOAD HERE

# load riverdata
load("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\riverdata.Rdata")
nid_dams <- read.csv("D:\\valencig\\Thesis\\Data\\NID_dataset\\NID_raw.csv")
nid_dams <- na.omit(nid_dams, cols=c('Longitude','Latitude'))
# Convert to sf object
nid_dams_sf <- st_as_sf(nid_dams, coords = c("Longitude", "Latitude")) %>% st_set_crs(4326)

# Get dams in the basin
basin_dams <- st_intersection(nid_dams_sf, riverdata$basin$geometry)
# Get their rows from nid_dams (geom_point needs lat and lon coordinates)
dams <- nid_dams %>% filter(NID.ID %in% basin_dams$NID.ID)
# One small dam (<25 feet) doesnt have a height category
dams[is.na(dams$NID.Height.Category),]$NID.Height.Category <- "Less than 25 feet"
dams$NID.Height.Category <- factor(dams$NID.Height.Category, levels=
                 c("Less than 25 feet",
                   "51-100 feet",
                   "Greater than 100 feet"))
erdc_dams <- dams %>% filter(Owner.Types == "Federal")
erdc_dams$V_m3 <- erdc_dams$Volume..Cubic.Yards. * 0.7646

# load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 # "erdc_dams.RData"))


# ch_basemap <- get_map(location=c(-85.68610, 29.54263, -83.10837, 35.92789),
#                       color="color",force=T) # manually input bbox
ch_basemap <- get_map(location=c(-85.68610, 29.54263, -83.10837, 34.9),
                      color="color",force=T) # manually input bbox

# All dams
ggmap(ch_basemap) +
  geom_sf(data = riverdata$basin, fill=NA, inherit.aes = FALSE) +
  geom_sf(data = riverdata$UM_flowlines, color="#1067b3", inherit.aes = FALSE) +
  geom_point(data = dams, aes(x=Longitude, y=Latitude, color=NID.Height.Category), 
             shape=17,
             inherit.aes = FALSE,
             size = 2) +
  geom_label_repel(data = erdc_dams, aes(x=Longitude, y=Latitude,
                label = Dam.Name),
                size=2,
                nudge_y = 0.25,
                nudge_x = 0.5,
                # family= "LM Roman 10"
                family = "Helvetica"
                ) +
  labs(color="NID Height Category",
       caption="The outline delineates the Chattahoochee\n-Apalachicola river basin. Labeled dams are \noperated by the US Army Corps of Engineers.") +
  # THEMES
  theme_bw() +
  # theme(text=element_text(family="LM Roman 10")) +
  theme(text=element_text(family="Helvetica")) +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # for no-axis labels
  theme(
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()
  ) +
  # scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  # scale_color_paletteer_d("nord::lumina",-1)
  scale_color_manual(values=paletteer_dynamic("cartography::wine.pal",3))

# Only ERDC Dams
ggmap(ch_basemap) +
  geom_sf(data = riverdata$basin, fill=NA, inherit.aes = FALSE) +
  geom_sf(data = riverdata$UM_flowlines, color="#1067b3", inherit.aes = FALSE) +
  geom_point(data = erdc_dams, aes(x=Longitude, y=Latitude, color=V_m3, size=V_m3), 
             shape=17,
             inherit.aes = FALSE) +
  geom_label_repel(data = erdc_dams, aes(x=Longitude, y=Latitude,
                                        label = Dam.Name),
                  size=3,
                  nudge_y = 0.2,
                  nudge_x = 0.5,
                  family= "LM Roman 10") +
  labs(color="Volume in Cubic Meters",
       caption="The outline delineates the Chattahoochee/Apalachicola river basin.\nThe upstream extent of the study transect starts at Buford Dam, GA.") +
  scale_size(guide="none") +
  # THEMES
  theme_bw() +
  theme(text=element_text(family="LM Roman 10")) +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # for no-axis labels
  theme(
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank()
  ) +
  # scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  # scale_color_paletteer_d("nord::lumina",-1)
  # scale_color_manual(values=paletteer_dynamic("cartography::wine.pal",3))
  scale_color_paletteer_c("pals::ocean.tempo")

############################### Data Density ###################################
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "ssc_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "landsat_mouth.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "sentinel_mouth.RData"))

# Add year to ssc_day_avg 
ssc_day_avg$year <- format(as.Date(ssc_day_avg$date), "%Y")

Samples <- data.frame(Source = c("USGS", "Landsat","Sentinel-2"),
                      Samples = c(nrow(ssc_day_avg), nrow(landsat_mouth), nrow(sentinel_mouth)))

ggplot(Samples, aes(x=Source, y=Samples, fill=Source)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label=Samples), vjust=-0.2) +
  theme_minimal() +
  labs(
    # title="Number of Samples at the Mouth of the Chattahoochee (1985-2023)",
       x = "",
       caption = "Samples are taken at site USGS-02359170: APALACHICOLA RIVER NR SUMATRA, FLA.\nFor Landsat/Sentinel, the mouth of the river is defined as the last 10km of the transect.") +
  geom_label(data = data.frame(mult=c("≈ 21.6x USGS","≈ 4x USGS"), 
                               Source = c("Landsat", "Sentinel-2")),
             aes(y=c(6000, 700),label=mult),
             label.size = 0,
             show.legend = F
             ) +
  # THEMES
  theme_bw() +
  # theme(text=element_text(family="LM Roman 10")) +
  theme(text=element_text(family="Helvetica")) +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  # scale_fill_paletteer_d("nord::lumina",-1)
  scale_fill_manual(values=paletteer_dynamic("cartography::harmo.pal",3))
  # scale_fill_paletteer_c("pals::ocean.tempo")


############################## Flux at Mouth ###################################

load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "ssc_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "landsat_mouth.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "sentinel_mouth.RData"))
# Increase flux by 10%
landsat_mouth$SSC_flux_MTyr <- landsat_mouth$SSC_flux_MTyr*1.1
sentinel_mouth$SSC_flux_MTyr <- sentinel_mouth$SSC_flux_MTyr*1.1
ssc_day_avg$SSC_flux_MTyr <- ssc_day_avg$SSC_flux_MTyr*1.1

# Total Data Density
ggplot(ssc_day_avg, aes(x=year_chunk, y=SSC_flux_MTyr, color="USGS")) +
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
  # scale_color_manual(values = c('#9CDB79', '#999999','#E69F00')) +
  labs(
    # title="Sediment Flux at the Mouth of the Chattahoochee",
       y="Sediment Flux (Mt/yr)",
       x="",
       color="Source",
       caption="USGS samples are taken at site USGS-02359170: APALACHICOLA RIVER NR SUMATRA, FLA.\nFor Landsat/Sentinel, the mouth of the river is defined as the last 10km of river.\nFlux increased by 10% to account for bed load (Trimble, 1977).") +
  # THEMES
  theme_bw() +
  # theme(text=element_text(family="LM Roman 10")) +
  theme(text=element_text(family="Helvetica")) +
  # Remove background grid
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  # scale_fill_paletteer_d("nord::lumina",-1)
  scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",3))

############################# Flux over Time ###################################
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "ssc_day_avg.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "landsat_mouth.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "sentinel_mouth.RData"))
load(file=paste0("D:\\valencig\\Thesis\\USGS_dataquery\\1985-01-01_2023-01-01\\variables\\",
                 "riverdata.RData"))
basin_area <- st_area(riverdata$basin) %>% set_units(km^2)
### PRE-ANTHROPOGENIC ###

# From octopusdata.org 
# Octupus ID: S156WTS010 Sample ID: SAP66
Be10E <- 9.15 %>% set_units(mm/kyr) %>% set_units(m/yr)
# Be10dE <- 2.15 %>% set_units(mm/kyr) %>% set_units(m/yr)

# Get total flux over one year (using density of quartz, 2.648 g/cm^3) --> switched to 1440 from Trimble 1977
Be10flux <- (Be10E * set_units(basin_area, m^2) * set_units(1440, kg/m^3)) %>% set_units(Mt/yr)
# Be10dflux <- (Be10dE * set_units(basin_area, m^2) * set_units(1440, kg/m^3)) %>% set_units(Mt/yr)

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
# USGSdflux <- sd(ssc_day_avg$SSC_flux_MTyr * 1.1) %>% set_units(Mt/yr)
# USGSQ <- mean(ssc_day_avg$flow_m3s) %>% set_units(m^3/s)
# USGSdQ <- sd(ssc_day_avg$flow_m3s) %>% set_units(m^3/s)
# USGSDisYield <- (USGSflux / basin_area) %>% set_units(t/(km^2*yr))
# USGSSSC <- mean(ssc_day_avg$sscmg_L) %>% set_units(mg/L)
# USGSdSSC <- sd(ssc_day_avg$sscmg_L) %>% set_units(mg/L)
# USGSsamp <- nrow(ssc_day_avg)

# Landsat Data (use last 10 km of the river)
# Increase by 10% to account for bed load
Landsatflux <- mean(landsat_mouth$SSC_flux_MTyr * 1.1, na.rm=T) %>% set_units(Mt/yr)
# Landsatdflux <- sd(landsat_mouth$SSC_flux_MTyr * 1.1, na.rm=T) %>% set_units(Mt/yr)
# LandsatQ <- mean(landsat_mouth$flow_m3s, na.rm=T) %>% set_units(m^3/s)
# LandsatdQ <- sd(landsat_mouth$flow_m3s, na.rm=T) %>% set_units(m^3/s)
# LandsatSSC <- mean(landsat_mouth$SSC_mgL, na.rm=T) %>% set_units(mg/L)
# LandsatdSSC <- mean(landsat_mouth$SSC_mgL, na.rm=T) %>% set_units(mg/L)
# Landsatsamp <- nrow(landsat_mouth)
# To use rating curve derived data uncomment below
# Landsatflux <- mean(mouth_q$rating_flux_landsat * 1.1, na.rm=T) %>% set_units(Mt/yr)
# Landsatflux5 <- mean(mouth_q[mouth_q$year_chunk %in% unique(sentinel_mouth$year_chunk),]$rating_flux_landsat * 1.1, na.rm=T) %>% set_units(Mt/yr)
# Landsatdflux <- sd(mouth_q$rating_flux_landsat * 1.1, na.rm=T) %>% set_units(Mt/yr)
# LandsatQ <- mean(mouth_q$flow_m3s, na.rm=T) %>% set_units(m^3/s)
# LandsatdQ <- sd(mouth_q$flow_m3s, na.rm=T) %>% set_units(m^3/s)


# Sentinel Data (use last 10 km of the river)
# Increase by 10% to account for bed load
Sentinelflux <- mean(sentinel_mouth$SSC_flux_MTyr * 1.1, na.rm=T) %>% set_units(Mt/yr)
# Sentineldflux <- sd(sentinel_mouth$SSC_flux_MTyr * 1.1, na.rm=T) %>% set_units(Mt/yr)
# SentinelQ <- mean(sentinel_mouth$flow_m3s, na.rm=T) %>% set_units(m^3/s)
# SentineldQ <- sd(sentinel_mouth$flow_m3s, na.rm=T) %>% set_units(m^3/s)
# SentinelSSC <- mean(sentinel_mouth$SSC_mgL, na.rm=T) %>% set_units(mg/L)
# SentineldSSC <- mean(sentinel_mouth$SSC_mgL, na.rm=T) %>% set_units(mg/L)
# Sentinelsamp <- nrow(sentinel_mouth)

# Data from Milliman et. al
# Increase by 10% to account for bed load
# MillimanQ <- (set_units(22, km^3/yr)) %>% set_units(m^3/s) # Average discharge
MillimanTDS <- set_units(1.1 * 1.1, Mt/yr) # TDS load
# MillimanDisYield <- set_units(0.021*1000, t/(km^2*yr)) # They use basin area of 52,000 km^2 but use 52 as divisor, need to multiply by 1,000
# MillimanSSC <- set_units(50, mg/L) # Average SSC concentration

# Flux Dataframe
Timeflux <- data.frame(source=c("Reusser (2015)", "Trimble (1977)", "USGS*", "Milliman (2013)*", "Landsat (Dethier, 2020)*","Sentinel-2**"),
                       flux=c(Be10flux, Tflux, USGSflux, MillimanTDS, Landsatflux, Sentinelflux),
                       period=c("Pre-Colonial", "Pre-Dam", "Post-Dam", "Post-Dam", "Post-Dam","Post-Dam"))
Timeflux$period <- factor(Timeflux$period, levels = c("Pre-Colonial", "Pre-Dam", "Post-Dam"))
Timeflux$source <- factor(Timeflux$source, levels = Timeflux$source)

ggplot(Timeflux, aes(x=source, y=flux, fill=period)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label=sprintf("%0.2f", round(flux, digits = 2))), vjust=-0.2) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    # title="Sediment Flux in the Chattahoochee Basin",
       x = "",
       y = "Sediment Flux",
       period  = "Time Period",
       caption = "Using bulk density of 1440 kg per cubic meter (Trimble, 1977) and assuming bed load accounts for an extra 10% across SSC derived flux.\n *Defined as 1985-present.\n**Data only available from 2018-2022.") +
  # THEMES
  theme_bw() +
  # theme(text=element_text(family="LM Roman 10")) +
  theme(text=element_text(family="Helvetica")) +
  # Remove background grid
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  # scale_fill_paletteer_d("nord::lumina",-1)
  # scale_fill_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  scale_fill_manual(values=c('#999999','#E69F00','#9CDB79'))
