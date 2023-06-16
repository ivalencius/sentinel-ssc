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
library(mapview)
library(dataRetrieval)
Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/Pandoc")
# Sys.setenv(WEBSHOTJS="C:/Users/ilanv/AppData/Roaming/PhantomJS")

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
  width = 9,
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

ggplot(df, aes(x = cluster)) +
  geom_histogram(stat="count", color = "black", fill="#3bbecc") +
  # coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  # xlim(c(2, 4)) + # Try to remove that to see how to make a pie chart +
  labs(fill="Cluster",
       title="Number of Training Samples Available For Each Cluster"
      #  caption="Across every training SSC sample collected for Sentinel-2 Regression."
       ) +
  # THEMES
  theme_bw() +
  theme(text=element_text(family="Avenir"))
  # Remove background grid
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        # axis.text.y=element_blank())
  # scale_fill_manual(values=paletteer_dynamic("cartography::multi.pal",6))

ggsave(
  filename = paste0(wd_figure, "CLUSTER_COUNT.png"),
  height = 15,
  width = 15,
  units = "cm"
)

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
  height = 8,
  width = 8,
  units = "in",
  dpi = 150
)
########################### SAMPLE BASINS  #####################################

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
# mapshot(hucBasins, file = paste0(wd_figure, "test_map.png"))

# Now generate number of samples per GNIS ID
gnisSamps <- stations %>% 
  group_by(GNIS_nm) %>%
  summarize(
    totalSamples = sum(samples),
    totalSites = n()
  ) %>% arrange(-totalSamples)

write.csv(gnisSamps, file = paste0(wd_exports, "SAMPLES_PER_GNIS.csv"), row.names = FALSE)

# Combine upper/lower GNIS IDs
# Remove "Lower Prefix"
lowSections <- gnisSamps %>% filter(startsWith(GNIS_nm, "Lower"))
lowSections$GNIS_nm <- sapply(
  sapply(strsplit(lowSections$GNIS_nm, split = ' '), 
  function(x) x[2:length(x)]), paste, collapse = ' '
  )
# Remove "Upper prefix"
highSections <- gnisSamps %>% filter(startsWith(GNIS_nm, "Upper"))
highSections$GNIS_nm <- sapply(
  sapply(strsplit(highSections$GNIS_nm, split = ' '), 
  function(x) x[2:length(x)]), paste, collapse = ' '
  )
# Remove "Middle prefix"
middleSections <- gnisSamps %>% filter(startsWith(GNIS_nm, "Middle"))
middleSections$GNIS_nm <- sapply(
  sapply(strsplit(middleSections$GNIS_nm, split = ' '), 
  function(x) x[2:length(x)]), paste, collapse = ' '
  )
mainSections <- gnisSamps %>% 
filter(!(startsWith(GNIS_nm, "Lower") | startsWith(GNIS_nm, "Upper") | startsWith(GNIS_nm, "Middle")))
# Combine
mainStems <- bind_rows(lowSections, highSections, middleSections, mainSections)
# Remove different parts of the same river (i.e. Colorado-Cummins vs. Colorado-Marble)
mainStems$GNIS_nm <- sub("-.*", "", mainStems$GNIS_nm)
mainStemSamps <- mainStems %>% 
  group_by(GNIS_nm) %>%
  summarize(
    totalSamples = sum(totalSamples),
    totalSites = sum(totalSites)
  ) %>% arrange(-totalSamples)

write.csv(mainStemSamps, file = paste0(wd_exports, "SAMPLES_PER_RIVER.csv"), row.names = FALSE)

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

######################### Dams and USGS Locations ##############################

# Load ERDC Dams
raw_dams <- read.csv(paste0(wd_exports, "ERDC_CHATTAHOOCHEE_DAMS.csv"))
dams <- raw_dams %>% st_as_sf(coords=c("Longitude", "Latitude")) %>% st_set_crs(4326)
# Load USGS sites
raw_sites <- read.csv(paste0(wd_exports, "CHATTAHOOCHEE_STATIONINFO.csv"))
pad0 <- function(x) {
  while (nchar(x) < 8) {
    x <- paste0("0", x)
  }
  return(x)
}
USGS_siteinfo <- readNWISsite(lapply(unique(raw_sites$site_no), pad0))
sites <- data.frame(
  site_no = USGS_siteinfo$site_no,
  lat = as.numeric(USGS_siteinfo$dec_lat_va),
  lon = as.numeric(USGS_siteinfo$dec_long_va)
)
sites <- sites %>% st_as_sf(coords=c("lon", "lat")) %>% st_set_crs(4326)
# Load Chattahoochee centerline
centerline <- st_read("C:/Users/ilanv/Desktop/sentinel-ssc/imports/chattahoochee_centerline.kml")
# Load basin
riverData <- findNLDI(
    location = c(-85.03114856911581, 29.84561518630128),
    nav = "UM",
    find = c("basin"),
    distance_km = 2
    )
basin <- riverData$basin

# ch_basemap <- get_map(location=c(-85.68610, 29.54263, -83.10837, 35.92789),
#                       color="color",force=T) # manually input bbox
register_google("AIzaSyBipmtLXjlhfoLYiIgQmJKryq6BcP8gyiY")
ch_basemap <- get_map(
  # location=c(-85.68610, 29.54263, -83.10837, 34.9),
  "chattahoochee",
  zoom = 7, scale = "auto",
  maptype = "satellite",
  source = "google") # manually input bbox

# All dams
ggmap(ch_basemap) +
  geom_sf(data = basin, color = "Red", fill=NA, inherit.aes = FALSE) +
  geom_sf(data = centerline, color="Blue", inherit.aes = FALSE) +
  geom_sf(data = dams, aes(size = V_m3, color = V_m3, shape="Dam"),
    inherit.aes = FALSE) +
  geom_sf(data = sites, aes(shape = "USGS Site"),
    color = "green", 
    inherit.aes = FALSE) +
  geom_label_repel(data = raw_dams, aes(x=Longitude, y=Latitude,
    label = Dam.Name),
    size=3,
    nudge_y = 0.25,
    nudge_x = 0.5,
    inherit.aes = FALSE
    # family= "LM Roman 10"
    # family = "Helvetica"
    ) +
  # labs(color="NID Height Category",
      #  caption="The outline delineates the Chattahoochee\n-Apalachicola river basin. Labeled dams are \noperated by the US Army Corps of Engineers.") +
  # THEMES
  theme_bw() +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.background = element_rect(color = "black", size=0.55),
    legend.box.margin = margin(6, 6, 6, 6)) +
  guides(size = "none") +
  # scale_color_distiller(palette = "YlOrBr",
    # limits = c(0, 1000000, 2000000, 3000000),
    # labels = c(TeX("0 x $10^6$"), TeX("1 x $10^6$"), TeX("2 x $10^6$"), TeX("3 x $10^6$"))) +
  scale_color_continuous(
    type = "gradient",
    low = "#FFE5BD",
    high = "#FB9800",
    breaks = c(0, 1000000, 2000000, 3000000),
    labels = c(TeX("0 x $10^6$"), TeX("1 x $10^6$"), TeX("2 x $10^6$"), TeX("3 x $10^6$"))) +
  # scale_color_manual(values=paletteer_dynamic("cartography::multi.pal",3))
  # scale_color_paletteer_d("nord::lumina",-1)
  # scale_color_manual(values=paletteer_dynamic("cartography::wine.pal",3))
  labs(
    x = "Lon", y = "Lat", title = "Chattahoochee-Apalachicola Basin",
    color = "Reservoir Volume (m^3)", shape = "Feature"
  )

ggsave(
  filename = paste0(wd_figure, "BASIN.png"),
  width = 8,
  height = 8,
  units = "in",
  dpi = 150
)
