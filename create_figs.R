# Code to make presentable plots

# Import packages
library(data.table)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(maps)
library(glmnet)
library(ggpubr)
library(scales)
library(hydroGOF)

# Various themes
library(paletteer) # color palettes, 
library(latex2exp) # use TeX() to render latex

# Latex fonts
library(extrafont)
# font_import()

# Working Directory
wd_root <- "/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc"
setwd(wd_root)
# Exports folder
wd_imports <- paste0(wd_root, "/figure-data/")
# Exports folder
wd_exports <- paste0(wd_root, "/figures/")

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- log10(l)
  return(parse(text = paste("10^", as.character(l), sep = "")))
}
##################  Significance of Regression  ################################
df <- read.csv('/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/RATING_regression/holdout_data.csv')

n = nrow(df)
names=c(rep("Holdout", n) , rep("Predicted", n))
# Draw samples from two different populations
value=c(df$SSC_MGL, 10^df$PRED_LOG10_SSC_MGL)
# Create a dataframe
data=data.frame(names,value)
t.test.lm = lm(value ~ names, data=data) 
summary(t.test.lm) # Means of two datasets are not different

##########################   Figure 1  #########################################
stations <- read.csv(paste0(wd_imports, "stations-over-time.csv"))

# Remove 2023 because 2023 not over yet 
stations <- stations %>% filter(year <= 2022)

coeff <- 1000
col1 <- "red"
col2 <- "blue"
ggplot(stations, aes(x = year, y = n_stations)) +
  geom_point(color = col1) +
  geom_line(color = col1) +
  geom_line(aes(y=n_samples/coeff), color=col2) +
  geom_point(aes(y=n_samples/coeff), color = col2) +
  labs(x = "Year", y = "Number of USGS-SSC Monitoring Sites") +
  theme_bw() +
  #theme(text=element_text(family="JetBrains Mono NL")) +
  scale_y_continuous(
     labels = scales::comma,
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Number of Samples Collected",
    labels = scales::label_number(scale_cut = scales::cut_short_scale()))) +
  theme(
    axis.title.y = element_text(color = col1, size=13),
    axis.title.y.right = element_text(color = col2, size=13)
  )

ggsave(
  filename = paste0(wd_exports, "stations-over-time.png"),
  height = 6,
  width = 6,
  units = "in")

##########################   Figure 2  #########################################
# Map of all training stations overlayed with distribution of data
# samples <- read.csv(paste0(wd_imports, "INSITU_TRAINING_CLUSTERED.csv"))
rating_samples <- read.csv(paste0(wd_imports, "lag_4_RATING_CLUSTERED.csv"))

# station_df <- read.csv(paste0(wd_imports, "STATIONINFO.csv"))
# rating_samples$LAT <- round(rating_samples$LAT, digits=3)
# rating_samples$LON <- round(rating_samples$LON, digits=3)
# Issue is lat lon are off due to floating point error
latlon <- rating_samples %>% select(SITE_NO, LAT, LON) %>% 
group_by(SITE_NO) %>% 
summarize_at(vars("LAT", "LON"), mean) %>% 
distinct()

# Summarize training samples by sites
samples_summed <- rating_samples %>% 
  group_by(SITE_NO) %>%
  summarize(
    n_samp = n(),
    mean_SSC = mean(SSC_MGL),
    cluster = mean(CLUSTER),
    mean_lag = mean(LAG_DAYS)
  )

stations <- left_join(samples_summed, latlon, by = "SITE_NO") %>%
  filter(is.finite(LAT))

# Get basemap
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
rivers50 <- ne_download(scale = "medium", type = "rivers_lake_centerlines", category = "physical", returnclass="sf")

mean_ssc_map <- ggplot(data = world) +
    geom_sf() +
    geom_sf(data=rivers50, col="Blue", alpha = 0.5) +
    geom_sf(data = states, fill = NA) +
    geom_point(data = stations, aes(x = LON, y = LAT, fill=mean_SSC), size = 4, 
        shape = 23) +
    scale_fill_gradientn(
      colors = c("#9DBF9E", "#FCB97D", "#A84268"),
      limits = c(0, 1500),
      # Set the out-of-bounds ("oob") rule to squish out-of-bounds values to the nearest limit
      oob = scales::squish,
      name = TeX("$\\bar{SSC}$ [mg/L]")) +
    coord_sf(xlim = c(-125, -67), ylim = c(25, 50), expand = FALSE) +
    theme_bw() +
    labs(x = "Longitude", y = "Latitude") +
    theme(legend.position = c(0.92, 0.2), legend.direction = "vertical") +
    theme(legend.box.background = element_rect(colour = "black", size=1))  +
    theme(legend.key.size = unit(0.2, "cm"))
    #theme(text = element_text(family = "JetBrains Mono NL"))

ggsave(filename = paste0(wd_exports, "sites-mean-ssc-main.png"))

# ggplot(data = world) +
#     geom_sf() +
#     geom_sf(data = states, fill = NA) +
#     geom_point(data = stations, aes(x = lon, y = lat, fill=mean_SSC), size = 27, 
#         shape = 23) + # size by approx. 6.6x to account for fact that plot will be 0.15 the size
#     scale_fill_gradientn(
#       colors = c("#9DBF9E", "#FCB97D", "#A84268"),
#       limits = c(0, 1500),
#       # Set the out-of-bounds ("oob") rule to squish out-of-bounds values to the nearest limit
#       oob = scales::squish) +
#     coord_sf(xlim = c(-155, -125), ylim = c(50, 71), expand = FALSE) +
#     theme_void() +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
#     theme(legend.position = "none") +
#     theme(text = element_text(family = "JetBrains Mono NL"))

# ggsave(filename = paste0(wd_exports, "sites-mean-ssc-inset.png"))


# Now get plot of number of samples

n_samp_map <- ggplot(data = world) +
    geom_sf() +
    geom_sf(data = rivers50, col="Blue", alpha = 0.5) +
    geom_sf(data = states, fill = NA) +
    geom_point(data = stations, aes(x = LON, y = LAT, fill=n_samp), size = 4, 
        shape = 23) +
    scale_fill_gradientn(
      colors = c("#9DBF9E", "#FCB97D", "#A84268"),
      limits = c(0, 300),
      # Set the out-of-bounds ("oob") rule to squish out-of-bounds values to the nearest limit
      oob = scales::squish,
      name = "# of Samples") +
    coord_sf(xlim = c(-125, -67), ylim = c(25, 50), expand = FALSE) +
    theme_bw() +
    labs(x = "Longitude", y = "Latitude") +
    theme(legend.position = c(0.92, 0.2), legend.direction = "vertical") +
    theme(legend.box.background = element_rect(colour = "black", size=1)) +
    theme(legend.key.size = unit(0.2, "cm"))
    # theme(legend.text = element_text(size = rel(2))) +
    # theme(legend.title = element_text(size = rel(3))) +
    # theme(axis.text = element_text(size = rel(2.5))) +
    # theme(axis.title = element_text(size = rel(4))) 
    #theme(text = element_text(family = "JetBrains Mono NL"))

ggsave(filename = paste0(wd_exports, "sites-n-samp-main.png"))

# ggplot(data = world) +
#     geom_sf() +
#     geom_sf(data = states, fill = NA) +
#     geom_point(data = stations, aes(x = LON, y = LAT, fill=n_samp), size = 27, 
#         shape = 23) + # size by approx. 6.6x to account for fact that plot will be 0.15 the size
#     scale_fill_gradientn(
#       colors = c("#9DBF9E", "#FCB97D", "#A84268"),
#       limits = c(0, 300),
#       # Set the out-of-bounds ("oob") rule to squish out-of-bounds values to the nearest limit
#       oob = scales::squish) +
#     coord_sf(xlim = c(-155, -125), ylim = c(50, 71), expand = FALSE) +
#     theme_void() +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
#     theme(legend.position = "none") +
#     theme(text = element_text(family = "JetBrains Mono NL"))

# ggsave(filename = paste0(wd_exports, "sites-n-samp-inset.png"))

# Now plot cluster locations

ggplot(data = world) +
    geom_sf() +
    geom_sf(data=rivers50, color="blue", size=0.5, alpha=0.5) +
    geom_sf(data = states, fill = NA) +
    geom_point(data = stations, aes(x = LON, y = LAT, fill=factor(cluster)), size = 4, 
        shape = 23) +
    coord_sf(xlim = c(-125, -67), ylim = c(25, 50), expand = FALSE) +
    theme_bw() +
    labs(x = "Longitude", y = "Latitude") +
    theme(legend.position = c(0.92, 0.2), legend.direction = "vertical") +
    theme(legend.box.background = element_rect(colour = "black", size=1)) +
    theme(legend.key.size = unit(0.1, "cm")) +
    #theme(text = element_text(family = "JetBrains Mono NL")) +
    labs(fill = "River Cluster")

ggsave(filename = paste0(wd_exports, "sites-cluster-main.png"))

# ggplot(data = world) +
#     geom_sf() +
#     geom_sf(data=rivers50, color="blue", size=0.5, alpha=0.5) +
#     geom_sf(data = states, fill = NA) +
#     geom_point(data = stations, aes(x = lon, y = lat, fill=factor(cluster)), size = 27, 
#         shape = 23) + # size by approx. 6.6x to account for fact that plot will be 0.15 the size
#     coord_sf(xlim = c(-155, -125), ylim = c(50, 71), expand = FALSE) +
#     theme_void() +
#     theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
#     theme(legend.position = "none") +
#     theme(text = element_text(family = "JetBrains Mono NL"))

# ggsave(filename = paste0(wd_exports, "sites-cluster-inset.png"))

#######################   Figure 2 - Hist  #####################################
samples <- read.csv(paste0(wd_imports, "lag_4_INSITU_CLUSTERED.csv"))
rating_samples <- read.csv(paste0(wd_imports, "lag_4_RATING_CLUSTERED.csv"))

sample_hist <- ggplot(rating_samples, aes(log10(SSC_MGL), fill="Rating Curve Derived")) +
  geom_histogram(col="black", bins=40, alpha=1) +
  geom_histogram(data=samples, aes(fill="In-Situ Samples"), col="black", bins=40, alpha=1) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.75)) +
  theme(legend.box.background = element_rect(colour = "black", size=1)) +
  labs(x = "Suspended Sediment Concentration [log10(mg/L)]", y = "Number of Samples", fill = "Data Source")
  #theme(text = element_text(family = "JetBrains Mono NL"))

ggsave(
   filename = paste0(wd_exports, "samples-hist.png"),
   height = 5,
   width = 5,
   unit = "in")

###################   Figure 2 - Rating Curve  #################################
ratingCurveData <- read.csv(paste0(wd_imports, "RATINGCURVE_SSC.csv"))
station_sel <- 10336610

rating_data_plot <- ratingCurveData %>% filter(SITE_NO == station_sel)
insitu_rating_data <- rating_data_plot %>% filter(RATINGCURVE == FALSE)

# ex_reg <- lm(log10_SSC_flux_MTyr~log10_Q_m3s, data=insitu_rating_data)
ex_R2 <- round(max(rating_data_plot$R2), digits=3)

rating_curve <- ggplot(insitu_rating_data, aes(x = LOG10_Q_M3S, y = LOG10_SSC_FLUX_MTYR)) +
  geom_point(color = "green", size=2) +
  geom_smooth(method="lm", color = "black", linetype = "dashed") +
  theme_bw() +
  labs(
    x = "River Discharge [log10(m^3/s)]",
    y = "Suspended Sediment Flux [log10(Mt/yr)]",
    title = paste0("USGS Station ", station_sel, " (R^2=", ex_R2,")"))
  #theme(text = element_text(family = "JetBrains Mono NL"))

ggsave(
   filename = paste0(wd_exports, "rating-curve.png"),
   height = 5,
   width = 5,
   units = "in")

################################################################################
# Arrange all figures
ggarrange(
  n_samp_map, rating_curve, mean_ssc_map, sample_hist,
  widths = c(2, 1), ncol = 2, nrow  = 2, 
  labels = c('(a)','(b)','(c)','(d)'))

ggsave(filename = paste0(wd_exports, "fig2.png"), height = 8, width = 12, units = "in", bg='white')

###################   Figure 3 - Color  #################################
df <- read.csv(paste0(wd_imports, "lag_4_RATING_CLUSTERED.csv"))
# ssc_categories <- c(0,50,100,250,500,750,1e6)
ssc_categories <- c(0,50,100,200,500,1e6)
# ssc_categories <- c(0,25,50,75,100,150,200,300,400,500,1e6)

# Generate SSC labels as 'low value' - 'high value'
ssc_category_labels <- paste0(ssc_categories[-length(ssc_categories)],'-',c(ssc_categories[-1]))
# Make highest SSC category "> highest value"
ssc_category_labels[length(ssc_category_labels)] <- paste0('> ', ssc_categories[length(ssc_category_labels)])

df <- as.data.table(df)[
  ,':='(cluster_sel = CLUSTER,
        # # Categorize SSC value as one of selected categories
        ssc_category = cut(SSC_MGL, 
                           breaks = ssc_categories,
                           labels = ssc_category_labels))][]

raster_color_types <- geom_raster(aes(fill = rgb(B4/2200,B3/2200,B2/2200))) # true color

# Investigate green for category 2!!
cat4green <- df %>% filter((SSC_MGL > 400) & (CLUSTER==2))
# cat4green <- df %>% filter(SITE_NO == 2175000)
cat4green_sites <- cat4green %>% group_by(SITE_NO) %>% 
  summarise_at(vars("B4","B3","B2"), mean)

ggplot(cat4green_sites, aes(x = factor(SITE_NO), y = 1)) +
  geom_raster(aes(fill = rgb(B4/2300,B3/2300,B2/2300))) +
  scale_fill_identity() +
  theme(legend.position = "none") 

# sites_of_interest <- c(2175000)

# Generate median B,G,R for each SSC category and each cluster or site
# ssc_category_color <- df %>% filter(!(SITE_NO %in% sites_of_interest)) %>% group_by(cluster_sel, ssc_category) %>% 
  # summarise_at(vars("B4","B3","B2"), mean)


ssc_category_color <- df %>% group_by(cluster_sel, ssc_category) %>% 
  summarise_at(vars("B4","B3","B2"), mean)

ggplot(ssc_category_color, aes(x = cluster_sel, y = ssc_category)) +
  raster_color_types +
  scale_fill_identity() +
  labs(
    y = 'SSC range [mg/L]',
    x = 'River Grouping') +
  # THEMES
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # Remove background grid
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #theme(text=element_text(family="JetBrains Mono NL"))

ggsave(
  filename = paste0(wd_exports, "cluster-groupings.png"),
  height = 5,
  width = 5,
  units = "in")

####################   Figure 4? -  REGRESSION #################################

# Step 1 - Need to generate predictions for all sites
# holdout <- read.csv("figure-data/RATING_regression/holdout_data.csv") %>% select(-starts_with("B"))
# df <- read.csv("figure-data/all_ssc_pred.csv")
df <- read.csv('/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/RATING_regression/holdout_data.csv')
df$LOG10_RESIDUAL <- df$PRED_LOG10_SSC_MGL - df$LOG10_SSC_MGL
# Residuals to do
# Width 
# Sample_depth_m
# Thin cirrus percentage
# Drainage_km2
# lag_days (?) - not sure if this is a good idea

r_width <- ggplot(df %>% filter(MERIT_WIDTH_M >= 60), aes(x = log10(MERIT_WIDTH_M), y = LOG10_RESIDUAL)) +
  geom_point(color = "black", alpha = 0.2) +
  stat_density_2d(
        aes(fill = ..level.., alpha = ..level..),
        bins = 10,
        geom = "polygon",
        colour = "black",
        alpha = 0.25,
        na.rm = TRUE
        ) +
  guides(fill = "none") +
  scale_fill_gradient(low = "black", high = "#00ffd971") +
  # geom_smooth(method="lm", color = "black", linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(-3,3)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "Red") +
  labs(
    x = "MERIT Hydro River Width [log10(m)]",
    y = "Model Residual\n[log10(predicted SSC) - log10(training sample SSC)]")
  #theme(text = element_text(family = "JetBrains Mono NL"))

# ggplot(df, aes(x = SAMPLE_DEPTH_M, y = LOG10_RESIDUAL)) +
#   geom_point(color = "black", alpha = 0.3) +
#   # geom_smooth(method="lm", color = "black", linetype = "dashed") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   coord_cartesian(ylim=c(-4,4)) +
#   geom_hline(yintercept = 0, linetype="dashed", color = "Red") +
#   labs(
#     x = "Sample Depth [m]",
#     y = "Model Residual\n[log10(predicted SSC) - log10(training sample SSC)]") +
#   theme(text = element_text(family = "JetBrains Mono NL"))

r_cirrus <- ggplot(df %>% filter(THIN_CIRRUS_PERCENTAGE > 10^-6), aes(x = log10(THIN_CIRRUS_PERCENTAGE), y = LOG10_RESIDUAL)) +
  geom_point(color = "black", alpha = 0.2) +
  stat_density_2d(
        aes(fill = ..level.., alpha = ..level..),
        bins = 10,
        geom = "polygon",
        colour = "black",
        alpha = 0.25,
        na.rm = TRUE
        ) +
  guides(fill = "none") +
  scale_fill_gradient(low = "black", high = "#00ffd971") +
  # geom_smooth(method="lm", color = "black", linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(-3,3)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "Red") +
  labs(
    x = "Thin Cirrus Percentage [log10(%)]",
    y = "Model Residual\n[log10(predicted SSC) - log10(training sample SSC)]")
  #heme(text = element_text(family = "JetBrains Mono NL"))


r_drainage <- ggplot(df, aes(x = log10(DRAINAGE_KM2), y = LOG10_RESIDUAL)) +
  geom_point(color = "black", alpha = 0.2) +
  stat_density_2d(
        aes(fill = ..level.., alpha = ..level..),
        bins = 10,
        geom = "polygon",
        colour = "black",
        alpha = 0.25,
        na.rm = TRUE
        ) +
  guides(fill = "none") +
  scale_fill_gradient(low = "black", high = "#00ffd971") +
  # geom_smooth(method="lm", color = "black", linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(-3,3)) +
  geom_hline(yintercept = 0, linetype="dashed", color = "Red") +
  labs(
    x = "Drainage Area [log10(km^2)]",
    y = "Model Residual\n[log10(predicted SSC) - log10(training sample SSC)]")
  #heme(text = element_text(family = "JetBrains Mono NL"))

df_factor <- df %>% mutate(LAG_DAYS = factor(abs(LAG_DAYS)),TYPE = factor(TYPE),
CLUSTER = factor(CLUSTER), SOURCE = factor(SOURCE),
SITE_TYPE = factor(SITE_TYPE))

# boxes <- c("LAG_DAYS", "TYPE", "CLUSTER", "SOURCE", "SITE_TYPE") 
# for (boxvar in boxes) {
#   foo <- ggplot(df_factor, aes(x = df_factor[[boxvar]], y = df_factor$LOG10_RESIDUAL)) +
#   geom_boxplot() +
#   # geom_point(color = "black", size=2, alpha = 0.3) +
#   geom_hline(yintercept = 0, linetype="dashed", color = "Red") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   coord_cartesian(ylim=c(-4,4)) +
#   geom_hline(yintercept = 0, linetype="dashed",  color = "Red") +
#   labs(
#     x = boxvar,
#     y = "Model Residual\n[log10(predicted SSC) - log10(training sample SSC)]") +
#   theme(text = element_text(family = "JetBrains Mono NL"))

#   ggsave(foo, filename = paste0(wd_exports, "boxplots/residuals-", boxvar, ".png"), 
#   height = 5, width = 5, units = "in")
# }

r_boxplot <- ggplot(df_factor, aes(x = SOURCE, y = LOG10_RESIDUAL)) +
  geom_boxplot() +
  # geom_point(color = "black", size=2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype="dashed", color = "Red") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(-3,3)) +
  geom_hline(yintercept = 0, linetype="dashed",  color = "Red") +
  labs(
    x = "Data Source",
    y = "Model Residual\n[log10(predicted SSC) - log10(training sample SSC)]")
  #theme(text = element_text(family = "JetBrains Mono NL"))

ggarrange(r_width, r_cirrus, r_drainage, r_boxplot, ncol = 2, nrow = 2, labels = c('(a)','(b)','(c)','(d)'))

ggsave(filename = paste0(wd_exports, "regression-residuals.png"), height = 10, width = 10, units = "in")

###################### Fig 5 - Model Bias ######################################
# df <- read.csv("figure-data/all_ssc_pred.csv")
df <- read.csv('/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/RATING_regression/holdout_data.csv')
# Investigate bias
totalBias <- function(df) {
  return(10^median(abs(log10(10^df$PRED_LOG10_SSC_MGL/10^df$LOG10_SSC_MGL)), na.rm = T)-1)
}
totalBias(df)

holdout <- df
holdout$BIAS <- (10^abs(log10(10^holdout$PRED_LOG10_SSC_MGL/10^holdout$LOG10_SSC_MGL)))-1
#abs(log10(10^holdout$PRED_LOG10_SSC_MGL/10^holdout$LOG10_SSC_MGL)-1)
# One invalid bias number - remove
holdout <- holdout %>% filter(is.finite(BIAS))

holdout_bias <- ggplot(holdout, aes(
  y = BIAS, 
  x="")) + 
  #geom_boxplot(color = "#1B6B1B") + 
  geom_boxplot(color='black') +
  # geom_jitter(width = 0.1, alpha = 0.3) +
  # coord_cartesian(ylim = c(0.001, 200)) +
  scale_y_log10(limits = c(0.001, 200), labels = fancy_scientific) + 
  theme_bw() +
  guides(color = "none") +
  geom_hline(yintercept = 0, linetype="dashed",  color = "Red") +
  labs(x = "Holdout Data", y = "Relative Bias")
  #theme(text = element_text(family = "JetBrains Mono NL"))

# cluster_bias <- ggplot(holdout, aes(
#   y = BIAS, 
#   x=factor(CLUSTER))) + 
#   geom_boxplot(aes(color = factor(SOURCE)), position = position_dodge(preserve = "single")) + 
#   # geom_jitter(width = 0.1, alpha = 0.3) +
#   # coord_cartesian(ylim = c(0.001, 200)) +
#   scale_y_log10(limits = c(0.001, 200), labels = fancy_scientific) + 
#   theme_bw() +
#   # guides(color = "none") +
#   scale_color_manual(values = c("#440154FF", "#2A788EFF", 
#   "#7AD151FF")) + 
#   theme(legend.position = c(0.8, 0.1)) +
#   geom_hline(yintercept = 0, linetype="dashed",  color = "Red") +
#   labs(x = "CLUSTER", y = "", color = "Data Source")
#   #theme(text = element_text(family = "JetBrains Mono NL"))

lag_bias <- ggplot(holdout, aes(y = BIAS, x=factor(abs(LAG_DAYS)))) +
    geom_boxplot(color='black', position = position_dodge(preserve = "single")) +
    # geom_jitter(width = 0.1, alpha = 0.3) +
    # coord_cartesian(ylim = c(0.001, 200)) +
    scale_y_log10(limits = c(0.001, 200), labels = fancy_scientific) +
    theme_bw() +
    # guides(color = "none") +
    theme(legend.position = "None") +
    geom_hline(yintercept = 0, linetype="dashed",  color = "Red") +
    labs(x = "|Lag Days|", y = "", color = "Data Source")


ggarrange(holdout_bias, lag_bias, ncol = 2, widths=c(1,3), labels = c('(a)','(b)'))

ggsave(filename=paste0(wd_exports, "lag-bias.png"), height = 7, width = 10, units = "in")

###################### Fig 6 - Chattahoochee ###################################
# Majority cluster
chatM <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/CHAT4.csv")
# Variable cluster
chatV <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/CHATVAR.csv")

chatM <- chatM %>% group_by(DISTANCE_KM) %>% 
  reframe(
    L = quantile(PRED_SSC_MGL, probs=0.10), 
    H = quantile(PRED_SSC_MGL, probs=0.90),
    M = mean(PRED_SSC_MGL),
    d2 = lead(DISTANCE_KM)
    ) %>% filter(d2 >= 0) %>% unique()

chatV <- chatV %>% group_by(DISTANCE_KM) %>% 
  reframe(
    L = quantile(PRED_SSC_MGL, probs=0.25), 
    H = quantile(PRED_SSC_MGL, probs=0.75),
    M = mean(PRED_SSC_MGL),
    d2 = lead(DISTANCE_KM)
    ) %>% filter(d2 >= 0) %>% unique()

m <- ggplot(chatM, aes(x=DISTANCE_KM, y=M)) +
  geom_ribbon(aes(xmin=DISTANCE_KM, xmax=d2, ymin=L, ymax=H, fill="Majority Cluster"), 
  alpha=0.3, color="#8fae8f") +
  geom_line() +
  scale_fill_manual(values = c("#064012")) +
  coord_cartesian(ylim = c(0, 200)) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.9)) +
  # guides(fill = guide_legend(title = "Cluster Assignment")) +
  guides(fill = "none") +
  labs(
    x = "Distance Downstream of Buford Dam [km]",
    y = "SSC [mg/L]",
    title="Majority Cluster")
  #theme(text = element_text(family = "JetBrains Mono NL"))

v <- ggplot(chatV, aes(x=DISTANCE_KM, y=M)) +
  geom_ribbon(aes(xmin=DISTANCE_KM, xmax=d2, ymin=L, ymax=H, fill="Variable Cluster"), 
  alpha=0.3, color="#3a3abb") +
  geom_line() +
  scale_fill_manual(values = c("#070457")) +
  coord_cartesian(ylim = c(0, 200)) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.9)) +
  # guides(fill = guide_legend(title = "Cluster Assignment")) +
  guides(fill = "none") +
  labs(
    x = "Distance Downstream of Buford Dam [km]",
    y = "SSC [mg/L]",
    title="Variable Cluster Assignment")
  #theme(text = element_text(family = "JetBrains Mono NL"))

ggarrange(v, m, ncol = 2, widths=c(1,1), labels = c('(a)','(b)'))
# ggsave(m, filename=paste0(wd_exports, "chatM.png"), height = 5, width = 7, units = "in", dpi=300)
ggsave(filename=paste0(wd_exports, "chat.png"), height = 5, width = 10, units = "in", dpi=300)

##################### Fig 6 - Chattahoochee 2 ##################################
# Majority cluster
chatM <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/CHAT4.csv")
# Add label for every 20 km
chatM$dist_group <- floor(chatM$DISTANCE_KM/10)*10
# USGS data
usgs <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/PROCESSED_CHATTAHOOCHEE_SSC.csv")
usgs$dist_group <- floor(usgs$dist_downstream_km/10) * 10
# Filter from data 2000 - now
usgs <- usgs %>% filter(sample_dt >= 1990)
# Load dams
dams <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/ERDC_CHATTAHOOCHEE_DAMS.csv")

# Active sites
active <- usgs %>% filter(sample_dt > 2022) %>% select(dist_downstream_km) %>% unique()
# chatM <- chatM %>% group_by(DISTANCE_KM) %>% 
#   reframe(
#     L = quantile(PRED_SSC_MGL, probs=0.10), 
#     H = quantile(PRED_SSC_MGL, probs=0.90),
#     M = mean(PRED_SSC_MGL),
#     d2 = lead(DISTANCE_KM)
#     ) %>% filter(d2 >= 0) %>% unique()

ggplot(chatM, aes(group=dist_group, x=dist_group, y=PRED_SSC_MGL)) +
  geom_boxplot(aes(color='Sentinel-2'), outlier.shape=NA) +
  geom_boxplot(data=usgs, aes(y=SSC_mgL, color='USGS'), width=20, outlier.shape = NA, alpha=0.5, linewidth=0.75) +
  geom_vline(xintercept = dams$distance_km, color="red", linewidth=1) +
  coord_cartesian(ylim = c(0, 150)) +
  geom_point(data=active, aes(x=dist_downstream_km, y=0, group=NA), fill="purple", size=5, shape=24) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.9)) +
  scale_color_manual(values = c("#064012", "#3a3abb")) +
  guides(color = guide_legend(title = "Data Source")) +
  labs(
    x = "Distance Downstream of Buford Dam [km]",
    y = "SSC [mg/L]",
    title="SSC Across the Chattahoocee and Apalachicola Rivers")
  #theme(text = element_text(family = "JetBrains Mono NL"))

ggsave(filename=paste0(wd_exports, "chat_boxplot.png"), height = 7, width = 10, units = "in", dpi=300)
