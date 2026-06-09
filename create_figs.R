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
  labs(x = "Year", y = "Active USGS SSC/Turbidity Monitoring Sites") +
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
latlon <- rating_samples %>% dplyr::select(SITE_NO, LAT, LON) %>% 
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
    theme(legend.position = c(0.9, 0.2), legend.direction = "vertical") +
    theme(legend.box.background = element_rect(colour = "black", linewidth=1))  +
    theme(legend.key.size = unit(0.2, "cm"))

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
    theme(legend.position = c(0.9, 0.2), legend.direction = "vertical") +
    theme(legend.box.background = element_rect(colour = "black", linewidth=1)) +
    theme(legend.key.size = unit(0.2, "cm"))

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
    theme(legend.position = c(0.9, 0.2), legend.direction = "vertical") +
    theme(legend.box.background = element_rect(colour = "black", size=1)) +
    theme(legend.key.size = unit(0.2, "cm")) +
    theme(
      # Set the legend title size (e.g., 8 is quite small)
      legend.title = element_text(size = 5)
    ) +
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
#samples <- read.csv(paste0(wd_imports, "lag_4_INSITU_CLUSTERED.csv"))
samples <- read.csv(paste0(wd_imports, "lag_4_RATING_CLUSTERED.csv")) %>% 
  mutate(
    SOURCE = factor(
      SOURCE, c( "LANDSAT", "USGS RATING CURVE", "USGS")
    )
  )

sample_hist <- ggplot(samples, aes(log10(SSC_MGL), fill=SOURCE)) +
  geom_histogram(col="black", bins=40, alpha=1, position='stack') +
  theme_bw() +
  theme(legend.position = c(0.25, 0.92)) +
  theme(legend.box.background = element_rect(colour = "black", size=0.5)) +
  labs(x = "Suspended Sediment Concentration [log10(mg/L)]", y = "Number of Samples", fill=NULL) +
  scale_fill_manual(
    values = c(
      "USGS"       = "orange",
      "LANDSAT"  = "steelblue",
      "USGS RATING CURVE"       = "forestgreen"
    ),
    labels = c(
        "USGS"       = "In-situ samples",
        "LANDSAT"  = TeX("Landsat (Dethier et al., 2020)"),
        "USGS RATING CURVE"       = TeX("Rating curve derived (R$^2$>0.8)")
    )
  ) +
  theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),   # shrink label text
    legend.key.size = unit(0.15, "cm"),      # shrink key boxes (you already had this)
    legend.spacing.y = unit(0.05, "cm"),     # reduce vertical spacing
    legend.box.margin = margin(1, 1, 1, 1)   # tighten box padding
  )

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
    title = bquote("USGS Station " * .(station_sel) * " (R"^2*" = " * .(ex_R2) * ")"))
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
chatM <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/CHAT2.csv")
# Variable cluster
#chatV <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/CHATVAR.csv")

chatM <- chatM %>% group_by(DISTANCE_KM) %>% 
  reframe(
    L = quantile(PRED_SSC_MGL, probs=0.10), 
    H = quantile(PRED_SSC_MGL, probs=0.90),
    M = mean(PRED_SSC_MGL),
    d2 = lead(DISTANCE_KM)
    ) %>% filter(d2 >= 0) %>% unique()

ggplot(chatM, aes(x=DISTANCE_KM, y=M)) +
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

ggsave(filename=paste0(wd_exports, "chat.png"), height = 5, width = 10, units = "in", dpi=300)

##################### Fig 6 - Chattahoochee 2 ##################################
# Majority cluster
chatM <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/CHAT2.csv")
# Add label for every 2 km (helps with harmonization)
chatM$dist_group <- floor(chatM$DISTANCE_KM/2)*2
# USGS data
usgs <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/PROCESSED_CHATTAHOOCHEE_SSC.csv")
usgs$dist_group <- floor(usgs$dist_downstream_km/2) * 2

usgs_Q <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/PROCESSED_CHATTAHOOCHEE_Q.csv")
usgs_Q$dist_group <- floor(usgs_Q$dist_downstream_km/2) * 2
usgs_Q$dam_discharge_m3s <- 0

dams <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/ERDC_CHATTAHOOCHEE_DAMS.csv")

# Filter from data 1990 - now
usgs <- usgs %>% filter(sample_dt >= 1990)
usgs_Q <- usgs_Q %>% filter(sample_dt >= 1990)
usgs_Q_mean <- usgs_Q %>%
  dplyr::group_by(dist_group, sample_dt) %>%
  dplyr::summarise(Q_m3s = mean(Q_m3s, na.rm = TRUE))

# Active sites
active <- usgs %>% filter(sample_dt > 2022) %>% dplyr::select(dist_group) %>% unique()
Q_stations <- usgs_Q %>% filter(sample_dt > 2022) %>% dplyr::select(dist_group) %>% unique()

# Get SSC flux
# Create flux
df_flux <- left_join(
  chatM,
  usgs_Q_mean,
  by = c("DATE" = "sample_dt", "dist_group" = "dist_group")
)
# Ensure dist_group is sorted within each DATE
setDT(df_flux)
setorder(df_flux, DATE, dist_group)

# Add annual discharge for dams (data in \figures\trapping-efficiencies\)
#df_flux$dam_discharge_m3s <- NA
#df_flux <- df_flux %>%
#  mutate(dam_discharge_m3s = case_when(
#    dist_group == 10  ~ 1914  * 0.0283168,  # Buford (dist = 0)
#    dist_group == 240 ~ 4712  * 0.0283168,  # West Point dam (dist = 240)
#    dist_group == 430 ~ 9155  * 0.0283168,  # Walter F. George (dist = 432)
#    dist_group == 480 ~ 10767 * 0.0283168,  # George W. Andrews (dist = 478)
#    dist_group == 550 ~ 21489 * 0.0283168,  # Jim Woodruff (dist = 552)
#    .default = dam_discharge_m3s            # keep existing values for all other rows
#  ))

# Forward fill across dist_group for each DATE
df_flux[, Q_m3s := zoo::na.locf(Q_m3s, na.rm = FALSE), by = DATE]
# Filter flux to be 0 after Flint shows up
df_flux <- df_flux %>%
  mutate(Q_m3s = ifelse(dist_group > 520, NA, Q_m3s))

df_flux$flux_kg_s <- df_flux$PRED_SSC_MGL * df_flux$Q_m3s
#df_flux$dam_flux_kg_s <- df_flux$PRED_SSC_MGL * df_flux$dam_discharge_m3s

scale_factor <- 1 / 100

# Compute quartiles if needed
sentinel_summary <- df_flux %>%
  dplyr::group_by(dist_group) %>%
  dplyr::summarise(
    q1 = quantile(PRED_SSC_MGL, 0.25, na.rm = TRUE),
    median = median(PRED_SSC_MGL, na.rm = TRUE),
    q3 = quantile(PRED_SSC_MGL, 0.75, na.rm = TRUE),
    flux_median = median(flux_kg_s, na.rm=TRUE)
  )

ggplot(sentinel_summary, aes(x = dist_group, y = median)) +
  geom_ribbon(
    data = sentinel_summary,
    aes(ymin = q1, ymax = q3, group = 1, fill = "Sentinel-2"),
    alpha = 0.25
  ) +
  geom_line(
    data = sentinel_summary,
    aes(y = median, group = 1, color = "Sentinel-2"),
    linewidth = 1
  ) +
  geom_boxplot(
    data = usgs,
    aes(group=dist_group, x=dist_group, y = SSC_mgL, color = "USGS"),
    width = 20, outlier.shape = NA, alpha = 0.5, linewidth = 0.75
  ) +
  geom_point(data = Q_stations, aes(x = dist_group, y = 1, group=NA), 
           fill = "green", size = 5, shape = 24) +
  geom_point(data = active, aes(x = dist_group, y = 1, group=NA), 
             fill = "purple", size = 5, shape = 24) +
  geom_vline(xintercept = dams$distance_km, color = "red", linewidth = 1) +
  coord_cartesian(ylim = c(0, 75), xlim=c(0, 710), expand=FALSE) +
  scale_y_continuous(
    name = "SSC [mg/L]",
    sec.axis = sec_axis(~ . / scale_factor, name = "Flux [kg/s]")
  ) +
  geom_line(
    data=sentinel_summary, aes(x=dist_group, y=flux_median * scale_factor),
    linewidth=1, color='black'
    ) +
  theme_bw() +
  theme(
    legend.position = c(0.9, 0.85),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
  ) +
  theme(legend.box.background = element_rect(colour = "black", size=2)) +
  scale_color_manual(values = c("Sentinel-2" = "#064012", "USGS" = "#3a3abb")) +
  scale_fill_manual(values = c("Sentinel-2" = "#064012")) +
  guides(color = guide_legend(title = "Data Source"),
         fill = "none") +
  labs(
    x = "Distance Downstream of Buford Dam [km]",
    title = "SSC across the Chattahoocee and Apalachicola Rivers"
  ) 

ggsave(filename=paste0(wd_exports, "chat_boxplot.png"), height = 5, width = 7, units = "in", dpi=300)


calc_trapping_efficiency <- function(medians) {
  
  # Define dams with their upstream/downstream dist_group boundaries
  dam_definitions <- list(
    list(name = "West Point",    before_min = 0,   before_max = 240, after_min = 240, after_max = 240),
    list(name = "Walter F. George", before_min = 240, before_max = 432, after_min = 432, after_max = 442),
    list(name = "George W. Andrews", before_min = 432, before_max = 478, after_min = 478, after_max = 488)
    #list(name = "Jim Woodruff",  before_min = 478, before_max = 552, after_min = 552, after_max = Inf)
  )
  
  results <- lapply(dam_definitions, function(dam) {
    before_vals <- medians$median_value[
      medians$dist_group > dam$before_min & medians$dist_group < dam$before_max
    ]
    after_vals <- medians$median_value[
      medians$dist_group > dam$after_min & medians$dist_group < dam$after_max
    ]

    median_before <- max(before_vals, na.rm = TRUE)
    median_after  <- min(after_vals,  na.rm = TRUE)
    trapping_eff  <- (median_before - median_after) / median_before * 100
    
    data.frame(
      dam_name      = dam$name,
      median_before = median_before,
      median_after  = median_after,
      trapping_efficiency_pct = trapping_eff
    )
  })
  
  do.call(rbind, results)
}
medians <- df_flux %>% group_by(dist_group) %>% summarise(median_value = median(flux_kg_s, na.rm=TRUE))
medians_ssc <- df_flux %>% group_by(dist_group) %>% summarise(median_value = median(PRED_SSC_MGL, na.rm=TRUE),)
# Usage
te_flux <- calc_trapping_efficiency(medians)
te_ssc <- calc_trapping_efficiency(medians_ssc)
write.csv(
  te_flux,
  paste0(wd_imports, 'TE_flux_kgs.csv')
)
write.csv(
  te_ssc,
  paste0(wd_imports, 'TE_SSC_mgL.csv')
)

####################### Plot of regression coefficients ########################
df <- read.csv('/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/exports/cluster_regress/lag_4/RATING_regression/regression_variables.csv')
library(tidyverse)

df_long <- df %>%
  filter(Variable != "Intercept") %>% 
  pivot_longer(
    cols = starts_with("CLUSTER_"),
    names_to = "Cluster",
    values_to = "Value"
  ) %>%
  mutate(
    Cluster = str_replace(Cluster, "CLUSTER_", "Cluster "),
    Variable = str_replace_all(Variable, "\\.", "/")
  )

# Select top 10 by magnitude *within each cluster*
df_top10 <- df_long %>%
  mutate(mag = abs(Value)) %>%
  group_by(Cluster) %>%
  slice_max(order_by = mag, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Variable = reorder(paste(Variable, Cluster, sep = "___"), mag)
  )

ggplot(df_top10, aes(x = Variable, y = mag, fill = Cluster)) +
  geom_col(color = "black", width = 0.75) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  facet_wrap(~ Cluster, scales = "free") +
  coord_flip() + 
  scale_x_discrete(labels = function(x) sub("___.*$", "", x)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) +
  labs(
    x = "Top 15 Coefficients",
    y = "Magnitude"
  )


ggsave(filename = paste0(wd_exports, "top10_coefficients.png"),
       height = 6, width = 6, units = "in", dpi = 300)



