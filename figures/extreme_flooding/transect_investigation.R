### LIBRARY IMPORTS ###

library(ggplot2)
library(data.table)
library(matrixStats)
library(glmnet)
library(hydroGOF)
library(clue)
library(dplyr)

# Set root directory
wd_root <- "/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figures/extreme_flooding"
setwd(wd_root)

# Create events dataframe
names <- c(
  "Hurricane Delta",
  "Hurricane Helene",
  "Extreme Event",
  "Extreme Event",
  "Extreme Event",
  "Extreme Event"
)
startdates <- c(
  "2024-10-11",
  "2024-09-27",
  "2015-12-24",
  "2020-03-06",
  "2024-03-12",
  "2018-12-28"
)
events <- data.frame(
  names=names,
  date=as.Date(startdates)
)


chatM <- read.csv("/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figure-data/CHAT2.csv")
chatM$dist_group <- floor(chatM$DISTANCE_KM/10)*10
chatM$date <- as.Date(chatM$DATE)

median_ssc <- chatM %>%
  dplyr::group_by(dist_group) %>%
  dplyr::summarise(median_ssc = median(PRED_SSC_MGL, na.rm = TRUE))


# TESTING, we have data for events 4,6
target <- events$date[2]

chatM_filtered <- chatM %>%
  dplyr::filter(date - target >= 0,
                date - target <= 14)

chatM %>% dplyr::filter(abs(date - target) < 14)


ggplot(chatM_filtered, aes(group=dist_group, x=dist_group, y=PRED_SSC_MGL)) +
  geom_boxplot(aes(color='Sentinel-2'), outlier.shape=NA) +
  geom_line(data=median_ssc, aes(group=NA, x=dist_group, y=median_ssc)) + 
  #geom_boxplot(data=chatM, aes(group=dist_group, x=dist_group, y=PRED_SSC_MGL, color='REF')) +
  coord_cartesian(ylim = c(0, 150)) +
  theme_bw() +
  theme(
    legend.position = c(0.12, 0.9),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
  ) +
  scale_color_manual(values = c("#064012", "#3a3abb")) +
  guides(color = guide_legend(title = "Data Source")) +
  labs(
    x = "Distance Downstream of Buford Dam [km]",
    title="SSC Across the Chattahoocee and Apalachicola Rivers",
  )
