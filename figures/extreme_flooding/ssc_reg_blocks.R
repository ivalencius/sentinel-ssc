### LIBRARY IMPORTS ###

library(ggplot2)
library(data.table)
library(matrixStats)
library(glmnet)
library(hydroGOF)
library(clue)
library(dplyr)
library(raster)

# Set root directory
wd_root <- "/Users/ilanvalencius/Documents/River-Sed-Manuscript/sentinel-ssc/figures/extreme_flooding"
setwd(wd_root)

# Get clusters
load("../../exports/cluster_regress/lag_4/lag_4_RATING_clusters_calculated.RData")
cluster <- 2 # same cluster as the rest of the Chattahoochee (from apply-regressions.R)


# Load events
events_df <- read.csv("events.csv")

for (id in events_df$image_id) {
  cat("Processing image:", id, "\n")
  
  file <- paste0("raw-imgs/", id, ".tif")
  if (!file.exists(file)) {
    cat("  File not found, skipping:", file, "\n")
    next
  }
  
  # Load image - generated using GEE file Extreme-Flooding
  gtif <- stack(file)
  
  # Load the proper cluster model
  load(
    file = paste0(
      "../../exports/cluster_regress/lag_4/RATING_regression/ssc_lm_", cluster, ".RData"
    )
  )
  
  # Prepare output raster
  save_name <- paste0(wd_root, "/ssc/", id, "_ssc_cluster", cluster, ".tif")
  tmp_raster <- raster(gtif, 1)
  
  bs <- blockSize(gtif)
  tmp_raster <- writeStart(tmp_raster, filename = save_name, format = "GTiff",
                           overwrite = TRUE)
  
  # ── BLOCKWISE PROCESSING ─────────────────────────────────────
  for (i in 1:bs$n) {
    cat("  Block", i, "of", bs$n, "\n")
    
    # Read block values (all bands)
    v <- getValuesBlock(gtif, row = bs$row[i], nrows = bs$nrows[i])
    
    # If block is completely empty, skip
    if (is.null(v)) {
      next
    }
    
    # Convert to data.table and keep only the bands you use
    band_df <- as.data.table(v)
    
    # Ensure band names match your TIFF (adjust if needed)
    # Here assuming layers are named B2, B3, ..., B12 already
    setnames(band_df, old = names(band_df),
             new = names(gtif))  # sync names with stack
    
    band_df <- band_df[, .(B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12)]
    
    # Filter out rows that have NaN / NA data
    valid_rows <- which(complete.cases(band_df))
    n_block <- nrow(band_df)
    
    # Initialize prediction vector for this block
    ssc_block <- rep(NA_real_, n_block)
    
    if (length(valid_rows) > 0) {
      band_valid <- band_df[valid_rows, ]
      
      # ---- YOUR ORIGINAL FEATURE ENGINEERING ------------------
      fullVars <- setDT(copy(band_valid))[, ":="(
        # Add squared columns
        B2.2 = B2^2,
        B3.2 = B3^2,
        B4.2 = B4^2,
        B5.2 = B5^2,
        B6.2 = B6^2,
        B7.2 = B7^2,
        B8.2 = B8^2,
        B8A.2 = B8A^2,
        B11.2 = B11^2,
        B12.2 = B12^2,
        # Add square root columns
        B2.0.5 = B2^0.5,
        B3.0.5 = B3^0.5,
        B4.0.5 = B4^0.5,
        B5.0.5 = B5^0.5,
        B6.0.5 = B6^0.5,
        B7.0.5 = B7^0.5,
        B8.0.5 = B8^0.5,
        B8A.0.5 = B8A^0.5,
        B11.0.5 = B11^0.5,
        B12.0.5 = B12^0.5,
        
        # Add band ratios
        B3.B2  = B3/B2,
        B4.B2  = B4/B2,
        B5.B2  = B5/B2,
        B6.B2  = B6/B2,
        B7.B2  = B7/B2,
        B8.B2  = B8/B2,
        B8A.B2 = B8A/B2,
        B11.B2 = B11/B2,
        B12.B2 = B12/B2,
        
        B4.B3  = B4/B3,
        B5.B3  = B5/B3,
        B6.B3  = B6/B3,
        B7.B3  = B7/B3,
        B8.B3  = B8/B3,
        B8A.B3 = B8A/B3,
        B11.B3 = B11/B3,
        B12.B3 = B12/B3,
        
        B5.B4  = B5/B4,
        B6.B4  = B6/B4,
        B7.B4  = B7/B4,
        B8.B4  = B8/B4,
        B8A.B4 = B8A/B4,
        B11.B4 = B11/B4,
        B12.B4 = B12/B4,
        
        B6.B5  = B6/B5,
        B7.B5  = B7/B5,
        B8.B5  = B8/B5,
        B8A.B5 = B8A/B5,
        B11.B5 = B11/B5,
        B12.B5 = B12/B5,
        
        B7.B6  = B7/B6,
        B8.B6  = B8/B6,
        B8A.B6 = B8A/B6,
        B11.B6 = B11/B6,
        B12.B6 = B12/B6,
        
        B8.B7  = B8/B7,
        B8A.B7 = B8A/B7,
        B11.B7 = B11/B7,
        B12.B7 = B12/B7,
        
        B8A.B8 = B8A/B8,
        B11.B8 = B11/B8,
        B12.B8 = B12/B8,
        
        B11.B8A = B11/B8A,
        B12.B8A = B12/B8A,
        
        B12.B11 = B12/B11
      )]
      # ---------------------------------------------------------
      
      # Predict on valid rows only
      log10_ssc_valid <- predict(
        ssc_lm,
        newx = as.matrix(fullVars),
        s = "lambda.1se"
      )
      ssc_valid <- 10^log10_ssc_valid
      
      # Slot valid predictions back into block vector
      ssc_block[valid_rows] <- ssc_valid[, 1]
    }
    
    # Write this block’s predictions to the output raster
    tmp_raster <- writeValues(tmp_raster, ssc_block, bs$row[i])
  }
  
  tmp_raster <- writeStop(tmp_raster)
  cat("  Saved to:", save_name, "\n")
  
  gc()
}
