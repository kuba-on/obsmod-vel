########################################
# The code below calculates a multiyear mean of the ratio data, estimating the dominant driver of seasonality in ice velocity (frontal vs fricitonal) and generates a multiyear heatmap visualisation.
########################################

##### MULTI-GLACIER MULTIYEAR MEAN PROCESSING & PLOTTING #####

# Parent directory containing all glacier folders
parent_directory <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3"

# Unified output directory (independent of friction law)
output_base_dir <- file.path(parent_directory, "Outputs", "fric_comp", "Multiyear Mean")
if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# Path to the glacier name mapping CSV
glacier_info_path <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/New Points v3/Input/Box Coordinates/box_sp_all_v3.csv"
glacier_info <- read.csv(glacier_info_path, check.names = FALSE, stringsAsFactors = FALSE)

# Expect EXACT column names in box_sp_all_v3.csv
# Required: feature_ID, glacier_ID, glacier_name, and coord_*_* columns
required_cols <- c("feature_ID","glacier_ID","glacier_name")
if (!all(required_cols %in% names(glacier_info))) {
  stop("box_sp_all_v3.csv must contain columns: feature_ID, glacier_ID, glacier_name, and coord_*_*.")
}

# Collect ALL coord_*_* columns once so rbind() works across mixed flowlines
all_coord_cols <- grep("^coord_\\d+_\\d+$", names(glacier_info), value = TRUE)

# Years to average (inclusive)
years_to_average <- 2016:2021

# Discover available friction-law directories under Outputs/, e.g., fric1, fric2, ...
outputs_dir <- file.path(parent_directory, "Outputs")
fric_dirs <- list.dirs(outputs_dir, full.names = TRUE, recursive = FALSE)
fric_dirs <- fric_dirs[grepl("^fric\\d+$", basename(fric_dirs))]
fric_ratio_dirs <- Filter(function(d) dir.exists(file.path(d, "Ratios")), fric_dirs)
if (length(fric_ratio_dirs) == 0) stop("No 'fric*/Ratios' directories found under Outputs/.")

# Colors and legend setup (keep style)
breaks <- seq(-1, 1, length.out = 1000)  # Log scale
color_palette <- colorRampPalette(c("#364B9A", "#EAECCC", "#A50026"))(999)

# Loop over each friction law directory
for (fric_path in fric_ratio_dirs) {
  
  fric_name <- basename(fric_path)                    # e.g., "fric2"
  fric_num  <- sub("^fric(\\d+)$", "\\1", fric_name)  # "2"
  ratio_dir <- file.path(fric_path, "Ratios")
  
  # Files named: gl_{featureID}_{flowlineID}_ratio_fric{n}.csv
  ratio_files <- list.files(
    ratio_dir,
    pattern = sprintf("^gl_\\d+_\\d+_ratio_fric%s\\.csv$", fric_num),
    full.names = TRUE
  )
  if (length(ratio_files) == 0) {
    cat("No ratio files found for", fric_name, "- skipping.\n")
    next
  }
  
  heatmap_data <- list()  # each row = one flowline
  meta_rows    <- list()  # per-row metadata + ratios
  
  for (file_path in ratio_files) {
    bname <- basename(file_path)
    
    # Parse feature_ID and flowline_ID from filename
    m <- regexec("^gl_(\\d+)_(\\d+)_ratio_fric(\\d+)\\.csv$", bname)
    mm <- regmatches(bname, m)[[1]]
    if (length(mm) != 4) { cat("Unexpected filename:", bname, "\n"); next }
    feature_id  <- as.integer(mm[2])
    flowline_id <- as.integer(mm[3])
    # fric_num already known
    
    # Read ratios
    ratio_data <- read.csv(file_path, row.names = 1, check.names = FALSE)
    
    # Compute 2016–2021 mean if missing
    if (!"mean_2016_21" %in% rownames(ratio_data)) {
      valid_years <- suppressWarnings(as.numeric(rownames(ratio_data)))
      ratio_subset <- ratio_data[valid_years %in% years_to_average, , drop = FALSE]
      if (nrow(ratio_subset) > 0) {
        mean_ratios <- colMeans(ratio_subset, na.rm = TRUE)
        ratio_data["mean_2016_21", ] <- mean_ratios
        write.csv(ratio_data, file = file_path, row.names = TRUE)
        cat("Updated 2016–2021 mean for:", bname, "\n")
      } else {
        cat("Skipping mean for:", bname, "- no valid data for 2016–2021.\n")
      }
    }
    if (!"mean_2016_21" %in% rownames(ratio_data)) next
    
    # Extract multiyear mean row
    mean_ratios <- as.numeric(ratio_data["mean_2016_21", ])
    names(mean_ratios) <- colnames(ratio_data)
    
    # Match metadata by feature_ID (one row per feature)
    gi_row <- glacier_info[glacier_info$feature_ID == feature_id, , drop = FALSE]
    if (nrow(gi_row) != 1) { cat("feature_ID", feature_id, "not uniquely found in box_sp_all_v3.csv. Skipping", bname, "\n"); next }
    
    glacier_ID   <- gi_row$glacier_ID[1]
    glacier_name <- gi_row$glacier_name[1]
    
    # ----- coords for THIS flowline (SECOND number = flowline_ID) -----
    # Columns like coord_1_3, coord_2_3 ... -> pick those with _{flowline_ID} at the end
    coord_cols_this <- grep(paste0("^coord_\\d+_", flowline_id, "$"), names(gi_row), value = TRUE)
    if (length(coord_cols_this) == 0) {
      cat("No coord columns found for flowline_ID:", flowline_id, "in file:", bname, "\n")
      next
    }
    
    # Build a 1-row data.frame with ALL coord_*_* columns (NA by default) so rbind names align
    coords_df <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = length(all_coord_cols)))
    names(coords_df) <- all_coord_cols
    # Fill only the columns for THIS flowline_ID
    coords_df[coord_cols_this] <- gi_row[coord_cols_this]
    
    # Unique key to keep stable ordering later
    storage_key <- sprintf("gID%08d_f%06d_fl%04d", as.integer(glacier_ID), as.integer(feature_id), as.integer(flowline_id))
    
    # Store for heatmap
    heatmap_data[[storage_key]] <- mean_ratios
    
    # Assemble row for CSV (metadata + ratios)
    meta_piece <- data.frame(
      glacier_ID   = glacier_ID,
      glacier_name = glacier_name,
      feature_ID   = feature_id,
      flowline_ID  = flowline_id,
      check.names  = FALSE,
      stringsAsFactors = FALSE
    )
    row_df <- cbind(meta_piece, coords_df, as.data.frame(as.list(mean_ratios), check.names = FALSE))
    row_df$..storage_key <- storage_key
    meta_rows[[storage_key]] <- row_df
  }
  
  if (length(heatmap_data) == 0) {
    cat("No valid glacier/flowline data found for", fric_name, ".\n")
    next
  }
  
  # Union of elevation-band columns across rows, ordered numerically when possible
  all_bands <- unique(unlist(lapply(heatmap_data, names)))
  suppressWarnings(bn <- as.numeric(all_bands))
  band_order <- if (all(!is.na(bn))) order(bn) else order(all_bands)
  eb_ordered <- all_bands[band_order]
  
  # Build matrix in consistent column order
  heatmap_matrix <- do.call(rbind, lapply(heatmap_data, function(v) v[match(eb_ordered, names(v))]))
  colnames(heatmap_matrix) <- eb_ordered
  rownames(heatmap_matrix) <- names(heatmap_data)
  
  # Build metadata frame matching row order
  meta_df <- do.call(rbind, meta_rows[rownames(heatmap_matrix)])
  rownames(meta_df) <- NULL
  
  # Order by ascending glacier_ID (then feature_ID, then flowline_ID for stability)
  ord_ix <- order(as.numeric(meta_df$glacier_ID),
                  as.numeric(meta_df$feature_ID),
                  as.numeric(meta_df$flowline_ID),
                  na.last = TRUE)
  meta_df <- meta_df[ord_ix, , drop = FALSE]
  heatmap_matrix <- heatmap_matrix[ord_ix, , drop = FALSE]
  
  # Labels: glacier_name, flowline_ID
  row_labels <- paste0(meta_df$glacier_name, ", ", meta_df$flowline_ID)
  
  # Log scale (clip to avoid -Inf)
  log_ratios <- log10(pmax(heatmap_matrix, 10^-1))
  
  # Reverse so smallest glacier_ID rows appear at the TOP (style kept consistent)
  log_ratios_plot <- log_ratios[rev(seq_len(nrow(log_ratios))), , drop = FALSE]
  row_labels_plot <- rev(row_labels)
  
  # Output filenames (independent of friction law directory)
  pdf_file_path   <- file.path(output_base_dir, sprintf("mean_hm_fric%s.pdf", fric_num))
  csv_output_path <- file.path(output_base_dir, sprintf("mean_ratio_fric%s.csv", fric_num))
  
  ##### PLOTTING THE MULTI-GLACIER HEATMAP #####
  # Open PDF device
  pdf(pdf_file_path, width = 8, height = max(4, nrow(log_ratios_plot) * 0.3))  # Adjust height dynamically
  
  # Set up plotting layout (left: heatmaps, right: legend)
  layout(matrix(c(1,2), nrow = 1), widths = c(4, 0.5)) 
  par(mar = c(1.75, 0.2, 0.5, 0.5), oma = c(3, 4, 0.5, 0.75))
  
  # Plot heatmap for all glaciers/flowlines
  image(
    1:length(colnames(log_ratios_plot)),
    1:length(row_labels_plot),
    t(log_ratios_plot),
    col = color_palette, breaks = breaks, axes = FALSE,
    xlab = "", ylab = ""
  )
  
  # Draw black borders around each row
  for (i in 1:length(row_labels_plot)) {
    rect(
      xleft = 0.5, xright = length(colnames(log_ratios_plot)) + 0.5,
      ybottom = i - 0.5, ytop = i + 0.5, 
      border = "black", lwd = 2
    )
  }
  
  # Axes (style unchanged)
  axis(1,
       at = 1:length(colnames(log_ratios_plot)),
       labels = paste0(as.numeric(colnames(log_ratios_plot)), "m"),
       las = 2, tck = -0.02, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  
  axis(2,
       at = 1:length(row_labels_plot),
       labels = row_labels_plot,
       las = 1, tck = -0.02, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  
  # Global labels
  mtext("glacier", side = 2, line = 3, outer = TRUE, cex = 1.1)
  mtext("elevation band", side = 1, line = 1.7, outer = TRUE, cex = 1.1)
  
  ### ADDING COLOR LEGEND ###
  par(mar = c(16, 0, 16, 2.5))  # Adjust margin for legend panel
  image(1, seq(-1, 1, length.out = 1000), matrix(seq(-1, 1, length.out = 1000), nrow = 1),
        col = color_palette, breaks = breaks, axes = FALSE, xlab = "", ylab = "")
  
  # Define log-scale labels and positions
  log_ticks <- seq(-1, 1, by = 1)  # Positions in log10 space
  log_labels <- parse(text = paste0("10^", log_ticks))  # Log-scale labels
  
  # Add log-scale ticks and labels
  axis(4, at = log_ticks, labels = log_labels, las = 1, cex.axis = 0.8, line = 0.1)
  mtext("velocity driver", side = 4, line = 2.25, cex = 0.8)
  mtext("friction                                     front", side = 4, line = -1.3, cex = 0.8)
  
  # Draw a black border around the legend
  rect(xleft = 0.6, xright = 1.4, ybottom = min(log_ticks), ytop = max(log_ticks),
       border = "black", lwd = 2)
  
  # Close PDF device
  dev.off()
  
  # ---- Save multiyear means with metadata ----
  # Keep metadata first (glacier_ID, glacier_name, feature_ID, flowline_ID, then ALL coord_*_*),
  # followed by elevation bands (already aligned in heatmap_matrix)
  out_df <- cbind(
    meta_df[, c("glacier_ID","glacier_name","feature_ID","flowline_ID"), drop = FALSE],
    as.data.frame(heatmap_matrix, check.names = FALSE)
  )
  
  # Already ordered by ascending glacier_ID; write CSV
  write.csv(out_df, file = csv_output_path, row.names = FALSE)
  
  cat("Heatmap saved at:", pdf_file_path, "\n")
  cat("Mean ratios (with metadata) saved at:", csv_output_path, "\n\n")
}

