########################################

# The code below calculates a multiyear mean of the ratio data, estimating the dominant driver of seasonality in ice velocity (frontal vs fricitonal and generates a multiyear heatmap visualisation.

########################################



##### MULTI-GLACIER 5-YEAR MEAN PROCESSING & PLOTTING #####

# Parent directory containing all glacier folders
parent_directory <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2"

# Define the directory containing ratio CSV files
ratio_dir <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Outputs/Ratios"

# Path to the glacier name mapping CSV
glacier_info_path <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/glacier_name_jess.csv"
glacier_info <- read.csv(glacier_info_path)

# List all CSV files in the directory
ratio_files <- list.files(ratio_dir, pattern = "_ratio.csv$", full.names = TRUE)

# Extract glacier IDs from filenames and convert "GlX" → X
glacier_ids <- gsub("_ratio.csv", "", basename(ratio_files))  
numeric_glacier_ids <- suppressWarnings(as.numeric(gsub("Gl", "", glacier_ids)))  # Convert to numeric

# Sort glaciers numerically (ascending order)
sorted_indices <- order(numeric_glacier_ids, na.last = TRUE)  
ratio_files <- ratio_files[sorted_indices]
glacier_ids <- glacier_ids[sorted_indices]
numeric_glacier_ids <- numeric_glacier_ids[sorted_indices]

# Process each glacier: Compute 5-year mean if missing
for (i in seq_along(ratio_files)) {
  file_path <- ratio_files[i]
  glacier_id <- glacier_ids[i]
  
  # Read the CSV file
  ratio_data <- read.csv(file_path, row.names = 1, check.names = FALSE)
  
  # Check if "mean_2016_20" already exists
  if (!"mean_2016_20" %in% rownames(ratio_data)) {
    
    # Filter only the years 2016-2020
    valid_years <- suppressWarnings(as.numeric(rownames(ratio_data)))  
    ratio_subset <- ratio_data[valid_years %in% 2016:2020, , drop = FALSE]
    
    # Ensure valid data exists before proceeding
    if (nrow(ratio_subset) > 0) {
      # Calculate the 5-year mean for each elevation band
      mean_ratios <- colMeans(ratio_subset, na.rm = TRUE)
      
      # Append the mean as a new row labeled "mean_2016_20"
      ratio_data["mean_2016_20", ] <- mean_ratios
      
      # Save the updated matrix back to CSV
      write.csv(ratio_data, file = file_path, row.names = TRUE)
      
      cat("Updated ratio matrix with 5-year mean for:", glacier_id, "\n")
    } else {
      cat("Skipping:", glacier_id, "- No valid data for 2016-2020.\n")
    }
  }
}

##### PLOTTING THE MULTI-GLACIER HEATMAP #####

# Prepare a storage list for heatmap data
heatmap_data <- list()

# Reload updated CSV files
for (i in seq_along(ratio_files)) {
  file_path <- ratio_files[i]
  glacier_id <- glacier_ids[i]
  
  # Read the CSV file again (now with computed mean)
  ratio_data <- read.csv(file_path, row.names = 1, check.names = FALSE)
  
  # Ensure "mean_2016_20" exists
  if (!"mean_2016_20" %in% rownames(ratio_data)) {
    next  # Skip if 5-year mean row is missing
  }
  
  # Extract 5-year mean row
  mean_ratios <- as.numeric(ratio_data["mean_2016_20", ])
  names(mean_ratios) <- colnames(ratio_data)
  
  # Store in the list
  heatmap_data[[glacier_id]] <- mean_ratios
}

# Check if any valid data was collected
if (length(heatmap_data) == 0) {
  stop("No valid glacier data found! Ensure that at least one CSV contains 'mean_2016_20'.")
}

# Convert list to matrix (rows = glaciers, cols = elevation bands)
heatmap_matrix <- do.call(rbind, heatmap_data)

# Extract glacier names from mapping
glacier_names <- glacier_info$name[match(rownames(heatmap_matrix), glacier_info$ID)]

# Replace missing names with glacier IDs
glacier_names[is.na(glacier_names)] <- rownames(heatmap_matrix)

# Convert to log scale for better visualization
log_ratios <- log10(pmax(heatmap_matrix, 10^-1))  # Avoid -Inf by clipping small values

# **REVERSE ORDER SO LOWEST GLACIER NUMBERS ARE AT THE TOP**
log_ratios <- log_ratios[rev(seq_len(nrow(log_ratios))), ]
glacier_names <- rev(glacier_names)
glacier_ids <- rev(glacier_ids)

# Define colors (same as before)
breaks <- seq(-1, 1, length.out = 1000)  # Log scale
color_palette <- colorRampPalette(c("#364B9A", "#EAECCC", "#A50026"))(999)

# Define output file
pdf_file_path <- file.path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Outputs/Heatmaps/Mean",
                           "2016_20_hmm.pdf")

# Open PDF device
pdf(pdf_file_path, width = 8, height = max(4, length(glacier_ids) * 0.3))  # Adjust height dynamically

# Set up plotting layout (left: heatmaps, right: legend)
layout(matrix(c(1,2), nrow = 1), widths = c(4, 0.5)) 
par(mar = c(1.75, 0.2, 0.5, 0.5), oma =c(3, 4, 0.5, 0.75))

# Plot heatmap for all glaciers
image(1:length(colnames(heatmap_matrix)), 1:length(glacier_names), t(log_ratios),
      col = color_palette, breaks = breaks, axes = FALSE,
      xlab = "", ylab = "")

# Draw black borders around each glacier row
for (i in 1:length(glacier_names)) {
  rect(xleft = 0.5, xright = length(colnames(heatmap_matrix)) + 0.5,
       ybottom = i - 0.5, ytop = i + 0.5, 
       border = "black", lwd = 2)
}

# Add axis labels
axis(1, at = 1:length(colnames(heatmap_matrix)), labels = paste0(as.numeric(colnames(heatmap_matrix)), "m"), 
     las = 2, tck = -0.02, lwd = 0, lwd.ticks = 1, cex.axis = 1)
axis(2, at = 1:length(glacier_names), labels = glacier_ids, las = 1, tck = -0.02, lwd = 0, lwd.ticks = 1, cex.axis = 1)

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

# Save the raw 5-year means as CSV
csv_output_path <- file.path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Outputs/Heatmaps/Mean", 
                             "2016_20_hmm.csv")
write.csv(heatmap_matrix, file = csv_output_path, row.names = TRUE)
cat("Heatmap data saved at:", csv_output_path, "\n")

# Close PDF device
dev.off()

cat("Multi-glacier heatmap saved at:", pdf_file_path, "\n")


