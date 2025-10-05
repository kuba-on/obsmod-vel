########################################

# The code below compares friction–frontal ratio matrices that use different friction laws (fric1, fric2, fric3) in a loop.
# The code computes the log10 fold-change between each pair of friction laws and generates heatmaps to visualize difference with regard to first friction law.
# The heatmaps use a diverging colour scale centered at 0 (no change), with blue denoting the first friction law being more friction-dominated (ratio < 0) and red denoting it being more frontal-dominated (ratio > 0).

########################################


options(stringsAsFactors = FALSE)

# Directories
dir_fric1 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric1/Ratios"
dir_fric2 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric2/Ratios"
dir_fric3 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric3/Ratios"

out_ratio_base   <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Ratios"
out_heatmap_base <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Heatmaps/Regular"
out_heatmap_masked_base <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Heatmaps/Masked"

# Glacier name lookup
csv_path <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/New Points v3/Input/Box Coordinates/box_sp_all_v3.csv"
glacier_info <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
gl_name_map <- setNames(glacier_info$glacier_name, as.character(glacier_info$feature_ID))

# Ensure output subfolders exist
pairs <- c("fric1_2","fric1_3","fric2_3")
for (p in pairs) {
  dir.create(file.path(out_ratio_base, p),    recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_heatmap_base, p),  recursive = TRUE, showWarnings = FALSE)
}

# Function to read a ratio CSV into a numeric matrix
read_ratio_matrix <- function(fp) {
  if (!file.exists(fp)) stop("Missing file: ", fp)
  df <- read.csv(fp, check.names = FALSE, row.names = 1)
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  return(mat)
}

# Function to parse glacier and flowline numbers from filename stem
parse_gl_fl <- function(stem) {
  m <- regexec("^gl_(\\d+)_(\\d+)_ratio$", stem)
  r <- regmatches(stem, m)[[1]]
  if (length(r) != 3) return(list(glacier_num = NA_integer_, flowline_num = NA_integer_))
  list(glacier_num = as.integer(r[2]), flowline_num = as.integer(r[3]))
}

# Helper to compute log10 fold-change matrix
logratio_matrix <- function(matA, matB, eps = .Machine$double.eps) {
  lr <- log10(matA / pmax(matB, eps))       # unified: log10(A/B), no pre-clipping
  lr[is.na(matA) | is.na(matB)] <- NA_real_ # preserve NA where inputs are NA
  lr
}

# Function to save log10 fold-change as CSV
save_difference_csv <- function(matA, matB, years, elevs, out_path) {
  stopifnot(identical(rownames(matA), years), identical(colnames(matA), elevs))
  stopifnot(identical(rownames(matB), years), identical(colnames(matB), elevs))
  lr_mat <- logratio_matrix(matA, matB)
  write.csv(lr_mat, out_path, row.names = TRUE)
}

# Heatmap plotting function (unmasked)
plot_heatmap_pair <- function(matA, matB, glacier_num, flowline_num, pair_label, out_pdf) {
  # unified values for plotting
  log_ratios <- logratio_matrix(matA, matB)
  years <- rownames(log_ratios)
  elevations <- colnames(log_ratios)
  
  # Set up PDF
  pdf(out_pdf, width = 9, height = 4)
  par(mar = c(0.2, 0.2, 0.2, 2), oma = c(4.5, 4, 1.5, 0.5))
  
  # Colour scale in log10 units, centred at 0
  breaks <- seq(-1, 1, length.out = 1000)
  color_palette <- colorRampPalette(c("#364B9A", "#EAECCC", "#A50026"))(999)
  
  # Reorder years top-to-bottom
  rev_years <- rev(years)
  log_ratios_reordered <- log_ratios[rev(seq_len(nrow(log_ratios))), , drop = FALSE]
  
  # Setup breaks (ensures consistent colour scale, plotting beyond range is clamped)
  lo <- min(breaks) + .Machine$double.eps
  hi <- max(breaks) - .Machine$double.eps
  lr_plot <- pmin(pmax(log_ratios_reordered, lo), hi)
  
  layout(matrix(c(1,2), nrow = 1), widths = c(4, 0.5))
  
  # Main panel
  image(1:length(elevations), 1:length(rev_years), t(lr_plot),
        col = color_palette, breaks = breaks, axes = FALSE,
        xlab = "", ylab = "")
  
  rect(xleft = 0.5, xright = length(elevations) + 0.5,
       ybottom = 0.5, ytop = length(rev_years) + 0.5,
       border = "black", lwd = 2)
  
  axis(1, at = 1:length(elevations), labels = paste0(elevations, "m"),
       las = 2, tck = -0.04, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  axis(2, at = 1:length(rev_years), labels = rev_years,
       las = 1, tck = -0.04, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  
  # Title
  glacier_name <- gl_name_map[as.character(glacier_num)]
  if (is.na(glacier_name) || length(glacier_name) == 0) glacier_name <- "Unknown glacier"
  mtext(paste0("Glacier ", glacier_num, " Flowline ", flowline_num,
               " - ", glacier_name, ", ", pair_label),
        side = 3, line = 0, outer = TRUE, cex = 1.3)

  # Global labels
  mtext("year of study", side = 2, line = 3, outer = TRUE, cex = 1.1)
  mtext("elevation band", side = 1, line = 3.5, outer = TRUE, cex = 1.1)      
  
  # Legend
  image(1, seq(-1, 1, length.out = 1000),
        matrix(seq(-1, 1, length.out = 1000), nrow = 1),
        col = color_palette, breaks = breaks, axes = FALSE,
        xlab = "", ylab = "")
  log_ticks <- seq(-1, 1, by = 1)
  axis(4, at = log_ticks, labels = log_ticks, las = 1, cex.axis = 0.8, line = 0.1)
  mtext("log-ratio", side = 4, line = 2.25, cex = 0.8)
  mtext("friction                                            front", side = 4, line = -1.3, cex = 0.8)
  
  rect(xleft = 0.6, xright = 1.4, ybottom = min(log_ticks), ytop = max(log_ticks),
       border = "black", lwd = 2)
  
  dev.off()
  cat("Heatmap saved:", out_pdf, "\n")
}

# Heatmap plotting function (masked)
plot_heatmap_pair_masked <- function(matA, matB, glacier_num, flowline_num, pair_label, out_pdf) {
  # unified values for plotting (same as CSV)
  log_ratios <- logratio_matrix(matA, matB)
  
  # Flip mask (cells where ratios A and B, obtained from two different friction laws, lie on different sides of 1)
  class_rel1 <- function(x) ifelse(is.na(x), NA_integer_,
                                   ifelse(x > 1, 1L, ifelse(x < 1, -1L, 0L)))
  clsA <- class_rel1(matA)
  clsB <- class_rel1(matB)
  flip_mask <- (clsA != clsB) & !is.na(clsA) & !is.na(clsB)
  
  # Valid data in both inputs
  valid_mask <- !is.na(matA) & !is.na(matB)
  # No-flip cells that still have data
  noflip_mask <- valid_mask & !flip_mask
  
  # Apply mask (keep only flips for the coloured layer; others NA)
  log_ratios[!flip_mask] <- NA_real_
  
  years <- rownames(log_ratios)
  elevations <- colnames(log_ratios)
  
  # Set up PDF
  pdf(out_pdf, width = 9, height = 4)
  par(mar = c(0.2, 0.2, 0.2, 2), oma = c(4.5, 4, 1.5, 0.5))
  
  # Colour scale in log10 units, centered at 0
  breaks <- seq(-1, 1, length.out = 1000)
  color_palette <- colorRampPalette(c("#364B9A", "#EAECCC", "#A50026"))(999)
  
  # Reorder years top-to-bottom (and reorder masks consistently)
  rev_years <- rev(years)
  rev_idx <- rev(seq_len(nrow(log_ratios)))
  log_ratios_reordered <- log_ratios[rev_idx, , drop = FALSE]
  noflip_mask_reordered <- noflip_mask[rev_idx, , drop = FALSE]
  
  # Setup breaks (ensures consistent colour scale, plotting beyond range is clamped)
  lo <- min(breaks) + .Machine$double.eps
  hi <- max(breaks) - .Machine$double.eps
  lr_plot <- pmin(pmax(log_ratios_reordered, lo), hi)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(4, 0.5))
  
  # Main panel
  # Draws only the "no-flip" cells as grey; leaves true NA (no data) cells white by not drawing them.
  noflip_grid <- matrix(NA_real_, nrow = nrow(noflip_mask_reordered), ncol = ncol(noflip_mask_reordered))
  noflip_grid[noflip_mask_reordered] <- 0
  image(1:length(elevations), 1:length(rev_years), t(noflip_grid),
        col = "grey85", breaks = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
  
  # Overlay the colored flips
  image(1:length(elevations), 1:length(rev_years), t(lr_plot),
        col = color_palette, breaks = breaks, axes = FALSE, add = TRUE)
  
  rect(xleft = 0.5, xright = length(elevations) + 0.5,
       ybottom = 0.5, ytop = length(rev_years) + 0.5, border = "black", lwd = 2)
  
  axis(1, at = 1:length(elevations), labels = paste0(elevations, "m"),
       las = 2, tck = -0.04, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  axis(2, at = 1:length(rev_years), labels = rev_years,
       las = 1, tck = -0.04, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  
  # Title
  glacier_name <- gl_name_map[as.character(glacier_num)]
  if (is.na(glacier_name) || length(glacier_name) == 0) glacier_name <- "Unknown glacier"
  mtext(paste0("Glacier ", glacier_num, " Flowline ", flowline_num,
               " — ", glacier_name, ", ", pair_label),
        side = 3, line = 0, outer = TRUE, cex = 1.1)
  
  # Global labels
  mtext("year of study", side = 2, line = 3, outer = TRUE, cex = 1.1)
  mtext("elevation band", side = 1, line = 3.5, outer = TRUE, cex = 1.1)
  
  # Legend
  par(mar = c(0.2, 0, 0.2, 3))
  image(1, seq(-1, 1, length.out = 1000),
        matrix(seq(-1, 1, length.out = 1000), nrow = 1),
        col = color_palette, breaks = breaks, axes = FALSE, xlab = "", ylab = "")
  log_ticks <- seq(-1, 1, by = 1)
  axis(4, at = log_ticks, labels = log_ticks, las = 1, cex.axis = 0.8, line = 0.1)
  mtext("log-ratio", side = 4, line = 2.25, cex = 0.8)
  mtext("friction                                            front", side = 4, line = -1.3, cex = 0.8)
  
  rect(xleft = 0.6, xright = 1.4, ybottom = min(log_ticks), ytop = max(log_ticks),
       border = "black", lwd = 2)
  
  dev.off()
  cat("Heatmap saved:", out_pdf, "\n")
}


# Main processing loop
all_f1 <- list.files(dir_fric1,
                     pattern = "^gl_\\d+_\\d+_ratio_fric1\\.csv$",
                     full.names = TRUE)

if (length(all_f1) == 0L) {
  stop("No fric1 CSVs found in: ", dir_fric1)
}

for (fp1 in all_f1) {
  base1 <- basename(fp1)
  stem  <- sub("_fric1\\.csv$", "", base1)   # e.g., "gl_0_1_ratio"
  ids   <- parse_gl_fl(stem)
  
  if (is.na(ids$glacier_num) || is.na(ids$flowline_num)) {
    message("Skipping unrecognized file name: ", base1)
    next
  }
  
  # Matching files in fric2 / fric3
  fp2 <- file.path(dir_fric2, paste0(stem, "_fric2.csv"))
  fp3 <- file.path(dir_fric3, paste0(stem, "_fric3.csv"))
  
  if (!file.exists(fp2) || !file.exists(fp3)) {
    message("Missing companion file(s) for ", stem,
            "  ->  fric2 exists? ", file.exists(fp2),
            " | fric3 exists? ", file.exists(fp3))
    next
  }
  
  # Read matrices
  m1 <- read_ratio_matrix(fp1)
  m2 <- read_ratio_matrix(fp2)
  m3 <- read_ratio_matrix(fp3)
  
  # Consistency check
  if (!identical(dim(m1), dim(m2)) || !identical(dim(m1), dim(m3)) ||
      !identical(rownames(m1), rownames(m2)) || !identical(rownames(m1), rownames(m3)) ||
      !identical(colnames(m1), colnames(m2)) || !identical(colnames(m1), colnames(m3))) {
    message("Dimension/name mismatch for ", stem, " — skipping.")
    next
  }
  
  years <- rownames(m1)
  elevs <- colnames(m1)
  
  # Save CSVs (now log10 fold-change via the updated function)
  out12_csv <- file.path(out_ratio_base, "fric1_2",
                         paste0(stem, "_fric1_2.csv"))
  save_difference_csv(m1, m2, years, elevs, out12_csv)
  
  out13_csv <- file.path(out_ratio_base, "fric1_3",
                         paste0(stem, "_fric1_3.csv"))
  save_difference_csv(m1, m3, years, elevs, out13_csv)
  
  out23_csv <- file.path(out_ratio_base, "fric2_3",
                         paste0(stem, "_fric2_3.csv"))
  save_difference_csv(m2, m3, years, elevs, out23_csv)
  
  # Generate heatmaps
  gl <- ids$glacier_num
  fl <- ids$flowline_num
  
  # Regular heatmaps
  out12_pdf <- file.path(out_heatmap_base, "fric1_2",
                         paste0("gl_", gl, "_", fl, "_hm_fric1_2.pdf"))
  plot_heatmap_pair(m1, m2, gl, fl, "fric1_2", out12_pdf)
  
  out13_pdf <- file.path(out_heatmap_base, "fric1_3",
                         paste0("gl_", gl, "_", fl, "_hm_fric1_3.pdf"))
  plot_heatmap_pair(m1, m3, gl, fl, "fric1_3", out13_pdf)
  
  out23_pdf <- file.path(out_heatmap_base, "fric2_3",
                         paste0("gl_", gl, "_", fl, "_hm_fric2_3.pdf"))
  plot_heatmap_pair(m2, m3, gl, fl, "fric2_3", out23_pdf)
  
  # Masked heatmaps
  out12_pdf <- file.path(out_heatmap_masked_base, "fric1_2",
                         paste0("gl_", gl, "_", fl, "_hmm_fric1_2.pdf"))
  plot_heatmap_pair_masked(m1, m2, gl, fl, "fric1_2", out12_pdf)
  
  out13_pdf <- file.path(out_heatmap_masked_base, "fric1_3",
                         paste0("gl_", gl, "_", fl, "_hmm_fric1_3.pdf"))
  plot_heatmap_pair_masked(m1, m3, gl, fl, "fric1_3", out13_pdf)
  
  out23_pdf <- file.path(out_heatmap_masked_base, "fric2_3",
                         paste0("gl_", gl, "_", fl, "_hmm_fric2_3.pdf"))
  plot_heatmap_pair_masked(m2, m3, gl, fl, "fric2_3", out23_pdf)
  
  cat("Processed: ", stem, " (Glacier ", gl, ", Flowline ", fl, ")\n", sep = "")
}

cat("Done.\n")