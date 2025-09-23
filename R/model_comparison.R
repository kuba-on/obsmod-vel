########################################

# The code below compares friction–frontal ratio matrices that use different friction laws (fric1, fric2, fric3) in a loop.
# It computes the raw differences between each pair of friction laws and generates heatmaps of the log10 fold-change.

########################################


options(stringsAsFactors = FALSE)

# Directories
dir_fric1 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric1/Ratios"
dir_fric2 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric2/Ratios"
dir_fric3 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric3/Ratios"

out_ratio_base   <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Ratios"
out_heatmap_base <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Heatmaps"

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

# Function to save raw difference CSV
save_difference_csv <- function(matA, matB, years, elevs, out_path) {
  stopifnot(identical(rownames(matA), years), identical(colnames(matA), elevs))
  stopifnot(identical(rownames(matB), years), identical(colnames(matB), elevs))
  diff_mat <- matA - matB
  write.csv(diff_mat, out_path, row.names = TRUE)
}

# Heatmap plotting function
plot_heatmap_pair <- function(matA, matB, glacier_num, flowline_num, pair_label, out_pdf) {
  ratio_matrix <- matA / pmax(matB, .Machine$double.eps)  # avoid /0
  years <- rownames(ratio_matrix)
  elevations <- colnames(ratio_matrix)
  
  # Set up PDF
  pdf(out_pdf, width = 9, height = 4)
  par(mar = c(0.2, 0.2, 0.2, 2), oma = c(4.5, 4, 1.5, 0.5))
  
  log_ratios <- log10(pmax(ratio_matrix, 10^-1))
  breaks <- seq(-1, 1, length.out = 1000)
  color_palette <- colorRampPalette(c("#364B9A", "#EAECCC", "#A50026"))(999)
  
  rev_years <- rev(years)
  log_ratios_reordered <- log_ratios[rev(seq_len(nrow(log_ratios))), ]
  
  layout(matrix(c(1,2), nrow = 1), widths = c(4, 0.5))
  
  image(1:length(elevations), 1:length(rev_years), t(log_ratios_reordered),
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
  
  par(mar = c(0.2, 0, 0.2, 3))
  image(1, seq(-1, 1, length.out = 1000),
        matrix(seq(-1, 1, length.out = 1000), nrow = 1),
        col = color_palette, breaks = breaks, axes = FALSE,
        xlab = "", ylab = "")
  
  log_ticks <- seq(-1, 1, by = 1)
  log_labels <- parse(text = paste0("10^", log_ticks))
  axis(4, at = log_ticks, labels = log_labels, las = 1, cex.axis = 0.8, line = 0.1)
  mtext("velocity driver", side = 4, line = 2.25, cex = 0.8)
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
  
  # Save raw difference CSVs
  # fric1 - fric2
  out12_csv <- file.path(out_ratio_base, "fric1_2",
                         paste0(stem, "_fric1_2.csv"))
  save_difference_csv(m1, m2, years, elevs, out12_csv)
  
  # fric1 - fric3
  out13_csv <- file.path(out_ratio_base, "fric1_3",
                         paste0(stem, "_fric1_3.csv"))
  save_difference_csv(m1, m3, years, elevs, out13_csv)
  
  # fric2 - fric3
  out23_csv <- file.path(out_ratio_base, "fric2_3",
                         paste0(stem, "_fric2_3.csv"))
  save_difference_csv(m2, m3, years, elevs, out23_csv)
  
  # Generate heatmaps
  gl <- ids$glacier_num
  fl <- ids$flowline_num
  
  out12_pdf <- file.path(out_heatmap_base, "fric1_2",
                         paste0("gl_", gl, "_", fl, "_hm_fric1_2.pdf"))
  plot_heatmap_pair(m1, m2, gl, fl, "fric1_2", out12_pdf)
  
  out13_pdf <- file.path(out_heatmap_base, "fric1_3",
                         paste0("gl_", gl, "_", fl, "_hm_fric1_3.pdf"))
  plot_heatmap_pair(m1, m3, gl, fl, "fric1_3", out13_pdf)
  
  out23_pdf <- file.path(out_heatmap_base, "fric2_3",
                         paste0("gl_", gl, "_", fl, "_hm_fric2_3.pdf"))
  plot_heatmap_pair(m2, m3, gl, fl, "fric2_3", out23_pdf)
  
  cat("Processed: ", stem, " (Glacier ", gl, ", Flowline ", fl, ")\n", sep = "")
}

cat("Done.\n")

cat("Done.\n")