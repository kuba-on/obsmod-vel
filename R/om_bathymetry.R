########################################
# Correlate velocity & latitude with ratio bands (fric1/2/3)
# BASE R ONLY
########################################

# ---- Shared inputs (reuse your ratio files setup) ----
bed_csv <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Bathymetry Analysis/Mean Depth/bed_mean.csv"
vel_summary_csv <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Correlation Analysis/Mean V/a_obs_vel_mean.csv"

ratio_dir   <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Multiyear Mean"
ratio_files <- c(fric1 = file.path(ratio_dir, "mean_ratio_fric1.csv"),
                 fric2 = file.path(ratio_dir, "mean_ratio_fric2.csv"),
                 fric3 = file.path(ratio_dir, "mean_ratio_fric3.csv"))

# ---- Outputs: Velocity ----
vel_csv_dir   <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Correlation Analysis/Mean V/CSV"
vel_pdf_base  <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Correlation Analysis/Mean V/Graphs/PDF"
vel_png_base  <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Correlation Analysis/Mean V/Graphs/PNG"

# ---- Outputs: Latitude ----
lat_csv_dir   <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Correlation Analysis/Latitude/CSV"
lat_pdf_base  <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Correlation Analysis/Latitude/Graphs/PDF"
lat_png_base  <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Correlation Analysis/Latitude/Graphs/PNG"

suppressWarnings(dir.create(vel_csv_dir, recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(vel_pdf_base, recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(vel_png_base, recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(lat_csv_dir, recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(lat_pdf_base, recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(lat_png_base, recursive = TRUE, showWarnings = FALSE))

# ---- Read base tables ----
bed_df <- read.csv(bed_csv, check.names = FALSE, stringsAsFactors = FALSE)
vel_df <- read.csv(vel_summary_csv, check.names = FALSE, stringsAsFactors = FALSE)

key_cols <- c("glacier_ID", "glacier_name", "feature_ID", "flowline_ID")
missing_bed <- setdiff(key_cols, names(bed_df))
if (length(missing_bed)) stop(paste("Missing key columns in bed_mean.csv:", paste(missing_bed, collapse = ", ")))
missing_vel <- setdiff(key_cols, names(vel_df))
if (length(missing_vel)) stop(paste("Missing key columns in a_obs_vel_mean.csv:", paste(missing_vel, collapse = ", ")))

# Exclude glacier_ID 103 & 104
bed_df$glacier_ID <- suppressWarnings(as.integer(bed_df$glacier_ID))
vel_df$glacier_ID <- suppressWarnings(as.integer(vel_df$glacier_ID))
bed_df <- subset(bed_df, !(glacier_ID %in% c(103, 104)))
vel_df <- subset(vel_df, !(glacier_ID %in% c(103, 104)))

# Latitude column robustness (handle earlier typo lat_4236 vs lat_4326)
lat_col <- if ("lat_4326" %in% names(bed_df)) "lat_4326" else if ("lat_4236" %in% names(bed_df)) "lat_4236" else NA
if (is.na(lat_col)) stop("Neither 'lat_4326' nor 'lat_4236' found in bed_mean.csv.")

# ---- Helpers (same style as bathymetry code) ----
get_elev_cols <- function(df) {
  cols <- setdiff(names(df), key_cols)
  elev_nums <- suppressWarnings(as.numeric(cols))
  cols[is.finite(elev_nums)]
}

kendall_theil <- function(xlog, y) {
  ok <- is.finite(xlog) & is.finite(y)
  x <- xlog[ok]; y <- y[ok]
  n <- length(x)
  if (n < 2) return(list(a = NA_real_, b = NA_real_))
  idx_i <- rep.int(1:(n-1), (n-1):1)
  idx_j <- sequence((n-1):1, from = 2, by = 1) + rep.int(0:(n-2), (n-1):1)
  dx <- x[idx_j] - x[idx_i]
  dy <- y[idx_j] - y[idx_i]
  keep <- is.finite(dx) & dx != 0 & is.finite(dy)
  if (!any(keep)) return(list(a = NA_real_, b = NA_real_))
  slopes <- dy[keep] / dx[keep]
  b <- stats::median(slopes, na.rm = TRUE)
  a <- stats::median(y - b * x, na.rm = TRUE)
  list(a = a, b = b)
}

# Plot builder (scatter with fixed axes, double x-axis, 3-decimal rho/p)
build_plotter <- function(global_breaks, xticks, global_yrange, y_label) {
  function(ratio_vals, y_vals, main_label, col_pts = "#648FFF") {
    ok <- is.finite(ratio_vals) & ratio_vals > 0 & is.finite(y_vals)
    x <- ratio_vals[ok]; y <- y_vals[ok]
    rho <- NA_real_; pval <- NA_real_
    if (length(x) >= 3) {
      ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
      rho <- unname(ct$estimate); pval <- ct$p.value
    }
    lx <- log10(x)
    par(mar = c(6, 6, 3, 2))
    plot(NA, xlim = range(global_breaks), ylim = global_yrange,
         xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
    box()
    points(lx, y, pch = 16, col = col_pts)
    kt <- kendall_theil(lx, y)
    if (is.finite(kt$a) && is.finite(kt$b)) abline(a = kt$a, b = kt$b, lty = 2, col = col_pts)
    axis(2, las = 1)
    axis(1, at = xticks, labels = parse(text = paste0("10^", xticks)))
    mtext(y_label,                     side = 2, line = 3, cex = 1.1)
    mtext("velocity driver, R",        side = 1, line = 4.5, cex = 1.1)
    abline(v = 0, col = "black", lwd = 1)  # 10^0
    usr <- par("usr"); xpd_prev <- par("xpd"); par(xpd = TRUE)
    axis_bottom_y <- usr[3] - 0.15 * diff(usr[3:4])
    tick_length   <- 0.03 * diff(usr[3:4])
    x_min_tick <- min(xticks, na.rm = TRUE); x_max_tick <- max(xticks, na.rm = TRUE)
    segments(x0 = x_min_tick, y0 = axis_bottom_y, x1 = x_max_tick, y1 = axis_bottom_y, col = "black")
    for (tx in c(x_min_tick, 0, x_max_tick)) {
      segments(x0 = tx, y0 = axis_bottom_y, x1 = tx, y1 = axis_bottom_y + tick_length, col = "black")
    }
    mid_left  <- (x_min_tick + 0) / 2
    mid_right <- (x_max_tick + 0) / 2
    text(x = mid_left,  y = axis_bottom_y - 0.035 * diff(usr[3:4]), labels = "friction", cex = 1)
    text(x = mid_right, y = axis_bottom_y - 0.035 * diff(usr[3:4]), labels = "front",    cex = 1)
    par(xpd = xpd_prev)
    dx <- diff(par("usr")[1:2]); dy <- diff(par("usr")[3:4])
    tlx <- par("usr")[1] + 0.02 * dx
    tly <- par("usr")[4] - 0.05 * dy
    text(tlx, tly,             labels = paste0("rho = ", if (is.finite(rho)) sprintf("%.3f", rho) else "NA"), adj = c(0, 1))
    text(tlx, tly - 0.06 * dy, labels = paste0("p-value = ", if (is.finite(pval)) sprintf("%.3f", pval) else "NA"), adj = c(0, 1))
    mtext(main_label, side = 3, line = 0.5, cex = 1.4)
    invisible(list(rho = rho, p = pval))
  }
}

# Compute global X (ratio) breaks and return per-law merged frames with elevation columns
prepare_merged_and_axes <- function(y_df, ratio_files) {
  merged_list <- list()
  all_ratio_vals <- numeric(0)
  for (fric in names(ratio_files)) {
    rdf <- read.csv(ratio_files[[fric]], check.names = FALSE, stringsAsFactors = FALSE)
    missing_keys2 <- setdiff(key_cols, names(rdf))
    if (length(missing_keys2)) { warning(paste("Skipping", fric, "- missing key columns:", paste(missing_keys2, collapse = ", "))); next }
    rdf$glacier_ID <- suppressWarnings(as.integer(rdf$glacier_ID))
    rdf <- subset(rdf, !(glacier_ID %in% c(103, 104)))
    merged <- merge(y_df, rdf, by = key_cols, all = FALSE)
    if (!nrow(merged)) { warning(paste("No matches after merge for", fric)); next }
    elev_cols <- get_elev_cols(rdf)
    merged_list[[fric]] <- list(df = merged, elev_cols = elev_cols)
    if (length(elev_cols)) {
      v <- unlist(merged[, elev_cols, drop = FALSE], use.names = FALSE)
      v <- as.numeric(v); v <- v[is.finite(v) & v > 0]
      all_ratio_vals <- c(all_ratio_vals, v)
    }
  }
  if (!length(all_ratio_vals)) stop("No positive ratio values after filtering & merging.")
  lx_all <- log10(all_ratio_vals)
  global_xrange <- range(lx_all, na.rm = TRUE)
  break_spacing <- 0.1
  raw_breaks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = break_spacing)
  if (!any(abs(raw_breaks) < 1e-8)) raw_breaks <- sort(unique(c(raw_breaks, 0)))
  global_breaks <- raw_breaks
  xticks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = 1)
  list(merged_list = merged_list, global_breaks = global_breaks, xticks = xticks)
}

# =========================
# 1) VELOCITY vs RATIO
# =========================
# Build per-law merges (on vel_df) and global fixed axes from ratio values
prep_vel <- prepare_merged_and_axes(vel_df, ratio_files)

# Fixed Y across all velocity plots
all_vel_vals <- numeric(0)
for (fric in names(prep_vel$merged_list)) {
  mdf <- prep_vel$merged_list[[fric]]$df
  # collect ALL vel means present (per-elev columns start with vel_mean_elev_)
  vel_cols <- grep("^vel_mean_elev_", names(mdf), value = TRUE)
  if (length(vel_cols)) {
    v <- unlist(mdf[, vel_cols, drop = FALSE], use.names = FALSE)
    v <- as.numeric(v); v <- v[is.finite(v)]
    all_vel_vals <- c(all_vel_vals, v)
  }
  # include 100–500 column if present
  if ("vel_mean_elev_100-500" %in% names(mdf)) {
    vv <- as.numeric(mdf$vel_mean_elev_100-500)
    vv <- vv[is.finite(vv)]
    all_vel_vals <- c(all_vel_vals, vv)
  }
}
if (!length(all_vel_vals)) stop("No velocity values found for plotting.")
yr <- range(all_vel_vals, na.rm = TRUE)
pad_y <- 0.05 * diff(yr); if (!is.finite(pad_y)) pad_y <- 1
vel_global_yrange <- c(yr[1] - pad_y, yr[2] + pad_y)

plot_vel <- build_plotter(prep_vel$global_breaks, prep_vel$xticks, vel_global_yrange,
                          y_label = "mean observed velocity")

for (fric in names(prep_vel$merged_list)) {
  merged    <- prep_vel$merged_list[[fric]]$df
  elev_cols <- prep_vel$merged_list[[fric]]$elev_cols
  if (!length(elev_cols)) { warning(paste("No elevation columns for", fric)); next }
  
  elev_nums <- as.numeric(elev_cols); ord <- order(elev_nums)
  elev_cols <- elev_cols[ord]; elev_nums <- elev_nums[ord]
  
  # Per-law output folders
  law_pdf_dir <- file.path(vel_pdf_base, fric)
  law_png_dir <- file.path(vel_png_base, fric)
  suppressWarnings(dir.create(law_pdf_dir, recursive = TRUE, showWarnings = FALSE))
  suppressWarnings(dir.create(law_png_dir, recursive = TRUE, showWarnings = FALSE))
  
  # Results container (add an extra row for "100-500")
  res_vel <- data.frame(elev = as.character(elev_nums),
                        Spearman_rho = NA_real_,
                        p_value = NA_real_,
                        stringsAsFactors = FALSE)
  
  # Per-elevation correlations and plots
  for (i in seq_along(elev_cols)) {
    col_elev <- elev_cols[i]
    elnum    <- elev_nums[i]
    x_ratio  <- as.numeric(merged[[col_elev]])
    y_velcol <- paste0("vel_mean_elev_", elnum)
    if (!(y_velcol %in% names(merged))) {
      warning(sprintf("Missing %s for %s; skipping elev %s", y_velcol, fric, elnum))
      next
    }
    y_vel <- as.numeric(merged[[y_velcol]])
    
    ok <- is.finite(x_ratio) & x_ratio > 0 & is.finite(y_vel)
    x_use <- x_ratio[ok]; y_use <- y_vel[ok]
    
    rho <- NA_real_; pval <- NA_real_
    if (length(x_use) >= 3) {
      ct <- suppressWarnings(cor.test(x_use, y_use, method = "spearman", exact = FALSE))
      rho <- unname(ct$estimate); pval <- ct$p.value
    }
    res_vel$Spearman_rho[res_vel$elev == as.character(elnum)] <- rho
    res_vel$p_value[res_vel$elev == as.character(elnum)]      <- pval
    
    # Save plots
    pdf(file.path(law_pdf_dir, sprintf("vel_ratio_corr_%d_%s.pdf", elnum, fric)), width = 9, height = 6)
    plot_vel(x_use, y_use, main_label = sprintf("Elevation %d, %s", elnum, fric))
    dev.off()
    
    png(file.path(law_png_dir, sprintf("vel_ratio_corr_%d_%s.png", elnum, fric)),
        width = 9, height = 6, units = "in", res = 300)
    plot_vel(x_use, y_use, main_label = sprintf("Elevation %d, %s", elnum, fric))
    dev.off()
  }
  
  # Extra: 100–500 mean correlation (ratio mean vs vel_mean_elev_100-500)
  use_elevs <- elev_nums[elev_nums >= 100 & elev_nums <= 500]
  use_cols  <- as.character(use_elevs)
  if (length(use_cols) && ("vel_mean_elev_100-500" %in% names(merged))) {
    ratio_mean_100_500 <- rowMeans(merged[, use_cols, drop = FALSE], na.rm = TRUE)
    x_use <- as.numeric(ratio_mean_100_500)
    y_use <- as.numeric(merged$vel_mean_elev_100-500)
    
    ok <- is.finite(x_use) & x_use > 0 & is.finite(y_use)
    x_use <- x_use[ok]; y_use <- y_use[ok]
    
    rho <- NA_real_; pval <- NA_real_
    if (length(x_use) >= 3) {
      ct <- suppressWarnings(cor.test(x_use, y_use, method = "spearman", exact = FALSE))
      rho <- unname(ct$estimate); pval <- ct$p.value
    }
    
    res_vel <- rbind(res_vel,
                     data.frame(elev = "100-500", Spearman_rho = rho, p_value = pval, stringsAsFactors = FALSE))
    
    pdf(file.path(law_pdf_dir, sprintf("vel_ratio_corr_100-500_%s.pdf", fric)), width = 9, height = 6)
    plot_vel(x_use, y_use, main_label = sprintf("Elevation 100–500 mean, %s", fric))
    dev.off()
    
    png(file.path(law_png_dir, sprintf("vel_ratio_corr_100-500_%s.png", fric)),
        width = 9, height = 6, units = "in", res = 300)
    plot_vel(x_use, y_use, main_label = sprintf("Elevation 100–500 mean, %s", fric))
    dev.off()
  }
  
  # Save CSV (per law)
  write.csv(res_vel, file.path(vel_csv_dir, sprintf("vel_ratio_corr_%s.csv", fric)), row.names = FALSE)
}

# =========================
# 2) LATITUDE vs RATIO
# =========================
# Build per-law merges (on bed_df for latitude) and global fixed axes from ratio values
prep_lat <- prepare_merged_and_axes(bed_df, ratio_files)

# Fixed Y across all latitude plots
all_lat_vals <- numeric(0)
for (fric in names(prep_lat$merged_list)) {
  mdf <- prep_lat$merged_list[[fric]]$df
  if (nrow(mdf)) {
    vv <- as.numeric(mdf[[lat_col]])
    vv <- vv[is.finite(vv)]
    all_lat_vals <- c(all_lat_vals, vv)
  }
}
if (!length(all_lat_vals)) stop("No latitude values found for plotting.")
yr <- range(all_lat_vals, na.rm = TRUE)
pad_y <- 0.05 * diff(yr); if (!is.finite(pad_y)) pad_y <- 1
lat_global_yrange <- c(yr[1] - pad_y, yr[2] + pad_y)

plot_lat <- build_plotter(prep_lat$global_breaks, prep_lat$xticks, lat_global_yrange,
                          y_label = "latitude (°)")

for (fric in names(prep_lat$merged_list)) {
  merged    <- prep_lat$merged_list[[fric]]$df
  elev_cols <- prep_lat$merged_list[[fric]]$elev_cols
  if (!length(elev_cols)) { warning(paste("No elevation columns for", fric)); next }
  
  elev_nums <- as.numeric(elev_cols); ord <- order(elev_nums)
  elev_cols <- elev_cols[ord]; elev_nums <- elev_nums[ord]
  
  law_pdf_dir <- file.path(lat_pdf_base, fric)
  law_png_dir <- file.path(lat_png_base, fric)
  suppressWarnings(dir.create(law_pdf_dir, recursive = TRUE, showWarnings = FALSE))
  suppressWarnings(dir.create(law_png_dir, recursive = TRUE, showWarnings = FALSE))
  
  res_lat <- data.frame(elev = as.character(elev_nums),
                        Spearman_rho = NA_real_,
                        p_value = NA_real_,
                        stringsAsFactors = FALSE)
  
  y_lat <- as.numeric(merged[[lat_col]])
  
  for (i in seq_along(elev_cols)) {
    col_elev <- elev_cols[i]
    elnum    <- elev_nums[i]
    x_ratio  <- as.numeric(merged[[col_elev]])
    
    ok <- is.finite(x_ratio) & x_ratio > 0 & is.finite(y_lat)
    x_use <- x_ratio[ok]; y_use <- y_lat[ok]
    
    rho <- NA_real_; pval <- NA_real_
    if (length(x_use) >= 3) {
      ct <- suppressWarnings(cor.test(x_use, y_use, method = "spearman", exact = FALSE))
      rho <- unname(ct$estimate); pval <- ct$p.value
    }
    res_lat$Spearman_rho[res_lat$elev == as.character(elnum)] <- rho
    res_lat$p_value[res_lat$elev == as.character(elnum)]      <- pval
    
    pdf(file.path(law_pdf_dir, sprintf("lat_ratio_corr_%d_%s.pdf", elnum, fric)), width = 9, height = 6)
    plot_lat(x_use, y_use, main_label = sprintf("Elevation %d, %s", elnum, fric))
    dev.off()
    
    png(file.path(law_png_dir, sprintf("lat_ratio_corr_%d_%s.png", elnum, fric)),
        width = 9, height = 6, units = "in", res = 300)
    plot_lat(x_use, y_use, main_label = sprintf("Elevation %d, %s", elnum, fric))
    dev.off()
  }
  
  # 100–500 mean correlation (ratio mean vs latitude)
  use_elevs <- elev_nums[elev_nums >= 100 & elev_nums <= 500]
  use_cols  <- as.character(use_elevs)
  if (length(use_cols)) {
    ratio_mean_100_500 <- rowMeans(merged[, use_cols, drop = FALSE], na.rm = TRUE)
    x_use <- as.numeric(ratio_mean_100_500)
    y_use <- as.numeric(merged[[lat_col]])
    
    ok <- is.finite(x_use) & x_use > 0 & is.finite(y_use)
    x_use <- x_use[ok]; y_use <- y_use[ok]
    
    rho <- NA_real_; pval <- NA_real_
    if (length(x_use) >= 3) {
      ct <- suppressWarnings(cor.test(x_use, y_use, method = "spearman", exact = FALSE))
      rho <- unname(ct$estimate); pval <- ct$p.value
    }
    
    res_lat <- rbind(res_lat,
                     data.frame(elev = "100-500", Spearman_rho = rho, p_value = pval, stringsAsFactors = FALSE))
    
    pdf(file.path(law_pdf_dir, sprintf("lat_ratio_corr_100-500_%s.pdf", fric)), width = 9, height = 6)
    plot_lat(x_use, y_use, main_label = sprintf("Elevation 100–500 mean, %s", fric))
    dev.off()
    
    png(file.path(law_png_dir, sprintf("lat_ratio_corr_100-500_%s.png", fric)),
        width = 9, height = 6, units = "in", res = 300)
    plot_lat(x_use, y_use, main_label = sprintf("Elevation 100–500 mean, %s", fric))
    dev.off()
  }
  
  write.csv(res_lat, file.path(lat_csv_dir, sprintf("lat_ratio_corr_%s.csv", fric)), row.names = FALSE)
}