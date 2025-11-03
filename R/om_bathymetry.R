########################################
# Correlate bed_mean with ratio bands (fric1/2/3) and plot results
# BASE R ONLY — fixed p_value column + PDF & PNG outputs
########################################

# ---- Paths ----
bed_csv <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Bathymetry Analysis/Mean Depth/bed_mean.csv"

ratio_dir      <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Multiyear Mean"
ratio_files    <- c(fric1 = file.path(ratio_dir, "mean_ratio_fric1.csv"),
                    fric2 = file.path(ratio_dir, "mean_ratio_fric2.csv"),
                    fric3 = file.path(ratio_dir, "mean_ratio_fric3.csv"))

base_corr_dir  <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Bathymetry Analysis/Correlation"
csv_out_dir    <- file.path(base_corr_dir, "CSV")
graphs_base    <- file.path(base_corr_dir, "Graphs")
graphs_pdf_dir <- file.path(graphs_base, "PDF")
graphs_png_dir <- file.path(graphs_base, "PNG")

suppressWarnings(dir.create(csv_out_dir,    recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(graphs_pdf_dir, recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(graphs_png_dir, recursive = TRUE, showWarnings = FALSE))

# ---- Read data ----
bed_df <- read.csv(bed_csv, check.names = FALSE, stringsAsFactors = FALSE)

# keys (first 4 columns identical across files)
key_cols <- c("glacier_ID", "glacier_name", "feature_ID", "flowline_ID")
missing_keys <- setdiff(key_cols, names(bed_df))
if (length(missing_keys) > 0) stop(paste("Missing key columns in bed_mean.csv:", paste(missing_keys, collapse = ", ")))

# EXCLUDE glacier_ID 103 and 104 (coerce safely)
bed_df$glacier_ID <- suppressWarnings(as.integer(bed_df$glacier_ID))
bed_df <- subset(bed_df, !(feature_ID %in% c(103, 104)))

# ---- Helpers ----
get_elev_cols <- function(df) {
  cols <- setdiff(names(df), key_cols)
  elev_nums <- suppressWarnings(as.numeric(cols))
  cols[is.finite(elev_nums)]
}

# Kendall–Theil (Theil–Sen) robust line on log10(x)
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

# ---- Pre-merge to compute fixed axes ----
merged_list <- list()
all_ratio_vals <- numeric(0)
all_bed_vals   <- numeric(0)

for (fric in names(ratio_files)) {
  ratio_df <- read.csv(ratio_files[[fric]], check.names = FALSE, stringsAsFactors = FALSE)
  
  missing_keys2 <- setdiff(key_cols, names(ratio_df))
  if (length(missing_keys2) > 0) { warning(paste("Skipping", fric, "- missing key columns:", paste(missing_keys2, collapse = ", "))); next }
  
  ratio_df$glacier_ID <- suppressWarnings(as.integer(ratio_df$glacier_ID))
  ratio_df <- subset(ratio_df, !(glacier_ID %in% c(103, 104)))
  
  merged <- merge(bed_df, ratio_df, by = key_cols, all = FALSE)
  if (!nrow(merged)) { warning(paste("No matching rows after merge for", fric)); next }
  
  elev_cols <- get_elev_cols(ratio_df)
  merged_list[[fric]] <- list(df = merged, elev_cols = elev_cols)
  
  if (length(elev_cols)) {
    v <- unlist(merged[ , elev_cols, drop = FALSE], use.names = FALSE)
    v <- as.numeric(v); v <- v[is.finite(v) & v > 0]
    all_ratio_vals <- c(all_ratio_vals, v)
  }
  yv <- as.numeric(merged$bed_mean)
  yv <- yv[is.finite(yv)]
  all_bed_vals <- c(all_bed_vals, yv)
}

if (!length(all_ratio_vals)) stop("No positive ratio values found after filtering & merging.")
if (!length(all_bed_vals))   stop("No bed_mean values found after filtering & merging.")

# ---- Fixed axes (X & Y) ----
lx_all <- log10(all_ratio_vals)
global_xrange <- range(lx_all, na.rm = TRUE)

break_spacing <- 0.1
raw_breaks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = break_spacing)
if (!any(abs(raw_breaks) < 1e-8)) raw_breaks <- sort(unique(c(raw_breaks, 0)))
global_breaks <- raw_breaks
xticks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = 1)

global_yrange <- range(all_bed_vals, na.rm = TRUE)
pad_y <- 0.05 * diff(global_yrange); if (!is.finite(pad_y)) pad_y <- 1
global_yrange <- c(global_yrange[1] - pad_y, global_yrange[2] + pad_y)

# ---- Plot function ----
plot_corr_scatter <- function(ratio_vals, bed_vals, main_label,
                              col_pts = "#648FFF",
                              y_lim = global_yrange) {
  ok <- is.finite(ratio_vals) & ratio_vals > 0 & is.finite(bed_vals)
  x <- ratio_vals[ok]; y <- bed_vals[ok]
  
  rho <- NA_real_; pval <- NA_real_
  if (length(x) >= 3) {
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    rho <- unname(ct$estimate); pval <- ct$p.value
  }
  
  lx <- log10(x)
  
  par(mar = c(6, 6, 3, 2))  # larger margins
  plot(NA, xlim = range(global_breaks), ylim = y_lim,
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  box()
  
  points(lx, y, pch = 16, col = col_pts)
  
  kt <- kendall_theil(lx, y)
  if (is.finite(kt$a) && is.finite(kt$b)) {
    abline(a = kt$a, b = kt$b, lty = 2, col = col_pts)
  }
  
  axis(2, las = 1)
  axis(1, at = xticks, labels = parse(text = paste0("10^", xticks)))
  mtext("bed elevation at the glacier front [in m a.s.l.]", side = 2, line = 3, cex = 1.1)
  mtext("velocity driver, R",                              side = 1, line = 4.5, cex = 1.1)
  
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

# ---- Main loop per friction law ----
for (fric in names(merged_list)) {
  merged    <- merged_list[[fric]]$df
  elev_cols <- merged_list[[fric]]$elev_cols
  if (!length(elev_cols)) { warning(paste("No elevation columns found for", fric)); next }
  
  # Sort elevation columns numerically
  elev_nums <- as.numeric(elev_cols); ord <- order(elev_nums)
  elev_cols <- elev_cols[ord]; elev_nums <- elev_nums[ord]
  
  # Per-law output folders
  law_pdf_dir <- file.path(graphs_pdf_dir, fric)
  law_png_dir <- file.path(graphs_png_dir, fric)
  suppressWarnings(dir.create(law_pdf_dir, recursive = TRUE, showWarnings = FALSE))
  suppressWarnings(dir.create(law_png_dir, recursive = TRUE, showWarnings = FALSE))
  
  # Results container — NOTE: p_value (underscore) to avoid name mangling
  res_mat <- data.frame(elev = elev_nums,
                        Spearman_rho = NA_real_,
                        p_value = NA_real_,
                        stringsAsFactors = FALSE)
  
  # Loop over elevations
  for (i in seq_along(elev_cols)) {
    col_elev <- elev_cols[i]
    elnum    <- elev_nums[i]
    
    x_ratio  <- as.numeric(merged[[col_elev]])
    y_bed    <- as.numeric(merged[["bed_mean"]])
    
    ok <- is.finite(x_ratio) & x_ratio > 0 & is.finite(y_bed)
    x_use <- x_ratio[ok]; y_use <- y_bed[ok]
    
    rho <- NA_real_; pval <- NA_real_
    if (length(x_use) >= 3) {
      ct <- suppressWarnings(cor.test(x_use, y_use, method = "spearman", exact = FALSE))
      rho <- unname(ct$estimate); pval <- ct$p.value
    }
    res_mat$Spearman_rho[i] <- rho
    res_mat$p_value[i]      <- pval
    
    # --- Save PDF ---
    pdf(file.path(law_pdf_dir, sprintf("bed_ratio_corr_%d_%s.pdf", elnum, fric)), width = 9, height = 6)
    plot_corr_scatter(x_use, y_use, main_label = sprintf("Elevation %d, %s", elnum, fric))
    dev.off()
    
    # --- Save PNG ---
    png(file.path(law_png_dir, sprintf("bed_ratio_corr_%d_%s.png", elnum, fric)),
        width = 9, height = 6, units = "in", res = 300)
    plot_corr_scatter(x_use, y_use, main_label = sprintf("Elevation %d, %s", elnum, fric))
    dev.off()
  }
  
  # CSV (per law) — only Spearman_rho and p_value columns (no extra p.value)
  out_csv <- file.path(csv_out_dir, sprintf("bed_ratio_corr_%s.csv", fric))
  write.csv(res_mat, out_csv, row.names = FALSE)
}
