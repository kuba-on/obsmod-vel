########################################

# The code below compares distributions of ratios of seasonality drivers obtained from different friction laws. The processing includes:
# - Reading CSV files containing all ratio data for different friction laws.
# - Plotting overlaid histograms with transparency, grouped by year, elevation, flowline, and globally with all data combined.
# - Plots contain vertical lines at the mean of each distribution.

########################################

# Directories
dir_fric1 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric1/Ratios"
dir_fric2 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric2/Ratios"
dir_fric3 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric3/Ratios"

out_global <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Histograms/Ratios/Global"
out_years  <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Histograms/Ratios/Years"
out_elevs  <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Histograms/Ratios/Elevations"
out_flows  <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Histograms/Ratios/Flowlines"

suppressWarnings(dir.create(out_global, recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(out_years,  recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(out_elevs,  recursive = TRUE, showWarnings = FALSE))
suppressWarnings(dir.create(out_flows,  recursive = TRUE, showWarnings = FALSE))

# Glacier name lookup
csv_lookup <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/New Points v3/Input/Box Coordinates/box_sp_all_v3.csv"
lookup_df <- tryCatch(read.csv(csv_lookup), error = function(e) NULL)

# Read all ratio CSVs for a friction law folder.
# Return a named list of numeric matrices; names like "gl_2_1"
read_ratio_folder <- function(folder_path, friction_tag) {
  patt <- paste0("^gl_\\d+_\\d+_ratio_", friction_tag, "\\.csv$")
  files <- list.files(folder_path, pattern = patt, full.names = TRUE)
  out <- list()
  for (f in files) {
    df <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) next
    
    # If row names were written as first column (from write.csv row.names=TRUE), move them into rownames and drop the first column
    if (ncol(df) > 1 && (is.character(df[[1]]) || is.numeric(df[[1]]))) {
      maybe_year <- suppressWarnings(as.numeric(df[[1]]))
      if (sum(!is.na(maybe_year)) > 0) {
        rownames(df) <- as.character(df[[1]])
        df <- df[, -1, drop = FALSE]
      }
    }
    
    # Coerce to numeric matrix, preserving column names (elevations)
    colnames(df) <- trimws(colnames(df))
    suppressWarnings(colnames(df) <- as.character(as.numeric(colnames(df))))
    m <- suppressWarnings(apply(df, 2, function(x) as.numeric(as.character(x))))
    if (is.null(dim(m))) {
      m <- matrix(m, ncol = 1)
      colnames(m) <- colnames(df)
      rownames(m) <- rownames(df)
    } else {
      rownames(m) <- rownames(df)
      colnames(m) <- colnames(df)
    }
    
    base <- basename(f)
    id <- sub(paste0("_ratio_", friction_tag, "\\.csv$"), "", base)
    out[[id]] <- m
  }
  out
}

# Flatten list of matrices to a single numeric vector (drop NA & non-positive)
flatten_vec <- function(mat_list) {
  if (length(mat_list) == 0) return(numeric(0))
  v <- unlist(mat_list, use.names = FALSE)
  v <- as.numeric(v)
  v <- v[is.finite(v) & !is.na(v) & v > 0]
  v
}

# Collect values by year
collect_by_year <- function(mat_list) {
  years <- unique(as.numeric(unlist(lapply(mat_list, function(m) rownames(m)))))
  years <- years[is.finite(years)]
  out <- setNames(vector("list", length(years)), as.character(sort(years)))
  for (yr in names(out)) {
    vals <- unlist(lapply(mat_list, function(m) {
      if (!is.null(rownames(m)) && yr %in% rownames(m)) as.numeric(m[yr, , drop = TRUE])
    }))
    vals <- vals[is.finite(vals) & vals > 0]
    out[[yr]] <- vals
  }
  out
}

# Collect values by elevation
collect_by_elev <- function(mat_list) {
  elevs <- unique(as.numeric(unlist(lapply(mat_list, function(m) colnames(m)))))
  elevs <- elevs[is.finite(elevs)]
  out <- setNames(vector("list", length(elevs)), as.character(sort(elevs)))
  for (el in names(out)) {
    vals <- unlist(lapply(mat_list, function(m) {
      if (!is.null(colnames(m)) && el %in% colnames(m)) as.numeric(m[, el, drop = TRUE])
    }))
    vals <- vals[is.finite(vals) & vals > 0]
    out[[el]] <- vals
  }
  out
}

# Glacier name lookup helper
get_glacier_name <- function(glacier_num) {
  if (is.null(lookup_df)) return("Unknown Glacier")
  nm <- lookup_df$glacier_name[lookup_df$feature_ID == as.integer(glacier_num)]
  if (length(nm) == 0) "Unknown Glacier" else as.character(nm[1])
}

#  Read all three friction sets
fric1_list <- read_ratio_folder(dir_fric1, "fric1")
fric2_list <- read_ratio_folder(dir_fric2, "fric2")
fric3_list <- read_ratio_folder(dir_fric3, "fric3")

# Define global axis range and breaks
x1_all <- flatten_vec(fric1_list)
x2_all <- flatten_vec(fric2_list)
x3_all <- flatten_vec(fric3_list)

x_all_global <- c(x1_all, x2_all, x3_all)
x_all_global <- x_all_global[is.finite(x_all_global) & x_all_global > 0]
lx_all_global <- log10(x_all_global)

global_xrange <- range(lx_all_global, na.rm = TRUE)
break_spacing <- 0.1
raw_breaks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = break_spacing)
if (!any(abs(raw_breaks) < 1e-8)) raw_breaks <- sort(unique(c(raw_breaks, 0)))
global_breaks <- raw_breaks

# Plot overlay function
plot_overlay_hist <- function(x1, x3, x2,
                              breaks = global_breaks,
                              col1 = "#648FFF", col2 = "#DC267F", col3 = "#FFB000",
                              alpha_fill = 0.35, alpha_border = 0.9,
                              main_label = "",
                              legend_pos = "topright") {
  allx <- c(x1, x2, x3)
  allx <- allx[is.finite(allx) & allx > 0]
  if (length(allx) == 0) {
    plot.new(); box(); mtext("No data", side = 3, line = 0.5, cex = 1)
    return(invisible(NULL))
  }
  lx1 <- log10(x1[x1 > 0])
  lx2 <- log10(x2[x2 > 0])
  lx3 <- log10(x3[x3 > 0])
  
  h1 <- hist(lx1, breaks = breaks, plot = FALSE)
  h2 <- hist(lx2, breaks = breaks, plot = FALSE)
  h3 <- hist(lx3, breaks = breaks, plot = FALSE)
  ylim_max <- max(c(h1$counts, h2$counts, h3$counts, 1))
  
  plot(NA, xlim = range(breaks), ylim = c(0, ylim_max),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  box()
  
  rect(xleft = h1$breaks[-length(h1$breaks)], ybottom = 0,
       xright = h1$breaks[-1], ytop = h1$counts,
       col = adjustcolor(col1, alpha.f = alpha_fill),
       border = adjustcolor(col1, alpha.f = alpha_border))
  rect(xleft = h2$breaks[-length(h2$breaks)], ybottom = 0,
       xright = h2$breaks[-1], ytop = h2$counts,
       col = adjustcolor(col2, alpha.f = alpha_fill),
       border = adjustcolor(col2, alpha.f = alpha_border))
  rect(xleft = h3$breaks[-length(h3$breaks)], ybottom = 0,
       xright = h3$breaks[-1], ytop = h3$counts,
       col = adjustcolor(col3, alpha.f = alpha_fill),
       border = adjustcolor(col3, alpha.f = alpha_border))
  
  axis(2, las = 1)
  xticks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = 1)
  axis(1, at = xticks, labels = parse(text = paste0("10^", xticks)))
  mtext("frequency", side = 2, line = 3, cex = 1.1)
  mtext("velocity driver, R", side = 1, line = 4.5, cex = 1.1)
  
  # --- Secondary inverted axis below main x-axis (aligned with tick range) ---
  usr <- par("usr")
  xpd_prev <- par("xpd")
  par(xpd = TRUE)
  
  axis_bottom_y <- usr[3] - 0.15 * diff(usr[3:4])  # vertical offset below x-axis
  tick_length   <- 0.03 * diff(usr[3:4])
  
  # Use tick limits instead of full x-range
  x_min_tick <- min(xticks, na.rm = TRUE)
  x_max_tick <- max(xticks, na.rm = TRUE)
  
  # Draw secondary axis line exactly between the first and last main ticks
  segments(x0 = x_min_tick, y0 = axis_bottom_y,
           x1 = x_max_tick, y1 = axis_bottom_y, col = "black")
  
  # Upward ticks at left end, 10^0, and right end
  ticks_x <- c(x_min_tick, 0, x_max_tick)
  for (tx in ticks_x) {
    segments(x0 = tx, y0 = axis_bottom_y,
             x1 = tx, y1 = axis_bottom_y + tick_length, col = "black")
  }
  
  # Compute midpoints between 10^0 and ends for labels
  mid_left  <- (x_min_tick + 0) / 2
  mid_right <- (x_max_tick + 0) / 2
  
  # Labels: "front" (left side) and "friction" (right side)
  text(x = mid_left,  y = axis_bottom_y - 0.035 * diff(usr[3:4]),
       labels = "friction", cex = 1, col = "black", font = 1)
  text(x = mid_right, y = axis_bottom_y - 0.035 * diff(usr[3:4]),
       labels = "front", cex = 1, col = "black", font = 1)
  
  par(xpd = xpd_prev)
  
  m1 <- if (length(x1)) mean(x1, na.rm = TRUE) else NA
  m2 <- if (length(x2)) mean(x2, na.rm = TRUE) else NA
  m3 <- if (length(x3)) mean(x3, na.rm = TRUE) else NA
  if (is.finite(m1) && m1 > 0) abline(v = log10(m1), col = col1, lty = 3, lwd = 2)
  if (is.finite(m2) && m2 > 0) abline(v = log10(m2), col = col2, lty = 3, lwd = 2)
  if (is.finite(m3) && m3 > 0) abline(v = log10(m3), col = col3, lty = 3, lwd = 2)
  
  legend(legend_pos,
         legend = c("Budd, m=1", "Budd, m=6", "Schoof, m=3"),
         fill = c(adjustcolor(col1, alpha.f = alpha_fill),
                  adjustcolor(col2, alpha.f = alpha_fill),
                  adjustcolor(col3, alpha.f = alpha_fill)),
         border = c(adjustcolor(col1, alpha.f = alpha_border),
                    adjustcolor(col2, alpha.f = alpha_border),
                    adjustcolor(col3, alpha.f = alpha_border)),
         bty = "n")
  
  mtext(main_label, side = 3, line = 0.5, outer = TRUE, cex = 1.4)
}

# 1) Global histogram
pdf(file.path(out_global, "all_hist.pdf"), width = 9, height = 6)
par(mar = c(4, 4, 1, 1), oma = c(3.5, 4.5, 2, 1))
plot_overlay_hist(x1_all, x2_all, x3_all,
                  main_label = "Distribution of velocity drivers, Global")
dev.off()

# 2) Year-specific histograms
fric1_by_year <- collect_by_year(fric1_list)
fric2_by_year <- collect_by_year(fric2_list)
fric3_by_year <- collect_by_year(fric3_list)

years_union <- sort(unique(as.numeric(c(names(fric1_by_year),
                                        names(fric2_by_year),
                                        names(fric3_by_year)))))
years_union <- years_union[is.finite(years_union)]

for (yr in years_union) {
  ychr <- as.character(yr)
  x1 <- fric1_by_year[[ychr]] %||% numeric(0)
  x2 <- fric2_by_year[[ychr]] %||% numeric(0)
  x3 <- fric3_by_year[[ychr]] %||% numeric(0)
  
  pdf(file.path(out_years, paste0("year_", ychr, "_hist.pdf")), width = 9, height = 6)
  par(mar = c(4, 4, 1, 1), oma = c(3.5, 4.5, 2, 1))
  plot_overlay_hist(x1, x2, x3,
                    main_label = paste0("Distribution of velocity drivers, Year ", ychr))
  dev.off()
}

# 3) Elevation-specific histograms
fric1_by_elev <- collect_by_elev(fric1_list)
fric2_by_elev <- collect_by_elev(fric2_list)
fric3_by_elev <- collect_by_elev(fric3_list)

elevs_union <- sort(unique(as.numeric(c(names(fric1_by_elev),
                                        names(fric2_by_elev),
                                        names(fric3_by_elev)))))
elevs_union <- elevs_union[is.finite(elevs_union)]

for (el in elevs_union) {
  echr <- as.character(el)
  x1 <- fric1_by_elev[[echr]] %||% numeric(0)
  x2 <- fric2_by_elev[[echr]] %||% numeric(0)
  x3 <- fric3_by_elev[[echr]] %||% numeric(0)
  
  pdf(file.path(out_elevs, paste0("elev_", echr, "_hist.pdf")), width = 9, height = 6)
  par(mar = c(4, 4, 1, 1), oma = c(3.25, 1, 2, 1))
  plot_overlay_hist(x1, x2, x3,
                    main_label = paste0("Distribution of velocity drivers, Elevation ", echr))
  dev.off()
  
  png(file.path(out_elevs, paste0("elev_", echr, "_hist.png")),
      width = 9, height = 6, units = "in", res = 600)
  par(mar = c(4, 4, 1, 1), oma = c(3.25, 1, 2, 1))
  plot_overlay_hist(x1, x2, x3,
                    main_label = paste0("Distribution of velocity drivers, Elevation ", echr))
  dev.off()
}

# 4) Flowline-specific histograms
flow_ids <- sort(unique(c(names(fric1_list), names(fric2_list), names(fric3_list))))

for (fid in flow_ids) {
  parts <- strsplit(fid, "_")[[1]]
  glacier_num <- suppressWarnings(as.integer(parts[2]))
  flowline_num <- suppressWarnings(as.integer(parts[3]))
  glacier_name <- if (is.finite(glacier_num)) get_glacier_name(glacier_num) else "Unknown Glacier"
  
  m1 <- fric1_list[[fid]]
  m2 <- fric2_list[[fid]]
  m3 <- fric3_list[[fid]]
  
  x1 <- if (!is.null(m1)) as.numeric(m1) else numeric(0)
  x2 <- if (!is.null(m2)) as.numeric(m2) else numeric(0)
  x3 <- if (!is.null(m3)) as.numeric(m3) else numeric(0)
  x1 <- x1[is.finite(x1) & x1 > 0]
  x2 <- x2[is.finite(x2) & x2 > 0]
  x3 <- x3[is.finite(x3) & x3 > 0]
  
  title_main <- paste0("Distribution of velocity drivers\n",
                       "Glacier ", glacier_num, ", Flowline ", flowline_num,
                       " \u2013 ", glacier_name)
  
  pdf(file.path(out_flows, paste0(fid, "_hist.pdf")), width = 9, height = 6)
  par(mar = c(4, 4, 1, 1), oma = c(3.5, 4.5, 2, 1))
  plot_overlay_hist(x1, x2, x3, main_label = title_main)
  dev.off()
}

cat("Done. Histograms saved to:\n",
    "- Global: ", out_global, "\n",
    "- Years:  ", out_years,  "\n",
    "- Elevs:  ", out_elevs,  "\n",
    "- Flows:  ", out_flows,  "\n", sep = "")
