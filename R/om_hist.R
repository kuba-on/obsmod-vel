########################################

# The code below compares distributions of ratios of seasonality drivers obtained from different friction laws. The processing includes:
# - Reading CSV files containing all ratio data for different friction laws.
# - Plotting overlaid histograms with transparency, grouped by year, elevation, flowline, and globally with all data combined.
# - Plots contain vertical lines at the mean of each distribution.

########################################
# Directories (unchanged)
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

# Manuscript output dir
out_manus <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Histograms/Ratios/Manuscript"
suppressWarnings(dir.create(out_manus, recursive = TRUE, showWarnings = FALSE))

# Glacier name lookup (unchanged)
csv_lookup <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/New Points v3/Input/Box Coordinates/box_sp_all_v3.csv"
lookup_df <- tryCatch(read.csv(csv_lookup), error = function(e) NULL)

########################################
# Helpers (unchanged except where noted)
read_ratio_folder <- function(folder_path, friction_tag) {
  patt <- paste0("^gl_\\d+_\\d+_ratio_", friction_tag, "\\.csv$")
  files <- list.files(folder_path, pattern = patt, full.names = TRUE)
  out <- list()
  for (f in files) {
    df <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) next
    if (ncol(df) > 1 && (is.character(df[[1]]) || is.numeric(df[[1]]))) {
      maybe_year <- suppressWarnings(as.numeric(df[[1]]))
      if (sum(!is.na(maybe_year)) > 0) {
        rownames(df) <- as.character(df[[1]])
        df <- df[, -1, drop = FALSE]
      }
    }
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

flatten_vec <- function(mat_list) {
  if (length(mat_list) == 0) return(numeric(0))
  v <- unlist(mat_list, use.names = FALSE)
  v <- as.numeric(v)
  v <- v[is.finite(v) & !is.na(v) & v > 0]
  v
}

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

get_glacier_name <- function(glacier_num) {
  if (is.null(lookup_df)) return("Unknown Glacier")
  nm <- lookup_df$glacier_name[lookup_df$feature_ID == as.integer(glacier_num)]
  if (length(nm) == 0) "Unknown Glacier" else as.character(nm[1])
}

########################################
# Read data (unchanged)
fric1_list <- read_ratio_folder(dir_fric1, "fric1")
fric2_list <- read_ratio_folder(dir_fric2, "fric2")
fric3_list <- read_ratio_folder(dir_fric3, "fric3")

# Global breaks (unchanged)
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

########################################
# Plot overlay function (only minimal edits)
plot_overlay_hist <- function(x1, x3, x2,
                              breaks = global_breaks,
                              col1 = "#648FFF", col2 = "#DC267F", col3 = "#FFB000",
                              alpha_fill = 0.35, alpha_border = 0.9,
                              main_label = "",
                              legend_pos = "topleft",
                              ylim_override = NULL) {
  allx <- c(x1, x2, x3)
  allx <- allx[is.finite(allx) & allx > 0]
  if (length(allx) == 0) {
    plot.new(); box(); mtext("No data", side = 3, line = 0.5, cex = 1)
    return(invisible(NULL))
  }
  lx1 <- log10(x1[x1 > 0]); h1 <- hist(lx1, breaks = breaks, plot = FALSE)
  lx2 <- log10(x2[x2 > 0]); h2 <- hist(lx2, breaks = breaks, plot = FALSE)
  lx3 <- log10(x3[x3 > 0]); h3 <- hist(lx3, breaks = breaks, plot = FALSE)
  ylim_max <- max(c(h1$counts, h2$counts, h3$counts, 1))
  if (!is.null(ylim_override) && is.finite(ylim_override) && ylim_override > 0) ylim_max <- ylim_override
  
  plot(NA, xlim = range(breaks), ylim = c(0, ylim_max),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  box()
  
  # Bars
  rect(h1$breaks[-length(h1$breaks)], 0, h1$breaks[-1], h1$counts,
       col = adjustcolor(col1, alpha.f = alpha_fill),
       border = adjustcolor(col1, alpha.f = alpha_border))
  rect(h2$breaks[-length(h2$breaks)], 0, h2$breaks[-1], h2$counts,
       col = adjustcolor(col2, alpha.f = alpha_fill),
       border = adjustcolor(col2, alpha.f = alpha_border))
  rect(h3$breaks[-length(h3$breaks)], 0, h3$breaks[-1], h3$counts,
       col = adjustcolor(col3, alpha.f = alpha_fill),
       border = adjustcolor(col3, alpha.f = alpha_border))
  
  # Axes
  axis(2, las = 1)
  xticks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = 1)
  axis(1, at = xticks, labels = parse(text = paste0("10^", xticks)))
  mtext("count", side = 2, line = 3, cex = 1.1)
  mtext("velocity driver, R", side = 1, line = 4.5, cex = 1.1)
  
  # Secondary axis (basal / front)
  usr <- par("usr"); xpd_prev <- par("xpd"); par(xpd = TRUE)
  axis_bottom_y <- usr[3] - 0.15 * diff(usr[3:4]); tick_length <- 0.03 * diff(usr[3:4])
  x_min_tick <- min(xticks, na.rm = TRUE); x_max_tick <- max(xticks, na.rm = TRUE)
  segments(x_min_tick, axis_bottom_y, x_max_tick, axis_bottom_y, col = "black")
  for (tx in c(x_min_tick, 0, x_max_tick)) {
    segments(tx, axis_bottom_y, tx, axis_bottom_y + tick_length, col = "black")
  }
  mid_left  <- (x_min_tick + 0) / 2; mid_right <- (x_max_tick + 0) / 2
  text(mid_left,  axis_bottom_y - 0.035 * diff(usr[3:4]), labels = "basal", cex = 1)
  text(mid_right, axis_bottom_y - 0.035 * diff(usr[3:4]), labels = "front", cex = 1)
  par(xpd = xpd_prev)
  
  # Medians + threshold
  m1 <- if (length(x1)) median(x1, na.rm = TRUE) else NA
  m2 <- if (length(x2)) median(x2, na.rm = TRUE) else NA
  m3 <- if (length(x3)) median(x3, na.rm = TRUE) else NA
  if (is.finite(m1) && m1 > 0) abline(v = log10(m1), col = col1, lty = 3, lwd = 2)
  if (is.finite(m2) && m2 > 0) abline(v = log10(m2), col = col2, lty = 3, lwd = 2)
  if (is.finite(m3) && m3 > 0) abline(v = log10(m3), col = col3, lty = 3, lwd = 2)
  abline(v = 0, col = "black", lwd = 1)
  
  legend(legend_pos,
         legend = c("Budd, m=1", "Budd, m=6", "Schoof, m=3"),
         fill = c(adjustcolor(col1, alpha.f = alpha_fill),
                  adjustcolor(col2, alpha.f = alpha_fill),
                  adjustcolor(col3, alpha.f = alpha_fill)),
         border = c(adjustcolor(col1, alpha.f = alpha_border),
                    adjustcolor(col2, alpha.f = alpha_border),
                    adjustcolor(col3, alpha.f = alpha_border)),
         bty = "n", ncol = 1)
  mtext(main_label, side = 3, line = 0.5, outer = TRUE, cex = 1.4)
}

########################################
# Existing outputs (unchanged use of function; labels fixed by function)
# 1) Global histogram
pdf(file.path(out_global, "all_hist.pdf"), width = 9, height = 6)
par(mar = c(4, 4, 1, 1), oma = c(3.5, 4.5, 2, 1))
plot_overlay_hist(x1_all, x2_all, x3_all,
                  main_label = "Distribution of velocity drivers, Global")
dev.off()

# 2) Year-specific histograms (now using fixed Y across all years)
fric1_by_year <- collect_by_year(fric1_list)
fric2_by_year <- collect_by_year(fric2_list)
fric3_by_year <- collect_by_year(fric3_list)
years_union <- sort(unique(as.numeric(c(names(fric1_by_year),
                                        names(fric2_by_year),
                                        names(fric3_by_year)))))
years_union <- years_union[is.finite(years_union)]

# Precompute global ymax for years
year_ylim_max <- {
  ymax <- 1
  for (yr in years_union) {
    ychr <- as.character(yr)
    x1 <- fric1_by_year[[ychr]] %||% numeric(0)
    x2 <- fric2_by_year[[ychr]] %||% numeric(0)
    x3 <- fric3_by_year[[ychr]] %||% numeric(0)
    if (length(c(x1,x2,x3))) {
      h1 <- hist(log10(x1[x1>0]), breaks = global_breaks, plot = FALSE)
      h2 <- hist(log10(x2[x2>0]), breaks = global_breaks, plot = FALSE)
      h3 <- hist(log10(x3[x3>0]), breaks = global_breaks, plot = FALSE)
      ymax <- max(ymax, h1$counts, h2$counts, h3$counts, na.rm = TRUE)
    }
  }
  ymax
}

for (yr in years_union) {
  ychr <- as.character(yr)
  x1 <- fric1_by_year[[ychr]] %||% numeric(0)
  x2 <- fric2_by_year[[ychr]] %||% numeric(0)
  x3 <- fric3_by_year[[ychr]] %||% numeric(0)
  
  pdf(file.path(out_years, paste0("year_", ychr, "_hist.pdf")), width = 9, height = 6)
  par(mar = c(4, 4, 1, 1), oma = c(3.5, 4.5, 2, 1))
  plot_overlay_hist(x1, x2, x3,
                    main_label = paste0("Distribution of velocity drivers, Year ", ychr),
                    ylim_override = year_ylim_max)
  dev.off()
}

# 3) Elevation-specific histograms (now using fixed Y across all elevations)
fric1_by_elev <- collect_by_elev(fric1_list)
fric2_by_elev <- collect_by_elev(fric2_list)
fric3_by_elev <- collect_by_elev(fric3_list)
elevs_union <- sort(unique(as.numeric(c(names(fric1_by_elev),
                                        names(fric2_by_elev),
                                        names(fric3_by_elev)))))
elevs_union <- elevs_union[is.finite(elevs_union)]

# Precompute global ymax for elevations
elev_ylim_max <- {
  ymax <- 1
  for (el in elevs_union) {
    echr <- as.character(el)
    x1 <- fric1_by_elev[[echr]] %||% numeric(0)
    x2 <- fric2_by_elev[[echr]] %||% numeric(0)
    x3 <- fric3_by_elev[[echr]] %||% numeric(0)
    if (length(c(x1,x2,x3))) {
      h1 <- hist(log10(x1[x1>0]), breaks = global_breaks, plot = FALSE)
      h2 <- hist(log10(x2[x2>0]), breaks = global_breaks, plot = FALSE)
      h3 <- hist(log10(x3[x3>0]), breaks = global_breaks, plot = FALSE)
      ymax <- max(ymax, h1$counts, h2$counts, h3$counts, na.rm = TRUE)
    }
  }
  ymax
}

for (el in elevs_union) {
  echr <- as.character(el)
  x1 <- fric1_by_elev[[echr]] %||% numeric(0)
  x2 <- fric2_by_elev[[echr]] %||% numeric(0)
  x3 <- fric3_by_elev[[echr]] %||% numeric(0)
  
  pdf(file.path(out_elevs, paste0("elev_", echr, "_hist.pdf")), width = 9, height = 6)
  par(mar = c(4, 4, 1, 1), oma = c(3.25, 1, 2, 1))
  plot_overlay_hist(x1, x2, x3,
                    main_label = paste0("Distribution of velocity drivers, Elevation ", echr),
                    ylim_override = elev_ylim_max)
  dev.off()
  
  png(file.path(out_elevs, paste0("elev_", echr, "_hist.png")),
      width = 9, height = 6, units = "in", res = 600)
  par(mar = c(4, 4, 1, 1), oma = c(3.25, 1, 2, 1))
  plot_overlay_hist(x1, x2, x3,
                    main_label = paste0("Distribution of velocity drivers, Elevation ", echr),
                    ylim_override = elev_ylim_max)
  dev.off()
}

# 4) Flowline-specific histograms (unchanged)
flow_ids <- sort(unique(c(names(fric1_list), names(fric2_list), names(fric3_list))))
for (fid in flow_ids) {
  parts <- strsplit(fid, "_")[[1]]
  glacier_num <- suppressWarnings(as.integer(parts[2]))
  flowline_num <- suppressWarnings(as.integer(parts[3]))
  glacier_name <- if (is.finite(glacier_num)) get_glacier_name(glacier_num) else "Unknown Glacier"
  m1 <- fric1_list[[fid]]; m2 <- fric2_list[[fid]]; m3 <- fric3_list[[fid]]
  x1 <- if (!is.null(m1)) as.numeric(m1) else numeric(0)
  x2 <- if (!is.null(m2)) as.numeric(m2) else numeric(0)
  x3 <- if (!is.null(m3)) as.numeric(m3) else numeric(0)
  x1 <- x1[is.finite(x1) & x1 > 0]; x2 <- x2[is.finite(x2) & x2 > 0]; x3 <- x3[is.finite(x3) & x3 > 0]
  title_main <- paste0("Distribution of velocity drivers\n",
                       "Glacier ", glacier_num, ", Flowline ", flowline_num,
                       " \u2013 ", glacier_name)
  pdf(file.path(out_flows, paste0(fid, "_hist.pdf")), width = 9, height = 6)
  par(mar = c(4, 4, 1, 1), oma = c(3.5, 4.5, 2, 1))
  plot_overlay_hist(x1, x2, x3, main_label = title_main)
  dev.off()
}

########################################
# 5) Friction-law-specific histograms (by elevation and by year) — minimal edits:
extend_palette <- function(base_cols, n) colorRampPalette(base_cols)(n)

plot_multi_hist <- function(values_list, group_label, friction_label, base_cols, ylim_override = NULL) {
  if (length(values_list) == 0) return(NULL)
  all_vals <- unlist(values_list, use.names = FALSE)
  all_vals <- all_vals[is.finite(all_vals) & all_vals > 0]
  if (length(all_vals) == 0) { plot.new(); box(); mtext("No data", side=3,line=0.5,cex=1); return(invisible(NULL)) }
  
  n_groups <- length(values_list)
  group_names <- names(values_list)
  cols <- extend_palette(base_cols, n_groups)
  breaks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = 0.1)
  
  # precompute
  hists <- lapply(values_list, function(v) hist(log10(v[v>0]), breaks = breaks, plot = FALSE))
  ylim_max <- max(unlist(lapply(hists, function(h) h$counts)), 1)
  if (!is.null(ylim_override) && is.finite(ylim_override) && ylim_override > 0) ylim_max <- ylim_override
  
  plot(NA, xlim = range(breaks), ylim = c(0, ylim_max),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  box(); axis(2, las = 1)
  xticks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = 1)
  axis(1, at = xticks, labels = parse(text = paste0("10^", xticks)))
  mtext("count", side = 2, line = 3, cex = 1.1)
  mtext("velocity driver, R", side = 1, line = 4.5, cex = 1.1)
  
  for (i in seq_along(hists)) {
    h <- hists[[i]]
    rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts,
         col = adjustcolor(cols[i], alpha.f = 0.35),
         border = adjustcolor(cols[i], alpha.f = 0.9))
  }
  
  # Secondary axis (basal/front)
  usr <- par("usr"); xpd_prev <- par("xpd"); par(xpd = TRUE)
  axis_bottom_y <- usr[3] - 0.15 * diff(usr[3:4]); tick_length <- 0.03 * diff(usr[3:4])
  x_min_tick <- min(xticks, na.rm = TRUE); x_max_tick <- max(xticks, na.rm = TRUE)
  segments(x_min_tick, axis_bottom_y, x_max_tick, axis_bottom_y, col = "black")
  for (tx in c(x_min_tick, 0, x_max_tick)) segments(tx, axis_bottom_y, tx, axis_bottom_y + tick_length, col = "black")
  mid_left  <- (x_min_tick + 0)/2; mid_right <- (x_max_tick + 0)/2
  text(mid_left,  axis_bottom_y - 0.035 * diff(usr[3:4]), labels = "basal", cex = 1)
  text(mid_right, axis_bottom_y - 0.035 * diff(usr[3:4]), labels = "front", cex = 1)
  par(xpd = xpd_prev)
  
  # Medians and threshold
  for (i in seq_along(values_list)) {
    m <- median(values_list[[i]], na.rm = TRUE)
    if (is.finite(m) && m > 0) abline(v = log10(m), col = cols[i], lty = 3, lwd = 1.8)
  }
  abline(v = 0, col = "black", lwd = 1)
  
  legend("topleft",
         legend = c("(m a.s.l.)", group_names),
         fill = c(NA, adjustcolor(cols, alpha.f = 0.35)),
         border = c(NA, adjustcolor(cols, alpha.f = 0.9)),
         bty = "n", ncol = 1,
         text.col = c("black", rep("black", length(group_names))))
  mtext(paste0("Distribution of velocity drivers by ", group_label, ", ", friction_label),
        side = 3, line = 0.5, outer = TRUE, cex = 1.4)
}

fric_sets <- list(fric1_list = fric1_list, fric2_list = fric2_list, fric3_list = fric3_list)
fric_labels <- c("Budd, m=1","Schoof, m=3","Budd, m=6")
fric_ids <- c("fric1","fric2","fric3")
base_cols <- c("#648FFF", "#DC267F", "#FFB000")

# Global ymax for the grouped (multi) plots: reuse elev_ylim_max / year_ylim_max
for (i in seq_along(fric_sets)) {
  flist <- fric_sets[[i]]; friction_label <- fric_labels[i]; friction_id <- fric_ids[i]
  by_elev <- collect_by_elev(flist); if (length(by_elev) == 0) next
  pdf(file.path(out_elevs, paste0("all_elevs_", friction_id, ".pdf")), width = 9, height = 6)
  par(mar = c(4,4,1,1), oma = c(3.5,4.5,2,1))
  plot_multi_hist(by_elev, "elevation", friction_label, base_cols, ylim_override = elev_ylim_max)
  dev.off()
}
for (i in seq_along(fric_sets)) {
  flist <- fric_sets[[i]]; friction_label <- fric_labels[i]; friction_id <- fric_ids[i]
  by_year <- collect_by_year(flist); if (length(by_year) == 0) next
  pdf(file.path(out_years, paste0("all_years_", friction_id, ".pdf")), width = 9, height = 6)
  par(mar = c(4,4,1,1), oma = c(3.5,4.5,2,1))
  plot_multi_hist(by_year, "year", friction_label, base_cols, ylim_override = year_ylim_max)
  dev.off()
}

########################################
# 6) Manuscript figure: 5 stacked elevation panels
#    - Only bottom shows x-axis
#    - Only middle shows y-axis label ("count")
#    - Legend with bold title "<elev> m a.s.l."
#    - No figure title
#    - Uses fixed Y (same elev_ylim_max)
########################################

## ============================
## Manuscript 5-panel elevation figure (PDF + PNG)
## ============================

####### --- PUBLICATION PANEL (stacked elevation histograms; equal heights) --- #######

# Output folder & filenames
out_manus <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Histograms/Ratios/Manuscript"
suppressWarnings(dir.create(out_manus, recursive = TRUE, showWarnings = FALSE))
pub_pdf <- file.path(out_manus, "hist_text.pdf")
pub_png <- file.path(out_manus, "hist_text.png")

# Elevations to plot (choose any; default 5 as requested)
sel_elevs <- c(200, 600, 1200, 1600)  # <-- edit this vector to pick which panels to show

# Precomputed per-elevation lists (from earlier code)
fric1_by_elev <- collect_by_elev(fric1_list)
fric2_by_elev <- collect_by_elev(fric2_list)
fric3_by_elev <- collect_by_elev(fric3_list)

# Colors and legend labels (unchanged)
col1 <- "#648FFF"; col2 <- "#DC267F"; col3 <- "#FFB000"
leg_labels <- c("Budd, m=1", "Budd, m=6", "Schoof, m=3")

# Helper: fetch values for an elevation
.get_vals_el <- function(el) {
  key <- as.character(el)
  x1 <- fric1_by_elev[[key]]; x2 <- fric2_by_elev[[key]]; x3 <- fric3_by_elev[[key]]
  x1 <- if (length(x1)) x1[is.finite(x1) & x1 > 0] else numeric(0)
  x2 <- if (length(x2)) x2[is.finite(x2) & x2 > 0] else numeric(0)
  x3 <- if (length(x3)) x3[is.finite(x3) & x3 > 0] else numeric(0)
  list(x1 = x1, x2 = x2, x3 = x3)
}

# Fixed y-maximum across ALL selected panels (as requested)
.panel_max <- (function() {
  ymax <- 0
  for (el in sel_elevs) {
    v <- .get_vals_el(el)
    if (length(v$x1)) ymax <- max(ymax, max(hist(log10(v$x1), breaks = global_breaks, plot = FALSE)$counts, na.rm = TRUE))
    if (length(v$x2)) ymax <- max(ymax, max(hist(log10(v$x2), breaks = global_breaks, plot = FALSE)$counts, na.rm = TRUE))
    if (length(v$x3)) ymax <- max(ymax, max(hist(log10(v$x3), breaks = global_breaks, plot = FALSE)$counts, na.rm = TRUE))
  }
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  ymax
})()

# One panel draw (no figure title)
plot_one_panel_pub <- function(el, show_xaxis = FALSE, show_ylabel = FALSE, ylim_max = .panel_max) {
  v  <- .get_vals_el(el)
  lx1 <- if (length(v$x1)) log10(v$x1) else numeric(0)
  lx2 <- if (length(v$x2)) log10(v$x2) else numeric(0)
  lx3 <- if (length(v$x3)) log10(v$x3) else numeric(0)
  
  h1 <- hist(lx1, breaks = global_breaks, plot = FALSE)
  h2 <- hist(lx2, breaks = global_breaks, plot = FALSE)
  h3 <- hist(lx3, breaks = global_breaks, plot = FALSE)
  
  # Identical margins for all panels (ensures equal visual heights)
  plot(NA, xlim = range(global_breaks), ylim = c(0, ylim_max),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", xaxs = "r", yaxs = "r")
  box()
  
  # Bars
  rect(h1$breaks[-length(h1$breaks)], 0, h1$breaks[-1], h1$counts,
       col = adjustcolor(col1, 0.35), border = adjustcolor(col1, 0.9))
  rect(h2$breaks[-length(h2$breaks)], 0, h2$breaks[-1], h2$counts,
       col = adjustcolor(col2, 0.35), border = adjustcolor(col2, 0.9))
  rect(h3$breaks[-length(h3$breaks)], 0, h3$breaks[-1], h3$counts,
       col = adjustcolor(col3, 0.35), border = adjustcolor(col3, 0.9))
  
  # Y axis + label (middle panel only)
  axis(2, las = 1)
  if (show_ylabel) mtext("count", side = 2, line = -2, cex = 1.1, outer = TRUE)
  
  # Medians (log domain)
  m1 <- if (length(v$x1)) median(v$x1, na.rm = TRUE) else NA
  m2 <- if (length(v$x2)) median(v$x2, na.rm = TRUE) else NA
  m3 <- if (length(v$x3)) median(v$x3, na.rm = TRUE) else NA
  if (is.finite(m1) && m1 > 0) abline(v = log10(m1), col = col1, lty = 3, lwd = 2)
  if (is.finite(m2) && m2 > 0) abline(v = log10(m2), col = col2, lty = 3, lwd = 2)
  if (is.finite(m3) && m3 > 0) abline(v = log10(m3), col = col3, lty = 3, lwd = 2)
  
  # Threshold at 10^0
  abline(v = 0, col = "black", lwd = 1)
  
  # Legend INSIDE (top-left) with bold elevation label
  usr <- par("usr")
  lg_x <- usr[1] + 0.02 * diff(usr[1:2])
  title_y <- usr[4] - 0.06 * diff(usr[3:4])  # slightly lower to avoid clipping
  text(lg_x, title_y, paste0(el, " m a.s.l."), adj = c(0, 1), font = 2, xpd = FALSE)
  legend(lg_x, title_y - 0.03 * diff(usr[3:4]),
         legend = leg_labels,
         fill   = c(adjustcolor(col1, 0.35), adjustcolor(col2, 0.35), adjustcolor(col3, 0.35)),
         border = c(adjustcolor(col1, 0.9),  adjustcolor(col2, 0.9),  adjustcolor(col3, 0.9)),
         bty = "n", ncol = 1, xpd = FALSE)
  
  # X-axis and secondary inverted axis only on bottom panel
  if (show_xaxis) {
    xticks <- seq(floor(global_xrange[1]), ceiling(global_xrange[2]), by = 1)
    axis(1, at = xticks, labels = parse(text = paste0("10^", xticks)))
    mtext("seasonality index, S", side = 1, line = 4.5, cex = 1.1)
    
    # Secondary inverted axis (basal / front)
    xpd_prev <- par("xpd"); par(xpd = NA)   # allow drawing into the outer margin
    usr <- par("usr")
    axis_bottom_y <- usr[3] - 0.18 * diff(usr[3:4])  # positioned into outer margin
    tick_len      <- 0.03 * diff(usr[3:4])
    x_min_tick <- min(xticks, na.rm = TRUE); x_max_tick <- max(xticks, na.rm = TRUE)
    segments(x_min_tick, axis_bottom_y, x_max_tick, axis_bottom_y)
    for (tx in c(x_min_tick, 0, x_max_tick)) {
      segments(tx, axis_bottom_y, tx, axis_bottom_y + tick_len)
    }
    mid_left  <- (x_min_tick + 0) / 2
    mid_right <- (x_max_tick + 0) / 2
    text(mid_left,  axis_bottom_y - 0.045 * diff(usr[3:4]), "basal",  cex = 1)
    text(mid_right, axis_bottom_y - 0.045 * diff(usr[3:4]), "front",  cex = 1)
    par(xpd = xpd_prev)
  }
}

# Panel-count helpers (to position y-label for even panel counts)
n_panels    <- length(sel_elevs)
even_panels <- (n_panels %% 2 == 0)

# Draw and save (NO cairo; identical margins in both devices to ensure equal heights)
draw_publication_panels <- function() {
  width_in   <- 9
  panel_h_in <- 2.5
  total_h_in <- panel_h_in * length(sel_elevs)  # stack height
  
  # Unified margins for BOTH devices:
  inner_mar <- c(1.0, 4.0, 0.6, 1.0)  # c(bottom, left, top, right)
  # Add left outer margin only if we place an outer y-label (even # of panels)
  outer_mar <- if (even_panels) c(7.0, 1.5, 0.0, 0.0) else c(7.0, 1.5, 0.0, 0.0)
  
  # PDF
  if (file.exists(pub_pdf)) unlink(pub_pdf, force = TRUE)
  pdf(pub_pdf, width = width_in, height = total_h_in)
  par(mfrow = c(length(sel_elevs), 1), mar = inner_mar, oma = outer_mar, xpd = FALSE)
  for (i in seq_along(sel_elevs)) {
    show_x <- (i == length(sel_elevs))
    show_y <- (!even_panels && i == ceiling(n_panels/2))
    plot_one_panel_pub(sel_elevs[i], show_xaxis = show_x, show_ylabel = show_y)
  }
  # Outer y-label centered on figure when even number of panels
  if (even_panels) mtext("count", side = 2, line = 0, cex = 1.1, outer = TRUE)
  dev.off()
  
  # PNG
  if (file.exists(pub_png)) unlink(pub_png, force = TRUE)
  png(pub_png, width = width_in, height = total_h_in, units = "in", res = 600)
  par(mfrow = c(length(sel_elevs), 1), mar = inner_mar, oma = outer_mar, xpd = FALSE)
  for (i in seq_along(sel_elevs)) {
    show_x <- (i == length(sel_elevs))
    show_y <- (!even_panels && i == ceiling(n_panels/2))
    plot_one_panel_pub(sel_elevs[i], show_xaxis = show_x, show_ylabel = show_y)
  }
  if (even_panels) mtext("count", side = 2, line = 2, cex = 1.1, outer = TRUE)
  dev.off()
}

# Build both files
draw_publication_panels()

########################################
cat("Done. Histograms saved to:\n",
    "- Global: ", out_global, "\n",
    "- Years:  ", out_years,  "\n",
    "- Elevs:  ", out_elevs,  "\n",
    "- Flows:  ", out_flows,  "\n",
    "- Manuscript: ", out_manus, "\n", sep = "")