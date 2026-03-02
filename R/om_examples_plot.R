##### CUSTOM OBS–MOD–RES PANEL GRID — 12 panels with margin labels #####

# --- Output ---
out_pdf <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Manuscript/Graphs/Examples/glacier_examples.pdf"

# --- Grid (3×4 = 12) preserving original per-panel size (25×10 over 21×6 panels) ---
n_rows <- 4
n_cols <- 4
cell_w <- 25/21
cell_h <- 7.5/6
pdf_width  <- n_cols * cell_w
pdf_height <- n_rows * cell_h

# --- Inputs on disk ---
res_base   <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/data/RDS/res_fric3"
ratio_base <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric3/Ratios"
csv_path   <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/New Points v3/Input/Box Coordinates/box_sp_all_v3.csv"

# --- Selections (your test set) ---
selections <- list(
  # ---- Row 1 (Inngia) ----
  list(idx=1,  gl=22,  fl=1, yr=2017, el=100),
  list(idx=2,  gl=22,  fl=1, yr=2017, el=300),
  list(idx=3,  gl=22,  fl=1, yr=2017, el=500),
  list(idx=4,  gl=22,  fl=1, yr=2017, el=1000),
  
  # ---- Row 2 (Kangilernata) ----
  list(idx=5,  gl=81, fl=1, yr=2019, el=100),
  list(idx=6,  gl=81, fl=1, yr=2019, el=600),
  list(idx=7,  gl=81, fl=1, yr=2019, el=1200),
  list(idx=8,  gl=81, fl=1, yr=2019, el=1900),
  
  # ---- Row 3 (Hayes N)) ----
  list(idx=9,  gl=213, fl=4, yr=2019, el=100),
  list(idx=10, gl=213, fl=4, yr=2019, el=200),
  list(idx=11, gl=213, fl=4, yr=2019, el=400),
  list(idx=12, gl=213, fl=4, yr=2019, el=600),
  
  # ---- Row 4 (Jakobshavn) ----
  list(idx=13, gl=217, fl=6, yr=2019, el=400),
  list(idx=14, gl=217, fl=6, yr=2019, el=700),
  list(idx=15, gl=217, fl=6, yr=2019, el=1100),
  list(idx=16, gl=217, fl=6, yr=2019, el=1600)
)

# =========================
# Helpers
# =========================
.letter_for_index <- function(i) letters[((i - 1) %% length(letters)) + 1]

.rds_cache   <- new.env(parent = emptyenv())
.ratio_cache <- new.env(parent = emptyenv())
.name_cache  <- new.env(parent = emptyenv())

.load_residuals <- function(gl, fl) {
  key <- paste(gl, fl, sep = "_")
  if (exists(key, envir = .rds_cache, inherits = FALSE)) return(get(key, envir = .rds_cache))
  f <- file.path(res_base, sprintf("gl_%d_%d_data_res_fric3.rds", gl, fl))
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f); assign(key, x, envir = .rds_cache); x
}

.load_ratio_matrix <- function(gl, fl) {
  key <- paste(gl, fl, sep = "_")
  if (exists(key, envir = .ratio_cache, inherits = FALSE)) return(get(key, envir = .ratio_cache))
  f <- file.path(ratio_base, sprintf("gl_%d_%d_ratio_fric3.csv", gl, fl))
  if (!file.exists(f)) return(NULL)
  m <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  assign(key, m, envir = .ratio_cache); m
}

.get_ratio_value_from_csv <- function(gl, fl, yr, el) {
  mat <- .load_ratio_matrix(gl, fl)
  if (is.null(mat)) return(NA_real_)
  yr <- as.character(yr); el <- as.character(el)
  if (!yr %in% rownames(mat) || !el %in% colnames(mat)) return(NA_real_)
  val <- suppressWarnings(as.numeric(mat[yr, el]))
  if (!is.finite(val)) NA_real_ else val
}

.get_glacier_name <- function(gl) {
  key <- as.character(gl)
  if (exists(key, envir = .name_cache, inherits = FALSE)) return(get(key, envir = .name_cache))
  info <- try(read.csv(csv_path, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(info, "try-error")) {
    assign(key, key, envir = .name_cache)
    return(key)
  }
  nm <- info$glacier_name[info$feature_ID == gl]
  nm <- if (length(nm) && !is.na(nm[1])) nm[1] else key
  assign(key, nm, envir = .name_cache)
  nm
}

.fit_spline <- function(dates, y, df = 20, spar = NULL) {
  ok <- is.finite(y) & !is.na(dates)
  if (sum(ok) < 4) return(NULL)
  x <- as.numeric(dates[ok]); y <- y[ok]
  fit <- try(
    if (is.null(spar)) smooth.spline(x, y, df = df) else smooth.spline(x, y, df = df, spar = spar),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NULL)
  lev <- fit$lev
  res <- (fit$yin - fit$y) / pmax(1 - lev, .Machine$double.eps)
  sigma <- sqrt(stats::var(res, na.rm = TRUE))
  upper <- fit$y + 2.0 * sigma * sqrt(lev)
  lower <- fit$y - 2.0 * sigma * sqrt(lev)
  data.frame(date = as.Date(fit$x, origin = "1970-01-01"),
             y = fit$y,
             upper = upper,
             lower = lower)
}

.compute_ratio_from_splines <- function(so, sm, sr, year) {
  if (is.null(so) || is.null(sm) || is.null(sr)) return(NA_real_)
  so <- subset(so, date >= as.Date(paste0(year, "-01-01")) & date <= as.Date(paste0(year, "-12-31")))
  sm <- subset(sm, date >= as.Date(paste0(year, "-01-01")) & date <= as.Date(paste0(year, "-12-31")))
  sr <- subset(sr, date >= as.Date(paste0(year, "-01-01")) & date <= as.Date(paste0(year, "-12-31")))
  if (nrow(so) == 0 || nrow(sm) == 0 || nrow(sr) == 0) return(NA_real_)
  common_dates <- Reduce(intersect, list(so$date, sm$date, sr$date))
  if (!length(common_dates)) return(NA_real_)
  so_v <- so$y[so$date %in% common_dates]
  sm_v <- sm$y[sm$date %in% common_dates]
  sr_v <- sr$y[sr$date %in% common_dates]
  valid <- complete.cases(so_v, sm_v, sr_v)
  if (!any(valid)) return(NA_real_)
  area_front <- sum(abs(so_v[valid] - sr_v[valid]), na.rm = TRUE)
  area_basal <- sum(abs(so_v[valid] - sm_v[valid]), na.rm = TRUE)
  if (!is.finite(area_basal) || area_basal == 0) return(NA_real_)
  area_front / area_basal
}

ratio_palette <- colorRampPalette(c("#364B9A", "#EAECCC", "#A50026"))(999)
.map_ratio_to_color <- function(val) {
  if (!is.finite(val)) return("#EAECCC")
  v <- log10(pmax(val, 1e-1))
  v <- max(-1, min(1, v))
  breaks <- seq(-1, 1, length.out = 1000)
  idx <- findInterval(v, breaks, all.inside = TRUE)
  ratio_palette[pmax(1, pmin(999, idx))]
}
.text_contrast <- function(col_hex) {
  rgb <- grDevices::col2rgb(col_hex)/255
  lum <- 0.2126*rgb[1] + 0.7152*rgb[2] + 0.0722*rgb[3]
  if (lum < 0.6) "white" else "black"
}

# ============================================================
# FIXED FULL-WIDTH MARGIN STRIP DRAWER
# ============================================================
.draw_margin_strip <- function(col_bg) {
  
  usr  <- par("usr")
  xrng <- usr[2] - usr[1]
  yrng <- usr[4] - usr[3]
  
  # same vertical scaling as your working version
  y0_user <- usr[4] + 0.065 * yrng
  y1_user <- usr[4] + 0.165 * yrng
  
  # extend slightly past panel borders
  x0_user <- usr[1] - 0.004 * xrng
  x1_user <- usr[2] + 0.004 * xrng
  
  # draw rectangle in user coordinates
  xpd_old <- par(xpd = NA)
  rect(x0_user, y0_user, x1_user, y1_user,
       col = col_bg, border = NA)
  par(xpd = xpd_old)
}

# =========================
# Panel plotter
# =========================
.plot_panel <- function(sel, row_i, col_i) {
  
  dat_all <- .load_residuals(sel$gl, sel$fl)
  key <- paste0("year", sel$yr, "_", sel$el)
  
  # ======================
  # EMPTY PANEL HANDLING
  # ======================
  
  if (is.null(dat_all) || is.null(dat_all[[key]]) || nrow(dat_all[[key]]) == 0) {
    
    plot(0, 0, type = "n", xlab = "", ylab = "",
         xaxt = "n", yaxt = "n", bty = "o",
         xlim = c(-1, 1), ylim = c(-1, 1))
    
    # Letter
    mtext(.letter_for_index(sel$idx), side = 3, line = -1.0,
          adj = 0.05, cex = 0.6, font = 1.5)
    
    # Left column y-labels
    if (col_i == 1) {
      mtext("", side = 2, line = 2.1, cex = 0.7)
      mtext(paste0(.get_glacier_name(sel$gl), ", ", sel$yr),
            side = 2, line = 3.0, cex = 0.7)
    }
    
    # Elevation text line
    if (col_i == 1) mtext("elevation:",
                          side = 3, line = 1.0, adj = 0, cex = 0.6)
    mtext(sprintf("%d", sel$el), side = 3, line = 1.0,
          cex = 0.6, adj = 1)
    
    # Colored S-index strip (neutral) BEHIND text
    .draw_margin_strip(col_bg = "#EAECCC")
    
    # S index text
    if (col_i == 1) mtext("S index:", side = 3, line = 0.3,
                          cex = 0.6, adj = 0)
    mtext("NA", side = 3, line = 0.3, cex = 0.6, adj = 1)
    
    return(invisible())
  }
  
  # ======================
  # NORMAL PANEL CONTENT
  # ======================
  
  merged <- dat_all[[key]]
  
  subset_data_obs <- data.frame(
    date      = merged$date,
    anomaly   = merged$anomaly_obs,
    vel       = merged$vel_obs,
    vel_error = if ("vel_error_obs" %in% names(merged)) merged$vel_error_obs else if ("vel_error" %in% names(merged)) merged$vel_error else NA_real_,
    weight    = if ("weight" %in% names(merged)) merged$weight else NA_real_
  )
  subset_data_model <- data.frame(
    date        = merged$date,
    anomaly     = merged$anomaly_model,
    vel         = merged$vel_model,
    is_original = if ("is_original" %in% names(merged)) merged$is_original else NA
  )
  subset_residual_data <- data.frame(
    date             = merged$date,
    residual_anomaly = merged$residual_anomaly,
    residual         = merged$residual,
    weight           = if ("weight" %in% names(merged)) merged$weight else NA_real_
  )
  
  spline_obs      <- .fit_spline(subset_data_obs$date,      subset_data_obs$anomaly, df = 20, spar = 0.7)
  spline_model    <- .fit_spline(subset_data_model$date,    subset_data_model$anomaly, df = 50, spar = NULL)
  spline_residual <- .fit_spline(subset_residual_data$date, subset_residual_data$residual_anomaly, df = 20, spar = 0.7)
  
  all_anoms <- c(subset_data_obs$anomaly, subset_data_model$anomaly, subset_residual_data$residual_anomaly)
  all_anoms <- all_anoms[is.finite(all_anoms)]
  y_max <- if (length(all_anoms)) max(abs(all_anoms), na.rm = TRUE) else 1
  y_min <- -y_max; y_max <- y_max
  
  # Base plot
  plot(subset_data_obs$date, subset_data_obs$anomaly, xlab = "", ylab = "",
       main = "", ylim = c(y_min, y_max),
       xlim = as.Date(c(paste0(sel$yr, "-01-01"), paste0(sel$yr, "-12-31"))),
       yaxt = "n", xaxt = "n",
       pch = 1, col = grDevices::adjustcolor("black", alpha.f = 0.5),
       cex = 0.5)
  
  # Left axis + glacier/year label
  if (col_i == 1) {
    axis(2, at = c(y_min, 0, y_max),
         labels = c("-1.0", "0.0", "1.0"), las = 1)
    mtext("", side = 2, line = 2.1, cex = 0.7)
    mtext(paste0(.get_glacier_name(sel$gl), ", ", sel$yr),
          side = 2, line = 3.0, cex = 0.7)
  }
  
  # Months only bottom row
  if (row_i == n_rows) {
    axis.Date(1,
              at = seq.Date(as.Date(paste0(sel$yr, "-03-01")),
                            as.Date(paste0(sel$yr, "-11-30")),
                            by = "4 months"),
              format = "%b", cex.axis = 1, las = 1)
  }
  
  # obs band
  if (!is.null(spline_obs) && nrow(spline_obs) > 0) {
    so <- subset(spline_obs, date >= as.Date(paste0(sel$yr, "-01-01")) &
                   date <= as.Date(paste0(sel$yr, "-12-31")))
    if (nrow(so) > 0) {
      polygon(c(so$date, rev(so$date)),
              c(so$upper, rev(so$lower)),
              col = grDevices::adjustcolor("black", alpha.f = 0.15),
              border = NA)
      lines(so$date, so$y, col = "black", lwd = 1)
    }
  }
  
  # Model
  original_points     <- subset(subset_data_model, is_original %in% TRUE)
  interpolated_points <- subset(subset_data_model, is_original %in% FALSE)
  
  if (nrow(interpolated_points) > 0) {
    points(interpolated_points$date, interpolated_points$anomaly,
           pch = 1, cex = 0.4,
           col = grDevices::adjustcolor("#ffb09c", alpha.f = 0.5))
  }
  
  if (!is.null(spline_model) && nrow(spline_model) > 0) {
    sm <- subset(spline_model, date >= as.Date(paste0(sel$yr, "-01-01")) &
                   date <= as.Date(paste0(sel$yr, "-12-31")))
    if (nrow(sm) > 0) {
      lines(sm$date, sm$y,
            col = grDevices::adjustcolor("red", alpha.f = 0.0),
            lwd = 1)
    }
  }
  
  if (nrow(original_points) > 0) {
    points(original_points$date, original_points$anomaly,
           pch = 1, col = grDevices::adjustcolor("red", alpha.f = 1.0),
           cex = 0.6)
  }
  
  # Residuals
  if (sum(is.finite(subset_residual_data$residual_anomaly)) > 0) {
    points(subset_residual_data$date,
           subset_residual_data$residual_anomaly,
           pch = 1,
           col = grDevices::adjustcolor("blue", alpha.f = 0.5),
           cex = 0.5)
  }
  
  if (!is.null(spline_residual) && nrow(spline_residual) > 0) {
    sr <- subset(spline_residual, date >= as.Date(paste0(sel$yr, "-01-01")) &
                   date <= as.Date(paste0(sel$yr, "-12-31")))
    if (nrow(sr) > 0) {
      polygon(c(sr$date, rev(sr$date)),
              c(sr$upper, rev(sr$lower)),
              col = grDevices::adjustcolor("blue", alpha.f = 0.15),
              border = NA)
      lines(sr$date, sr$y, col = "blue")
      abline(h = 0, lty = 3, lwd = 1, col = "black")
    }
  }
  
  # Letter
  mtext(.letter_for_index(sel$idx),
        side = 3, line = -1.0,
        adj = 0.05, cex = 0.6, font = 1.5)
  
  # ========================================================
  #  TOP MARGIN LABELS: ELEVATION AND S-INDEX WITH STRIP
  # ========================================================
  
  # --- elevation line ---
  if (col_i == 1) {
    mtext("elevation:",
          side = 3, line = 1.0,
          adj = 0, cex = 0.6)
  }
  mtext(sprintf("%d", sel$el),
        side = 3, line = 1.0,
        adj = 1, cex = 0.6)
  
  # --- S index value ---
  r_val <- .get_ratio_value_from_csv(sel$gl, sel$fl, sel$yr, sel$el)
  if (!is.finite(r_val)) {
    so <- if (!is.null(spline_obs)) subset(spline_obs, date >= as.Date(paste0(sel$yr, "-01-01")) &
                                             date <= as.Date(paste0(sel$yr, "-12-31"))) else NULL
    sm <- if (!is.null(spline_model)) subset(spline_model, date >= as.Date(paste0(sel$yr, "-01-01")) &
                                               date <= as.Date(paste0(sel$yr, "-12-31"))) else NULL
    sr <- if (!is.null(spline_residual)) subset(spline_residual, date >= as.Date(paste0(sel$yr, "-01-01")) &
                                                  date <= as.Date(paste0(sel$yr, "-12-31"))) else NULL
    r_val <- .compute_ratio_from_splines(so, sm, sr, sel$yr)
  }
  lab <- if (is.finite(r_val)) formatC(r_val, digits = 2, format = "f") else "NA"
  col_bg <- .map_ratio_to_color(r_val)
  text_col <- .text_contrast(col_bg)
  
  # --- full-width colored strip BEHIND S-index text ---
  .draw_margin_strip(col_bg = col_bg)
  
  # --- S index text ---
  if (col_i == 1) {
    mtext("S index:", side = 3, line = 0.3,
          adj = 0, cex = 0.6, col = text_col)
  }
  mtext(lab, side = 3, line = 0.3,
        adj = 1, cex = 0.6, col = text_col)
}

# =========================
# DRAW
# =========================
pdf(out_pdf, width = pdf_width, height = pdf_height)

par(mfrow = c(n_rows, n_cols),
    mar = c(0.2, 0.5, 2.0, 0.2),   # top margin ↑ enough for 2 lines
    oma = c(2, 4, 1, 0))

for (k in seq_along(selections)) {
  sel <- selections[[k]]
  r <- ((sel$idx - 1) %/% n_cols) + 1
  c <- ((sel$idx - 1) %%  n_cols) + 1
  .plot_panel(sel, row_i = r, col_i = c)
}

mtext("", side = 2, line = 3.8, outer = TRUE, cex = 1.2)
mtext("",     side = 1, line = 3.2, outer = TRUE, cex = 1.2)

dev.off()
cat("Saved:", out_pdf, "\n")