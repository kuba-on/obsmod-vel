########################################
#  Ratio statistics for friction-law ratio matrices
#  - Reads all CSV matrices per friction-law folder (rows=years, cols=elevations).
#  - Skips any files whose basename starts with gl_103 or gl_104.
#  - Excludes the summary row named "mean_2016_21" from all statistics.
#  - Computes statistics on finite, non-NA values (zeros/negatives allowed):
#      0)  num_files
#      1)  num_results
#      2)  below_1
#      3)  below_1_%
#      4)  over_1
#      5)  over_1_%
#      6)  mixed                     (>= 0.1 & <= 10)
#      7)  mixed_%
#      8)  below_01                  (< 0.1)
#      9)  below_01_%
#     10)  over_10                   (>= 10)
#     11)  over_10_%
#     12)  max_elev_mixed            (highest elevation with >=0.1 & <=10) [not in per-elev tables]
#     13)  max_elev_front            (highest elevation with >=1)          [not in per-elev tables]
#     14)  fric_group                (triples with all three frics present)
#     15)  num_change                (triples with cross-1 inconsistency)
#     16)  num_change_%
#     17)  max_value
#     18)  min_value
#     19)  median_value
#
#  - Saves:
#    * Global (all vs each fric): ratio_stats_all.csv
#    * Elevations: ratio_stats_elevs_{fric1,fric2,fric3,all}.csv
#    * Years:      ratio_stats_years_{fric1,fric2,fric3,all}.csv
#    * Flowlines:  ratio_stats_flowlines_{fric1,fric2,fric3,all}.csv
#    * FRONTS (100–500 m only):
#         - Global:     ratio_stats_fronts_all.csv
#         - Flowlines:  ratio_stats_fronts_flowlines_{fric1,fric2,fric3,all}.csv
########################################

# ---------- Directories ----------
dir_fric1 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric1/Ratios"
dir_fric2 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric2/Ratios"
dir_fric3 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric3/Ratios"

out_stats <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Outputs/fric_comp/Ratio Statistics"
suppressWarnings(dir.create(out_stats, recursive = TRUE, showWarnings = FALSE))

# ---------- Small helpers ----------
`%||%` <- function(x, y) if (!is.null(x)) x else y

count_files_with_valid_values <- function(mat_list) {
  if (length(mat_list) == 0) return(0L)
  sum(vapply(mat_list, function(m) {
    v <- as.numeric(m)
    any(is.finite(v) & !is.na(v))
  }, logical(1)))
}

count_files_with_valid_values_in_elev_range <- function(mat_list, emin = 100, emax = 500) {
  if (length(mat_list) == 0) return(0L)
  sum(vapply(mat_list, function(m) {
    if (is.null(m)) return(FALSE)
    if (!is.null(rownames(m))) m <- m[rownames(m) != "mean_2016_21", , drop = FALSE]
    elevs <- suppressWarnings(as.numeric(colnames(m)))
    keep_cols <- is.finite(elevs) & elevs >= emin & elevs <= emax
    if (!any(keep_cols)) return(FALSE)
    v <- as.numeric(m[, keep_cols, drop = FALSE])
    any(is.finite(v) & !is.na(v))
  }, logical(1)))
}

# ---------- I/O helpers ----------
read_ratio_folder <- function(folder_path, friction_tag) {
  patt <- paste0("^gl_\\d+_\\d+_ratio_", friction_tag, "\\.csv$")
  files <- list.files(folder_path, pattern = patt, full.names = TRUE)
  
  base_names <- basename(files)
  keep <- !grepl("^gl_(103|104)_", base_names)
  files <- files[keep]
  base_names <- base_names[keep]
  
  out <- list()
  for (i in seq_along(files)) {
    f <- files[i]
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
    
    if (!is.null(rownames(df))) {
      df <- df[rownames(df) != "mean_2016_21", , drop = FALSE]
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
    
    base <- base_names[i]
    id <- sub(paste0("_ratio_", friction_tag, "\\.csv$"), "", base)
    if (grepl("^gl_(103|104)_", id)) next
    
    out[[id]] <- m
  }
  out
}

# ---------- Flatteners ----------
flatten_vals_elevs <- function(mat_list) {
  vals_out <- numeric(0)
  elevs_out <- numeric(0)
  if (length(mat_list) == 0) return(list(vals = vals_out, elevs = elevs_out))
  for (m in mat_list) {
    if (is.null(m)) next
    if (!is.null(rownames(m))) m <- m[rownames(m) != "mean_2016_21", , drop = FALSE]
    vals <- as.numeric(m)
    elevs <- suppressWarnings(as.numeric(colnames(m)))
    elevs_rep <- rep(elevs, each = nrow(m))
    ok <- is.finite(vals) & !is.na(vals) & is.finite(elevs_rep) & !is.na(elevs_rep)
    vals_out  <- c(vals_out,  vals[ok])
    elevs_out <- c(elevs_out, elevs_rep[ok])
  }
  list(vals = vals_out, elevs = elevs_out)
}

flatten_vals_elevs_range <- function(mat_list, emin = 100, emax = 500) {
  vals_out <- numeric(0)
  elevs_out <- numeric(0)
  if (length(mat_list) == 0) return(list(vals = vals_out, elevs = elevs_out))
  for (m in mat_list) {
    if (is.null(m)) next
    if (!is.null(rownames(m))) m <- m[rownames(m) != "mean_2016_21", , drop = FALSE]
    elevs <- suppressWarnings(as.numeric(colnames(m)))
    keep_cols <- is.finite(elevs) & elevs >= emin & elevs <= emax
    if (!any(keep_cols)) next
    vals <- as.numeric(m[, keep_cols, drop = FALSE])
    elevs_rep <- rep(elevs[keep_cols], each = nrow(m))
    ok <- is.finite(vals) & !is.na(vals) & is.finite(elevs_rep) & !is.na(elevs_rep)
    vals_out  <- c(vals_out,  vals[ok])
    elevs_out <- c(elevs_out, elevs_rep[ok])
  }
  list(vals = vals_out, elevs = elevs_out)
}

# ---------- Collectors (per-year / per-elevation) ----------
collect_by_year_detailed <- function(mat_list) {
  out <- list()
  if (length(mat_list) == 0) return(out)
  years <- unique(as.character(unlist(lapply(mat_list, function(m) rownames(m)))))
  years <- setdiff(years, "mean_2016_21")
  years <- as.character(sort(suppressWarnings(as.numeric(years[is.finite(suppressWarnings(as.numeric(years)))]))))
  for (yr in years) {
    vals <- numeric(0); elevs <- numeric(0); n_files <- 0L
    for (fid in names(mat_list)) {
      m <- mat_list[[fid]]
      if (!is.null(rownames(m)) && yr %in% rownames(m)) {
        row_sel <- rownames(m) == yr
        v_row <- as.numeric(m[row_sel, , drop = TRUE])
        e_row <- suppressWarnings(as.numeric(colnames(m)))
        ok <- is.finite(v_row) & !is.na(v_row) & is.finite(e_row) & !is.na(e_row)
        if (any(ok)) { vals <- c(vals, v_row[ok]); elevs <- c(elevs, e_row[ok]); n_files <- n_files + 1L }
      }
    }
    out[[yr]] <- list(vals = vals, elevs = elevs, n_files = n_files)
  }
  out
}

collect_by_elev_detailed <- function(mat_list) {
  out <- list()
  if (length(mat_list) == 0) return(out)
  elevs <- unique(as.character(unlist(lapply(mat_list, function(m) colnames(m)))))
  elevs_num <- suppressWarnings(as.numeric(elevs))
  elevs <- as.character(sort(elevs_num[is.finite(elevs_num)]))
  for (el in elevs) {
    vals <- numeric(0); n_files <- 0L
    for (fid in names(mat_list)) {
      m <- mat_list[[fid]]
      if (!is.null(colnames(m)) && el %in% colnames(m)) {
        if (!is.null(rownames(m))) m <- m[rownames(m) != "mean_2016_21", , drop = FALSE]
        v <- as.numeric(m[, el, drop = TRUE]); v <- v[is.finite(v) & !is.na(v)]
        if (length(v) > 0) { vals <- c(vals, v); n_files <- n_files + 1L }
      }
    }
    out[[el]] <- list(vals = vals, n_files = n_files)
  }
  out
}

# ---------- Cross-friction change counts ----------
# Count groups/changes for triples (f1,f2,f3) at same (flowline, year, elevation).
# Optional filters: years (character), elevs (numeric or character), fids (character), elevation range (emin/emax).
cross_fric_change_counts <- function(f1_list, f2_list, f3_list,
                                     years = NULL, elevs = NULL, fids = NULL,
                                     emin = NULL, emax = NULL) {
  groups <- 0L
  changes <- 0L
  # Which flowlines to consider
  fids_all <- Reduce(intersect, list(names(f1_list), names(f2_list), names(f3_list)))
  if (!is.null(fids)) fids_all <- intersect(fids_all, fids)
  if (length(fids_all) == 0) return(list(groups = groups, changes = changes))
  
  for (fid in fids_all) {
    m1 <- f1_list[[fid]]; m2 <- f2_list[[fid]]; m3 <- f3_list[[fid]]
    # common rows/cols
    rows <- Reduce(intersect, list(rownames(m1), rownames(m2), rownames(m3)))
    rows <- setdiff(rows, "mean_2016_21")
    cols <- Reduce(intersect, list(colnames(m1), colnames(m2), colnames(m3)))
    
    # filter by years
    if (!is.null(years)) rows <- intersect(rows, years)
    # filter by elev values or range
    if (!is.null(elevs)) {
      cols <- intersect(cols, as.character(suppressWarnings(as.numeric(elevs))))
    }
    if (!is.null(emin) || !is.null(emax)) {
      en <- suppressWarnings(as.numeric(cols))
      keep <- is.finite(en)
      if (!is.null(emin)) keep <- keep & (en >= emin)
      if (!is.null(emax)) keep <- keep & (en <= emax)
      cols <- as.character(en[keep])
    }
    
    if (length(rows) == 0 || length(cols) == 0) next
    
    for (yr in rows) {
      for (el in cols) {
        v1 <- suppressWarnings(as.numeric(m1[yr, el]))
        v2 <- suppressWarnings(as.numeric(m2[yr, el]))
        v3 <- suppressWarnings(as.numeric(m3[yr, el]))
        if (all(is.finite(c(v1, v2, v3)) & !is.na(c(v1, v2, v3)))) {
          groups <- groups + 1L
          sides <- c(v1 < 1, v2 < 1, v3 < 1)  # TRUE = below, FALSE = >= 1
          if (!(all(sides) || all(!sides))) changes <- changes + 1L
        }
      }
    }
  }
  list(groups = groups, changes = changes)
}

# ---------- Statistics labels & core ----------
stat_labels <- c(
  "num_files",
  "num_results",
  "below_1",
  "below_1_%",
  "over_1",
  "over_1_%",
  "mixed",
  "mixed_%",
  "below_01",
  "below_01_%",
  "over_10",
  "over_10_%",
  "max_elev_mixed",
  "max_elev_front",
  "fric_group",
  "num_change",
  "num_change_%",
  "max_value",
  "min_value",
  "median_value"
)

compute_stats_vec <- function(x, n_files = NA_integer_, elevs = NULL, include_elev_metrics = TRUE,
                              fric_group = NA_integer_, num_change = NA_integer_) {
  ok <- is.finite(x) & !is.na(x)
  x <- x[ok]
  if (!is.null(elevs)) elevs <- elevs[ok]
  n <- length(x)
  
  c_below_1   <- sum(x < 1,  na.rm = TRUE)
  c_overeq_1  <- sum(x >= 1, na.rm = TRUE)
  c_below_01  <- sum(x < 0.1, na.rm = TRUE)
  c_overeq_10 <- sum(x >= 10,  na.rm = TRUE)
  c_mixed     <- sum(x >= 0.1 & x <= 10, na.rm = TRUE)
  
  pct <- function(count, total) if (total > 0) round(100 * count / total, 2) else NA_real_
  
  if (include_elev_metrics && !is.null(elevs) && length(elevs) == length(x) && length(x) > 0) {
    mix_idx <- which(x >= 0.1 & x <= 10)
    front_idx <- which(x >= 1)
    max_elev_mixed <- if (length(mix_idx) > 0) max(elevs[mix_idx], na.rm = TRUE) else NA_real_
    max_elev_front <- if (length(front_idx) > 0) max(elevs[front_idx], na.rm = TRUE) else NA_real_
  } else {
    max_elev_mixed <- NA_real_
    max_elev_front <- NA_real_
  }
  
  out <- c(
    num_files       = as.numeric(n_files),
    num_results     = n,
    below_1         = c_below_1,
    `below_1_%`     = pct(c_below_1, n),
    over_1          = c_overeq_1,
    `over_1_%`      = pct(c_overeq_1, n),
    mixed           = c_mixed,
    `mixed_%`       = pct(c_mixed, n),
    below_01        = c_below_01,
    `below_01_%`    = pct(c_below_01, n),
    over_10         = c_overeq_10,
    `over_10_%`     = pct(c_overeq_10, n),
    max_elev_mixed  = max_elev_mixed,
    max_elev_front  = max_elev_front,
    fric_group      = as.numeric(fric_group),
    num_change      = as.numeric(num_change),
    `num_change_%`  = if (is.finite(fric_group) && fric_group > 0) round(100 * num_change / fric_group, 2) else NA_real_,
    max_value       = if (n > 0) max(x, na.rm = TRUE) else NA_real_,
    min_value       = if (n > 0) min(x, na.rm = TRUE) else NA_real_,
    median_value    = if (n > 0) median(x, na.rm = TRUE) else NA_real_
  )
  out[stat_labels]
}

build_stats_matrix_detailed <- function(detail_list, col_order, include_elev_metrics = TRUE) {
  if (length(col_order) == 0) {
    m <- matrix(numeric(0), nrow = length(stat_labels), ncol = 0)
    rownames(m) <- stat_labels
    return(m)
  }
  M <- matrix(NA_real_, nrow = length(stat_labels), ncol = length(col_order))
  rownames(M) <- stat_labels
  colnames(M) <- col_order
  for (nm in col_order) {
    info <- detail_list[[nm]]
    if (is.null(info)) {
      M[, nm] <- compute_stats_vec(numeric(0), n_files = 0L, elevs = NULL, include_elev_metrics = include_elev_metrics)
    } else {
      v  <- info$vals %||% numeric(0)
      e  <- info$elevs %||% NULL
      nf <- info$n_files %||% 0L
      fg <- info$fric_group %||% NA_integer_
      ch <- info$num_change %||% NA_integer_
      if (!include_elev_metrics) e <- NULL
      M[, nm] <- compute_stats_vec(v, n_files = nf, elevs = e, include_elev_metrics = include_elev_metrics,
                                   fric_group = fg, num_change = ch)
    }
  }
  M
}

# ---------- Read data ----------
fric1_list <- read_ratio_folder(dir_fric1, "fric1")
fric2_list <- read_ratio_folder(dir_fric2, "fric2")
fric3_list <- read_ratio_folder(dir_fric3, "fric3")

# ---------- 1) GLOBAL ----------
fe1 <- flatten_vals_elevs(fric1_list); fe2 <- flatten_vals_elevs(fric2_list); fe3 <- flatten_vals_elevs(fric3_list)
x1_all <- fe1$vals; e1_all <- fe1$elevs
x2_all <- fe2$vals; e2_all <- fe2$elevs
x3_all <- fe3$vals; e3_all <- fe3$elevs
x_all <- c(x1_all, x2_all, x3_all); e_all <- c(e1_all, e2_all, e3_all)

nfiles_f1  <- count_files_with_valid_values(fric1_list)
nfiles_f2  <- count_files_with_valid_values(fric2_list)
nfiles_f3  <- count_files_with_valid_values(fric3_list)
nfiles_all <- nfiles_f1 + nfiles_f2 + nfiles_f3

cc_global <- cross_fric_change_counts(fric1_list, fric2_list, fric3_list)

stats_all   <- compute_stats_vec(x_all,  n_files = nfiles_all, elevs = e_all, include_elev_metrics = TRUE,
                                 fric_group = cc_global$groups, num_change = cc_global$changes)
stats_fric1 <- compute_stats_vec(x1_all, n_files = nfiles_f1,  elevs = e1_all, include_elev_metrics = TRUE)
stats_fric2 <- compute_stats_vec(x2_all, n_files = nfiles_f2,  elevs = e2_all, include_elev_metrics = TRUE)
stats_fric3 <- compute_stats_vec(x3_all, n_files = nfiles_f3,  elevs = e3_all, include_elev_metrics = TRUE)

df_global <- data.frame(
  stat  = stat_labels,
  all   = as.numeric(stats_all),
  fric1 = as.numeric(stats_fric1),
  fric2 = as.numeric(stats_fric2),
  fric3 = as.numeric(stats_fric3),
  check.names = FALSE,
  row.names = NULL
)
write.csv(df_global, file = file.path(out_stats, "ratio_stats_all.csv"), row.names = FALSE)

# ---------- 2) BY ELEVATION (no max_elev_* here) ----------
f1_by_elev_det <- collect_by_elev_detailed(fric1_list)
f2_by_elev_det <- collect_by_elev_detailed(fric2_list)
f3_by_elev_det <- collect_by_elev_detailed(fric3_list)

elev_union <- sort(unique(as.character(c(names(f1_by_elev_det), names(f2_by_elev_det), names(f3_by_elev_det)))))
elev_union_num <- suppressWarnings(as.numeric(elev_union))
elev_union <- as.character(sort(elev_union_num[is.finite(elev_union_num)]))

# per-fric matrices (fric_group/num_change not applicable -> NA internally)
M_elev_f1  <- build_stats_matrix_detailed(f1_by_elev_det, elev_union, include_elev_metrics = FALSE)
M_elev_f2  <- build_stats_matrix_detailed(f2_by_elev_det, elev_union, include_elev_metrics = FALSE)
M_elev_f3  <- build_stats_matrix_detailed(f3_by_elev_det, elev_union, include_elev_metrics = FALSE)

# "All" by elevation: add cross-fric change counts per elevation
f_all_by_elev_det <- setNames(vector("list", length(elev_union)), elev_union)
for (el in elev_union) {
  v <- c((f1_by_elev_det[[el]]$vals %||% numeric(0)),
         (f2_by_elev_det[[el]]$vals %||% numeric(0)),
         (f3_by_elev_det[[el]]$vals %||% numeric(0)))
  nf <- (f1_by_elev_det[[el]]$n_files %||% 0L) +
    (f2_by_elev_det[[el]]$n_files %||% 0L) +
    (f3_by_elev_det[[el]]$n_files %||% 0L)
  cc <- cross_fric_change_counts(fric1_list, fric2_list, fric3_list, elevs = el)
  f_all_by_elev_det[[el]] <- list(vals = v, n_files = nf, fric_group = cc$groups, num_change = cc$changes)
}
M_elev_all <- build_stats_matrix_detailed(f_all_by_elev_det, elev_union, include_elev_metrics = FALSE)

write.csv(data.frame(stat = rownames(M_elev_f1), M_elev_f1, check.names = FALSE, row.names = NULL),
          file.path(out_stats, "ratio_stats_elevs_fric1.csv"), row.names = FALSE)
write.csv(data.frame(stat = rownames(M_elev_f2), M_elev_f2, check.names = FALSE, row.names = NULL),
          file.path(out_stats, "ratio_stats_elevs_fric2.csv"), row.names = FALSE)
write.csv(data.frame(stat = rownames(M_elev_f3), M_elev_f3, check.names = FALSE, row.names = NULL),
          file.path(out_stats, "ratio_stats_elevs_fric3.csv"), row.names = FALSE)
write.csv(data.frame(stat = rownames(M_elev_all), M_elev_all, check.names = FALSE, row.names = NULL),
          file.path(out_stats, "ratio_stats_elevs_all.csv"), row.names = FALSE)

# ---------- 3) BY YEAR ----------
f1_by_year_det <- collect_by_year_detailed(fric1_list)
f2_by_year_det <- collect_by_year_detailed(fric2_list)
f3_by_year_det <- collect_by_year_detailed(fric3_list)

year_union <- sort(unique(as.character(c(names(f1_by_year_det), names(f2_by_year_det), names(f3_by_year_det)))))
year_union_num <- suppressWarnings(as.numeric(year_union))
year_union <- as.character(sort(year_union_num[is.finite(year_union_num)]))

M_year_f1 <- build_stats_matrix_detailed(f1_by_year_det, year_union, include_elev_metrics = TRUE)
M_year_f2 <- build_stats_matrix_detailed(f2_by_year_det, year_union, include_elev_metrics = TRUE)
M_year_f3 <- build_stats_matrix_detailed(f3_by_year_det, year_union, include_elev_metrics = TRUE)

f_all_by_year_det <- setNames(vector("list", length(year_union)), year_union)
for (yr in year_union) {
  v <- c((f1_by_year_det[[yr]]$vals %||% numeric(0)),
         (f2_by_year_det[[yr]]$vals %||% numeric(0)),
         (f3_by_year_det[[yr]]$vals %||% numeric(0)))
  e <- c((f1_by_year_det[[yr]]$elevs %||% numeric(0)),
         (f2_by_year_det[[yr]]$elevs %||% numeric(0)),
         (f3_by_year_det[[yr]]$elevs %||% numeric(0)))
  nf <- (f1_by_year_det[[yr]]$n_files %||% 0L) +
    (f2_by_year_det[[yr]]$n_files %||% 0L) +
    (f3_by_year_det[[yr]]$n_files %||% 0L)
  cc <- cross_fric_change_counts(fric1_list, fric2_list, fric3_list, years = yr)
  f_all_by_year_det[[yr]] <- list(vals = v, elevs = e, n_files = nf, fric_group = cc$groups, num_change = cc$changes)
}
M_year_all <- build_stats_matrix_detailed(f_all_by_year_det, year_union, include_elev_metrics = TRUE)

write.csv(data.frame(stat = rownames(M_year_f1), M_year_f1, check.names = FALSE, row.names = NULL),
          file.path(out_stats, "ratio_stats_years_fric1.csv"), row.names = FALSE)
write.csv(data.frame(stat = rownames(M_year_f2), M_year_f2, check.names = FALSE, row.names = NULL),
          file.path(out_stats, "ratio_stats_years_fric2.csv"), row.names = FALSE)
write.csv(data.frame(stat = rownames(M_year_f3), M_year_f3, check.names = FALSE, row.names = NULL),
          file.path(out_stats, "ratio_stats_years_fric3.csv"), row.names = FALSE)
write.csv(data.frame(stat = rownames(M_year_all), M_year_all, check.names = FALSE, row.names = NULL),
          file.path(out_stats, "ratio_stats_years_all.csv"), row.names = FALSE)

# ---------- 4) BY FLOWLINE ----------
flow_ids <- sort(unique(c(names(fric1_list), names(fric2_list), names(fric3_list))))

flowline_stats_df <- function(flist) {
  if (length(flow_ids) == 0) return(data.frame())
  M <- matrix(NA_real_, nrow = length(flow_ids), ncol = length(stat_labels))
  colnames(M) <- stat_labels; rownames(M) <- flow_ids
  for (fid in flow_ids) {
    m <- flist[[fid]]
    if (!is.null(m) && !is.null(rownames(m))) m <- m[rownames(m) != "mean_2016_21", , drop = FALSE]
    if (!is.null(m)) {
      vals <- as.numeric(m)
      elevs <- suppressWarnings(as.numeric(colnames(m)))
      elevs_rep <- rep(elevs, each = nrow(m))
      ok <- is.finite(vals) & !is.na(vals) & is.finite(elevs_rep) & !is.na(elevs_rep)
      vals <- vals[ok]; elevs_rep <- elevs_rep[ok]
    } else {
      vals <- numeric(0); elevs_rep <- NULL
    }
    n_files <- if (!is.null(m) && length(vals) > 0) 1L else 0L
    M[fid, ] <- compute_stats_vec(vals, n_files = n_files, elevs = elevs_rep, include_elev_metrics = TRUE)
  }
  as.data.frame(M, check.names = FALSE)
}

df_flow_f1 <- flowline_stats_df(fric1_list)
df_flow_f2 <- flowline_stats_df(fric2_list)
df_flow_f3 <- flowline_stats_df(fric3_list)

# Flowlines "all": add cross-fric change counts per flowline
df_flow_all <- {
  M <- matrix(NA_real_, nrow = length(flow_ids), ncol = length(stat_labels))
  colnames(M) <- stat_labels; rownames(M) <- flow_ids
  for (fid in flow_ids) {
    m1 <- fric1_list[[fid]]; m2 <- fric2_list[[fid]]; m3 <- fric3_list[[fid]]
    if (!is.null(m1) && !is.null(rownames(m1))) m1 <- m1[rownames(m1) != "mean_2016_21", , drop = FALSE]
    if (!is.null(m2) && !is.null(rownames(m2))) m2 <- m2[rownames(m2) != "mean_2016_21", , drop = FALSE]
    if (!is.null(m3) && !is.null(rownames(m3))) m3 <- m3[rownames(m3) != "mean_2016_21", , drop = FALSE]
    
    add_me <- function(mat) {
      if (is.null(mat)) return(list(v = numeric(0), e = numeric(0), nf = 0L))
      vals <- as.numeric(mat)
      elevs <- suppressWarnings(as.numeric(colnames(mat)))
      elevs_rep <- rep(elevs, each = nrow(mat))
      ok <- is.finite(vals) & !is.na(vals) & is.finite(elevs_rep) & !is.na(elevs_rep)
      list(v = vals[ok], e = elevs_rep[ok], nf = if (any(ok)) 1L else 0L)
    }
    a1 <- add_me(m1); a2 <- add_me(m2); a3 <- add_me(m3)
    vals_all  <- c(a1$v, a2$v, a3$v)
    elevs_all <- c(a1$e, a2$e, a3$e)
    n_files   <- a1$nf + a2$nf + a3$nf
    
    cc <- cross_fric_change_counts(fric1_list, fric2_list, fric3_list, fids = fid)
    M[fid, ] <- compute_stats_vec(vals_all, n_files = n_files, elevs = elevs_all, include_elev_metrics = TRUE,
                                  fric_group = cc$groups, num_change = cc$changes)
  }
  as.data.frame(M, check.names = FALSE)
}

write.csv(data.frame(flowline = rownames(df_flow_f1), df_flow_f1, row.names = NULL, check.names = FALSE),
          file.path(out_stats, "ratio_stats_flowlines_fric1.csv"), row.names = FALSE)
write.csv(data.frame(flowline = rownames(df_flow_f2), df_flow_f2, row.names = NULL, check.names = FALSE),
          file.path(out_stats, "ratio_stats_flowlines_fric2.csv"), row.names = FALSE)
write.csv(data.frame(flowline = rownames(df_flow_f3), df_flow_f3, row.names = NULL, check.names = FALSE),
          file.path(out_stats, "ratio_stats_flowlines_fric3.csv"), row.names = FALSE)
write.csv(data.frame(flowline = rownames(df_flow_all), df_flow_all, row.names = NULL, check.names = FALSE),
          file.path(out_stats, "ratio_stats_flowlines_all.csv"), row.names = FALSE)

# ---------- 5) FRONTS (100–500 m) ----------
FRONT_MIN_E <- 100; FRONT_MAX_E <- 500

fe1_front <- flatten_vals_elevs_range(fric1_list, FRONT_MIN_E, FRONT_MAX_E)
fe2_front <- flatten_vals_elevs_range(fric2_list, FRONT_MIN_E, FRONT_MAX_E)
fe3_front <- flatten_vals_elevs_range(fric3_list, FRONT_MIN_E, FRONT_MAX_E)

x1_front <- fe1_front$vals; e1_front <- fe1_front$elevs
x2_front <- fe2_front$vals; e2_front <- fe2_front$elevs
x3_front <- fe3_front$vals; e3_front <- fe3_front$elevs
x_front_all <- c(x1_front, x2_front, x3_front)
e_front_all <- c(e1_front, e2_front, e3_front)

nf1_front <- count_files_with_valid_values_in_elev_range(fric1_list, FRONT_MIN_E, FRONT_MAX_E)
nf2_front <- count_files_with_valid_values_in_elev_range(fric2_list, FRONT_MIN_E, FRONT_MAX_E)
nf3_front <- count_files_with_valid_values_in_elev_range(fric3_list, FRONT_MIN_E, FRONT_MAX_E)
nf_all_front <- nf1_front + nf2_front + nf3_front

cc_front_global <- cross_fric_change_counts(fric1_list, fric2_list, fric3_list, emin = FRONT_MIN_E, emax = FRONT_MAX_E)

stats_front_all   <- compute_stats_vec(x_front_all, n_files = nf_all_front, elevs = e_front_all, include_elev_metrics = TRUE,
                                       fric_group = cc_front_global$groups, num_change = cc_front_global$changes)
stats_front_fric1 <- compute_stats_vec(x1_front,    n_files = nf1_front,    elevs = e1_front, include_elev_metrics = TRUE)
stats_front_fric2 <- compute_stats_vec(x2_front,    n_files = nf2_front,    elevs = e2_front, include_elev_metrics = TRUE)
stats_front_fric3 <- compute_stats_vec(x3_front,    n_files = nf3_front,    elevs = e3_front, include_elev_metrics = TRUE)

df_front_global <- data.frame(
  stat  = stat_labels,
  all   = as.numeric(stats_front_all),
  fric1 = as.numeric(stats_front_fric1),
  fric2 = as.numeric(stats_front_fric2),
  fric3 = as.numeric(stats_front_fric3),
  check.names = FALSE,
  row.names = NULL
)
write.csv(df_front_global, file = file.path(out_stats, "ratio_stats_fronts_all.csv"), row.names = FALSE)

# Flowline fronts per friction
flowline_fronts_df <- function(flist, emin = 100, emax = 500) {
  if (length(flow_ids) == 0) return(data.frame())
  M <- matrix(NA_real_, nrow = length(flow_ids), ncol = length(stat_labels))
  colnames(M) <- stat_labels; rownames(M) <- flow_ids
  for (fid in flow_ids) {
    m <- flist[[fid]]
    if (!is.null(m) && !is.null(rownames(m))) m <- m[rownames(m) != "mean_2016_21", , drop = FALSE]
    if (!is.null(m)) {
      elevs <- suppressWarnings(as.numeric(colnames(m)))
      keep_cols <- is.finite(elevs) & elevs >= emin & elevs <= emax
      if (any(keep_cols)) {
        vals <- as.numeric(m[, keep_cols, drop = FALSE])
        elevs_rep <- rep(elevs[keep_cols], each = nrow(m))
        ok <- is.finite(vals) & !is.na(vals) & is.finite(elevs_rep) & !is.na(elevs_rep)
        vals <- vals[ok]; elevs_rep <- elevs_rep[ok]
      } else { vals <- numeric(0); elevs_rep <- NULL }
    } else { vals <- numeric(0); elevs_rep <- NULL }
    n_files <- if (length(vals) > 0) 1L else 0L
    M[fid, ] <- compute_stats_vec(vals, n_files = n_files, elevs = elevs_rep, include_elev_metrics = TRUE)
  }
  as.data.frame(M, check.names = FALSE)
}

df_front_flow_f1 <- flowline_fronts_df(fric1_list, FRONT_MIN_E, FRONT_MAX_E)
df_front_flow_f2 <- flowline_fronts_df(fric2_list, FRONT_MIN_E, FRONT_MAX_E)
df_front_flow_f3 <- flowline_fronts_df(fric3_list, FRONT_MIN_E, FRONT_MAX_E)

# Flowline fronts "all" (combine frictions per flowline, 100–500 m) + cross-fric change counts
df_front_flow_all <- {
  M <- matrix(NA_real_, nrow = length(flow_ids), ncol = length(stat_labels))
  colnames(M) <- stat_labels; rownames(M) <- flow_ids
  for (fid in flow_ids) {
    m1 <- fric1_list[[fid]]; m2 <- fric2_list[[fid]]; m3 <- fric3_list[[fid]]
    if (!is.null(m1) && !is.null(rownames(m1))) m1 <- m1[rownames(m1) != "mean_2016_21", , drop = FALSE]
    if (!is.null(m2) && !is.null(rownames(m2))) m2 <- m2[rownames(m2) != "mean_2016_21", , drop = FALSE]
    if (!is.null(m3) && !is.null(rownames(m3))) m3 <- m3[rownames(m3) != "mean_2016_21", , drop = FALSE]
    
    pick_range <- function(mat, emin, emax) {
      if (is.null(mat)) return(list(v = numeric(0), e = numeric(0), nf = 0L))
      elevs <- suppressWarnings(as.numeric(colnames(mat)))
      keep_cols <- is.finite(elevs) & elevs >= emin & elevs <= emax
      if (!any(keep_cols)) return(list(v = numeric(0), e = numeric(0), nf = 0L))
      vals <- as.numeric(mat[, keep_cols, drop = FALSE])
      elevs_rep <- rep(elevs[keep_cols], each = nrow(mat))
      ok <- is.finite(vals) & !is.na(vals) & is.finite(elevs_rep) & !is.na(elevs_rep)
      list(v = vals[ok], e = elevs_rep[ok], nf = if (any(ok)) 1L else 0L)
    }
    
    a1 <- pick_range(m1, FRONT_MIN_E, FRONT_MAX_E)
    a2 <- pick_range(m2, FRONT_MIN_E, FRONT_MAX_E)
    a3 <- pick_range(m3, FRONT_MIN_E, FRONT_MAX_E)
    
    vals_all  <- c(a1$v, a2$v, a3$v)
    elevs_all <- c(a1$e, a2$e, a3$e)
    n_files   <- a1$nf + a2$nf + a3$nf
    
    cc <- cross_fric_change_counts(fric1_list, fric2_list, fric3_list, fids = fid,
                                   emin = FRONT_MIN_E, emax = FRONT_MAX_E)
    M[fid, ] <- compute_stats_vec(vals_all, n_files = n_files, elevs = elevs_all, include_elev_metrics = TRUE,
                                  fric_group = cc$groups, num_change = cc$changes)
  }
  as.data.frame(M, check.names = FALSE)
}

# Save FRONTS outputs
write.csv(data.frame(stat = stat_labels,
                     all   = df_front_global$all,
                     fric1 = df_front_global$fric1,
                     fric2 = df_front_global$fric2,
                     fric3 = df_front_global$fric3),
          file.path(out_stats, "ratio_stats_fronts_all.csv"), row.names = FALSE)

write.csv(data.frame(flowline = rownames(df_front_flow_f1), df_front_flow_f1, row.names = NULL, check.names = FALSE),
          file.path(out_stats, "ratio_stats_fronts_flowlines_fric1.csv"), row.names = FALSE)
write.csv(data.frame(flowline = rownames(df_front_flow_f2), df_front_flow_f2, row.names = NULL, check.names = FALSE),
          file.path(out_stats, "ratio_stats_fronts_flowlines_fric2.csv"), row.names = FALSE)
write.csv(data.frame(flowline = rownames(df_front_flow_f3), df_front_flow_f3, row.names = NULL, check.names = FALSE),
          file.path(out_stats, "ratio_stats_fronts_flowlines_fric3.csv"), row.names = FALSE)
write.csv(data.frame(flowline = rownames(df_front_flow_all), df_front_flow_all, row.names = NULL, check.names = FALSE),
          file.path(out_stats, "ratio_stats_fronts_flowlines_all.csv"), row.names = FALSE)

cat("Done. Ratio statistics saved to:\n",
    "- Global:     ", file.path(out_stats, "ratio_stats_all.csv"), "\n",
    "- Elevations: ", paste0(out_stats, c("/ratio_stats_elevs_fric1.csv",
                                          "/ratio_stats_elevs_fric2.csv",
                                          "/ratio_stats_elevs_fric3.csv",
                                          "/ratio_stats_elevs_all.csv")), "\n",
    "- Years:      ", paste0(out_stats, c("/ratio_stats_years_fric1.csv",
                                          "/ratio_stats_years_fric2.csv",
                                          "/ratio_stats_years_fric3.csv",
                                          "/ratio_stats_years_all.csv")), "\n",
    "- Flowlines:  ", paste0(out_stats, c("/ratio_stats_flowlines_fric1.csv",
                                          "/ratio_stats_flowlines_fric2.csv",
                                          "/ratio_stats_flowlines_fric3.csv",
                                          "/ratio_stats_flowlines_all.csv")), "\n",
    "- FRONTS:     ", paste0(out_stats, c("/ratio_stats_fronts_all.csv",
                                          "/ratio_stats_fronts_flowlines_fric1.csv",
                                          "/ratio_stats_fronts_flowlines_fric2.csv",
                                          "/ratio_stats_fronts_flowlines_fric3.csv",
                                          "/ratio_stats_fronts_flowlines_all.csv")), "\n",
    sep = "")