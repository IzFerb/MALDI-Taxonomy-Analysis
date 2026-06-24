#' Peak picking and alignment (MATLAB findpeaks + peak screening)
#'
#' @param x Sample matrix after preprocessing (samples x m/z points).
#' @param mz m/z axis values.
#' @param reference_idx Integer indices of reference samples for width calculation.
#' @param prominence_factor Multiplier for MinPeakProminence (default 10, as in MATLAB).
find_and_align_peaks <- function(
    x,
    mz,
    reference_idx,
    prominence_factor = 10) {

  x <- as.matrix(x)
  n_samples <- nrow(x)
  peaks <- vector("list", n_samples)

  for (i in seq_len(n_samples)) {
    y <- x[i, ]
    prom <- stats::sd(y, na.rm = TRUE) * prominence_factor
    peaks[[i]] <- findpeaks_with_props(y, mz, min_prominence = prom)
  }

  ref_widths <- unlist(lapply(peaks[reference_idx], function(p) p$width))
  ref_widths <- ref_widths[is.finite(ref_widths) & ref_widths > 0]
  if (length(ref_widths) == 0) {
    maximum_width <- 1
  } else if (length(ref_widths) == 1) {
    maximum_width <- ref_widths[1]
  } else {
    maximum_width <- mean(ref_widths, na.rm = TRUE) + stats::sd(ref_widths, na.rm = TRUE) * 3
  }
  if (!is.finite(maximum_width) || maximum_width <= 0) maximum_width <- 1

  ref_locs <- unique(unlist(lapply(peaks[reference_idx], function(p) p$location)))
  ref_locs <- sort(ref_locs)

  feature_names <- numeric()
  feature_matrix <- matrix(0, nrow = n_samples, ncol = 0)

  while (length(ref_locs) > 0L) {
    center <- ref_locs[1]
    in_window <- ref_locs >= center - maximum_width & ref_locs <= center + maximum_width
    if (!any(in_window, na.rm = TRUE)) {
      ref_locs <- ref_locs[-1]
      next
    }
    feature_center <- stats::median(ref_locs[in_window])
    ref_locs <- ref_locs[!in_window]

    col <- numeric(n_samples)
    for (i in seq_len(n_samples)) {
      p <- peaks[[i]]
      hits <- which(p$location > feature_center - maximum_width &
                      p$location < feature_center + maximum_width)
      if (length(hits) == 0L) {
        col[i] <- 0
      } else if (length(hits) == 1L) {
        col[i] <- p$height[hits]
      } else {
        col[i] <- max(p$height[hits])
      }
    }

    feature_names <- c(feature_names, feature_center)
    feature_matrix <- cbind(feature_matrix, col)
  }

  list(
    matrix = feature_matrix,
    peak_mz = feature_names,
    maximum_width = maximum_width
  )
}

findpeaks_with_props <- function(y, x, min_prominence, search_window = 500L) {
  n <- length(y)
  if (n < 3L) {
    return(list(location = numeric(), height = numeric(), width = numeric(), prominence = numeric()))
  }

  peak_idx <- which(y[2:(n - 1L)] > y[1:(n - 2L)] & y[2:(n - 1L)] > y[3:n]) + 1L
  if (length(peak_idx) == 0L) {
    return(list(location = numeric(), height = numeric(), width = numeric(), prominence = numeric()))
  }

  height_floor <- stats::median(y, na.rm = TRUE) + min_prominence / 2
  peak_idx <- peak_idx[y[peak_idx] >= height_floor]
  if (length(peak_idx) == 0L) {
    return(list(location = numeric(), height = numeric(), width = numeric(), prominence = numeric()))
  }

  heights <- y[peak_idx]
  prominences <- numeric(length(peak_idx))
  widths <- rep(NA_real_, length(peak_idx))

  for (j in seq_along(peak_idx)) {
    i <- peak_idx[j]
    left <- max(1L, i - search_window)
    right <- min(n, i + search_window)
    left_min <- min(y[left:i])
    right_min <- min(y[i:right])
    prominences[j] <- y[i] - max(left_min, right_min)
    if (prominences[j] < min_prominence) next

    half <- prominences[j] / 2
    li <- i
    while (li > left && (y[li] - min(y[left:li])) < half) li <- li - 1L
    ri <- i
    while (ri < right && (y[ri] - min(y[ri:right])) < half) ri <- ri + 1L
    widths[j] <- x[ri] - x[li]
  }

  keep <- prominences >= min_prominence
  list(
    location = x[peak_idx[keep]],
    height = heights[keep],
    width = widths[keep],
    prominence = prominences[keep]
  )
}

#' Remove peaks present in blank samples (PQN step in Analysis.m)
remove_blank_peaks <- function(feature_matrix, n_blanks = 4L) {
  if (n_blanks <= 0L || nrow(feature_matrix) < n_blanks) {
    return(list(matrix = feature_matrix, keep_idx = seq_len(ncol(feature_matrix))))
  }
  blank_rows <- feature_matrix[seq_len(n_blanks), , drop = FALSE]
  keep <- !apply(blank_rows, 2, function(col) any(col > 0))
  keep_idx <- which(keep)
  list(
    matrix = feature_matrix[, keep, drop = FALSE],
    keep_idx = keep_idx
  )
}
