#' Interval Correlation Optimized Shifting (icoshift)
#'
#' R implementation of the icoshift alignment algorithm used in the original
#' MATLAB pipeline. Supports median/average targets and whole-spectrum alignment.
#'
#' @param target Reference mode: "median", "average", or a numeric row vector.
#' @param spectra Sample matrix (samples x variables).
#' @param mode Alignment mode: "whole" or integer number of intervals.
#' @param max_shift Maximum allowed shift in data points, or "f" for fast auto-search.
#' @return List with aligned spectra, interval info, shift indices, and target.
icoshift <- function(
    target = "median",
    spectra,
    mode = "whole",
    max_shift = 34) {

  spectra <- as.matrix(spectra)
  target_mode <- target

  if (is.character(target)) {
    if (target == "median") {
      ref <- matrix(apply(spectra, 2, stats::median, na.rm = TRUE), nrow = 1)
    } else if (target == "average") {
      ref <- matrix(colMeans(spectra, na.rm = TRUE), nrow = 1)
    } else {
      stop('target must be "median", "average", or a numeric vector')
    }
  } else {
    ref <- matrix(as.numeric(target), nrow = 1)
    target_mode <- "custom"
  }

  np <- nrow(spectra)
  mp <- ncol(spectra)

  if (identical(mode, "whole")) {
    intervals <- matrix(c(1L, mp), nrow = 1)
    whole <- TRUE
  } else if (is.numeric(mode) && length(mode) == 1L) {
    intervals <- build_intervals(mp, as.integer(mode))
    whole <- FALSE
  } else {
    stop('mode must be "whole" or an integer number of intervals')
  }

  n_intervals <- nrow(intervals)
  shift_matrix <- matrix(0, nrow = np, ncol = n_intervals)
  aligned <- spectra

  for (i in seq_len(n_intervals)) {
    start <- intervals[i, 1L]
    end <- intervals[i, 2L]
    seg_target <- ref[, start:end, drop = FALSE]
    seg_data <- aligned[, start:end, drop = FALSE]

    shifted <- coshift_segment(
      target = seg_target,
      spectra = seg_data,
      max_shift = max_shift,
      fill_previous = TRUE
    )

    aligned[, start:end] <- shifted$spectra
    shift_matrix[, i] <- shifted$shifts

    if (!whole && target_mode == "median") {
      ref[, start:end] <- matrix(apply(shifted$spectra, 2, stats::median, na.rm = TRUE), nrow = 1)
    }
  }

  ints <- cbind(seq_len(n_intervals), intervals, intervals[, 2] - intervals[, 1] + 1L)

  list(
    spectra = aligned,
    intervals = ints,
    shifts = shift_matrix,
    target = ref
  )
}

build_intervals <- function(n_points, n_intervals) {
  step <- n_points %/% n_intervals
  remain <- n_points %% n_intervals
  bounds <- integer()
  pos <- 1L
  for (i in seq_len(n_intervals)) {
    extra <- if (i <= remain) 1L else 0L
    end <- pos + step - 1L + extra
    if (i == n_intervals) end <- n_points
    bounds <- c(bounds, pos, end)
    pos <- end + 1L
  }
  matrix(bounds, ncol = 2, byrow = TRUE)
}

coshift_segment <- function(target, spectra, max_shift = 34, fill_previous = TRUE) {
  spectra <- as.matrix(spectra)
  target <- as.numeric(target)
  np <- nrow(spectra)
  mp <- ncol(spectra)

  if (identical(max_shift, "f")) {
    max_shift <- mp - 1L
    step <- max(1L, round(mp / 2) - 1L)
    repeat {
      result <- coshift_fixed(target, spectra, max_shift, fill_previous)
      if (max(abs(result$shifts)) < max_shift || max_shift + step >= mp) break
      max_shift <- max_shift + step
    }
    result
  } else if (identical(max_shift, "b")) {
    max_shift <- max(10L, round(0.05 * mp))
    step <- max(1L, round(0.05 * mp))
    repeat {
      result <- coshift_fixed(target, spectra, max_shift, fill_previous)
      if (max(abs(result$shifts)) < max_shift || max_shift + step >= mp) break
      max_shift <- max_shift + step
    }
    result
  } else {
    coshift_fixed(target, spectra, as.integer(max_shift), fill_previous)
  }
}

coshift_fixed <- function(target, spectra, max_shift, fill_previous = TRUE) {
  np <- nrow(spectra)
  mp <- ncol(spectra)
  shifts <- integer(np)
  aligned <- matrix(0, nrow = np, ncol = mp)

  for (i in seq_len(np)) {
    sig <- spectra[i, ]
    shift <- best_fft_shift(target, sig, max_shift)
    shifts[i] <- shift
    aligned[i, ] <- shift_vector(sig, shift, fill_previous)
  }

  list(spectra = aligned, shifts = shifts)
}

best_fft_shift <- function(target, signal, max_shift) {
  cc <- stats::ccf(
    target,
    signal,
    lag.max = max_shift,
    plot = FALSE,
    type = "correlation",
    na.action = stats::na.pass
  )
  best_lag <- as.integer(cc$lag[which.max(cc$acf)])
  as.integer(-best_lag)
}

shift_vector <- function(x, shift, fill_previous = TRUE) {
  n <- length(x)
  if (shift == 0L) return(x)
  if (shift > 0L) {
    if (fill_previous) {
      c(rep(x[1], shift), x[seq_len(n - shift)])
    } else {
      c(rep(NA_real_, shift), x[seq_len(n - shift)])
    }
  } else {
    s <- abs(shift)
    if (fill_previous) {
      c(x[(s + 1L):n], rep(x[n], s))
    } else {
      c(x[(s + 1L):n], rep(NA_real_, s))
    }
  }
}

#' Row-wise SNV-style standardization used in Analysis.m
snv_rows <- function(x) {
  x <- as.matrix(x)
  mu <- rowMeans(x, na.rm = TRUE)
  sd <- apply(x, 1, stats::sd, na.rm = TRUE)
  sd[sd == 0] <- 1
  sweep(sweep(x, 1, mu, "-"), 1, sd, "/")
}
