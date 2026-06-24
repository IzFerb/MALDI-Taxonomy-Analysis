#' Load MALDI spectra from folder structure
#'
#' Mirrors the MATLAB Analysis.m directory scan: each subfolder is treated as
#' one sample group. Blank .txt files in the project root are loaded separately.
#'
#' @param data_dir Root directory containing subfolders with .txt spectra.
#' @param exclude_dirs Directory names to skip (e.g. resources, .git).
#' @param blank_pattern Regex for blank files in data_dir root.
#' @return List with spectra matrix metadata.
load_spectra <- function(
    data_dir = ".",
    exclude_dirs = c("resources", "icoshift3p0", ".git", "r", "app", "output"),
    blank_pattern = "^Blank") {

  data_dir <- normalizePath(data_dir, mustWork = TRUE)

  entries <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
  entries <- sort(entries[!entries %in% exclude_dirs])

  filenames <- character()
  group_ids <- integer()
  group_names <- character()
  group_id <- 0L

  for (folder in entries) {
    folder_path <- file.path(data_dir, folder)
    txt_files <- sort(list.files(folder_path, pattern = "\\.txt$", ignore.case = TRUE))
    if (length(txt_files) == 0L) next
    group_id <- group_id + 1L
    group_names <- c(group_names, folder)
    filenames <- c(filenames, file.path(folder, txt_files))
    group_ids <- c(group_ids, rep(group_id, length(txt_files)))
  }

  blank_files <- list.files(
    data_dir,
    pattern = paste0(blank_pattern, ".*\\.txt$"),
    ignore.case = TRUE,
    full.names = FALSE
  )

  list(
    data_dir = data_dir,
    filenames = filenames,
    group_ids = group_ids,
    group_names = group_names,
    blank_files = blank_files
  )
}

#' Interpolate spectra onto a common m/z grid
#'
#' @param data_dir Project directory.
#' @param filenames Character vector of relative file paths.
#' @param mz_grid Numeric m/z grid (default matches MATLAB: 800:0.01:1300).
interpolate_spectra <- function(data_dir, filenames, mz_grid = seq(800, 1300, by = 0.01)) {
  n <- length(filenames)
  mat <- matrix(0, nrow = n, ncol = length(mz_grid))
  for (i in seq_len(n)) {
    spec <- read_spectrum_txt(file.path(data_dir, filenames[i]))
    mat[i, ] <- stats::spline(spec$mz, spec$intensity, xout = mz_grid)$y
  }
  list(
    mz = mz_grid,
    matrix = mat,
    filenames = filenames
  )
}
