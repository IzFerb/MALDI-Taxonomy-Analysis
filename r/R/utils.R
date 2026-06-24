# Shared utilities for the MALDI gum taxonomy pipeline.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

message_step <- function(step, text) {
  message(sprintf("[%d/%d] %s", step$current, step$total, text))
  step$current <- step$current + 1L
  invisible(step)
}

new_step_counter <- function(total) {
  list(current = 1L, total = as.integer(total))
}

read_spectrum_txt <- function(path) {
  dat <- utils::read.table(
    path,
    header = FALSE,
    sep = "",
    dec = ".",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  if (ncol(dat) < 2L) {
    stop("Spectrum file must contain at least two columns (m/z and intensity): ", path)
  }
  list(mz = dat[[1]], intensity = dat[[2]])
}

parse_class_from_filename <- function(filename, anchor = "Acacia") {
  pos <- regexpr(anchor, filename, fixed = TRUE)
  if (pos[1] == -1L) return(NA_character_)
  k <- pos[1]
  parts <- gregexpr("_", filename, fixed = TRUE)[[1]]
  parts <- parts[parts > k]
  if (length(parts) < 2L) return(NA_character_)
  substr(filename, parts[1] + 1L, parts[2] - 1L)
}
