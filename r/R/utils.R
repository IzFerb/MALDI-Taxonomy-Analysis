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

# Species label = the word right after "Acacia" (separated by "_" or a space),
# case-insensitive, returned in lowercase. E.g.
#   "SEN_A_Acacia_senegal_F_123.txt"      -> "senegal"
#   "SEY_C_Acacia seyal_K_78685.txt"      -> "seyal"
#   "STEN_B_Acacia_Stenocarpa_F_269.txt"  -> "stenocarpa"
parse_class_from_filename <- function(filename, anchor = "Acacia") {
  pattern <- paste0(anchor, "[_ ]+([^_ .]+)")
  m <- regmatches(filename, regexec(pattern, filename, ignore.case = TRUE))[[1]]
  if (length(m) < 2L) return(NA_character_)
  tolower(m[2])
}
