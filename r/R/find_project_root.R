#' Locate the project root (folder containing r/R/)
find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = FALSE)
  for (i in 0:6) {
    cand <- if (i == 0) {
      path
    } else {
      normalizePath(file.path(path, paste(rep("..", i), collapse = "/")), mustWork = FALSE)
    }
    if (dir.exists(file.path(cand, "r", "R"))) return(cand)
  }
  stop(
    "Could not find the project root (folder containing r/R/).\n",
    "In RStudio: Session → Set Working Directory → Choose Directory → select art_institute_code"
  )
}
