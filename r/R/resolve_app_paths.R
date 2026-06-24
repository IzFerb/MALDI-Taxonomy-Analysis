#' Resolve app, R code, and project root directories reliably
resolve_app_paths <- function() {
  # Method 1: --file= from commandArgs (works when Rscript/shiny sources app.R)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    app_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    return(list(
      app_dir = app_dir,
      r_code_dir = normalizePath(file.path(app_dir, "..", "R")),
      project_root = normalizePath(file.path(app_dir, "../.."))
    ))
  }

  # Method 2: already inside r/app
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(cwd) == "app" && basename(dirname(cwd)) == "r") {
    return(list(
      app_dir = cwd,
      r_code_dir = normalizePath(file.path(cwd, "..", "R")),
      project_root = normalizePath(file.path(cwd, "../.."))
    ))
  }

  # Method 3: project root (contains r/app/)
  if (dir.exists(file.path(cwd, "r", "app"))) {
    app_dir <- normalizePath(file.path(cwd, "r", "app"))
    return(list(
      app_dir = app_dir,
      r_code_dir = normalizePath(file.path(cwd, "r", "R")),
      project_root = cwd
    ))
  }

  # Method 4: inside r/ folder
  if (basename(cwd) == "r" && dir.exists(file.path(cwd, "app"))) {
    return(list(
      app_dir = normalizePath(file.path(cwd, "app")),
      r_code_dir = normalizePath(file.path(cwd, "R")),
      project_root = normalizePath(file.path(cwd, ".."))
    ))
  }

  # Method 5: walk up the directory tree
  path <- cwd
  for (i in 0:6) {
    cand <- if (i == 0) path else normalizePath(file.path(path, paste(rep("..", i), collapse = "/")), mustWork = FALSE)
    if (dir.exists(file.path(cand, "r", "app")) && dir.exists(file.path(cand, "r", "R"))) {
      return(list(
        app_dir = normalizePath(file.path(cand, "r", "app")),
        r_code_dir = normalizePath(file.path(cand, "r", "R")),
        project_root = cand
      ))
    }
  }

  stop(
    "Could not locate the project folders (r/app and r/R).\n",
    "Current working directory: ", cwd, "\n",
    "Try: setwd('/Users/isaac/Desktop/art_institute_code')"
  )
}
