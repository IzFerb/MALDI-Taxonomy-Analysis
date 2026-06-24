#!/usr/bin/env Rscript
# Launch the Shiny analysis app
#
# In RStudio (from project root):
#   source("r/launch_app.R")
#
# Or use the simpler wrapper:
#   source("START_APP.R")

if (!requireNamespace("shiny", quietly = TRUE)) {
  source(file.path(getwd(), "r", "install_packages.R"))
}

# --- Find project root ---
start_dir <- normalizePath(getwd(), mustWork = FALSE)

if (basename(start_dir) == "app" && basename(dirname(start_dir)) == "r") {
  project_root <- normalizePath(file.path(start_dir, "../.."))
} else if (dir.exists(file.path(start_dir, "r", "R"))) {
  project_root <- start_dir
} else if (basename(start_dir) == "r" && dir.exists(file.path(start_dir, "R"))) {
  project_root <- normalizePath(file.path(start_dir, ".."))
} else {
  stop(
    "Cannot find project root.\n",
    "Current directory: ", start_dir, "\n",
    "Run: setwd('/Users/isaac/Desktop/art_institute_code')"
  )
}

app_dir <- file.path(project_root, "r", "app")
if (!file.exists(file.path(app_dir, "app.R"))) {
  stop("App not found: ", file.path(app_dir, "app.R"))
}

options(
  gums.project_root = project_root,
  gums.r_code_dir   = file.path(project_root, "r", "R")
)

message("Project root: ", project_root)
message("Starting Shiny app...")

# Open app in RStudio Viewer, system browser, or print URL as fallback
launch_browser <- function(url) {
  message("\n========================================")
  message("App running at: ", url)
  message("========================================\n")

  opened <- FALSE

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    tryCatch({
      rstudioapi::viewer(url)
      message("Opened in RStudio Viewer pane (right side of IDE).")
      message("If you don't see it: click the 'Show in new window' icon in the Viewer.\n")
      opened <- TRUE
    }, error = function(e) NULL)
  }

  if (!opened) {
    tryCatch({
      utils::browseURL(url)
      message("Opened in your default web browser.\n")
      opened <- TRUE
    }, error = function(e) NULL)
  }

  if (!opened) {
    message("Could not open automatically.")
    message("Copy this URL into Chrome or Safari: ", url, "\n")
  }
}

# In RStudio console, launch.browser must be TRUE or a function; FALSE opens nothing
shiny::runApp(
  appDir = app_dir,
  launch.browser = launch_browser,
  host = "127.0.0.1"
)
