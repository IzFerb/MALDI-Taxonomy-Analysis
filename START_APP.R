# One-line launcher — paste into the RStudio Console:
#
#   source("START_APP.R")
#
# The app opens in RStudio's Viewer pane (right side).
# To open in an external browser instead, run:
#   options(shiny.launch.browser = TRUE); source("START_APP.R")

here <- "/Users/isaac/Desktop/art_institute_code"
if (dir.exists(here)) {
  setwd(here)
} else if (file.exists("r/launch_app.R")) {
  # already in project root
} else {
  stop("Could not find project. Edit the 'here' path in START_APP.R.")
}

source("r/launch_app.R")
