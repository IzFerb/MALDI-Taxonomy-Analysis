# Extract sample data and print folder-to-group mapping
#
# Usage:
#   source("r/setup_sample_data.R")

project_root <- getwd()
zip_path <- file.path(project_root, "Gums files.zip")
dest <- file.path(project_root, "sample_data")

if (!file.exists(zip_path)) {
  stop("Could not find 'Gums files.zip' in the project root.")
}

if (!dir.exists(dest)) {
  dir.create(dest, recursive = TRUE)
  message("Extracting ", zip_path, " ...")
  utils::unzip(zip_path, exdir = dest)
} else {
  message("Using existing folder: ", dest)
}

source(file.path(project_root, "r", "R", "utils.R"))
source(file.path(project_root, "r", "R", "load_spectra.R"))

meta <- load_spectra(dest, exclude_dirs = c("validation"))

cat("\nSample data layout (folder order = group ID):\n")
for (i in seq_along(meta$group_names)) {
  n <- sum(meta$group_ids == i)
  cat(sprintf("  Group %d: %s (%d files)\n", i, meta$group_names[i], n))
}
cat("\nTotal samples:", length(meta$filenames), "\n")
cat("Run validation with:\n")
cat('  Rscript r/validate_pipeline.R\n')
