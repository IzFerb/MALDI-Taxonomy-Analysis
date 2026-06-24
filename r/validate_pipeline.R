#!/usr/bin/env Rscript
# Run the R pipeline on sample data and compare against MATLAB (if available)
#
# Usage:
#   Rscript r/validate_pipeline.R              # run R pipeline + compare
#   Rscript r/validate_pipeline.R --compare-only  # compare existing exports

args <- commandArgs(trailingOnly = TRUE)
compare_only <- "--compare-only" %in% args

project_root <- getwd()
config_path <- file.path(project_root, "r", "config", "config_sample_data.yaml")

source(file.path(project_root, "r", "R", "run_pipeline.R"), local = FALSE)

if (!compare_only) {
  if (!dir.exists(file.path(project_root, "sample_data"))) {
    source(file.path(project_root, "r", "setup_sample_data.R"))
  }
  cat("Running R pipeline on sample_data/ (this may take several minutes)...\n\n")
  results <- run_maldi_pipeline(config = config_path, verbose = TRUE)
  cat("\n")
}

source(file.path(project_root, "r", "R", "validate.R"), local = FALSE)

matlab_dir <- file.path(project_root, "validation", "matlab")
if (file.exists(file.path(matlab_dir, "similarity_matrix.csv"))) {
  compare_validation(
    r_dir = file.path(project_root, "validation", "r"),
    matlab_dir = matlab_dir
  )
} else {
  cat("R validation exports are ready in validation/r/\n\n")
  cat("Next step — run MATLAB to generate reference outputs:\n")
  cat("  1. Open MATLAB with the icoshift toolbox on your path\n")
  cat("  2. cd to this project folder\n")
  cat("  3. Run: Analysis_export\n")
  cat("  4. Re-run: Rscript r/validate_pipeline.R --compare-only\n")
}
