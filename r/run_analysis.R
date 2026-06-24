#!/usr/bin/env Rscript
# One-click analysis runner for researchers
#
# Usage (from project root in Terminal or RStudio):
#   Rscript r/run_analysis.R
#
# Or in R:
#   source("r/run_analysis.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path("r", "config", "config.yaml")

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}

source(file.path("r", "R", "run_pipeline.R"), local = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  source(file.path("r", "install_packages.R"), local = TRUE)
}

results <- run_maldi_pipeline(config = config_path, verbose = TRUE)
cat("\nAnalysis complete.\n")
cat("Output folder:", normalizePath(results$output_dir), "\n")
