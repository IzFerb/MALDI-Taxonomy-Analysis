required_packages <- c(
  "ggplot2",
  "yaml",
  "shiny",
  "DT",
  "colorspace"
)

missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0L) {
  message("Installing required packages: ", paste(missing, collapse = ", "))
  utils::install.packages(missing, repos = "https://cloud.r-project.org")
}

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Failed to install package: ", pkg)
  }
}

message("All dependencies are ready.")
