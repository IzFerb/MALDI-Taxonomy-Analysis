#' Export matrices for cross-validation against MATLAB
export_validation_artifacts <- function(results, validation_dir, meta = NULL) {
  ensure_dir(validation_dir)

  if (!is.null(results$similarity_matrix)) {
    utils::write.csv(
      results$similarity_matrix,
      file.path(validation_dir, "similarity_matrix.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(results$acacia_features)) {
    utils::write.csv(
      results$acacia_features,
      file.path(validation_dir, "acacia_feature_matrix.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(results$feature_matrix)) {
    utils::write.csv(
      results$feature_matrix,
      file.path(validation_dir, "feature_matrix.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(results$acacia_classes)) {
    utils::write.csv(
      data.frame(class = results$acacia_classes),
      file.path(validation_dir, "acacia_classes.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(meta)) {
    utils::write.csv(
      data.frame(
        index = seq_along(meta$filenames),
        filename = meta$filenames,
        group_id = meta$group_ids,
        group_name = meta$group_names[meta$group_ids]
      ),
      file.path(validation_dir, "sample_manifest.csv"),
      row.names = FALSE
    )
  }

  summary_lines <- c(
    paste("n_samples_total", length(meta$filenames %||% NA), sep = "="),
    paste("n_features", ncol(results$feature_matrix %||% matrix(0, 0, 0)), sep = "="),
    paste("n_acacia_samples", length(results$acacia_classes %||% character()), sep = "="),
    paste("n_acacia_features", ncol(results$acacia_features %||% matrix(0, 0, 0)), sep = "=")
  )
  if (!is.null(results$similarity_matrix)) {
    d <- results$similarity_matrix
    summary_lines <- c(
      summary_lines,
      paste("similarity_dim", paste(dim(d), collapse = "x"), sep = "="),
      paste("similarity_mean", mean(d, na.rm = TRUE), sep = "="),
      paste("similarity_sd", stats::sd(d, na.rm = TRUE), sep = "=")
    )
  }
  writeLines(summary_lines, file.path(validation_dir, "summary.txt"))
  invisible(validation_dir)
}

#' Compare R and MATLAB validation exports
compare_validation <- function(
    r_dir = "validation/r",
    matlab_dir = "validation/matlab",
    tolerance = 1e-4) {

  r_sim <- file.path(r_dir, "similarity_matrix.csv")
  m_sim <- file.path(matlab_dir, "similarity_matrix.csv")

  if (!file.exists(r_sim)) {
    stop("R similarity matrix not found. Run validate_pipeline.R first.")
  }
  if (!file.exists(m_sim)) {
    stop(
      "MATLAB reference not found at ", m_sim, ".\n",
      "Run Analysis_export.m in MATLAB first (see validation/README.md)."
    )
  }

  r_mat <- as.matrix(utils::read.csv(r_sim, check.names = FALSE))
  m_mat <- as.matrix(utils::read.csv(m_sim, check.names = FALSE))

  if (!all(dim(r_mat) == dim(m_mat))) {
    stop(
      "Dimension mismatch: R is ", paste(dim(r_mat), collapse = "x"),
      " but MATLAB is ", paste(dim(m_mat), collapse = "x"),
      ". Check sample order in sample_manifest.csv."
    )
  }

  diff_mat <- r_mat - m_mat
  max_abs <- max(abs(diff_mat), na.rm = TRUE)
  rmse <- sqrt(mean(diff_mat^2, na.rm = TRUE))
  corr <- stats::cor(c(r_mat), c(m_mat), use = "complete.obs")

  report <- list(
    dimensions = dim(r_mat),
    max_absolute_difference = max_abs,
    rmse = rmse,
    correlation = corr,
    within_tolerance = max_abs <= tolerance,
    tolerance = tolerance
  )

  report_lines <- c(
    "MALDI pipeline validation report",
    "=============================",
    paste("Matrix dimensions:", paste(report$dimensions, collapse = " x ")),
    paste("Max absolute difference:", signif(report$max_absolute_difference, 6)),
    paste("RMSE:", signif(report$rmse, 6)),
    paste("Pearson correlation:", signif(report$correlation, 6)),
    paste("Tolerance:", report$tolerance),
    if (report$within_tolerance) {
      "PASS: R and MATLAB similarity matrices match within tolerance."
    } else {
      "REVIEW: Differences exceed tolerance. See validation/r/difference_matrix.csv"
    }
  )

  ensure_dir(r_dir)
  utils::write.csv(diff_mat, file.path(r_dir, "difference_matrix.csv"), row.names = FALSE)
  writeLines(report_lines, file.path(r_dir, "comparison_report.txt"))

  cat(paste(report_lines, collapse = "\n"), "\n")
  invisible(report)
}
