# Auto-load sibling R files when run_pipeline.R is sourced
.local_r_dir <- local({
  if (exists(".r_pipeline_dir", envir = .GlobalEnv)) {
    return(get(".r_pipeline_dir", envir = .GlobalEnv))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  NA_character_
})
if (!is.na(.local_r_dir) && dir.exists(.local_r_dir)) {
  .r_siblings <- list.files(.local_r_dir, pattern = "\\.R$", full.names = TRUE)
  .r_siblings <- .r_siblings[!grepl("run_pipeline\\.R$", basename(.r_siblings))]
  invisible(lapply(.r_siblings, source, local = FALSE))
  assign(".r_pipeline_dir", .local_r_dir, envir = .GlobalEnv)
}

source_r_dir <- function() {
  r_dir <- if (exists(".r_pipeline_dir", envir = .GlobalEnv)) {
    get(".r_pipeline_dir", envir = .GlobalEnv)
  } else {
    path <- normalizePath(getwd(), mustWork = FALSE)
    found <- NA_character_
    for (i in 0:6) {
      cand <- if (i == 0) path else normalizePath(file.path(path, paste(rep("..", i), collapse = "/")), mustWork = FALSE)
      if (dir.exists(file.path(cand, "r", "R"))) {
        found <- file.path(cand, "r", "R")
        break
      }
    }
    found
  }
  if (is.na(r_dir) || !dir.exists(r_dir)) {
    stop("Could not locate r/R source directory.")
  }
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  r_files <- r_files[!grepl("run_pipeline\\.R$", basename(r_files))]
  invisible(lapply(r_files, source, local = FALSE))
}

#' Run the full MALDI analysis pipeline
#'
#' @param config Path to YAML config file or a named list of settings.
#' @param data_dir Override data directory from config.
#' @param verbose Print progress messages.
run_maldi_pipeline <- function(
    config = file.path("r", "config", "config.yaml"),
    data_dir = NULL,
    verbose = TRUE) {

  if (is.character(config)) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Install the yaml package: install.packages('yaml')")
    }
    cfg <- yaml::read_yaml(config)
  } else {
    cfg <- config
  }

  source_r_dir()

  data_dir <- data_dir %||% cfg$data_dir %||% "."
  output_dir <- cfg$output_dir %||% file.path("r", "output")
  ensure_dir(output_dir)

  step <- if (verbose) new_step_counter(8) else NULL

  if (verbose) step <- message_step(step, "Loading spectra from folders...")
  meta <- load_spectra(
    data_dir = data_dir,
    exclude_dirs = cfg$exclude_dirs %||% c("resources", "icoshift3p0", ".git", "r", "app", "output"),
    blank_pattern = cfg$blank_pattern %||% "^Blank"
  )

  if (length(meta$filenames) == 0L) {
    stop(
      "No sample .txt files found in subfolders of: ", data_dir, "\n",
      "Create one folder per sample group (e.g. Acacia/, Historical/) ",
      "and place .txt spectra inside each folder."
    )
  }

  drop_idx <- cfg$drop_sample_index %||% 5L
  if (drop_idx <= length(meta$filenames)) {
    keep <- seq_along(meta$filenames) != drop_idx
    meta$filenames <- meta$filenames[keep]
    meta$group_ids <- meta$group_ids[keep]
  }

  ref_groups <- cfg$reference_groups %||% c(3L, 4L, 6L)
  reference_idx <- which(meta$group_ids %in% ref_groups)

  mz_grid <- seq(
    cfg$mz_min %||% 800,
    cfg$mz_max %||% 1300,
    by = cfg$mz_step %||% 0.01
  )

  if (verbose) step <- message_step(step, "Interpolating spectra...")
  interp <- interpolate_spectra(data_dir, meta$filenames, mz_grid)

  if (verbose) step <- message_step(step, "Aligning spectra (icoshift)...")
  aligned <- icoshift(
    target = cfg$alignment_target %||% "median",
    spectra = interp$matrix,
    mode = cfg$alignment_mode %||% "whole",
    max_shift = cfg$max_shift %||% 34L
  )

  if (isTRUE(cfg$save_alignment_plot)) {
    png(file.path(output_dir, "01_aligned_spectra.png"), width = 1200, height = 700)
    matplot(mz_grid, t(aligned$spectra), type = "l", lty = 1,
            xlab = "m/z", ylab = "Intensity", main = "Aligned spectra")
    dev.off()
  }

  if (verbose) step <- message_step(step, "Standardizing spectra (row-wise SNV)...")
  x_norm <- snv_rows(aligned$spectra)

  if (isTRUE(cfg$save_snv_plot)) {
    png(file.path(output_dir, "02_snv_spectra.png"), width = 1200, height = 700)
    matplot(mz_grid, t(x_norm), type = "l", lty = 1,
            xlab = "m/z", ylab = "SNV intensity", main = "Standardized spectra")
    dev.off()
  }

  if (verbose) step <- message_step(step, "Peak picking and alignment...")
  peaks <- find_and_align_peaks(
    x = x_norm,
    mz = mz_grid,
    reference_idx = reference_idx,
    prominence_factor = cfg$peak_prominence_factor %||% 10
  )

  n_blanks <- cfg$n_blank_samples %||% 4L
  blanked <- remove_blank_peaks(peaks$matrix, n_blanks = n_blanks)
  z <- blanked$matrix

  exclude_group <- cfg$exclude_group_id %||% 7L
  keep_samples <- apply(z, 1, function(row) any(row > 0)) & meta$group_ids != exclude_group
  z2 <- z[keep_samples, , drop = FALSE]
  filenames_kept <- meta$filenames[keep_samples]
  group_kept <- meta$group_ids[keep_samples]

  acacia_group <- cfg$acacia_group_id %||% 3L
  acacia_idx <- which(group_kept == acacia_group)

  if (length(acacia_idx) == 0L) {
    warning("No Acacia samples found (group id ", acacia_group, "). Skipping taxonomy PCA.")
    return(list(
      config = cfg,
      meta = meta,
      z = z,
      z2 = z2,
      output_dir = output_dir
    ))
  }

  if (verbose) step <- message_step(step, "Building Acacia similarity matrix...")
  z_acacia <- z2[acacia_idx, , drop = FALSE]
  fn_acacia <- filenames_kept[acacia_idx]

  class_labels <- vapply(fn_acacia, parse_class_from_filename, character(1), anchor = "Acacia")
  manual <- cfg$manual_class_overrides %||% list()
  for (nm in names(manual)) {
    idx <- as.integer(nm)
    if (idx <= length(class_labels)) class_labels[idx] <- manual[[nm]]
  }

  drop_acacia <- cfg$drop_acacia_index %||% 20L
  if (drop_acacia <= length(class_labels)) {
    keep_a <- seq_along(class_labels) != drop_acacia
    z_acacia <- z_acacia[keep_a, , drop = FALSE]
    class_labels <- class_labels[keep_a]
    fn_acacia <- fn_acacia[keep_a]
  }

  d <- similarity_matrix(z_acacia, type = cfg$similarity_type %||% "Pei1")
  class_labels[is.na(class_labels)] <- "unknown"
  colors <- distinguishable_colors(nlevels(factor(class_labels)))

  if (verbose) step <- message_step(step, "Generating PCA plots...")
  pc_pairs <- cfg$pca_pairs %||% list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4))
  n_pc_max <- min(nrow(d), ncol(d))
  plots <- list()
  for (pair in pc_pairs) {
    if (any(pair > n_pc_max)) next
    res <- plot_pca_metabolomics(
      x = d,
      groups = class_labels,
      colors = colors,
      text_labels = seq_along(class_labels),
      pc_x = pair[[1]],
      pc_y = pair[[2]],
      title = sprintf("Acacia PCA (PC %d vs PC %d)", pair[[1]], pair[[2]])
    )
    plots[[length(plots) + 1L]] <- res$plot
  }
  save_pca_plots(plots, output_dir, prefix = "acacia_pca")

  results <- list(
    config = cfg,
    meta = meta,
    mz = mz_grid,
    aligned = aligned$spectra,
    normalized = x_norm,
    feature_matrix = z,
    feature_matrix_filtered = z2,
    acacia_features = z_acacia,
    acacia_classes = class_labels,
    similarity_matrix = d,
    output_dir = output_dir,
    plots = plots
  )

  historical_group <- cfg$historical_group_id %||% 7L
  hist_idx <- which(meta$group_ids == historical_group)
  if (length(hist_idx) > 0L) {
    z_hist <- z[hist_idx, , drop = FALSE]
    drop_hist <- cfg$drop_historical_index %||% c(6L, 9L, 10L)
    keep_h <- seq_len(nrow(z_hist)) %notin% drop_hist
    z_hist <- z_hist[keep_h, , drop = FALSE]

    d_hist <- cross_similarity_matrix(z_hist, z_acacia, type = cfg$similarity_type %||% "Pei1")
    class2 <- c(class_labels, rep("historical", nrow(d_hist)))
    combined <- rbind(d, d_hist)
    trainset <- c(rep(TRUE, nrow(d)), rep(FALSE, nrow(d_hist)))

    hist_plot <- plot_pca_metabolomics(
      x = combined,
      groups = class2,
      colors = distinguishable_colors(nlevels(factor(class2))),
      pc_x = 1L,
      pc_y = 2L,
      trainset = trainset,
      title = "Acacia + historical samples PCA"
    )
    ggplot2::ggsave(
      file.path(output_dir, "historical_pca.png"),
      plot = hist_plot$plot,
      width = 8,
      height = 7,
      dpi = 150
    )
    results$historical_similarity <- d_hist
    results$historical_plot <- hist_plot$plot
  }

  pigmented_group <- cfg$pigmented_group_id %||% 5L
  pig_idx <- which(meta$group_ids == pigmented_group)
  if (length(pig_idx) > 0L) {
    z_pig <- z[pig_idx, , drop = FALSE]
    d_pig <- cross_similarity_matrix(z_pig, z_acacia, type = cfg$similarity_type %||% "Pei1")
    class3 <- c(class_labels, rep("pigmented", nrow(d_pig)))
    combined <- rbind(d, d_pig)
    trainset <- c(rep(TRUE, nrow(d)), rep(FALSE, nrow(d_pig)))

    pig_plot <- plot_pca_metabolomics(
      x = combined,
      groups = class3,
      colors = distinguishable_colors(nlevels(factor(class3))),
      pc_x = 1L,
      pc_y = 2L,
      trainset = trainset,
      title = "Acacia + pigmented references PCA"
    )
    ggplot2::ggsave(
      file.path(output_dir, "pigmented_pca.png"),
      plot = pig_plot$plot,
      width = 8,
      height = 7,
      dpi = 150
    )
    results$pigmented_similarity <- d_pig
    results$pigmented_plot <- pig_plot$plot
  }

  utils::write.csv(d, file.path(output_dir, "acacia_similarity_matrix.csv"))
  utils::write.csv(
    data.frame(sample = fn_acacia, class = class_labels),
    file.path(output_dir, "acacia_sample_classes.csv"),
    row.names = FALSE
  )

  if (isTRUE(cfg$export_validation %||% FALSE)) {
    val_dir <- cfg$validation_dir %||% "validation/r"
    export_validation_artifacts(results, val_dir, meta = meta)
    if (verbose) message("Validation exports saved to: ", normalizePath(val_dir))
  }

  if (verbose) message("Done. Results saved to: ", normalizePath(output_dir))
  invisible(results)
}

`%notin%` <- function(x, y) !(x %in% y)
