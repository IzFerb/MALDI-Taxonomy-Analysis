#' PCA utilities matching plot_pca_metabolomics.m

pca_function <- function(x) {
  x <- as.matrix(x)
  svd_result <- svd(x, nu = min(nrow(x), ncol(x)), nv = ncol(x))
  scores <- svd_result$u %*% diag(svd_result$d)
  loadings <- svd_result$v
  variance <- round(100 * (svd_result$d^2) / sum(svd_result$d^2), 1)
  list(scores = scores, loadings = loadings, variance = variance)
}

prepare_pca_input <- function(x, trainset = NULL) {
  if (is.null(trainset)) trainset <- rep(TRUE, nrow(x))
  x_mean <- colMeans(x[trainset, , drop = FALSE], na.rm = TRUE)
  x_asc <- sweep(x, 2, x_mean, "-")
  x_asc <- x_asc[, colSums(is.nan(x_asc)) == 0, drop = FALSE]
  x_asc
}

plot_pca_metabolomics <- function(
    x,
    groups,
    colors = NULL,
    legend_labels = NULL,
    text_labels = NULL,
    pc_x = 1L,
    pc_y = 2L,
    trainset = NULL,
    var_select = NULL,
    title = NULL) {

  x_asc <- prepare_pca_input(x, trainset)
  if (!is.null(var_select)) x_asc <- x_asc[, var_select, drop = FALSE]

  pca <- pca_function(x_asc[trainset %||% rep(TRUE, nrow(x_asc)), , drop = FALSE])
  scores <- x_asc %*% pca$loadings

  groups <- factor(groups)
  n_groups <- nlevels(groups)
  if (is.null(colors)) {
    colors <- distinguishable_colors(n_groups)
  } else {
    colors <- rep(colors, length.out = n_groups)
  }
  color_map <- stats::setNames(colors[seq_len(n_groups)], levels(groups))

  df <- data.frame(
    PCx = scores[, pc_x],
    PCy = scores[, pc_y],
    Group = groups
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PCx, y = PCy, color = Group)) +
    ggplot2::geom_point(size = 3, stroke = 1) +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::labs(
      title = title,
      x = sprintf("PC %d (%s%% variance explained)", pc_x, pca$variance[pc_x]),
      y = sprintf("PC %d (%s%% variance explained)", pc_y, pca$variance[pc_y])
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(aspect.ratio = 1, legend.title = ggplot2::element_blank())

  if (!is.null(text_labels) && length(text_labels) > 0L) {
    df$Label <- as.character(text_labels)
    p <- p + ggplot2::geom_text(
      data = df,
      ggplot2::aes(x = .data$PCx, y = .data$PCy, label = .data$Label),
      hjust = -0.1,
      size = 3,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  list(plot = p, scores = scores, loadings = pca$loadings, variance = pca$variance)
}

distinguishable_colors <- function(n, bg = matrix(c(0, 0, 0, 1, 1, 1), nrow = 2, byrow = TRUE)) {
  grDevices::hcl(h = seq(15, 375, length.out = n + 1)[seq_len(n)], c = 100, l = 65)
}

save_pca_plots <- function(plot_objects, output_dir, prefix = "pca") {
  ensure_dir(output_dir)
  paths <- character()
  for (i in seq_along(plot_objects)) {
    path <- file.path(output_dir, sprintf("%s_%02d.png", prefix, i))
    ggplot2::ggsave(path, plot = plot_objects[[i]], width = 8, height = 7, dpi = 150)
    paths <- c(paths, path)
  }
  paths
}
