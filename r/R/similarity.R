#' Binary similarity measures from simbin.m (Jeffrey M. Spielberg)
#'
#' Only the measure used in the publication (Pei1) is required for the
#' default pipeline, but several common measures are included.

simbin <- function(mat1, mat2, type = "Pei1", mask = NULL) {
  mat1 <- as.numeric(mat1 > 0)
  mat2 <- as.numeric(mat2 > 0)

  if (is.null(mask)) {
    mask <- rep(1, length(mat1))
  } else {
    mask <- as.numeric(mask > 0)
  }

  mat1 <- mat1 * mask
  mat2 <- mat2 * mask
  nmat1 <- (1 - mat1) * mask
  nmat2 <- (1 - mat2) * mask

  a <- sum(mat1 * mat2)
  b <- sum(nmat1 * mat2)
  c <- sum(mat1 * nmat2)
  d <- sum(nmat1 * nmat2)

  switch(
    type,
    Pei1 = ((a * b) - (b * c)) / ((a + c) * (b + d)),
    dice = (2 * a) / ((2 * a) + b + c),
    jaccard = a / (a + b + c),
    ochiai = a / sqrt((a + b) * (a + c)),
    stop("Unsupported simbin type: ", type)
  )
}

similarity_matrix <- function(data, type = "Pei1") {
  n <- nrow(data)
  d <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      d[i, j] <- simbin(data[i, ], data[j, ], type = type)
    }
  }
  d
}

cross_similarity_matrix <- function(query, reference, type = "Pei1") {
  n_q <- nrow(query)
  n_r <- nrow(reference)
  d <- matrix(0, nrow = n_q, ncol = n_r)
  for (i in seq_len(n_q)) {
    for (j in seq_len(n_r)) {
      d[i, j] <- simbin(query[i, ], reference[j, ], type = type)
    }
  }
  d
}
