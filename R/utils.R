.logit_inv <- function(x) {
  stats::plogis(x)
}

.is_positive_definite <- function(x, tol = 1e-8) {
  values <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  all(is.finite(values)) && all(values > tol)
}

.make_positive_definite <- function(x, eps = 1e-6) {
  x <- 0.5 * (x + t(x))
  eig <- eigen(x, symmetric = TRUE)
  eig$values <- pmax(eig$values, eps)
  eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*%
    t(eig$vectors)
}

.safe_solve <- function(x) {
  tryCatch(solve(x), error = function(e) MASS::ginv(x))
}

.matrix_power <- function(x, power = 0.5, tol = 1e-10) {
  x <- 0.5 * (x + t(x))
  eig <- eigen(x, symmetric = TRUE)
  if (any(eig$values < -tol)) {
    stop("Matrix has negative eigenvalues.", call. = FALSE)
  }
  values <- pmax(eig$values, tol)^power
  eig$vectors %*% diag(values, nrow = length(values)) %*%
    t(eig$vectors)
}

.hz_beta <- function(k, d = 2L) {
  (((2 * d + 1) * k / 4)^(1 / (d + 4))) / sqrt(2)
}

.apply_continuity_correction <- function(data, cc = 0.5) {
  corrected <- data
  cells <- c("TP", "FN", "FP", "TN")
  affected <- apply(corrected[, cells, drop = FALSE] == 0, 1L, any)
  corrected[affected, cells] <- corrected[affected, cells, drop = FALSE] + cc
  corrected
}

.validate_dta_table <- function(data) {
  required <- c("study_id", "TP", "FN", "FP", "TN")
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  data <- as.data.frame(data[, required], stringsAsFactors = FALSE)
  if (nrow(data) < 2L) {
    stop("At least two studies are required.", call. = FALSE)
  }
  if (anyDuplicated(data$study_id)) {
    stop("study_id must uniquely identify each study.", call. = FALSE)
  }
  counts <- as.matrix(data[, c("TP", "FN", "FP", "TN")])
  storage.mode(counts) <- "double"
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Cell counts must be finite and nonnegative.", call. = FALSE)
  }
  if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    stop("Cell counts must be integers.", call. = FALSE)
  }
  if (any(rowSums(counts) == 0)) {
    stop("Every study must contain at least one participant.", call. = FALSE)
  }
  data[, c("TP", "FN", "FP", "TN")] <- counts
  data
}

.align_table_to_effects <- function(data, effects) {
  effect_ids <- rownames(effects)
  if (is.null(effect_ids)) {
    if (nrow(data) != nrow(effects)) {
      stop("The study table and random effects cannot be aligned.", call. = FALSE)
    }
    return(data)
  }
  index <- match(effect_ids, as.character(data$study_id))
  if (anyNA(index)) {
    stop("The study table and random effects cannot be aligned by study_id.",
         call. = FALSE)
  }
  data[index, , drop = FALSE]
}

