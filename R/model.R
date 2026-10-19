.fit_dta_glmm <- function(ipd, nAGQ = 7L,
                          fixed_formula = cbind(correct_classification, misclassification) ~ group - 1) {
  GLMMadaptive::mixed_model(
    fixed = fixed_formula,
    random = ~ group - 1 | study_id,
    data = ipd,
    family = stats::binomial(),
    nAGQ = nAGQ,
    control = list(iter_EM = 0L, verbose = FALSE)
  )
}

.extract_random_effects <- function(model) {
  effects <- as.matrix(nlme::ranef(model))
  if (ncol(effects) != 2L || any(!is.finite(effects))) {
    stop("The fitted model did not produce two finite random effects per study.",
         call. = FALSE)
  }
  colnames(effects) <- c("Sensitivity", "Specificity")
  effects
}

.weighted_effects <- function(model, corrected_data) {
  effects <- .extract_random_effects(model)
  corrected_data <- .align_table_to_effects(corrected_data, effects)
  k <- nrow(effects)

  Sigma <- 0.5 * (model$D + t(model$D))
  if (!.is_positive_definite(Sigma)) {
    Sigma <- .make_positive_definite(Sigma)
  }

  var_se <- 1 / corrected_data$TP + 1 / corrected_data$FN
  var_sp <- 1 / corrected_data$FP + 1 / corrected_data$TN
  variance <- lapply(seq_len(k), function(i) {
    value <- Sigma + diag(c(var_se[i], var_sp[i]), nrow = 2L)
    value <- 0.5 * (value + t(value))
    if (!.is_positive_definite(value)) value <- .make_positive_definite(value)
    value
  })
  precision <- lapply(variance, .safe_solve)

  total_precision <- Reduce(`+`, precision)
  total_precision <- 0.5 * (total_precision + t(total_precision))
  if (!.is_positive_definite(total_precision)) {
    total_precision <- .make_positive_definite(total_precision)
  }
  total_precision_inv <- .safe_solve(total_precision)
  numerator <- Reduce(`+`, lapply(seq_len(k), function(i) {
    precision[[i]] %*% matrix(effects[i, ], ncol = 1L)
  }))
  weighted_mean <- as.vector(total_precision_inv %*% numerator)

  middle <- matrix(0, nrow = 2L, ncol = 2L)
  for (i in seq_len(k)) {
    difference <- matrix(effects[i, ] - weighted_mean, ncol = 1L)
    precision_sqrt <- .matrix_power(precision[[i]], 0.5)
    middle <- middle + precision_sqrt %*%
      (difference %*% t(difference)) %*% precision_sqrt
  }
  total_precision_inv_sqrt <- .matrix_power(total_precision, -0.5)
  weighted_covariance <- total_precision_inv_sqrt %*% middle %*%
    total_precision_inv_sqrt
  weighted_covariance <- 0.5 * (weighted_covariance + t(weighted_covariance))
  if (!.is_positive_definite(weighted_covariance)) {
    weighted_covariance <- .make_positive_definite(weighted_covariance)
  }

  covariance_inv_sqrt <- .matrix_power(weighted_covariance, -0.5)
  standardized <- t(apply(effects, 1L, function(value) {
    as.vector(covariance_inv_sqrt %*% matrix(value - weighted_mean, ncol = 1L))
  }))
  rownames(standardized) <- rownames(effects)
  colnames(standardized) <- c("Sensitivity", "Specificity")

  list(
    effects = effects,
    standardized = standardized,
    weighted_mean = weighted_mean,
    weighted_covariance = weighted_covariance,
    Sigma = Sigma,
    variance = variance,
    precision = precision
  )
}
