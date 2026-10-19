#' Weighted Henze-Zirkler Statistic
#'
#' Computes the proposed precision-weighted HZ statistic from a fitted
#' bivariate diagnostic-accuracy GLMM and the corresponding study tables.
#'
#' @param model A model returned by `GLMMadaptive::mixed_model()` using the
#'   package's bivariate diagnostic-accuracy specification.
#' @param data A data frame with raw `TP`, `FN`, `FP`, and `TN` counts.
#' @param cc Continuity correction added to all four cells of a study when at
#'   least one cell is zero.
#' @param beta Optional HZ smoothing parameter. The standard choice is used
#'   when `NULL`.
#' @return A list containing the statistic, beta, standardized effects, and
#'   weighting components.
#' @export
weighted_hz_statistic <- function(model, data, cc = 0.5, beta = NULL) {
  data <- .validate_dta_table(data)
  if (!is.numeric(cc) || length(cc) != 1L || !is.finite(cc) || cc <= 0) {
    stop("cc must be one positive finite number.", call. = FALSE)
  }
  corrected <- .apply_continuity_correction(data, cc)
  weighted <- .weighted_effects(model, corrected)
  z <- weighted$standardized
  k <- nrow(z)
  d <- ncol(z)
  if (is.null(beta)) beta <- .hz_beta(k, d)
  if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta) || beta <= 0) {
    stop("beta must be one positive finite number.", call. = FALSE)
  }

  squared_distance <- as.matrix(stats::dist(z))^2
  term1 <- sum(exp(-0.5 * beta^2 * squared_distance)) / k^2
  term2 <- 2 * (1 + beta^2)^(-d / 2) *
    mean(exp(-(beta^2 / (2 * (1 + beta^2))) * rowSums(z^2)))
  term3 <- (1 + 2 * beta^2)^(-d / 2)
  statistic <- max(as.numeric(k * (term1 - term2 + term3)), 0)

  c(list(statistic = statistic, beta = beta), weighted)
}

#' Classical Henze-Zirkler Test for Empirical Bayes Effects
#'
#' @param model A fitted bivariate diagnostic-accuracy GLMM.
#' @return A list containing the HZ statistic and p-value.
#' @export
classical_hz_test <- function(model) {
  effects <- .extract_random_effects(model)
  result <- mvnTest::HZ.test(effects, qqplot = FALSE)
  list(
    statistic = as.numeric(result@HZ),
    p_value = as.numeric(result@p.value),
    effects = effects
  )
}

#' Weighted HZ Test with Parametric Bootstrap Calibration
#'
#' @param data A data frame with one row per study and columns `study_id`,
#'   `TP`, `FN`, `FP`, and `TN`.
#' @param B Number of parametric bootstrap samples.
#' @param cc Continuity correction for within-study variance calculations.
#' @param nAGQ Number of adaptive Gaussian quadrature points used by the GLMM.
#' @param seed Optional random-number seed.
#' @param keep_boot Logical; retain bootstrap statistics in the result?
#' @param parallel Run bootstrap replications in parallel? This is supported on
#'   Windows, macOS, and Linux.
#' @param n_cores Number of worker processes when `parallel = TRUE`.
#' @return An object of class `dtaIPDMA_test`.
#' @export
weighted_hz_test <- function(data, B = 499L, cc = 0.5, nAGQ = 7L,
                             seed = NULL, keep_boot = FALSE,
                             parallel = FALSE, n_cores = 1L) {
  data <- .validate_dta_table(data)
  B <- as.integer(B)
  nAGQ <- as.integer(nAGQ)
  if (length(B) != 1L || is.na(B) || B < 1L) {
    stop("B must be a positive integer.", call. = FALSE)
  }
  if (length(nAGQ) != 1L || is.na(nAGQ) || nAGQ < 1L) {
    stop("nAGQ must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be either TRUE or FALSE.", call. = FALSE)
  }
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    stop("n_cores must be a positive integer.", call. = FALSE)
  }
  if (isTRUE(parallel)) {
    available_cores <- parallel::detectCores(logical = TRUE)
    if (is.na(available_cores) || available_cores < 2L) {
      stop("Parallel processing requires at least two available cores.", call. = FALSE)
    }
    n_cores <- min(n_cores, available_cores, B)
  } else {
    n_cores <- 1L
  }
  if (!is.null(seed)) set.seed(seed)

  ipd <- table2x2_to_ipd(data)
  fixed_formula <- cbind(correct_classification, misclassification) ~ group - 1
  model <- .fit_dta_glmm(ipd, nAGQ, fixed_formula)
  observed_details <- weighted_hz_statistic(model, data, cc)
  observed <- observed_details$statistic
  beta_hat <- as.numeric(nlme::fixef(model))
  D_hat <- 0.5 * (model$D + t(model$D))
  if (!.is_positive_definite(D_hat)) D_hat <- .make_positive_definite(D_hat)

  studies <- unique(ipd$study_id)
  k <- length(studies)
  bootstrap_indices <- seq_len(B)
  if (isTRUE(parallel)) {
    cluster <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, library(dtaIPDMA))
    parallel::clusterExport(
      cluster,
      varlist = c("ipd", "studies", "beta_hat", "D_hat", "nAGQ", "cc", "fixed_formula"),
      envir = environment()
    )
    if (!is.null(seed)) parallel::clusterSetRNGStream(cluster, iseed = seed)
    bootstrap_statistics <- unlist(parallel::parLapply(cluster, bootstrap_indices, function(i) {
      utils::getFromNamespace(".bootstrap_hz_once", "dtaIPDMA")(
        ipd = ipd,
        studies = studies,
        beta_hat = beta_hat,
        D_hat = D_hat,
        nAGQ = nAGQ,
        cc = cc,
        fixed_formula = fixed_formula
      )
    }), use.names = FALSE)
  } else {
    bootstrap_statistics <- vapply(bootstrap_indices, function(i) {
      .bootstrap_hz_once(ipd, studies, beta_hat, D_hat, nAGQ, cc, fixed_formula)
    }, numeric(1))
  }

  valid <- is.finite(bootstrap_statistics)
  effective_B <- sum(valid)
  p_value <- if (effective_B) {
    (sum(bootstrap_statistics[valid] >= observed) + 1) / (effective_B + 1)
  } else {
    NA_real_
  }
  result <- list(
    method = "Weighted Henze-Zirkler test with parametric bootstrap",
    statistic = observed,
    p_value = p_value,
    B = B,
    B_effective = effective_B,
    failed_bootstraps = B - effective_B,
    parallel = parallel,
    n_cores = n_cores,
    beta = observed_details$beta,
    model = model,
    effects = observed_details$effects,
    standardized_effects = observed_details$standardized,
    weighted_mean = observed_details$weighted_mean,
    weighted_covariance = observed_details$weighted_covariance,
    data = data,
    call = match.call()
  )
  if (isTRUE(keep_boot)) result$bootstrap_statistics <- bootstrap_statistics
  class(result) <- "dtaIPDMA_test"
  result
}

#' Compare Classical and Weighted HZ Tests
#'
#' @inheritParams weighted_hz_test
#' @return A list containing both tests and a two-row results table.
#' @export
dta_hz_test <- function(data, B = 499L, cc = 0.5, nAGQ = 7L,
                        seed = NULL, keep_boot = FALSE,
                        parallel = FALSE, n_cores = 1L) {
  weighted <- weighted_hz_test(
    data, B, cc, nAGQ, seed, keep_boot, parallel, n_cores
  )
  classical <- classical_hz_test(weighted$model)
  summary <- data.frame(
    method = c("Classical HZ", "Weighted HZ (parametric bootstrap)"),
    statistic = c(classical$statistic, weighted$statistic),
    p_value = c(classical$p_value, weighted$p_value),
    B_effective = c(NA_integer_, weighted$B_effective)
  )
  list(classical = classical, weighted = weighted, summary = summary)
}

#' HZ Tests for Individual Participant Data Meta-Analysis
#'
#' Fits a covariate-adjusted bivariate diagnostic-accuracy GLMM directly to
#' individual participant data. Covariates are centred within study and their
#' effects are allowed to differ between diseased and non-diseased participants.
#'
#' @param data A participant-level data frame.
#' @param study_id Name of the study-identifier column.
#' @param disease Name of a numeric disease-status column coded 0 and 1.
#' @param test Name of a numeric index-test-result column coded 0 and 1.
#' @param covariates Optional character vector of numeric or binary categorical
#'   covariate column names.
#' @inheritParams weighted_hz_test
#' @return A list containing classical and weighted tests and a two-row summary.
#' @export
ipd_hz_test <- function(data, study_id, disease, test, covariates = NULL,
                        B = 499L, cc = 0.5, nAGQ = 7L, seed = NULL,
                        keep_boot = FALSE, parallel = FALSE, n_cores = 1L) {
  prepared <- .prepare_ipd_data(data, study_id, disease, test, covariates)
  ipd <- prepared$ipd
  fixed_formula <- prepared$fixed_formula
  B <- as.integer(B)
  nAGQ <- as.integer(nAGQ)
  if (length(B) != 1L || is.na(B) || B < 1L) {
    stop("B must be a positive integer.", call. = FALSE)
  }
  if (length(nAGQ) != 1L || is.na(nAGQ) || nAGQ < 1L) {
    stop("nAGQ must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be either TRUE or FALSE.", call. = FALSE)
  }
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    stop("n_cores must be a positive integer.", call. = FALSE)
  }
  if (isTRUE(parallel)) {
    available_cores <- parallel::detectCores(logical = TRUE)
    if (is.na(available_cores) || available_cores < 2L) {
      stop("Parallel processing requires at least two available cores.", call. = FALSE)
    }
    n_cores <- min(n_cores, available_cores, B)
  } else {
    n_cores <- 1L
  }
  if (!is.null(seed)) set.seed(seed)

  raw_tables <- .ipd_to_2x2(ipd)
  model <- .fit_dta_glmm(ipd, nAGQ, fixed_formula)
  observed_details <- weighted_hz_statistic(model, raw_tables, cc)
  beta_hat <- as.numeric(nlme::fixef(model))
  D_hat <- 0.5 * (model$D + t(model$D))
  if (!.is_positive_definite(D_hat)) D_hat <- .make_positive_definite(D_hat)
  studies <- unique(ipd$study_id)
  bootstrap_indices <- seq_len(B)

  if (isTRUE(parallel)) {
    cluster <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, library(dtaIPDMA))
    parallel::clusterExport(
      cluster,
      varlist = c("ipd", "studies", "beta_hat", "D_hat", "nAGQ", "cc", "fixed_formula"),
      envir = environment()
    )
    if (!is.null(seed)) parallel::clusterSetRNGStream(cluster, iseed = seed)
    bootstrap_statistics <- unlist(parallel::parLapply(cluster, bootstrap_indices, function(i) {
      utils::getFromNamespace(".bootstrap_hz_once", "dtaIPDMA")(
        ipd = ipd,
        studies = studies,
        beta_hat = beta_hat,
        D_hat = D_hat,
        nAGQ = nAGQ,
        cc = cc,
        fixed_formula = fixed_formula
      )
    }), use.names = FALSE)
  } else {
    bootstrap_statistics <- vapply(bootstrap_indices, function(i) {
      .bootstrap_hz_once(ipd, studies, beta_hat, D_hat, nAGQ, cc, fixed_formula)
    }, numeric(1))
  }

  valid <- is.finite(bootstrap_statistics)
  effective_B <- sum(valid)
  weighted <- list(
    method = "Weighted Henze-Zirkler test with parametric bootstrap (IPD)",
    statistic = observed_details$statistic,
    p_value = if (effective_B) {
      (sum(bootstrap_statistics[valid] >= observed_details$statistic) + 1) /
        (effective_B + 1)
    } else NA_real_,
    B = B,
    B_effective = effective_B,
    failed_bootstraps = B - effective_B,
    parallel = parallel,
    n_cores = n_cores,
    beta = observed_details$beta,
    model = model,
    effects = observed_details$effects,
    standardized_effects = observed_details$standardized,
    weighted_mean = observed_details$weighted_mean,
    weighted_covariance = observed_details$weighted_covariance,
    data = data,
    covariates = prepared$covariates,
    call = match.call()
  )
  if (isTRUE(keep_boot)) weighted$bootstrap_statistics <- bootstrap_statistics
  class(weighted) <- "dtaIPDMA_test"
  classical <- classical_hz_test(model)
  summary <- data.frame(
    method = c("Classical HZ", "Weighted HZ (parametric bootstrap)"),
    statistic = c(classical$statistic, weighted$statistic),
    p_value = c(classical$p_value, weighted$p_value),
    B_effective = c(NA_integer_, weighted$B_effective)
  )
  list(classical = classical, weighted = weighted, summary = summary)
}

.bootstrap_hz_once <- function(ipd, studies, beta_hat, D_hat, nAGQ, cc,
                                fixed_formula = cbind(correct_classification, misclassification) ~ group - 1) {
  k <- length(studies)
  random_effects <- tryCatch(
    MASS::mvrnorm(k, mu = c(0, 0), Sigma = D_hat),
    error = function(e) {
      adjusted <- as.matrix(Matrix::nearPD(D_hat, corr = FALSE)$mat)
      MASS::mvrnorm(k, mu = c(0, 0), Sigma = adjusted)
    }
  )
  random_effects <- matrix(random_effects, nrow = k, ncol = 2L)
  bootstrap_ipd <- ipd
  bootstrap_ipd$correct_classification <- NA_integer_

  for (i in seq_len(k)) {
    rows <- which(bootstrap_ipd$study_id == studies[i])
    fixed_matrix <- stats::model.matrix(
      stats::delete.response(stats::terms(fixed_formula)),
      bootstrap_ipd[rows, , drop = FALSE]
    )
    random_matrix <- stats::model.matrix(
      ~ group - 1,
      bootstrap_ipd[rows, , drop = FALSE]
    )
    if (ncol(fixed_matrix) != length(beta_hat) ||
        ncol(random_matrix) != ncol(random_effects)) {
      stop("Bootstrap design matrices do not match fitted coefficients.",
           call. = FALSE)
    }
    probability <- .logit_inv(as.numeric(
      fixed_matrix %*% beta_hat +
        random_matrix %*% random_effects[i, ]
    ))
    bootstrap_ipd$correct_classification[rows] <-
      stats::rbinom(length(rows), 1L, probability)
  }
  bootstrap_ipd$misclassification <- 1L - bootstrap_ipd$correct_classification
  bootstrap_table <- .ipd_to_2x2(bootstrap_ipd)
  bootstrap_model <- tryCatch(
    .fit_dta_glmm(bootstrap_ipd, nAGQ, fixed_formula),
    error = function(e) NULL
  )
  if (is.null(bootstrap_model)) return(NA_real_)
  tryCatch(
    weighted_hz_statistic(bootstrap_model, bootstrap_table, cc)$statistic,
    error = function(e) NA_real_
  )
}
