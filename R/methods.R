#' @export
print.dtaIPDMA_test <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("HZ statistic:", formatC(x$statistic, digits = digits), "\n")
  cat("Bootstrap p-value:", base::format.pval(x$p_value, digits = digits), "\n")
  cat("Effective bootstrap samples:", x$B_effective, "of", x$B, "\n")
  
  if (isTRUE(x$parallel)) {
    cat("Parallel workers:", x$n_cores, "\n")
  }
  
  invisible(x)
}

#' @export
summary.dtaIPDMA_test <- function(object, ...) {
  result <- list(
    method = object$method,
    results = data.frame(
      statistic = object$statistic,
      p_value = object$p_value,
      requested_B = object$B,
      effective_B = object$B_effective,
      failed_B = object$failed_bootstraps
    ),
    weighted_mean = object$weighted_mean,
    weighted_covariance = object$weighted_covariance
  )
  class(result) <- "summary.dtaIPDMA_test"
  result
}

#' @export
print.summary.dtaIPDMA_test <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  print(x$results, row.names = FALSE)
  cat("\nWeighted mean of empirical Bayes effects:\n")
  print(x$weighted_mean)
  cat("\nWeighted covariance matrix:\n")
  print(x$weighted_covariance)
  invisible(x)
}

#' Plot Empirical Bayes Effects
#'
#' Produces marginal density or normal Q-Q plots for raw or weighted,
#' standardized empirical Bayes effects.
#'
#' @param x An object returned by `weighted_hz_test()`.
#' @param type Either `"density"` or `"qq"`.
#' @param weighted Plot weighted standardized effects when `TRUE`.
#' @param ... Additional arguments passed to the primary plotting function.
#' @return `x`, invisibly.
#' @export
plot_dta_effects <- function(x, type = c("density", "qq"), weighted = TRUE, ...) {
  if (!inherits(x, "dtaIPDMA_test")) {
    stop("x must be a dtaIPDMA_test object.", call. = FALSE)
  }
  type <- match.arg(type)
  values <- if (isTRUE(weighted)) x$standardized_effects else x$effects
  labels <- c("Sensitivity", "Specificity")
  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters))
  graphics::par(mfrow = c(1L, 2L))
  for (j in seq_len(2L)) {
    if (type == "qq") {
      stats::qqnorm(values[, j], main = paste("Q-Q Plot:", labels[j]), ...)
      stats::qqline(values[, j], col = "red", lwd = 2)
    } else {
      estimate <- stats::density(values[, j])
      grid <- seq(min(estimate$x), max(estimate$x), length.out = 200L)
      normal <- if (isTRUE(weighted)) {
        stats::dnorm(grid)
      } else {
        stats::dnorm(grid, mean(values[, j]), stats::sd(values[, j]))
      }
      graphics::plot(estimate, ylim = c(0, 1.1 * max(estimate$y, normal)),
                     main = labels[j], xlab = "Empirical Bayes effect", ...)
      graphics::lines(grid, normal, lty = 2, lwd = 2)
      graphics::rug(values[, j])
    }
  }
  invisible(x)
}
