#' SCOFF Diagnostic Accuracy Data
#'
#' Study-level 2-by-2 diagnostic accuracy data from nine
#' studies evaluating the SCOFF screening questionnaire.
#'
#' @format A data frame with nine rows and five variables:
#' \describe{
#'   \item{study_id}{Study identifier.}
#'   \item{TP}{Number of true positives.}
#'   \item{FN}{Number of false negatives.}
#'   \item{FP}{Number of false positives.}
#'   \item{TN}{Number of true negatives.}
#' }
#'
#' @details
#' This dataset is provided to illustrate diagnostic test
#' accuracy meta-analysis and statistical methods for
#' assessing the bivariate normality assumption of
#' study-specific random effects.
#'
#' @source
#' Botella, J., Huang, H., and Suero, M. (2015).
#' Meta-analysis of the accuracy of tools used for binary
#' classification when the primary studies employ different
#' references.
#' \emph{Psychological Methods}, 20(3), 331.
#'
#' @examples
#' data("scoff", package = "dtaIPDMA")
#' head(scoff)
#'
"scoff"