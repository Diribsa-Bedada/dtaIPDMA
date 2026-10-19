#' Expand Diagnostic Accuracy Tables to Participant-Level Rows
#'
#' Converts one row per study containing a diagnostic 2-by-2 table into the
#' participant-level representation used to fit the bivariate GLMM. The
#' expansion is exact; continuity-corrected counts must not be supplied.
#'
#' @param data A data frame with `study_id`, `TP`, `FN`, `FP`, and `TN`.
#' @return A data frame with one row per participant.
#' @export
table2x2_to_ipd <- function(data) {
  data <- .validate_dta_table(data)

  pieces <- lapply(seq_len(nrow(data)), function(i) {
    counts <- unlist(data[i, c("TP", "FN", "FP", "TN")], use.names = TRUE)
    study <- data$study_id[i]
    definitions <- data.frame(
      cell = c("TP", "FN", "FP", "TN"),
      disease = c(1L, 1L, 0L, 0L),
      correct = c(1L, 0L, 0L, 1L)
    )
    rows <- lapply(seq_len(nrow(definitions)), function(j) {
      n <- as.integer(counts[definitions$cell[j]])
      if (n == 0L) return(NULL)
      data.frame(
        study_id = rep(study, n),
        disease = rep(definitions$disease[j], n),
        correct_classification = rep(definitions$correct[j], n)
      )
    })
    do.call(rbind, rows)
  })

  ipd <- do.call(rbind, pieces)
  ipd$group <- factor(
    ifelse(ipd$disease == 1L, "diseased", "non-diseased"),
    levels = c("diseased", "non-diseased")
  )
  ipd$misclassification <- 1L - ipd$correct_classification
  rownames(ipd) <- NULL
  ipd[, c("study_id", "disease", "group", "correct_classification",
          "misclassification")]
}

.ipd_to_2x2 <- function(ipd) {
  studies <- unique(ipd$study_id)
  result <- lapply(studies, function(study) {
    x <- ipd[ipd$study_id == study, , drop = FALSE]
    data.frame(
      study_id = study,
      TP = sum(x$disease == 1L & x$correct_classification == 1L),
      FN = sum(x$disease == 1L & x$correct_classification == 0L),
      FP = sum(x$disease == 0L & x$correct_classification == 0L),
      TN = sum(x$disease == 0L & x$correct_classification == 1L)
    )
  })
  do.call(rbind, result)
}

.prepare_ipd_data <- function(data, study_id, disease, test, covariates = NULL) {
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  required <- c(study_id, disease, test, covariates)
  if (!is.character(required) || any(lengths(list(study_id, disease, test)) != 1L)) {
    stop("study_id, disease, and test must each be one column name.", call. = FALSE)
  }
  if (length(covariates) && !is.character(covariates)) {
    stop("covariates must be NULL or a character vector of column names.", call. = FALSE)
  }
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }
  if (anyNA(data[, required, drop = FALSE])) {
    stop("IPD input cannot contain missing values in analysis columns.", call. = FALSE)
  }

  disease_value <- data[[disease]]
  test_value <- data[[test]]
  if (!is.numeric(disease_value) || !is.numeric(test_value) ||
      !all(disease_value %in% c(0, 1)) || !all(test_value %in% c(0, 1))) {
    stop("disease and test must be numeric binary variables coded 0 and 1.",
         call. = FALSE)
  }
  if (length(unique(data[[study_id]])) < 2L) {
    stop("At least two studies are required.", call. = FALSE)
  }

  ipd <- data.frame(
    study_id = data[[study_id]],
    disease = as.integer(disease_value),
    correct_classification = ifelse(disease_value == 1, test_value, 1 - test_value)
  )
  ipd$group <- factor(
    ifelse(ipd$disease == 1L, "diseased", "non-diseased"),
    levels = c("diseased", "non-diseased")
  )
  ipd$misclassification <- 1L - ipd$correct_classification

  centered_covariates <- character()
  for (j in seq_along(covariates)) {
    value <- data[[covariates[j]]]
    if (is.logical(value)) value <- as.integer(value)
    if (is.factor(value) || is.character(value)) {
      level <- unique(as.character(value))
      if (length(level) != 2L) {
        stop("Categorical covariates must have exactly two levels.", call. = FALSE)
      }
      value <- as.integer(as.character(value) == level[2L])
    }
    if (!is.numeric(value) || any(!is.finite(value))) {
      stop("Covariates must be finite numeric variables or binary factors.",
           call. = FALSE)
    }
    name <- paste0(".whz_cov_", j)
    ipd[[name]] <- value - stats::ave(value, ipd$study_id, FUN = mean)
    centered_covariates <- c(centered_covariates, name)
  }
  fixed_formula <- if (length(centered_covariates)) {
    stats::as.formula(paste(
      "cbind(correct_classification, misclassification) ~ group +",
      paste(paste0("group:", centered_covariates), collapse = " + "), "- 1"
    ))
  } else {
    cbind(correct_classification, misclassification) ~ group - 1
  }
  list(ipd = ipd, fixed_formula = fixed_formula, covariates = covariates)
}
