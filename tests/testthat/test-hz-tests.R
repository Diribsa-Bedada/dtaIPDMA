test_that("weighted HZ test returns the documented structure", {
  skip_if_not_installed("GLMMadaptive")
  skip_if_not_installed("mvnTest")

  fit <- weighted_hz_test(scoff, B = 1, nAGQ = 3, seed = 17)
  expect_s3_class(fit, "dtaIPDMA_test")
  expect_true(is.finite(fit$statistic))
  expect_true(is.na(fit$p_value) || fit$p_value >= 0 && fit$p_value <= 1)
  expect_equal(fit$B, 1L)
})

test_that("combined interface reports two methods", {
  skip_if_not_installed("GLMMadaptive")
  skip_if_not_installed("mvnTest")

  fit <- dta_hz_test(scoff, B = 1, nAGQ = 3, seed = 19)
  expect_equal(nrow(fit$summary), 2L)
  expect_equal(
    fit$summary$method,
    c("Classical HZ", "Weighted HZ (parametric bootstrap)")
  )
})

test_that("IPD interface accepts covariates", {
  skip_if_not_installed("GLMMadaptive")
  skip_if_not_installed("mvnTest")

  ipd <- table2x2_to_ipd(scoff)
  ipd$test_result <- ifelse(
    ipd$disease == 1L,
    ipd$correct_classification,
    1L - ipd$correct_classification
  )
  ipd$age <- rep(seq(40, 60, length.out = 10), length.out = nrow(ipd))
  ipd$sex <- factor(rep(c("F", "M"), length.out = nrow(ipd)))

  fit <- ipd_hz_test(
    ipd,
    study_id = "study_id",
    disease = "disease",
    test = "test_result",
    covariates = c("age", "sex"),
    B = 1,
    nAGQ = 3,
    seed = 31
  )
  expect_equal(nrow(fit$summary), 2L)
  expect_equal(fit$weighted$covariates, c("age", "sex"))
})


test_that("covariate-adjusted bootstrap uses separate random design", {
  skip_if_not_installed("GLMMadaptive")
  ipd <- table2x2_to_ipd(scoff)
  ipd$test_result <- ifelse(ipd$disease == 1L,
                            ipd$correct_classification,
                            1L - ipd$correct_classification)
  ipd$age <- rep(seq(40, 60, length.out = 10), length.out = nrow(ipd))
  ipd$sex <- factor(rep(c("F", "M"), length.out = nrow(ipd)))
  prepared <- dtaIPDMA:::.prepare_ipd_data(
    ipd, "study_id", "disease", "test_result", c("age", "sex")
  )
  fit <- dtaIPDMA:::.fit_dta_glmm(prepared$ipd, nAGQ = 3L,
                                fixed_formula = prepared$fixed_formula)
  beta_hat <- as.numeric(nlme::fixef(fit))
  expect_gt(length(beta_hat), 2L)
  expect_equal(ncol(model.matrix(
    delete.response(terms(prepared$fixed_formula)), prepared$ipd)),
    length(beta_hat))
  expect_equal(ncol(model.matrix(~ group - 1, prepared$ipd)), 2L)
  set.seed(117)
  stat <- dtaIPDMA:::.bootstrap_hz_once(
    prepared$ipd, unique(prepared$ipd$study_id), beta_hat,
    fit$D, nAGQ = 3L, cc = 0.5,
    fixed_formula = prepared$fixed_formula
  )
  expect_true(is.numeric(stat) && length(stat) == 1L)
  expect_true(is.na(stat) || is.finite(stat))
})
