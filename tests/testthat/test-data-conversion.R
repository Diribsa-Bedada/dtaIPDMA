test_that("2x2 tables expand to the correct number of participants", {
  dat <- scoff[1:2, ]
  ipd <- table2x2_to_ipd(dat)

  expect_equal(nrow(ipd), sum(dat$TP + dat$FN + dat$FP + dat$TN))
  expect_equal(sort(unique(ipd$study_id)), dat$study_id)
  expect_equal(ipd$misclassification, 1L - ipd$correct_classification)
})

test_that("invalid tables are rejected", {
  dat <- scoff
  expect_error(table2x2_to_ipd(dat[, -1]), "Missing required columns")

  dat$TP[1] <- -1
  expect_error(table2x2_to_ipd(dat), "nonnegative")

  dat <- scoff
  dat$FN[1] <- 1.5
  expect_error(table2x2_to_ipd(dat), "integers")
})

