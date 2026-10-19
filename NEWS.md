# dtaIPDMA 0.1.0 (development)

* Adopted package title: "Robust Statistical Methods for Individual Participant Data Meta-Analysis of Diagnostic Test Accuracy".
* Updated citation and README to reflect the broader package scope and the
  original Project I covariate-adjusted simulation formulation.
* Retained the existing Project I implementation; no changes to test statistics
  or bootstrap algorithms in this metadata update.

* Initial package under a broader name for diagnostic test accuracy IPD meta-analysis.
* Includes Project I weighted HZ test, bootstrap calibration, and covariate support.
* Renamed package namespace, S3 result classes, tests, help files, and vignette.
* Previous development history was maintained under the former wHZDTA name.

# wHZDTA 0.3.1 (development)

- Corrected parametric-bootstrap linear predictor for covariate-adjusted IPD:
  fixed effects use the full fixed-effects design matrix, while study-level
  bivariate random intercepts use a separate two-column group design matrix.
- Added targeted regression tests for covariate-adjusted bootstrap generation.
- Clarified that the covariate-adjusted weighting scheme still needs
  methodological validation before being claimed as validated.

# wHZDTA 0.3.0

- Added `ipd_hz_test()` for direct individual participant data analysis with
  within-study-centred covariates and parallel parametric bootstrap.

# wHZDTA 0.2.0

- Added Windows-compatible parallel parametric bootstrap processing through
  `parallel = TRUE` and `n_cores`.

# wHZDTA 0.1.1

- Corrected the namespace for normal Q-Q plotting functions.

# wHZDTA 0.1.0

- Initial implementation of the weighted HZ statistic.
- Parametric-bootstrap p-value calibration under the fitted BVN GLMM.
- Classical HZ comparison.
- Diagnostic plots and SCOFF example data.
