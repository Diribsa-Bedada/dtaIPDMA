# dtaIPDMA

**Robust Statistical Methods for Individual Participant Data Meta-Analysis of Diagnostic Test Accuracy**

The initial release provides the Project I distributional-diagnostics methods. Future projects (flexible random-effects models, verification-bias correction, and imperfect reference standards) are planned but **not yet implemented**.


`dtaIPDMA` (version 0.1.0) implements the classical Henze–Zirkler (HZ) test and a
precision-weighted HZ test for empirical Bayes study effects from bivariate
diagnostic test accuracy meta-analysis.

The weighted procedure combines the fitted between-study covariance matrix
with study-specific within-study uncertainty. Its p-value is calibrated using
a parametric bootstrap under the fitted bivariate normal GLMM. The bootstrap
preserves the original study sizes and disease-group structure and generates
new outcomes under the fitted null model.

## Installation

```r
# After pushing this repository to GitHub:
# install.packages("remotes")
# remotes::install_github("YOUR_GITHUB_USERNAME/dtaIPDMA")
```

## Quick example

```r
library(dtaIPDMA)

dat <- scoff_data()
result <- dta_hz_test(dat, B = 499, nAGQ = 7, seed = 2026)
result$summary

print(result$weighted)
summary(result$weighted)
plot_dta_effects(result$weighted, type = "qq", weighted = TRUE)
```

On a Windows laptop with eight logical processors, a faster bootstrap run can
use six worker processes:

```r
result <- dta_hz_test(
  dat, B = 499, nAGQ = 5, seed = 2026,
  parallel = TRUE, n_cores = 6
)
```

The required input contains one row per study and the columns `study_id`,
`TP`, `FN`, `FP`, and `TN`. Raw integer counts must be supplied. A continuity
correction is used only for within-study variance calculations when a study has
at least one zero cell.

## Individual participant data with covariates

For IPD, use one row per participant with a study identifier, disease status
(`0/1`), index-test result (`0/1`), and optional covariates. Covariates are
centred within study and fitted with separate effects for diseased and
non-diseased participants.

```r
ipd_result <- ipd_hz_test(
  data = ipd_data,
  study_id = "study_id",
  disease = "disease",
  test = "test_result",
  covariates = c("age", "sex"),
  B = 499,
  nAGQ = 7,
  parallel = TRUE,
  n_cores = 6,
  seed = 2026
)

ipd_result$summary
```

## Methods currently implemented

- Classical HZ test applied to empirical Bayes bivariate study effects.
- Weighted HZ statistic incorporating study-specific precision.
- Parametric-bootstrap calibration of the weighted statistic.
- Optional parallel bootstrap processing on Windows, macOS, and Linux.
- Direct IPD analysis with within-study covariate centring.

The nonparametric bootstrap is intentionally not included.

## Development status and covariate caveat

This initial release includes a correction that fixes the covariate-adjusted parametric-bootstrap linear
predictor. The weighted statistic uses within-study variance approximations derived
from collapsed 2x2 cell counts, even when participant-level covariates are
included in the BGLMM. This is intentional and matches the original Project I
simulation formulation. Covariates are centred within study before fitting;
the bootstrap retains the original study membership, disease status and
covariate design, regenerates outcomes under the fitted bivariate-normal
null, and refits the same model. The original simulation script includes age
and sex, but numerical equivalence between the package and simulation
implementation remains to be verified. Other covariate configurations have
not been separately validated.

Run `devtools::test()` and `devtools::check()` in an R environment before
publishing a release. These checks have not been run in the preparation
environment.
