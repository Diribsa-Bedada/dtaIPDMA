# dtaIPDMA

**Robust Statistical Methods for Individual Participant Data Meta-Analysis of Diagnostic Test Accuracy**

`dtaIPDMA` is an R package for methodological analysis of diagnostic test accuracy (DTA) studies using individual participant data (IPD) and study-level diagnostic 2 × 2 tables. The package is designed to support assessment of distributional assumptions and, in future versions, extensions to flexible hierarchical models and methods for addressing verification bias.

**Version:** 0.1.0  
**Development status:** Research software; initial methodological release

## Motivation

Meta-analysis of diagnostic test accuracy studies commonly models sensitivity and specificity jointly using a bivariate generalized linear mixed model (GLMM). Conventional formulations assume that study-specific random effects follow a bivariate normal distribution. Departures from this assumption may affect model adequacy and motivate alternative random-effects specifications.

Assessing bivariate normality is challenging when random effects are unobserved and must be estimated from studies of differing sizes and precision. In particular, empirical Bayes estimates have heterogeneous uncertainty and shrinkage. Applying a conventional multivariate normality test directly to these estimates may therefore yield misleading evidence about the underlying random-effects distribution.

`dtaIPDMA` provides a precision-weighted Henze–Zirkler (HZ) testing procedure, with inference calibrated using a parametric bootstrap that refits the diagnostic accuracy GLMM. This approach is intended to account for differences in the precision of estimated study effects and uncertainty introduced by model fitting.

## Implemented methods

Version 0.1.0 includes:

- **Precision-weighted HZ test for study-level DTA data:** `weighted_hz_test()` accepts one diagnostic 2 × 2 table per study (true positives, false negatives, false positives and true negatives).
- **Precision-weighted HZ test for IPD:** `ipd_hz_test()` accepts participant-level disease and index-test results, study identifiers and optional participant-level covariates.
- **Classical HZ test:** `classical_hz_test()` provides a conventional comparison based on fitted model effects.
- **Data conversion:** `table2x2_to_ipd()` converts aggregated diagnostic 2 × 2 tables to participant-level records for supported workflows.
- **Example datasets:** `scoff` and `hcc` are included for illustrating DTA meta-analysis workflows.

The precision-weighted procedure uses a parametric bootstrap for hypothesis testing. Bootstrap replication counts and quadrature settings can be adjusted to balance computational demands and numerical accuracy.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("Diribsa-Bedada/dtaIPDMA")
```

Load the package:

```r
library(dtaIPDMA)
```

> **Note:** The package is currently distributed through GitHub; it is not described here as a CRAN release.

## Examples

### Study-level diagnostic accuracy data

The package includes a study-level SCOFF dataset:

```r
library(dtaIPDMA)

data("scoff", package = "dtaIPDMA")

set.seed(2026)
result <- weighted_hz_test(
  data = scoff,
  B = 499L,
  nAGQ = 7L,
  seed = 2026
)

print(result)
```

The input should contain one row per study with the diagnostic classification counts required by the function. Consult `?weighted_hz_test` for the exact column requirements and return-value definitions.

### Individual participant data

For participant-level data, identify the study, reference-standard disease status and index-test result. Optional covariates may be supplied by name.

```r
# Replace `my_ipd` and column names with those in your dataset.
result_ipd <- ipd_hz_test(
  data = my_ipd,
  study_id = "study_id",
  disease = "disease",
  test = "test",
  covariates = c("age", "sex"),
  B = 499L,
  nAGQ = 7L,
  seed = 2026
)

print(result_ipd)
```

This second example is a template and requires a user-supplied dataset. Check `?ipd_hz_test` for coding requirements, particularly for binary outcomes and categorical covariates.

## Statistical considerations

**Null hypothesis.** The testing framework evaluates compatibility of the study-level random-effects distribution with the bivariate normal assumption under the fitted DTA GLMM.

**Bootstrap inference.** Because study effects are estimated rather than observed, and their uncertainty differs across studies, the weighted test uses a model-based parametric bootstrap with GLMM refitting rather than relying solely on the classical HZ reference distribution.

**Interpretation.** A small bootstrap p-value provides evidence against the assumed bivariate normal random-effects distribution. A non-significant result does not establish normality. Test performance depends on the number of studies, data configuration, model specification and numerical estimation.

**Computation.** Bootstrap refitting can be computationally intensive, particularly for larger meta-analyses or models with covariates. For substantive analyses, select a sufficiently large `B`, inspect convergence diagnostics and assess sensitivity to numerical integration settings.

## Scope and development roadmap

The current release focuses on normality assessment for bivariate DTA random effects. Planned research directions include:

- Flexible bivariate random-effects distributions, including skew-normal specifications.
- Likelihood-based approaches for selective reference-standard verification in IPD DTA meta-analysis.
- Further model diagnostics and robustness assessments.

These directions are **not implemented features of version 0.1.0**.

## Reproducibility and documentation

Function-level documentation is available through the R help system:

```r
?weighted_hz_test
?ipd_hz_test
?classical_hz_test
?table2x2_to_ipd
```

The GitHub repository contains the package source code and associated tests. Results should be reported together with the package version, bootstrap replication count, random seed, GLMM specification and numerical integration settings.

## Citation

If you use `dtaIPDMA` in research, please cite the package repository and identify the version used. A formal software citation and methodological manuscript reference will be added when available.

**Repository:** https://github.com/Diribsa-Bedada/dtaIPDMA

## License

A software license has not yet been specified in this README. Users should consult the repository for current licensing information before redistribution or reuse.

## Contact and contributions

Questions, reproducible bug reports and methodological feedback can be submitted through the [GitHub issue tracker](https://github.com/Diribsa-Bedada/dtaIPDMA/issues).

---

*`dtaIPDMA` is research software under active development. Users are encouraged to validate numerical results and report any unexpected behavior.*
