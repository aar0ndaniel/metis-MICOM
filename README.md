# metis

MICOM measurement invariance and multigroup analysis for [seminr](https://github.com/sem-in-r/seminr) partial least squares structural equation models.

`metis` implements the Measurement Invariance of Composite Models (MICOM) procedure of Henseler, Ringle and Sarstedt (2016) on top of the `seminr` estimation engine, and provides thin wrappers for bootstrap PLS multigroup analysis (PLS-MGA) and an optional permutation multigroup test. Constructs may use correlation weights (`seminr::mode_A`, reflective-style measurement) or regression weights (`seminr::mode_B`, formative measurement). The Step 2 compositional invariance statistic uses deterministic composite sign alignment and drops inadmissible permutation re-estimations, and the result reports permutation admissibility explicitly.

Step 2 and Step 3 results are cross-verified against `cSEM::testMICOM`. `cSEM` is used only for validation and is never a runtime dependency.

## Installation

```r
# install.packages("remotes")
remotes::install_github("aar0ndaniel/metis-MICOM")
```

Runtime dependencies are `seminr` and base R.

## Usage

```r
library(metis)
library(seminr)

# Data and a two-group variable
data <- mobi
data$group <- rep(c("A", "B"), length.out = nrow(data))

# A seminr composite model
mm <- constructs(
  composite("Image",       multi_items("IMAG", 1:5)),
  composite("Expectation", multi_items("CUEX", 1:3)),
  composite("Quality",     multi_items("PERQ", 1:7))
)
sm <- relationships(
  paths(from = "Image",       to = "Expectation"),
  paths(from = "Expectation", to = "Quality")
)
model <- estimate_pls(data, mm, sm)

# MICOM: three-step measurement invariance
mi <- metis_micom(model, data, group_var = "group", permutations = 1000)
print(mi)
mi$step2          # compositional invariance, per construct
mi$admissibility  # permutation admissibility report

# Bootstrap PLS-MGA, cautioned by the MICOM result
mga <- metis_pls_mga(model, data, group_var = "group", nboot = 1000,
                     micom_result = mi)
print(mga)
```

Formative constructs use regression weights:

```r
mm <- constructs(
  composite("Image",       multi_items("IMAG", 1:5), weights = mode_B),
  composite("Expectation", multi_items("CUEX", 1:3), weights = mode_B)
)
```

## Scope

`metis` targets composite models (Mode A and Mode B). It does not target common-factor (consistent PLS) constructs, for which MICOM is not the appropriate invariance procedure.

## Reference

Henseler, J., Ringle, C. M., & Sarstedt, M. (2016). Testing measurement invariance of composites using partial least squares. *International Marketing Review*, 33(3), 405-431. https://doi.org/10.1108/IMR-09-2014-0304

## License

MIT. See `LICENSE`.
