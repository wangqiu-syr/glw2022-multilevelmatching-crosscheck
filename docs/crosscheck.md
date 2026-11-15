# Cross-check with multilevelMatching (Yang et al., 2016)

Both simulations, 500 replications each. Sources: `crosscheck/results/` (CSV names below) and the scripts `crosscheck/R/04_crosscheck_all.R` and `crosscheck/R/05_summarize_crosscheck.R`.

## What is checked

Guo, Liu and Wang (2022) state that without the school-level covariates Z their methods reduce to those of Yang et al. (2016). Yang and Barkley's R package `multilevelMatching` (version 1.0.0, CRAN 2019) implements those estimators. The rebuild matches on X and Z together, so the package, given the same seven covariates, should compute the same CM and GPS estimates as the rebuild. This check asks whether it does, replication by replication.

For each replication of each simulation:

1. **Same sample.** The sample is regenerated from the main run's seed stream (master seed 20220000 + simulation, L'Ecuyer stream r) and checked against the main run's stored diagnostics (`results/sim<k>_diagnostics.csv`).
2. **Our estimates.** DIF, CM and GPS are computed exactly as in the main run (`estimate_all()`).
3. **Three package fits,** each after `set.seed(123)`, with `trimming = FALSE` and `M_matches = 1`:
   - `covariates`: Mahalanobis matching on the seven covariates, compared with our CM;
   - `multinom`: the package fits its own multinomial logit (`nnet::multinom` at its default settings) and matches on it, compared with our GPS;
   - `existing`: the package matches on our fitted GPS matrix, compared with our GPS.
4. **Recorded:** all estimates (package signs flipped to our tau(t, t') = E[Y(t) - Y(t')]), the package's Abadie–Imbens variances, and unit-level diagnostics: how many imputed potential outcomes differ from ours, ties, the fitted-probability difference for `multinom`.

Two input workarounds are applied; the package itself is not modified. On R 4.3 or later its `class(X)` check fails for any matrix, so `X` is passed with `oldClass(X) <- "matrix"` (values unchanged). `existing` requires each row of the GPS matrix to sum to exactly 1, so rows off by rounding were adjusted by at most 3.3e-16. Details are in `BRIDGE.md`.

## Verification: the regenerated samples are the main run's

SIM_VERIFICATION_TABLE

## Agreement

SIM_AGREEMENT_TABLE

## Simulation 1

SIM1_SECTION

## Simulation 2

SIM2_SECTION

## Why the fits differ when they do

WHY_SECTION

## What this does and does not show

LIMITS_SECTION
