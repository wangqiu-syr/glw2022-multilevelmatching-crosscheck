# Cross-check with multilevelMatching (Yang et al., 2016)

Both simulations, 500 replications each. Sources: the CSVs in `crosscheck/results/` named below, and the scripts `crosscheck/R/04_crosscheck_all.R` (all replications), `crosscheck/R/05_summarize_crosscheck.R` (summaries) and `crosscheck/R/06_existing_mismatches.R` (unit-level check of the `existing` differences).

## What is checked

Guo, Liu and Wang (2022) state that without the school-level covariates Z their methods reduce to those of Yang et al. (2016). Yang and Barkley's R package `multilevelMatching` (version 1.0.0, CRAN 2019) implements those estimators. The rebuild matches on X and Z together, so the package, given the same seven covariates, should compute the same CM and GPS estimates as the rebuild. This check asks whether it does, replication by replication.

For each replication of each simulation:

1. **Same sample.** The sample is regenerated from the main run's seed stream (master seed 20220000 + simulation, L'Ecuyer stream r) and checked against the main run's stored diagnostics (`results/sim<k>_diagnostics.csv`).
2. **Our estimates.** DIF, CM and GPS are computed exactly as in the main run (`estimate_all()`).
3. **Three package fits,** each after `set.seed(123)`, with `trimming = FALSE` and `M_matches = 1`:
   - `covariates`: Mahalanobis matching on the seven covariates, compared with our CM;
   - `existing`: the package matches on our fitted GPS matrix, compared with our GPS;
   - `multinom`: the package fits its own multinomial logit (`nnet::multinom` at its default settings) and matches on it, also compared with our GPS.
4. **Recorded:** all estimates (package signs flipped to our tau(t, t') = E[Y(t) - Y(t')]), the package's Abadie–Imbens variances, and unit-level diagnostics: how many imputed potential outcomes differ from ours, ties and near-ties in our matching, and the fitted-probability difference for `multinom`.

Two input workarounds are applied; the package itself is not modified. On R 4.3 or later its `class(X)` check fails for any matrix, so `X` is passed with `oldClass(X) <- "matrix"` (values unchanged). `existing` requires each row of the GPS matrix to sum to exactly 1, so rows off by rounding (on average 75 of 400 rows in Simulation 1 and 602 of 1,800 in Simulation 2) were adjusted by at most 4.2e-16. Details are in `BRIDGE.md`.

The cross-check ran on Linux with R 4.3.3, Matching 4.10-14 and nnet 7.3-19; the main run used R 4.5.1 on Windows with nnet 7.3.20. Run time on 4 workers: about 1 minute for Simulation 1 and 18 minutes for Simulation 2.

## Verification: the regenerated samples are the main run's

From `crosscheck/results/crosscheck_verification.csv`.

| | Simulation 1 | Simulation 2 |
|---|---|---|
| Replications | 500 | 500 |
| Stored diagnostics (pool size, CM and GPS matching distances), largest difference | 2.9e-9 | 4.8e-9 |
| Mean of the regenerated DIF estimates vs the main run's stored means, largest difference | 9.6e-14 | 2.0e-13 |
| Same for CM | 1.0e-13 | 9.2e-14 |
| Same for GPS | 2.1e-13 | **5.6e-5** |
| Exact ties in our own matching | 0 | 0 |
| Package errors / warnings | 0 / 0 | 0 / 0 |

Every regenerated sample is the main run's: pool sizes are identical and the CM matching distances agree to rounding. Our DIF and CM estimates reproduce the main run's in every replication of both simulations, and our GPS estimates do so in Simulation 1.

In Simulation 2 our GPS estimates average to the main run's to within 5.6e-5, not to rounding. The GPS matching distances differ by up to 4.8e-9 because the fitted probabilities differ slightly across the two setups (R 4.3.3 / nnet 7.3-19 / Linux against R 4.5.1 / nnet 7.3.20 / Windows), even at reltol = 1e-14. Simulation 2 is full of near-ties: the smallest gap between a unit's best and second-best GPS match is below 1e-8 in 151 of 500 replications (median 1.9e-8). A few such near-ties are decided differently on the two setups, which moves some replications' GPS estimates. The main run's per-replication estimates were not saved, so the affected replications cannot be listed. The effect on the published summaries is at most 5.6e-5 in a mean estimate, against GPS biases of about 0.1. It does not affect the cross-check, which compares the package with our estimates on the same setup and the same samples.

## Agreement

From `crosscheck/results/crosscheck_agreement.csv`. "Identical" means every imputed potential outcome is the same as ours (400 x 4 = 1,600 per replication in Simulation 1; 1,800 x 9 = 16,200 in Simulation 2).

| Package fit | Compared with | Sim. 1: identical replications | Sim. 1: largest difference in an estimate | Sim. 2: identical replications | Sim. 2: largest difference in an estimate |
|---|---|---|---|---|---|
| `covariates` | our CM | 500 of 500 | 0 | 500 of 500 | 1e-13 (rounding) |
| `existing` | our GPS | 499 of 500 | 0.014 | 476 of 500 | 0.019 |
| `multinom` | our GPS | 405 of 500 | 0.25 | 0 of 500 | 2.72 |

- **`covariates` reproduces our CM in every replication.**
- **`existing` reproduces our GPS except at near-ties,** 30 imputed outcomes in 25 of 1,000 replications. In each of the 30 our match is the nearest candidate, and the package's is behind by less than 1e-10 in squared distance, which the package treats as a tie (see below).
- **`multinom` differs because of how the package fits the GPS model,** not because of how it matches. In Simulation 1, 95 of 500 replications differ, in at most 20 of 1,600 imputed outcomes. In Simulation 2 every replication differs, in on average 1,003 of 16,200 (6.2%).

## Simulation 1

Absolute bias against τ\_pop by pair, from `crosscheck/results/sim1_crosscheck_summary.csv` (τ\_bal is in the same file). The paper's values are from its Table 1.

| Pair | Ours CM | Pkg CM | Ours GPS | Pkg GPS (existing) | Pkg GPS (multinom) | Paper CM | Paper GPS |
|---|---|---|---|---|---|---|---|
| (1,2) | 0.110 | 0.110 | 0.024 | 0.024 | 0.024 | 0.026 | 0.060 |
| (1,3) | 0.345 | 0.345 | 0.037 | 0.037 | 0.036 | 0.166 | 0.374 |
| (1,4) | 0.199 | 0.199 | 0.048 | 0.048 | 0.048 | 0.240 | 0.481 |
| (2,3) | 0.235 | 0.235 | 0.013 | 0.013 | 0.012 | 0.140 | 0.315 |
| (2,4) | 0.088 | 0.088 | 0.024 | 0.024 | 0.024 | 0.213 | 0.421 |
| (3,4) | 0.146 | 0.146 | 0.011 | 0.011 | 0.012 | 0.074 | 0.107 |
| Mean | 0.187 | 0.187 | 0.026 | 0.026 | 0.026 | 0.143 | 0.293 |

Coverage (%) of 95% intervals against τ\_pop. Ours are the main run's percentile bootstrap intervals (B = 200, D8). The package's are estimate ± 1.96 SE with its Abadie–Imbens (2006) variance, and for `multinom` also its Abadie–Imbens (2016) variance, which accounts for estimating the GPS.

| Pair | Ours CM (bootstrap) | Pkg CM (AI 2006) | Ours GPS (bootstrap) | Pkg GPS (AI 2006) | Pkg GPS, multinom (AI 2016) | Paper CM | Paper GPS |
|---|---|---|---|---|---|---|---|
| (1,2) | 91.6 | 94.8 | 100 | 99.2 | 93.6 | 92 | 100 |
| (1,3) | 83.4 | 89.0 | 100 | 99.6 | 96.6 | 92 | 99 |
| (1,4) | 93.6 | 96.0 | 100 | 99.4 | 92.7 | 96 | 100 |
| (2,3) | 91.4 | 95.6 | 100 | 99.8 | 97.8 | 95 | 99 |
| (2,4) | 95.4 | 97.2 | 100 | 100.0 | 97.8 | 97 | 100 |
| (3,4) | 93.2 | 96.4 | 100 | 99.8 | 96.8 | 98 | 100 |
| Mean | 91.4 | 94.8 | 100 | 99.6 | 95.9 | 95 | 99.7 |

The "Pkg GPS (AI 2006)" column is the `existing` fit; the `multinom` fit gives the same coverage. Six of the 3,000 AI 2016 variances are negative (five in pair (1,4), one in (2,4)); those replications have no interval and are left out of that column.

## Simulation 2

Means over the 36 pairs, from `crosscheck/results/sim2_crosscheck_summary.csv` and `crosscheck/results/crosscheck_overview.csv`. Pair-by-pair values are in the CSV.

| Estimator | Mean abs. bias, τ\_pop | Mean abs. bias, τ\_bal | Mean MSE, τ\_pop | Mean coverage, τ\_pop (%) |
|---|---|---|---|---|
| Ours CM | 0.295 | 0.262 | 0.288 | 78.4 (bootstrap) |
| Pkg CM | 0.295 | 0.262 | 0.288 | 90.1 (AI 2006) |
| Ours GPS | 0.099 | 0.096 | 0.626 | 100.0 (bootstrap) |
| Pkg GPS (existing) | 0.099 | 0.096 | 0.626 | 99.6 (AI 2006) |
| Pkg GPS (multinom) | 0.100 | 0.092 | 0.627 | 99.6 (AI 2006); 94.7 (AI 2016) |

The paper's Table 3 gives mean absolute biases of 0.139 (CM) and 0.307 (GPS); its signs and its target are discussed in `docs/spec.md` and `docs/table3_checks.md`. Three of the 18,000 AI 2016 variances are negative and are left out.

Coverage by pair ranges from 39.6% to 93.4% for our bootstrap CM intervals and from 62.6% to 97.8% for the package's AI 2006 CM intervals; for GPS, from 87.6% to 97.8% with the package's AI 2016 variance.

## Why the fits differ when they do

**`covariates`: no differences.** The package calls `Matching::Match` with `Weight = 2`, which standardizes each covariate and weights by the inverse of their correlation matrix over the whole sample. That is Mahalanobis distance with the pooled sample covariance up to a constant factor, the same metric as our CM (D6), with one nearest neighbour, with replacement. Mahalanobis distances on seven covariates have no near-ties at the scale that matters here.

**`existing`: identical up to the package's notion of a tie.** The package imputes Y(s) by matching on column s of the GPS matrix, so the t' term uses p_hat_t'. This is our D7 reading. The 30 differing imputed outcomes are all near-ties (`crosscheck/results/existing_mismatches.csv`). `Matching::Match` treats candidates whose squared standardized distances differ by less than about 1e-10 as tied, even with `distance.tolerance = 0` and `tolerance = 1e-20`. This was confirmed directly on replication 351 of Simulation 1 and with constructed examples. The package uses `ties = FALSE`, so `Match` then picks one of the tied candidates at random. Our code picks the exact nearest neighbour and breaks only exact ties at random, and no exact ties occurred. In all 30 cases ours is the nearest candidate; the package's choice is farther by at most 2.3e-7 in |p_hat difference| and 9.3e-11 in squared standardized distance. The largest effect on an estimate is 0.019.

**`multinom`: the GPS model is fitted at nnet's default tolerance.** The package fits `nnet::multinom` with its defaults (maxit = 100, reltol = 1e-8); the rebuild uses reltol = 1e-14 and maxit = 10000 (D5). The resulting fitted probabilities differ from ours by up to 1.7e-4 (median per replication 1.2e-6) in Simulation 1 and up to 2.2e-3 (median 2.7e-4) in Simulation 2. In 2 of the 500 Simulation 2 replications the package's fit stops at the iteration limit without converging. GPS matching distances in Simulation 2 are of the same order as these differences, so 6.2% of imputed outcomes change in every replication, and single estimates move by up to 2.7. This is the problem described under D5 and in the README's Lesson, now seen with the package as shipped.

Averaged over replications, the `multinom` differences are small:

- **Simulation 1:** the mean package-minus-ours difference is at most 0.0007 in any pair, with |z| < 1.7 in all 6 pairs.
- **Simulation 2:** at most 0.027, with |z| > 2 in 9 of 36 pairs and none above 3 (`crosscheck/results/sim2_crosscheck_multinom_minus_ours.csv`). The pairs share treatment levels, so these z-values are not independent.
- **Bias barely changes:** mean absolute bias against τ\_pop is 0.0995 against our 0.0991.

**Intervals.** Our bootstrap intervals and the package's analytic ones answer the same question differently, and coverage differs accordingly:

- **CM:** the AI 2006 variance gives coverage closer to 95% than our percentile bootstrap (94.8% against 91.4% in Simulation 1; 90.1% against 78.4% in Simulation 2).
- **GPS:** the AI 2006 variance, which treats the GPS as known, gives near-100% coverage, as our bootstrap does. The AI 2016 variance, which accounts for estimating the GPS, gives 95.9% and 94.7% with intervals about 30% narrower.

The bootstrap is not formally valid for matching estimators (Abadie and Imbens, 2008; see D8), so this is consistent with theory. The AI 2016 variance is computed only for the package's own `multinom` fit.

## What this does and does not show

- **Shows:** given the same seven covariates, Yang et al.'s package computes the rebuild's CM estimates in all 1,000 replications. Given the same fitted GPS, it computes the rebuild's GPS estimates except for 30 near-ties that its matching routine treats as ties. D7 (the t' term matched on p_hat_t') is the estimator as the package implements it.
- **Shows:** fitting the GPS model at nnet's default tolerance, as the package does, changes individual GPS estimates noticeably in Simulation 2, though not the average bias. Anyone using the package on data like Simulation 2 should expect replication-level results to depend on the optimizer tolerance.
- **Shows:** the rebuild's own GPS estimates in Simulation 2 are reproducible across setups only up to near-ties, with effects of at most 5.6e-5 on the summary means.
- **Does not show** what Guo, Liu and Wang's own code did; the paper reports "Code Availability N/A".
- **Does not test** anything specific to clustering: like the rebuild, the package treats z1–z3 as ordinary covariates.
- **Version caveat:** multilevelMatching 1.0.0 does not run unmodified on R 4.3 or later; the input workaround is described in `BRIDGE.md`.

## References

- Abadie, A., & Imbens, G. W. (2006). Large sample properties of matching estimators for average treatment effects. *Econometrica, 74*(1), 235–267.
- Abadie, A., & Imbens, G. W. (2008). On the failure of the bootstrap for matching estimators. *Econometrica, 76*(6), 1537–1557.
- Abadie, A., & Imbens, G. W. (2016). Matching on the estimated propensity score. *Econometrica, 84*(2), 781–807.
- Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19: Multilevel covariates matching and propensity score matching. *Annals of Data Science, 9*(5), 967–982.
- Yang, S., & Barkley, B. G. (2019). *multilevelMatching* (R package version 1.0.0).
- Yang, S., Imbens, G., Cui, Z., Faries, D., & Kadziola, Z. (2016). Propensity score matching and subclassification in observational studies with multi-level treatments. *Biometrics, 72*, 1055–1065.
