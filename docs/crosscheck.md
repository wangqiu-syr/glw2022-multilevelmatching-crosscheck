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

   The estimator and the three GPS versions are written out in the README's Background ("The GPS estimator in the cross-check"); [crosscheck_math.md](crosscheck_math.md) derives step by step why all three compute it and when they can differ.
4. **Recorded:** all estimates (package signs flipped to our tau(t, t') = E[Y(t) - Y(t')]), the package's Abadie–Imbens variances, and unit-level diagnostics: how many imputed potential outcomes differ from ours, ties and near-ties in our matching, and the fitted-probability difference for `multinom`.

Two input workarounds are applied; the package itself is not modified. On R 4.3 or later its `class(X)` check fails for any matrix, so `X` is passed with `oldClass(X) <- "matrix"` (values unchanged). `existing` requires each row of the GPS matrix to sum to exactly 1, so rows off by rounding (on average 75 of 400 rows in Simulation 1 and 602 of 1,800 in Simulation 2) were adjusted by at most 4.2e-16. Details are in Appendix B.

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

### Summary across the three fits

This summary was in the main README until 2026-11-17 and moved here unchanged, apart from link paths.

| Package fit | Matches on | Compared with | Replications with every matched outcome identical, Sim. 1 / Sim. 2 | Matched outcomes that differ, average per replication, Sim. 1 / Sim. 2 | Mean absolute bias vs τ\_pop, Sim. 1 / Sim. 2: package (ours) |
|---|---|---|---|---|---|
| `covariates` | the 7 covariates (Mahalanobis) | our CM | 500 / 500 of 500 | 0 / 0 | 0.187 / 0.295 (0.187 / 0.295) |
| `existing` | our fitted GPS, p̂ | our GPS | 499 / 476 of 500 | 0.002 / 0.06 | 0.026 / 0.099 (0.026 / 0.099) |
| `multinom` | its own fitted GPS, p̃ | our GPS | 405 / 0 of 500 | 0.4 of 1,600 / 1,003 of 16,200 (6.2%) | 0.026 / 0.100 (0.026 / 0.099) |

A matched outcome is one imputed potential outcome Ŷ\_i(s); a replication has 400 x 4 = 1,600 of them in Simulation 1 and 1,800 x 9 = 16,200 in Simulation 2. "0 of 500" for `multinom` in Simulation 2 means that every replication has at least one differing matched outcome (6.2% on average), not that the estimates disagree broadly: the bias is the same to two decimals. The three fits are defined in [rebuild.md](rebuild.md#background) ("The GPS estimator in the cross-check").

- **The package computes the rebuild's estimators.** Covariate matching agrees in every replication. Given our fitted GPS, GPS matching agrees except for 30 imputed outcomes in 25 of 1,000 replications. All 30 are near-ties that the package's matching routine (`Matching::Match`) treats as ties and breaks at random; in each, ours is the nearest candidate. [docs/crosscheck_math.md](crosscheck_math.md) shows, from the code of both, why the three GPS versions compute the same estimator.
- **The GPS model's fitting tolerance matters for single estimates (D5).** The package fits the GPS model at nnet's default tolerance. In Simulation 2 that changes about 6% of the matched outcomes in every replication and moves single estimates by up to 2.7 ([side-by-side note, Figure 4](../crosscheck/side_by_side/README.md#agreement-replication-by-replication)), but not the bias or MSE.
- **Reproducibility across setups.** The regenerated samples are the main run's. Our DIF and CM estimates, and our GPS estimates in Simulation 1, reproduce the main run's exactly. In Simulation 2 a few GPS near-ties fall differently on this setup (Linux, R 4.3.3) than in the main run (Windows, R 4.5.1), moving GPS summary means by at most 5.6e-5.

<img src="crosscheck/side_by_side/figures/fig4_gps_per_replication.png" width="100%" alt="Every GPS estimate, the package's own GPS fit against the rebuild, Simulations 1 and 2">

**Figure J. Every GPS estimate, package against rebuild.** Each point is one pair in one replication (3,000 in Simulation 1, 18,000 in Simulation 2). In Simulation 2 the package's own GPS fit moves estimates off the diagonal in both directions; the averages, and hence bias and MSE, are unchanged.

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
- **Version caveat:** multilevelMatching 1.0.0 does not run unmodified on R 4.3 or later; the input workaround is described in Appendix B.

## Appendix A. The first check, and the full run as first recorded

Moved here from BRIDGE.md on 2026-11-17, unchanged. A.1 is the check of replication 1 of Simulation 1 that started the cross-check; A.2 is the summary written when all replications had run, which the sections above report in full.

### A.1 Replication 1 of Simulation 1

Same sample as the main run: master seed 20220001, L'Ecuyer stream 1. Per-replication estimates of the main run were not kept in git (`results/**/*.rds` is ignored), so the regenerated sample was checked against the stored rep-1 diagnostics in `results/sim1_diagnostics.csv`: pool size identical, CM matching distances within 4.7e-15, GPS matching distances within 2e-17. No ties occurred in our own matching.

Three `multiMatch()` fits, each after `set.seed(123)`, with `trimming = FALSE` and `M_matches = 1`. Package signs flipped to our convention tau(t, t') = E[Y(t) - Y(t')]:

| Pair | Ours CM | Pkg CM | Ours GPS | Pkg GPS (multinom) | Pkg GPS (existing) |
|---|---|---|---|---|---|
| (1,2) | -2.017504 | -2.017504 | -1.226368 | -1.226368 | -1.226368 |
| (1,3) | -1.374542 | -1.374542 | -0.961368 | -0.961368 | -0.961368 |
| (1,4) | -2.146265 | -2.146265 | -1.360152 | -1.360152 | -1.360152 |
| (2,3) | 0.642962 | 0.642962 | 0.265000 | 0.265000 | 0.265000 |
| (2,4) | -0.128762 | -0.128762 | -0.133784 | -0.133784 | -0.133784 |
| (3,4) | -0.771724 | -0.771724 | -0.398784 | -0.398784 | -0.398784 |

All 18 package-minus-ours differences are exactly 0, and all 1,600 imputed potential outcomes (400 units x 4 levels) are identical in all three fits.

- **`existing` reproduces our GPS exactly.** The package imputes Y(s) by matching on column s of the GPS matrix, so the t' term uses p_hat_t'. This is our D7 reading: D7 is the estimator as implemented in Yang et al.'s own package. It does not show what Guo, Liu & Wang's code did.
- **`multinom` agreement is not guaranteed in other replications.** The package fits `nnet::multinom` at its defaults (maxit 100, reltol 1e-8); ours uses reltol 1e-14 (D5). Fitted probabilities differ by up to 8.2e-7, while the smallest best-versus-second-best GPS match gap in this sample is 1.4e-7. No match changed here, but one could elsewhere, more likely in Simulation 2. (Confirmed by the full run below.)

### A.2 All replications of both simulations

Run on request on 2026-11-15 with `crosscheck/R/04_crosscheck_all.R` (500 replications each; about 1 min and 18 min on 4 workers), summarized by `05_summarize_crosscheck.R` and `06_existing_mismatches.R`. Full report: `docs/crosscheck.md`.

- **Same samples as the main run.** Every regenerated sample matches the stored diagnostics. Our DIF and CM estimates average to the main run's stored means to within 2e-13 in both simulations, and our GPS estimates do in Simulation 1. In Simulation 2, our GPS means differ from the main run's by up to 5.6e-5: near-ties are common there (best-vs-second-best gap below 1e-8 in 151 of 500 replications), and the fitted probabilities differ at about the 1e-9 level between this setup (Linux, R 4.3.3, nnet 7.3-19) and the main run's (Windows, R 4.5.1, nnet 7.3.20), so a few matches fall differently.
- **`covariates`:** identical imputations in 500 of 500 replications in both simulations.
- **`existing`:** identical in 499 of 500 (Simulation 1) and 476 of 500 (Simulation 2). The 30 differing imputed outcomes are all near-ties: ours is the nearest candidate in every case, and the package's choice is behind by at most 9.3e-11 in squared standardized distance. `Matching::Match` treats such candidates as tied, even with `distance.tolerance = 0`, and picks one at random (`ties = FALSE`). Largest effect on an estimate: 0.019.
- **`multinom`:** identical in 405 of 500 (Simulation 1) and 0 of 500 (Simulation 2). The package's default-tolerance fit differs from ours by up to 2.2e-3 in the probabilities in Simulation 2, changing about 6.2% of imputed outcomes per replication and single estimates by up to 2.7; 2 fits stopped at maxit = 100. Averaged over replications the shift is small (mean absolute bias 0.0995 against our 0.0991, τ_pop).
- **Coverage with the package's variances:** GPS with AI 2016 variance 95.9% (Sim. 1) and 94.7% (Sim. 2), against 100% for our bootstrap; CM with AI 2006 variance 94.8% and 90.1%, against 91.4% and 78.4%.

## Appendix B. multilevelMatching facts and workarounds

Moved here from BRIDGE.md on 2026-11-17, unchanged.

### B.1 multilevelMatching facts (version 1.0.0, CRAN 2019-05-08)

From `?multiMatch` and the package source (`R/prepareData.R`, `R/matchAllTreatments.R`, `R/estimateTau.R`):

- **Signature:** `multiMatch(Y, W, X, trimming = NULL, match_on, model_options = list(reference_level = sort(W)[1]), M_matches = 1, J_var_matches = 1)`.
- **Trimming:** default `NULL` = no trimming for `covariates` and `existing`; an error for `multinom` (must be set). `TRUE` is supported only for `multinom` (Crump et al. 2009 rule).
- **Matching:** `Matching::Match(..., Weight = 2, distance.tolerance = 0, ties = FALSE, M = M_matches)`, with replacement. `Weight = 2` is Mahalanobis on the full pooled sample, as in D6; there is no option to change the metric. Ties are broken at random, which is why a seed is set. `Match` treats candidates whose squared standardized distances differ by less than about 1e-10 as tied, even with `distance.tolerance = 0` and `tolerance = 1e-20` (confirmed on replication 351 of Simulation 1 and with constructed examples).
- **Output fields:** `results` (`Param`, `Trt1`, `Trt2`, `Estimate`, `Variance`, `VarianceAI2016`), `mu` (`Param`, `Trt`, `Estimate`), `analysis_idx`, `impute_mat`, `impute_mat_sorted`, `propensity_scores`, `model`, `estimate_args`, `impute_match_data`, plus `AI2016_args` for `multinom`. The help page's Value section lists older names (`tauestimate`, `varestimate`) that `multiMatch()` does not return.
- **Sign:** `Param` is `EY(Trt2)-EY(Trt1)` with Trt1 < Trt2, and `Estimate` = mu(Trt2) - mu(Trt1). Our tau(t, t') = -Estimate.
- **Treatment coding:** 0-based coding gives the same estimates with shifted labels. For `existing`, column k of the GPS matrix must belong to the k-th treatment level in sorted order.
- **GitHub master** (`shuyang1987/multilevelMatching`) has the same code as the CRAN release for every function `multiMatch()` uses.

### B.2 Workarounds needed (package not modified)

- **R >= 4.3:** `prepareData()` runs `(length(class(X)) != 1 && class(X) != "matrix") || ...`. A plain matrix has class `c("matrix", "array")`, and R 4.3 made that an error (`'length = 2' in coercion to 'logical(1)'`). It would fail on the main run's R 4.5.1 too. Fix: pass `X` after `oldClass(X) <- "matrix"`; values are unchanged.
- **`existing` needs exact row sums:** every row of the GPS matrix must satisfy `sum == 1` in floating point. In replication 1, 75 of our 400 rows were off by rounding (on average 75 of 400 in Simulation 1 and 602 of 1,800 in Simulation 2); the largest entry in each was adjusted by at most 2.2e-16 (4.2e-16 over all replications). Our own matching on the adjusted matrix gives identical estimates.

## Appendix C. Files of the cross-check

Moved here from BRIDGE.md on 2026-11-17; the per-replication CSVs were gzip-compressed the same day (`.csv.gz`; `read.csv()` reads them as they are). Scripts 09–13 and the notes they support are also described in [crosscheck/intervals/README.md](../crosscheck/intervals/README.md) and [crosscheck/rebuild2/README.md](../crosscheck/rebuild2/README.md).

```
crosscheck/.gitignore              ignores the regenerable *.rds written by 01/02
crosscheck/R/01_export_rep001.R    regenerates rep 1, checks it against the stored diagnostics, exports CSVs
crosscheck/R/02_pkg_fits.R         three multiMatch() fits, comparison table, unit-level diagnostics
crosscheck/R/03_check_pasted.R     tests a cross-check script suggested by another assistant
crosscheck/results/sim1_rep001.csv            T, Y, x1-x4, z1-z3 (17 significant digits)
crosscheck/results/sim1_rep001_ours.csv       our CM, GPS (and DIF) estimates for the 6 pairs
crosscheck/results/sim1_rep001_gps.csv        our fitted GPS matrix, 400 x 4
crosscheck/results/sim1_rep001_crosscheck.csv the comparison table above, full precision
crosscheck/results/pkg_session_info.txt       sessionInfo() of the package run (written by 02)
crosscheck/logs/                   console logs of 01-03
crosscheck/drafts/d7_d13_wording.md  draft wording for D7, D13 and related sentences (applied; kept for the record)
crosscheck/R/04_crosscheck_all.R   all replications of Simulation 1 or 2: regenerate, verify, our estimates, three package fits
crosscheck/R/05_summarize_crosscheck.R  verification, agreement, bias, MSE and coverage summaries
crosscheck/R/06_existing_mismatches.R  unit-level explanation of every `existing` difference
crosscheck/results/sim{1,2}_crosscheck_reps.csv.gz     every estimate, ours and the package's, with package variances
crosscheck/results/sim{1,2}_crosscheck_info.csv.gz     per-replication checks (diagnostics, ties, differing cells, ...)
crosscheck/results/sim{1,2}_crosscheck_summary.csv  bias, MSE and coverage by estimator, pair and target
crosscheck/results/sim{1,2}_crosscheck_multinom_minus_ours.csv  paired differences by pair
crosscheck/results/crosscheck_{agreement,verification,overview}.csv  the tables in docs/crosscheck.md
crosscheck/results/existing_mismatches.csv  the 30 near-ties
crosscheck/reps/                  per-replication checkpoints (*.rds, git-ignored, local only)
docs/crosscheck.md                 the full report
docs/crosscheck_math.md            why the three GPS versions compute the same estimator (math, from the code of both)
crosscheck/R/07_d8_analytic_intervals.R  D8: AI 2006 and AI 2016 variances for the rebuild's own GPS fit
crosscheck/R/08_side_by_side.R    side-by-side tables and figures, D8 summaries
crosscheck/results/sim{1,2}_d8_gps_variances.csv.gz  per-replication analytic GPS variances
crosscheck/side_by_side/          README.md (discussion: package rerun, CM vs GPS, D7, the other D's, D8), tables/, figures/
crosscheck/R/09_cm_bias_decomposition.R  exact split of each CM estimate (sample mix, matching gap, noise); regression-adjusted check; AI 2006 variance parts
crosscheck/R/10_cm_bootstrap_centring.R  CM-only rerun of the percentile bootstrap, to see where its intervals are centred
crosscheck/R/11_intervals_summary.R     tables and figure for crosscheck/intervals/
crosscheck/results/sim{1,2}_cm_{gap,terms,bc,boot}_reps.csv.gz  per-replication results of 09 and 10
crosscheck/intervals/             README.md (interval methods: sources, math, coverage; why CM coverage is low), tables/, figures/
crosscheck/R/12_cm_biasadjust.R   CM through Matching::Match, without and with BiasAdjust = TRUE, on the rebuild's samples
crosscheck/R/13_rebuild2_comparison.R  tables and figure for crosscheck/rebuild2/
crosscheck/results/sim{1,2}_cm_biasadjust_reps.csv.gz  per-replication results of 12
crosscheck/rebuild2/              README.md (rebuild-2 vs this rebuild vs the paper; bias-adjustment test), data/ (rebuild-2's summaries, copied), tables/, figures/
```

The script suggested by another assistant fails as written: a plain matrix `X` errors on R >= 4.3, `fit$estimate` partially matches `estimate_args`, `summary()` returns `NULL`, it does not flip signs, and it omits `existing`. `03_check_pasted.R` records this.

## References

- Abadie, A., & Imbens, G. W. (2006). Large sample properties of matching estimators for average treatment effects. *Econometrica, 74*(1), 235–267.
- Abadie, A., & Imbens, G. W. (2008). On the failure of the bootstrap for matching estimators. *Econometrica, 76*(6), 1537–1557.
- Abadie, A., & Imbens, G. W. (2016). Matching on the estimated propensity score. *Econometrica, 84*(2), 781–807.
- Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19: Multilevel covariates matching and propensity score matching. *Annals of Data Science, 9*(5), 967–982.
- Yang, S., & Barkley, B. G. (2019). *multilevelMatching* (R package version 1.0.0).
- Yang, S., Imbens, G., Cui, Z., Faries, D., & Kadziola, Z. (2016). Propensity score matching and subclassification in observational studies with multi-level treatments. *Biometrics, 72*, 1055–1065.
