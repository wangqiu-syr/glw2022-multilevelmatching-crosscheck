# Simulation spec: multilevel covariate matching vs. GPS matching

Source: Guo, Liu & Wang (2022), *Annals of Data Science* 9(5):967-982, https://doi.org/10.1007/s40745-022-00392-x. This file restates the set-up in our own words and fixes the points the paper leaves open (marked **D1–D13**, listed at the end). D1–D8 were fixed before the full run; D9–D13 were added during the work. Any reader who seeks to validate the proposed methods, or to reproduce the reported results, from the published description alone and without access to the original code faces the same thirteen decisions. For each, this rebuild adopts a default that either follows directly from the paper's description of its simulations or, where the text is ambiguous, is a reasoned reading of it; each choice and its rationale are recorded below. Under these defaults, the rebuilt simulations support the paper's central claims: covariate matching (CM) and generalized propensity score (GPS) matching both remove most of the bias of the unadjusted comparison, and GPS matching shows the smallest bias of the three estimators. The rebuild therefore offers no evidence against the use of generalized propensity score matching to estimate the effects of multilevel treatments in the setting the paper considers.

## Goal

Replicate the two simulation studies (paper Section 5; Tables 1 and 3) in R and compare three estimators of the pairwise average treatment effect across treatment combinations:

- **DIF**: raw difference in group means
- **CM**: covariate matching
- **GPS**: generalized propensity score matching

## Notation

- X = (X1, X2, X3, X4): individual-level covariates. Z = (Z1, Z2, Z3): school-level covariates.
- T in {1, ..., A}: treatment combination. A = 4 in Simulation 1, A = 9 in Simulation 2.
- Estimand: tau(t, t') = E[Y(t) - Y(t')] for every pair t < t'. This gives 6 effects in Simulation 1 and 36 in Simulation 2.

## Covariates (both simulations)

X ~ MVN(0, Sigma), with

```
Sigma = [  2.0   1.0  -1.0   0.5
           1.0   1.0  -0.5   0.2
          -1.0  -0.5   1.0  -0.3
           0.5   0.2  -0.3   2.0 ]
```

(positive definite; smallest eigenvalue about 0.31). Independent of X and of each other:

- Z1 ~ Uniform(-3, 3)  **(D2)**
- Z2 ~ chi-square(1)
- Z3 ~ Bernoulli(0.5)

## Treatment assignment

Multinomial logit, with eta_1 = 0 as the baseline:

p_t(x, z) = exp(eta_t) / sum over s of exp(eta_s)  **(D1)**

**Simulation 1 (A = 4):**

| t | eta_t / 0.1 |
|---|---|
| 1 | 0 (eta_1 = 0) |
| 2 | x1^2 + x2 + x3 + x4 + z1 + z2 + z3 |
| 3 | x1 + x2^2 + x3^2 + x4 + z1 + z2 + z3 |
| 4 | x1 + x2^2 + x3^2 + x4^2 + z1 + z2 + z3 |

**Simulation 2 (A = 9):** eta_1 = 0, and eta_t = 0.1 times:

| t | expression |
|---|---|
| 2 | x1^2 + x2 + x3 + x4 + z1 + z2 + z3 |
| 3 | x1 + x2^2 + x3^2 + x4 + z1 + z2 + z3 |
| 4 | x1 + x2^2 + x3^2 + x4^2 + z1 + z2 + z3 |
| 5 | x1 + x2^2 + x3^2 + x4^2 + z1^2 + z2 + z3 |
| 6 | x1 + x2^2 + x3^2 + x4^2 + z1^2 + z2^2 + z3 |
| 7 | x1 + x2^2 + x3^2 + x4^2 + z1^2 + z2^2 + z3^2 |
| 8 | x1^2 + x2^2 + x3^2 + x4^2 + z1^2 + z2^2 + z3^2 |
| 9 | x1 + x2 + x3 + x4 + z1 + z2 + z3 |

Only the index t matters. Ignore the cell labels (Aa, Ba, ...) in the paper's Table 2, which contain numbering typos.

**Note (D10).** eta_6 and eta_7 differ only in z3 versus z3^2. z3 is binary, so z3^2 = z3, which makes eta_6 = eta_7 and p_6 = p_7 for every unit. We keep the paper's specification as written. The two groups still have different outcome coefficients, so tau(6, 7) is well defined. Expect the fitted p_hat_6 to be close to p_hat_7.

**Note (D11).** Groups t1-t5 and t9 have thin lower tails in p_t. In a 1e5-unit check, 3.4% of units have some p_t below 0.01, and the 1st percentile of p_1 is about 0.0016. This means overlap is poor for those levels. No trimming is applied. Matching diagnostics (mean matching distance per target group, and the share of units with p_hat below 0.01) are recorded for every replication.

## Potential outcomes

epsilon ~ N(0, 1), drawn independently for each unit. Observed outcome Y = Y(T).

**Simulation 1:** Y(t) = beta_t' (x1, x2, x3, x4, z1, z2, z3) + epsilon, with

| t | beta_t |
|---|---|
| 1 | (1, 1, 1, 1, 1, 1, 1) |
| 2 | (2, 3, 1, 2, 2, 2, 2) |
| 3 | (3, 1, 2, 1, 3, 1, 3) |
| 4 | (4, 1, 3, 4, 4, 1, 4) |

**Simulation 2:** Y(t) = c_t' (x1, x2, x3, x4, z1, z2 - 1, z3 - 0.5) + epsilon, with

| t | c_t |
|---|---|
| 1 | (1, 1, 1, 1, 1, 1, 1) |
| 2 | (2, 3, 1, 2, 2, 2, 2) |
| 3 | (3, 1, 2, 1, 3, 1, 3) |
| 4 | (4, 1, 3, 4, 4, 1, 4) |
| 5 | (5, 1, 1, 1, 1, 1, 1) |
| 6 | (6, 3, 1, 2, 2, 2, 9) |
| 7 | (7, 1, 8, 1, 3, 1, 8) |
| 8 | (8, 1, 5, 2.5, 4, 1, 7) |
| 9 | (9, 1, 7, 1.5, 4, 1, 7) |

The paper prints no epsilon in Y(1) for Simulation 2. Treat that as a typo and add it.

## Sample design  (D3)

| | Simulation 1 | Simulation 2 |
|---|---|---|
| Units per treatment combination | 100 | 200 |
| Total n | 400 | 1800 |
| Replications | 500 | 500 |

The paper fixes the group sizes but also assigns treatment by a multinomial model. Default: draw a large pool of units, assign T from p_t, then keep the first n_g units in each group. The pool is large enough that every group fills; R/dgp.R stops with an error if any group falls short.

## True effects  (D4, D9)

D-numbers (D1–D13) label the decisions this rebuild had to make where the paper is silent or ambiguous; the full list is in the [decisions table](#decisions-d1d13).

Write W for the outcome regressors: W = (x1, x2, x3, x4, z1, z2, z3) in Simulation 1 and W = (x1, x2, x3, x4, z1, z2 - 1, z3 - 0.5) in Simulation 2. Because outcomes are linear in W, every target has the form tau(t, t') = (beta_t - beta_t')' m for some covariate mean m. Two targets are computed and every result is reported against both.

D4 (target: τ_pop). Decision 4 sets which true effect the bias is measured against. **Primary: population ATE, τ_pop.** m = mu = E[W] over the superpopulation. With Z1 mean 0:

- **Simulation 1:** mu = (0, 0, 0, 0, 0, 1, 0.5), so τ_pop(1, 2) = -1.5, and so on.
- **Simulation 2:** mu = 0, so all 36 τ_pop are 0.

Verify these closed forms against a Monte Carlo average of Y(t) - Y(t') on a large population (about one million (1,000,000) units, large enough that Monte Carlo error is negligible next to the biases reported). `R/check_truth.R` uses exactly one million (1,000,000) units for this τ_pop check (`N_POP <- 1e6`).

D9 (target: τ_bal). Decision 9 records what the matching estimators actually estimate when every group has a fixed size (see D3). **Secondary: balanced-design ATE, τ_bal.** m = mu_bal, the covariate mean of the sampling design in D3 (equal group sizes):

mu_bal = (1/A) * sum over t of E[W | T = t],  with  E[W | T = t] = E[W p_t(X, Z)] / E[p_t(X, Z)].

Compute mu_bal by Monte Carlo on a pool of two million (2,000,000) units in `R/check_truth.R` (`N_BAL <- 2e6`), weighting by p_t rather than by the sampled T, which lowers the variance. Check it against a direct average of Y_i(t) - Y_i(t') over many simulated balanced samples. τ_bal is generally not 0 in Simulation 2. See D9 for why this target exists.

**Why Simulation 2 centres z2 and z3.** Simulation 2's outcome uses z2 − 1 and z3 − 0.5 instead of z2 and z3. Because Z2 averages 1 and Z3 averages 0.5 (and X and Z1 average 0), every term in the outcome averages zero over the population, so all 36 population effects are exactly 0. Effects still differ from unit to unit, since each treatment has its own coefficients, but they cancel on average. Treatment assignment depends on the same covariates, so the treatment groups differ in X and Z, and any raw difference between groups is bias rather than effect. In a real study, this is the case where two learning modes have the same average effect across a district, but schools with more digital devices (a school-level covariate, one of Z1–Z3) choose one mode more often: a raw comparison shows a difference that matching on X and Z should remove. It also explains why τ_bal is not 0: the balanced sample over-represents some groups' covariate mix, over which the varying effects no longer cancel.

## Estimators

For each replication and each pair (t, t'):

**DIF.** mean(Y | T = t) minus mean(Y | T = t').

**CM.** For every unit i in the sample and every level s, impute Y_i(s) = Y_i if T_i = s. Otherwise use the Y of the nearest unit j with T_j = s, where distance is Mahalanobis on all 7 covariates (X and Z together, V = pooled sample covariance)  **(D6)**. Use 1 match, with replacement, and break ties at random. Then tau_hat_CM(t, t') = mean over all i of [Y_i(t) - Y_i(t')].

**GPS.** Fit a multinomial logit of T on x1..x4, z1, z2, z3 (linear terms only, e.g. `nnet::multinom`)  **(D5)**. This is deliberately misspecified, because the true eta contains squares. For every unit i and level s, impute Y_i(s) = Y_i if T_i = s. Otherwise use the Y of the unit j with T_j = s that minimizes |p_hat_s(x_j, z_j) - p_hat_s(x_i, z_i)|. Then tau_hat_GPS(t, t') = mean over all i of [Y_i(t) - Y_i(t')]. Use p_hat_t' for the t' term. The paper's Eq. (4) prints p_t in both terms, which looks like a typo  **(D7)**. In the implementation, column s of the n x A imputed-outcome matrix is always matched on p_hat_s, so the t' term automatically uses p_hat_t'. The R package of Yang et al. (2016), `multilevelMatching`, matches the same way; on a test replication it reproduces our GPS estimates exactly (D7).

## Inference and coverage  (D8)

The paper does not say how its confidence intervals were built. Default: nonparametric bootstrap, percentile 95% CI, B = 200 for the final run (B = 20 for tests).

- **Stratified resampling.** In each bootstrap draw, resample n_g units with replacement within each treatment group, so the group sizes stay fixed as in the design (D3).
- **Full re-estimation.** Recompute every estimator from scratch on each bootstrap sample. That means refitting the GPS multinomial logit, and recomputing the pooled covariance V for CM.
- **One CI, two targets.** Each replication gives one CI per estimator and pair. Coverage is the share of replications whose CI contains τ_pop, and separately the share whose CI contains τ_bal.
- **Timing.** Time the 20-replication test first. If B = 200 is too slow, lower B and record the value used.

Bootstrapping matching estimators is not formally valid (Abadie and Imbens, 2008), so coverage may differ from the paper. The main run computes no Abadie–Imbens standard errors and does not use the `Matching` package.

**Sensitivity analysis (added after the main run).** The cross-check with `multilevelMatching` adds analytic Abadie–Imbens intervals as a sensitivity analysis of D8. It keeps every point estimate of the main run, and therefore its bias and MSE, and replaces only the intervals, with estimate ± 1.96 SE: the AI 2006 variance (Abadie and Imbens, 2006) for CM and GPS, and for GPS also the AI 2016 variance (Abadie and Imbens, 2016), which allows for the estimated GPS. The bootstrap remains the D8 default. Results: `crosscheck/side_by_side/README.md`, section 5.2; code: `crosscheck/R/07_d8_analytic_intervals.R`.

## Outputs

Per estimator, pair, and target (τ_pop and τ_bal), over the replications:

- signed bias = mean(tau_hat) - tau
- absolute bias = |mean(tau_hat) - tau| (Table 1 reports this)
- MSE = mean((tau_hat - tau)^2)
- coverage of the 95% CI (%)

The summary tables use long format, with a `target` column taking the values `pop` or `bal`. Write the true effects to `results/true_effects_sim1.csv` and `results/true_effects_sim2.csv`. Write tidy summaries to `results/sim1_summary.csv` and `results/sim2_summary.csv`. Write per-replication estimates to `results/*.rds` (git-ignored). Reproduce the bias plot of the 36 effects (paper Fig. 1) for Simulation 2.

## Comparison with the paper: sign of the bias

In this project, bias = mean(tau_hat) - tau, with tau(t, t') = E[Y(t) - Y(t')]. The paper's Simulation 2 bias values (`results/paper_table3.csv`, transcribed as published) correlate negatively with ours. `R/plots.R` therefore draws the overlay twice, keeping both versions:

- `results/fig1_bias_sim2_pop_vs_paper.png`: the paper's bias as published.
- `results/fig1_bias_sim2_pop_vs_paper_flipped.png`: the paper's bias multiplied by -1.

Observed correlations across the 36 pairs, with our bias against τ_pop (500 replications, B = 200). They are saved in `results/sim2_bias_correlations.csv`, together with the τ_bal versions.

| Estimator | r, paper sign as published | r, paper sign flipped | Same sign (published / flipped) |
|---|---|---|---|
| DIF | -0.954 | 0.954 | 5 / 31 of 36 |
| CM | -0.764 | 0.764 | 7 / 29 of 36 |
| GPS | -0.874 | 0.874 | 7 / 29 of 36 |

The flipped version correlates at 0.95 (DIF), 0.76 (CM) and 0.87 (GPS) with our bias. Flipping the sign only negates r, so the two columns differ only in sign. Against τ_bal, the flipped correlations are 0.94, 0.59 and 0.47.

This is an observation about the published numbers. It does not establish which convention the paper used, for example whether it reports tau - E(tau_hat) or defines the effect as Y(t') - Y(t). The paper's Simulation 1 Table 1 reports absolute bias, so it cannot show the sign. Comparisons of absolute bias, MSE and coverage do not depend on the sign.

## Implementation notes

- Folder plan: `R/dgp.R`, `R/estimators.R`, `R/run_sim1.R`, `R/run_sim2.R`, `R/summarize.R`.
- Stage 1 is `R/dgp.R` plus `R/check_truth.R`, which runs the sanity checks on the data-generating process and both true-effect targets.
- Per-replication seeds (pre-generated L'Ecuyer-CMRG streams), parallelized with `parallel`. Test with 20 replications before 500.
- Package versions: the scripts use only base R, `parallel` and `nnet`, plus `MASS` for data generation. The R version and session details of the full run are in `results/run_session_info.txt`, and the `nnet` and `MASS` versions are logged in `results/run_log.txt`. No `renv` lockfile exists. Adding one with `renv` is an option for later.
- Sanity checks: p_t sums to 1; every group reaches its target size; closed-form τ_pop matches Monte Carlo; the Monte Carlo τ_bal matches a direct average over balanced samples; GPS overlap looks sensible.
- Stop after each stage for review. `.Rhistory` is never committed.

## Out of scope

Paper Section 6 (the learning-mode application). Its data-generating process is not described, so Table 5 cannot be reproduced. For reference, the cell sizes in the paper's Table 4 are 8, 56, 118, 93, 75 and 70.

## Decisions D1–D13

| ID | Issue | Default |
|---|---|---|
| D1 | Paper's denominator is 1 + sum(exp(eta)), which makes probabilities sum below 1 | Standard multinomial logit, no extra 1 |
| D2 | "Z1 ~ N[-3, 3]" is ambiguous | Uniform(-3, 3). Its mean of 0 is consistent with the centred z2, z3 terms, which make all true effects 0 in Simulation 2 |
| D3 | Fixed group size versus random assignment | Large pool, assign by p_t, keep n_g per group |
| D4 | True tau not stated | Two targets, both reported: population ATE τ_pop (primary; closed form checked by Monte Carlo) and balanced-design ATE τ_bal (secondary; Monte Carlo) |
| D5 | GPS model not stated | Linear multinomial logit (misspecified on purpose), fitted with `nnet::multinom` at reltol = 1e-14 and maxit = 10000. The defaults (reltol = 1e-8, maxit = 100) left the fitted probabilities too uncertain for the very small matching distances. Cross-check (`docs/crosscheck.md`): `multilevelMatching` 1.0.0 fits the same multinomial logit at those defaults. Its fitted probabilities differ from ours by up to 1.7e-4 (Simulation 1) and 2.2e-3 (Simulation 2), which changes about 6% of the matched outcomes in every Simulation 2 replication and single estimates by up to 2.7, while mean absolute bias hardly changes (0.0995 against 0.0991, τ_pop); 2 of its 500 Simulation 2 fits stop at maxit = 100. Even at reltol = 1e-14, fitted probabilities differ at about the 1e-9 level across setups (the cross-check's R 4.3.3, nnet 7.3-19, Linux against the main run's R 4.5.1, nnet 7.3.20, Windows). Near-ties are common in Simulation 2, so a few GPS matches fall differently, and the Simulation 2 GPS summary means reproduce across setups to within 5.6e-5 rather than exactly; DIF, CM and the Simulation 1 GPS reproduce exactly |
| D6 | Paper defines V from X only but matches on X and Z | Mahalanobis on (X, Z) jointly |
| D7 | Eq. (4) uses p_t in both terms | Use p_t' for the t' term: column s of the imputed-outcome matrix is always matched on p_hat_s. This is how the R package of Yang et al. (2016), to whose methods the paper's reduce without Z, implements GPS matching: `multilevelMatching::multiMatch()` (version 1.0.0, CRAN 2019) imputes Y(s) by matching on column s of the GPS matrix. Cross-check on replication 1 of Simulation 1 (master seed 20220001, same sample as the main run): given our fitted GPS matrix (`match_on = "existing"`), `multiMatch()` reproduces our six GPS estimates and all 1,600 imputed potential outcomes exactly; with its own multinomial logit (`"multinom"`) and with covariate matching (`"covariates"`, Mahalanobis, pooled covariance) it also reproduces our GPS and CM estimates exactly. This shows that D7 is the estimator as implemented in the reference software for the method; it does not show what the paper's own code did. Run notes: R 4.3.3, Matching 4.10-14, nnet 7.3-19, no trimming, M_matches = 1, `set.seed(123)` (no ties occurred, so the seed did not affect the estimates). Version 1.0.0 stops on R >= 4.3 at its `class(X)` check, so `X` was passed with `oldClass(X) <- "matrix"` (values unchanged); `"existing"` requires every row to sum to exactly 1, so 75 of 400 rows were adjusted by at most 2.2e-16, which leaves our own estimates unchanged. Code and results: `crosscheck/R/01_export_rep001.R`, `crosscheck/R/02_pkg_fits.R`, `crosscheck/results/sim1_rep001_crosscheck.csv` |
| D8 | CI method not stated | Percentile bootstrap, B = 200, stratified by treatment group, full re-estimation (GPS refit, V recomputed) in each draw. Sensitivity analysis, added after the main run (`crosscheck/side_by_side/README.md`, section 5.2): the same estimates with analytic Abadie–Imbens intervals, estimate ± 1.96 SE. Mean coverage of τ_pop, Simulation 1 / 2: CM 94.8% / 90.1% with AI 2006, against 91.4% / 78.4% with the bootstrap (paper 95.0% / 97.1%); GPS 95.9% / 94.6% with AI 2016, against 100% / 100% with the bootstrap (paper 99.7% / 98.9%). For CM the paper's coverage is closer to AI 2006 than to the bootstrap (mean distance per pair 1.4 against 3.6 points in Simulation 1, 7.4 against 18.6 in Simulation 2); for GPS the bootstrap and AI 2006 are equally close. This is consistent with, but not proof of, the paper having used analytic variances. The bootstrap stays the default; bias and MSE are unchanged |
| D9 | The D3 design does not sample from the population: keeping n_g units per group gives an equal-weight mixture of f(w \| T = t), not f(w). CM and GPS average over the sample's units, so they estimate the ATE for that mixture rather than τ_pop | Report bias, MSE and coverage against both τ_pop and τ_bal |
| D10 | eta_6 = eta_7 in Simulation 2 (z3 is binary, so z3^2 = z3), so p_6 = p_7 | Keep the paper's specification; expect p_hat_6 close to p_hat_7 |
| D11 | Thin p_t tails in t1-t5 and t9 (poor overlap) | No trimming; record matching diagnostics per replication |
| D12 | How to record the software environment | Scripts use only base R, `parallel` and `nnet` (plus `MASS` for data generation). `R/run_all.R` writes `results/run_session_info.txt`. In the completed full run, that file lists only R and the base packages, and the `nnet` and `MASS` versions (7.3.20, 7.3.65) are in `results/run_log.txt`; later runs also list `nnet` and `MASS` in the file. No `renv` lockfile; adding one is an option for later |
| D13 | Sensitivity check of the D7 reading (Simulation 1 only) | Variant `GPS_lit` follows Eq. (4) as printed: for a unit in neither group t nor t', the t' term matches within group t' on \|p_hat_t'(j) - p_hat_t(i)\|. Units in group t or t' keep their observed Y on their own level, and their other term is imputed as in the main GPS. It is a variant that tests the D7 reading, not a replacement: the main GPS estimator and all main results are unchanged. Same seeds, samples and bootstrap draws as the main run (paired; checked: 500/500 replications reproduce the main results exactly). Code: `R/gps_literal.R`, `R/run_gps_literal.R`; results: `results/gps_literal_sim1.csv`. Second variant `GPS_lit2`, the fully literal reading: units in group t also get their t' term by matching within group t' on \|p_hat_t'(j) - p_hat_t(i)\|; every unit keeps its observed Y for its own group's term, and nothing else changes. It is computed in the same paired run as GPS and GPS_lit (same seeds and datasets), also as a sensitivity check, not a replacement. Code: `R/gps_literal2.R`, `R/run_gps_literal2.R`; results: `results/gps_literal2_sim1.csv`. Pairing checked: 500/500 replications reproduce the main DIF/CM/GPS results and the GPS_lit results exactly. Observation (absolute bias vs τ_pop, mean over the 6 pairs): paper GPS 0.293, main GPS 0.026, GPS_lit 0.044, GPS_lit2 0.159. GPS_lit2 is nearer the paper's value in 5 of 6 pairs, its bias is beyond Monte Carlo noise in 3 of 6 pairs (|z| > 3; 5 of 6 at 2 SE), and its MSE is 50-134% higher than the main GPS. This shows that the reading of Eq. (4) changes GPS bias in the direction of the paper's values; it does not establish which reading the paper used. Neither variant is what the package of Yang et al. (2016), `multilevelMatching`, computes: that package matches as in D7 and reproduces the main GPS exactly (see D7). If the paper's GPS results came from that package, the literal reading cannot explain the gap; it can do so only if the paper's own code followed Eq. (4) as printed. The paper names no software ("Code Availability N/A"). |

Exact agreement with the published tables should not be expected, because of these open points. The qualitative pattern should hold: DIF is clearly biased, while CM and GPS reduce the bias substantially.
