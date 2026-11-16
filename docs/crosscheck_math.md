# Why the three GPS versions compute the same estimator

The cross-check ([docs/crosscheck.md](crosscheck.md)) compares three GPS matching estimates for every replication: the rebuild's own, the package's `existing` fit and the package's `multinom` fit. This note shows that all three compute one and the same estimator, a fixed function of two inputs. It then explains why `existing` agrees with the rebuild except at near-ties, and why `multinom` differs in many Simulation 2 replications. The same argument covers covariate matching (CM).

Throughout, "the package" is `multilevelMatching` 1.0.0 (Yang and Barkley, 2019), which calls `Matching::Match`; "the rebuild" is `R/estimators.R`.

## Summary

The estimator is a function $`\Phi(r, \varepsilon)`$ of a score $`r`$ (what units are matched on) and a tie tolerance $`\varepsilon`$ (how close two candidates must be to count as tied):

| Version | Score r | Tie tolerance ε | Estimate |
|---|---|---|---|
| Our GPS (rebuild) | our fitted probabilities p̂ | 0 (exact) | Φ(p̂, 0) |
| `existing` (package, given our scores) | the same p̂ | about 1e-10 | Φ(p̂, ε) |
| `multinom` (package's own fit) | its own fitted probabilities p̃ | about 1e-10 | Φ(p̃, ε) |

`existing` changes only the tie rule, so it can differ from ours only at near-ties. `multinom` also changes the score, so it differs wherever p̃ and p̂ rank candidates differently.

## 1. The estimator as a function of its inputs

The data are $`(Y_i, T_i, w_i)`$ for $`i = 1, \dots, n`$, with treatment levels $`T_i \in \{1, \dots, A\}`$ and $`w_i = (1, x_{1}, \dots, x_{4}, z_{1}, z_{2}, z_{3})_i`$. A score assigns each unit a number for each level, $`r = (r_1, \dots, r_A)`$; for GPS, $`r_s(i)`$ is unit i's fitted probability of level s.

**Step 1: nearest neighbours.** For a unit i outside level s ($`T_i \ne s`$), the distance to a level-s unit j is $`d(i,j) = \lvert r_s(i) - r_s(j) \rvert`$. The set of candidates within ε of the nearest is

```math
N_\varepsilon(i, s) = \Bigl\{\, j : T_j = s,\;\; d(i, j) \le \min_{k:\,T_k = s} d(i, k) + \varepsilon \,\Bigr\}.
```

If $`N_\varepsilon(i, s)`$ has more than one member, one is drawn uniformly at random; call the chosen unit $`j(i, s)`$. Matching is with replacement: a level-s unit can serve as the match for many units.

**Step 2: imputation.** Units in level s keep their observed outcome; all others take their match's:

```math
\hat Y_i(s) =
\begin{cases}
Y_i & \text{if } T_i = s, \\
Y_{j(i,s)} & \text{if } T_i \ne s.
\end{cases}
```

**Step 3: averaging.** The mean potential outcome and the pairwise effect are

```math
\hat\mu_s = \frac{1}{n} \sum_{i=1}^{n} \hat Y_i(s), \qquad \hat\tau(t, t') = \hat\mu_t - \hat\mu_{t'} .
```

This is the matching estimator of Yang et al. (2016), and the paper's Eqs. (3)–(4) read as in D7. Call its value $`\Phi(r, \varepsilon)`$.

## 2. The rebuild computes Φ(p̂, 0)

In `R/estimators.R`:

- **Step 1.** `impute_from_scores()` builds `abs(outer(P[, s], P[g, s], "-"))`, the distances $`\lvert \hat p_s(i) - \hat p_s(j) \rvert`$ from every unit to the units g in level s. `nn_pick()` takes the exact minimum (ε = 0) and breaks exact ties at random. No exact ties occurred in either simulation.
- **Step 2.** `impute_by_matching()` sets `Yimp[own, s] <- Y[own]` for units in level s and the matched outcome for the rest.
- **Step 3.** `tau_from_imputed()` takes `colMeans(Yimp)`, which gives $`\hat\mu_s`$, and then the differences $`\hat\mu_t - \hat\mu_{t'}`$.

## 3. The package computes Φ(r, ε) for whatever score it is given

In `multilevelMatching` (`matchAllTreatments()`, `matchImputePO()`, `wrangleImputations()`, `estimateTau()`):

**Step 1: its distance ranks candidates in the same order.** For level s the package takes the score column `X_kk <- X[, kk]` (our p̂ for `existing`, its own p̃ for `multinom`) and calls

```
Matching::Match(Y, Tr = (W != s), X = X_kk, Weight = 2, M = 1,
                distance.tolerance = 0, ties = FALSE)
```

With the defaults `estimand = "ATT"` and `replace = TRUE`, every unit outside s (the "treated", `Tr = 1`) is matched, with replacement, to one unit in s (the "controls"). `Match` first standardizes the score with its sample mean m and sample standard deviation $`\sigma_s`$, $`\tilde x = (r_s - m)/\sigma_s`$. It then weights by $`W = (\tilde x^\top \tilde x / N)^{-1}`$ (the Mahalanobis option, `Weight = 2`). Since $`\tilde x^\top \tilde x = N - 1`$, this gives $`W = N/(N-1)`$ and the squared distance

```math
d^2_{\text{pkg}}(i, j) = W \,\bigl(\tilde x_i - \tilde x_j\bigr)^2 = \underbrace{\frac{N}{(N-1)\,\sigma_s^2}}_{c\,>\,0}\,\bigl(r_s(i) - r_s(j)\bigr)^2 = c\, d(i,j)^2 .
```

Squaring and multiplying by a positive constant do not change which candidate is nearest, so the package's nearest neighbours are ours. The one difference is the tie tolerance. Even with `distance.tolerance = 0`, the compiled matching loop treats candidates whose squared standardized distances lie within about $`10^{-10}`$ of the minimum as tied. With `ties = FALSE` it then picks one of them at random. This was confirmed directly on replication 351 of Simulation 1, where two candidates 8.3e-11 apart in squared standardized distance were treated as tied. Constructed examples gave the same threshold at every distance scale tried: between about 1e-10 and 5e-10 in squared standardized distance. That is the package's ε.

**Step 2: same imputation.** `wrangleImputations()` sets `Yiw[same] <- Y[same]` and `Yiw[diff] <-` the matched outcome, which is $`\hat Y_i(s)`$.

**Step 3: same averages, opposite sign.** `mean_Yiw[s] <- mean(Yiw_s)` is $`\hat\mu_s`$. `estimateTau()` reports, for $`t < t'`$,

```math
\text{Estimate}(t, t') = \hat\mu_{t'} - \hat\mu_t = -\,\hat\tau(t, t'),
```

labelled `EY(t')-EY(t)`, so the cross-check flips its sign. The package also sorts units by treatment before matching. That changes their order but not the matches or the averages; the simulated samples are already sorted.

Steps 1–3 compute $`\Phi(r, \varepsilon)`$ whatever score is passed in. With our p̂ (`existing`) the result is $`\Phi(\hat p, \varepsilon)`$; with the package's own p̃ (`multinom`) it is $`\Phi(\tilde p, \varepsilon)`$.

## 4. Covariate matching: the same distance up to a constant

For `match_on = "covariates"` the score is replaced by the seven covariates $`v = (x_1, \dots, x_4, z_1, z_2, z_3)`$, and the rebuild uses the Mahalanobis distance with the pooled sample covariance S (D6):

```math
d(i, j)^2 = (v_i - v_j)^\top S^{-1} (v_i - v_j).
```

`Match` standardizes each column, $`\tilde v = D^{-1}(v - m)`$, where D is the diagonal matrix of sample standard deviations. With R the sample correlation matrix, $`\tilde X^\top \tilde X / N = \tfrac{N-1}{N} R`$, so

```math
d^2_{\text{pkg}}(i, j) = (\tilde v_i - \tilde v_j)^\top \Bigl(\tfrac{N-1}{N} R\Bigr)^{-1} (\tilde v_i - \tilde v_j)
= \tfrac{N}{N-1}\,(v_i - v_j)^\top D^{-1} R^{-1} D^{-1} (v_i - v_j)
= \tfrac{N}{N-1}\, d(i, j)^2 ,
```

because $`S = D R D`$. (`Match` would add a small ridge to the matrix if it were nearly singular; with these covariates it is not.) Again the nearest neighbours are the same, and the rest of the estimator follows Steps 2–3.

## 5. Where the two scores come from

Both scores come from the same multinomial logit, with level 1 as the baseline ($`\beta_1 = 0`$):

```math
p_s(w; \beta) = \frac{\exp(w^\top \beta_s)}{\sum_{k=1}^{A} \exp(w^\top \beta_k)} .
```

Both are fitted by maximizing the same log-likelihood

```math
\ell(\beta) = \sum_{i=1}^{n} \log p_{T_i}(w_i; \beta)
```

with the same optimizer (`nnet::multinom`, BFGS). The optimizer stops when an iteration improves the fit by less than a relative `reltol`, or after `maxit` iterations. Only these two settings differ:

| | Coefficients | reltol | maxit | Score |
|---|---|---|---|---|
| Rebuild (D5) | β̂ | 1e-14 | 10,000 | p̂ = p(w; β̂) |
| Package `multinom` | β̃ | 1e-8 (nnet default) | 100 (nnet default) | p̃ = p(w; β̃) |

At the rebuild's settings the fitted probabilities do not depend on the optimizer's starting point to about 5e-8 (D5). At nnet's defaults the optimizer stops earlier, so β̃ is a coarser approximation to the same maximum, and in 2 of the 500 Simulation 2 replications it stops at maxit = 100 before converging. Both p̂ and p̃ approximate the same fitted probabilities, but p̃ less precisely.

For `existing`, the package also requires each row of the score matrix to sum to exactly 1 in floating point, so rows off by rounding were adjusted by at most 4.2e-16. That moves any distance by at most 8.4e-16, far below the smallest gap between competing candidates in either simulation (about 2e-10), so it cannot change a match.

## 6. When the outputs differ

The estimate depends on the data only through the matches $`j(i, s)`$: if two versions choose the same match for every unit and level, every imputed outcome $`\hat Y_i(s)`$ is the same, and so is every $`\hat\tau(t, t')`$. Differences can arise only where some match changes.

**`existing` against the rebuild: same score, different tie rule.** A match can differ only if $`N_\varepsilon(i, s)`$ has two or more members and the random draw is not the exact nearest neighbour. This happened for 30 imputed outcomes in 25 of the 1,000 replications. In every case the rebuild's match was the exact nearest, and the package's choice was farther by at most 9.3e-11 in squared standardized distance ([existing_mismatches.csv](../crosscheck/results/existing_mismatches.csv)).

**`multinom` against the rebuild: different score.** Let $`\delta = \max_{i,s} \lvert \tilde p_s(i) - \hat p_s(i) \rvert`$ be the largest probability difference in a replication. Every distance then moves by at most $`2\delta`$:

```math
\bigl|\, \lvert \tilde p_s(i) - \tilde p_s(j) \rvert - \lvert \hat p_s(i) - \hat p_s(j) \rvert \,\bigr| \;\le\; \lvert \tilde p_s(i) - \hat p_s(i) \rvert + \lvert \tilde p_s(j) - \hat p_s(j) \rvert \;\le\; 2\delta .
```

Let $`j^*`$ be unit i's nearest level-s candidate under p̂, and $`g = \min_{j \ne j^*,\, T_j = s} d(i, j) - d(i, j^*)`$ its margin over the runner-up. Under p̃, every other candidate is still farther away by at least $`g - 4\delta`$. So

```math
g > 4\delta \;\;\Longrightarrow\;\; \text{unit } i \text{ keeps the same match } j^* \text{ for level } s
```

(ignoring the package's ε, which is negligible here). Matches can switch only for units whose margin is at most $`4\delta`$:

| | Simulation 1 | Simulation 2 |
|---|---|---|
| δ per replication, median | 1.2e-6 | 2.7e-4 |
| δ per replication, largest | 1.7e-4 | 2.2e-3 |
| Imputed outcomes that switch, average per replication | 0.4 of 1,600 (0.03%) | 1,003 of 16,200 (6.2%) |
| Replications with no switch | 405 of 500 | 0 of 500 |

GPS matching distances in Simulation 2 average 6.5e-4, the same order as δ, so many margins fall below 4δ. Single estimates move by up to 2.7, but over replications the switches largely average out. The mean package-minus-rebuild difference is at most 0.027 in any pair (|z| > 2 in 9 of 36 pairs, none above 3; the pairs share levels, so these are not independent), and mean absolute bias against τ\_pop is 0.0995 for `multinom` against 0.0991 for the rebuild (Simulation 2).

**CM.** The covariate distances agree up to the constant N/(N − 1), and no two candidates came within the package's tie tolerance: the smallest margin between a unit's nearest and second-nearest CM candidate was 6.5e-7 in Simulation 1 and 9.4e-9 in Simulation 2 (about 3e-8 in squared distance). So `covariates` reproduces the rebuild's CM in all 1,000 replications. The only differences are floating-point rounding of the averages (at most 1e-13), because the package uses `mean()` and the rebuild `colMeans()`.

## Sources

- Rebuild: `R/estimators.R` (`impute_from_scores`, `nn_pick`, `impute_by_matching`, `tau_from_imputed`, `fit_gps`, `mahalanobis_d2`).
- Package: `multilevelMatching` 1.0.0 (`R/matchAllTreatments.R`, `R/matchImputePO.R`, `R/estimateTau.R`) and `Matching` 4.10-14 (`Match`, `RmatchLoop`, the compiled `MatchLoopC`).
- Numbers: `crosscheck/results/crosscheck_agreement.csv`, `sim{1,2}_crosscheck_info.csv`, `sim{1,2}_crosscheck_multinom_minus_ours.csv` and `existing_mismatches.csv`; see [docs/crosscheck.md](crosscheck.md).
- Yang, S., Imbens, G., Cui, Z., Faries, D., & Kadziola, Z. (2016). Propensity score matching and subclassification in observational studies with multi-level treatments. *Biometrics, 72*, 1055–1065. https://doi.org/10.1111/biom.12505
- Yang, S., & Barkley, B. G. (2019). *multilevelMatching* (R package version 1.0.0). https://CRAN.R-project.org/package=multilevelMatching
