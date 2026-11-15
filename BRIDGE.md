# Bridge: from glw2022-matching-replication to glw2022-multilevelmatching-crosscheck

This file carries the work and decisions from the last session in the original repo into this one, so the next session (or reader) can pick up without the chat history. Written 2026-11-15.

## What this repo is

This repo starts as an import of [glw2022-matching-replication](https://github.com/wangqiu-syr/glw2022-matching-replication), as of commit `aea8268` ("Add exploratory Table 3 checks; relabel figures; clarify Simulation 1 scope"). The import is a single initial commit with the same 48 files and identical contents, without the earlier history, which stays in the original repo. The original repo is left unchanged.

This repo adds a cross-check: the rebuild's covariate matching (CM) and generalized propensity score (GPS) matching estimators are compared with the R package `multilevelMatching` of Yang et al. (2016). Guo, Liu & Wang (2022) state that without the school-level covariates Z their methods reduce to Yang et al. (2016), their reference [44]. In Yang et al., "multilevel" means a treatment with several levels, not nested data.

## Status

| Item | Status |
|---|---|
| Cross-check, replication 1 of Simulation 1 | Done: exact agreement (below) |
| Cross-check, all 500 replications of Simulation 1 | **On hold.** Not run; wait for the go-ahead |
| D7 and D13 wording in `docs/spec.md` and related docs | Applied (branch `claude/apply-d7-d13-wording`); drafts kept in `crosscheck/drafts/d7_d13_wording.md` |
| Repo-specific text inherited from the import (README, CITATION.cff) | Still describes the original repo; not yet updated |

## Result: replication 1 of Simulation 1

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
- **`multinom` agreement is not guaranteed in other replications.** The package fits `nnet::multinom` at its defaults (maxit 100, reltol 1e-8); ours uses reltol 1e-14 (D5). Fitted probabilities differ by up to 8.2e-7, while the smallest best-versus-second-best GPS match gap in this sample is 1.4e-7. No match changed here, but one could elsewhere, more likely in Simulation 2.

## multilevelMatching facts (version 1.0.0, CRAN 2019-05-08)

From `?multiMatch` and the package source (`R/prepareData.R`, `R/matchAllTreatments.R`, `R/estimateTau.R`):

- **Signature:** `multiMatch(Y, W, X, trimming = NULL, match_on, model_options = list(reference_level = sort(W)[1]), M_matches = 1, J_var_matches = 1)`.
- **Trimming:** default `NULL` = no trimming for `covariates` and `existing`; an error for `multinom` (must be set). `TRUE` is supported only for `multinom` (Crump et al. 2009 rule).
- **Matching:** `Matching::Match(..., Weight = 2, distance.tolerance = 0, ties = FALSE, M = M_matches)`, with replacement. `Weight = 2` is Mahalanobis on the full pooled sample, as in D6; there is no option to change the metric. Ties are broken at random, which is why a seed is set.
- **Output fields:** `results` (`Param`, `Trt1`, `Trt2`, `Estimate`, `Variance`, `VarianceAI2016`), `mu` (`Param`, `Trt`, `Estimate`), `analysis_idx`, `impute_mat`, `impute_mat_sorted`, `propensity_scores`, `model`, `estimate_args`, `impute_match_data`, plus `AI2016_args` for `multinom`. The help page's Value section lists older names (`tauestimate`, `varestimate`) that `multiMatch()` does not return.
- **Sign:** `Param` is `EY(Trt2)-EY(Trt1)` with Trt1 < Trt2, and `Estimate` = mu(Trt2) - mu(Trt1). Our tau(t, t') = -Estimate.
- **Treatment coding:** 0-based coding gives the same estimates with shifted labels. For `existing`, column k of the GPS matrix must belong to the k-th treatment level in sorted order.
- **GitHub master** (`shuyang1987/multilevelMatching`) has the same code as the CRAN release for every function `multiMatch()` uses.

## Workarounds needed (package not modified)

- **R >= 4.3:** `prepareData()` runs `(length(class(X)) != 1 && class(X) != "matrix") || ...`. A plain matrix has class `c("matrix", "array")`, and R 4.3 made that an error (`'length = 2' in coercion to 'logical(1)'`). It would fail on the main run's R 4.5.1 too. Fix: pass `X` after `oldClass(X) <- "matrix"`; values are unchanged.
- **`existing` needs exact row sums:** every row of the GPS matrix must satisfy `sum == 1` in floating point. 75 of our 400 rows were off by rounding; the largest entry in each was adjusted by at most 2.2e-16. Our own matching on the adjusted matrix gives identical estimates.

## Facts from the paper (checked against the PDF, 2026-11-15)

- Eq. (4) prints Pr(t | x_i, z_i) in both terms. Read with Eq. (3), the literal reading is GPS_lit2 (D13).
- "Code Availability N/A"; the paper names no software.
- Yang et al. (2016) is reference [44]; the methods "are reduced to that in [44]" without Z.
- The PDF is under copyright. Keep it out of the repo (`.gitignore` excludes `*.pdf`).

## Files added by this cross-check

```
BRIDGE.md                          this file
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
```

The script suggested by another assistant fails as written: a plain matrix `X` errors on R >= 4.3, `fit$estimate` partially matches `estimate_args`, `summary()` returns `NULL`, it does not flip signs, and it omits `existing`. `03_check_pasted.R` records this.

## How to reproduce

The cross-check ran on Linux with R 4.3.3, Matching 4.10-14, nnet 7.3-19, MASS 7.3-60.0.1 and reference BLAS/LAPACK 3.12.0. The main run used R 4.5.1 on Windows with nnet 7.3.20 and MASS 7.3.65; the regenerated sample still matched as described above.

In a Claude Code cloud session R is not pre-installed and CRAN was blocked by the network policy, so R came from the Ubuntu archive and the package from GitHub's read-only CRAN mirror. From the repo root:

```sh
apt-get install -y --no-install-recommends r-base-core r-cran-matching r-cran-mass r-cran-nnet r-cran-boot
git clone https://github.com/cran/multilevelMatching.git /tmp/multilevelMatching   # v1.0.0
mkdir -p ~/R/mlm-lib && R CMD INSTALL -l ~/R/mlm-lib /tmp/multilevelMatching
Rscript crosscheck/R/01_export_rep001.R . crosscheck/results > crosscheck/logs/01_export_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/02_pkg_fits.R crosscheck/results > crosscheck/logs/02_pkg_fits_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/03_check_pasted.R crosscheck/results > crosscheck/logs/03_check_pasted_log.txt 2>&1
```

On a laptop: `install.packages("Matching")`, then install `multilevelMatching` from CRAN if it is available there, or from source (`remotes::install_github("cran/multilevelMatching")`). Either way, on R >= 4.3 the `oldClass` workaround in `02_pkg_fits.R` is needed.

## Working rules carried over

- Leave the original repo, glw2022-matching-replication, unchanged.
- Claude does not commit or push unless asked; the spec says the owner commits from GitHub Desktop.
- Do not run the 500-replication cross-check until asked.
- Keep the paper's PDF out of the repo.

## Next steps

1. Done: the D7 and D13 wording is applied in `docs/spec.md`, `docs/replication_verdict.md`, `docs/sensitivity.md`, `docs/table3_checks.md` and `README.md`.
2. Update the inherited README and `CITATION.cff` (title, description, repository URL) for this repo.
3. When approved, extend the cross-check to all 500 replications of Simulation 1. Expect `existing` to match exactly (up to ties). Any `multinom` differences should be traced to the fitted probabilities (D5), not to D7.
4. Optional: add a note to D5 that `multilevelMatching` fits `nnet::multinom` at its default tolerance.
5. Optional: a Simulation 2 (9 levels) cross-check.

## References

- Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19: Multilevel covariates matching and propensity score matching. *Annals of Data Science, 9*(5), 967–982. https://doi.org/10.1007/s40745-022-00392-x
- Yang, S., Imbens, G., Cui, Z., Faries, D., & Kadziola, Z. (2016). Propensity score matching and subclassification in observational studies with multi-level treatments. *Biometrics, 72*, 1055–1065. https://doi.org/10.1111/biom.12505
- Yang, S., & Barkley, B. G. (2019). *multilevelMatching: Propensity score matching and subclassification in observational studies with multi-level treatments* (R package version 1.0.0). https://CRAN.R-project.org/package=multilevelMatching
