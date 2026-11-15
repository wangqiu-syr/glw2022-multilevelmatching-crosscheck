# Draft wording: D7 and D13 after the multilevelMatching cross-check

Drafts only; not yet applied to `docs/` or `README.md`. Paths refer to the glw2022-multilevelmatching-crosscheck layout (see BRIDGE.md). The software fact in D13 was checked against the paper's PDF on 2026-11-15.

## 1. docs/spec.md — decisions table

### D7 (replace the whole row)

| D7 | Eq. (4) uses p_t in both terms | Use p_t' for the t' term: column s of the imputed-outcome matrix is always matched on p_hat_s. This is how the R package of Yang et al. (2016), to whose methods the paper's reduce without Z, implements GPS matching: `multilevelMatching::multiMatch()` (version 1.0.0, CRAN 2019) imputes Y(s) by matching on column s of the GPS matrix. Cross-check on replication 1 of Simulation 1 (master seed 20220001, same sample as the main run): given our fitted GPS matrix (`match_on = "existing"`), `multiMatch()` reproduces our six GPS estimates and all 1,600 imputed potential outcomes exactly; with its own multinomial logit (`"multinom"`) and with covariate matching (`"covariates"`, Mahalanobis, pooled covariance) it also reproduces our GPS and CM estimates exactly. This shows that D7 is the estimator as implemented in the reference software for the method; it does not show what the paper's own code did. Run notes: R 4.3.3, Matching 4.10-14, nnet 7.3-19, no trimming, M_matches = 1, `set.seed(123)` (no ties occurred, so the seed did not affect the estimates). Version 1.0.0 stops on R >= 4.3 at its `class(X)` check, so `X` was passed with `oldClass(X) <- "matrix"` (values unchanged); `"existing"` requires every row to sum to exactly 1, so 75 of 400 rows were adjusted by at most 2.2e-16, which leaves our own estimates unchanged. Code and results: `crosscheck/R/01_export_rep001.R`, `crosscheck/R/02_pkg_fits.R`, `crosscheck/results/sim1_rep001_crosscheck.csv` |

### D13 (replace only the closing sentence of the row)

Old:

> This shows that the reading of Eq. (4) changes GPS bias in the direction of the paper's values; it does not establish which reading the paper used

New:

> This shows that the reading of Eq. (4) changes GPS bias in the direction of the paper's values; it does not establish which reading the paper used. Neither variant is what the package of Yang et al. (2016), `multilevelMatching`, computes: that package matches as in D7 and reproduces the main GPS exactly (see D7). If the paper's GPS results came from that package, the literal reading cannot explain the gap; it can do so only if the paper's own code followed Eq. (4) as printed. The paper names no software ("Code Availability N/A").

### Optional: GPS estimator paragraph ("Estimators" section)

After "The paper's Eq. (4) prints p_t in both terms, which looks like a typo **(D7)**." add:

> The R package of Yang et al. (2016), `multilevelMatching`, matches the same way; on a test replication it reproduces our GPS estimates exactly (D7).

## 2. docs/replication_verdict.md — decisions table and limits

D7 row:

| D7 | Eq. (4) uses p\_t in both terms of the effect | Use p\_t′ for the t′ term, treated as a typo. This is how Yang et al.'s (2016) R package multilevelMatching computes the estimator; given our fitted scores, it reproduces our GPS estimates exactly on a test replication. The literal reading is tested in D13 |

D13 row:

| D13 | Our GPS bias is smaller than the paper's | Sensitivity runs GPS\_lit and GPS\_lit2 (Simulation 1 only) that read Eq. (4) literally; the main GPS is unchanged. Neither variant is what multilevelMatching computes |

"What this does not show", first bullet:

> - **Which reading of Eq. (4) the paper used.** The literal reading explains part of the GPS gap, but it is not how Yang et al.'s own package computes the estimator, and we do not know whether the paper used that package or its own code.

## 3. Knock-on sentences (same claim elsewhere)

These sentences say the D7 reading "accounts for about half of the GPS gap". After the cross-check, that holds only if the paper's code read Eq. (4) literally. Suggested replacement, adapted to each spot:

> Reading Eq. (4) literally would close about half of the Simulation 1 GPS gap, but that is not how Yang et al.'s own package computes the estimator, so it explains the gap only if the paper's code followed Eq. (4) as printed.

Locations: README.md lines 11, 83 (replace "This is consistent with the published figures reflecting, at least in part, the literal reading."), 149, 166, 241; docs/replication_verdict.md lines 60 and 99; docs/sensitivity.md, "What this does and does not show".
