# Cross-checking a rebuild of Guo, Liu, and Wang (2022) with Yang et al.'s (2016) multilevelMatching package

Qiu Wang, Syracuse University · v1.0 · November 2026

This repository continues [glw2022-matching-replication](https://github.com/wangqiu-syr/glw2022-matching-replication), an R rebuild of the two simulation studies in Guo, Liu and Wang (2022), and keeps that rebuild unchanged. It checks the rebuild three ways:
- **Against Yang et al.'s package.** `multilevelMatching` implements the matching estimators of Yang et al. (2016), to which the paper's methods reduce without the school-level covariates.
- **Against an independent second rebuild,** written from the paper alone.
- **Against analytic intervals,** to see what the rebuild's percentile bootstrap intervals do and do not show.

As in the paper, the main criteria are bias and mean squared error (MSE); interval coverage comes last.

## At a glance

Mean absolute bias and mean MSE against τ\_pop, averaged over the 6 pairs of Simulation 1 and the 36 pairs of Simulation 2.
- **Paper:** its Tables 1 and 3. Table 1 reports no MSE.
- **Rebuild:** the main run.
- **Package:** `multilevelMatching` 1.0.0, run as a user would run it, on the rebuild's samples.
- **Rebuild-2:** the independent rebuild, on its own samples.

| Simulation | Estimator | Abs. bias: paper | rebuild | package | rebuild-2 | MSE: paper | rebuild | package | rebuild-2 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | CM | 0.143 | 0.187 | 0.187 | 0.166 |  | 0.363 | 0.363 | 0.364 |
| 1 | GPS | 0.293 | 0.026 | 0.026 | 0.043 |  | 0.676 | 0.676 | 0.656 |
| 2 | CM | 0.139 | 0.295 | 0.295 | 0.299 | 0.368 | 0.288 | 0.288 | 0.294 |
| 2 | GPS | 0.307 | 0.099 | 0.100 | 0.108 | 1.916 | 0.626 | 0.627 | 0.612 |

- **The package reproduces the rebuild.** With the same samples, Yang et al.'s own software gives the rebuild's bias and MSE: exactly for CM, and to within 0.001 for GPS.
- **An independent rebuild agrees, pair by pair.** Written from the paper alone, rebuild-2 matches the rebuild's Simulation 2 biases with correlations of 0.996 (CM) and 0.93 (GPS), differing by about the Monte Carlo error. It made the same choices on most of the points the paper leaves open.
- **GPS is less biased, CM has the lower MSE.** This bias–variance trade-off holds in the rebuild, the package and rebuild-2. The paper agrees on MSE but not on bias: its GPS bias is 3 to 11 times the rebuild's, while its CM bias is smaller.
- **So the gap with the paper is not in how the estimators are computed.** It lies in what the paper does not report: how its data were generated, which true effect its bias refers to, the GPS model, or its own code. The paper's CM column remains unexplained. It is not bias-adjusted matching.
- **Coverage depends on the interval method.**
  - **CM:** coverage is low (78% in Simulation 2) because of the bias from the covariate gaps between units and their matches, which no interval accounts for.
  - **GPS:** coverage is about 100% because the intervals ignore that the GPS was estimated; correcting for it (AI 2016) gives about 95%.
  - **The paper:** it does not say how its intervals were built. Its application's limits are exactly estimate ± 1.96 sd.

## Where to read more

| Topic | Document |
|---|---|
| The rebuild itself: background, how the paper handles the hierarchy, design, what differed from the paper, the two true effects τ\_pop and τ\_bal, results (Figures A–F) | [docs/rebuild.md](docs/rebuild.md) |
| The simulation specification and the thirteen decisions D1–D13 | [docs/spec.md](docs/spec.md) |
| What replicates and what does not; pair by pair against the paper's Tables 1 and 3 | [docs/replication_verdict.md](docs/replication_verdict.md), [docs/two_way.md](docs/two_way.md) |
| The literal reading of Eq. (4) (D13); exploratory checks of the paper's Table 3 | [docs/sensitivity.md](docs/sensitivity.md), [docs/table3_checks.md](docs/table3_checks.md) |
| The cross-check with `multilevelMatching`, all 500 replications of both simulations, with the package's facts and the workarounds it needs (Appendix B); why the three GPS versions compute the same estimator | [docs/crosscheck.md](docs/crosscheck.md), [docs/crosscheck_math.md](docs/crosscheck_math.md) |
| The paper, the rebuild and the package side by side: bias, MSE, the decisions D1–D13, analytic intervals (D8), what "multilevel" means | [crosscheck/side_by_side/README.md](crosscheck/side_by_side/README.md) |
| Intervals: where each method comes from, their math, why coverage is not accuracy, why CM coverage is low, with a math appendix | [crosscheck/intervals/README.md](crosscheck/intervals/README.md) |
| The independent rebuild-2, and the test of bias-adjusted CM | [crosscheck/rebuild2/README.md](crosscheck/rebuild2/README.md) |
| How this repository was set up: status, open questions, how to install R and the package in a cloud session, working rules, log | [BRIDGE.md](BRIDGE.md) |

## Reproduce

From the repository root (R 4.5; packages MASS and nnet):

| Step | Command | Time |
|---|---|---|
| 0. True effects (needed by the summaries) | `Rscript R/check_truth.R` | about 1 min |
| (optional) Tests | `Rscript R/test_estimators.R` (`quick` for unit tests only) | about 6 min |
| 1. Simulations 1 and 2, then summaries | `Rscript R/run_all.R` (options: `R=500 B=200 workers=10`) | On 10 workers: Simulation 1 about 15 min, Simulation 2 about 4 h 48 min, full run about 5 h (log: 14m28s, 4h48m19s, 5h03m08s) |
| 2. Rebuild summaries only | `Rscript R/summarize.R` | seconds |
| 3. Figures and comparisons | `Rscript R/plots.R` | seconds |
| 4. D13 sensitivity (Simulation 1) | `Rscript R/run_gps_literal.R`, then `Rscript R/run_gps_literal2.R` (`check` runs the self-tests) | |
| 5. Cross-check with multilevelMatching | `R_LIBS=<lib> Rscript crosscheck/R/04_crosscheck_all.R sim=1`, then `sim=2` (options: `workers=4`) | On 4 workers: Simulation 1 about 1 min, Simulation 2 about 18 min (log: 55 s, 18 min 16 s) |
| 6. Cross-check summaries | `Rscript crosscheck/R/05_summarize_crosscheck.R`, then `R_LIBS=<lib> Rscript crosscheck/R/06_existing_mismatches.R` | seconds; about 2 min |
| 7. Analytic variances for the rebuild's GPS (D8) | `R_LIBS=<lib> Rscript crosscheck/R/07_d8_analytic_intervals.R sim=1`, then `sim=2` (options: `workers=4`) | On 4 workers: Simulation 1 under 1 min, Simulation 2 about 12 min |
| 8. Side-by-side tables and figures | `Rscript crosscheck/R/08_side_by_side.R` | seconds |
| 9. CM bias decomposition and AI 2006 variance parts | `Rscript crosscheck/R/09_cm_bias_decomposition.R sim=1`, then `sim=2` | On 4 workers: seconds; about 1 min |
| 10. CM-only bootstrap rerun (interval centring) | `Rscript crosscheck/R/10_cm_bootstrap_centring.R sim=2 R=120`, then `sim=1 R=500` | On 4 workers: about 36 min; about 7 min |
| 11. Intervals tables and figure | `Rscript crosscheck/R/11_intervals_summary.R` | seconds |
| 12. Bias-adjusted CM (`Matching::Match(BiasAdjust = TRUE)`) | `Rscript crosscheck/R/12_cm_biasadjust.R sim=1`, then `sim=2` | On 4 workers: about 15 s; about 6 min |
| 13. Rebuild-2 comparison tables and figure | `Rscript crosscheck/R/13_rebuild2_comparison.R` | seconds |

**Checkpoints.** Each replication is saved to `results/reps/simK/rep_###.rds` when it finishes. Re-running `run_all.R` skips replications that are already done, so an interrupted run resumes where it stopped. It refuses to mix checkpoints made with different settings (R, B, seed, optimizer tolerance). Progress goes to standard output; redirect it to `results/run_log.txt` for a background run.

**Packages needed:**
- **Steps 9–11:** only R with MASS and nnet.
- **Step 12:** also `Matching`.
- **Steps 5–7:** `multilevelMatching` 1.0.0 and `Matching`. On R 4.3 or later the package needs an input workaround, which the scripts apply. Installation is in [BRIDGE.md](BRIDGE.md), and the workaround in [docs/crosscheck.md](docs/crosscheck.md), Appendix B. The cross-check ran on R 4.3.3; its replications regenerate the main run's samples from the same seeds.

**Random numbers and versions.** Each replication uses its own pre-generated L'Ecuyer-CMRG random-number stream, so results do not depend on the number of workers or the order in which replications run. The scripts use only base R, `parallel` and `nnet`, plus `MASS` for data generation. The R version and session details of the full run are in `results/run_session_info.txt`; the `nnet` and `MASS` versions used (7.3.20 and 7.3.65) are on the second line of `results/run_log.txt`. No `renv` lockfile is included.

## Folders

- **[R/](R/):** the rebuild's scripts, listed by role in [docs/rebuild.md](docs/rebuild.md#files-of-the-rebuild).
- **[results/](results/):** the main run's summaries and figures, and the values transcribed from the paper's Tables 1, 3 and 5 (`paper_table*.csv`).
- **[docs/](docs/):** the rebuild and cross-check documents in the table above.
- **[crosscheck/R/](crosscheck/R/):** cross-check scripts 01–13, run in the order of the Reproduce table.
- **[crosscheck/results/](crosscheck/results/)** and **[crosscheck/logs/](crosscheck/logs/):** their results and logs. The 16 per-replication result files are gzip-compressed CSVs (`.csv.gz`), which R's `read.csv()` reads as they are.
- **[crosscheck/side_by_side/](crosscheck/side_by_side/), [crosscheck/intervals/](crosscheck/intervals/), [crosscheck/rebuild2/](crosscheck/rebuild2/):** one note each, with its own tables/ and figures/.
- **Project files:** [LICENSE](LICENSE) (MIT), [CITATION.cff](CITATION.cff) (citation metadata) and [.gitignore](.gitignore), which keeps checkpoints, PDFs and local files out of the repo.

Per-replication checkpoints (`results/reps*/`, `crosscheck/reps/`) are local only, not in the repo. The paper's PDF is not included.

## Conclusions

In these simulations, covariate matching and GPS matching both remove most of the bias in the raw group difference, as the paper reports. GPS had the smallest bias and CM the smallest mean squared error, so the better choice depends on whether bias or variance matters more.

The published numbers cannot be reproduced exactly from the paper alone, because thirteen points are left open. Reading Eq. (4) literally would close about half of the GPS discrepancy in Simulation 1. But that is not how Yang et al.'s own package computes the estimator, so it explains the gap only if the paper's code followed Eq. (4) as printed.

These conclusions do not rest on one implementation or one reading of the paper:
- Yang et al.'s own package, run on the rebuild's samples, gives the rebuild's bias and MSE.
- An independent rebuild, written from the paper alone, gives the same results pair by pair.

So the remaining differences from the paper lie in what the paper does not report, chiefly how its data were generated and which true effect its bias refers to. The paper's reported coverage is closer to analytic Abadie–Imbens intervals than to the rebuild's bootstrap intervals. CM's low coverage comes from the covariate gaps between units and their matches, which no interval accounts for. These results are a baseline for extending the comparison to clustered, multilevel data.

## Acknowledgment

This is a re-implementation, written from the published description of Guo, Liu & Wang (2022). The original authors' code was not used.

Disclosure: the author of this repository is a coauthor of Guo, Liu and Wang (2022).

## AI assistance

The code was written with assistance from Claude Code (Anthropic).

## References

Each document has its own references. Cited here:

- Abadie, A., & Imbens, G. W. (2006). Large sample properties of matching estimators for average treatment effects. *Econometrica, 74*(1), 235–267.
- Abadie, A., & Imbens, G. W. (2016). Matching on the estimated propensity score. *Econometrica, 84*(2), 781–807.
- Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19: Multilevel covariates matching and propensity score matching. *Annals of Data Science, 9*(5), 967–982. https://doi.org/10.1007/s40745-022-00392-x
- Wang, Q. (2026). *Generalized propensity score and covariate matching for multilevel treatments in hierarchical data: Rebuilding Guo, Liu, and Wang (2022) in R* (Version 1.0) [Computer software]. https://github.com/wangqiu-syr/glw2022-matching-replication
- Yang, S., & Barkley, B. G. (2019). *multilevelMatching: Propensity score matching and subclassification in observational studies with multi-level treatments* (R package version 1.0.0). https://CRAN.R-project.org/package=multilevelMatching
- Yang, S., Imbens, G., Cui, Z., Faries, D., & Kadziola, Z. (2016). Propensity score matching and subclassification in observational studies with multi-level treatments. *Biometrics, 72*, 1055–1065. https://doi.org/10.1111/biom.12505

## How to cite

Wang, Q. (2026). *Cross-checking a rebuild of Guo, Liu, and Wang (2022) with Yang et al.'s (2016) multilevelMatching package* (Version 1.0) [Computer software]. https://github.com/wangqiu-syr/glw2022-multilevelmatching-crosscheck

DOI: to be added after the Zenodo release. Citation metadata are in [`CITATION.cff`](CITATION.cff); please also cite the original paper, Yang et al. (2016) and the rebuild it builds on.

## License

MIT (see `LICENSE`).
