# Bridge: from glw2022-matching-replication to glw2022-multilevelmatching-crosscheck

This file lets the next session, or reader, pick up without the chat history. It records:
- what this repository is;
- where each piece of work stands;
- what is still open;
- how to set up R in a cloud session;
- the working rules;
- a short log.

Written 2026-11-15. Condensed 2026-11-17: the detailed cross-check results, the package facts and workarounds, and the file list moved, unchanged, to [docs/crosscheck.md](docs/crosscheck.md) (Appendices A–C).

## What this repo is

This repo starts as an import of [glw2022-matching-replication](https://github.com/wangqiu-syr/glw2022-matching-replication), as of commit `aea8268` ("Add exploratory Table 3 checks; relabel figures; clarify Simulation 1 scope"). The import is a single initial commit, `4bfdd50`, with the same 48 files and identical contents, without the earlier history, which stays in the original repo. The original repo is left unchanged.

This repo checks the rebuild three ways:
- **Against Yang et al.'s package.** The rebuild's covariate matching (CM) and generalized propensity score (GPS) matching are compared with the R package `multilevelMatching` of Yang et al. (2016). Guo, Liu & Wang (2022) state that without the school-level covariates Z their methods reduce to Yang et al. (2016), their reference [44]. In Yang et al., "multilevel" means a treatment with several levels, not nested data.
- **Against analytic intervals,** in a study of the interval methods.
- **Against an independent second rebuild** written from the paper alone.

The [README](README.md) has the summary and a map of all documents.

## Where things stand

| Work | Document | Merged in |
|---|---|---|
| The rebuild (imported, unchanged) | [docs/rebuild.md](docs/rebuild.md), [docs/spec.md](docs/spec.md) and the other docs/ files | import `4bfdd50` |
| Cross-check of replication 1, then all 500 replications of both simulations | [docs/crosscheck.md](docs/crosscheck.md) | PRs #1, #3 |
| D7 and D13 wording after the cross-check | [docs/spec.md](docs/spec.md), [docs/replication_verdict.md](docs/replication_verdict.md), [docs/sensitivity.md](docs/sensitivity.md), [docs/table3_checks.md](docs/table3_checks.md); drafts in [crosscheck/drafts/](crosscheck/drafts/) | PR #2 |
| Why the three GPS versions compute the same estimator | [docs/crosscheck_math.md](docs/crosscheck_math.md) | PR #4 |
| Paper, rebuild and package side by side; D8 analytic-interval sensitivity | [crosscheck/side_by_side/README.md](crosscheck/side_by_side/README.md) | PRs #5, #6, #7 |
| Intervals: sources, math, why CM coverage is low; math appendix; the paper's interval facts from its PDF | [crosscheck/intervals/README.md](crosscheck/intervals/README.md) | PRs #8, #9 |
| Rebuild-2 comparison; bias-adjusted CM test | [crosscheck/rebuild2/README.md](crosscheck/rebuild2/README.md) | PRs #10, #11 |
| README slimmed to a summary and map; this file condensed | [README.md](README.md), [docs/rebuild.md](docs/rebuild.md), this file | branch `claude/slim-readme-bridge` |

Rebuild-2 lives in its own private repository, `wangqiu-syr/glw2022-rebuild2`: its stages 1–4 are in its PR #1, and its figures commit in its PR #2 (merge commit `0a4a038`). This repository only reads it. Its four summary CSVs are copied into [crosscheck/rebuild2/data/](crosscheck/rebuild2/data/) from `7b25857`, and they are unchanged at `0a4a038`.

## Open questions

- **What the paper's CM column is.**
  - **Bias:** it is less biased than plain CM. Its biases fit 1.38 × the sample-mix term + 0.22 × the matching-gap term (R² = 0.93).
  - **Variance:** it is about twice as variable as plain CM.
  - **What doesn't fit:** neither plain CM, nor bias-adjusted CM (`Matching::Match(BiasAdjust = TRUE)`), nor an exchange of the CM and GPS labels fits both its bias and its MSE ([crosscheck/rebuild2/README.md](crosscheck/rebuild2/README.md), section 5).
- **How the paper built its intervals.** It does not say. Its application's 95% limits are exactly estimate ± 1.96 sd, so they are normal intervals, not percentile ones. Its simulation coverage is closer to analytic Abadie–Imbens (2006) intervals than to the rebuild's bootstrap (side-by-side note, section 5.2).
- **The DIF column.** The rebuild's DIF bias tracks the paper's (with the sign flipped) but does not match it value for value. Rebuild-2 found the same and could not trace the cause (its D20).

## Facts from the paper (checked against the PDF, 2026-11-15 and 2026-11-16)

- Eq. (4) prints Pr(t | x_i, z_i) in both terms. Read with Eq. (3), the literal reading is GPS_lit2 (D13).
- **Nothing on how intervals were computed:** the words "bootstrap" and "Abadie" do not occur, no software is named, and the paper never says how its intervals or standard deviations were computed (D8).
- **Table 5's limits are symmetric normal intervals:** each 95% limit equals the estimate ∓ qnorm(0.975) × "sd" ("the sample standard deviation") to within rounding, for all 30 CM and GPS intervals. These are not percentile intervals (`crosscheck/intervals/tables/paper_table5_interval_form.md`).
- "Data Availability N/A" and "Code Availability N/A"; the paper names no software.
- Yang et al. (2016) is reference [44]; the methods "are reduced to that in [44]" without Z.
- The PDF is under copyright. Keep it out of the repo (`.gitignore` excludes `*.pdf`).

## Setting up R and the package in a cloud session

The cross-check ran on Linux with R 4.3.3, Matching 4.10-14, nnet 7.3-19, MASS 7.3-60.0.1 and reference BLAS/LAPACK 3.12.0. The main run used R 4.5.1 on Windows with nnet 7.3.20 and MASS 7.3.65; the regenerated samples still match ([docs/crosscheck.md](docs/crosscheck.md)).

In a Claude Code cloud session R is not pre-installed and CRAN was blocked by the network policy, so R came from the Ubuntu archive and the package from GitHub's read-only CRAN mirror. From the repo root:

```sh
apt-get install -y --no-install-recommends r-base-core r-cran-matching r-cran-mass r-cran-nnet r-cran-boot
git clone https://github.com/cran/multilevelMatching.git /tmp/multilevelMatching   # v1.0.0
mkdir -p ~/R/mlm-lib && R CMD INSTALL -l ~/R/mlm-lib /tmp/multilevelMatching
Rscript crosscheck/R/01_export_rep001.R . crosscheck/results > crosscheck/logs/01_export_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/02_pkg_fits.R crosscheck/results > crosscheck/logs/02_pkg_fits_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/03_check_pasted.R crosscheck/results > crosscheck/logs/03_check_pasted_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/04_crosscheck_all.R sim=1 workers=4 > crosscheck/logs/04_crosscheck_sim1_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/04_crosscheck_all.R sim=2 workers=4 > crosscheck/logs/04_crosscheck_sim2_log.txt 2>&1
Rscript crosscheck/R/05_summarize_crosscheck.R > crosscheck/logs/05_summarize_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/06_existing_mismatches.R > crosscheck/logs/06_existing_mismatches_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/07_d8_analytic_intervals.R sim=1 workers=4 > crosscheck/logs/07_d8_sim1_log.txt 2>&1
R_LIBS=~/R/mlm-lib Rscript crosscheck/R/07_d8_analytic_intervals.R sim=2 workers=4 > crosscheck/logs/07_d8_sim2_log.txt 2>&1
Rscript crosscheck/R/08_side_by_side.R > crosscheck/logs/08_side_by_side_log.txt 2>&1
for s in 1 2; do Rscript crosscheck/R/09_cm_bias_decomposition.R sim=$s workers=4 > crosscheck/logs/09_cm_bias_decomposition_sim${s}_log.txt 2>&1; done
Rscript crosscheck/R/10_cm_bootstrap_centring.R sim=2 R=120 workers=4 > crosscheck/logs/10_cm_bootstrap_centring_sim2_log.txt 2>&1
Rscript crosscheck/R/10_cm_bootstrap_centring.R sim=1 R=500 workers=4 > crosscheck/logs/10_cm_bootstrap_centring_sim1_log.txt 2>&1
Rscript crosscheck/R/11_intervals_summary.R
for s in 1 2; do Rscript crosscheck/R/12_cm_biasadjust.R sim=$s workers=4 > crosscheck/logs/12_cm_biasadjust_sim${s}_log.txt 2>&1; done
Rscript crosscheck/R/13_rebuild2_comparison.R
```

On a laptop: `install.packages("Matching")`, then install `multilevelMatching` from CRAN if it is available there, or from source (`remotes::install_github("cran/multilevelMatching")`). Either way, on R >= 4.3 the `oldClass` workaround in `02_pkg_fits.R` is needed ([docs/crosscheck.md](docs/crosscheck.md), Appendix B.2). Other hosts the sessions needed and could not reach: arxiv.org, link.springer.com, sticerd.lse.ac.uk (full texts of papers).

## Working rules carried over

- Leave the original repo, glw2022-matching-replication, unchanged.
- Claude does not commit or push unless asked; the spec says the owner commits from GitHub Desktop. In practice, since 2026-11-15, the owner has asked Claude to commit on `claude/...` branches, open pull requests, and merge them when told to.
- Re-running a script resumes from its checkpoints in `results/reps*/` or `crosscheck/reps/`, which are git-ignored and local only.
- Keep the paper's PDF out of the repo.
- Rebuild-2 is read here, never written to from this repository's work. While it ran, nothing was passed back to it.

## Log

Dates are UTC, as GitHub records the merges.

- **2026-11-15:**
  - Repository imported (`4bfdd50`).
  - Cross-check of replication 1 of Simulation 1 (PR #1).
- **2026-11-16:**
  - D7 and D13 wording (PR #2).
  - Cross-check of all 500 replications of both simulations, with the README, CITATION.cff and D5 updated (PR #3).
  - Why the three GPS versions agree (PR #4).
  - Side-by-side note and D8 analytic intervals (PR #5).
  - README reorganized around the side-by-side results (PR #6).
  - D8 in the spec and the verdict (PR #7).
  - Intervals note (PR #8) and its math appendix, with the paper's interval facts from the PDF (PR #9).
  - Rebuild-2 comparison and the bias-adjusted CM test (PR #10). The blind prediction planned for rebuild-2 was not made: it had finished, and its log already quoted results, before it was read here.
  - In rebuild-2's repository, its figures commit was merged as its PR #2.
- **2026-11-17:**
  - BRIDGE updated for those merges (PR #11).
  - README slimmed to a summary and map. Its rebuild sections moved unchanged to `docs/rebuild.md`; its side-by-side, intervals and rebuild-2 sections are covered by their notes.
  - This file condensed, with its detailed sections moved unchanged to `docs/crosscheck.md` (Appendices A–C).
  - The 16 per-replication CSVs in `crosscheck/results/` gzip-compressed (`.csv.gz`), taking them from 12.7 MB to 4.7 MB and the repository from 15.1 MB to 7.1 MB; the writer scripts (04, 07, 09, 10, 12) now write them that way.
    - **Content:** each file decompresses to exactly the committed CSV.
    - **Readers:** scripts 05, 06, 08, 11 and 13, rerun on the compressed files, reproduce all their outputs unchanged.
    - **Writers:** 09, 10 and 12, rerun, rewrite identical content. 04 and 07 were tested on 2 replications in a copy of the repo.
    - **Logs:** the earlier run logs keep the old file names.

## References

- Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19: Multilevel covariates matching and propensity score matching. *Annals of Data Science, 9*(5), 967–982. https://doi.org/10.1007/s40745-022-00392-x
- Yang, S., Imbens, G., Cui, Z., Faries, D., & Kadziola, Z. (2016). Propensity score matching and subclassification in observational studies with multi-level treatments. *Biometrics, 72*, 1055–1065. https://doi.org/10.1111/biom.12505
- Yang, S., & Barkley, B. G. (2019). *multilevelMatching: Propensity score matching and subclassification in observational studies with multi-level treatments* (R package version 1.0.0). https://CRAN.R-project.org/package=multilevelMatching
