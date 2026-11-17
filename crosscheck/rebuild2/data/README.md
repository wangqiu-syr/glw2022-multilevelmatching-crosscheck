# Rebuild 2 summary results and table transcriptions (copied)

These four files are copied unchanged from the independent rebuild, the repository [`wangqiu-syr/glw2022-rebuild2`](https://github.com/wangqiu-syr/glw2022-rebuild2), at commit `7b25857` on `main` (Merge pull request #1, 2026-11-16). They are unchanged at `706d514`, its `main` on 2026-11-17: its later PRs #2–#4 add figures, a provenance note and a new long-format summary file only.

| File here | Source in rebuild 2 |
|---|---|
| `rebuild2_sim1_performance.csv`, `rebuild2_sim2_performance.csv` | `results/sim1/performance.csv`, `results/sim2/performance.csv`: bias, absolute bias, MSE, SD and Monte Carlo SEs by estimator, pair and truth (`tau_super`, `tau_mix`) |
| `rebuild2_sim1_coverage.csv`, `rebuild2_sim2_coverage.csv` | `results/sim1/coverage.csv`, `results/sim2/coverage.csv`: coverage of its analytic 95% intervals (its D26) against `tau_super`, `tau_mix` and the sample-average effect; oracle coverage; mean SE and Monte Carlo SD |

The two transcriptions of the paper's tables are copied unchanged from rebuild 2's `main` at `f8050d7` (2026-11-17). They have not changed since its Stage 1 commit `76c84de`, so they are also what `7b25857` holds:

| File here | Source in rebuild 2 |
|---|---|
| `rebuild2_paper_table1_sim1.csv` | `data/paper/table1_sim1.csv`: the paper's Table 1 (Simulation 1), as rebuild 2 transcribed it |
| `rebuild2_paper_table3_sim2.csv` | `data/paper/table3_sim2.csv`: the paper's Table 3 (Simulation 2), as rebuild 2 transcribed it |

`crosscheck/R/13_rebuild2_comparison.R` compares them cell by cell with this repository's own transcription, `results/paper_table1.csv` and `results/paper_table3.csv`.

Rebuild 2's `tau_super` and `tau_mix` are this repository's τ\_pop and τ\_bal. Its pairs are labelled `t`, `tprime`. Estimators: DIF, CM, GPS, GPS\_cs (correctly specified logit) and GPS\_true (population scores: the propensity scores of the data-generating model).
