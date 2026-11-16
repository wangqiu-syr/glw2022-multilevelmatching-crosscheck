# Rebuild-2 summary results (copied)

These four files are copied unchanged from the independent rebuild, the private repository `wangqiu-syr/glw2022-rebuild2`, at commit `7b25857` on `main` (Merge pull request #1, 2026-11-16). They are the same at `8690904` on its branch `claude/new-session-2wp49a`, which adds only figures.

| File here | Source in rebuild-2 |
|---|---|
| `rebuild2_sim1_performance.csv`, `rebuild2_sim2_performance.csv` | `results/sim1/performance.csv`, `results/sim2/performance.csv`: bias, absolute bias, MSE, SD and Monte Carlo SEs by estimator, pair and truth (`tau_super`, `tau_mix`) |
| `rebuild2_sim1_coverage.csv`, `rebuild2_sim2_coverage.csv` | `results/sim1/coverage.csv`, `results/sim2/coverage.csv`: coverage of its analytic 95% intervals (its D26) against `tau_super`, `tau_mix` and the sample-average effect; oracle coverage; mean SE and Monte Carlo SD |

Rebuild-2's `tau_super` and `tau_mix` are this repository's τ\_pop and τ\_bal. Its pairs are labelled `t`, `tprime`. Estimators: DIF, CM, GPS, GPS\_cs (correctly specified logit) and GPS\_true (true scores).
