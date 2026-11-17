# Three readings of Eq. (4), Simulation 1

The paper's GPS estimator, Eq. (4) (Guo, Liu and Wang, 2022, Section 4, p. 972), prints the score for t, Pr(t | x_i, z_i), in both of its terms:

```math
\hat\tau_{GPS}(t,t') = \frac{1}{n} \sum_{i=1}^{n} \left\lbrace Y^{obs}_{GPS}\bigl(t, \Pr(t \mid x_i, z_i)\bigr) - Y^{obs}_{GPS}\bigl(t', \Pr(t \mid x_i, z_i)\bigr) \right\rbrace \qquad (4)
```

Here $`Y^{obs}_{GPS}(t, p)`$ is the observed outcome of the unit picked by Eq. (3):

```math
GPS(t,p) = \mathop{\arg\min}_{j:\ T_j = t} \left\lVert \Pr(t \mid x_j, z_j) - p \right\rVert \qquad (3)
```

The first term is clear: each unit's outcome under t comes from its closest match in group t on the score for t. The second term, the outcome under t′, can be read three ways. The rebuild's main GPS and its D13 check cover two of them ([docs/sensitivity.md](../../docs/sensitivity.md)). This note adds the third: one score, the score for t, used for both groups, as in binary propensity-score matching.

| Reading | Unit i's match in group t′ is the unit j that minimizes | Where |
|---|---|---|
| **Main GPS** (D7) | \|p̂\_t′(j) − p̂\_t′(i)\|: the score for t′ on both sides | the rebuild; also what Yang et al.'s package computes |
| **Literal, as printed** (GPS\_lit2) | \|p̂\_t′(j) − p̂\_t(i)\|: Eq. (3) for group t′, given unit i's score for t | D13 |
| **Literal, outside units only** (GPS\_lit) | as GPS\_lit2, but only for units in neither group t nor t′ | D13 |
| **Same score for t on both sides** (GPS\_same) | \|p̂\_t(j) − p̂\_t(i)\|: the score for t on both sides | new, this note |

In every reading the t term is the main GPS's, and each unit keeps its observed outcome for its own group. GPS\_same is not symmetric in t and t′, because it uses the score for t, so pairs are taken as in the paper's Table 1, with t < t′.

## Results

Simulation 1, the main run's 500 samples, point estimates only: there is no bootstrap, so there is no coverage. Bias is against τ\_pop unless stated. "–": not reported or not applicable.

| Reading | Mean abs. bias vs τ\_pop | Mean abs. bias vs τ\_bal | Mean MSE | Mean SD | Pairs nearer the paper than the main GPS | Correlation with the paper, pair by pair | Pairs with bias beyond 3 MC SE |
|---|---|---|---|---|---|---|---|
| Paper (Table 1) | 0.293 | – | – | – | – | – | – |
| Main GPS (D7) | 0.026 | 0.072 | 0.676 | 0.809 | – | 0.620 | 0 |
| Literal, outside units only (GPS\_lit) | 0.044 | 0.080 | 0.941 | 0.954 | 4 | −0.126 | 0 |
| Literal, as printed (GPS\_lit2) | 0.159 | 0.205 | 1.150 | 1.045 | 5 | 0.405 | 3 |
| Same score for t on both sides (GPS\_same) | 0.313 | 0.310 | 1.168 | 1.003 | 4 | 0.157 | 6 |

Absolute bias by pair, against τ\_pop:

| Pair | Paper | Main GPS | GPS\_lit | GPS\_lit2 | GPS\_same |
|---|---|---|---|---|---|
| (1,2) | 0.060 | 0.024 | 0.036 | 0.161 | 0.161 |
| (1,3) | 0.374 | 0.037 | 0.025 | 0.255 | 0.156 |
| (1,4) | 0.481 | 0.048 | 0.005 | 0.142 | 0.326 |
| (2,3) | 0.315 | 0.013 | 0.027 | 0.224 | 0.560 |
| (2,4) | 0.421 | 0.024 | 0.103 | 0.141 | 0.308 |
| (3,4) | 0.107 | 0.011 | 0.066 | 0.030 | 0.366 |

## What this shows

- **The same-score reading matches the paper's average.** Its mean absolute bias, 0.313, is the closest of the four readings to the paper's 0.293.
- **It does not match the paper pair by pair.** Its correlation with the paper's six values is 0.16, against 0.62 for the main GPS, and it is nearer the paper than the main GPS in only 4 of 6 pairs. Its largest misses are (2,3), 0.56 against 0.31, and (3,4), 0.37 against 0.11.
- **It is the least accurate reading.** Its mean MSE is 1.17, against 0.68 for the main GPS, and all six of its biases are beyond Monte Carlo noise.
- **Why it is biased.** With more than two groups, units with the same score for t can have very different scores for t′. So matching group t′ on the score for t does not balance the covariates that drive assignment to t′. The paper's own Section 2 (p. 970) says the outcome under t is identified within subpopulations "defined by the value of a single score, Pr(t | x, z)"; for the outcome under t′, that is Pr(t′ | x, z). With two groups, Pr(t′ | x, z) = 1 − Pr(t | x, z), so one score serves both groups; with four or nine it does not.

## What this does not show

- **Which reading the paper used.** A mean level over six pairs is weak evidence, and the pair-by-pair pattern points away from this reading. The paper names no software ("Code Availability N/A"), and Yang et al.'s package computes the main GPS.
- **Simulation 2 or coverage.** Only Simulation 1 was run, and without the bootstrap.

## Reproduce

`Rscript crosscheck/R/14_gps_same_score.R` (options: `R=500 workers=4`), from the repository root, with R, MASS and nnet. It takes about 15 seconds on 4 workers.

Each replication regenerates the main run's sample from its stored seed, then computes CM and GPS in the main run's order. Checks against the stored results:
- **DIF, CM and the main GPS:** they reproduce exactly (largest difference in a pair mean, 4e-15), and the script stops if they do not.
- **GPS\_lit and GPS\_lit2:** they reproduce exactly in 5 of 6 pairs. In pair (1,4), both means differ from the stored D13 run by 3.1e-5. That run is from another setup (R 4.5.1 on Windows), whose fitted scores differ from this one's at about the 1e-9 level ([docs/spec.md](../../docs/spec.md), D5). At a near-tie that can change a literal-reading match; the main GPS is not affected.

## Files

- **Script:** [crosscheck/R/14_gps_same_score.R](../R/14_gps_same_score.R). **Log:** `crosscheck/logs/14_gps_same_score_log.txt`.
- **Per-replication estimates:** `crosscheck/results/sim1_gps_same_score_reps.csv.gz` (DIF, CM, GPS, GPS\_lit, GPS\_lit2 and GPS\_same, by replication and pair).
- **Tables** ([tables/](tables/)): `sim1_eq4_readings_summary`, `sim1_eq4_readings_abs_bias_by_pair` (CSV and Markdown), and `sim1_eq4_readings_by_pair.csv`, with bias against both truths, SD, Monte Carlo SE and MSE for each reading and pair.

## References

- Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19: Multilevel covariates matching and propensity score matching. *Annals of Data Science, 9*(5), 967–982. https://doi.org/10.1007/s40745-022-00392-x
- Yang, S., Imbens, G., Cui, Z., Faries, D., & Kadziola, Z. (2016). Propensity score matching and subclassification in observational studies with multi-level treatments. *Biometrics, 72*, 1055–1065. https://doi.org/10.1111/biom.12505
