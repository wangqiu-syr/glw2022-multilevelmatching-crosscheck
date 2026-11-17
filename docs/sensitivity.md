# Sensitivity check D13: reading Eq. (4) literally

Simulation 1 only (4 treatment combinations, 6 pairwise effects, 500 replications). Sources: `results/gps_literal_sim1.csv`, `results/gps_literal2_sim1.csv` and the D13 row of `docs/spec.md`.

## What is tested

The paper's Eq. (4) (Guo, Liu and Wang, 2022, Section 4, p. 972) prints the same score, Pr(t | x, z), in both terms of the GPS effect estimate:

```math
\hat\tau_{GPS}(t,t') = \frac{1}{n} \sum_{i=1}^{n} \left\lbrace Y^{obs}_{GPS}\bigl(t, \Pr(t \mid x_i, z_i)\bigr) - Y^{obs}_{GPS}\bigl(t', \Pr(t \mid x_i, z_i)\bigr) \right\rbrace \qquad (4)
```

Here $`Y^{obs}_{GPS}(t, p)`$ is the observed outcome of the unit in group t whose score for t is closest to p, as defined by the paper's Eq. (3):

```math
GPS(t,p) = \mathop{\arg\min}_{j:\ T_j = t} \left\lVert \Pr(t \mid x_j, z_j) - p \right\rVert \qquad (3)
```

In words: each unit's outcome under t and under t′ is imputed from its closest match in each group, and the differences are averaged over all n units. The main GPS estimator reads the second Pr(t | x_i, z_i) as a typo for Pr(t′ | x_i, z_i) and matches the t′ term on the score for t′ (D7). This check tests the literal reading instead:

- **GPS_lit.** For a unit in neither group t nor t′, the t′ term is matched within group t′ on |p̂_t′(j) − p̂_t(i)|. Units in group t or t′ keep their observed outcome for their own level; their other term is imputed as in the main GPS.
- **GPS_lit2.** The fully literal reading: units in group t also get their t′ term by matching within group t′ on |p̂_t′(j) − p̂_t(i)|. Every unit keeps its observed outcome for its own group's term.

Both variants use the same seeds, datasets and bootstrap draws as the main run (paired). All 500 replications reproduce the main DIF, CM and GPS results exactly, and the GPS_lit results are reproduced exactly in the GPS_lit2 run. The main GPS estimator and all main results are unchanged.

## Results

Absolute bias against τ\_pop, by pair. The paper's column is its Table 1 GPS value.

| Pair | Paper GPS | GPS | GPS_lit | GPS_lit2 |
| --- | --- | --- | --- | --- |
| (1,2) | 0.0597 | 0.024 | 0.036 | 0.161 |
| (1,3) | 0.3742 | 0.037 | 0.025 | 0.255 |
| (1,4) | 0.4811 | 0.048 | 0.005 | 0.142 |
| (2,3) | 0.3145 | 0.013 | 0.027 | 0.224 |
| (2,4) | 0.4214 | 0.024 | 0.103 | 0.141 |
| (3,4) | 0.1069 | 0.011 | 0.066 | 0.030 |

Mean over the 6 pairs (D13 row of `docs/spec.md`): paper GPS 0.293, main GPS 0.026, GPS_lit 0.044, GPS_lit2 0.159.

MSE and 95% interval coverage (%) against τ\_pop. The "GPS_lit2 vs GPS" column is the MSE increase computed from the two MSE columns.

| Pair | MSE GPS | MSE GPS_lit | MSE GPS_lit2 | GPS_lit2 vs GPS | Coverage GPS_lit | Coverage GPS_lit2 | Paper coverage |
| --- | --- | --- | --- | --- | --- | --- | --- |
| (1,2) | 0.476 | 0.650 | 0.826 | +74% | 100.0 | 99.6 | 100 |
| (1,3) | 0.297 | 0.468 | 0.694 | +134% | 100.0 | 99.2 |  99 |
| (1,4) | 0.832 | 1.323 | 1.593 | +92% | 100.0 | 99.8 | 100 |
| (2,3) | 0.630 | 0.814 | 0.942 | +50% | 100.0 | 100.0 |  99 |
| (2,4) | 0.902 | 1.212 | 1.352 | +50% | 100.0 | 100.0 | 100 |
| (3,4) | 0.922 | 1.179 | 1.489 | +62% | 100.0 | 100.0 | 100 |

- **GPS_lit** changes the bias little (mean 0.044 against 0.026).
- **GPS_lit2** raises the mean absolute bias to 0.159, about half the way from the main GPS to the paper's 0.293. It is nearer the paper's value in 5 of 6 pairs, Its bias is beyond Monte Carlo noise in 3 of 6 pairs (|z| > 3; 5 of 6 at 2 SE).
- **Cost.** GPS_lit2's MSE is 50–134% higher than the main GPS, and its coverage stays at 99.2–100.0%.

## What this does and does not show

This tests one reading of Eq. (4). It shows that the literal reading moves the GPS bias toward the paper's values and would account for part of the gap. It does not show which reading the paper used, and the rest of the gap remains unexplained. Neither variant is what Yang et al.'s (2016) R package `multilevelMatching` computes: given our fitted scores, that package reproduces the main GPS exactly on a test replication (D7 in `docs/spec.md`). The literal reading explains the gap only if the paper's code followed Eq. (4) as printed; the paper names no software ("Code Availability N/A").

Code: `R/gps_literal.R`, `R/run_gps_literal.R`, `R/gps_literal2.R`, `R/run_gps_literal2.R`.
