# Replication of Guo, Liu & Wang (2022): results and decisions D1–D13

Nov 4, 2026

Our R replication reproduces the paper's main pattern: both matching methods remove most of the bias in the raw group difference. It does not reproduce the paper's exact numbers, and our GPS estimator is far less biased than the paper's, so the thirteen choices below, where the paper is silent, matter.

## Why replicate this paper

Guo, Liu & Wang (2022) propose covariate matching and generalized propensity score (GPS) matching for hierarchical data with multi-level treatments, a setting common in education research. The paper provides no code or data, and several simulation details are not specified. We rebuilt the simulations from the published description to check which reported results hold up, and which depend on choices the paper leaves open (D1–D13). Code and results are provided so that each choice can be checked.

## Main findings

- **The estimators are the paper's, not modified.** CM and GPS are implemented as described. The one departure from the printed text is Eq. (4), where the main GPS uses the score for t′ in the second term (D7); the literal reading is reported as a sensitivity check (D13).
- **Matching removes most of the raw-difference bias.** Mean absolute bias is 0.6–0.9 for DIF, 0.2–0.3 for CM and 0.03–0.10 for GPS.
- **There is a bias–variance trade-off.** GPS has the smallest bias but about twice CM's MSE and intervals about 30% too wide; CM is more efficient but more biased.
- **Our GPS bias is much smaller than the published values.** In Simulation 1 it averages 0.026, against 0.29 in the paper's Table 1, and the Simulation 2 bias pattern tracks the paper's with the opposite sign. This is a difference between implementations, not an improvement of the method, and we cannot tell which unreported choices explain it.

## What we did

We re-ran both simulations from Guo, Liu & Wang (2022, *Annals of Data Science* 9(5):967–982, [doi:10.1007/s40745-022-00392-x](https://doi.org/10.1007/s40745-022-00392-x)) in R, with 500 replications each.

- **Simulation 1:** 4 treatment combinations, 100 units each, 6 pairwise effects.
- **Simulation 2:** 9 treatment combinations, 200 units each, 36 pairwise effects.
- **Estimators:** the raw difference in group means (DIF); covariate matching (CM: one nearest neighbour by Mahalanobis distance, which scales the covariates by their variances and allows for their correlations, with replacement); generalized propensity score matching (GPS: linear multinomial logit, nearest neighbour on the fitted score). The GPS model is fitted much more precisely than the software's default (D5), because the default left the fitted probabilities too imprecise for the very small matching distances.
- **Intervals:** percentile bootstrap, B = 200, resampling within treatment groups.
- **Run:** about 5 hours on 10 workers. Code, results and the full specification (docs/spec.md) are in the repository wangqiu-syr/glw2022-matching-replication.

## Results

Averaged over pairs, DIF is biased by 0.6–0.9, CM by 0.2–0.3 and GPS by 0.03–0.10. Both matching methods remove most of the raw-difference bias.

Mean absolute bias (95% CI coverage), averaged over pairs. Truth is the population ATE (τ\_pop) or the balanced-design ATE (τ\_bal).

| Estimator | Sim 1, vs τ\_pop | Sim 1, vs τ\_bal | Sim 2, vs τ\_pop | Sim 2, vs τ\_bal |
| --- | --- | --- | --- | --- |
| DIF | 0.588 (85%) | 0.612 (84%) | 0.941 (78%) | 0.927 (79%) |
| CM | 0.187 (91%) | 0.223 (90%) | 0.295 (78%) | 0.262 (81%) |
| GPS | 0.026 (100%) | 0.072 (100%) | 0.099 (100%) | 0.096 (100%) |

- **The two truths differ.** Fixed group sizes change the covariate mix. In Simulation 2 all 36 population effects are 0, while the balanced-design effects range from −0.18 to 0.17; in Simulation 1 the gap is up to about 0.09.
- **CM is more efficient.** GPS has about double CM's MSE (against τ\_pop; Simulation 1: 0.68 vs 0.36; Simulation 2: 0.63 vs 0.29), because its estimates vary more (SD 0.81 vs 0.56, and 0.77 vs 0.40).
- **Interval widths.** DIF and CM intervals are about the right width (0.93–0.98 of the width implied by their spread on average, and 0.82–1.03 for CM pair by pair, across both simulations), so their undercoverage comes from bias. GPS intervals are about 30% too wide on average (1.20–1.41 pair by pair, across both simulations), which gives its 100% coverage.
- **Which biases are real.** CM and DIF are biased beyond Monte Carlo noise (the chance variation from running a finite number of replications) in nearly every pair (|z| > 3 in 34 of 36 pairs each in Simulation 2). GPS is beyond noise in 19 of 36 pairs in Simulation 2 and in 0 of 6 in Simulation 1, against τ\_pop.

## Against the paper

The qualitative pattern replicates; the GPS bias in Simulation 1 and some pair-level numbers do not.

**Replicated**

- **Pattern.** DIF is clearly biased, CM and GPS reduce the bias, CM is more efficient, and GPS coverage is 99–100% in both the paper's Table 1 and ours.
- **Simulation 2 bias.** Our bias for the 36 effects tracks the paper's Table 3 closely but with the opposite sign: across pairs, against τ\_pop, the correlation with the paper's published values is −0.95 (DIF), −0.76 (CM) and −0.87 (GPS). It fits best against τ\_pop (flipped correlations 0.95, 0.76, 0.87, against 0.94, 0.59, 0.47 for τ\_bal). The paper's sign convention is unknown, because its Table 1 reports only absolute bias.

**Not replicated**

- **GPS bias in Simulation 1.** The paper's is 0.06–0.48 (mean 0.29) and larger than CM's in all 6 pairs. Ours is 0.01–0.05 (mean 0.026) and smaller than CM's in every pair.
- **DIF in Simulation 1.** Similar size, but different pair by pair (pair (1,3): ours 1.16, paper 0.47).
- **CM coverage in Simulation 1.** Ours 83–95% against τ\_pop (79–95% against τ\_bal), the paper's 92–98%.

**Sensitivity check (D13).** The paper's Eq. (4) prints the same score in both terms of the effect. We read it as a typo in the main run, then tested the literal reading on Simulation 1. Applied only to units outside both groups (GPS\_lit), mean absolute bias is 0.044. Applied to units in group t as well (GPS\_lit2), it rises to 0.159, about half the way to the paper's 0.293; MSE rises 50–134% against the main GPS and coverage stays at 99–100% (both against τ\_pop). This would explain part of the gap, not all of it, and only if the paper's code followed Eq. (4) as printed: Yang et al.'s (2016) own package computes the main GPS, not the literal reading (D7).

## Decisions D1–D13

Thirteen points where the paper is silent, ambiguous or inconsistent, and what we chose. The full text is in docs/spec.md.

| ID | Open point | What we did |
| --- | --- | --- |
| D1 | The paper's denominator, 1 + sum of exp(η), makes the probabilities sum to less than 1 | Standard multinomial logit with η₁ = 0 and no extra 1 |
| D2 | Z1 \~ "N\[−3, 3\]" is ambiguous | Uniform(−3, 3); its mean of 0 fits the centred z2 and z3 terms, which make all Simulation 2 true effects 0 |
| D3 | Fixed group sizes versus random treatment assignment | Draw a large pool, assign by the multinomial probabilities, keep 100 (Sim 1) or 200 (Sim 2) per group |
| D4 | True effects are not stated | Two targets, both reported: population ATE (primary, the paper's Eq. 1) and balanced-design ATE |
| D5 | The GPS model is not stated | Linear multinomial logit in x1–x4 and z1–z3, misspecified on purpose |
| D6 | The paper builds the covariance matrix from X only but matches on X and Z | Mahalanobis distance on all 7 covariates, pooled sample covariance |
| D7 | Eq. (4) uses p\_t in both terms of the effect | Use p\_t′ for the t′ term, treated as a typo. This is how Yang et al.'s (2016) R package multilevelMatching computes the estimator; given our fitted scores, it reproduces our GPS estimates exactly on a test replication. The literal reading is tested in D13 |
| D8 | The confidence interval method is not stated | Percentile bootstrap, B = 200, resampling within treatment groups, refitting the GPS model and covariance matrix each time. A later sensitivity analysis keeps every estimate and replaces only the intervals with analytic Abadie–Imbens ones. With AI 2006 intervals, CM coverage moves closer to the paper's (94.8% against 91.4% in Simulation 1, 90.1% against 78.4% in Simulation 2; paper 95.0% and 97.1%), consistent with, but not proof of, the paper having used analytic variances. AI 2016 intervals, which allow for the estimated GPS, bring GPS coverage to about 95%. See crosscheck/side\_by\_side/README.md, section 5.2 |
| D9 | The D3 design does not sample from the population: it gives an equal-weight mixture of the group-specific covariate distributions. CM and GPS average over the sample's units, so they estimate the ATE for that mixture, not τ\_pop | Report bias, MSE and coverage against both τ\_pop and τ\_bal |
| D10 | η₆ and η₇ differ only by z3 versus z3², and z3 is binary | Kept as printed; p₆ = p₇ for every unit, so the two groups cannot be told apart by covariates |
| D11 | Thin tails in groups 1–5 and 9: poor overlap, some fitted probabilities near 0 | No trimming; the paper's simulation section describes none |
| D12 | Recording the software environment | Base R, parallel, nnet and MASS; no renv lockfile. R and base-package versions are in results/run\_session\_info.txt; for the full run, nnet 7.3.20 and MASS 7.3.65 appear only in results/run\_log.txt (later runs list them in the session file too) |
| D13 | Our GPS bias is smaller than the paper's | Sensitivity runs GPS\_lit and GPS\_lit2 (Simulation 1 only) that read Eq. (4) literally; the main GPS is unchanged. Neither variant is what multilevelMatching computes |

## What this does not show

These results support the pattern, not the paper's exact numbers, and several claims stay out of reach.

- **Which reading of Eq. (4) the paper used.** The literal reading explains part of the GPS gap, but it is not how Yang et al.'s own package computes the estimator, and we do not know whether the paper used that package or its own code.
- **Why some numbers differ.** Our Simulation 1 DIF differences do not depend on any estimator choice, so the data design probably differs somewhere (D1–D3). Simulation 2's DIF pattern matches closely, which argues against a large mismatch.
- **Which sign convention the paper used.** Its Table 1 reports only absolute bias.
- **Anything beyond this data-generating process.** As in the paper (its Sections 2 and 4), the two treatment levels are combined into one multi-group treatment, so the results say little about other multilevel designs.
- **The learning-mode application.** The paper does not describe how its Section 6 data were generated, so Table 5 was not replicated; its Fig. 2 was redrawn from the published table.
- **A pattern match in Simulation 1.** Six pairs are too few: the pattern correlation with the paper is 0.62 for the main GPS and 0.40 for GPS\_lit2.

## Acknowledgments

This replication was written from the simulation procedures described in the published paper (Guo, Liu & Wang, 2022). The authors' original simulation code was not available to us. The code was written with assistance from Claude Code.

## Conclusions

In these simulations, covariate matching and GPS matching both remove most of the bias in the raw group difference, as the paper reports. GPS had the smallest bias and CM the smallest mean squared error, so the better choice depends on whether bias or variance matters more. The published numbers cannot be reproduced exactly from the paper alone, because thirteen points are left open; reading Eq. (4) literally would close about half of the GPS discrepancy in Simulation 1, but that is not how Yang et al.'s own package computes the estimator, so it explains the gap only if the paper's code followed Eq. (4) as printed. Like the paper's, these simulations have no school clusters: the two treatments are combined into one, and matching ignores which school a student attends. How well the methods work on genuinely clustered data is a separate question that these results do not answer.
