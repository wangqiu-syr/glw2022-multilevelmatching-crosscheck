# Exploratory checks of the paper's Table 3

These checks compare the published Simulation 2 results (Table 3) with the rebuild, using only existing results; no simulations were run. They describe patterns in the published numbers and cannot establish their cause without the original code or output.

## Exchanging the CM and GPS labels

With the paper's sign reversed (see the note on sign in [two_way.md](two_way.md)), the published GPS bias matches our CM bias closely (correlation 0.946, mean absolute difference 0.097), and the published CM bias matches our GPS bias (mean absolute difference 0.119). Compared as printed, the correlations are 0.874 (GPS) and 0.764 (CM), with mean absolute differences of 0.232 and 0.202. After the exchange, the mean absolute biases also align (0.307 vs 0.295; 0.139 vs 0.099), and GPS becomes the more biased estimator in 4 of 36 pairs in the paper against 5 of 36 in the rebuild; as printed, GPS is the more biased estimator in 32 of 36 pairs in the paper but in 5 of 36 in the rebuild.

## The MSE column does not fit the exchange

The published MSE column matches the rebuild without any exchange: published and rebuilt GPS MSE correlate at 0.963, and GPS has the larger MSE in 36 of 36 pairs in the paper and 33 of 36 in the rebuild. Exchanging only the bias labels would leave three CM pairs with a squared bias larger than their MSE, which cannot occur within one set of runs. A simple exchange of two columns therefore does not explain Table 3; the bias and MSE columns would have to come from different runs or be computed differently.

## Raw-difference MSE

The published MSE for the raw difference (DIF) is about 100 times ours (median ratio 100.7), while the published CM and GPS MSEs are on the rebuild's scale. This suggests a difference in scaling or transcription in that column; its source cannot be determined here.

## Simulation 1

With six pairs, the evidence is weak. An exchange fits the mean levels better, but the correlations are mixed, and the literal reading of Eq. (4) (D7, tested in D13) already accounts for about half of the Simulation 1 GPS gap without an exchange.

## Reproduce

`Rscript R/table3_checks.R` (reads results/*.csv; writes results/table3_checks.csv).
