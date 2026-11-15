# table3_checks.R -- exploratory checks of the paper's Table 3 (and Table 1)
# against the rebuild (docs/table3_checks.md). Uses existing results only;
# no simulations are run. Run from the repository root:
#   Rscript R/table3_checks.R
# Reads results/sim1_vs_paper.csv, results/sim2_vs_paper.csv and
# results/gps_literal2_sim1.csv. Writes results/table3_checks.csv.

s1  <- read.csv("results/sim1_vs_paper.csv")
s2  <- read.csv("results/sim2_vs_paper.csv")
lit <- read.csv("results/gps_literal2_sim1.csv")
g <- function(d, e, col) d[d$estimator == e, col]

out <- list()
add <- function(section, simulation, comparison, statistic, value) {
  out[[length(out) + 1]] <<- data.frame(section, simulation, comparison,
                                        statistic, value)
}

# Simulation 2 signed bias; the paper's sign is reversed (see docs/two_way.md).
pb <- list(CM = -g(s2, "CM", "paper_bias"), GPS = -g(s2, "GPS", "paper_bias"))
ob <- list(CM = g(s2, "CM", "ours_bias_pop"), GPS = g(s2, "GPS", "ours_bias_pop"))
for (pe in c("CM", "GPS")) for (oe in c("CM", "GPS")) {
  lab <- if (pe == oe) "labels as printed" else "labels exchanged"
  cmp <- sprintf("paper %s vs rebuild %s", pe, oe)
  add(lab, 2, cmp, "bias correlation", cor(pb[[pe]], ob[[oe]]))
  add(lab, 2, cmp, "bias mean absolute difference", mean(abs(pb[[pe]] - ob[[oe]])))
}
for (e in c("CM", "GPS")) {
  add("mean absolute bias", 2, paste("paper", e), "mean absolute bias", mean(abs(pb[[e]])))
  add("mean absolute bias", 2, paste("rebuild", e), "mean absolute bias", mean(abs(ob[[e]])))
}
paper_gps_more <- abs(pb$GPS) > abs(pb$CM)
ours_gps_more  <- abs(ob$GPS) > abs(ob$CM)
add("ordering", 2, "paper, labels as printed", "pairs with GPS more biased", sum(paper_gps_more))
add("ordering", 2, "rebuild", "pairs with GPS more biased", sum(ours_gps_more))
add("ordering", 2, "paper, labels exchanged", "pairs with GPS more biased", sum(!paper_gps_more))
add("ordering", 2, "paper exchanged vs rebuild", "pairs with the same ordering", sum(!paper_gps_more == ours_gps_more))
add("ordering", 2, "all", "pairs", length(ours_gps_more))

# Simulation 2 MSE.
pm <- list(CM = g(s2, "CM", "paper_mse"), GPS = g(s2, "GPS", "paper_mse"))
om <- list(CM = g(s2, "CM", "ours_mse_pop"), GPS = g(s2, "GPS", "ours_mse_pop"))
for (pe in c("CM", "GPS")) for (oe in c("CM", "GPS"))
  add("MSE", 2, sprintf("paper %s vs rebuild %s", pe, oe), "MSE correlation", cor(pm[[pe]], om[[oe]]))
add("MSE", 2, "paper", "pairs with GPS larger MSE", sum(pm$GPS > pm$CM))
add("MSE", 2, "rebuild", "pairs with GPS larger MSE", sum(om$GPS > om$CM))
add("MSE", 2, "paper CM MSE minus paper GPS bias squared",
    "pairs with negative implied variance", sum(pm$CM - pb$GPS^2 < 0))
add("MSE", 2, "paper GPS MSE minus paper CM bias squared",
    "pairs with negative implied variance", sum(pm$GPS - pb$CM^2 < 0))
add("MSE", 2, "paper, labels as printed",
    "pairs with negative implied variance", sum(pm$CM - pb$CM^2 < 0) + sum(pm$GPS - pb$GPS^2 < 0))

# Raw-difference (DIF) MSE scale.
for (e in c("DIF", "CM", "GPS"))
  add("MSE scale", 2, paste(e, "paper / rebuild"), "median MSE ratio",
      median(g(s2, e, "paper_mse") / g(s2, e, "ours_mse_pop")))

# Simulation 1 absolute bias (Table 1).
for (pe in c("CM", "GPS")) for (oe in c("CM", "GPS")) {
  p <- g(s1, pe, "paper_abs_bias"); o <- g(s1, oe, "ours_abs_bias_pop")
  lab <- if (pe == oe) "labels as printed" else "labels exchanged"
  add(lab, 1, sprintf("paper %s vs rebuild %s", pe, oe), "absolute bias correlation", cor(p, o))
}
for (e in c("CM", "GPS")) {
  add("mean absolute bias", 1, paste("paper", e), "mean absolute bias", mean(g(s1, e, "paper_abs_bias")))
  add("mean absolute bias", 1, paste("rebuild", e), "mean absolute bias", mean(g(s1, e, "ours_abs_bias_pop")))
}
gps <- mean(lit$gps_abs_bias_pop); lit2 <- mean(lit$gpslit2_abs_bias_pop)
paper <- mean(lit$paper_gps_abs_bias)
add("Eq. (4) literal reading", 1, "GPS_lit2", "mean absolute bias", lit2)
add("Eq. (4) literal reading", 1, "GPS_lit2", "share of GPS gap to paper closed",
    (lit2 - gps) / (paper - gps))

res <- do.call(rbind, out)
write.csv(res, "results/table3_checks.csv", row.names = FALSE)
res$value <- signif(res$value, 4)
print(res, right = FALSE)
