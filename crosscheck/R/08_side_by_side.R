# 08_side_by_side.R -- side-by-side tables and figures: the paper (Guo, Liu &
# Wang, 2022), the rebuild (main run) and Yang et al.'s package
# (multilevelMatching), plus the D8 sensitivity analysis (analytic intervals).
#
#   Rscript crosscheck/R/08_side_by_side.R      (from the repo root; base R only)
#
# Reads  results/paper_table{1,3}.csv, results/sim{1,2}_summary.csv,
#        results/true_effects_sim{1,2}.csv,
#        crosscheck/results/sim{1,2}_crosscheck_{reps,summary}.csv,
#        crosscheck/results/sim{1,2}_d8_gps_variances.csv.gz (from 07_d8_analytic_intervals.R).
# Writes crosscheck/side_by_side/tables/*.csv and *.md, crosscheck/side_by_side/figures/*.png.
#
# "Rebuild" values are the main run's (results/sim<k>_summary.csv). "Package"
# is multiMatch() as a user would run it: CM from match_on = "covariates", GPS
# from match_on = "multinom" (its own GPS fit). D8 intervals are centred on the
# rebuild's own estimates: CM +/- 1.96 SE with the AI 2006 variance (the
# "covariates" fit reproduces our CM exactly); GPS +/- 1.96 SE with the AI 2006
# and AI 2016 variances computed for our fitted GPS (07_d8_analytic_intervals.R).

OUT <- file.path("crosscheck", "side_by_side")
TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
COL <- c(CM = "green3", GPS = "blue")          # the paper's colours, as in R/plots.R
z <- qnorm(0.975)
key <- function(t, tp) paste0("(", t, ",", tp, ")")

md_table <- function(df, digits = 3) {
  # digits: one number, or a named vector by column; whole-number columns print as integers
  dg <- function(j) if (!is.null(names(digits)) && j %in% names(digits)) digits[[j]] else if (is.null(names(digits))) digits else 3
  f <- function(v, j) {
    if (!is.numeric(v)) return(as.character(v))
    w <- v[!is.na(v)]; d <- if (length(w) && all(w == round(w))) 0 else dg(j)
    ifelse(is.na(v), "", formatC(v, format = "f", digits = d))
  }
  cells <- vapply(names(df), function(j) f(df[[j]], j), character(nrow(df)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = 1)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(cells, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
save_table <- function(df, name, digits = 3) {
  write.csv(df, file.path(TAB, paste0(name, ".csv")), row.names = FALSE)
  writeLines(md_table(df, digits), file.path(TAB, paste0(name, ".md")))
}

# ---- assemble per-pair results ------------------------------------------------
per_pair <- list(); d8_long <- list()
for (sim in 1:2) {
  te   <- read.csv(file.path("results", sprintf("true_effects_sim%d.csv", sim)))
  main <- read.csv(file.path("results", sprintf("sim%d_summary.csv", sim)))
  xs   <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_crosscheck_summary.csv", sim)))
  reps <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_crosscheck_reps.csv.gz", sim)))
  d8   <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_d8_gps_variances.csv.gz", sim)))
  pairs <- te[c("t", "tp")]; pairs$pair <- key(pairs$t, pairs$tp)

  if (sim == 1) {
    p <- read.csv(file.path("results", "paper_table1.csv"))
    paper <- data.frame(pair = paste0("(", p$pair, ")"), CM_bias = NA, CM_abs = p$CM_abs_bias, CM_mse = NA,
                        CM_cov = p$CM_coverage, GPS_bias = NA, GPS_abs = p$GPS_abs_bias, GPS_mse = NA, GPS_cov = p$GPS_coverage)
  } else {
    p <- read.csv(file.path("results", "paper_table3.csv"))
    paper <- data.frame(pair = paste0("(", p$pair, ")"), CM_bias = p$CM_bias, CM_abs = abs(p$CM_bias), CM_mse = p$CM_mse,
                        CM_cov = p$CM_coverage, GPS_bias = p$GPS_bias, GPS_abs = abs(p$GPS_bias), GPS_mse = p$GPS_mse, GPS_cov = p$GPS_coverage)
  }
  paper <- paper[match(pairs$pair, paper$pair), ]

  m_get <- function(e, col, tg = "pop") { x <- main[main$estimator == e & main$target == tg, ]; x[[col]][match(pairs$pair, key(x$t, x$tp))] }
  x_get <- function(e, col, tg = "pop") { x <- xs[xs$estimator == e & xs$target == tg, ];   x[[col]][match(pairs$pair, key(x$t, x$tp))] }

  # ---- D8: analytic intervals centred on the rebuild's estimates ----
  stopifnot(identical(signif(d8$GPS, 10), signif(reps$GPS[match(paste(d8$rep, d8$t, d8$tp), paste(reps$rep, reps$t, reps$tp))], 10)))
  r8 <- merge(reps[c("rep", "t", "tp", "CM", "GPS", "var_pkg_CM")], d8[c("rep", "t", "tp", "var_AI2006", "var_AI2016")])
  for (tg in c("pop", "bal")) {
    tau <- te[[paste0("tau_", tg)]][match(key(r8$t, r8$tp), pairs$pair)]
    for (spec in list(c("CM", "CM", "var_pkg_CM", "AI 2006"), c("GPS", "GPS", "var_AI2006", "AI 2006"),
                      c("GPS", "GPS", "var_AI2016", "AI 2016"))) {
      v <- r8[[spec[3]]]; se <- sqrt(ifelse(v < 0, NA, v)); est <- r8[[spec[2]]]
      cover <- abs(est - tau) <= z * se
      a <- aggregate(data.frame(coverage = 100 * cover, width = 2 * z * se, neg = v < 0),
                     by = list(t = r8$t, tp = r8$tp), FUN = function(x) mean(x, na.rm = TRUE))
      a$n_unavailable <- aggregate(v < 0, by = list(r8$t, r8$tp), FUN = sum)$x
      d8_long[[length(d8_long) + 1]] <- data.frame(sim = sim, estimator = spec[1], method = spec[4], target = tg,
        t = a$t, tp = a$tp, coverage = a$coverage, mean_width = a$width, n_unavailable = a$n_unavailable)
    }
    for (e in c("CM", "GPS")) {
      x <- main[main$estimator == e & main$target == tg, ]
      d8_long[[length(d8_long) + 1]] <- data.frame(sim = sim, estimator = e, method = "bootstrap (main run)", target = tg,
        t = x$t, tp = x$tp, coverage = x$coverage, mean_width = x$mean_ci_width, n_unavailable = 0)
    }
  }
  d8s <- do.call(rbind, d8_long[vapply(d8_long, function(x) x$sim[1] == sim, logical(1))])
  dget <- function(e, m, col = "coverage", tg = "pop") { x <- d8s[d8s$estimator == e & d8s$method == m & d8s$target == tg, ]; x[[col]][match(pairs$pair, key(x$t, x$tp))] }

  per_pair[[sim]] <- data.frame(
    sim = sim, pair = pairs$pair, tau_pop = te$tau_pop,
    CM_bias_paper = paper$CM_bias, CM_bias_rebuild = m_get("CM", "bias"), CM_bias_package = x_get("pkg_CM", "bias"),
    CM_absbias_paper = paper$CM_abs, CM_absbias_rebuild = m_get("CM", "abs_bias"), CM_absbias_package = x_get("pkg_CM", "abs_bias"),
    CM_mse_paper = paper$CM_mse, CM_mse_rebuild = m_get("CM", "mse"), CM_mse_package = x_get("pkg_CM", "mse"),
    CM_cov_paper = paper$CM_cov, CM_cov_bootstrap = m_get("CM", "coverage"), CM_cov_AI2006 = dget("CM", "AI 2006"),
    GPS_bias_paper = paper$GPS_bias, GPS_bias_rebuild = m_get("GPS", "bias"), GPS_bias_package = x_get("pkg_GPS_multinom", "bias"),
    GPS_absbias_paper = paper$GPS_abs, GPS_absbias_rebuild = m_get("GPS", "abs_bias"), GPS_absbias_package = x_get("pkg_GPS_multinom", "abs_bias"),
    GPS_absbias_package_existing = x_get("pkg_GPS_existing", "abs_bias"),
    GPS_mse_paper = paper$GPS_mse, GPS_mse_rebuild = m_get("GPS", "mse"), GPS_mse_package = x_get("pkg_GPS_multinom", "mse"),
    GPS_cov_paper = paper$GPS_cov, GPS_cov_bootstrap = m_get("GPS", "coverage"),
    GPS_cov_AI2006 = dget("GPS", "AI 2006"), GPS_cov_AI2016 = dget("GPS", "AI 2016"))
}
d8_long <- do.call(rbind, d8_long)

# ---- tables -------------------------------------------------------------------
for (sim in 1:2) {
  pp <- per_pair[[sim]]
  write.csv(pp, file.path(TAB, sprintf("sim%d_side_by_side.csv", sim)), row.names = FALSE)
  mean_row <- function(df) { m <- df[1, ]; m[] <- NA; m$Pair <- "Mean"; for (j in names(df)[-1]) m[[j]] <- mean(df[[j]]); rbind(df, m) }
  bias <- data.frame(Pair = pp$pair,
    `CM paper` = pp$CM_absbias_paper, `CM rebuild` = pp$CM_absbias_rebuild, `CM package` = pp$CM_absbias_package,
    `GPS paper` = pp$GPS_absbias_paper, `GPS rebuild` = pp$GPS_absbias_rebuild, `GPS package` = pp$GPS_absbias_package, check.names = FALSE)
  save_table(mean_row(bias), sprintf("sim%d_abs_bias", sim))
  mse <- data.frame(Pair = pp$pair, `CM paper` = pp$CM_mse_paper, `CM rebuild` = pp$CM_mse_rebuild, `CM package` = pp$CM_mse_package,
    `GPS paper` = pp$GPS_mse_paper, `GPS rebuild` = pp$GPS_mse_rebuild, `GPS package` = pp$GPS_mse_package, check.names = FALSE)
  if (sim == 1) mse <- mse[, !grepl("paper", names(mse))]
  save_table(mean_row(mse), sprintf("sim%d_mse", sim))
  cov <- data.frame(Pair = pp$pair, `CM paper` = pp$CM_cov_paper, `CM bootstrap` = pp$CM_cov_bootstrap, `CM AI 2006` = pp$CM_cov_AI2006,
    `GPS paper` = pp$GPS_cov_paper, `GPS bootstrap` = pp$GPS_cov_bootstrap, `GPS AI 2006` = pp$GPS_cov_AI2006,
    `GPS AI 2016` = pp$GPS_cov_AI2016, check.names = FALSE)
  save_table(mean_row(cov), sprintf("sim%d_coverage", sim), digits = 1)
}

# Orderings: which estimator is less biased / more efficient, per pair.
ord <- do.call(rbind, lapply(1:2, function(sim) {
  pp <- per_pair[[sim]]; n <- nrow(pp)
  rbind(
    data.frame(Simulation = sim, Comparison = "GPS less biased than CM (pairs)",
      Paper = sum(pp$GPS_absbias_paper < pp$CM_absbias_paper), Rebuild = sum(pp$GPS_absbias_rebuild < pp$CM_absbias_rebuild),
      Package = sum(pp$GPS_absbias_package < pp$CM_absbias_package), Pairs = n),
    data.frame(Simulation = sim, Comparison = "CM lower MSE than GPS (pairs)",
      Paper = if (sim == 1) NA else sum(pp$CM_mse_paper < pp$GPS_mse_paper), Rebuild = sum(pp$CM_mse_rebuild < pp$GPS_mse_rebuild),
      Package = sum(pp$CM_mse_package < pp$GPS_mse_package), Pairs = n),
    data.frame(Simulation = sim, Comparison = "Mean absolute bias, CM",
      Paper = mean(pp$CM_absbias_paper), Rebuild = mean(pp$CM_absbias_rebuild), Package = mean(pp$CM_absbias_package), Pairs = n),
    data.frame(Simulation = sim, Comparison = "Mean absolute bias, GPS",
      Paper = mean(pp$GPS_absbias_paper), Rebuild = mean(pp$GPS_absbias_rebuild), Package = mean(pp$GPS_absbias_package), Pairs = n))
}))
save_table(ord[grepl("pairs", ord$Comparison), ], "orderings_counts")
mb <- ord[grepl("Mean absolute", ord$Comparison), ]; names(mb)[2] <- "Mean absolute bias vs τ_pop"
mb[[2]] <- sub("Mean absolute bias, ", "", mb[[2]])
save_table(mb[c("Simulation", "Mean absolute bias vs τ_pop", "Paper", "Rebuild", "Package")], "mean_abs_bias")
write.csv(ord, file.path(TAB, "orderings.csv"), row.names = FALSE)

# Headline: mean absolute bias and mean MSE over pairs, paper vs rebuild vs package.
head_tab <- do.call(rbind, lapply(1:2, function(sim) do.call(rbind, lapply(c("CM", "GPS"), function(e) {
  pp <- per_pair[[sim]]; g <- function(x) mean(pp[[paste0(e, "_", x)]])
  data.frame(Simulation = sim, Estimator = e,
             `Abs. bias: paper` = g("absbias_paper"), `Abs. bias: rebuild` = g("absbias_rebuild"), `Abs. bias: package` = g("absbias_package"),
             `MSE: paper` = if (sim == 1) NA else g("mse_paper"), `MSE: rebuild` = g("mse_rebuild"), `MSE: package` = g("mse_package"),
             check.names = FALSE) }))))
save_table(head_tab, "bias_mse_summary")

# D8 summary: coverage by interval method, and distance from the paper's coverage.
d8_sum <- do.call(rbind, lapply(split(d8_long, list(d8_long$sim, d8_long$estimator, d8_long$method, d8_long$target), drop = TRUE), function(g) {
  pp <- per_pair[[g$sim[1]]]
  pc <- pp[[paste0(g$estimator[1], "_cov_paper")]][match(key(g$t, g$tp), pp$pair)]
  data.frame(Simulation = g$sim[1], Estimator = g$estimator[1], Method = g$method[1], Target = g$target[1],
             `Mean coverage` = mean(g$coverage), `Min` = min(g$coverage), `Max` = max(g$coverage),
             `Mean width` = mean(g$mean_width),
             `Mean abs. diff. from paper` = if (g$target[1] == "pop") mean(abs(g$coverage - pc)) else NA,
             `Unavailable intervals` = sum(g$n_unavailable), check.names = FALSE)
}))
d8_sum <- d8_sum[order(d8_sum$Target != "pop", d8_sum$Simulation, d8_sum$Estimator,
                       match(d8_sum$Method, c("bootstrap (main run)", "AI 2006", "AI 2016"))), ]
paper_rows <- do.call(rbind, lapply(1:2, function(sim) do.call(rbind, lapply(c("CM", "GPS"), function(e) {
  pc <- per_pair[[sim]][[paste0(e, "_cov_paper")]]
  data.frame(Simulation = sim, Estimator = e, Method = "paper", Target = "pop", `Mean coverage` = mean(pc),
             Min = min(pc), Max = max(pc), `Mean width` = NA, `Mean abs. diff. from paper` = NA,
             `Unavailable intervals` = NA, check.names = FALSE) }))))
d8_sum <- rbind(paper_rows, d8_sum)
d8_sum <- d8_sum[order(d8_sum$Target != "pop", d8_sum$Simulation, d8_sum$Estimator,
                       match(d8_sum$Method, c("paper", "bootstrap (main run)", "AI 2006", "AI 2016"))), ]
d8dig <- c(`Mean coverage` = 1, Min = 1, Max = 1, `Mean width` = 2, `Mean abs. diff. from paper` = 1)
write.csv(d8_sum, file.path(TAB, "d8_coverage_summary.csv"), row.names = FALSE)
writeLines(md_table(d8_sum[d8_sum$Target == "pop", names(d8_sum) != "Target"], d8dig), file.path(TAB, "d8_coverage_summary_pop.md"))
writeLines(md_table(d8_sum[d8_sum$Target == "bal", !names(d8_sum) %in% c("Target", "Mean abs. diff. from paper")], d8dig),
           file.path(TAB, "d8_coverage_summary_bal.md"))
write.csv(d8_long, file.path(TAB, "d8_coverage_by_pair.csv"), row.names = FALSE)

# Interval width against the actual spread of the estimates: implied standard
# error (mean width / 3.92; for the percentile bootstrap an approximation)
# versus the Monte Carlo SD of the estimates, averaged over pairs (tau_pop).
spread <- do.call(rbind, lapply(1:2, function(sim) do.call(rbind, lapply(c("CM", "GPS"), function(e) {
  m <- read.csv(file.path("results", sprintf("sim%d_summary.csv", sim))); m <- m[m$estimator == e & m$target == "pop", ]
  w <- function(meth) { x <- d8_long[d8_long$sim == sim & d8_long$estimator == e & d8_long$method == meth & d8_long$target == "pop", ]
                        if (nrow(x)) mean(x$mean_width) / (2 * z) else NA }
  data.frame(Simulation = sim, Estimator = e, `Monte Carlo SD of the estimates` = mean(m$sd_est),
             `Bootstrap` = w("bootstrap (main run)"), `AI 2006` = w("AI 2006"), `AI 2016` = w("AI 2016"),
             `Mean abs. bias / SD` = mean(m$abs_bias / m$sd_est), check.names = FALSE) }))))
save_table(spread, "d8_interval_width_vs_spread")

# Is low bootstrap coverage explained by bias? For each pair, a normal
# approximation predicts coverage from the estimator's bias b, its Monte Carlo
# SD s and the interval's half-width h = mean width / 2:
#   P(|b + s Z| <= h),  Z ~ N(0, 1);  "if unbiased" sets b = 0.
# Main-run bootstrap intervals, tau_pop.
cov_bias <- do.call(rbind, lapply(1:2, function(sim) do.call(rbind, lapply(c("CM", "GPS"), function(e) {
  m <- read.csv(file.path("results", sprintf("sim%d_summary.csv", sim))); x <- m[m$estimator == e & m$target == "pop", ]
  h <- x$mean_ci_width / 2; b <- x$bias; sg <- x$sd_est
  pred <- 100 * (pnorm((h - b) / sg) - pnorm((-h - b) / sg)); nob <- 100 * (pnorm(h / sg) - pnorm(-h / sg))
  i <- which.min(x$coverage)
  data.frame(Simulation = sim, Estimator = e, `Observed mean coverage` = mean(x$coverage),
             `Predicted from bias, SD and width` = mean(pred), `Predicted if unbiased` = mean(nob),
             `Correlation, observed vs predicted (pairs)` = if (sd(x$coverage) > 0) cor(pred, x$coverage) else NA,
             `Lowest pair` = key(x$t[i], x$tp[i]), `Its coverage` = x$coverage[i], `Its abs. bias / SD` = abs(b[i]) / sg[i],
             check.names = FALSE) }))))
save_table(cov_bias, "coverage_vs_bias", digits = 2)

# ---- figures --------------------------------------------------------------------
pch_src <- c(paper = 1, rebuild = 16, package = 17)
fig_bias <- function(sim, file, w, h) {
  pp <- per_pair[[sim]]; k <- nrow(pp)
  png(file.path(FIG, file), width = w, height = h, res = 150)
  op <- par(mfrow = if (sim == 1) c(1, 2) else c(2, 1), mar = c(if (sim == 1) 4.5 else 5.5, 4.5, 2.5, 1), oma = c(0, 0, 2.2, 0))
  for (e in c("CM", "GPS")) {
    v <- cbind(paper = pp[[paste0(e, "_absbias_paper")]], rebuild = pp[[paste0(e, "_absbias_rebuild")]],
               package = pp[[paste0(e, "_absbias_package")]])
    off <- c(paper = 0, rebuild = -0.18, package = 0.18) * (if (sim == 1) 1 else 0)
    plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0, max(v) * (if (sim == 1) 1.55 else 1.3)), xaxt = "n",
         xlab = if (sim == 1) "Pair (t, t')" else "", ylab = expression("Absolute bias vs " * tau[pop]), main = e, las = 1)
    abline(h = axTicks(2), col = "grey92"); abline(h = 0, col = "grey60")
    axis(1, at = 1:k, labels = pp$pair, las = if (sim == 1) 1 else 2, cex.axis = if (sim == 1) 0.9 else 0.7)
    if (sim == 2) {
      lines(1:k, v[, "paper"], lty = 3, lwd = 1.3, col = COL[[e]])
      lines(1:k, v[, "rebuild"], lty = 1, lwd = 1.8, col = COL[[e]])
      lines(1:k, v[, "package"], lty = 2, lwd = 1.3, col = "grey20")
    }
    for (s in names(pch_src)) points((1:k) + off[[s]], v[, s], pch = pch_src[[s]],
                                     col = if (s == "package") "grey20" else COL[[e]], cex = if (sim == 1) 1.4 else 0.9, lwd = 1.5)
    legend(if (sim == 1) "topleft" else "top", horiz = sim == 2, bty = "n", cex = 0.8, pch = pch_src, pt.lwd = 1.5,
           lty = if (sim == 2) c(3, 1, 2) else c(0, 0, 0), col = c(COL[[e]], COL[[e]], "grey20"),
           legend = c(sprintf("paper (Table %d)", if (sim == 1) 1 else 3), "rebuild (main run)", "package (multilevelMatching)"))
  }
  mtext(sprintf("Simulation %d: absolute bias by pair, paper vs rebuild vs Yang et al.'s package", sim), outer = TRUE, cex = 1.05, font = 2)
  par(op); dev.off()
}
fig_bias(1, "fig1_sim1_abs_bias.png", 1600, 750)
fig_bias(2, "fig2_sim2_abs_bias.png", 1600, 1250)

# MSE by pair: Simulation 1 (rebuild, package; the paper's Table 1 has no MSE)
# and Simulation 2 (the paper's Table 3, rebuild, package).
png(file.path(FIG, "fig3_mse_by_pair.png"), width = 1600, height = 1300, res = 150)
op <- par(mfrow = c(2, 2), mar = c(5.5, 4.5, 2.5, 1), oma = c(0, 0, 2.2, 0))
for (sim in 1:2) for (e in c("CM", "GPS")) {
  pp <- per_pair[[sim]]; k <- nrow(pp)
  src <- if (sim == 1) c("rebuild", "package") else c("paper", "rebuild", "package")
  v <- sapply(src, function(s) pp[[paste0(e, "_mse_", s)]])
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0, max(v) * 1.3), xaxt = "n", xlab = "", ylab = "MSE",
       main = sprintf("Simulation %d, %s", sim, e), las = 1)
  abline(h = axTicks(2), col = "grey92"); abline(h = 0, col = "grey60")
  axis(1, at = 1:k, labels = pp$pair, las = 2, cex.axis = if (sim == 1) 0.9 else 0.6)
  off <- if (sim == 1) c(rebuild = -0.15, package = 0.15) else c(paper = 0, rebuild = 0, package = 0)
  if (sim == 2) {
    lines(1:k, v[, "paper"], lty = 3, lwd = 1.3, col = COL[[e]])
    lines(1:k, v[, "rebuild"], lty = 1, lwd = 1.8, col = COL[[e]])
    lines(1:k, v[, "package"], lty = 2, lwd = 1.3, col = "grey20")
  }
  for (s in src) points((1:k) + off[[s]], v[, s], pch = pch_src[[s]], col = if (s == "package") "grey20" else COL[[e]],
                        cex = if (sim == 1) 1.3 else 0.8, lwd = 1.5)
  legend("top", horiz = TRUE, bty = "n", cex = 0.75, pch = pch_src[src], pt.lwd = 1.5,
         lty = if (sim == 2) c(3, 1, 2) else c(0, 0), col = ifelse(src == "package", "grey20", COL[[e]]),
         legend = c(paper = "paper (Table 3)", rebuild = "rebuild (main run)", package = "package")[src])
}
mtext(expression(bold("Mean squared error by pair, against ") * bold(tau[pop]) * bold(": paper vs rebuild vs Yang et al.'s package")),
      outer = TRUE, cex = 1.05)
par(op); dev.off()

# Coverage by pair and interval method.
png(file.path(FIG, "fig5_coverage_by_method.png"), width = 1600, height = 1300, res = 150)
op <- par(mfrow = c(2, 2), mar = c(5.5, 4.5, 2.5, 1), oma = c(3, 0, 2.2, 0))
meth <- c("paper", "bootstrap", "AI2006", "AI2016"); pch_m <- c(paper = 1, bootstrap = 16, AI2006 = 2, AI2016 = 18)
lab_m <- c(paper = "paper", bootstrap = "rebuild: bootstrap (main run)", AI2006 = "rebuild: AI 2006", AI2016 = "rebuild: AI 2016 (GPS only)")
for (sim in 1:2) for (e in c("CM", "GPS")) {
  pp <- per_pair[[sim]]; k <- nrow(pp)
  ms <- if (e == "CM") meth[1:3] else meth
  v <- sapply(ms, function(m) pp[[paste0(e, "_cov_", m)]])
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(min(v, 90, na.rm = TRUE) - 2, 101), xaxt = "n", xlab = "",
       ylab = expression("Coverage of " * tau[pop] * " (%)"),
       main = sprintf("Simulation %d, %s", sim, e), las = 1)
  abline(h = axTicks(2), col = "grey92"); abline(h = 95, lty = 2, col = "grey40")
  axis(1, at = 1:k, labels = pp$pair, las = 2, cex.axis = if (sim == 1) 0.9 else 0.6)
  for (m in ms) points(1:k, pp[[paste0(e, "_cov_", m)]], pch = pch_m[[m]], col = if (m == "paper") "black" else COL[[e]],
                       cex = if (sim == 1) 1.3 else 0.85, lwd = 1.3)
}
mtext("Coverage of 95% intervals by pair: the paper, and the rebuild's estimates with bootstrap or analytic (Abadie-Imbens) intervals (D8)",
      outer = TRUE, cex = 1, font = 2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
legend("bottom", horiz = TRUE, bty = "n", cex = 0.85, pch = c(pch_m, NA), pt.lwd = 1.3, col = c("black", rep("grey25", 3), NA),
       legend = c(lab_m, "colour: green = CM, blue = GPS"), inset = 0.005)
par(op); dev.off()

# Per-replication GPS estimates: package's own fit vs the rebuild (D5).
png(file.path(FIG, "fig4_gps_per_replication.png"), width = 1600, height = 800, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 2.2, 0))
for (sim in 1:2) {
  reps <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_crosscheck_reps.csv.gz", sim)))
  lim <- range(reps$GPS, reps$pkg_GPS_multinom)
  plot(reps$GPS, reps$pkg_GPS_multinom, pch = 16, cex = 0.35, col = adjustcolor(COL[["GPS"]], 0.25), xlim = lim, ylim = lim,
       xlab = "Rebuild GPS estimate", ylab = "Package GPS estimate (own fit, nnet defaults)", las = 1,
       main = sprintf("Simulation %d: %s of %s identical", sim,
                      format(sum(abs(reps$GPS - reps$pkg_GPS_multinom) < 1e-12), big.mark = ","), format(nrow(reps), big.mark = ",")))
  abline(0, 1, col = "grey30", lwd = 1)
}
mtext("Every GPS estimate (all pairs, all replications): the package's own GPS fit against the rebuild's (D5)", outer = TRUE, cex = 1, font = 2)
par(op); dev.off()

cat("Wrote tables:\n"); print(list.files(TAB)); cat("Wrote figures:\n"); print(list.files(FIG))
