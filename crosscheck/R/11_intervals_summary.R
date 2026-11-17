# 11_intervals_summary.R -- tables and figure for crosscheck/intervals/README.md:
# why CM coverage is low, pair by pair, and where the percentile bootstrap
# interval is centred.
#
# Inputs (all in the repo): results/sim<k>_summary.csv and true_effects_sim<k>.csv
# (main run); crosscheck/results/sim<k>_cm_gap_reps.csv.gz and sim<k>_cm_bc_reps.csv.gz
# (09); crosscheck/results/sim<k>_crosscheck_reps.csv.gz (04: the package's AI 2006
# variance for CM, whose estimates equal ours); crosscheck/side_by_side/tables/
# d8_coverage_by_pair.csv (08); and, when present, crosscheck/results/
# sim<k>_cm_boot_reps.csv.gz (10).
#
#   Rscript crosscheck/R/11_intervals_summary.R
# Writes crosscheck/intervals/tables/*.{csv,md} and crosscheck/intervals/figures/*.png.

source("R/dgp.R")
OUT <- file.path("crosscheck", "intervals")
TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
COL <- c(CM = "green3", GPS = "blue")          # the paper's colours, as in R/plots.R
z <- qnorm(0.975)
key <- function(t, tp) paste0("(", t, ",", tp, ")")

md_table <- function(df, digits = 3) {        # as in 08_side_by_side.R
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
# Normal approximation: coverage (%) of est +- h when est ~ N(tau + b, s^2).
pcov <- function(b, s, h) 100 * (pnorm((h - b) / s) - pnorm((-h - b) / s))

d8 <- read.csv(file.path("crosscheck", "side_by_side", "tables", "d8_coverage_by_pair.csv"))
per_pair <- list(); summ_rows <- list(); boot_rows <- list(); var_rows <- list(); shift_rows <- list()

for (sim in 1:2) {
  setup <- sim_setup(sim); A <- setup$A; B <- setup$beta
  pr <- t(combn(A, 2)); pk <- key(pr[, 1], pr[, 2])
  te <- read.csv(file.path("results", sprintf("true_effects_sim%d.csv", sim)))
  te <- te[match(pk, key(te$t, te$tp)), ]
  su <- read.csv(file.path("results", sprintf("sim%d_summary.csv", sim)))
  ms <- function(e, col) { x <- su[su$estimator == e & su$target == "pop", ]; x[[col]][match(pk, key(x$t, x$tp))] }
  dd <- function(e, m, col = "coverage") { x <- d8[d8$sim == sim & d8$estimator == e & d8$method == m & d8$target == "pop", ]; x[[col]][match(pk, key(x$t, x$tp))] }

  # 09: mean matching gaps, noise and sample covariate mean over replications
  g <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_cm_gap_reps.csv.gz", sim)), check.names = FALSE)
  gc <- paste0("gap_", COV_NAMES); mc <- paste0("m_", COV_NAMES)
  Dbar <- t(vapply(1:A, function(s) colMeans(g[g$s == s, gc]), numeric(7))); dimnames(Dbar) <- list(paste0("t", 1:A), COV_NAMES)
  Nbar <- vapply(1:A, function(s) mean(g$N[g$s == s]), numeric(1))
  mbar <- colMeans(g[g$s == 1, mc])
  contrib <- B * Dbar                          # beta_s,k * Dbar_s,k
  gap_lvl <- rowSums(contrib)                  # beta_s' Dbar_s

  bias <- ms("CM", "bias"); sdv <- ms("CM", "sd_est"); h_boot <- ms("CM", "mean_ci_width") / 2
  mix <- as.vector((B[pr[, 1], ] - B[pr[, 2], ]) %*% mbar) - te$tau_pop
  gap_t <- gap_lvl[pr[, 1]]; gap_tp <- -gap_lvl[pr[, 2]]
  noise <- Nbar[pr[, 1]] - Nbar[pr[, 2]]
  stopifnot(max(abs(mix + gap_t + gap_tp + noise - bias)) < 1e-10)

  # Regression-adjusted CM (diagnostic). Its own interval would need its own
  # variance; here only an oracle-width check: |est - tau_pop| <= 1.96 * SD of
  # its estimates over replications.
  bc <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_cm_bc_reps.csv.gz", sim)))
  bc$tau <- te$tau_pop[match(key(bc$t, bc$tp), pk)]
  agg <- function(x, f) { a <- tapply(x, key(bc$t, bc$tp), f); unname(a[pk]) }
  bc_bias <- agg(bc$CM_regadj - bc$tau, mean); bc_sd <- agg(bc$CM_regadj, sd)
  bc$sd <- bc_sd[match(key(bc$t, bc$tp), pk)]
  bc_cov <- 100 * agg(abs(bc$CM_regadj - bc$tau) <= z * bc$sd, mean)
  h_ai <- dd("CM", "AI 2006", "mean_width") / 2

  # 09: the three terms of each estimate, and the AI 2006 variance parts
  tm <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_cm_terms_reps.csv.gz", sim)))
  ta <- function(x, f = mean) { a <- tapply(x, key(tm$t, tm$tp), f); unname(a[pk]) }
  vt <- data.frame(
    Pair = pk,
    `Var. of estimates` = ta(tm$est, var),
    `from sample mix` = ta(tm$mix, var),
    `from matching gap` = ta(tm$gap_term, var),
    `from noise` = ta(tm$noise_term, var),
    `covariances` = ta(tm$est, var) - ta(tm$mix, var) - ta(tm$gap_term, var) - ta(tm$noise_term, var),
    `AI 2006: V_tau / n` = ta(tm$V_tau),
    `AI 2006: V_E / n` = ta(tm$V_E),
    `AI 2006: V_E / n with true sigma2 = 1` = ta(tm$V_E_true),
    `AI 2006: total` = ta(tm$var_AI2006),
    `Bootstrap (main run), approx.` = (ms("CM", "mean_ci_width") / (2 * z))^2,
    check.names = FALSE)
  save_table(vt, sprintf("sim%d_cm_variance_parts", sim), digits = 4)
  vt_mean <- data.frame(Simulation = sim, vt[1, -1, drop = FALSE], check.names = FALSE)
  vt_mean[1, -1] <- colMeans(vt[, -1])

  # sigma2_hat of AI 2006 by level, and its parts (true sigma^2 = 1)
  sg <- aggregate(g[, c("sigma2_hat", "sigma2_covariate_part", "sigma2_noise_part", "sigma2_cross")], list(Level = paste0("t", g$s)), mean)
  names(sg) <- c("Level", "sigma2_hat", "from covariate gap", "from noise", "cross term")
  sg <- rbind(sg, data.frame(Level = "Mean", t(colMeans(sg[, -1])), check.names = FALSE))
  save_table(sg, sprintf("sim%d_ai2006_sigma2", sim), digits = 3)

  pp <- data.frame(
    Pair = pk,
    `tau_pop` = te$tau_pop,
    Bias = bias,
    `Sample mix` = mix,
    `Gap, level t` = gap_t,
    `Gap, level t'` = gap_tp,
    Noise = noise,
    SD = sdv,
    `Abs. bias / SD` = abs(bias) / sdv,
    `Coverage: bootstrap` = ms("CM", "coverage"),
    `Predicted: bootstrap width` = pcov(bias, sdv, h_boot),
    `Coverage: AI 2006` = dd("CM", "AI 2006"),
    `Predicted: AI 2006 width` = pcov(bias, sdv, h_ai),
    `Reg.-adjusted: bias` = bc_bias,
    `Reg.-adjusted: bias vs tau_bal` = bc_bias - te$bal_minus_pop,
    `Reg.-adjusted: SD` = bc_sd,
    `Reg.-adjusted: coverage, oracle width` = bc_cov,
    check.names = FALSE)
  per_pair[[sim]] <- pp
  save_table(pp, sprintf("sim%d_cm_bias_coverage", sim),
             digits = c(tau_pop = 3, Bias = 3, `Sample mix` = 3, `Gap, level t` = 3, `Gap, level t'` = 3, Noise = 3, SD = 3,
                        `Abs. bias / SD` = 2, `Coverage: bootstrap` = 1, `Predicted: bootstrap width` = 1,
                        `Coverage: AI 2006` = 1, `Predicted: AI 2006 width` = 1, `Reg.-adjusted: bias` = 3,
                        `Reg.-adjusted: bias vs tau_bal` = 3, `Reg.-adjusted: SD` = 3, `Reg.-adjusted: coverage, oracle width` = 1))

  # Mean matching gap by level and covariate, and its contribution beta * gap
  gl <- data.frame(Level = paste0("t", 1:A), round(Dbar, 6), check.names = FALSE)
  save_table(gl, sprintf("sim%d_cm_mean_gap", sim), digits = 3)
  cl <- data.frame(Level = paste0("t", 1:A), contrib, Total = gap_lvl, check.names = FALSE)
  save_table(cl, sprintf("sim%d_cm_gap_contribution", sim), digits = 3)

  # The lowest-coverage pair, covariate by covariate
  w <- which.min(pp$`Coverage: bootstrap`); t1 <- pr[w, 1]; t2 <- pr[w, 2]
  lowp <- data.frame(c(COV_NAMES, "Total"), c(B[t1, ], NA), c(Dbar[t1, ], NA), c(contrib[t1, ], gap_lvl[t1]),
                     c(B[t2, ], NA), c(Dbar[t2, ], NA), c(-contrib[t2, ], -gap_lvl[t2]),
                     c(contrib[t1, ] - contrib[t2, ], gap_lvl[t1] - gap_lvl[t2]))
  names(lowp) <- c("Covariate", sprintf("beta_%d", t1), sprintf("Gap D_%d", t1), sprintf("beta_%d x D_%d", t1, t1),
                   sprintf("beta_%d", t2), sprintf("Gap D_%d", t2), sprintf("-beta_%d x D_%d", t2, t2), "Net")
  save_table(lowp, sprintf("sim%d_cm_lowest_pair_by_covariate", sim), digits = 3)

  summ_rows[[sim]] <- data.frame(
    Simulation = sim, Pairs = nrow(pp),
    `Mean abs. bias` = mean(abs(bias)),
    `Mean abs. gap term` = mean(abs(gap_t + gap_tp)),
    `Mean abs. sample-mix term` = mean(abs(mix)),
    `Mean abs. noise term` = mean(abs(noise)),
    `Mean SD` = mean(sdv),
    `Mean abs. bias / SD` = mean(abs(bias) / sdv),
    `Coverage: bootstrap` = mean(pp$`Coverage: bootstrap`),
    `Coverage: AI 2006` = mean(pp$`Coverage: AI 2006`),
    `Reg.-adjusted: mean abs. bias` = mean(abs(bc_bias)),
    `Reg.-adjusted: mean abs. bias vs tau_bal` = mean(abs(bc_bias - te$bal_minus_pop)),
    `Reg.-adjusted: mean SD` = mean(bc_sd),
    `Reg.-adjusted: coverage, oracle width` = mean(bc_cov),
    `Lowest pair` = pk[w], check.names = FALSE)
  var_rows[[sim]] <- vt_mean

  # 10: centring of the percentile interval (CM-only rerun of the bootstrap)
  fb <- file.path("crosscheck", "results", sprintf("sim%d_cm_boot_reps.csv.gz", sim))
  if (file.exists(fb)) {
    bt <- read.csv(fb); bt$tau <- te$tau_pop[match(key(bt$t, bt$tp), pk)]
    bt$shift <- bt$boot_mean - bt$est
    bt$mid <- (bt$pct_lo + bt$pct_hi) / 2 - bt$est
    bt$cov_pct <- bt$pct_lo <= bt$tau & bt$tau <= bt$pct_hi
    bt$cov_basic <- 2 * bt$est - bt$pct_hi <= bt$tau & bt$tau <= 2 * bt$est - bt$pct_lo
    bt$cov_norm <- abs(bt$est - bt$tau) <= z * bt$boot_sd
    ag <- function(x, f = mean) { a <- tapply(x, key(bt$t, bt$tp), f); unname(a[pk]) }
    bias_sub <- ag(bt$est - bt$tau)
    br <- data.frame(Pair = pk, Replications = ag(bt$rep, length),
                     `Bias (these replications)` = bias_sub,
                     `Bootstrap mean - estimate` = ag(bt$shift),
                     `Percentile midpoint - estimate` = ag(bt$mid),
                     `SD of estimates` = ag(bt$est, sd),
                     `Mean bootstrap SD` = ag(bt$boot_sd),
                     `Coverage: percentile` = 100 * ag(bt$cov_pct),
                     `Coverage: basic` = 100 * ag(bt$cov_basic),
                     `Coverage: normal` = 100 * ag(bt$cov_norm),
                     `Coverage: main run (500 replications)` = ms("CM", "coverage"),
                     check.names = FALSE)
    save_table(br, sprintf("sim%d_cm_bootstrap_centring", sim),
               digits = c(`Bias (these replications)` = 3, `Bootstrap mean - estimate` = 3, `Percentile midpoint - estimate` = 3,
                          `SD of estimates` = 3, `Mean bootstrap SD` = 3, `Coverage: percentile` = 1, `Coverage: basic` = 1,
                          `Coverage: normal` = 1, `Coverage: main run (500 replications)` = 1))
    boot_rows[[sim]] <- data.frame(
      Simulation = sim, Replications = length(unique(bt$rep)),
      `Mean abs. bias` = mean(abs(bias_sub)),
      `Mean (bootstrap mean - estimate) / bias` = mean(br$`Bootstrap mean - estimate` / bias_sub),
      `Correlation over pairs: shift vs bias` = cor(br$`Bootstrap mean - estimate`, bias_sub),
      `Mean SD of estimates` = mean(br$`SD of estimates`),
      `Mean bootstrap SD` = mean(br$`Mean bootstrap SD`),
      `Coverage: percentile` = mean(br$`Coverage: percentile`),
      `Coverage: basic` = mean(br$`Coverage: basic`),
      `Coverage: normal` = mean(br$`Coverage: normal`),
      `Coverage: main run (500 replications)` = mean(br$`Coverage: main run (500 replications)`),
      check.names = FALSE)

    # Appendix A.3: predicted coverage of the main run's percentile intervals,
    # centred on the estimate (b) and centred on the bootstrap mean (b + shift),
    # with the main run's bias, SD and mean half-width and the rerun's shift.
    p0 <- pcov(bias, sdv, h_boot); p1 <- pcov(bias + br$`Bootstrap mean - estimate`, sdv, h_boot)
    obs <- ms("CM", "coverage")
    sp <- data.frame(Pair = pk, Bias = bias, `Shift (bootstrap mean - estimate)` = br$`Bootstrap mean - estimate`,
                     SD = sdv, `Half-width` = h_boot, `Coverage: main run` = obs,
                     `Predicted: centred on estimate` = p0, `Predicted: centred on bootstrap mean` = p1, check.names = FALSE)
    save_table(sp, sprintf("sim%d_cm_bootstrap_shift_prediction", sim),
               digits = c(Bias = 3, `Shift (bootstrap mean - estimate)` = 3, SD = 3, `Half-width` = 3, `Coverage: main run` = 1,
                          `Predicted: centred on estimate` = 1, `Predicted: centred on bootstrap mean` = 1))
    shift_rows[[sim]] <- data.frame(
      Simulation = sim, `Coverage: main run` = mean(obs),
      `Predicted: centred on estimate` = mean(p0), `Predicted: centred on bootstrap mean` = mean(p1),
      `Mean abs. error: centred on estimate` = mean(abs(obs - p0)), `Mean abs. error: centred on bootstrap mean` = mean(abs(obs - p1)),
      `Correlation: centred on estimate` = cor(obs, p0), `Correlation: centred on bootstrap mean` = cor(obs, p1),
      check.names = FALSE)
  }
}

summ <- do.call(rbind, summ_rows)
save_table(summ, "cm_coverage_summary", digits = c(`Mean abs. bias` = 3, `Mean abs. gap term` = 3, `Mean abs. sample-mix term` = 3,
           `Mean abs. noise term` = 3, `Mean SD` = 3, `Mean abs. bias / SD` = 2, `Coverage: bootstrap` = 1, `Coverage: AI 2006` = 1,
           `Reg.-adjusted: mean abs. bias` = 3, `Reg.-adjusted: mean abs. bias vs tau_bal` = 3, `Reg.-adjusted: mean SD` = 3,
           `Reg.-adjusted: coverage, oracle width` = 1))
save_table(do.call(rbind, var_rows), "cm_variance_parts_summary", digits = 4)
if (length(boot_rows)) save_table(do.call(rbind, boot_rows), "cm_bootstrap_centring_summary",
                                  digits = c(`Mean abs. bias` = 3, `Mean (bootstrap mean - estimate) / bias` = 2,
                                             `Correlation over pairs: shift vs bias` = 2, `Mean SD of estimates` = 3,
                                             `Mean bootstrap SD` = 3, `Coverage: percentile` = 1, `Coverage: basic` = 1,
                                             `Coverage: normal` = 1, `Coverage: main run (500 replications)` = 1))

# ---- figure: coverage against |bias| / SD -----------------------------------
# Curves: coverage of est +- 1.96 * r * SD when est ~ N(tau + b, SD^2), for
# interval widths r = SE / SD of 1 (correct width) and 1.3 (30% too wide).
png(file.path(FIG, "fig1_coverage_vs_bias_ratio.png"), width = 1600, height = 760, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(3.2, 0, 0, 0))
xs <- seq(0, 2, by = 0.01)
for (sim in 1:2) {
  pp <- per_pair[[sim]]
  su <- read.csv(file.path("results", sprintf("sim%d_summary.csv", sim)))
  pr <- t(combn(sim_setup(sim)$A, 2)); pk <- key(pr[, 1], pr[, 2])
  gp <- su[su$estimator == "GPS" & su$target == "pop", ]; gp <- gp[match(pk, key(gp$t, gp$tp)), ]
  g16 <- d8[d8$sim == sim & d8$estimator == "GPS" & d8$method == "AI 2016" & d8$target == "pop", ]
  g16 <- g16$coverage[match(pk, key(g16$t, g16$tp))]
  plot(NA, xlim = c(0, 2), ylim = c(30, 100), xlab = "|bias| / SD of the estimates",
       ylab = "Coverage of 95% intervals (%)", main = sprintf("Simulation %d", sim), las = 1)
  abline(h = seq(30, 100, 10), col = "grey92"); abline(h = 95, lty = 2, col = "grey40")
  lines(xs, pcov(xs, 1, z), lwd = 2, col = "grey30")
  lines(xs, pcov(xs, 1, 1.3 * z), lwd = 1.5, lty = 3, col = "grey30")
  points(abs(gp$bias) / gp$sd_est, gp$coverage, pch = 1, col = COL[["GPS"]], cex = 1.1)
  points(abs(gp$bias) / gp$sd_est, g16, pch = 5, col = COL[["GPS"]], cex = 1.0)
  points(pp$`Abs. bias / SD`, pp$`Coverage: bootstrap`, pch = 1, col = COL[["CM"]], cex = 1.1, lwd = 1.3)
  points(pp$`Abs. bias / SD`, pp$`Coverage: AI 2006`, pch = 2, col = COL[["CM"]], cex = 1.0, lwd = 1.3)
  w <- which.min(pp$`Coverage: bootstrap`)
  text(pp$`Abs. bias / SD`[w], pp$`Coverage: bootstrap`[w], pp$Pair[w], pos = 2, cex = 0.85)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
legend("bottom", horiz = FALSE, ncol = 3, bty = "n", cex = 0.8, inset = c(0, 0.005),
       legend = c("CM, bootstrap (main run)", "CM, AI 2006", "GPS, bootstrap (main run)",
                  "GPS, AI 2016", "Theory: correct width", "Theory: 30% too wide"),
       pch = c(1, 2, 1, 5, NA, NA), lty = c(NA, NA, NA, NA, 1, 3), lwd = c(1.3, 1.3, 1, 1, 2, 1.5),
       col = c(COL[["CM"]], COL[["CM"]], COL[["GPS"]], COL[["GPS"]], "grey30", "grey30"))
par(op); dev.off()

if (length(shift_rows)) save_table(do.call(rbind, shift_rows), "cm_bootstrap_shift_prediction_summary",
                                   digits = c(`Coverage: main run` = 1, `Predicted: centred on estimate` = 1,
                                              `Predicted: centred on bootstrap mean` = 1, `Mean abs. error: centred on estimate` = 1,
                                              `Mean abs. error: centred on bootstrap mean` = 1, `Correlation: centred on estimate` = 3,
                                              `Correlation: centred on bootstrap mean` = 3))

# ---- the paper's application intervals (Table 5): what form are they? ------
# Checks whether each 95% limit equals estimate -/+ z * sd, with z = 1.96 and
# with z = qnorm(0.975), using the values transcribed in results/paper_table5.csv.
p5 <- read.csv(file.path("results", "paper_table5.csv"))
t5 <- do.call(rbind, lapply(c("CM", "GPS"), function(e) {
  est <- p5[[paste0(e, "_est")]]; s <- p5[[paste0(e, "_sd")]]; lo <- p5[[paste0(e, "_LL")]]; hi <- p5[[paste0(e, "_UL")]]
  err <- function(zz) max(abs(c(lo - (est - zz * s), hi - (est + zz * s))))
  data.frame(Estimator = e, Pairs = nrow(p5),
             `Max abs. (midpoint - estimate)` = max(abs((lo + hi) / 2 - est)),
             `Max abs. error, est +- 1.96 sd` = err(1.96),
             `Max abs. error, est +- qnorm(0.975) sd` = err(qnorm(0.975)),
             `Implied z: min` = min((hi - lo) / (2 * s)), `Implied z: max` = max((hi - lo) / (2 * s)), check.names = FALSE)
}))
save_table(t5, "paper_table5_interval_form", digits = 6)

cat("CM coverage summary:\n"); print(summ, digits = 3)
if (length(boot_rows)) { cat("\nBootstrap centring:\n"); print(do.call(rbind, boot_rows), digits = 3) }
