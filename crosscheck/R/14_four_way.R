# 14_four_way.R -- four-way comparison for crosscheck/four_way/README.md:
# the paper (Guo, Liu & Wang, 2022), this rebuild (R1), the independent
# rebuild-2 (R2) and Yang et al.'s package (multilevelMatching, run on R1's
# samples by 04_crosscheck_all.R).
#
# Inputs (all already in this repository):
#   results/paper_table{1,3}.csv            the paper's tables (R1's transcription)
#   crosscheck/four_way/data/               R2's own transcription of the same tables (copied)
#   results/sim{1,2}_summary.csv            R1, main run
#   crosscheck/rebuild2/data/               R2's performance and coverage (copied; see its README)
#   crosscheck/results/sim{1,2}_crosscheck_summary.csv   the package (05_summarize_crosscheck.R)
#
#   Rscript crosscheck/R/14_four_way.R      (from the repo root; base R only)
# Writes crosscheck/four_way/tables/*.{csv,md} and crosscheck/four_way/figures/*.png.
#
# Conventions: bias = E(tau_hat) - tau, tau(t,t') = E[Y(t) - Y(t')], against the
# population effect (R1's tau_pop = R2's tau_super). The paper's Table 3 signs
# are reversed relative to this (docs/two_way.md; R2's D5), so its signed biases
# are multiplied by -1 ("sign-aligned"). Table 1 gives absolute bias only.

OUT <- file.path("crosscheck", "four_way")
TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
COL <- c(DIF = "red", CM = "green3", GPS = "blue")      # the paper's colours, as in R/plots.R
SRC <- c("paper", "R1", "R2", "package")
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

# ---- 1. the paper's tables, transcribed twice -----------------------------------
p1 <- read.csv(file.path("results", "paper_table1.csv")); q1 <- read.csv(file.path(DAT, "rebuild2_paper_table1_sim1.csv"))
p3 <- read.csv(file.path("results", "paper_table3.csv")); q3 <- read.csv(file.path(DAT, "rebuild2_paper_table3_sim2.csv"))
m1 <- match(p1$pair, paste(q1$t, q1$tprime, sep = ",")); m3 <- match(p3$pair, paste(q3$t, q3$tprime, sep = ","))
stopifnot(!anyNA(m1), !anyNA(m3))
d1 <- cbind(p1$DIF_abs_bias - q1$dif_absbias[m1], p1$CM_abs_bias - q1$cm_absbias[m1], p1$CM_coverage - q1$cm_coverage[m1],
            p1$GPS_abs_bias - q1$gps_absbias[m1], p1$GPS_coverage - q1$gps_coverage[m1])
d3 <- cbind(p3$DIF_bias - q3$dif_bias[m3], p3$DIF_mse - q3$dif_mse[m3], p3$CM_bias - q3$cm_bias[m3], p3$CM_mse - q3$cm_mse[m3],
            p3$CM_coverage - q3$cm_coverage[m3], p3$GPS_bias - q3$gps_bias[m3], p3$GPS_mse - q3$gps_mse[m3],
            p3$GPS_coverage - q3$gps_coverage[m3])
save_table(data.frame(Table = c("Table 1 (Simulation 1)", "Table 3 (Simulation 2)"), Pairs = c(nrow(d1), nrow(d3)),
                      `Cells compared` = c(length(d1), length(d3)), `Cells that differ` = c(sum(d1 != 0), sum(d3 != 0)),
                      `Largest difference` = c(max(abs(d1)), max(abs(d3))), check.names = FALSE), "transcription_check")

# ---- 2. per-pair results for the four sources ----------------------------------
load_sim <- function(sim) {
  r1 <- read.csv(file.path("results", sprintf("sim%d_summary.csv", sim))); r1 <- r1[r1$target == "pop", ]
  r2 <- read.csv(file.path("crosscheck", "rebuild2", "data", sprintf("rebuild2_sim%d_performance.csv", sim))); r2 <- r2[r2$truth == "tau_super", ]
  r2c <- read.csv(file.path("crosscheck", "rebuild2", "data", sprintf("rebuild2_sim%d_coverage.csv", sim)))
  pk <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_crosscheck_summary.csv", sim))); pk <- pk[pk$target == "pop", ]
  pp <- if (sim == 1) p1 else p3
  pairs <- unique(r1[c("t", "tp")]); pairs <- pairs[order(pairs$t, pairs$tp), ]; k <- key(pairs$t, pairs$tp)
  g <- function(d, est, col, tcol = "tp") { x <- d[d$estimator == est, ]; x[[col]][match(k, key(x$t, x[[tcol]]))] }
  pm <- match(k, paste0("(", pp$pair, ")"))
  out <- data.frame(sim = sim, pair = k)
  for (e in c("DIF", "CM", "GPS")) {
    out[[paste0(e, "_bias_paper")]] <- if (sim == 2) -pp[[paste0(e, "_bias")]][pm] else NA       # sign-aligned
    out[[paste0(e, "_absbias_paper")]] <- if (sim == 2) abs(pp[[paste0(e, "_bias")]][pm]) else pp[[paste0(e, "_abs_bias")]][pm]
    out[[paste0(e, "_mse_paper")]] <- if (sim == 2 && e != "DIF") pp[[paste0(e, "_mse")]][pm] else NA   # DIF "MSE" is not an MSE
    out[[paste0(e, "_cov_paper")]] <- if (e != "DIF") pp[[paste0(e, "_coverage")]][pm] else NA
    out[[paste0(e, "_bias_R1")]] <- g(r1, e, "bias"); out[[paste0(e, "_absbias_R1")]] <- g(r1, e, "abs_bias")
    out[[paste0(e, "_mse_R1")]] <- g(r1, e, "mse"); out[[paste0(e, "_cov_R1")]] <- g(r1, e, "coverage")
    out[[paste0(e, "_mcse_R1")]] <- g(r1, e, "sd_est") / sqrt(500)
    out[[paste0(e, "_bias_R2")]] <- g(r2, e, "bias", "tprime"); out[[paste0(e, "_absbias_R2")]] <- g(r2, e, "abs_bias", "tprime")
    out[[paste0(e, "_mse_R2")]] <- g(r2, e, "mse", "tprime"); out[[paste0(e, "_cov_R2")]] <- 100 * g(r2c, e, "cov_super", "tprime")
    out[[paste0(e, "_mcse_R2")]] <- g(r2, e, "mcse_bias", "tprime")
    pe <- c(DIF = NA, CM = "pkg_CM", GPS = "pkg_GPS_multinom")[[e]]
    for (m in c("bias", "absbias", "mse", "cov", "mcse")) out[[paste0(e, "_", m, "_package")]] <- NA
    if (!is.na(pe)) {
      out[[paste0(e, "_bias_package")]] <- g(pk, pe, "bias"); out[[paste0(e, "_absbias_package")]] <- g(pk, pe, "abs_bias")
      out[[paste0(e, "_mse_package")]] <- g(pk, pe, "mse"); out[[paste0(e, "_cov_package")]] <- g(pk, pe, "coverage")
      out[[paste0(e, "_mcse_package")]] <- g(pk, pe, "mc_se")
      if (e == "GPS") out$GPS_cov_package_AI2016 <- g(pk, "pkg_GPS_multinom_AI2016", "coverage")
    }
  }
  out
}
S <- list(load_sim(1), load_sim(2))
for (sim in 1:2) write.csv(S[[sim]], file.path(TAB, sprintf("sim%d_four_way_by_pair.csv", sim)), row.names = FALSE)

# ---- 3. headline: mean absolute bias and MSE -----------------------------------
mn <- function(v) if (all(is.na(v))) NA else mean(v)
head_tab <- do.call(rbind, lapply(1:2, function(sim) do.call(rbind, lapply(c("DIF", "CM", "GPS"), function(e) {
  x <- S[[sim]]; v <- function(m, s) mn(x[[paste0(e, "_", m, "_", s)]])
  data.frame(Simulation = sim, Estimator = e,
             `Abs. bias: paper` = v("absbias", "paper"), `Abs. bias: R1` = v("absbias", "R1"),
             `Abs. bias: R2` = v("absbias", "R2"), `Abs. bias: package` = v("absbias", "package"),
             `MSE: paper` = v("mse", "paper"), `MSE: R1` = v("mse", "R1"), `MSE: R2` = v("mse", "R2"),
             `MSE: package` = v("mse", "package"), check.names = FALSE) }))))
save_table(head_tab, "bias_mse_four_way")

# ---- 4. which estimator is less biased, and which more efficient ----------------
ord <- do.call(rbind, lapply(1:2, function(sim) {
  x <- S[[sim]]
  cnt <- function(a, b, s) { u <- x[[paste0(a, "_", s)]]; w <- x[[paste0(b, "_", s)]]; if (all(is.na(u)) || all(is.na(w))) NA else sum(u < w) }
  rbind(data.frame(Simulation = sim, Comparison = "GPS less biased than CM",
                   Paper = cnt("GPS_absbias", "CM_absbias", "paper"), R1 = cnt("GPS_absbias", "CM_absbias", "R1"),
                   R2 = cnt("GPS_absbias", "CM_absbias", "R2"), Package = cnt("GPS_absbias", "CM_absbias", "package"), Pairs = nrow(x)),
        data.frame(Simulation = sim, Comparison = "CM lower MSE than GPS",
                   Paper = cnt("CM_mse", "GPS_mse", "paper"), R1 = cnt("CM_mse", "GPS_mse", "R1"),
                   R2 = cnt("CM_mse", "GPS_mse", "R2"), Package = cnt("CM_mse", "GPS_mse", "package"), Pairs = nrow(x)))
}))
save_table(ord, "orderings_four_way")

# ---- 5. agreement between sources, pair by pair --------------------------------
# Simulation 2: signed bias (paper sign-aligned). Simulation 1: absolute bias
# (the only form Table 1 gives). For R1 vs R2 and R2 vs package (independent
# samples), the difference is also scaled by its Monte Carlo SE.
agree <- list()
for (sim in 1:2) for (e in c("CM", "GPS")) {
  x <- S[[sim]]; m <- if (sim == 2) "bias" else "absbias"
  for (i in 1:3) for (j in (i + 1):4) {
    a <- SRC[i]; b <- SRC[j]; u <- x[[paste0(e, "_", m, "_", a)]]; w <- x[[paste0(e, "_", m, "_", b)]]
    indep <- (a == "R1" && b == "R2") || (a == "R2" && b == "package")
    zz <- if (indep) abs(x[[paste0(e, "_bias_", b)]] - x[[paste0(e, "_bias_", a)]]) /
                     sqrt(x[[paste0(e, "_mcse_", a)]]^2 + x[[paste0(e, "_mcse_", b)]]^2) else NA
    agree[[length(agree) + 1]] <- data.frame(Simulation = sim, Estimator = e,
      Compared = if (sim == 2) "signed bias" else "absolute bias", Sources = paste(a, "vs", b),
      Correlation = cor(u, w), `Mean abs. difference` = mean(abs(u - w)),
      `Pairs within 2 MC SE` = if (indep) sum(zz <= 2) else NA, Pairs = nrow(x), check.names = FALSE)
  }
}
agree <- do.call(rbind, agree)
save_table(agree, "agreement_by_source_pair")
# Compact: correlations only, one row per pair of sources.
agree$col <- paste0("Sim. ", agree$Simulation, " ", agree$Estimator, ifelse(agree$Simulation == 1, " (abs. bias)", " (signed bias)"))
cw <- reshape(agree[c("Sources", "col", "Correlation")], idvar = "Sources", timevar = "col", direction = "wide")
names(cw) <- sub("^Correlation\\.", "", names(cw)); rownames(cw) <- NULL
save_table(cw, "agreement_correlations")

# ---- 6. Simulation 1 by pair --------------------------------------------------
x <- S[[1]]
save_table(data.frame(Pair = x$pair,
  `CM paper` = x$CM_absbias_paper, `CM R1` = x$CM_absbias_R1, `CM R2` = x$CM_absbias_R2, `CM package` = x$CM_absbias_package,
  `GPS paper` = x$GPS_absbias_paper, `GPS R1` = x$GPS_absbias_R1, `GPS R2` = x$GPS_absbias_R2, `GPS package` = x$GPS_absbias_package,
  check.names = FALSE), "sim1_abs_bias_four_way")

# ---- 7. coverage (secondary) --------------------------------------------------
cov_tab <- do.call(rbind, lapply(1:2, function(sim) do.call(rbind, lapply(c("CM", "GPS"), function(e) {
  x <- S[[sim]]
  data.frame(Simulation = sim, Estimator = e, Paper = mn(x[[paste0(e, "_cov_paper")]]),
             `R1: bootstrap` = mn(x[[paste0(e, "_cov_R1")]]), `R2: analytic` = mn(x[[paste0(e, "_cov_R2")]]),
             `Package: AI 2006` = mn(x[[paste0(e, "_cov_package")]]),
             `Package: AI 2016` = if (e == "GPS") mn(x$GPS_cov_package_AI2016) else NA, check.names = FALSE) }))))
save_table(cov_tab, "coverage_four_way", digits = 1)

# ---- figures --------------------------------------------------------------------
sty <- list(paper = list(pch = 1, lty = 3, col = function(e) "black", lwd = 1.2),
            R1 = list(pch = 16, lty = 1, col = function(e) COL[[e]], lwd = 1.8),
            R2 = list(pch = 0, lty = 2, col = function(e) COL[[e]], lwd = 1.3),
            package = list(pch = 17, lty = 4, col = function(e) "grey25", lwd = 1.2))
lab <- c(paper = "paper (sign-aligned)", R1 = "R1 (this rebuild)", R2 = "R2 (rebuild-2)", package = "package (multilevelMatching)")

png(file.path(FIG, "fig1_sim2_bias_four_way.png"), width = 1600, height = 1250, res = 150)
op <- par(mfrow = c(2, 1), mar = c(5.5, 4.5, 2.5, 1), oma = c(0, 0, 2.2, 0))
for (e in c("CM", "GPS")) {
  x <- S[[2]]; k <- nrow(x); v <- sapply(SRC, function(s) x[[paste0(e, "_bias_", s)]])
  plot(NA, xlim = c(0.5, k + 0.5), ylim = range(v) * c(1.1, 1.45), xaxt = "n", xlab = "",
       ylab = expression("Bias vs " * tau[pop]), main = e, las = 1)
  abline(h = axTicks(2), col = "grey92"); abline(h = 0, col = "grey50")
  axis(1, at = 1:k, labels = x$pair, las = 2, cex.axis = 0.7)
  for (s in SRC) { st <- sty[[s]]; lines(1:k, v[, s], lty = st$lty, lwd = st$lwd, col = st$col(e))
                   points(1:k, v[, s], pch = st$pch, col = st$col(e), cex = 0.8, lwd = 1.2) }
  legend("top", horiz = TRUE, bty = "n", cex = 0.75, pch = sapply(sty, `[[`, "pch"), lty = sapply(sty, `[[`, "lty"),
         col = sapply(SRC, function(s) sty[[s]]$col(e)), legend = lab[SRC])
}
mtext(expression(bold("Simulation 2: bias of the 36 effects against ") * bold(tau[pop]) * bold(", paper vs R1 vs R2 vs package")),
      outer = TRUE, cex = 1.05)
par(op); dev.off()

png(file.path(FIG, "fig2_sim1_abs_bias_four_way.png"), width = 1600, height = 750, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 2.2, 0))
off <- c(paper = -0.27, R1 = -0.09, R2 = 0.09, package = 0.27)
for (e in c("CM", "GPS")) {
  x <- S[[1]]; k <- nrow(x); v <- sapply(SRC, function(s) x[[paste0(e, "_absbias_", s)]])
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0, max(v) * 1.6), xaxt = "n", xlab = "Pair (t, t')",
       ylab = expression("Absolute bias vs " * tau[pop]), main = e, las = 1)
  abline(h = axTicks(2), col = "grey92"); abline(h = 0, col = "grey50")
  axis(1, at = 1:k, labels = x$pair, cex.axis = 0.9)
  for (s in SRC) points((1:k) + off[[s]], v[, s], pch = sty[[s]]$pch, col = sty[[s]]$col(e), cex = 1.3, lwd = 1.4)
  legend("topleft", bty = "n", cex = 0.8, pch = sapply(sty, `[[`, "pch"), pt.lwd = 1.4,
         col = sapply(SRC, function(s) sty[[s]]$col(e)), legend = sub(" \\(sign-aligned\\)", " (Table 1)", lab[SRC]))
}
mtext(expression(bold("Simulation 1: absolute bias by pair against ") * bold(tau[pop]) * bold(", paper vs R1 vs R2 vs package")),
      outer = TRUE, cex = 1.05)
par(op); dev.off()

png(file.path(FIG, "fig3_mse_four_way.png"), width = 1600, height = 1300, res = 150)
op <- par(mfrow = c(2, 2), mar = c(5.5, 4.5, 2.5, 1), oma = c(0, 0, 2.2, 0))
for (sim in 1:2) for (e in c("CM", "GPS")) {
  x <- S[[sim]]; k <- nrow(x); src <- if (sim == 1) SRC[-1] else SRC
  v <- sapply(src, function(s) x[[paste0(e, "_mse_", s)]])
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0, max(v) * (if (sim == 1) 1.3 else 1.45)), xaxt = "n", xlab = "", ylab = "MSE",
       main = sprintf("Simulation %d, %s", sim, e), las = 1)
  abline(h = axTicks(2), col = "grey92"); abline(h = 0, col = "grey50")
  axis(1, at = 1:k, labels = x$pair, las = 2, cex.axis = if (sim == 1) 0.9 else 0.6)
  o <- if (sim == 1) c(R1 = -0.2, R2 = 0, package = 0.2) else c(paper = 0, R1 = 0, R2 = 0, package = 0)
  for (s in src) { st <- sty[[s]]
    if (sim == 2) lines(1:k, v[, s], lty = st$lty, lwd = st$lwd, col = st$col(e))
    points((1:k) + o[[s]], v[, s], pch = st$pch, col = st$col(e), cex = if (sim == 1) 1.2 else 0.7, lwd = 1.2) }
  legend(if (sim == 1) "top" else "topleft", horiz = sim == 1, ncol = if (sim == 1) 1 else 2, bty = "n",
         cex = if (sim == 1) 0.62 else 0.72, pch = sapply(sty[src], `[[`, "pch"),
         lty = if (sim == 2) sapply(sty[src], `[[`, "lty") else rep(0, length(src)),
         col = sapply(src, function(s) sty[[s]]$col(e)), legend = sub(" \\(sign-aligned\\)", " (Table 3)", lab[src]))
}
mtext(expression(bold("MSE by pair against ") * bold(tau[pop]) * bold(": paper (Simulation 2 only) vs R1 vs R2 vs package")),
      outer = TRUE, cex = 1.05)
par(op); dev.off()

cat("Wrote tables:\n"); print(list.files(TAB)); cat("Wrote figures:\n"); print(list.files(FIG))
