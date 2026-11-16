# 13_rebuild2_comparison.R -- tables and figure for crosscheck/rebuild2/README.md:
# the paper, this rebuild (R1) and the independent rebuild-2 (R2), and the
# bias-adjusted CM of 12_cm_biasadjust.R.
#
# Inputs: results/ (main run, paper tables), crosscheck/results/ (04, 09, 12),
# crosscheck/side_by_side/tables/d8_coverage_by_pair.csv (08),
# crosscheck/intervals/tables/sim<k>_cm_bias_coverage.csv (11) and
# crosscheck/rebuild2/data/ (R2's summaries, copied; see its README).
#
#   Rscript crosscheck/R/13_rebuild2_comparison.R
# Writes crosscheck/rebuild2/tables/*.{csv,md} and crosscheck/rebuild2/figures/*.png.

source("R/dgp.R")
OUT <- file.path("crosscheck", "rebuild2")
TAB <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures"); DAT <- file.path(OUT, "data")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
COL <- c(DIF = "red", CM = "green3", GPS = "blue")      # the paper's colours, as in R/plots.R
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

d8 <- read.csv(file.path("crosscheck", "side_by_side", "tables", "d8_coverage_by_pair.csv"))
load_sim <- function(sim) {
  A <- sim_setup(sim)$A; pr <- t(combn(A, 2)); pk <- key(pr[, 1], pr[, 2])
  te <- read.csv(file.path("results", sprintf("true_effects_sim%d.csv", sim))); te <- te[match(pk, key(te$t, te$tp)), ]
  su <- read.csv(file.path("results", sprintf("sim%d_summary.csv", sim))); su <- su[su$target == "pop", ]
  r1 <- function(e, col) { x <- su[su$estimator == e, ]; x[[col]][match(pk, key(x$t, x$tp))] }
  pf <- read.csv(file.path(DAT, sprintf("rebuild2_sim%d_performance.csv", sim))); pf <- pf[pf$truth == "tau_super", ]
  r2 <- function(e, col) { x <- pf[pf$estimator == e, ]; x[[col]][match(pk, key(x$t, x$tprime))] }
  cv <- read.csv(file.path(DAT, sprintf("rebuild2_sim%d_coverage.csv", sim)))
  r2c <- function(e, col = "cov_super") { x <- cv[cv$estimator == e, ]; x[[col]][match(pk, key(x$t, x$tprime))] }
  d8c <- function(e, m) { x <- d8[d8$sim == sim & d8$estimator == e & d8$method == m & d8$target == "pop", ]; x$coverage[match(pk, key(x$t, x$tp))] }
  b <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_cm_biasadjust_reps.csv", sim)))
  kb <- key(b$t, b$tp); tau <- te$tau_pop[match(kb, pk)]; bal <- te$tau_bal[match(kb, pk)]
  agb <- function(x, f = mean) unname(tapply(x, kb, f)[pk])
  ba <- list(bias = agb(b$CM_BiasAdjust - tau), bias_bal = agb(b$CM_BiasAdjust - bal),
             mse = agb((b$CM_BiasAdjust - tau)^2), sd = agb(b$CM_BiasAdjust, sd))
  if (sim == 1) {
    p <- read.csv(file.path("results", "paper_table1.csv")); kp <- match(pk, paste0("(", p$pair, ")"))
    paper <- list(CM_abs = p$CM_abs_bias[kp], GPS_abs = p$GPS_abs_bias[kp], DIF_abs = p$DIF_abs_bias[kp],
                  CM_cov = p$CM_coverage[kp], GPS_cov = p$GPS_coverage[kp])
  } else {
    p <- read.csv(file.path("results", "paper_table3.csv")); kp <- match(pk, paste0("(", p$pair, ")"))
    paper <- list(CM_bias = -p$CM_bias[kp], GPS_bias = -p$GPS_bias[kp], DIF_bias = -p$DIF_bias[kp],   # sign-aligned (D5 / two_way.md)
                  CM_abs = abs(p$CM_bias[kp]), GPS_abs = abs(p$GPS_bias[kp]), DIF_abs = abs(p$DIF_bias[kp]),
                  CM_mse = p$CM_mse[kp], GPS_mse = p$GPS_mse[kp], CM_cov = p$CM_coverage[kp], GPS_cov = p$GPS_coverage[kp])
  }
  list(sim = sim, pk = pk, te = te, r1 = r1, r2 = r2, r2c = r2c, d8c = d8c, ba = ba, paper = paper)
}
S <- lapply(1:2, load_sim)

# ---- 1. Bias and MSE, means over pairs ------------------------------------
t1 <- do.call(rbind, lapply(S, function(s) do.call(rbind, lapply(c("DIF", "CM", "GPS", "GPS_cs", "GPS_true", "CM, bias-adjusted"), function(e) {
  pa <- if (e %in% c("DIF", "CM", "GPS")) mean(s$paper[[paste0(e, "_abs")]]) else NA
  pm <- if (s$sim == 2 && e %in% c("CM", "GPS")) mean(s$paper[[paste0(e, "_mse")]]) else NA
  r1a <- if (e %in% c("DIF", "CM", "GPS")) mean(abs(s$r1(e, "bias"))) else if (e == "CM, bias-adjusted") mean(abs(s$ba$bias)) else NA
  r1m <- if (e %in% c("DIF", "CM", "GPS")) mean(s$r1(e, "mse")) else if (e == "CM, bias-adjusted") mean(s$ba$mse) else NA
  r2a <- if (e != "CM, bias-adjusted") mean(abs(s$r2(e, "bias"))) else NA
  r2m <- if (e != "CM, bias-adjusted") mean(s$r2(e, "mse")) else NA
  data.frame(Simulation = s$sim, Estimator = e, `Abs. bias: paper` = pa, `Abs. bias: R1` = r1a, `Abs. bias: R2` = r2a,
             `MSE: paper` = pm, `MSE: R1` = r1m, `MSE: R2` = r2m, check.names = FALSE)
}))))
save_table(t1, "bias_mse_summary", digits = 3)

# ---- 2. R1 against R2, pair by pair ----------------------------------------
t2 <- do.call(rbind, lapply(S, function(s) do.call(rbind, lapply(c("DIF", "CM", "GPS"), function(e) {
  b1 <- s$r1(e, "bias"); b2 <- s$r2(e, "bias")
  se <- sqrt((s$r1(e, "sd_est") / sqrt(500))^2 + s$r2(e, "mcse_bias")^2)   # MC SE of the difference
  m1 <- s$r1(e, "mse"); m2 <- s$r2(e, "mse")
  data.frame(Simulation = s$sim, Estimator = e, Pairs = length(b1),
             `Bias: correlation` = cor(b1, b2), `Bias: mean abs. difference` = mean(abs(b1 - b2)),
             `Bias: mean MC SE of difference` = mean(se), `Pairs within 2 MC SE` = sum(abs(b1 - b2) <= 2 * se),
             `MSE: correlation` = cor(m1, m2), `MSE: R2 / R1, ratio of means` = mean(m2) / mean(m1), check.names = FALSE)
}))))
save_table(t2, "r1_vs_r2_by_pair", digits = c(`Bias: correlation` = 3, `Bias: mean abs. difference` = 3,
           `Bias: mean MC SE of difference` = 3, `MSE: correlation` = 3, `MSE: R2 / R1, ratio of means` = 3))

# ---- 3. Coverage ------------------------------------------------------------
t3 <- do.call(rbind, lapply(S, function(s) do.call(rbind, lapply(c("CM", "GPS"), function(e) {
  boot <- s$r1(e, "coverage"); ai06 <- s$d8c(e, "AI 2006"); ai16 <- if (e == "GPS") s$d8c(e, "AI 2016") else NA
  r2 <- 100 * s$r2c(e); r2t <- if (e == "GPS") 100 * s$r2c("GPS_true") else NA
  data.frame(Simulation = s$sim, Estimator = e, Paper = mean(s$paper[[paste0(e, "_cov")]]),
             `R1: bootstrap` = mean(boot), `R1: AI 2006` = mean(ai06), `R1: AI 2016` = mean(ai16),
             `R2: analytic (its D26)` = mean(r2), `R2: GPS with true scores` = mean(r2t),
             `R1 AI 2006 vs R2, per pair: correlation` = cor(ai06, r2), `R1 AI 2006 vs R2: mean abs. difference` = mean(abs(ai06 - r2)),
             `R2: SE / SD` = mean(s$r2c(e, "se_ratio")), check.names = FALSE)
}))))
save_table(t3, "coverage_summary", digits = c(Paper = 1, `R1: bootstrap` = 1, `R1: AI 2006` = 1, `R1: AI 2016` = 1,
           `R2: analytic (its D26)` = 1, `R2: GPS with true scores` = 1, `R1 AI 2006 vs R2, per pair: correlation` = 3,
           `R1 AI 2006 vs R2: mean abs. difference` = 1, `R2: SE / SD` = 2))
low <- do.call(rbind, lapply(S, function(s) {
  ai06 <- s$d8c("CM", "AI 2006"); r2 <- 100 * s$r2c("CM"); o1 <- order(ai06)[1:4]; o2 <- order(r2)[1:4]
  data.frame(Simulation = s$sim, Rank = 1:4, `R1 (AI 2006): pair` = s$pk[o1], `R1: coverage` = ai06[o1],
             `R2 (analytic): pair` = s$pk[o2], `R2: coverage` = r2[o2], check.names = FALSE)
}))
save_table(low, "cm_lowest_coverage_pairs", digits = 1)

# ---- 4. Simulation 2: which bias column tracks which -------------------------
s <- S[[2]]
ia <- read.csv(file.path("crosscheck", "intervals", "tables", "sim2_cm_bias_coverage.csv"), check.names = FALSE)
ia <- ia[match(s$pk, ia$Pair), ]
mix <- ia[["Sample mix"]]; gap <- ia[["Gap, level t"]] + ia[["Gap, level t'"]]
M <- cbind(`Paper CM` = s$paper$CM_bias, `Paper GPS` = s$paper$GPS_bias, `R1 CM` = s$r1("CM", "bias"), `R1 GPS` = s$r1("GPS", "bias"),
           `R2 CM` = s$r2("CM", "bias"), `R2 GPS` = s$r2("GPS", "bias"), `R1 CM, bias-adjusted` = s$ba$bias,
           `Sample-mix term (R1)` = mix)
cm_mat <- data.frame(Column = colnames(M), round(cor(M), 3), check.names = FALSE)
save_table(cm_mat, "sim2_bias_correlations", digits = 2)

fit <- function(y, lab) { f <- lm(y ~ 0 + mix + gap); cf <- coef(summary(f))
  data.frame(`Bias column` = lab, `Coefficient on sample-mix term` = cf["mix", 1], SE = cf["mix", 2],
             `Coefficient on matching-gap term` = cf["gap", 1], `SE ` = cf["gap", 2], `R^2` = summary(f)$r.squared, check.names = FALSE) }
t5 <- rbind(fit(s$paper$CM_bias, "Paper CM (sign-aligned)"), fit(s$r1("CM", "bias"), "R1 CM"),
            fit(s$r2("CM", "bias"), "R2 CM"), fit(s$ba$bias, "R1 CM, bias-adjusted"))
save_table(t5, "sim2_cm_bias_on_terms", digits = 3)

# ---- 5. The bias-adjusted CM against the paper's CM ---------------------------
t6 <- do.call(rbind, lapply(S, function(s) {
  r1b <- s$r1("CM", "bias"); r1sd <- s$r1("CM", "sd_est")
  pv <- if (s$sim == 2) s$paper$CM_mse - s$paper$CM_bias^2 else NA
  data.frame(Simulation = s$sim,
             `Abs. bias: paper` = mean(s$paper$CM_abs), `Abs. bias: CM` = mean(abs(r1b)), `Abs. bias: bias-adjusted` = mean(abs(s$ba$bias)),
             `Bias-adjusted, abs. bias vs tau_bal` = mean(abs(s$ba$bias_bal)),
             `Corr. with paper's abs. bias: CM` = cor(s$paper$CM_abs, abs(r1b)), `Corr.: bias-adjusted` = cor(s$paper$CM_abs, abs(s$ba$bias)),
             `MSE: paper` = if (s$sim == 2) mean(s$paper$CM_mse) else NA, `MSE: CM` = mean(s$r1("CM", "mse")), `MSE: bias-adjusted` = mean(s$ba$mse),
             `Variance: paper (MSE - bias^2)` = mean(pv), `Variance: CM` = mean(r1sd^2), `Variance: bias-adjusted` = mean(s$ba$sd^2),
             check.names = FALSE)
}))
save_table(t6, "cm_biasadjust_vs_paper", digits = 3)
sim1_pairs <- data.frame(Pair = S[[1]]$pk, `Paper` = S[[1]]$paper$CM_abs, `R1 CM` = abs(S[[1]]$r1("CM", "bias")),
                         `R2 CM` = abs(S[[1]]$r2("CM", "bias")), `R1 CM, bias-adjusted` = abs(S[[1]]$ba$bias), check.names = FALSE)
save_table(sim1_pairs, "sim1_cm_abs_bias_by_pair", digits = 3)

# ---- figure ------------------------------------------------------------------
png(file.path(FIG, "fig1_r1_r2_paper_bias.png"), width = 1600, height = 780, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.6, 4.6, 2.6, 1), oma = c(3, 0, 0, 0))
lim <- range(c(s$r1("DIF", "bias"), s$r2("DIF", "bias")))
plot(NA, xlim = lim, ylim = lim, xlab = "R1: bias (36 pairs, Simulation 2)", ylab = "R2: bias", las = 1,
     main = "Two independent rebuilds, pair by pair")
abline(0, 1, col = "grey60"); abline(h = 0, v = 0, col = "grey90")
for (e in c("DIF", "CM", "GPS")) points(s$r1(e, "bias"), s$r2(e, "bias"), pch = c(DIF = 0, CM = 1, GPS = 2)[[e]], col = COL[[e]], lwd = 1.3)
lim2 <- range(c(s$paper$CM_bias, s$r1("CM", "bias"), s$ba$bias))
plot(NA, xlim = lim2, ylim = lim2, xlab = "R1: bias, CM and bias-adjusted CM", ylab = "Paper: CM bias (sign-aligned)", las = 1,
     main = "The paper's CM column, Simulation 2")
abline(0, 1, col = "grey60"); abline(h = 0, v = 0, col = "grey90")
points(s$r1("CM", "bias"), s$paper$CM_bias, pch = 1, col = COL[["CM"]], lwd = 1.3)
points(s$ba$bias, s$paper$CM_bias, pch = 17, col = "grey25")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
legend("bottom", ncol = 5, bty = "n", cex = 0.8, inset = c(0, 0.005),
       legend = c("DIF", "CM", "GPS", "R1 CM, bias-adjusted (right)", "45-degree line"),
       pch = c(0, 1, 2, 17, NA), lty = c(NA, NA, NA, NA, 1), col = c(COL[["DIF"]], COL[["CM"]], COL[["GPS"]], "grey25", "grey60"))
par(op); dev.off()

print(t1, digits = 3); print(t2, digits = 3); print(t3, digits = 3); print(t5, digits = 3); print(t6, digits = 3)
