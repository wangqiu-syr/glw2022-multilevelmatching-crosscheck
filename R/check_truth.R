# check_truth.R -- Stage 1 sanity checks for the DGP and both true-effect
# targets (docs/spec.md, D4 and D9). Run from the repository root:
#   Rscript R/check_truth.R
# Writes results/true_effects_sim1.csv and results/true_effects_sim2.csv.

source("R/dgp.R")
set.seed(20221104, kind = "L'Ecuyer-CMRG")
dir.create("results", showWarnings = FALSE)

N_POP     <- 1e6   # population size for the tau_pop Monte Carlo check
N_BAL     <- 2e6   # pool size for the p_t-weighted mu_bal estimate
R_BAL     <- 300   # number of balanced samples for the direct tau_bal check
R_POOL    <- 50    # number of balanced samples used to profile pool sizes
Z_TOL     <- 4     # |z| threshold for Monte Carlo agreement

all_ok <- TRUE
report <- function(label, ok) {
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) all_ok <<- FALSE
}

pair_diff <- function(Ypo, pr) Ypo[, pr[, 1], drop = FALSE] - Ypo[, pr[, 2], drop = FALSE]

for (sim in 1:2) {
  s  <- sim_setup(sim)
  pr <- t(combn(s$A, 2))
  cat(sprintf("\n================ Simulation %d (A = %d, n_g = %d, %d pairs) ================\n",
              sim, s$A, s$n_g, nrow(pr)))

  # 1. Treatment probabilities ------------------------------------------------
  V <- gen_covariates(1e5)
  P <- treatment_probs(V, s)
  cat("\n1. Treatment probabilities (1e5 units)\n")
  cat(sprintf("  max |sum_t p_t - 1| = %.2e;  range of p_t = [%.4f, %.4f]\n",
              max(abs(rowSums(P) - 1)), min(P), max(P)))
  report("p_t sums to 1 and lies in (0, 1)",
         max(abs(rowSums(P) - 1)) < 1e-12 && min(P) > 0 && max(P) < 1)
  if (sim == 2) {
    cat(sprintf("  share of units with any p_t < 0.01: %.4f  (any p_t < 0.001: %.4f)\n",
                mean(apply(P, 1, min) < 0.01), mean(apply(P, 1, min) < 0.001)))
    q <- apply(P, 2, quantile, probs = c(0.01, 0.05, 0.50))
    cat("  percentiles of p_t by group:\n")
    print(data.frame(group = colnames(P), p01 = signif(q[1, ], 3),
                     p05 = signif(q[2, ], 3), median = signif(q[3, ], 3),
                     share_below_0.01 = round(colMeans(P < 0.01), 4)),
          row.names = FALSE)
  }

  # 2. Balanced samples fill ----------------------------------------------------
  cat(sprintf("\n2. Balanced-sample construction (%d samples)\n", R_POOL))
  t0 <- proc.time()[["elapsed"]]
  pools <- vapply(seq_len(R_POOL), function(r) attr(gen_sample(s), "pool_size"), numeric(1))
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  pool units drawn per sample: min %d, median %d, max %d  (%.3f s per sample)\n",
              min(pools), as.integer(median(pools)), max(pools), el / R_POOL))
  report("every group reached n_g in every sample (checked inside gen_sample)", TRUE)

  # 3. tau_pop: closed form vs Monte Carlo ------------------------------------
  cat(sprintf("\n3. tau_pop: closed form vs Monte Carlo (N = %.0e)\n", N_POP))
  pop  <- gen_population(N_POP, s)
  D    <- pair_diff(pop$Ypo, pr)
  tp   <- true_tau_pop(s)
  mc   <- colMeans(D); mc_se <- apply(D, 2, sd) / sqrt(N_POP)
  z    <- (mc - tp$tau) / mc_se
  cat(sprintf("  group shares in population: %s\n",
              paste(sprintf("%.3f", tabulate(pop$T, s$A) / N_POP), collapse = " ")))
  cat(sprintf("  max |MC - closed form| = %.4f;  max |z| = %.2f\n",
              max(abs(mc - tp$tau)), max(abs(z))))
  report(sprintf("closed-form tau_pop within %g MC SEs for all pairs", Z_TOL),
         all(abs(z) < Z_TOL))
  rm(pop, D)

  # 4. tau_bal: p_t-weighted MC vs direct average over balanced samples -------
  cat(sprintf("\n4. tau_bal: p_t-weighted mu_bal (N = %.0e) vs %d balanced samples\n",
              N_BAL, R_BAL))
  mb <- mu_bal_mc(s, n = N_BAL)
  cat("  mu_pop: ", sprintf("%7.4f", mu_pop(s)), "\n")
  cat("  mu_bal: ", sprintf("%7.4f", mb$mu), "\n")
  cat("  (SE)    ", sprintf("%7.4f", mb$se), "\n")
  tb <- true_tau(s, mb$mu)
  tb_se <- true_tau_se(s, mb$vcov)   # MC SE of tau_bal, full covariance of mu_bal

  direct <- t(vapply(seq_len(R_BAL), function(r) {
    d <- gen_sample(s, keep_po = TRUE)
    colMeans(pair_diff(as.matrix(d[, paste0("Y", 1:s$A)]), pr))
  }, numeric(nrow(pr))))
  dir_mean <- colMeans(direct); dir_se <- apply(direct, 2, sd) / sqrt(R_BAL)
  zb <- (dir_mean - tb$tau) / sqrt(dir_se^2 + tb_se^2)
  cat(sprintf("  max |direct - weighted MC| = %.4f;  max |z| = %.2f\n",
              max(abs(dir_mean - tb$tau)), max(abs(zb))))
  report(sprintf("weighted-MC tau_bal within %g SEs of direct balanced-sample average", Z_TOL),
         all(abs(zb) < Z_TOL))

  # 5. Table and output ---------------------------------------------------------
  out <- data.frame(sim = sim, t = tp$t, tp = tp$tp,
                    tau_pop = tp$tau, tau_bal = tb$tau, tau_bal_mc_se = tb_se,
                    tau_bal_direct = dir_mean, tau_bal_direct_se = dir_se)
  out$bal_minus_pop <- out$tau_bal - out$tau_pop
  out$z_bal_vs_pop  <- out$bal_minus_pop / out$tau_bal_mc_se
  f <- sprintf("results/true_effects_sim%d.csv", sim)
  write.csv(out, f, row.names = FALSE)

  cat("\n5. True effects (mc_se = Monte Carlo SE of tau_bal; z = (tau_bal - tau_pop) / mc_se)\n")
  tab <- data.frame(t = out$t, tp = out$tp,
                    tau_pop = out$tau_pop,
                    tau_bal = round(out$tau_bal, 4),
                    mc_se   = round(out$tau_bal_mc_se, 4),
                    direct  = round(out$tau_bal_direct, 4),
                    bal_minus_pop = round(out$bal_minus_pop, 4),
                    z = round(out$z_bal_vs_pop, 1))
  print(tab, row.names = FALSE)
  cat(sprintf("  pairs with |z| > 3 (tau_bal clearly differs from tau_pop): %d of %d\n",
              sum(abs(out$z_bal_vs_pop) > 3), nrow(out)))
  cat(sprintf("  written to %s\n", f))
}

cat(sprintf("\nOverall: %s\n", if (all_ok) "ALL CHECKS PASSED" else "SOME CHECKS FAILED"))
cat(sprintf("R %s; MASS %s\n", getRversion(), packageVersion("MASS")))
