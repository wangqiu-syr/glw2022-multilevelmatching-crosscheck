# test_estimators.R -- Stage 2 tests for R/estimators.R. Run from repo root:
#   Rscript R/test_estimators.R          (unit tests + timing + reproducibility)
#   Rscript R/test_estimators.R quick    (unit tests only)

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

n_fail <- 0L
cur_test <- NA_character_
test_status <- c()
begin_test <- function(id, title) {
  cur_test <<- id
  test_status[id] <<- "PASS"
  cat(sprintf("\nTest %s. %s\n", id, title))
}
check <- function(label, ok) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (!isTRUE(ok)) { n_fail <<- n_fail + 1L; test_status[cur_test] <<- "FAIL" }
}
finish <- function() {
  cat("\n================ Test summary ================\n")
  for (id in names(test_status)) cat(sprintf("  Test %-3s %s\n", id, test_status[[id]]))
  cat(sprintf("  Failures: %d\n", n_fail))
  quit(status = if (n_fail > 0) 1 else 0)
}

# ---------------------------------------------------------------------------
begin_test("1", "Hand-made data with obvious nearest neighbours")

# 1a. CM: 8 well-separated clusters in 7-D, one unit per group (A = 3) in each,
# tiny within-cluster jitter. The NN in group s of any unit is the group-s unit
# of the same cluster. Y = 100 * cluster + group, so imputed Y is predictable.
set.seed(1)
A <- 3; K <- 8
centers <- matrix(rnorm(K * 7, sd = 100), K, 7)
grid <- expand.grid(s = 1:A, k = 1:K)
M <- centers[grid$k, ] + matrix(rnorm(nrow(grid) * 7, sd = 0.01), nrow(grid), 7)
colnames(M) <- COV_NAMES
d1 <- data.frame(M, T = grid$s, Y = 100 * grid$k + grid$s)
expected <- outer(100 * grid$k, 1:A, "+")
cm1 <- impute_cm(d1, A)
check("CM imputes the same-cluster unit's Y for every unit and level",
      identical(unname(cm1$Yimp), expected + 0))

D2 <- mahalanobis_d2(M)
ref <- mahalanobis(M, M[5, ], cov(M))
check("mahalanobis_d2 matches stats::mahalanobis (row 5)",
      isTRUE(all.equal(D2[5, ], unname(ref), tolerance = 1e-8)))

# 1b. GPS matching core with hand-set scores (binary fractions -> exact ties).
P <- cbind(c(0.125, 0.5, 0.875, 0.140625, 0.46875, 0.6875), NA)
P[, 2] <- 1 - P[, 1]
Tg <- c(1, 1, 1, 2, 2, 2); Yg <- c(10, 20, 30, 40, 50, 60)
g1 <- impute_from_scores(P, Tg, Yg, 2)
cat("  imputed matrix from hand-set scores:\n")
print(cbind(T = Tg, p1 = P[, 1], p2 = P[, 2], Y = Yg, g1$Yimp))
check("GPS column 1: units 4, 5 -> Y 10, 20; unit 6 is an exact tie (20 or 30)",
      all(g1$Yimp[4:5, 1] == c(10, 20)) && g1$Yimp[6, 1] %in% c(20, 30))
check("GPS column 2 (matched on p2): units 1, 2, 3 -> Y 40, 50, 60",
      all(g1$Yimp[1:3, 2] == c(40, 50, 60)))

# 1c. Random tie-breaking: two equidistant targets must each be chosen ~half
# the time, and no other target ever.
set.seed(2)
picks <- replicate(2000, impute_from_scores(P, Tg, Yg, 2)$Yimp[6, 1])
tab <- table(picks)
cat(sprintf("  unit 6 tie over 2000 draws: Y=20 %d times, Y=30 %d times\n",
            sum(picks == 20), sum(picks == 30)))
check("ties broken at random among the tied units only (binomial |z| < 4)",
      setequal(names(tab), c("20", "30")) &&
        abs(sum(picks == 20) - 1000) / sqrt(500) < 4)

# ---------------------------------------------------------------------------
begin_test("2", "Own-group units keep their observed Y")
set.seed(3)
for (sim in 1:2) {
  s <- sim_setup(sim)
  d <- gen_sample(s)
  idx <- unlist(lapply(split(seq_len(nrow(d)), d$T),
                       function(g) g[sample.int(length(g), replace = TRUE)]))
  for (lab in c("sample", "bootstrap sample")) {
    dd <- if (lab == "sample") d else d[idx, ]
    own <- cbind(seq_len(nrow(dd)), dd$T)
    cm  <- impute_cm(dd, s$A); gp <- impute_gps(dd, s$A)
    check(sprintf("Sim %d %s: CM and GPS keep observed Y on own level; no NA",
                  sim, lab),
          all(cm$Yimp[own] == dd$Y) && all(gp$Yimp[own] == dd$Y) &&
            !anyNA(cm$Yimp) && !anyNA(gp$Yimp))
  }
}

# ---------------------------------------------------------------------------
begin_test("3", "DIF equals plain group-mean differences")
for (sim in 1:2) {
  s <- sim_setup(sim)
  d <- gen_sample(s)
  gm <- tapply(d$Y, d$T, mean)
  manual <- unlist(lapply(seq_len(nrow(pairs_of(s$A))), function(k) {
    p <- pairs_of(s$A)[k, ]; gm[[p[1]]] - gm[[p[2]]] }))
  check(sprintf("Sim %d: tau_dif equals tapply mean differences (%d pairs)",
                sim, length(manual)),
        isTRUE(all.equal(unname(tau_dif(d, s$A)), manual)))
}

# ---------------------------------------------------------------------------
begin_test("4", "GPS probabilities")
for (sim in 1:2) {
  s <- sim_setup(sim)
  d <- gen_sample(s)
  gf <- fit_gps(d, s$A)
  check(sprintf("Sim %d: multinom converged; rows of p_hat sum to 1 (max err %.1e)",
                sim, max(abs(rowSums(gf$P) - 1))),
        gf$converged && max(abs(rowSums(gf$P) - 1)) < 1e-10)
  if (sim == 2) {
    P <- gf$P
    adj <- function(a, b) median(abs(P[, a] - P[, b]))
    cat(sprintf("  median |p_hat_6 - p_hat_7| = %.4f   (cor = %.4f)\n", adj(6, 7), cor(P[, 6], P[, 7])))
    cat(sprintf("  for comparison: median |p_hat_5 - p_hat_6| = %.4f, |p_hat_7 - p_hat_8| = %.4f\n",
                adj(5, 6), adj(7, 8)))
    cf <- coef(gf$fit)
    cat("  multinom coefficients (vs t1) for t6 and t7:\n")
    print(round(cf[c("6", "7"), ], 3))
    Ptrue <- treatment_probs(as.matrix(d[, COV_NAMES]), s)
    check(sprintf("true p_6 == p_7 for every unit, up to rounding (max diff %.1e) (D10)",
                  max(abs(Ptrue[, 6] - Ptrue[, 7]))),
          max(abs(Ptrue[, 6] - Ptrue[, 7])) < 1e-12)
    cat(sprintf("  estimation noise in this sample: median |p_hat_6 - p_6| = %.4f, |p_hat_7 - p_7| = %.4f\n",
                median(abs(P[, 6] - Ptrue[, 6])), median(abs(P[, 7] - Ptrue[, 7]))))

    # Over 40 samples: p_hat_6 - p_hat_7 should average to 0 (only noise),
    # whereas p_hat_7 - p_hat_8 should average to the true gap. In a balanced
    # sample multinom estimates the design probabilities
    #   q_t(x) = (p_t(x) / pi_t) / sum_r (p_r(x) / pi_r),  pi_t = E[p_t],
    # not p_t(x) itself; q_6 = q_7 still holds because pi_6 = pi_7.
    set.seed(4)
    pi_t <- colMeans(treatment_probs(gen_covariates(2e5), s))
    Vref <- gen_covariates(2000)
    Pref <- sweep(treatment_probs(Vref, s), 2, pi_t, "/")
    Pref <- Pref / rowSums(Pref)
    ref_df <- as.data.frame(Vref)
    gaps <- t(replicate(40, {
      fit <- fit_gps(gen_sample(s), s$A)$fit
      Ph <- predict(fit, newdata = ref_df, type = "probs")
      c(g67 = mean(Ph[, 6] - Ph[, 7]), g78 = mean(Ph[, 7] - Ph[, 8]))
    }))
    z67 <- mean(gaps[, "g67"]) / (sd(gaps[, "g67"]) / sqrt(40))
    cat(sprintf("  40 samples, mean over 2000 reference units:\n"))
    cat(sprintf("    p_hat_6 - p_hat_7: mean %+.4f (SE %.4f, z = %.2f); true gap %+.4f\n",
                mean(gaps[, "g67"]), sd(gaps[, "g67"]) / sqrt(40), z67, 0))
    cat(sprintf("    p_hat_7 - p_hat_8: mean %+.4f (SE %.4f);            true design gap %+.4f\n",
                mean(gaps[, "g78"]), sd(gaps[, "g78"]) / sqrt(40),
                mean(Pref[, 7] - Pref[, 8])))
    cat(sprintf("    (true design gap q_6 - q_7 = %+.1e)\n", mean(Pref[, 6] - Pref[, 7])))
    check("p_hat_6 - p_hat_7 averages to 0 across samples (|z| < 4): gap is estimation noise",
          abs(z67) < 4)
  }
}

# ---------------------------------------------------------------------------
begin_test("4b", sprintf("GPS fit is start-independent at reltol = %g (cold vs warm start, Sim 2)", GPS_RELTOL))
# The pipeline always cold-starts. Here a warm-started fit (from the original
# sample's weights, same reltol/maxit) must give the same p_hat and tau-hat.
fit_gps_warm <- function(d, A, start) {
  dd <- d; dd$T_f <- factor(d$T, levels = seq_len(A))
  fit <- nnet::multinom(gps_formula, data = dd, trace = FALSE, maxit = GPS_MAXIT,
                        reltol = GPS_RELTOL, Wts = start)
  list(P = unname(fitted(fit)), converged = fit$convergence == 0, fit = fit)
}
set.seed(6)
s <- sim_setup(2); d <- gen_sample(s)
orig <- fit_gps(d, s$A)
groups <- split(seq_len(nrow(d)), d$T)
n_bs <- 30
cmp <- t(vapply(seq_len(n_bs), function(b) {
  idx <- unlist(lapply(groups, function(g) g[sample.int(length(g), replace = TRUE)]))
  db <- d[idx, ]
  tc <- system.time(fc <- fit_gps(db, s$A))[["elapsed"]]
  tw <- system.time(fw <- fit_gps_warm(db, s$A, orig$fit$wts))[["elapsed"]]
  set.seed(1000 + b); ic <- impute_from_scores(fc$P, db$T, db$Y, s$A)
  set.seed(1000 + b); iw <- impute_from_scores(fw$P, db$T, db$Y, s$A)
  dtau <- tau_from_imputed(ic$Yimp) - tau_from_imputed(iw$Yimp)
  c(max_dP = max(abs(fc$P - fw$P)),
    rel_dev = abs(fc$fit$deviance - fw$fit$deviance) / fc$fit$deviance,
    max_dtau = max(abs(dtau)), n_tau_diff = sum(dtau != 0),
    conv = fc$converged && fw$converged, t_cold = tc, t_warm = tw)
}, numeric(7)))
cat(sprintf("  %d bootstrap samples; max |p_hat cold - warm| = %.2e; max relative deviance diff = %.2e\n",
            n_bs, max(cmp[, "max_dP"]), max(cmp[, "rel_dev"])))
cat(sprintf("  GPS tau-hat (36 pairs x %d samples): %d of %d differ; max |diff| = %.2e\n",
            n_bs, sum(cmp[, "n_tau_diff"]), 36 * n_bs, max(cmp[, "max_dtau"])))
cat(sprintf("  multinom time per fit: cold %.3f s, warm %.3f s\n",
            mean(cmp[, "t_cold"]), mean(cmp[, "t_warm"])))
check("all cold and warm fits converged", all(cmp[, "conv"] == 1))
check("same optimum: max |p_hat diff| < 1e-6 and relative deviance diff < 1e-10",
      max(cmp[, "max_dP"]) < 1e-6 && max(cmp[, "rel_dev"]) < 1e-10)
check("GPS tau-hat identical for all pairs and samples", sum(cmp[, "n_tau_diff"]) == 0)

if (quick) finish()

# ---------------------------------------------------------------------------
begin_test("5", "Timing: 20 replications with B = 20 (serial)")
R_T <- 20; B_T <- 20; R_FULL <- 500; B_FULL <- 200; WORKERS <- 10
truth <- list(read.csv("results/true_effects_sim1.csv"),
              read.csv("results/true_effects_sim2.csv"))
timing <- list()
for (sim in 1:2) {
  s <- sim_setup(sim)
  seeds <- make_seeds(R_T, master_seed = 1000 + sim)
  t0 <- proc.time()[["elapsed"]]
  res <- lapply(seq_len(R_T), function(r) run_replication(r, s, B_T, seeds[[r]]))
  el <- proc.time()[["elapsed"]] - t0
  saveRDS(res, sprintf("results/test_sim%d_R%d_B%d.rds", sim, R_T, B_T))
  per_fit <- el / (R_T * (B_T + 1))
  full_serial <- per_fit * R_FULL * (B_FULL + 1)
  timing[[sim]] <- c(elapsed = el, per_rep = el / R_T, per_fit = per_fit,
                     full_serial_h = full_serial / 3600)
  cat(sprintf("  Sim %d: %.1f s total, %.2f s per replication, %.4f s per full estimate (DIF+CM+GPS)\n",
              sim, el, el / R_T, per_fit))
  cat(sprintf("         extrapolated 500 reps x B = 200: %.1f h serial, ~%.1f h on %d workers (assuming 80%% efficiency)\n",
              full_serial / 3600, full_serial / 3600 / (WORKERS * 0.8), WORKERS))
  check(sprintf("Sim %d: all %d replications completed", sim, R_T), length(res) == R_T)

  # Smoke summary (20 reps, B = 20: noisy, not a result).
  est <- do.call(rbind, lapply(res, `[[`, "est"))
  tr <- truth[[sim]][, c("t", "tp", "tau_pop", "tau_bal")]
  est <- merge(est, tr, by = c("t", "tp"))
  sm <- do.call(rbind, lapply(split(est, est$estimator), function(e) {
    pb <- tapply(e$est - e$tau_pop, paste(e$t, e$tp), mean)
    bb <- tapply(e$est - e$tau_bal, paste(e$t, e$tp), mean)
    data.frame(estimator = e$estimator[1],
               mean_abs_bias_pop = mean(abs(pb)), mean_abs_bias_bal = mean(abs(bb)),
               cover_pop = mean(e$lo <= e$tau_pop & e$tau_pop <= e$hi),
               cover_bal = mean(e$lo <= e$tau_bal & e$tau_bal <= e$hi))
  }))
  cat("  smoke summary (averaged over pairs; noisy):\n")
  print(format(sm, digits = 3), row.names = FALSE)
  dg <- lapply(res, `[[`, "diag")
  cat(sprintf("  diagnostics: GPS converged in %d/%d samples; bootstrap non-converged fits: %d; ",
              sum(sapply(dg, `[[`, "gps_converged")), R_T, sum(sapply(dg, `[[`, "boot_gps_nonconverged"))))
  cat(sprintf("share of units with any p_hat < 0.01: mean %.3f\n",
              mean(sapply(dg, `[[`, "p_hat_any_below_01"))))
  cat("  mean CM matching distance by target group: ",
      sprintf("%.2f", rowMeans(sapply(dg, `[[`, "cm_match_dist"))), "\n")
  cat("  mean GPS matching distance by target group:",
      sprintf("%.4f", rowMeans(sapply(dg, `[[`, "gps_match_dist"))), "\n")
}

# ---------------------------------------------------------------------------
begin_test("6", "Reproducibility: serial vs parallel (PSOCK) give identical results")
s <- sim_setup(1); seeds <- make_seeds(4, master_seed = 99)
ser <- lapply(1:4, function(r) run_replication(r, s, 5, seeds[[r]])$est)
library(parallel)
cl <- makeCluster(4)
root <- getwd()
invisible(clusterCall(cl, function(p) { setwd(p); source("R/dgp.R"); source("R/estimators.R");
                                        source("R/replicate.R"); NULL }, root))
par <- parLapplyLB(cl, 4:1, function(r, s, seeds) run_replication(r, s, 5, seeds[[r]])$est,
                   s = s, seeds = seeds)
stopCluster(cl)
par <- rev(par)
check("4 replications (run in reverse order on 4 workers) identical to serial",
      identical(ser, par))
cat(sprintf("  cores detected: %d\n", detectCores()))

finish()
