# 09_cm_bias_decomposition.R -- where the CM bias comes from, replication by
# replication (docs/intervals.md, section 4).
#
# The outcomes are linear without intercept: Y_i(s) = beta_s' w_i + eps_i, with
# w = the 7 covariates minus the design's shift (R/dgp.R). CM imputes Y_i(s) by
# the observed Y of the nearest unit j(i, s) in group s (Mahalanobis distance on
# the 7 covariates, with replacement), and keeps Y_i for i's own group. Then,
# exactly, in every sample,
#
#   tau_hat(t, t') = (beta_t - beta_t')' m          sample-mix term
#                  + beta_t' D_t - beta_t'' D_t'     matching-gap term
#                  + N_t - N_t'                      noise term
#
#   m   = (1/n) sum_i w_i                     (covariate mean of the sample)
#   D_s = (1/n) sum_{i: T_i != s} (w_j(i,s) - w_i)   (mean matching gap, 7-vector)
#   N_s = (1/n) sum_i eps_j(i,s)              (j(i,s) = i for i in group s)
#
# The sample-mix term averages to tau_bal over replications; the noise term to 0
# (matching uses covariates only). So bias against tau_pop = (tau_bal - tau_pop)
# + E[matching-gap term] + Monte Carlo noise, and the matching-gap term splits
# by level (t or t') and by covariate.
#
#   Rscript crosscheck/R/09_cm_bias_decomposition.R [sim=2] [R=500] [workers=4]
# Run from the repo root. Writes crosscheck/results/sim<k>_cm_gap_reps.csv (per
# replication: D_s, m, N_s, tau_hat) and checks that tau_hat equals the main
# run's CM estimates (mean over replications) and the decomposition identity.
#
# Diagnostic (not an estimator of the rebuild): a regression-adjusted CM, in the
# spirit of the bias-corrected matching estimator of Abadie and Imbens (2011).
# For i not in group s it imputes Y_j(i,s) + mu_hat_s(x_i) - mu_hat_s(x_j(i,s)),
# with mu_hat_s an OLS fit of Y on the 7 covariates within group s, which removes
# the matching gap (the outcome model here is linear, so the adjustment is
# correctly specified: a favourable case). Same matches as CM. Written to
# crosscheck/results/sim<k>_cm_bc_reps.csv.

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")

opt <- list(sim = "2", R = "500", workers = "4")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2] else stop("unknown argument: ", a)
}
SIM <- as.integer(opt$sim); R_REPS <- as.integer(opt$R); WORKERS <- as.integer(opt$workers)
log_msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...))); flush(stdout()) }

one_rep <- function(r, sim, seeds) {
  setup <- sim_setup(sim); A <- setup$A
  assign(".Random.seed", seeds[[r]], envir = globalenv())
  d <- gen_sample(setup)                     # same stream use as the main run
  n <- nrow(d); T <- d$T
  W <- as.matrix(d[, COV_NAMES])
  w <- sweep(W, 2, setup$w_shift)            # outcome regressors
  eps <- d$Y - rowSums(w * setup$beta[T, ])  # eps_i = Y_i - beta_{T_i}' w_i
  Dm <- sqrt(mahalanobis_d2(W))              # as impute_cm(); first RNG use after gen_sample
  gap <- matrix(0, A, 7, dimnames = list(paste0("t", 1:A), COV_NAMES))
  N <- numeric(A); Yimp <- Ybc <- matrix(NA_real_, n, A)
  sig <- matrix(NA_real_, A, 4, dimnames = list(NULL, c("sigma2_hat", "sigma2_covariate_part", "sigma2_noise_part", "sigma2_cross")))
  K <- numeric(n); sig2_i <- numeric(n)       # times used as a match; sigma2_hat per unit
  n_tie <- 0L
  for (s in seq_len(A)) {
    g <- which(T == s)
    nn <- nn_pick(Dm[, g, drop = FALSE])
    j <- g[nn$j]; own <- T == s
    j[own] <- which(own)                     # own group keeps its observed Y
    n_tie <- n_tie + sum(nn$n_tie[!own] > 1)
    Yimp[, s] <- d$Y[j]
    cf <- lm.fit(cbind(1, W[g, , drop = FALSE]), d$Y[g])$coefficients
    mu_s <- as.vector(cbind(1, W) %*% cf)
    Ybc[, s] <- Yimp[, s] + mu_s - mu_s[j]   # zero adjustment for own group (j = i)
    gap[s, ] <- colSums(W[j[!own], , drop = FALSE] - W[!own, , drop = FALSE]) / n
    N[s] <- mean(eps[j])
    # AI 2006 conditional variance, as in multilevelMatching (J = 1): match each
    # unit of group s to its nearest other unit l in group s (Mahalanobis, with
    # group s's own covariance, as Matching::Match does on the group's data), and
    # sigma2_i = (Y_i - Y_l)^2 / 2. Split it into the part from mu_s(x_i) - mu_s(x_l)
    # (covariate gap), the noise part (eps_i - eps_l)^2 / 2, and the cross term.
    Wg <- W[g, , drop = FALSE]; Dg <- mahalanobis_d2(Wg); diag(Dg) <- Inf
    l <- apply(Dg, 1, which.min)
    dmu <- as.vector((w[g, , drop = FALSE] - w[g[l], , drop = FALSE]) %*% setup$beta[s, ])
    deps <- eps[g] - eps[g[l]]
    sig[s, ] <- c(mean((d$Y[g] - d$Y[g[l]])^2) / 2, mean(dmu^2) / 2, mean(deps^2) / 2, mean(dmu * deps))
    sig2_i[g] <- (d$Y[g] - d$Y[g[l]])^2 / 2
    K <- K + tabulate(j[!own], nbins = n)
  }
  # AI 2006 variance as in multilevelMatching::estVarAI2006 (M = 1), per pair:
  #   (1/n) [ mean_i (Yhat_i(t) - Yhat_i(t') - tau_hat)^2                  "V_tau"
  #         + mean_i (K_i^2 + K_i) sigma2_hat_i 1{T_i in {t, t'}} ]        "V_E"
  # and V_E with the true sigma^2 = 1 in place of sigma2_hat ("V_E_true").
  m <- colMeans(w)
  pr <- pairs_of(A)
  tau_hat <- colMeans(Yimp)[pr[, 1]] - colMeans(Yimp)[pr[, 2]]
  mix   <- as.vector((setup$beta[pr[, 1], ] - setup$beta[pr[, 2], ]) %*% m)
  gapt  <- rowSums(setup$beta * gap)         # beta_s' D_s, per level s
  recon <- mix + gapt[pr[, 1]] - gapt[pr[, 2]] + N[pr[, 1]] - N[pr[, 2]]
  tau_bc <- colMeans(Ybc)[pr[, 1]] - colMeans(Ybc)[pr[, 2]]
  kf <- K^2 + K
  terms <- data.frame(rep = r, t = pr[, 1], tp = pr[, 2], est = unname(tau_hat), mix = mix,
                      gap_term = gapt[pr[, 1]] - gapt[pr[, 2]], noise_term = N[pr[, 1]] - N[pr[, 2]])
  terms$V_tau <- vapply(seq_len(nrow(pr)), function(k) mean((Yimp[, pr[k, 1]] - Yimp[, pr[k, 2]] - tau_hat[k])^2), numeric(1)) / n
  terms$V_E <- vapply(seq_len(nrow(pr)), function(k) mean(kf * sig2_i * (T %in% pr[k, ])), numeric(1)) / n
  terms$V_E_true <- vapply(seq_len(nrow(pr)), function(k) mean(kf * (T %in% pr[k, ])), numeric(1)) / n
  list(rep = r, gap = gap, m = m, N = N, sig = sig, tau_hat = tau_hat, tau_bc = tau_bc, terms = terms,
       max_identity_err = max(abs(recon - tau_hat)), n_tie = n_tie)
}

seeds <- make_seeds(R_REPS, 20220000 + SIM)
t0 <- proc.time()[["elapsed"]]
cl <- parallel::makeCluster(WORKERS)
invisible(parallel::clusterCall(cl, function(p) {
  setwd(p); source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R"); NULL }, getwd()))
res <- parallel::parLapplyLB(cl, seq_len(R_REPS), one_rep, sim = SIM, seeds = seeds)
parallel::stopCluster(cl)
log_msg("Simulation %d: %d replications in %.0fs", SIM, R_REPS, proc.time()[["elapsed"]] - t0)

setup <- sim_setup(SIM); A <- setup$A; pr <- pairs_of(A)
log_msg("Decomposition identity: max |reconstructed - tau_hat| over all replications and pairs = %.2e",
        max(vapply(res, `[[`, numeric(1), "max_identity_err")))
log_msg("Ties among CM nearest neighbours (all replications): %d", sum(vapply(res, `[[`, integer(1), "n_tie")))

# Check against the main run's CM estimates (mean over replications).
summ <- read.csv(file.path("results", sprintf("sim%d_summary.csv", SIM)))
cm_main <- summ[summ$estimator == "CM" & summ$target == "pop", ]
tau_mean <- Reduce(`+`, lapply(res, `[[`, "tau_hat")) / R_REPS
key <- paste(pr[, 1], pr[, 2]); mk <- match(key, paste(cm_main$t, cm_main$tp))
log_msg("Mean CM estimate vs the main run's (results/sim%d_summary.csv): max |difference| = %.2e",
        SIM, max(abs(tau_mean - cm_main$mean_est[mk])))

# Per-replication output, long by level s.
out <- do.call(rbind, lapply(res, function(x) {
  data.frame(rep = x$rep, s = seq_len(A), x$gap, N = x$N,
             matrix(x$m, A, 7, byrow = TRUE, dimnames = list(NULL, paste0("m_", COV_NAMES))),
             x$sig, check.names = FALSE)
}))
names(out)[3:9] <- paste0("gap_", COV_NAMES)
num <- vapply(out, is.double, logical(1)); out[num] <- lapply(out[num], signif, digits = 12)
f <- file.path("crosscheck", "results", sprintf("sim%d_cm_gap_reps.csv", SIM))
write.csv(out, f, row.names = FALSE)
log_msg("Wrote %s (%d rows)", f, nrow(out))

bc <- do.call(rbind, lapply(res, function(x)
  data.frame(rep = x$rep, t = pr[, 1], tp = pr[, 2], CM = unname(x$tau_hat), CM_regadj = unname(x$tau_bc))))
num <- vapply(bc, is.double, logical(1)); bc[num] <- lapply(bc[num], signif, digits = 12)
f <- file.path("crosscheck", "results", sprintf("sim%d_cm_bc_reps.csv", SIM))
write.csv(bc, f, row.names = FALSE)
log_msg("Wrote %s (%d rows)", f, nrow(bc))

# Per replication and pair: the three terms of tau_hat and the AI 2006 variance
# parts; check the variance against the package's (crosscheck/results, 04).
tm <- do.call(rbind, lapply(res, `[[`, "terms"))
tm$var_AI2006 <- tm$V_tau + tm$V_E
cx <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_crosscheck_reps.csv", SIM)))
k <- match(paste(tm$rep, tm$t, tm$tp), paste(cx$rep, cx$t, cx$tp))
rel <- abs(tm$var_AI2006 / cx$var_pkg_CM[k] - 1)
log_msg("AI 2006 variance rebuilt here vs the package's (04): max relative difference %.2e; %d of %d (replication, pair) cells differ by more than 1e-8",
        max(rel), sum(rel > 1e-8), length(rel))
num <- vapply(tm, is.double, logical(1)); tm[num] <- lapply(tm[num], signif, digits = 12)
f <- file.path("crosscheck", "results", sprintf("sim%d_cm_terms_reps.csv", SIM))
write.csv(tm, f, row.names = FALSE)
log_msg("Wrote %s (%d rows)", f, nrow(tm))
