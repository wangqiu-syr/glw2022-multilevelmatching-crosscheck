# 12_cm_biasadjust.R -- covariate matching with Matching::Match(BiasAdjust = TRUE),
# on the rebuild's samples (crosscheck/rebuild2/README.md, section 5).
#
# Question: does a bias-adjusted CM reproduce the pattern of the paper's CM
# biases? The paper's Simulation 2 CM biases fit 1.38 x the sample-mix term +
# 0.22 x the matching-gap term of the rebuild's CM (R^2 = 0.93), as if most of
# the matching gap had been removed.
#
# For each level s, multilevelMatching calls Matching::Match with the units
# outside group s as "treated" and group s as the pool (estimand ATT, M = 1,
# Mahalanobis on the full sample, Weight = 2). Then the mean of the imputed
# Y(s) over all n units is
#   mu_hat(s) = mean(Y) - ((n - n_s) / n) * ATT_hat_s,
# because ATT_hat_s = mean over i outside s of (Y_i - Yhat_i(s)). With
# BiasAdjust = FALSE this is the rebuild's CM (checked below). With
# BiasAdjust = TRUE, Match adjusts each matched outcome by a regression on the
# covariates of the matched pool units (Abadie and Imbens, 2011); same matches.
#
#   Rscript crosscheck/R/12_cm_biasadjust.R sim=2 [R=500] [workers=4]
# Run from the repo root. Needs Matching. Writes
# crosscheck/results/sim<k>_cm_biasadjust_reps.csv.

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
suppressPackageStartupMessages(library(Matching))

opt <- list(sim = "2", R = "500", workers = "4")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2] else stop("unknown argument: ", a)
}
SIM <- as.integer(opt$sim); R_REPS <- as.integer(opt$R); WORKERS <- as.integer(opt$workers)
log_msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...))); flush(stdout()) }

# Match() only adjusts with ties = TRUE (it otherwise warns and silently sets
# BiasAdjust = FALSE), so the adjusted call uses ties = TRUE and returns both its
# adjusted estimate and, from the same matches, its unadjusted one (est.noadj).
# Warnings are turned into errors so that a silent switch-off cannot recur.
mu_match <- function(Y, T, X, s, bias_adjust) {
  m <- withCallingHandlers(
    Matching::Match(Y = Y, Tr = as.integer(T != s), X = X, estimand = "ATT", M = 1,
                    Weight = 2, distance.tolerance = 0, ties = bias_adjust,
                    BiasAdjust = bias_adjust),
    warning = function(w) stop("Match() warned: ", conditionMessage(w)))
  share <- mean(T != s)
  c(adj = mean(Y) - share * as.numeric(m$est)[1], noadj = mean(Y) - share * as.numeric(m$est.noadj)[1])
}

one_rep <- function(r, sim, seeds) {
  setup <- sim_setup(sim); A <- setup$A
  assign(".Random.seed", seeds[[r]], envir = globalenv())
  d <- gen_sample(setup)
  cm <- tau_from_imputed(impute_cm(d, A)$Yimp)        # the rebuild's CM
  X <- as.matrix(d[, COV_NAMES])
  set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  mu0 <- vapply(seq_len(A), function(s) mu_match(d$Y, d$T, X, s, FALSE)[["adj"]], numeric(1))
  m1 <- vapply(seq_len(A), function(s) mu_match(d$Y, d$T, X, s, TRUE), numeric(2))
  pr <- pairs_of(A)
  data.frame(rep = r, t = pr[, 1], tp = pr[, 2], CM = unname(cm),
             CM_Match = mu0[pr[, 1]] - mu0[pr[, 2]],                       # package settings, ties = FALSE
             CM_Match_ties = m1["noadj", pr[, 1]] - m1["noadj", pr[, 2]],  # ties = TRUE, unadjusted
             CM_BiasAdjust = m1["adj", pr[, 1]] - m1["adj", pr[, 2]])      # ties = TRUE, adjusted
}

seeds <- make_seeds(R_REPS, 20220000 + SIM)
t0 <- proc.time()[["elapsed"]]
cl <- parallel::makeCluster(WORKERS)
invisible(parallel::clusterCall(cl, function(p) {
  setwd(p); source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
  suppressPackageStartupMessages(library(Matching)); NULL }, getwd()))
parallel::clusterExport(cl, "mu_match")
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(R_REPS), one_rep, sim = SIM, seeds = seeds))
parallel::stopCluster(cl)
log_msg("Simulation %d: %d replications in %.0fs (Matching %s)", SIM, R_REPS,
        proc.time()[["elapsed"]] - t0, as.character(packageVersion("Matching")))
dif <- abs(res$CM_Match - res$CM)
log_msg("Check: Match() without bias adjustment vs the rebuild's CM: max |difference| = %.2e; cells differing by > 1e-10: %d of %d",
        max(dif), sum(dif > 1e-10), nrow(res))
dif2 <- abs(res$CM_Match_ties - res$CM)
log_msg("Check: the same with ties = TRUE (the adjusted call's own unadjusted estimate): max |difference| = %.2e; cells differing by > 1e-10: %d of %d",
        max(dif2), sum(dif2 > 1e-10), nrow(res))
log_msg("Bias adjustment: mean |adjusted - unadjusted| over all cells = %.4f", mean(abs(res$CM_BiasAdjust - res$CM_Match_ties)))
num <- vapply(res, is.double, logical(1)); res[num] <- lapply(res[num], signif, digits = 12)
f <- file.path("crosscheck", "results", sprintf("sim%d_cm_biasadjust_reps.csv", SIM))
write.csv(res, f, row.names = FALSE)
log_msg("Wrote %s (%d rows)", f, nrow(res))
