# 07_d8_analytic_intervals.R -- D8 sensitivity: analytic (Abadie-Imbens)
# variances for the rebuild's own GPS estimates, every replication.
#
# The main run's intervals are percentile bootstrap intervals (D8). This script
# provides the analytic alternative for the GPS:
#   - AI 2006: treats the GPS as known (also available from the "existing" fit);
#   - AI 2016: accounts for estimating the GPS by multinomial logit.
# multiMatch() computes the AI 2016 variance only for the package's own GPS fit
# (nnet defaults). Here the package's internal steps are run with OUR fitted
# model instead (reltol = 1e-14; its probabilities and coefficient covariance),
# so the variance belongs to the rebuild's estimator. The CM analytic variance
# needs no new run: the package's "covariates" fit reproduces our CM exactly, and
# its AI 2006 variance is in crosscheck/results/sim<k>_crosscheck_reps.csv.gz.
#
# Check built in: run once per simulation with the package's own fit in place of
# ours, the same steps must reproduce multiMatch(match_on = "multinom") exactly.
#
#   R_LIBS=~/R/mlm-lib Rscript crosscheck/R/07_d8_analytic_intervals.R sim=1 [R=500] [workers=4]
# Run from the repo root. Writes crosscheck/results/sim<k>_d8_gps_variances.csv.gz;
# per-replication checkpoints go to crosscheck/reps/d8_sim<k>/ (git-ignored).

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
suppressPackageStartupMessages(library(multilevelMatching))
# Per-replication results are written gzip-compressed; read.csv() reads them as they are.
write_csv_gz <- function(x, f) { con <- gzfile(f, "w"); on.exit(close(con)); write.csv(x, con, row.names = FALSE) }

opt <- list(sim = "1", R = "500", workers = "4", rep_root = file.path("crosscheck", "reps"))
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2] else stop("unknown argument: ", a)
}
SIM <- as.integer(opt$sim); R_REPS <- as.integer(opt$R); WORKERS <- as.integer(opt$workers)
rep_dir <- file.path(opt$rep_root, paste0("d8_sim", SIM))
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...))); flush(stdout()) }

# AI 2006 and AI 2016 variances from the package's internal steps, for a given
# fitted multinomial model (probabilities P, coefficient covariance V).
ai_variances <- function(Y, W, X, P, V) {
  Xc <- X; oldClass(Xc) <- "matrix"
  pd <- suppressMessages(multilevelMatching:::prepareData(
    Y = Y, W = W, X = Xc, match_on = "multinom", trimming = FALSE,
    model_options = list(reference_level = sort(W)[1]), M_matches = 1, J_var_matches = 1))
  X_covars <- pd$X
  pd <- append(pd, list(M_matches = 1, J_var_matches = 1))
  pd$X <- P
  set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  me <- do.call(multilevelMatching:::matchAllTreatments, pd)
  args <- append(pd, me)
  tau <- do.call(multilevelMatching:::estimateTau, args)$tau_dfm
  a16 <- append(args, list(prop_score_ests = P, vcov_coeff = V))
  a16$tau_dfm <- tau; a16$X_covars <- X_covars
  v16 <- do.call(multilevelMatching:::calcSigSqAI2016, a16)$est_sigsq_AI2016
  data.frame(t = tau$Trt1, tp = tau$Trt2, est_pkg_rule = -tau$Estimate,
             var_AI2006 = tau$Variance, var_AI2016 = v16)
}

one_rep <- function(r, sim, seeds) {
  setup <- sim_setup(sim); A <- setup$A
  assign(".Random.seed", seeds[[r]], envir = globalenv())
  d <- gen_sample(setup)
  dd <- d; dd$T_f <- factor(d$T, levels = seq_len(A))
  fit <- nnet::multinom(gps_formula, data = dd, trace = FALSE, maxit = GPS_MAXIT,
                        reltol = GPS_RELTOL, Hess = TRUE)
  P <- unname(fitted(fit)); colnames(P) <- paste0("p", 1:A)
  P_main <- fit_gps(d, A)$P                        # the main run's GPS fit (no Hessian)
  same_P <- isTRUE(all.equal(P, P_main, tolerance = 0, check.attributes = FALSE))
  X <- as.matrix(d[, COV_NAMES])
  out <- ai_variances(d$Y, d$T, X, P, vcov(fit))
  out$GPS <- unname(tau_from_imputed(impute_from_scores(P, d$T, d$Y, A)$Yimp))
  out$rep <- r; out$same_P_as_main_fit <- same_P; out$converged <- fit$convergence == 0
  out
}

seeds <- make_seeds(R_REPS, 20220000 + SIM)

# ---- check: with the package's own fit, the same steps reproduce multiMatch ----
local({
  setup <- sim_setup(SIM); A <- setup$A
  assign(".Random.seed", seeds[[1]], envir = globalenv())
  d <- gen_sample(setup); X <- as.matrix(d[, COV_NAMES])
  Xc <- X; oldClass(Xc) <- "matrix"
  set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  mm <- suppressMessages(multiMatch(d$Y, d$T, Xc, trimming = FALSE, match_on = "multinom"))
  P_pkg <- unname(mm$propensity_scores)
  mine <- ai_variances(d$Y, d$T, X, P_pkg, vcov(mm$model))
  ok <- identical(mine$var_AI2016, mm$results$VarianceAI2016) && identical(mine$var_AI2006, mm$results$Variance)
  log_msg("Check (Simulation %d, replication 1): internal steps with the package's own fit reproduce multiMatch's AI 2006 and AI 2016 variances exactly: %s", SIM, ok)
  if (!ok) stop("check failed: max diff AI2016 ", max(abs(mine$var_AI2016 - mm$results$VarianceAI2016)))
})

rep_file <- function(r) file.path(rep_dir, sprintf("rep_%03d.rds", r))
todo <- which(!file.exists(vapply(seq_len(R_REPS), rep_file, character(1))))
log_msg("Simulation %d: %d replications to run on %d workers", SIM, length(todo), WORKERS)
if (length(todo)) {
  cl <- parallel::makeCluster(WORKERS)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  invisible(parallel::clusterCall(cl, function(p) {
    setwd(p); source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
    suppressPackageStartupMessages(library(multilevelMatching)); NULL }, getwd()))
  parallel::clusterExport(cl, c("one_rep", "ai_variances"))
  t0 <- proc.time()[["elapsed"]]
  for (ch in split(todo, ceiling(seq_along(todo) / (25 * WORKERS)))) {
    ok <- parallel::parLapplyLB(cl, ch, function(r, sim, seeds, rep_dir) {
      res <- tryCatch(one_rep(r, sim, seeds), error = function(e) e)
      f <- file.path(rep_dir, sprintf("rep_%03d.rds", r))
      if (inherits(res, "error")) { writeLines(conditionMessage(res), sub("\\.rds$", "_error.txt", f)); return(FALSE) }
      saveRDS(res, f); TRUE
    }, sim = SIM, seeds = seeds, rep_dir = rep_dir)
    log_msg("Simulation %d: %d/%d done, elapsed %.0fs", SIM,
            sum(file.exists(vapply(seq_len(R_REPS), rep_file, character(1)))), R_REPS, proc.time()[["elapsed"]] - t0)
  }
}
files <- vapply(seq_len(R_REPS), rep_file, character(1))
if (all(file.exists(files))) {
  out <- do.call(rbind, lapply(files, readRDS))
  out <- out[c("rep", "t", "tp", "GPS", "est_pkg_rule", "var_AI2006", "var_AI2016", "same_P_as_main_fit", "converged")]
  num <- vapply(out, is.double, logical(1)); out[num] <- lapply(out[num], signif, digits = 12)
  write_csv_gz(out, file.path("crosscheck", "results", sprintf("sim%d_d8_gps_variances.csv.gz", SIM)))
  log_msg("Wrote sim%d_d8_gps_variances.csv.gz (%d rows); fits identical to the main GPS fit in %d of %d replications",
          SIM, nrow(out), sum(tapply(out$same_P_as_main_fit, out$rep, all)), R_REPS)
}
