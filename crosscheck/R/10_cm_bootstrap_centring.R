# 10_cm_bootstrap_centring.R -- where the percentile bootstrap interval for CM
# is centred, relative to the estimate (docs/intervals.md, section 4.4).
#
# The main run's percentile intervals were not kept per replication (only their
# coverage and mean width, in results/sim<k>_summary.csv). This script reruns the
# same bootstrap for CM alone, on the main run's samples: resample n_g units with
# replacement within each treatment group, recompute CM (the covariance matrix
# is re-estimated in each draw), B = 200 draws, 2.5% and 97.5% quantiles. Same
# samples and same procedure as the main run, but not the same draws: the main
# run also refitted the GPS in each draw, and random tie-breaking among duplicated
# units used the random-number stream differently, so its intervals cannot be
# regenerated draw for draw without rerunning the GPS bootstrap.
#
# Recorded per replication and pair: the estimate, the mean and SD of the
# bootstrap estimates, and the percentile limits. From these the summary builds
# percentile, basic (2*est - quantiles) and normal (est +- 1.96 SD*) intervals.
#
#   Rscript crosscheck/R/10_cm_bootstrap_centring.R sim=2 [R=200] [B=200] [workers=4]
# Run from the repo root. Writes crosscheck/results/sim<k>_cm_boot_reps.csv;
# per-replication checkpoints go to crosscheck/reps/cmboot_sim<k>/ (git-ignored).

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")

opt <- list(sim = "2", R = "200", B = "200", workers = "4", rep_root = file.path("crosscheck", "reps"))
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2] else stop("unknown argument: ", a)
}
SIM <- as.integer(opt$sim); R_REPS <- as.integer(opt$R); B <- as.integer(opt$B); WORKERS <- as.integer(opt$workers)
rep_dir <- file.path(opt$rep_root, paste0("cmboot_sim", SIM))
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...))); flush(stdout()) }

one_rep <- function(r, sim, seeds, B) {
  setup <- sim_setup(sim); A <- setup$A
  assign(".Random.seed", seeds[[r]], envir = globalenv())
  d <- gen_sample(setup)
  est <- tau_from_imputed(impute_cm(d, A)$Yimp)
  groups <- split(seq_len(nrow(d)), d$T)
  draws <- matrix(NA_real_, B, length(est))
  for (b in seq_len(B)) {
    idx <- unlist(lapply(groups, function(g) g[sample.int(length(g), replace = TRUE)]), use.names = FALSE)
    draws[b, ] <- tau_from_imputed(impute_cm(d[idx, , drop = FALSE], A)$Yimp)
  }
  q <- apply(draws, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
  pr <- pairs_of(A)
  data.frame(rep = r, t = pr[, 1], tp = pr[, 2], est = unname(est),
             boot_mean = colMeans(draws), boot_sd = apply(draws, 2, sd),
             pct_lo = q[1, ], pct_hi = q[2, ])
}

seeds <- make_seeds(R_REPS, 20220000 + SIM)
rep_file <- function(r) file.path(rep_dir, sprintf("rep_%03d.rds", r))
todo <- which(!file.exists(vapply(seq_len(R_REPS), rep_file, character(1))))
log_msg("Simulation %d: %d of %d replications to run, B = %d, %d workers", SIM, length(todo), R_REPS, B, WORKERS)
if (length(todo)) {
  cl <- parallel::makeCluster(WORKERS)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  invisible(parallel::clusterCall(cl, function(p) {
    setwd(p); source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R"); NULL }, getwd()))
  parallel::clusterExport(cl, "one_rep")
  t0 <- proc.time()[["elapsed"]]
  for (ch in split(todo, ceiling(seq_along(todo) / (10 * WORKERS)))) {
    parallel::parLapplyLB(cl, ch, function(r, sim, seeds, B, rep_dir) {
      saveRDS(one_rep(r, sim, seeds, B), file.path(rep_dir, sprintf("rep_%03d.rds", r))); TRUE
    }, sim = SIM, seeds = seeds, B = B, rep_dir = rep_dir)
    log_msg("Simulation %d: %d/%d done, elapsed %.0fs", SIM,
            sum(file.exists(vapply(seq_len(R_REPS), rep_file, character(1)))), R_REPS, proc.time()[["elapsed"]] - t0)
  }
}
files <- vapply(seq_len(R_REPS), rep_file, character(1))
if (all(file.exists(files))) {
  out <- do.call(rbind, lapply(files, readRDS))
  # The point estimates must equal the main run's CM estimates; check via the
  # mean over these replications against crosscheck/results/sim<k>_crosscheck_reps.csv.
  cx <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_crosscheck_reps.csv", SIM)))
  k <- match(paste(out$rep, out$t, out$tp), paste(cx$rep, cx$t, cx$tp))
  log_msg("CM estimates vs the main run's (per replication and pair): max |difference| = %.2e",
          max(abs(out$est - cx$CM[k])))
  num <- vapply(out, is.double, logical(1)); out[num] <- lapply(out[num], signif, digits = 12)
  f <- file.path("crosscheck", "results", sprintf("sim%d_cm_boot_reps.csv", SIM))
  write.csv(out, f, row.names = FALSE)
  log_msg("Wrote %s (%d rows)", f, nrow(out))
}
