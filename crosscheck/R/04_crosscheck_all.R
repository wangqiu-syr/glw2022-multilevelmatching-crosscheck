# 04_crosscheck_all.R -- cross-check every replication of Simulation 1 or 2
# against multilevelMatching::multiMatch (match_on = "covariates", "multinom",
# "existing"). For each replication r:
#   1. regenerate the sample from the main run's seed stream (master seed
#      20220000 + sim, L'Ecuyer stream r) and check it against the stored
#      diagnostics in results/sim<k>_diagnostics.csv;
#   2. compute our DIF, CM and GPS exactly as estimate_all() does;
#   3. fit multiMatch() three ways, each after set.seed(123) (Mersenne-Twister),
#      with trimming = FALSE and M_matches = 1; "existing" gets our fitted GPS
#      matrix (rows nudged by a few ulp so they sum to exactly 1);
#   4. record estimates (package signs flipped to tau(t,t') = E[Y(t) - Y(t')]),
#      the package's Abadie-Imbens variances, and unit-level diagnostics.
# One checkpoint per replication in <rep_root>/sim<k>/rep_###.rds (git-ignored);
# re-running skips finished replications. When all are done, the combined
# per-replication tables are written to crosscheck/results/.
#
# Run from the repo root, with multilevelMatching on the library path:
#   R_LIBS=~/R/mlm-lib Rscript crosscheck/R/04_crosscheck_all.R sim=1 [R=500] [workers=4] [reps=1:20]
# Package workarounds (R >= 4.3 class check; exact row sums) are described in docs/crosscheck.md, Appendix B.2.

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")

opt <- list(sim = "1", R = "500", workers = "4", reps = "", rep_root = file.path("crosscheck", "reps"),
            out = file.path("crosscheck", "results"))
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2] else stop("unknown argument: ", a)
}
SIM <- as.integer(opt$sim); R_REPS <- as.integer(opt$R); WORKERS <- as.integer(opt$workers)
REPS <- if (nzchar(opt$reps)) eval(parse(text = opt$reps)) else seq_len(R_REPS)
rep_dir <- file.path(opt$rep_root, paste0("sim", SIM))
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...))); flush(stdout())
}

# ---- one replication (runs on a worker) -------------------------------------

crosscheck_rep <- function(r, sim, seeds, diag_ref) {
  t0 <- proc.time()[["elapsed"]]
  setup <- sim_setup(sim); A <- setup$A

  # 1. Sample, exactly as run_replication() draws it.
  assign(".Random.seed", seeds[[r]], envir = globalenv())
  d <- gen_sample(setup)

  # 2. Our estimators: estimate_all() calls impute_cm() then impute_gps(), and
  #    nothing else draws random numbers, so this reproduces its estimates.
  cm  <- impute_cm(d, A)
  gps <- impute_gps(d, A)
  ours <- data.frame(t = pairs_of(A)[, 1], tp = pairs_of(A)[, 2],
                     DIF = unname(tau_dif(d, A)),
                     CM  = unname(tau_from_imputed(cm$Yimp)),
                     GPS = unname(tau_from_imputed(gps$Yimp)))

  ref <- diag_ref[diag_ref$rep == r, ]
  mine <- c(attr(d, "pool_size"), cm$match_dist, gps$match_dist)
  refv <- unlist(ref[c("pool_size", paste0("cm_dist_t", 1:A), paste0("gps_dist_t", 1:A))])
  diag_maxdiff <- max(abs(mine - refv))

  # Ties and near-ties in our own matching.
  M <- as.matrix(d[, COV_NAMES]); Dm <- sqrt(mahalanobis_d2(M))
  P <- gps$P
  tie_info <- function(dist_to) {
    n_tie <- 0; gap <- Inf
    for (s in seq_len(A)) {
      g <- which(d$T == s); D <- dist_to(g, s)[d$T != s, , drop = FALSE]
      dmin <- apply(D, 1, min)
      n_tie <- n_tie + sum(rowSums(D == dmin) > 1)
      if (ncol(D) > 1) gap <- min(gap, apply(D, 1, function(x) { o <- sort(x, partial = 1:2); o[2] - o[1] }))
    }
    c(n_tie = n_tie, min_gap = gap)
  }
  ti_cm  <- tie_info(function(g, s) Dm[, g, drop = FALSE])
  ti_gps <- tie_info(function(g, s) abs(outer(P[, s], P[g, s], "-")))

  # 3. Package fits.
  Y <- d$Y; W <- d$T
  X <- M; oldClass(X) <- "matrix"                # R >= 4.3 class-check workaround
  P_ex <- P
  bad <- which(rowSums(P_ex) != 1)
  for (i in bad) {                               # exact row sums for "existing"
    j <- which.max(P_ex[i, ])
    P_ex[i, j] <- P_ex[i, j] + (1 - sum(P_ex[i, ]))
    k <- 0
    while (sum(P_ex[i, ]) != 1 && k < 50) {
      P_ex[i, j] <- P_ex[i, j] + sign(1 - sum(P_ex[i, ])) * .Machine$double.eps * P_ex[i, j] / 2
      k <- k + 1
    }
  }
  rows_exact <- all(rowSums(P_ex) == 1)
  Xe <- P_ex; colnames(Xe) <- paste0("p", 1:A); oldClass(Xe) <- "matrix"

  fit_one <- function(method, XX) {
    warn <- character(0)
    set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    f <- withCallingHandlers(
      tryCatch(suppressMessages(multiMatch(Y, W, XX, trimming = FALSE, match_on = method,
                                           M_matches = 1)),
               error = function(e) e),
      warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
    list(fit = f, warn = warn)
  }
  fits <- list(cov = fit_one("covariates", X), mn = fit_one("multinom", X),
               ex = if (rows_exact) fit_one("existing", Xe) else list(fit = simpleError("rows not exact"), warn = character(0)))

  key <- paste0(ours$t, "-", ours$tp)
  get_col <- function(f, col, flip = FALSE) {
    if (inherits(f, "error")) return(rep(NA_real_, length(key)))
    r <- f$results
    v <- setNames(r[[col]], paste0(r$Trt1, "-", r$Trt2))[key]
    unname(if (flip) -v else v)
  }
  est <- cbind(rep = r, ours,
               pkg_CM = get_col(fits$cov$fit, "Estimate", TRUE),
               pkg_GPS_multinom = get_col(fits$mn$fit, "Estimate", TRUE),
               pkg_GPS_existing = get_col(fits$ex$fit, "Estimate", TRUE),
               var_pkg_CM = get_col(fits$cov$fit, "Variance"),
               var_pkg_GPS_multinom = get_col(fits$mn$fit, "Variance"),
               varAI2016_pkg_GPS_multinom = get_col(fits$mn$fit, "VarianceAI2016"),
               var_pkg_GPS_existing = get_col(fits$ex$fit, "Variance"))

  cells_diff <- function(f, Yimp) {
    if (inherits(f, "error")) return(NA_integer_)
    sum(unname(as.matrix(f$impute_mat)) != unname(Yimp))
  }
  mn_ok <- !inherits(fits$mn$fit, "error")
  info <- data.frame(
    rep = r, seconds = NA_real_, pool_size = attr(d, "pool_size"),
    diag_maxdiff = diag_maxdiff, ours_gps_converged = gps$converged,
    ties_cm = ti_cm[["n_tie"]], ties_gps = ti_gps[["n_tie"]],
    min_gap_cm = ti_cm[["min_gap"]], min_gap_gps = ti_gps[["min_gap"]],
    rows_nudged = length(bad), max_nudge = max(abs(P_ex - P)), rows_exact = rows_exact,
    cells_diff_cm = cells_diff(fits$cov$fit, cm$Yimp),
    cells_diff_mn = cells_diff(fits$mn$fit, gps$Yimp),
    cells_diff_ex = cells_diff(fits$ex$fit, gps$Yimp),
    mn_convergence = if (mn_ok) fits$mn$fit$model$convergence else NA_integer_,
    mn_max_dP = if (mn_ok) max(abs(unname(fits$mn$fit$propensity_scores) - P)) else NA_real_,
    err_cm = if (inherits(fits$cov$fit, "error")) conditionMessage(fits$cov$fit) else "",
    err_mn = if (inherits(fits$mn$fit, "error")) conditionMessage(fits$mn$fit) else "",
    err_ex = if (inherits(fits$ex$fit, "error")) conditionMessage(fits$ex$fit) else "",
    warn_cm = paste(unique(fits$cov$warn), collapse = " | "),
    warn_mn = paste(unique(fits$mn$warn), collapse = " | "),
    warn_ex = paste(unique(fits$ex$warn), collapse = " | "),
    stringsAsFactors = FALSE)
  info$seconds <- proc.time()[["elapsed"]] - t0
  list(est = est, info = info)
}

# ---- driver -----------------------------------------------------------------

seeds <- make_seeds(R_REPS, 20220000 + SIM)          # as run_sim(): master_seed = 20220000 + sim
diag_ref <- read.csv(file.path("results", sprintf("sim%d_diagnostics.csv", SIM)))
rep_file <- function(r) file.path(rep_dir, sprintf("rep_%03d.rds", r))
todo <- REPS[!file.exists(vapply(REPS, rep_file, character(1)))]
log_msg("Cross-check, Simulation %d: %d replications requested, %d to run, %d workers",
        SIM, length(REPS), length(todo), WORKERS)
log_msg("%s; multilevelMatching %s; Matching %s; nnet %s; MASS %s", R.version.string,
        packageVersion("multilevelMatching"), packageVersion("Matching"),
        packageVersion("nnet"), packageVersion("MASS"))

if (length(todo)) {
  cl <- parallel::makeCluster(WORKERS)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  invisible(parallel::clusterCall(cl, function(p) {
    setwd(p); source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
    suppressPackageStartupMessages(library(multilevelMatching)); NULL
  }, getwd()))
  parallel::clusterExport(cl, "crosscheck_rep")
  worker_fun <- function(r, sim, seeds, diag_ref, rep_dir) {
    f <- file.path(rep_dir, sprintf("rep_%03d.rds", r))
    res <- tryCatch(crosscheck_rep(r, sim, seeds, diag_ref), error = function(e) e)
    if (inherits(res, "error")) {
      writeLines(conditionMessage(res), sub("\\.rds$", "_error.txt", f)); return(FALSE)
    }
    tmp <- paste0(f, ".tmp"); saveRDS(res, tmp); file.rename(tmp, f)
  }
  chunks <- split(todo, ceiling(seq_along(todo) / (5 * WORKERS)))
  t0 <- proc.time()[["elapsed"]]; n_done <- 0L; n_err <- 0L
  for (ch in chunks) {
    ok <- unlist(parallel::parLapplyLB(cl, ch, worker_fun, sim = SIM, seeds = seeds,
                                       diag_ref = diag_ref, rep_dir = rep_dir))
    n_done <- n_done + sum(ok); n_err <- n_err + sum(!ok)
    el <- proc.time()[["elapsed"]] - t0
    log_msg("Simulation %d: %d/%d done (%d errors), elapsed %.0fs, ETA %.0fs", SIM,
            n_done + n_err, length(todo), n_err, el, el / (n_done + n_err) * (length(todo) - n_done - n_err))
  }
}

# ---- combine when every requested replication is done ----------------------

files <- vapply(REPS, rep_file, character(1))
if (all(file.exists(files)) && identical(sort(REPS), seq_len(R_REPS))) {
  res <- lapply(files, readRDS)
  est  <- do.call(rbind, lapply(res, `[[`, "est"))
  info <- do.call(rbind, lapply(res, `[[`, "info"))
  num <- vapply(est, is.double, logical(1))
  est[num] <- lapply(est[num], signif, digits = 12)
  write.csv(est,  file.path(opt$out, sprintf("sim%d_crosscheck_reps.csv", SIM)), row.names = FALSE)
  write.csv(info, file.path(opt$out, sprintf("sim%d_crosscheck_info.csv", SIM)), row.names = FALSE)
  log_msg("Wrote %s and %s", sprintf("sim%d_crosscheck_reps.csv", SIM), sprintf("sim%d_crosscheck_info.csv", SIM))
} else {
  log_msg("Not combining: %d of %d replications have checkpoints", sum(file.exists(files)), R_REPS)
}
