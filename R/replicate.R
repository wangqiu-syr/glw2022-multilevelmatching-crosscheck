# replicate.R -- one simulation replication and reproducible per-replication
# seeds. Used by R/run_all.R (via R/run_sim.R), R/test_estimators.R and the D13 runners.
# Requires R/dgp.R and R/estimators.R to be sourced.

# One L'Ecuyer-CMRG stream per replication, derived from a master seed.
# Replication r always uses seeds[[r]], so results do not depend on how
# replications are split across workers or in what order they run.
make_seeds <- function(R, master_seed) {
  old_kind <- RNGkind()[1]
  on.exit(RNGkind(old_kind))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(master_seed)
  s <- .Random.seed
  seeds <- vector("list", R)
  for (r in seq_len(R)) {
    seeds[[r]] <- s
    s <- parallel::nextRNGStream(s)
  }
  seeds
}

# Draw one balanced sample, compute the three estimators, bootstrap CIs and
# diagnostics. Returns list(est = long data.frame, diag = list).
run_replication <- function(r, setup, B, seed, maxit = GPS_MAXIT) {
  assign(".Random.seed", seed, envir = globalenv())
  t0 <- proc.time()[["elapsed"]]
  A <- setup$A
  d <- gen_sample(setup)
  pt <- estimate_all(d, A, diagnostics = TRUE, maxit = maxit)
  ci <- boot_ci(d, A, B, maxit = maxit)
  nm <- names(pt$est)
  est <- data.frame(
    rep       = r,
    estimator = sub("\\..*$", "", nm),
    pair      = sub("^[^.]*\\.", "", nm),
    est       = unname(pt$est),
    lo        = unname(ci$lo),
    hi        = unname(ci$hi)
  )
  pr <- do.call(rbind, strsplit(est$pair, "-"))
  est$t  <- as.integer(pr[, 1]); est$tp <- as.integer(pr[, 2])
  est$pair <- NULL
  diag <- c(pt$diag, list(rep = r, B = B, pool_size = attr(d, "pool_size"),
                          boot_gps_nonconverged = ci$n_nonconverged,
                          seconds = proc.time()[["elapsed"]] - t0))
  list(est = est, diag = diag)
}
