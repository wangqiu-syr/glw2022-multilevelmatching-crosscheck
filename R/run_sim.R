# run_sim.R -- parallel, checkpointed driver for one simulation.
# Requires R/dgp.R, R/estimators.R and R/replicate.R to be sourced.
#
# Each replication r uses the pre-generated L'Ecuyer stream seeds[[r]] and is
# saved to <rep_dir>/rep_###.rds as soon as it finishes. Re-running skips
# replications whose file already exists, so an interrupted run can resume.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
  flush(stdout())
}

fmt_dur <- function(sec) {
  if (!is.finite(sec)) return("?")
  sprintf("%dh%02dm%02ds", as.integer(sec %/% 3600), as.integer(sec %% 3600 %/% 60),
          as.integer(round(sec %% 60)))
}

run_sim <- function(sim, R = 500, B = 200, workers = 10,
                    master_seed = 20220000 + sim,
                    rep_dir = file.path("results", "reps", paste0("sim", sim))) {
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  s <- sim_setup(sim)

  # Refuse to mix checkpoints made with different settings.
  meta <- list(sim = sim, R = R, B = B, master_seed = master_seed,
               gps_reltol = GPS_RELTOL, gps_maxit = GPS_MAXIT)
  meta_file <- file.path(rep_dir, "meta.rds")
  if (file.exists(meta_file)) {
    if (!identical(readRDS(meta_file), meta))
      stop("run_sim: ", rep_dir, " holds checkpoints made with different settings")
  } else saveRDS(meta, meta_file)

  seeds <- make_seeds(R, master_seed)
  rep_file <- function(r) file.path(rep_dir, sprintf("rep_%03d.rds", r))
  todo <- which(!file.exists(vapply(seq_len(R), rep_file, character(1))))
  log_msg("Simulation %d: A = %d, n_g = %d, R = %d, B = %d, workers = %d, master seed %d",
          sim, s$A, s$n_g, R, B, workers, master_seed)
  log_msg("Simulation %d: %d replications already done, %d to run", sim, R - length(todo), length(todo))
  if (!length(todo)) return(invisible(TRUE))

  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  invisible(parallel::clusterCall(cl, function(p) {
    setwd(p); source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R"); NULL
  }, getwd()))

  worker_fun <- function(r, s, B, seeds, rep_dir) {
    f <- file.path(rep_dir, sprintf("rep_%03d.rds", r))
    res <- tryCatch(run_replication(r, s, B, seeds[[r]]), error = function(e) e)
    if (inherits(res, "error")) {
      writeLines(conditionMessage(res), sub("\\.rds$", "_error.txt", f))
      return(FALSE)
    }
    tmp <- paste0(f, ".tmp")
    saveRDS(res, tmp)
    file.rename(tmp, f)
  }

  chunks <- split(todo, ceiling(seq_along(todo) / (2 * workers)))
  t0 <- proc.time()[["elapsed"]]
  n_done <- 0L; n_err <- 0L
  for (ch in chunks) {
    ok <- unlist(parallel::parLapplyLB(cl, ch, worker_fun,
                                       s = s, B = B, seeds = seeds, rep_dir = rep_dir))
    n_done <- n_done + sum(ok); n_err <- n_err + sum(!ok)
    el <- proc.time()[["elapsed"]] - t0
    left <- length(todo) - n_done - n_err
    log_msg("Simulation %d: %d/%d done (%d errors), elapsed %s, ETA %s",
            sim, R - length(todo) + n_done, R, n_err, fmt_dur(el),
            fmt_dur(el / (n_done + n_err) * left))
  }
  if (n_err) log_msg("Simulation %d: %d replications failed; see %s/*_error.txt", sim, n_err, rep_dir)
  invisible(n_err == 0)
}
