# run_gps_literal2.R -- D13 sensitivity run: GPS_lit2 (fully literal Eq. (4)) on
# the Simulation 1 datasets of the main run, alongside GPS and GPS_lit (paired).
#   Rscript R/run_gps_literal2.R check       self-tests only (about 30 s)
#   Rscript R/run_gps_literal2.R             full run (500 reps, B from main run, 10 workers)
#   options: workers=10
# Reads (never writes): results/reps/sim1/ (main checkpoints, meta.rds),
#   results/reps_gps_literal/sim1/ (GPS_lit checkpoints), results/true_effects_sim1.csv,
#   results/paper_table1.csv.
# Writes: results/reps_gps_literal2/sim1/rep_###.rds (git-ignored checkpoints),
#   results/gps_literal2_sim1.csv.

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
source("R/run_sim.R"); source("R/gps_literal.R"); source("R/gps_literal2.R")

args <- commandArgs(trailingOnly = TRUE)
check_only <- "check" %in% args
WORKERS <- 10L
for (a in args[grepl("^workers=", args)]) WORKERS <- as.integer(sub("workers=", "", a))

MAIN_DIR <- file.path("results", "reps", "sim1")
LIT_DIR  <- file.path("results", "reps_gps_literal", "sim1")
LIT2_DIR <- file.path("results", "reps_gps_literal2", "sim1")
OUT_CSV  <- file.path("results", "gps_literal2_sim1.csv")
meta <- readRDS(file.path(MAIN_DIR, "meta.rds"))
stopifnot(meta$sim == 1, meta$gps_reltol == GPS_RELTOL, meta$gps_maxit == GPS_MAXIT)
s <- sim_setup(1)
seeds <- make_seeds(meta$R, meta$master_seed)

# Max |difference| of est, lo, hi between the new run and an earlier checkpoint,
# over the estimators that the earlier checkpoint contains.
diff_vs <- function(new, old, estimators) {
  old <- old[old$estimator %in% estimators, c("rep", "estimator", "t", "tp", "est", "lo", "hi")]
  x <- merge(old, new, by = c("rep", "estimator", "t", "tp"), suffixes = c(".old", ".new"))
  stopifnot(nrow(x) == nrow(old))
  max(abs(c(x$est.old - x$est.new, x$lo.old - x$lo.new, x$hi.old - x$hi.new)))
}
rep_file <- function(dir, r) file.path(dir, sprintf("rep_%03d.rds", r))
pairing <- function(new) {
  r <- new$rep[1]
  c(main = diff_vs(new, readRDS(rep_file(MAIN_DIR, r))$est, c("DIF", "CM", "GPS")),
    lit  = diff_vs(new, readRDS(rep_file(LIT_DIR, r)), "GPS_lit"))
}

# ---- self-tests ------------------------------------------------------------------
log_msg("GPS_lit2 self-tests")
set.seed(11)
d0 <- gen_sample(s); g0 <- impute_gps(d0, s$A)
# (a) If p_hat_1 == p_hat_2 for every unit, GPS_lit2 equals GPS for pair (1,2).
P_eq <- g0$P; P_eq[, 2] <- P_eq[, 1]
g_eq <- impute_from_scores(P_eq, d0$T, d0$Y, s$A)
same <- isTRUE(all.equal(unname(tau_gps_lit2(P_eq, d0$T, d0$Y, g_eq$Yimp, s$A)[1]),
                         unname(tau_from_imputed(g_eq$Yimp)[1])))
log_msg("  (a) pair (1,2) with p_hat_1 == p_hat_2: GPS_lit2 equals GPS: %s", same)
# (b) GPS_lit2 differs from GPS_lit only through units in group t.
col_lit  <- gps_lit_tp_column(g0$P, d0$T, d0$Y, g0$Yimp, 1, 2)
col_lit2 <- gps_lit2_tp_column(g0$P, d0$T, d0$Y, g0$Yimp, 1, 2)
only_t <- all(col_lit[d0$T != 1] == col_lit2[d0$T != 1])
log_msg("  (b) pair (1,2): t' column differs from GPS_lit only for units in group 1: %s (%d of %d group-1 units changed)",
        only_t, sum(col_lit[d0$T == 1] != col_lit2[d0$T == 1]), sum(d0$T == 1))
# (c) GPS_lit2 leaves the RNG state unchanged.
st <- .Random.seed; invisible(tau_gps_lit2(g0$P, d0$T, d0$Y, g0$Yimp, s$A))
rng_ok <- identical(st, .Random.seed)
log_msg("  (c) RNG state unchanged by GPS_lit2: %s", rng_ok)
# (d) Replications 1-2 reproduce the main and GPS_lit checkpoints exactly.
pd <- sapply(1:2, function(r) pairing(run_replication_lit2(r, s, meta$B, seeds[[r]])))
log_msg("  (d) reps 1-2, max |diff| vs main checkpoints %g, vs GPS_lit checkpoints %g",
        max(pd["main", ]), max(pd["lit", ]))
ok <- same && only_t && rng_ok && all(pd == 0)
log_msg("Self-tests: %s", if (ok) "PASS" else "FAIL")
if (!ok) quit(status = 1)
if (check_only) quit(status = 0)

# ---- full run --------------------------------------------------------------------
dir.create(LIT2_DIR, recursive = TRUE, showWarnings = FALSE)
todo <- which(!file.exists(vapply(seq_len(meta$R), function(r) rep_file(LIT2_DIR, r), character(1))))
log_msg("GPS_lit2 run: R = %d, B = %d, master seed %d, workers = %d; %d to run",
        meta$R, meta$B, meta$master_seed, WORKERS, length(todo))
if (length(todo)) {
  cl <- parallel::makeCluster(WORKERS)
  invisible(parallel::clusterCall(cl, function(p) {
    setwd(p); source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
    source("R/gps_literal.R"); source("R/gps_literal2.R"); NULL
  }, getwd()))
  t0 <- proc.time()[["elapsed"]]; n_done <- 0L
  for (ch in split(todo, ceiling(seq_along(todo) / (2 * WORKERS)))) {
    parallel::parLapplyLB(cl, ch, function(r, s, B, seeds, dir) {
      res <- run_replication_lit2(r, s, B, seeds[[r]])
      f <- file.path(dir, sprintf("rep_%03d.rds", r)); tmp <- paste0(f, ".tmp")
      saveRDS(res, tmp); file.rename(tmp, f)
    }, s = s, B = meta$B, seeds = seeds, dir = LIT2_DIR)
    n_done <- n_done + length(ch); el <- proc.time()[["elapsed"]] - t0
    log_msg("GPS_lit2: %d/%d done, elapsed %s, ETA %s", meta$R - length(todo) + n_done, meta$R,
            fmt_dur(el), fmt_dur(el / n_done * (length(todo) - n_done)))
  }
  parallel::stopCluster(cl)
}

# ---- pairing check and summary -----------------------------------------------------
reps <- lapply(seq_len(meta$R), function(r) readRDS(rep_file(LIT2_DIR, r)))
pd <- sapply(reps, pairing)
log_msg("Pairing check: %d/%d reps reproduce main DIF/CM/GPS exactly (max |diff| %g); %d/%d reproduce GPS_lit exactly (max |diff| %g)",
        sum(pd["main", ] == 0), ncol(pd), max(pd["main", ]),
        sum(pd["lit", ] == 0), ncol(pd), max(pd["lit", ]))

EST <- c(GPS = "gps_", GPS_lit = "gpslit_", GPS_lit2 = "gpslit2_")
est <- do.call(rbind, reps)
est <- est[est$estimator %in% names(EST), ]
truth <- read.csv(file.path("results", "true_effects_sim1.csv"))
paper <- read.csv(file.path("results", "paper_table1.csv"))
pp <- do.call(rbind, strsplit(paper$pair, ","))
paper$t <- as.integer(pp[, 1]); paper$tp <- as.integer(pp[, 2])

metrics <- function(e, tau) c(bias = mean(e$est) - tau, abs_bias = abs(mean(e$est) - tau),
                               mse = mean((e$est - tau)^2),
                               coverage = 100 * mean(e$lo <= tau & tau <= e$hi))
out <- do.call(rbind, lapply(seq_len(nrow(truth)), function(k) {
  t <- truth$t[k]; tp <- truth$tp[k]
  row <- data.frame(pair = sprintf("(%d,%d)", t, tp), t = t, tp = tp,
                    tau_pop = truth$tau_pop[k], tau_bal = truth$tau_bal[k])
  j <- which(paper$t == t & paper$tp == tp)
  row$paper_gps_abs_bias <- paper$GPS_abs_bias[j]
  row$paper_gps_coverage <- paper$GPS_coverage[j]
  for (e in names(EST)) {
    x <- est[est$estimator == e & est$t == t & est$tp == tp, ]
    mp <- metrics(x, truth$tau_pop[k]); mb <- metrics(x, truth$tau_bal[k]); pre <- EST[[e]]
    row[[paste0(pre, "abs_bias_pop")]] <- mp[["abs_bias"]]
    row[[paste0(pre, "abs_bias_bal")]] <- mb[["abs_bias"]]
    row[[paste0(pre, "bias_pop")]]     <- mp[["bias"]]
    row[[paste0(pre, "bias_bal")]]     <- mb[["bias"]]
    row[[paste0(pre, "mse_pop")]]      <- mp[["mse"]]
    row[[paste0(pre, "mse_bal")]]      <- mb[["mse"]]
    row[[paste0(pre, "coverage_pop")]] <- mp[["coverage"]]
    row[[paste0(pre, "coverage_bal")]] <- mb[["coverage"]]
    row[[paste0(pre, "mc_se")]]        <- sd(x$est) / sqrt(nrow(x))
  }
  row$n_reps <- length(unique(est$rep[est$t == t & est$tp == tp]))
  row
}))
write.csv(out, OUT_CSV, row.names = FALSE)
log_msg("wrote %s", OUT_CSV)
print(format(out[, c("pair", "paper_gps_abs_bias", "gps_abs_bias_pop", "gpslit_abs_bias_pop",
                     "gpslit2_abs_bias_pop", "gpslit2_bias_pop", "gpslit2_mc_se",
                     "gpslit2_abs_bias_bal", "gps_mse_pop", "gpslit2_mse_pop",
                     "paper_gps_coverage", "gpslit2_coverage_pop")], digits = 3), row.names = FALSE)
