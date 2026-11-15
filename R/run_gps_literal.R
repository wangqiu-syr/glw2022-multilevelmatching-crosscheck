# run_gps_literal.R -- D13 sensitivity run: GPS_lit on the Simulation 1 datasets
# of the main run (same seeds, samples and bootstrap draws; paired comparison).
#   Rscript R/run_gps_literal.R check        self-tests only (about 30 s)
#   Rscript R/run_gps_literal.R              full run (500 reps, B from main run, 10 workers)
#   options: workers=10
# Reads (never writes): results/reps/sim1/ (main checkpoints and meta.rds),
#   results/true_effects_sim1.csv, results/paper_table1.csv.
# Writes: results/reps_gps_literal/sim1/rep_###.rds (git-ignored checkpoints),
#   results/gps_literal_sim1.csv.

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
source("R/run_sim.R"); source("R/gps_literal.R")

args <- commandArgs(trailingOnly = TRUE)
check_only <- "check" %in% args
WORKERS <- 10L
for (a in args[grepl("^workers=", args)]) WORKERS <- as.integer(sub("workers=", "", a))

MAIN_DIR <- file.path("results", "reps", "sim1")
LIT_DIR  <- file.path("results", "reps_gps_literal", "sim1")
OUT_CSV  <- file.path("results", "gps_literal_sim1.csv")
meta <- readRDS(file.path(MAIN_DIR, "meta.rds"))
stopifnot(meta$sim == 1, meta$gps_reltol == GPS_RELTOL, meta$gps_maxit == GPS_MAXIT)
s <- sim_setup(1)
seeds <- make_seeds(meta$R, meta$master_seed)
main_file <- function(r) file.path(MAIN_DIR, sprintf("rep_%03d.rds", r))

# Compare the recomputed DIF/CM/GPS rows with the main checkpoint of the same rep.
pairing_diff <- function(lit_est) {
  r <- lit_est$rep[1]
  main <- readRDS(main_file(r))$est
  x <- merge(main, lit_est[lit_est$estimator != "GPS_lit", ],
             by = c("rep", "estimator", "t", "tp"), suffixes = c(".main", ".new"))
  stopifnot(nrow(x) == nrow(main))
  max(abs(c(x$est.main - x$est.new, x$lo.main - x$lo.new, x$hi.main - x$hi.new)))
}

# ---- self-tests ----------------------------------------------------------------
log_msg("GPS_lit self-tests")
# (a) If p_hat_t == p_hat_t' for every unit, GPS_lit must equal GPS.
set.seed(11)
d0 <- gen_sample(s); g0 <- impute_gps(d0, s$A)
P_eq <- g0$P; P_eq[, 2] <- P_eq[, 1]
g_eq <- impute_from_scores(P_eq, d0$T, d0$Y, s$A)
same <- all.equal(unname(tau_gps_lit(P_eq, d0$T, d0$Y, g_eq$Yimp, s$A)[1]),
                  unname(tau_from_imputed(g_eq$Yimp)[1]))
log_msg("  (a) pair (1,2) with p_hat_1 == p_hat_2: GPS_lit equals GPS: %s", isTRUE(same))
# (b) tau_gps_lit leaves the RNG state unchanged.
st <- .Random.seed; invisible(tau_gps_lit(g0$P, d0$T, d0$Y, g0$Yimp, s$A))
rng_ok <- identical(st, .Random.seed)
log_msg("  (b) RNG state unchanged by GPS_lit: %s", rng_ok)
# (c) Replications 1 and 2 reproduce the main checkpoints exactly.
dd <- vapply(1:2, function(r) pairing_diff(run_replication_lit(r, s, meta$B, seeds[[r]])), numeric(1))
log_msg("  (c) reps 1-2, max |difference| from main checkpoints (DIF/CM/GPS est, lo, hi): %g", max(dd))
ok <- isTRUE(same) && rng_ok && all(dd == 0)
log_msg("Self-tests: %s", if (ok) "PASS" else "FAIL")
if (!ok) quit(status = 1)
if (check_only) quit(status = 0)

# ---- full run --------------------------------------------------------------------
dir.create(LIT_DIR, recursive = TRUE, showWarnings = FALSE)
lit_file <- function(r) file.path(LIT_DIR, sprintf("rep_%03d.rds", r))
todo <- which(!file.exists(vapply(seq_len(meta$R), lit_file, character(1))))
log_msg("GPS_lit run: R = %d, B = %d, master seed %d, workers = %d; %d to run",
        meta$R, meta$B, meta$master_seed, WORKERS, length(todo))
if (length(todo)) {
  cl <- parallel::makeCluster(WORKERS)
  invisible(parallel::clusterCall(cl, function(p) {
    setwd(p); source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
    source("R/gps_literal.R"); NULL
  }, getwd()))
  t0 <- proc.time()[["elapsed"]]; n_done <- 0L
  for (ch in split(todo, ceiling(seq_along(todo) / (2 * WORKERS)))) {
    parallel::parLapplyLB(cl, ch, function(r, s, B, seeds, dir) {
      res <- run_replication_lit(r, s, B, seeds[[r]])
      f <- file.path(dir, sprintf("rep_%03d.rds", r)); tmp <- paste0(f, ".tmp")
      saveRDS(res, tmp); file.rename(tmp, f)
    }, s = s, B = meta$B, seeds = seeds, dir = LIT_DIR)
    n_done <- n_done + length(ch); el <- proc.time()[["elapsed"]] - t0
    log_msg("GPS_lit: %d/%d done, elapsed %s, ETA %s", meta$R - length(todo) + n_done, meta$R,
            fmt_dur(el), fmt_dur(el / n_done * (length(todo) - n_done)))
  }
  parallel::stopCluster(cl)
}

# ---- pairing check and summary ---------------------------------------------------------
reps <- lapply(seq_len(meta$R), function(r) readRDS(lit_file(r)))
pd <- vapply(reps, pairing_diff, numeric(1))
log_msg("Pairing check: %d/%d replications reproduce the main DIF/CM/GPS results exactly (max |diff| %g)",
        sum(pd == 0), length(pd), max(pd))

est <- do.call(rbind, reps)
est <- est[est$estimator %in% c("GPS", "GPS_lit"), ]
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
  for (e in c("GPS", "GPS_lit")) {
    x <- est[est$estimator == e & est$t == t & est$tp == tp, ]
    mp <- metrics(x, truth$tau_pop[k]); mb <- metrics(x, truth$tau_bal[k])
    pre <- if (e == "GPS") "gps_" else "gpslit_"
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
  j <- which(paper$t == t & paper$tp == tp)
  row$paper_gps_abs_bias <- paper$GPS_abs_bias[j]
  row$paper_gps_coverage <- paper$GPS_coverage[j]
  row$n_reps <- length(unique(est$rep[est$t == t & est$tp == tp]))
  row
}))
write.csv(out, OUT_CSV, row.names = FALSE)
log_msg("wrote %s", OUT_CSV)
print(format(out[, c("pair", "paper_gps_abs_bias", "gps_abs_bias_pop", "gpslit_abs_bias_pop",
                     "gps_abs_bias_bal", "gpslit_abs_bias_bal", "gps_bias_pop", "gpslit_bias_pop",
                     "gps_mse_pop", "gpslit_mse_pop", "paper_gps_coverage",
                     "gps_coverage_pop", "gpslit_coverage_pop")], digits = 3), row.names = FALSE)
