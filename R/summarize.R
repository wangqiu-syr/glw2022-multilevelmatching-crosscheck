# summarize.R -- combine checkpointed replications into summary tables
# (docs/spec.md, "Outputs"). Can be sourced, or run on its own after the
# simulation has finished:
#   Rscript R/summarize.R
# Writes, in results/:
#   sim1_summary.csv, sim2_summary.csv          per estimator x pair x target
#   sim1_diagnostics.csv, sim2_diagnostics.csv  per replication
#   sim1_vs_paper.csv                           Simulation 1 vs paper Table 1

load_reps <- function(rep_dir) {
  files <- list.files(rep_dir, pattern = "^rep_[0-9]+\\.rds$", full.names = TRUE)
  if (!length(files)) stop("load_reps: no replication files in ", rep_dir)
  lapply(files, readRDS)
}

diag_table <- function(reps) {
  do.call(rbind, lapply(reps, function(x) {
    g <- x$diag
    A <- length(g$cm_match_dist)
    data.frame(rep = g$rep, B = g$B, seconds = g$seconds, pool_size = g$pool_size,
               gps_converged = g$gps_converged,
               boot_gps_nonconverged = g$boot_gps_nonconverged,
               p_hat_any_below_01 = g$p_hat_any_below_01,
               setNames(as.list(g$cm_match_dist),  paste0("cm_dist_t", 1:A)),
               setNames(as.list(g$gps_match_dist), paste0("gps_dist_t", 1:A)),
               setNames(as.list(g$p_hat_below_01), paste0("p_hat_below_01_t", 1:A)))
  }))
}

summarize_sim <- function(sim, rep_dir = file.path("results", "reps", paste0("sim", sim)),
                          out_dir = "results", truth_dir = "results") {
  reps <- load_reps(rep_dir)
  est <- do.call(rbind, lapply(reps, `[[`, "est"))
  truth <- read.csv(file.path(truth_dir, sprintf("true_effects_sim%d.csv", sim)))
  long <- rbind(
    merge(est, data.frame(t = truth$t, tp = truth$tp, target = "pop", tau = truth$tau_pop)),
    merge(est, data.frame(t = truth$t, tp = truth$tp, target = "bal", tau = truth$tau_bal)))
  key <- interaction(long$estimator, long$t, long$tp, long$target, drop = TRUE, lex.order = TRUE)
  summ <- do.call(rbind, lapply(split(long, key), function(e) {
    m <- mean(e$est)
    data.frame(estimator = e$estimator[1], t = e$t[1], tp = e$tp[1], target = e$target[1],
               n_reps = nrow(e), tau = e$tau[1], mean_est = m,
               bias = m - e$tau[1], abs_bias = abs(m - e$tau[1]),
               sd_est = sd(e$est), mse = mean((e$est - e$tau[1])^2),
               coverage = 100 * mean(e$lo <= e$tau[1] & e$tau[1] <= e$hi),
               mean_ci_width = mean(e$hi - e$lo))
  }))
  summ$target <- factor(summ$target, levels = c("pop", "bal"))
  summ <- summ[order(summ$target, summ$estimator, summ$t, summ$tp), ]
  summ$target <- as.character(summ$target)
  rownames(summ) <- NULL
  f <- file.path(out_dir, sprintf("sim%d_summary.csv", sim))
  write.csv(summ, f, row.names = FALSE)

  dg <- diag_table(reps)
  write.csv(dg, file.path(out_dir, sprintf("sim%d_diagnostics.csv", sim)), row.names = FALSE)

  cat(sprintf("\nSimulation %d summary (%d replications), averaged over pairs:\n",
              sim, length(reps)))
  ov <- aggregate(cbind(abs_bias, mse, coverage, mean_ci_width) ~ estimator + target,
                  data = summ, FUN = mean)
  print(format(ov, digits = 3), row.names = FALSE)
  cat(sprintf("  GPS converged in %d/%d samples; non-converged bootstrap fits: %d; mean seconds/rep %.1f\n",
              sum(dg$gps_converged), nrow(dg), sum(dg$boot_gps_nonconverged), mean(dg$seconds)))
  cat(sprintf("  written: %s and sim%d_diagnostics.csv\n", f, sim))
  flush(stdout())
  invisible(summ)
}

# Simulation 1 vs the paper's Table 1 (results/paper_table1.csv; pair "t,tp").
compare_paper <- function(summary_file = file.path("results", "sim1_summary.csv"),
                          paper_file = file.path("results", "paper_table1.csv"),
                          out_file = file.path("results", "sim1_vs_paper.csv")) {
  summ <- read.csv(summary_file)
  paper <- read.csv(paper_file)
  pr <- do.call(rbind, strsplit(paper$pair, ","))
  paper$t <- as.integer(pr[, 1]); paper$tp <- as.integer(pr[, 2])
  ours <- function(estim, target, col) {
    x <- summ[summ$estimator == estim & summ$target == target, ]
    x[match(paste(paper$t, paper$tp), paste(x$t, x$tp)), col]
  }
  out <- do.call(rbind, lapply(c("DIF", "CM", "GPS"), function(estim) {
    pcov <- paste0(estim, "_coverage")
    data.frame(
      pair = paste0("(", paper$t, ",", paper$tp, ")"), t = paper$t, tp = paper$tp,
      estimator = estim,
      paper_abs_bias = paper[[paste0(estim, "_abs_bias")]],
      ours_abs_bias_pop = ours(estim, "pop", "abs_bias"),
      ours_abs_bias_bal = ours(estim, "bal", "abs_bias"),
      paper_coverage = if (pcov %in% names(paper)) paper[[pcov]] else NA_real_,
      ours_coverage_pop = ours(estim, "pop", "coverage"),
      ours_coverage_bal = ours(estim, "bal", "coverage"),
      n_reps = ours(estim, "pop", "n_reps"))
  }))
  out$abs_bias_diff_pop <- out$ours_abs_bias_pop - out$paper_abs_bias
  out$coverage_diff_pop <- out$ours_coverage_pop - out$paper_coverage
  write.csv(out, out_file, row.names = FALSE)
  cat("\nSimulation 1 vs paper Table 1 (abs bias = |mean(tau_hat) - tau|; coverage %):\n")
  print(format(out[, c("pair", "estimator", "paper_abs_bias", "ours_abs_bias_pop",
                       "ours_abs_bias_bal", "paper_coverage", "ours_coverage_pop",
                       "ours_coverage_bal")], digits = 3), row.names = FALSE)
  cat(sprintf("  written: %s\n", out_file))
  flush(stdout())
  invisible(out)
}

if (sys.nframe() == 0) {
  summarize_sim(1); summarize_sim(2); compare_paper()
}
