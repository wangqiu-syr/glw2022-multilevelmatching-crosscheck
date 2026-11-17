# 05_summarize_crosscheck.R -- summarize the 500-replication cross-check
# written by 04_crosscheck_all.R, for Simulations 1 and 2.
#   Rscript crosscheck/R/05_summarize_crosscheck.R [sims=1,2]      (from the repo root)
# Reads  crosscheck/results/sim<k>_crosscheck_{reps,info}.csv and, from the
#        main run, results/true_effects_sim<k>.csv and results/sim<k>_summary.csv.
# Writes crosscheck/results/crosscheck_agreement.csv   (how often each fit agrees with ours)
#        crosscheck/results/crosscheck_verification.csv (regenerated samples vs the main run)
#        crosscheck/results/sim<k>_crosscheck_summary.csv (bias, MSE, coverage by pair)
#        crosscheck/results/crosscheck_overview.csv     (means over pairs)

OUT <- file.path("crosscheck", "results")
TOL <- 1e-12                                     # "same estimate" up to floating-point rounding
z <- qnorm(0.975)
SIMS <- 1:2
for (a in commandArgs(trailingOnly = TRUE)) if (startsWith(a, "sims="))
  SIMS <- as.integer(strsplit(sub("sims=", "", a, fixed = TRUE), ",")[[1]])

agree <- list(); verif <- list(); overview <- list()
nz <- function(x) !is.na(x) & nzchar(x)       # read.csv turns "" into NA
for (sim in SIMS) {
  est  <- read.csv(file.path(OUT, sprintf("sim%d_crosscheck_reps.csv.gz", sim)))
  info <- read.csv(file.path(OUT, sprintf("sim%d_crosscheck_info.csv.gz", sim)))
  te   <- read.csv(file.path("results", sprintf("true_effects_sim%d.csv", sim)))
  main <- read.csv(file.path("results", sprintf("sim%d_summary.csv", sim)))
  R <- nrow(info)

  # ---- 1. regenerated samples and estimates vs the main run ------------------
  mean_check <- sapply(c("DIF", "CM", "GPS"), function(e) {
    m <- aggregate(est[[e]], est[c("t", "tp")], mean)
    ref <- main[main$estimator == e & main$target == "pop", c("t", "tp", "mean_est")]
    mm <- merge(m, ref)
    max(abs(mm$x - mm$mean_est))
  })
  verif[[sim]] <- data.frame(
    sim = sim, n_reps = R,
    # GPS matching distances differ by up to ~3e-9 from the main run because the
    # fitted probabilities differ slightly across platforms (Linux nnet 7.3-19 vs
    # Windows nnet 7.3.20); the estimates themselves are checked via the means below.
    reps_diagnostics_within_1e8 = sum(info$diag_maxdiff < 1e-8),
    max_diag_diff = max(info$diag_maxdiff),
    max_diff_mean_DIF = mean_check[["DIF"]], max_diff_mean_CM = mean_check[["CM"]],
    max_diff_mean_GPS = mean_check[["GPS"]],
    ours_gps_converged = sum(info$ours_gps_converged),
    reps_with_ties_cm = sum(info$ties_cm > 0), reps_with_ties_gps = sum(info$ties_gps > 0),
    rows_nudged_mean = mean(info$rows_nudged), max_nudge = max(info$max_nudge),
    package_errors = sum(nz(info$err_cm) | nz(info$err_mn) | nz(info$err_ex)),
    package_warnings = sum(nz(info$warn_cm) | nz(info$warn_mn) | nz(info$warn_ex)))

  # ---- 2. agreement, by package fit -----------------------------------------
  fits <- list(covariates = c("pkg_CM", "CM", "cells_diff_cm"),
               multinom   = c("pkg_GPS_multinom", "GPS", "cells_diff_mn"),
               existing   = c("pkg_GPS_existing", "GPS", "cells_diff_ex"))
  n_cells <- nrow(est) / R            # pairs per rep; cells are N x A
  for (f in names(fits)) {
    v <- fits[[f]]
    dif <- est[[v[1]]] - est[[v[2]]]
    rep_max <- tapply(abs(dif), est$rep, max)
    agree[[length(agree) + 1]] <- data.frame(
      sim = sim, fit = f, compared_with = paste("our", v[2]), n_reps = R,
      reps_identical_imputations = sum(info[[v[3]]] == 0, na.rm = TRUE),
      reps_same_estimates = sum(rep_max < TOL),
      max_abs_diff = max(abs(dif)),
      median_rep_max_abs_diff = median(rep_max),
      mean_cells_differing = mean(info[[v[3]]]),
      max_cells_differing = max(info[[v[3]]]),
      max_abs_dP = if (f == "multinom") max(info$mn_max_dP) else NA,
      median_abs_dP = if (f == "multinom") median(info$mn_max_dP) else NA,
      pkg_multinom_not_converged = if (f == "multinom") sum(info$mn_convergence != 0) else NA)
  }

  # ---- 3. bias, MSE and coverage by pair ------------------------------------
  long <- list()
  add <- function(name, e, v = NULL, var_label = NA) {
    for (tg in c("pop", "bal")) {
      tau <- setNames(te[[paste0("tau_", tg)]], paste(te$t, te$tp))[paste(est$t, est$tp)]
      df <- data.frame(t = est$t, tp = est$tp, rep = est$rep, est = est[[e]], tau = tau)
      if (!is.null(v)) {
        vv <- est[[v]]; se <- sqrt(ifelse(vv < 0, NA, vv))   # negative variance -> no interval
        df$cover <- abs(df$est - df$tau) <= z * se
        df$se <- se; df$neg_var <- est[[v]] < 0
      } else { df$cover <- NA; df$se <- NA; df$neg_var <- NA }
      a <- do.call(rbind, lapply(split(df, list(df$t, df$tp), drop = TRUE), function(g) data.frame(
        estimator = name, t = g$t[1], tp = g$tp[1], target = tg, tau = g$tau[1],
        mean_est = mean(g$est), bias = mean(g$est) - g$tau[1], abs_bias = abs(mean(g$est) - g$tau[1]),
        mc_se = sd(g$est) / sqrt(nrow(g)), sd_est = sd(g$est), mse = mean((g$est - g$tau[1])^2),
        variance = var_label,
        coverage = if (all(is.na(g$cover))) NA else 100 * mean(g$cover, na.rm = TRUE),
        mean_se = if (all(is.na(g$se))) NA else mean(g$se, na.rm = TRUE),
        n_neg_var = if (all(is.na(g$neg_var))) NA else sum(g$neg_var))))
      long[[length(long) + 1]] <<- a
    }
  }
  add("ours_CM", "CM"); add("ours_GPS", "GPS")
  add("pkg_CM", "pkg_CM", "var_pkg_CM", "AI2006")
  add("pkg_GPS_multinom", "pkg_GPS_multinom", "var_pkg_GPS_multinom", "AI2006")
  add("pkg_GPS_multinom_AI2016", "pkg_GPS_multinom", "varAI2016_pkg_GPS_multinom", "AI2016")
  add("pkg_GPS_existing", "pkg_GPS_existing", "var_pkg_GPS_existing", "AI2006")
  summ <- do.call(rbind, long)
  # Our bootstrap coverage from the main run, for comparison.
  boot <- main[main$estimator %in% c("CM", "GPS"), c("estimator", "t", "tp", "target", "coverage")]
  boot$estimator <- paste0("ours_", boot$estimator)
  names(boot)[5] <- "coverage_boot"
  summ <- merge(summ, boot, all.x = TRUE)
  summ$variance[summ$estimator %in% c("ours_CM", "ours_GPS")] <- "bootstrap (main run)"
  summ$coverage[summ$estimator %in% c("ours_CM", "ours_GPS")] <-
    summ$coverage_boot[summ$estimator %in% c("ours_CM", "ours_GPS")]
  summ$coverage_boot <- NULL
  summ <- summ[order(summ$target, summ$estimator, summ$t, summ$tp), ]
  write.csv(summ, file.path(OUT, sprintf("sim%d_crosscheck_summary.csv", sim)), row.names = FALSE)

  # Paired difference: package multinom GPS minus our GPS (does the default
  # tolerance shift the estimates systematically?).
  dmn <- est$pkg_GPS_multinom - est$GPS
  pd <- aggregate(dmn, est[c("t", "tp")], function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
  pd <- data.frame(t = pd$t, tp = pd$tp, mean_diff = pd$x[, "mean"], se = pd$x[, "se"])
  pd$z <- pd$mean_diff / pd$se
  write.csv(pd, file.path(OUT, sprintf("sim%d_crosscheck_multinom_minus_ours.csv", sim)), row.names = FALSE)

  ov <- aggregate(cbind(abs_bias, mse, coverage, mean_se) ~ estimator + target, data = summ,
                  FUN = mean, na.action = na.pass)
  ov$sim <- sim
  ov$pairs_abs_z_gt_3 <- NA
  ov$pairs_abs_z_gt_3[ov$estimator == "pkg_GPS_multinom"] <- sum(abs(pd$z) > 3)
  overview[[sim]] <- ov
}

agree <- do.call(rbind, agree); verif <- do.call(rbind, verif); overview <- do.call(rbind, overview)
overview <- overview[c("sim", setdiff(names(overview), "sim"))]
write.csv(agree, file.path(OUT, "crosscheck_agreement.csv"), row.names = FALSE)
write.csv(verif, file.path(OUT, "crosscheck_verification.csv"), row.names = FALSE)
write.csv(overview, file.path(OUT, "crosscheck_overview.csv"), row.names = FALSE)

options(width = 200)
cat("== Verification\n"); print(verif)
cat("\n== Agreement\n"); print(agree, digits = 4)
cat("\n== Overview (means over pairs)\n"); print(overview, digits = 4)
for (sim in SIMS) {
  cat(sprintf("\n== Sim %d: package multinom GPS minus our GPS, by pair\n", sim))
  print(read.csv(file.path(OUT, sprintf("sim%d_crosscheck_multinom_minus_ours.csv", sim))), digits = 3)
}
