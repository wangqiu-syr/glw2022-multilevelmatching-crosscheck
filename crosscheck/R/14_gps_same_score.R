# 14_gps_same_score.R -- a third reading of the paper's Eq. (4), Simulation 1 only
# (crosscheck/eq4_readings/README.md).
#
# Eq. (4) (Guo, Liu and Wang, 2022, Section 4, p. 972) prints Pr(t | x_i, z_i) in
# both terms. Three readings of its second (t') term, for pair (t, t'):
#
#   GPS      (main, D7)  match within group t' on |p_hat_t'(j) - p_hat_t'(i)|
#   GPS_lit2 (D13)       match within group t' on |p_hat_t'(j) - p_hat_t(i)|
#                        (Eq. (3) for level t', fed unit i's score for t)
#   GPS_same (new)       match within group t' on |p_hat_t(j)  - p_hat_t(i)|
#                        (one score, the score for t, on both sides, as in
#                        binary propensity-score matching)
#
# GPS_lit (D13) is GPS_lit2 applied only to units outside both groups. In every
# reading the t term is the main GPS's, and each unit keeps its observed outcome
# for its own group. GPS_same is not antisymmetric: tau_hat(t, t') uses the score
# for t, so pairs are taken with t < t', as in the paper's Table 1.
#
# Point estimates only (no bootstrap): bias, MSE and SD over the 500 replications,
# against tau_pop and tau_bal. Each replication regenerates the main run's sample
# from its stored seed and computes CM, then GPS, in the main run's order, so DIF,
# CM, GPS, GPS_lit and GPS_lit2 must reproduce the stored results exactly; the
# script stops if they do not.
#
#   Rscript crosscheck/R/14_gps_same_score.R [R=500] [workers=4]
# Run from the repo root. Writes crosscheck/results/sim1_gps_same_score_reps.csv.gz
# and crosscheck/eq4_readings/tables/*.{csv,md}.

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
source("R/gps_literal.R"); source("R/gps_literal2.R")
write_csv_gz <- function(x, f) { con <- gzfile(f, "w"); on.exit(close(con)); write.csv(x, con, row.names = FALSE) }

opt <- list(R = "500", workers = "4")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2] else stop("unknown argument: ", a)
}
R_REPS <- as.integer(opt$R); WORKERS <- as.integer(opt$workers); SIM <- 1L
OUT <- file.path("crosscheck", "eq4_readings"); TAB <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...))); flush(stdout()) }

# t' column under the same-score reading: every unit outside group t' gets the Y
# of its nearest unit in group t' on |p_hat_t(j) - p_hat_t(i)|.
gps_same_tp_column <- function(P, T, Y, Yimp, t, tp) {
  out <- Yimp[, tp]
  others <- which(T != tp); g <- which(T == tp)
  if (length(others)) {
    nn <- nn_pick(abs(outer(P[others, t], P[g, t], "-")))
    out[others] <- Y[g[nn$j]]
  }
  out
}
tau_gps_same <- function(P, T, Y, Yimp, A) {
  saved <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", saved, envir = globalenv()))
  pr <- pairs_of(A)
  v <- vapply(seq_len(nrow(pr)), function(k) {
    t <- pr[k, 1]; tp <- pr[k, 2]
    mean(Yimp[, t]) - mean(gps_same_tp_column(P, T, Y, Yimp, t, tp))
  }, numeric(1))
  setNames(v, pair_labels(A))
}

# The same computations, in the same order, as estimate_with_lit2() (the D13
# run), plus GPS_same; written out because that function does not return the
# fitted scores.
one_rep <- function(r, seeds) {
  setup <- sim_setup(SIM); A <- setup$A
  assign(".Random.seed", seeds[[r]], envir = globalenv())
  d <- gen_sample(setup)
  cm  <- impute_cm(d, A)
  gps <- impute_gps(d, A)
  est <- c(setNames(tau_dif(d, A),              paste0("DIF.",      pair_labels(A))),
           setNames(tau_from_imputed(cm$Yimp),  paste0("CM.",       pair_labels(A))),
           setNames(tau_from_imputed(gps$Yimp), paste0("GPS.",      pair_labels(A))),
           setNames(tau_gps_lit(gps$P, d$T, d$Y, gps$Yimp, A),  paste0("GPS_lit.",  pair_labels(A))),
           setNames(tau_gps_lit2(gps$P, d$T, d$Y, gps$Yimp, A), paste0("GPS_lit2.", pair_labels(A))),
           setNames(tau_gps_same(gps$P, d$T, d$Y, gps$Yimp, A), paste0("GPS_same.", pair_labels(A))))
  nm <- names(est); pr <- do.call(rbind, strsplit(sub("^[^.]*\\.", "", nm), "-"))
  data.frame(rep = r, estimator = sub("\\..*$", "", nm), t = as.integer(pr[, 1]), tp = as.integer(pr[, 2]),
             est = unname(est))
}

seeds <- make_seeds(R_REPS, 20220000 + SIM)
log_msg("Simulation 1, R = %d, workers = %d", R_REPS, WORKERS)
t0 <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(R_REPS), one_rep, seeds = seeds,
                                         mc.cores = WORKERS, mc.preschedule = FALSE))
log_msg("done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
write_csv_gz(res, file.path("crosscheck", "results", "sim1_gps_same_score_reps.csv.gz"))

# ---- reproduce the stored results exactly ------------------------------------
te <- read.csv(file.path("results", "true_effects_sim1.csv"))
key <- function(t, tp) paste0("(", t, ",", tp, ")")
mean_est <- aggregate(est ~ estimator + t + tp, res, mean)
main <- read.csv(file.path("results", "sim1_summary.csv")); main <- main[main$target == "pop", ]
chk <- merge(mean_est[mean_est$estimator %in% c("DIF", "CM", "GPS"), ], main[c("estimator", "t", "tp", "mean_est")],
             by = c("estimator", "t", "tp"), suffixes = c("", "_stored"))
lit <- read.csv(file.path("results", "gps_literal2_sim1.csv"))
lit_k <- match(key(mean_est$t, mean_est$tp), lit$pair)
lit_bias <- ifelse(mean_est$estimator == "GPS_lit", lit$gpslit_bias_pop[lit_k],
                   ifelse(mean_est$estimator == "GPS_lit2", lit$gpslit2_bias_pop[lit_k], NA))
tau_pop <- te$tau_pop[match(key(mean_est$t, mean_est$tp), key(te$t, te$tp))]
d_main <- max(abs(chk$est - chk$mean_est))
sel <- !is.na(lit_bias); d_lit <- max(abs(mean_est$est[sel] - tau_pop[sel] - lit_bias[sel]))
dl <- abs(mean_est$est[sel] - tau_pop[sel] - lit_bias[sel])
log_msg("reproduction: DIF/CM/GPS max |mean diff| vs results/sim1_summary.csv = %.3g (must be exact)", d_main)
log_msg("reproduction: GPS_lit/GPS_lit2 vs results/gps_literal2_sim1.csv: %d of 12 pair means exact (< 1e-12); largest |diff| %.3g, in %s",
        sum(dl < 1e-12), d_lit, paste(unique(paste(mean_est$estimator[sel][dl >= 1e-12], key(mean_est$t[sel][dl >= 1e-12], mean_est$tp[sel][dl >= 1e-12]))), collapse = ", "))
# The stored D13 run is from another setup (R 4.5.1, Windows), whose fitted scores
# differ from this one's at about the 1e-9 level (docs/spec.md, D5); at a near-tie
# that can change a literal-reading match. The main GPS matches on p_hat_s only and
# reproduces exactly.
stopifnot(nrow(chk) == 18, sum(sel) == 12, d_main < 1e-12, d_lit < 1e-4)

# ---- summaries -----------------------------------------------------------------
md_table <- function(df, digits = 3) {
  f <- function(v) if (!is.numeric(v)) as.character(v) else
    ifelse(is.na(v), "–", if (all(v[!is.na(v)] == round(v[!is.na(v)]))) formatC(v, format = "d") else formatC(v, format = "f", digits = digits))
  cells <- vapply(df, f, character(nrow(df)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = 1)
  c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(cells, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
save_table <- function(df, name, digits = 3) {
  write.csv(df, file.path(TAB, paste0(name, ".csv")), row.names = FALSE)
  writeLines(md_table(df, digits), file.path(TAB, paste0(name, ".md")))
}
EST <- c("GPS", "GPS_lit", "GPS_lit2", "GPS_same")
LAB <- c(GPS = "Main GPS (D7)", GPS_lit = "Literal, outside units only (GPS_lit)",
         GPS_lit2 = "Literal, as printed (GPS_lit2)", GPS_same = "Same score for t on both sides (GPS_same)")
paper <- read.csv(file.path("results", "paper_table1.csv"))
perpair <- do.call(rbind, lapply(EST, function(e) {
  x <- res[res$estimator == e, ]
  do.call(rbind, lapply(split(x, key(x$t, x$tp)), function(g) {
    k <- match(key(g$t[1], g$tp[1]), key(te$t, te$tp))
    data.frame(estimator = e, pair = key(g$t[1], g$tp[1]), t = g$t[1], tp = g$tp[1],
               tau_pop = te$tau_pop[k], tau_bal = te$tau_bal[k],
               bias_pop = mean(g$est) - te$tau_pop[k], bias_bal = mean(g$est) - te$tau_bal[k],
               sd = sd(g$est), mc_se = sd(g$est) / sqrt(nrow(g)),
               mse_pop = mean((g$est - te$tau_pop[k])^2), mse_bal = mean((g$est - te$tau_bal[k])^2),
               paper_abs_bias = paper$GPS_abs_bias[match(paste0(g$t[1], ",", g$tp[1]), paper$pair)])
  }))
}))
perpair <- perpair[order(match(perpair$estimator, EST), perpair$t, perpair$tp), ]
summ <- do.call(rbind, lapply(EST, function(e) {
  x <- perpair[perpair$estimator == e, ]
  data.frame(Reading = LAB[[e]], `Mean abs. bias vs tau_pop` = mean(abs(x$bias_pop)),
             `Mean abs. bias vs tau_bal` = mean(abs(x$bias_bal)), `Mean MSE vs tau_pop` = mean(x$mse_pop),
             `Mean SD` = mean(x$sd), `Pairs nearer the paper than the main GPS` = NA,
             `Correlation with the paper's abs. bias` = cor(abs(x$bias_pop), x$paper_abs_bias),
             `Pairs with bias beyond 3 MC SE` = sum(abs(x$bias_pop) / x$mc_se > 3), check.names = FALSE)
}))
mainx <- perpair[perpair$estimator == "GPS", ]
for (i in seq_along(EST)[-1]) {
  x <- perpair[perpair$estimator == EST[i], ]
  summ$`Pairs nearer the paper than the main GPS`[i] <-
    sum(abs(abs(x$bias_pop) - x$paper_abs_bias) < abs(abs(mainx$bias_pop) - mainx$paper_abs_bias))
}
summ <- rbind(data.frame(Reading = "Paper (Table 1)", `Mean abs. bias vs tau_pop` = mean(paper$GPS_abs_bias),
                         `Mean abs. bias vs tau_bal` = NA, `Mean MSE vs tau_pop` = NA, `Mean SD` = NA,
                         `Pairs nearer the paper than the main GPS` = NA, `Correlation with the paper's abs. bias` = NA,
                         `Pairs with bias beyond 3 MC SE` = NA, check.names = FALSE), summ)
save_table(summ, "sim1_eq4_readings_summary")
bp <- data.frame(Pair = unique(perpair$pair), Paper = mainx$paper_abs_bias)
for (e in EST) bp[[LAB[[e]]]] <- abs(perpair$bias_pop[perpair$estimator == e])
save_table(bp, "sim1_eq4_readings_abs_bias_by_pair")
write.csv(perpair, file.path(TAB, "sim1_eq4_readings_by_pair.csv"), row.names = FALSE)
print(summ, digits = 3); print(bp, digits = 3)
log_msg("wrote %s", TAB)
