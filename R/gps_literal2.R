# gps_literal2.R -- D13 sensitivity variant "GPS_lit2" (docs/spec.md, D7 and D13).
# Requires R/dgp.R, R/estimators.R, R/replicate.R and R/gps_literal.R to be sourced.
#
# GPS_lit2 is the fully literal reading of the paper's Eq. (4). For pair (t, t'):
#   t  term: units in group t keep their observed Y; every other unit gets the Y
#            of its nearest unit in group t on |p_hat_t(j) - p_hat_t(i)|  (as in GPS)
#   t' term: units in group t' keep their observed Y; every other unit, including
#            units in group t, gets the Y of its nearest unit in group t' on
#            |p_hat_t'(j) - p_hat_t(i)|  (i's score for t, not for t')
# Compared with GPS_lit (R/gps_literal.R), the only change is that units in
# group t now also use p_hat_t(i) for their t' term. Sensitivity check only.
#
# Pairing: estimate_with_lit2() consumes random numbers in the same order as
# estimate_all(); GPS_lit and GPS_lit2 each run with the RNG state saved and
# restored, so DIF/CM/GPS and GPS_lit reproduce the earlier runs exactly.

gps_lit2_tp_column <- function(P, T, Y, Yimp, t, tp) {
  out <- Yimp[, tp]
  others <- which(T != tp)                 # includes units in group t
  g <- which(T == tp)
  if (length(others)) {
    nn <- nn_pick(abs(outer(P[others, t], P[g, tp], "-")))
    out[others] <- Y[g[nn$j]]
  }
  out
}

# GPS_lit2 tau-hat for all pairs. Leaves the global RNG state unchanged.
tau_gps_lit2 <- function(P, T, Y, Yimp, A) {
  saved <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", saved, envir = globalenv()))
  pr <- pairs_of(A)
  v <- vapply(seq_len(nrow(pr)), function(k) {
    t <- pr[k, 1]; tp <- pr[k, 2]
    mean(Yimp[, t]) - mean(gps_lit2_tp_column(P, T, Y, Yimp, t, tp))
  }, numeric(1))
  setNames(v, pair_labels(A))
}

# Same computations, in the same order, as estimate_with_lit(), plus GPS_lit2.
estimate_with_lit2 <- function(d, A, maxit = GPS_MAXIT) {
  cm  <- impute_cm(d, A)
  gps <- impute_gps(d, A, maxit)
  c(setNames(tau_dif(d, A),               paste0("DIF.",      pair_labels(A))),
    setNames(tau_from_imputed(cm$Yimp),   paste0("CM.",       pair_labels(A))),
    setNames(tau_from_imputed(gps$Yimp),  paste0("GPS.",      pair_labels(A))),
    setNames(tau_gps_lit(gps$P, d$T, d$Y, gps$Yimp, A),  paste0("GPS_lit.",  pair_labels(A))),
    setNames(tau_gps_lit2(gps$P, d$T, d$Y, gps$Yimp, A), paste0("GPS_lit2.", pair_labels(A))))
}

# Mirrors run_replication_lit() (stratified bootstrap, percentile CI).
run_replication_lit2 <- function(r, setup, B, seed, level = 0.95) {
  assign(".Random.seed", seed, envir = globalenv())
  A <- setup$A
  d <- gen_sample(setup)
  est <- estimate_with_lit2(d, A)
  groups <- split(seq_len(nrow(d)), d$T)
  draws <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
  for (b in seq_len(B)) {
    idx <- unlist(lapply(groups, function(g) g[sample.int(length(g), replace = TRUE)]),
                  use.names = FALSE)
    draws[b, ] <- estimate_with_lit2(d[idx, , drop = FALSE], A)
  }
  a <- (1 - level) / 2
  q <- apply(draws, 2, quantile, probs = c(a, 1 - a), names = FALSE)
  nm <- names(est)
  pr <- do.call(rbind, strsplit(sub("^[^.]*\\.", "", nm), "-"))
  data.frame(rep = r, estimator = sub("\\..*$", "", nm),
             t = as.integer(pr[, 1]), tp = as.integer(pr[, 2]),
             est = unname(est), lo = q[1, ], hi = q[2, ])
}
