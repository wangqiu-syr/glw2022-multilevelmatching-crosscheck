# gps_literal.R -- D13 sensitivity variant "GPS_lit" (docs/spec.md, D7 and D13).
# Requires R/dgp.R, R/estimators.R and R/replicate.R to be sourced.
#
# GPS_lit follows the paper's Eq. (4) as printed. For pair (t, t') and a unit i
# in neither group t nor group t':
#   t  term: Y of the nearest unit j in group t  on |p_hat_t(j)  - p_hat_t(i)|  (as in GPS)
#   t' term: Y of the nearest unit j in group t' on |p_hat_t'(j) - p_hat_t(i)|  (i's score for t)
# Units in group t or t' keep their observed Y on their own level; their other
# term is imputed exactly as in the main GPS estimator. The variant is
# pair-specific. It is a sensitivity check of the D7 reading, not a replacement.
#
# Pairing with the main run: estimate_with_lit() consumes random numbers in
# exactly the same order as estimate_all(), and GPS_lit's own tie-breaking runs
# with the RNG state saved and restored, so the main pipeline's samples,
# bootstrap draws and estimates are reproduced bit for bit.

# Imputed t' column for pair (t, tp) under the literal reading.
gps_lit_tp_column <- function(P, T, Y, Yimp, t, tp) {
  out <- Yimp[, tp]
  others <- which(T != t & T != tp)
  g <- which(T == tp)
  if (length(others)) {
    nn <- nn_pick(abs(outer(P[others, t], P[g, tp], "-")))
    out[others] <- Y[g[nn$j]]
  }
  out
}

# GPS_lit tau-hat for all pairs. Leaves the global RNG state unchanged.
tau_gps_lit <- function(P, T, Y, Yimp, A) {
  saved <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", saved, envir = globalenv()))
  pr <- pairs_of(A)
  v <- vapply(seq_len(nrow(pr)), function(k) {
    t <- pr[k, 1]; tp <- pr[k, 2]
    mean(Yimp[, t]) - mean(gps_lit_tp_column(P, T, Y, Yimp, t, tp))
  }, numeric(1))
  setNames(v, pair_labels(A))
}

# Same computations, in the same order, as estimate_all(), plus GPS_lit.
estimate_with_lit <- function(d, A, maxit = GPS_MAXIT) {
  cm  <- impute_cm(d, A)
  gps <- impute_gps(d, A, maxit)
  c(setNames(tau_dif(d, A),               paste0("DIF.",     pair_labels(A))),
    setNames(tau_from_imputed(cm$Yimp),   paste0("CM.",      pair_labels(A))),
    setNames(tau_from_imputed(gps$Yimp),  paste0("GPS.",     pair_labels(A))),
    setNames(tau_gps_lit(gps$P, d$T, d$Y, gps$Yimp, A), paste0("GPS_lit.", pair_labels(A))))
}

# Mirrors run_replication() + boot_ci() (stratified bootstrap, percentile CI).
run_replication_lit <- function(r, setup, B, seed, level = 0.95) {
  assign(".Random.seed", seed, envir = globalenv())
  A <- setup$A
  d <- gen_sample(setup)
  est <- estimate_with_lit(d, A)
  groups <- split(seq_len(nrow(d)), d$T)
  draws <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
  for (b in seq_len(B)) {
    idx <- unlist(lapply(groups, function(g) g[sample.int(length(g), replace = TRUE)]),
                  use.names = FALSE)
    draws[b, ] <- estimate_with_lit(d[idx, , drop = FALSE], A)
  }
  a <- (1 - level) / 2
  q <- apply(draws, 2, quantile, probs = c(a, 1 - a), names = FALSE)
  nm <- names(est)
  pr <- do.call(rbind, strsplit(sub("^[^.]*\\.", "", nm), "-"))
  data.frame(rep = r, estimator = sub("\\..*$", "", nm),
             t = as.integer(pr[, 1]), tp = as.integer(pr[, 2]),
             est = unname(est), lo = q[1, ], hi = q[2, ])
}
