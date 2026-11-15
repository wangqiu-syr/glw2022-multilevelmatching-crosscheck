# estimators.R -- DIF, CM and GPS estimators, stratified bootstrap and
# per-replication diagnostics (docs/spec.md, "Estimators" and D5-D8, D10-D11).
# Base R plus nnet::multinom (GPS model) only.
#
# A sample `d` is a data.frame with columns x1..x4, z1..z3 (COV_NAMES), T in
# 1..A and Y, as returned by gen_sample() in R/dgp.R.
#
# Main entry points:
#   impute_cm(d, A)            n x A imputed outcomes, Mahalanobis NN matching
#   impute_gps(d, A)           n x A imputed outcomes, matching on p_hat_s
#   tau_dif(d, A)              pairwise differences in group means
#   tau_from_imputed(Yimp)     pairwise tau-hat from column means
#   estimate_all(d, A)         all three estimators (+ diagnostics)
#   boot_ci(d, A, B)           stratified bootstrap percentile CIs

if (!exists("COV_NAMES")) COV_NAMES <- c("x1", "x2", "x3", "x4", "z1", "z2", "z3")

pairs_of <- function(A) t(combn(A, 2))

pair_labels <- function(A) {
  pr <- pairs_of(A)
  paste0(pr[, 1], "-", pr[, 2])
}

# ---- nearest-neighbour core -------------------------------------------------

# Dsub: n_query x n_target distance matrix. Returns, for each query row, the
# column index of a nearest target (ties broken uniformly at random) and the
# minimum distance.
nn_pick <- function(Dsub) {
  n <- nrow(Dsub)
  dmin <- apply(Dsub, 1, min)
  is_min <- Dsub == dmin                    # dmin recycles down columns = by row
  j <- max.col(is_min * 1, ties.method = "first")
  n_tie <- rowSums(is_min)
  for (r in which(n_tie > 1)) {
    cand <- which(is_min[r, ])
    j[r] <- cand[sample.int(length(cand), 1L)]
  }
  list(j = j, d = dmin, n_tie = n_tie)
}

# Generic imputation given a distance function dist_to(g, s) that returns the
# n x |g| distance matrix from every unit to the units g in group s.
# Own-group units keep their observed Y. Returns imputed matrix and mean
# matching distance (over units not in s) per target group.
impute_by_matching <- function(T, Y, A, dist_to) {
  n <- length(Y)
  Yimp <- matrix(NA_real_, n, A, dimnames = list(NULL, paste0("Y", 1:A)))
  mdist <- numeric(A)
  for (s in seq_len(A)) {
    g <- which(T == s)
    if (!length(g)) stop("impute_by_matching: group ", s, " is empty")
    nn <- nn_pick(dist_to(g, s))
    Yimp[, s] <- Y[g[nn$j]]
    own <- T == s
    Yimp[own, s] <- Y[own]
    mdist[s] <- mean(nn$d[!own])
  }
  list(Yimp = Yimp, match_dist = mdist)
}

# ---- CM: Mahalanobis matching on (X, Z) -------------------------------------

# Mahalanobis squared distances between all rows of M, with V = cov(M).
mahalanobis_d2 <- function(M) {
  R <- chol(solve(cov(M)))                  # V^-1 = R'R
  Mt <- M %*% t(R)                          # whitened rows
  sq <- rowSums(Mt^2)
  D2 <- outer(sq, sq, "+") - 2 * tcrossprod(Mt)
  D2[D2 < 0] <- 0
  D2
}

impute_cm <- function(d, A) {
  M <- as.matrix(d[, COV_NAMES])
  D <- sqrt(mahalanobis_d2(M))
  impute_by_matching(d$T, d$Y, A, function(g, s) D[, g, drop = FALSE])
}

# ---- GPS: multinomial logit + matching on p_hat_s ---------------------------

gps_formula <- T_f ~ x1 + x2 + x3 + x4 + z1 + z2 + z3

# Optimizer settings. GPS matching distances in Simulation 2 are ~7e-4, while
# nnet's default reltol = 1e-8 leaves p_hat uncertain at ~1e-3, so matches
# (and tau-hat) depended on the starting point. At reltol = 1e-14 the fit is
# start-independent: p_hat agrees to ~5e-8 and tau-hat exactly (see Test 4b).
# multinom always starts from zero weights (rang = 0), so fitting uses no RNG.
GPS_RELTOL <- 1e-14
GPS_MAXIT  <- 10000

fit_gps <- function(d, A, maxit = GPS_MAXIT, reltol = GPS_RELTOL) {
  dd <- d
  dd$T_f <- factor(d$T, levels = seq_len(A))
  fit <- nnet::multinom(gps_formula, data = dd, trace = FALSE,
                        maxit = maxit, reltol = reltol)
  P <- fitted(fit)
  if (is.null(dim(P))) P <- cbind(1 - P, P) # A = 2 safeguard
  colnames(P) <- paste0("t", seq_len(A))
  list(P = unname(P), converged = fit$convergence == 0, fit = fit)
}

# Matching on the scores in P (n x A): column s uses p_hat_s (so the t' term
# of each pair uses p_hat_t', D7).
impute_from_scores <- function(P, T, Y, A) {
  impute_by_matching(T, Y, A, function(g, s) abs(outer(P[, s], P[g, s], "-")))
}

impute_gps <- function(d, A, maxit = GPS_MAXIT) {
  gf <- fit_gps(d, A, maxit)
  out <- impute_from_scores(gf$P, d$T, d$Y, A)
  out$P <- gf$P
  out$converged <- gf$converged
  out
}

# ---- point estimates ----------------------------------------------------------

tau_from_imputed <- function(Yimp) {
  A <- ncol(Yimp); pr <- pairs_of(A)
  m <- colMeans(Yimp)
  setNames(m[pr[, 1]] - m[pr[, 2]], pair_labels(A))
}

tau_dif <- function(d, A) {
  pr <- pairs_of(A)
  m <- vapply(seq_len(A), function(s) mean(d$Y[d$T == s]), numeric(1))
  setNames(m[pr[, 1]] - m[pr[, 2]], pair_labels(A))
}

# All three estimators. Returns a named vector est (DIF.1-2, ..., CM.., GPS..)
# and, if diagnostics = TRUE, a list of diagnostics.
estimate_all <- function(d, A, diagnostics = FALSE, maxit = GPS_MAXIT) {
  cm  <- impute_cm(d, A)
  gps <- impute_gps(d, A, maxit)
  est <- c(setNames(tau_dif(d, A),           paste0("DIF.", pair_labels(A))),
           setNames(tau_from_imputed(cm$Yimp),  paste0("CM.",  pair_labels(A))),
           setNames(tau_from_imputed(gps$Yimp), paste0("GPS.", pair_labels(A))))
  out <- list(est = est, gps_converged = gps$converged)
  if (diagnostics) {
    out$diag <- list(
      cm_match_dist   = cm$match_dist,              # mean Mahalanobis distance, per target group
      gps_match_dist  = gps$match_dist,             # mean |p_hat_s difference|, per target group
      p_hat_below_01  = colMeans(gps$P < 0.01),     # share of units with p_hat_s < 0.01, per group
      p_hat_any_below_01 = mean(apply(gps$P, 1, min) < 0.01),
      gps_converged   = gps$converged
    )
  }
  out
}

# ---- stratified bootstrap -----------------------------------------------------

# Resample n_t units with replacement within each treatment group; recompute
# every estimator (V and the GPS model are re-estimated on each draw, D8).
boot_draws <- function(d, A, B, maxit = GPS_MAXIT) {
  groups <- split(seq_len(nrow(d)), d$T)
  draws <- NULL; n_nonconv <- 0L
  for (b in seq_len(B)) {
    idx <- unlist(lapply(groups, function(g) g[sample.int(length(g), replace = TRUE)]),
                  use.names = FALSE)
    e <- estimate_all(d[idx, , drop = FALSE], A, maxit = maxit)
    if (is.null(draws)) draws <- matrix(NA_real_, B, length(e$est),
                                        dimnames = list(NULL, names(e$est)))
    draws[b, ] <- e$est
    n_nonconv <- n_nonconv + !e$gps_converged
  }
  attr(draws, "n_nonconverged") <- n_nonconv
  draws
}

boot_ci <- function(d, A, B, level = 0.95, maxit = GPS_MAXIT) {
  draws <- boot_draws(d, A, B, maxit)
  a <- (1 - level) / 2
  q <- apply(draws, 2, quantile, probs = c(a, 1 - a), names = FALSE)
  list(lo = q[1, ], hi = q[2, ], n_nonconverged = attr(draws, "n_nonconverged"))
}
