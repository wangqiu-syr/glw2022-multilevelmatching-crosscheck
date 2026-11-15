# dgp.R -- data-generating process for Simulations 1 and 2 (docs/spec.md).
#
# Main entry points:
#   sim_setup(sim)                 parameters for sim = 1 or 2
#   gen_population(n, setup)       n i.i.d. units with T and all potential outcomes
#   gen_sample(setup, n_g)         one balanced sample (D3): n_g units per group
#   true_tau_pop(setup)            closed-form population ATE (D4, primary)
#   mu_bal_mc(setup, n)            Monte Carlo covariate mean of the balanced design (D9)
#   true_tau(setup, m)             pairwise effects (beta_t - beta_t')' m, t < t'

COV_NAMES <- c("x1", "x2", "x3", "x4", "z1", "z2", "z3")

SIGMA_X <- matrix(c( 2.0,  1.0, -1.0,  0.5,
                     1.0,  1.0, -0.5,  0.2,
                    -1.0, -0.5,  1.0, -0.3,
                     0.5,  0.2, -0.3,  2.0), 4, 4, byrow = TRUE)

# ---- parameters -------------------------------------------------------------

# Treatment-score specification: for each t, which of the 7 covariates enter
# squared (1) or linearly (0). Row 1 is the baseline (eta_1 = 0) and is unused.
SQ_SIM2 <- rbind(
  c(NA, NA, NA, NA, NA, NA, NA),  # t = 1 baseline
  c(1, 0, 0, 0, 0, 0, 0),         # t = 2
  c(0, 1, 1, 0, 0, 0, 0),         # t = 3
  c(0, 1, 1, 1, 0, 0, 0),         # t = 4
  c(0, 1, 1, 1, 1, 0, 0),         # t = 5
  c(0, 1, 1, 1, 1, 1, 0),         # t = 6
  c(0, 1, 1, 1, 1, 1, 1),         # t = 7
  c(1, 1, 1, 1, 1, 1, 1),         # t = 8
  c(0, 0, 0, 0, 0, 0, 0)          # t = 9
)

BETA_SIM2 <- rbind(
  c(1, 1, 1, 1,   1, 1, 1),
  c(2, 3, 1, 2,   2, 2, 2),
  c(3, 1, 2, 1,   3, 1, 3),
  c(4, 1, 3, 4,   4, 1, 4),
  c(5, 1, 1, 1,   1, 1, 1),
  c(6, 3, 1, 2,   2, 2, 9),
  c(7, 1, 8, 1,   3, 1, 8),
  c(8, 1, 5, 2.5, 4, 1, 7),
  c(9, 1, 7, 1.5, 4, 1, 7)
)

sim_setup <- function(sim) {
  stopifnot(sim %in% c(1, 2))
  A <- if (sim == 1) 4L else 9L
  sq <- SQ_SIM2[1:A, , drop = FALSE]
  beta <- BETA_SIM2[1:A, , drop = FALSE]
  dimnames(sq) <- dimnames(beta) <- list(paste0("t", 1:A), COV_NAMES)
  list(
    sim       = sim,
    A         = A,
    sq        = sq,
    beta      = beta,
    eta_scale = 0.1,
    n_g       = if (sim == 1) 100L else 200L,
    # Outcome regressors W: Simulation 2 centres z2 and z3.
    w_shift   = if (sim == 1) rep(0, 7) else c(0, 0, 0, 0, 0, 1, 0.5)
  )
}

# ---- covariates, treatment, outcomes ---------------------------------------

gen_covariates <- function(n) {
  X <- MASS::mvrnorm(n, mu = rep(0, 4), Sigma = SIGMA_X)
  if (n == 1) X <- matrix(X, nrow = 1)
  V <- cbind(X,
             runif(n, -3, 3),        # z1 (D2)
             rchisq(n, df = 1),      # z2
             rbinom(n, 1, 0.5))      # z3
  colnames(V) <- COV_NAMES
  V
}

# n x A matrix of linear predictors eta_t (eta_1 = 0).
treatment_eta <- function(V, setup) {
  eta <- matrix(0, nrow(V), setup$A)
  V2 <- V^2
  for (t in seq_len(setup$A)[-1]) {
    s <- setup$sq[t, ]
    eta[, t] <- setup$eta_scale * (V2 %*% s + V %*% (1 - s))
  }
  eta
}

# n x A matrix of p_t(x, z): standard multinomial logit (D1), stabilised.
treatment_probs <- function(V, setup) {
  eta <- treatment_eta(V, setup)
  eta <- eta - apply(eta, 1, max)
  e <- exp(eta)
  P <- e / rowSums(e)
  colnames(P) <- paste0("t", seq_len(setup$A))
  P
}

# One categorical draw per row of P.
draw_treatment <- function(P) {
  cp <- P %*% upper.tri(diag(ncol(P)), diag = TRUE)  # row-wise cumsum
  u <- runif(nrow(P))
  pmin(rowSums(u > cp) + 1L, ncol(P))
}

outcome_regressors <- function(V, setup) sweep(V, 2, setup$w_shift)

# n x A matrix of potential outcomes; one epsilon per unit (spec).
potential_outcomes <- function(V, setup, eps = rnorm(nrow(V))) {
  Y <- outcome_regressors(V, setup) %*% t(setup$beta) + eps
  colnames(Y) <- paste0("Y", seq_len(setup$A))
  Y
}

# ---- populations and balanced samples --------------------------------------

gen_population <- function(n, setup, keep_probs = FALSE) {
  V <- gen_covariates(n)
  P <- treatment_probs(V, setup)
  T <- draw_treatment(P)
  Ypo <- potential_outcomes(V, setup)
  out <- list(V = V, T = T, Ypo = Ypo, Y = Ypo[cbind(seq_len(n), T)])
  if (keep_probs) out$P <- P
  out
}

# Balanced sample (D3): draw pool units in chunks, assign T from p_t, and keep
# the first n_g units of each group in draw order. Stops with an error if the
# groups are not filled within max_pool units.
gen_sample <- function(setup, n_g = setup$n_g, chunk = 20L * n_g * setup$A,
                       max_pool = 1e7, keep_po = FALSE) {
  A <- setup$A
  kept <- vector("list", A)
  have <- integer(A)
  pool <- 0
  while (any(have < n_g)) {
    if (pool >= max_pool) stop("gen_sample: groups not filled within max_pool")
    pop <- gen_population(chunk, setup)
    pool <- pool + chunk
    for (t in which(have < n_g)) {
      idx <- which(pop$T == t)
      idx <- idx[seq_len(min(length(idx), n_g - have[t]))]
      if (length(idx)) {
        kept[[t]] <- rbind(kept[[t]],
                           cbind(pop$V[idx, , drop = FALSE], T = t,
                                 Y = pop$Y[idx],
                                 pop$Ypo[idx, , drop = FALSE]))
        have[t] <- have[t] + length(idx)
      }
    }
  }
  d <- as.data.frame(do.call(rbind, kept))
  d$T <- as.integer(d$T)
  stopifnot(all(tabulate(d$T, A) == n_g))
  if (!keep_po) d <- d[, c(COV_NAMES, "T", "Y")]
  attr(d, "pool_size") <- pool
  d
}

# ---- true effects -----------------------------------------------------------

# Pairwise effects tau(t, t') = (beta_t - beta_t')' m for all t < t'.
true_tau <- function(setup, m) {
  pr <- t(combn(setup$A, 2))
  data.frame(t = pr[, 1], tp = pr[, 2],
             tau = as.vector((setup$beta[pr[, 1], ] - setup$beta[pr[, 2], ]) %*% m))
}

# Population mean of W in closed form: E[X] = 0, E[z1] = 0, E[z2] = 1, E[z3] = 0.5.
mu_pop <- function(setup) c(0, 0, 0, 0, 0, 1, 0.5) - setup$w_shift

true_tau_pop <- function(setup) true_tau(setup, mu_pop(setup))

# mu_bal = (1/A) sum_t E[W p_t] / E[p_t], estimated on n units drawn in chunks.
# Returns the estimate, its Monte Carlo SE, and its 7 x 7 Monte Carlo covariance
# (`vcov`), all from the spread of the chunk-level estimates (batch means).
mu_bal_mc <- function(setup, n = 1e6, n_chunks = 20) {
  m <- n / n_chunks
  A <- setup$A
  est <- matrix(NA, n_chunks, 7)
  num_all <- matrix(0, A, 7); den_all <- numeric(A)
  for (k in seq_len(n_chunks)) {
    V <- gen_covariates(m)
    W <- outcome_regressors(V, setup)
    P <- treatment_probs(V, setup)
    num <- crossprod(P, W) / m          # A x 7: E[p_t W]
    den <- colMeans(P)                  # E[p_t]
    est[k, ] <- colMeans(num / den)
    num_all <- num_all + num; den_all <- den_all + den
  }
  mu <- colMeans((num_all / n_chunks) / (den_all / n_chunks))
  names(mu) <- COV_NAMES
  vc <- cov(est) / n_chunks
  dimnames(vc) <- list(COV_NAMES, COV_NAMES)
  list(mu = mu, se = sqrt(diag(vc)), vcov = vc,
       group_prob = den_all / n_chunks)
}

# Monte Carlo SE of each pairwise effect (beta_t - beta_t')' m, given the
# covariance of the estimated m (full covariance, not just the diagonal).
true_tau_se <- function(setup, vcov_m) {
  pr <- t(combn(setup$A, 2))
  Dm <- setup$beta[pr[, 1], , drop = FALSE] - setup$beta[pr[, 2], , drop = FALSE]
  sqrt(rowSums((Dm %*% vcov_m) * Dm))
}
