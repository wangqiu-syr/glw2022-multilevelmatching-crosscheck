# 02_pkg_fits.R -- cross-check replication 1 of Simulation 1 against
# multilevelMatching::multiMatch. Reads only the exported CSVs (plus our
# internal RDS for diagnostics). Writes only to <out_dir>.
#   R_LIBS=<lib> Rscript crosscheck/R/02_pkg_fits.R crosscheck/results   (from the repo root)

args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]
suppressPackageStartupMessages(library(multilevelMatching))

dat  <- read.csv(file.path(OUT, "sim1_rep001.csv"))
ours <- read.csv(file.path(OUT, "sim1_rep001_ours.csv"))
P    <- as.matrix(read.csv(file.path(OUT, "sim1_rep001_gps.csv")))
int  <- readRDS(file.path(OUT, "rep001_ours_internal.rds"))
COV  <- c("x1", "x2", "x3", "x4", "z1", "z2", "z3")

Y <- dat$Y; W <- dat$T
X <- as.matrix(dat[, COV]); dimnames(P) <- list(NULL, colnames(P))
stopifnot(!is.unsorted(W))                     # data are sorted by T: no reordering inside multiMatch

# --- workaround 1: R >= 4.3 vs prepareData()'s class check --------------------
# prepareData() runs `(length(class(X)) != 1 && class(X) != "matrix") || ...`.
# For a plain matrix class(X) is c("matrix", "array") (R >= 4.0), and since
# R 4.3 `&&` on a length-2 operand is an error. Setting the class attribute
# explicitly makes class(X) == "matrix" (length 1) without changing any value.
as_mat1 <- function(M) { oldClass(M) <- "matrix"; M }

# --- workaround 2: "existing" requires rowSums(X) == 1 exactly ---------------
# Our fitted P sums to 1 only up to rounding in some rows. Absorb the residual
# into the largest entry of each failing row (a change of a few ulp), retrying
# with one-ulp nudges until the floating-point row sum is exactly 1.
exact_rows <- function(P) {
  for (i in which(rowSums(P) != 1)) {
    j <- which.max(P[i, ])
    P[i, j] <- P[i, j] + (1 - sum(P[i, ]))
    k <- 0
    while (sum(P[i, ]) != 1 && k < 50) {
      P[i, j] <- P[i, j] + sign(1 - sum(P[i, ])) * .Machine$double.eps * P[i, j] / 2
      k <- k + 1
    }
  }
  P
}
P_ex <- exact_rows(P)
cat(sprintf("GPS rows summing to exactly 1: before %d/%d, after %d/%d; max |change| = %.2e\n",
            sum(rowSums(P) == 1), nrow(P), sum(rowSums(P_ex) == 1), nrow(P),
            max(abs(P_ex - P))))
stopifnot(all(rowSums(P_ex) == 1))

# Does passing P unadjusted fail? (record the message)
msg_raw <- tryCatch({ set.seed(123)
  multiMatch(Y, W, as_mat1(P), trimming = FALSE, match_on = "existing"); "no error" },
  error = function(e) conditionMessage(e))
cat("existing with unadjusted P:", gsub("\\s+", " ", msg_raw), "\n")
msg_cls <- tryCatch({ set.seed(123)
  multiMatch(Y, W, X, trimming = FALSE, match_on = "covariates"); "no error" },
  error = function(e) conditionMessage(e))
cat("covariates with a plain matrix X:", msg_cls, "\n")

# --- the three fits -----------------------------------------------------------
fit <- list()
set.seed(123); fit$cov <- multiMatch(Y, W, as_mat1(X), trimming = FALSE,
                                     match_on = "covariates", M_matches = 1)
set.seed(123); fit$mn  <- multiMatch(Y, W, as_mat1(X), trimming = FALSE,
                                     match_on = "multinom", M_matches = 1)
set.seed(123); fit$ex  <- multiMatch(Y, W, as_mat1(P_ex), trimming = FALSE,
                                     match_on = "existing", M_matches = 1)

cat("\nOutput fields of a multiMatch object:\n"); print(names(fit$cov))
cat("Columns of $results:\n"); print(names(fit$cov$results))
cat("\n$results (covariates):\n"); print(fit$cov$results, digits = 10)
cat("\n$mu (covariates):\n"); print(fit$cov$mu, digits = 10)
cat("\nanalysis_idx (NULL = no trimming):", is.null(fit$cov$analysis_idx),
    is.null(fit$mn$analysis_idx), is.null(fit$ex$analysis_idx), "\n")

# Sign convention: Estimate = mu[Trt2] - mu[Trt1]
mu_of <- function(f) setNames(f$mu$Estimate, f$mu$Trt)
for (nm in names(fit)) {
  r <- fit[[nm]]$results; mu <- mu_of(fit[[nm]])
  stopifnot(isTRUE(all.equal(r$Estimate, unname(mu[as.character(r$Trt2)] - mu[as.character(r$Trt1)]),
                             tolerance = 0)))
}
cat("Checked: Estimate == EY(Trt2) - EY(Trt1) in all three fits\n")

# --- table in our convention: tau(t, t') = E[Y(t) - Y(t')] = -pkg -------------
flip <- function(f) {
  r <- f$results
  key <- paste0(r$Trt1, "-", r$Trt2)
  setNames(-r$Estimate, key)[ours$pair]
}
tab <- data.frame(pair = paste0("(", ours$t, ",", ours$tp, ")"),
                  ours_CM = ours$CM, pkg_CM = unname(flip(fit$cov)),
                  ours_GPS = ours$GPS, pkg_GPS_multinom = unname(flip(fit$mn)),
                  pkg_GPS_existing = unname(flip(fit$ex)))
tab$d_CM <- tab$pkg_CM - tab$ours_CM
tab$d_GPS_multinom <- tab$pkg_GPS_multinom - tab$ours_GPS
tab$d_GPS_existing <- tab$pkg_GPS_existing - tab$ours_GPS
cat("\nComparison, our sign convention tau(t,t') = E[Y(t) - Y(t')]:\n")
print(tab, digits = 10)
write.csv(tab, file.path(OUT, "sim1_rep001_crosscheck.csv"), row.names = FALSE)

# --- diagnostics: imputed potential outcomes, unit by unit --------------------
own <- cbind(seq_along(W), W)
cmp_imp <- function(pkg_mat, our_mat) {
  pkg_mat <- unname(as.matrix(pkg_mat)); our_mat <- unname(our_mat)
  stopifnot(all(pkg_mat[own] == Y), all(our_mat[own] == Y))
  colSums(pkg_mat != our_mat)
}
cat("\nUnits whose imputed Y(s) differs from ours, by level s = 1..4:\n")
cat("  covariates :", cmp_imp(fit$cov$impute_mat, int$Yimp_cm), "\n")
cat("  multinom   :", cmp_imp(fit$mn$impute_mat,  int$Yimp_gps), "\n")
cat("  existing   :", cmp_imp(fit$ex$impute_mat,  int$Yimp_gps), "\n")

# Package's multinom fit vs ours.
Pm <- unname(fit$mn$propensity_scores)
cat(sprintf("\nmultinom: package fit convergence code %d (0 = converged), iterations cap maxit = %s\n",
            fit$mn$model$convergence, "100 (nnet default)"))
cat(sprintf("  max |p_hat_pkg - p_hat_ours| = %.3e; median = %.3e\n",
            max(abs(Pm - unname(P))), median(abs(Pm - unname(P)))))
cat(sprintf("  deviance: package %.10f\n", fit$mn$model$deviance))
# Our deviance for reference (refit with our settings).
dd <- data.frame(X, T_f = factor(W))
f_ours <- nnet::multinom(T_f ~ x1 + x2 + x3 + x4 + z1 + z2 + z3, data = dd, trace = FALSE,
                         maxit = 10000, reltol = 1e-14)
f_def  <- nnet::multinom(T_f ~ x1 + x2 + x3 + x4 + z1 + z2 + z3, data = dd, trace = FALSE)
cat(sprintf("  deviance: ours (reltol 1e-14) %.10f; nnet defaults %.10f (convergence %d)\n",
            f_ours$deviance, f_def$deviance, f_def$convergence))
cat(sprintf("  package P == nnet-defaults P: max diff %.3e\n",
            max(abs(Pm - unname(fitted(f_def))))))
cat(sprintf("  our refit P == exported P: max diff %.3e\n",
            max(abs(unname(fitted(f_ours)) - unname(P)))))

# Our matching rule applied to the package's multinom P reproduces the package?
nn_imp <- function(PP) {
  Yi <- matrix(NA_real_, length(Y), 4)
  for (s in 1:4) {
    g <- which(W == s)
    D <- abs(outer(PP[, s], PP[g, s], "-"))
    Yi[, s] <- Y[g[max.col(-D, ties.method = "first")]]
    Yi[W == s, s] <- Y[W == s]
  }
  Yi
}
tau_of <- function(Yi) { m <- colMeans(Yi); pr <- t(combn(4, 2)); m[pr[, 1]] - m[pr[, 2]] }
cat(sprintf("  our NN rule on package P vs package multinom GPS: max |diff| = %.3e\n",
            max(abs(tau_of(nn_imp(Pm)) - tab$pkg_GPS_multinom))))
cat(sprintf("  our NN rule on P_ex vs our GPS: max |diff| = %.3e\n",
            max(abs(tau_of(nn_imp(P_ex)) - tab$ours_GPS))))

# Distance gap to the runner-up match in GPS matching (how close to a tie?)
gap <- unlist(lapply(1:4, function(s) {
  g <- which(W == s); D <- abs(outer(P[W != s, s], P[g, s], "-"))
  apply(D, 1, function(r) { o <- sort(r, partial = 1:2)[1:2]; o[2] - o[1] })
}))
cat(sprintf("  GPS matching, min gap between best and second-best |dp|: %.3e (median %.3e)\n",
            min(gap), median(gap)))

saveRDS(fit, file.path(OUT, "rep001_pkg_fits.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "pkg_session_info.txt"))
cat("\nmultilevelMatching", format(packageVersion("multilevelMatching")),
    "| Matching", format(packageVersion("Matching")),
    "| nnet", format(packageVersion("nnet")), "|", R.version.string, "\n")
