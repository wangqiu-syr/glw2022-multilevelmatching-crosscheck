# 01_export_rep001.R -- regenerate replication 1 of Simulation 1 exactly as the
# main run did (master seed 20220001, L'Ecuyer stream 1), export data, our
# CM/GPS estimates and our fitted GPS matrix. Reads the repo; writes only to
# <out_dir>.
#   Rscript crosscheck/R/01_export_rep001.R . crosscheck/results   (from the repo root)

args <- commandArgs(trailingOnly = TRUE)
REPO <- normalizePath(args[1]); OUT <- args[2]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

src <- function(f) source(file.path(REPO, "R", f), local = globalenv())
src("dgp.R"); src("estimators.R"); src("replicate.R")

setup <- sim_setup(1); A <- setup$A
seeds <- make_seeds(500, 20220000 + 1)          # as run_sim(): master_seed = 20220000 + sim

# --- exactly the start of run_replication(r = 1, ...) -------------------------
assign(".Random.seed", seeds[[1]], envir = globalenv())
d <- gen_sample(setup)
state_after_sample <- .Random.seed
pt <- estimate_all(d, A, diagnostics = TRUE)    # same call as run_replication()

# Re-run the two imputations from the same RNG state to keep Yimp and P
# (estimate_all calls impute_cm, then impute_gps; nothing else draws RNG).
assign(".Random.seed", state_after_sample, envir = globalenv())
cm  <- impute_cm(d, A)
gps <- impute_gps(d, A)
stopifnot(identical(unname(tau_from_imputed(cm$Yimp)),
                    unname(pt$est[paste0("CM.", pair_labels(A))])),
          identical(unname(tau_from_imputed(gps$Yimp)),
                    unname(pt$est[paste0("GPS.", pair_labels(A))])))

# Ties in our own matching (would make results RNG-dependent)?
M <- as.matrix(d[, COV_NAMES]); Dm <- sqrt(mahalanobis_d2(M))
n_ties <- function(dist_to) sum(vapply(seq_len(A), function(s) {
  g <- which(d$T == s); nn <- nn_pick(dist_to(g, s)); sum(nn$n_tie[d$T != s] > 1)
}, numeric(1)))
cat(sprintf("Units with tied nearest neighbours: CM %d, GPS %d\n",
            n_ties(function(g, s) Dm[, g, drop = FALSE]),
            n_ties(function(g, s) abs(outer(gps$P[, s], gps$P[g, s], "-")))))

# --- check against the main run's stored rep-1 diagnostics --------------------
dg <- read.csv(file.path(REPO, "results", "sim1_diagnostics.csv"))
ref <- dg[dg$rep == 1, ]
mine <- c(pool_size = attr(d, "pool_size"),
          setNames(pt$diag$cm_match_dist,  paste0("cm_dist_t", 1:A)),
          setNames(pt$diag$gps_match_dist, paste0("gps_dist_t", 1:A)),
          setNames(pt$diag$p_hat_below_01, paste0("p_hat_below_01_t", 1:A)),
          p_hat_any_below_01 = pt$diag$p_hat_any_below_01,
          gps_converged = as.numeric(pt$diag$gps_converged))
refv <- unlist(ref[names(mine)]) + 0
chk <- data.frame(item = names(mine), main_run = refv, this_run = mine,
                  abs_diff = abs(refv - mine), row.names = NULL)
print(chk, digits = 15)
cat(sprintf("Max abs diff vs stored rep-1 diagnostics: %.3g\n", max(chk$abs_diff)))
stopifnot(max(chk$abs_diff) < 1e-12)

# --- exports (17 significant digits so doubles round-trip) --------------------
f17 <- function(x) sprintf("%.17g", x)
write_exact <- function(df, file) {
  out <- df
  for (j in names(out)) if (is.double(out[[j]])) out[[j]] <- f17(out[[j]])
  write.csv(out, file, row.names = FALSE, quote = FALSE)
}

dat <- data.frame(T = d$T, Y = d$Y, d[, COV_NAMES])
write_exact(dat, file.path(OUT, "sim1_rep001.csv"))

pr <- pairs_of(A)
ours <- data.frame(pair = pair_labels(A), t = pr[, 1], tp = pr[, 2],
                   CM  = unname(pt$est[paste0("CM.",  pair_labels(A))]),
                   GPS = unname(pt$est[paste0("GPS.", pair_labels(A))]),
                   DIF = unname(pt$est[paste0("DIF.", pair_labels(A))]))
write_exact(ours, file.path(OUT, "sim1_rep001_ours.csv"))

P <- gps$P; colnames(P) <- paste0("p", 1:A)
write_exact(as.data.frame(P), file.path(OUT, "sim1_rep001_gps.csv"))

# Round-trip check: values read back are bit-identical.
rt <- read.csv(file.path(OUT, "sim1_rep001.csv"))
rp <- as.matrix(read.csv(file.path(OUT, "sim1_rep001_gps.csv")))
ro <- read.csv(file.path(OUT, "sim1_rep001_ours.csv"))
cat("Round trip identical: data", identical(unname(as.matrix(rt)), unname(as.matrix(dat))),
    "| gps", identical(unname(rp), unname(P)),
    "| ours", identical(ro$CM, ours$CM) && identical(ro$GPS, ours$GPS), "\n")

# Extra objects for the comparison step (scratch only).
saveRDS(list(d = d, est = pt$est, diag = pt$diag, Yimp_cm = cm$Yimp,
             Yimp_gps = gps$Yimp, P = gps$P, gps_fit = gps$converged,
             coef = coef(fit_gps(d, A)$fit)),
        file.path(OUT, "rep001_ours_internal.rds"))

cat("multinom converged:", gps$converged, "\n")
cat("rowSums(P) == 1 exactly:", sum(rowSums(P) == 1), "of", nrow(P), "rows\n")
print(ours, digits = 10)
cat(R.version.string, "| nnet", format(packageVersion("nnet")),
    "| MASS", format(packageVersion("MASS")), "| LAPACK", La_version(), "\n")
