# 06_existing_mismatches.R -- explain every imputed outcome where
# multiMatch(match_on = "existing") differs from our GPS matching.
# For each replication flagged in crosscheck/results/sim<k>_crosscheck_info.csv
# (cells_diff_ex > 0), regenerate the sample, refit "existing" exactly as
# 04_crosscheck_all.R does, and for each differing cell record our match, the
# package's match, both distances |p_hat_s(i) - p_hat_s(j)|, and the gap in
# squared standardized distance, the scale on which Matching::Match compares
# candidates. Matching treats candidates as tied when their squared distances
# differ by less than about 1e-10, even with distance.tolerance = 0, and with
# ties = FALSE it then picks one at random.
#   R_LIBS=~/R/mlm-lib Rscript crosscheck/R/06_existing_mismatches.R      (from the repo root)
# Writes crosscheck/results/existing_mismatches.csv.

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
suppressPackageStartupMessages(library(multilevelMatching))

rows <- list()
for (sim in 1:2) {
  info <- read.csv(file.path("crosscheck", "results", sprintf("sim%d_crosscheck_info.csv", sim)))
  reps <- info$rep[info$cells_diff_ex > 0]
  if (!length(reps)) next
  setup <- sim_setup(sim); A <- setup$A
  seeds <- make_seeds(nrow(info), 20220000 + sim)
  for (r in reps) {
    assign(".Random.seed", seeds[[r]], envir = globalenv())
    d <- gen_sample(setup)
    cm <- impute_cm(d, A); gps <- impute_gps(d, A); P <- gps$P
    P_ex <- P
    for (i in which(rowSums(P_ex) != 1)) {
      j <- which.max(P_ex[i, ]); P_ex[i, j] <- P_ex[i, j] + (1 - sum(P_ex[i, ])); k <- 0
      while (sum(P_ex[i, ]) != 1 && k < 50) {
        P_ex[i, j] <- P_ex[i, j] + sign(1 - sum(P_ex[i, ])) * .Machine$double.eps * P_ex[i, j] / 2; k <- k + 1
      }
    }
    Xe <- P_ex; colnames(Xe) <- paste0("p", 1:A); oldClass(Xe) <- "matrix"
    set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    f <- suppressMessages(multiMatch(d$Y, d$T, Xe, trimming = FALSE, match_on = "existing", M_matches = 1))
    pm <- unname(as.matrix(f$impute_mat)); om <- unname(gps$Yimp)
    w <- which(pm != om, arr.ind = TRUE)
    for (k in seq_len(nrow(w))) {
      i <- w[k, 1]; s <- w[k, 2]; g <- which(d$T == s)
      jo <- g[d$Y[g] == om[i, s]]; jp <- g[d$Y[g] == pm[i, s]]
      x <- P_ex[, s]; sdx <- sd(x)
      dpo <- abs(P_ex[i, s] - P_ex[jo, s]); dpp <- abs(P_ex[i, s] - P_ex[jp, s])
      rows[[length(rows) + 1]] <- data.frame(
        sim = sim, rep = r, unit = i, unit_T = d$T[i], level = s, our_match = jo, pkg_match = jp,
        dp_ours = dpo, dp_pkg = dpp, dp_gap = dpp - dpo,
        std_sq_gap = (dpp / sdx)^2 - (dpo / sdx)^2,
        ours_is_nearest = dpo == min(abs(P_ex[i, s] - P_ex[g, s])))
    }
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path("crosscheck", "results", "existing_mismatches.csv"), row.names = FALSE)
options(width = 200)
print(out, digits = 4)
cat(sprintf("\n%d differing cells in %d replications. Ours is the nearest candidate in all: %s.\n",
            nrow(out), nrow(unique(out[c("sim", "rep")])), all(out$ours_is_nearest)))
cat(sprintf("Gap in |dp|: max %.3g. Gap in squared standardized distance: max %.3g.\n",
            max(out$dp_gap), max(out$std_sq_gap)))
