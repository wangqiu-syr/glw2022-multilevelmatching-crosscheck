# 03_check_pasted.R -- tests an earlier draft of a cross-check script
# (see docs/crosscheck.md, Appendix A.1): plain-matrix X fails on R >= 4.3, fit$estimate partially
# matches estimate_args, and 0-based treatment coding gives the same estimates.
#   R_LIBS=<lib> Rscript crosscheck/R/03_check_pasted.R crosscheck/results   (from the repo root)
suppressPackageStartupMessages(library(multilevelMatching))
OUT <- commandArgs(trailingOnly = TRUE)[1]
d <- read.csv(file.path(OUT, "sim1_rep001.csv"))
Y <- d$Y; W <- as.integer(d$T)
X <- as.matrix(d[, c("x1", "x2", "x3", "x4", "z1", "z2", "z3")])

cat("== A. Pasted script as written (plain matrix X):\n")
r <- tryCatch({ set.seed(123); multiMatch(Y = Y, W = W, X = X, match_on = "covariates", trimming = FALSE) },
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)

Xc <- X; oldClass(Xc) <- "matrix"
set.seed(123); fit_cm <- suppressMessages(multiMatch(Y = Y, W = W, X = Xc, match_on = "covariates", trimming = FALSE))

cat("\n== B. Pasted extraction lines (with the class workaround so the fit runs):\n")
x <- fit_cm$estimate
cat("fit_cm$estimate partially matches:", if (is.list(x)) "estimate_args (a list of inputs)" else class(x), "\n")
pe <- tryCatch(dim(as.data.frame(fit_cm$estimate)), error = function(e) paste("ERROR:", conditionMessage(e)))
cat("as.data.frame(fit_cm$estimate):", if (is.character(pe)) pe else paste("data frame", paste(pe, collapse = " x ")), "\n")
s <- capture.output(sv <- summary(fit_cm))
cat("summary(fit_cm) returns:", if (is.null(sv)) "NULL (prints only)" else class(sv), "-> $estimates is",
    if (is.null(sv$estimates)) "NULL" else "non-NULL", "\n")
cat("Correct field: fit_cm$results\n"); print(fit_cm$results[, c("Param", "Estimate")], digits = 10)

cat("\n== C. Treatment coding 0..3 instead of 1..4:\n")
set.seed(123); f0 <- suppressMessages(multiMatch(Y = Y, W = W - 1L, X = Xc, match_on = "covariates", trimming = FALSE))
print(f0$results[, c("Param", "Estimate")], digits = 10)
cat("Same estimates as 1..4 coding:", identical(f0$results$Estimate, fit_cm$results$Estimate), "\n")
