# run_all.R -- full run: Simulation 1, then Simulation 2, then summaries.
# Run from the repository root, e.g. in the background from PowerShell:
#   Start-Process Rscript.exe -ArgumentList "R/run_all.R" `
#     -RedirectStandardOutput results/run_log.txt -RedirectStandardError results/run_log_stderr.txt
# Optional key=value arguments (defaults in brackets):
#   R=[500] B=[200] workers=[10] rep_root=[results/reps] out=[results]
# Re-running resumes from the checkpoint files under rep_root.

source("R/dgp.R"); source("R/estimators.R"); source("R/replicate.R")
source("R/run_sim.R"); source("R/summarize.R")

opt <- list(R = "500", B = "200", workers = "10", rep_root = file.path("results", "reps"),
            out = "results")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2 && kv[1] %in% names(opt)) opt[[kv[1]]] <- kv[2] else stop("unknown argument: ", a)
}
R_REPS <- as.integer(opt$R); B <- as.integer(opt$B); WORKERS <- as.integer(opt$workers)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

log_msg("run_all start: R = %d, B = %d, workers = %d, rep_root = %s, out = %s",
        R_REPS, B, WORKERS, opt$rep_root, opt$out)
log_msg("%s; nnet %s; MASS %s; GPS reltol %g, maxit %d; %d cores detected",
        R.version.string, packageVersion("nnet"), packageVersion("MASS"),
        GPS_RELTOL, GPS_MAXIT, parallel::detectCores())

t_all <- proc.time()[["elapsed"]]
ok <- TRUE
for (sim in 1:2) {
  ok <- run_sim(sim, R = R_REPS, B = B, workers = WORKERS,
                rep_dir = file.path(opt$rep_root, paste0("sim", sim))) && ok
}
for (sim in 1:2) {
  summarize_sim(sim, rep_dir = file.path(opt$rep_root, paste0("sim", sim)), out_dir = opt$out)
}
compare_paper(summary_file = file.path(opt$out, "sim1_summary.csv"),
              out_file = file.path(opt$out, "sim1_vs_paper.csv"))
# The workers load nnet and MASS, not this session, so load their namespaces
# here (without attaching) so that sessionInfo() records their versions.
invisible(lapply(c("nnet", "MASS"), loadNamespace))
writeLines(capture.output(sessionInfo()), file.path(opt$out, "run_session_info.txt"))
log_msg("run_all finished in %s%s", fmt_dur(proc.time()[["elapsed"]] - t_all),
        if (ok) "" else " (with failed replications; see *_error.txt)")
