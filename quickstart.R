#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# malariafwd quickstart
#
# Installs the package, runs one simulation, and writes the diagnostic and model-output
# plots. Use it to check a fresh install end to end, or as a worked example to copy from.
#
#   Rscript quickstart.R                 # run from a clone: installs from this directory
#   Rscript quickstart.R --github        # install from GitHub instead
#   Rscript quickstart.R --outdir path   # where to put the run and the figures
#
# Needs SLiM (>= 5.2) on your PATH; the IBD and ARG steps also need Python with tskit (>= 1.0),
# pandas and numpy, through reticulate (with no Python chosen, reticulate sets one up on first
# use). check_setup() reports what is found. Everything else degrades gracefully.
# ---------------------------------------------------------------------------------------

args    <- commandArgs(trailingOnly = TRUE)
use_gh  <- "--github" %in% args
outdir  <- if (any(args == "--outdir")) args[which(args == "--outdir") + 1] else "quickstart-output"
figs    <- file.path(outdir, "figures")
rundir  <- file.path(outdir, "run")
dir.create(figs, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat("\n", paste0(...), "\n", sep = "")
hr  <- function() cat(strrep("-", 78), "\n")

# -- 1. install ---------------------------------------------------------------------------
say("1. Installing malariafwd")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
if (use_gh) {
  remotes::install_github("shaziaruybal/malariafwd", upgrade = "never", dependencies = TRUE)
} else {
  here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
  if (!file.exists(file.path(here, "DESCRIPTION")))
    stop("no DESCRIPTION next to this script - run it from a clone, or pass --github", call. = FALSE)
  remotes::install_local(here, force = TRUE, upgrade = "never", dependencies = TRUE)
}
suppressPackageStartupMessages(library(malariafwd))
cat("   installed version:", as.character(utils::packageVersion("malariafwd")), "\n")

# -- 2. check the toolchain ---------------------------------------------------------------
say("2. Checking the toolchain")
print(check_setup())
if (!slim_available())
  stop("SLiM is not on your PATH. Install SLiM >= 5.2, or set MALARIAFWD_SLIM to the binary.",
       call. = FALSE)
has_py <- tskit_available()
if (!has_py) cat("   the IBD and ARG steps will be skipped (no Python with tskit: see above)\n")

# -- 3. run one simulation ----------------------------------------------------------------
say("3. Running a simulation")
p <- vivax_params(
  N_HOSTS = 100, N_MOSQ = 200,     # a small population, so this finishes quickly
  BITING_RATE = 0.11,              # moderate transmission
  SIM_DAYS = 548,                  # 1.5 years
  N_FOUNDER_HOSTS = 16,
  N_FOUNDER_LINEAGES = 10,         # 10 distinct founding genomes among 16 infected hosts
  LAMBDA_SPOROZOITES = 2           # >1 genotype can be co-transmitted in a bite
)
print(p)
unlink(rundir, recursive = TRUE)
t0  <- Sys.time()
run <- run_sim(p, seed = 1, outdir = rundir)
cat("   done in", round(as.numeric(Sys.time() - t0, units = "secs")), "s ->", rundir, "\n")

# -- 4. what the run produced --------------------------------------------------------------
say("4. Run summary")
ev <- load_events(run)
st <- load_state(run)
hr(); print(summarise_run(run)); hr()

cat("\nevents logged:\n")
print(ev[, .N, by = event_type][order(-N)])
cat("\ncases by true origin:\n")
cs <- classify_cases(ev)
print(cs[, .N, by = .(origin, imported)][order(origin)])
cat("\nforce of blood-stage infection (first windows):\n")
print(head(force_of_infection(run), 6))

# -- 5. plots -------------------------------------------------------------------------------
say("5. Writing plots to ", figs)
sv <- function(name, plot, w = 9, h = 5) {
  f <- file.path(figs, paste0(name, ".png"))
  ggplot2::ggsave(f, plot, width = w, height = h, dpi = 140)
  cat("   ", basename(f), "\n")
}

# 5a. the run-health dashboard: is this simulation behaving sensibly?
sv("01-diagnostics", plot_diagnostics(run), w = 11, h = 8)

# 5b. epidemiological output
host <- pick_host(ev)
sv("02-host-trajectory", plot_host_trajectory(st, ev, host_id = host), w = 9, h = 4)
sv("03-host-states",     plot_host_state_prevalence(host_state_prevalence(st)))
sv("04-time-to-relapse", plot_relapse_km(time_to_relapse(run)))
sv("05-moi-vs-relapse",  plot_moi_vs_relapse(moi_vs_relapse(st, ev)))

# 5c. genetic output (needs tskit)
if (has_py) {
  say("   computing ground-truth IBD (this is the slow step)")
  pairs <- ibd_truth_pairs(run)
  cat("\nIBD by known relationship:\n")
  print(pairs[, .(n = .N, median_ibd = round(median(ibd), 3)), by = cls][order(cls)])

  sv("06-relatedness-ladder", plot_relatedness_ladder(pairs), w = 9, h = 5.5)

  ex <- pairs[, .SD[which.min(abs(ibd - median(ibd)))], by = cls]
  sv("07-ibd-karyotype",
     plot_ibd_pair(run, ex$a, ex$b, labels = as.character(ex$cls)), w = 12, h = 4.5)

  cat("\nARG structure:\n"); print(arg_summary(run, cs$pedigree_id))
} else {
  cat("   skipped the IBD and ARG steps (no Python with tskit: see check_setup())\n")
}

say("Done. Run: ", rundir, "   Figures: ", figs)
cat("\nNext: ?plot_diagnostics, ?ibd_pairs, ?arg_summary, or vignette(\"anatomy\", \"malariafwd\")\n")
