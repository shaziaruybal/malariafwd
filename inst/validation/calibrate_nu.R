## calibrate_nu.R -- fit NU_HYPNOZOITE_BATCH, the mean hypnozoite batch a bite seeds, against four
## calibration targets, ranges taken from published P. vivax studies, at moderate transmission.
## Section 13 of the model description gives the targets, their sources and the table this script
## makes.
##   Rscript -e 'source(system.file("validation", "calibrate_nu.R", package = "malariafwd"))'
## Each value of nu runs 400 hosts and 800 mosquitoes for 5 years at BITING_RATE 0.08, with 20
## founder hosts carrying 10 lineages, for seeds 1 to 3: 15 runs, about 6 minutes on one core, or
## less with MALARIAFWD_WORKERS set to the cores to use. The other parameters are the defaults, so a
## bite can co-transmit several genotypes. The first 2 years are burn-in. Each figure is the mean
## over the seeds of summarise_run()'s, counting every relapse and, where marked, with masked
## relapses (into a patent clone of their own lineage) left out. A run can die out below the fitted
## value: it counts with no infection after its extinction, and the MOI and the relapse fraction,
## which need cases, are averaged over the runs that have them.
suppressMessages({
  library(malariafwd)
  library(data.table)
})
## MALARIAFWD_CALIBRATE_OUT puts the tables somewhere they outlive the run, for
## dev/calibration-table.R to generate section 13's table from.
OUT <- Sys.getenv("MALARIAFWD_CALIBRATE_OUT", unset = file.path(tempdir(), "mfw_calibrate_nu"))
unlink(OUT, recursive = TRUE)
PASS <- 0; FAIL <- 0
ok <- function(lbl, cond, detail = "") {
  if (isTRUE(cond)) { PASS <<- PASS + 1; cat(sprintf("  PASS  %-52s %s\n", lbl, detail)) }
  else              { FAIL <<- FAIL + 1; cat(sprintf("  FAIL  %-52s %s\n", lbl, detail)) }
}

NU <- c(2, 3, 5, 8, 12)
SEEDS <- 1:3
BURN_FRAC <- 730 / 1825
## the calibration targets, as section 13 of the model description gives them with their sources
TARGETS <- data.table(metric = c("relapses_per_py", "molfob_per_py", "moi_median", "relapse_frac"),
                      label = c("relapses per person-year", "molFOB", "median MOI", "relapse fraction"),
                      lo = c(0.4, 1, 1, 0.45), hi = c(4.5, 15, 2, 0.9))

sets <- lapply(NU, function(nu)
  vivax_params(N_HOSTS = 400L, N_MOSQ = 800L, SIM_DAYS = 1825L, BITING_RATE = 0.08,
               N_FOUNDER_HOSTS = 20L, N_FOUNDER_LINEAGES = 10L, NU_HYPNOZOITE_BATCH = nu))
names(sets) <- paste0("nu", NU)
workers <- suppressWarnings(as.integer(Sys.getenv("MALARIAFWD_WORKERS", "1")))
cat(sprintf("\n================ NU_HYPNOZOITE_BATCH against the targets (%d runs) ================\n",
            length(NU) * length(SEEDS)))
t0 <- Sys.time()
g <- run_grid(sets, seeds = SEEDS, outdir_base = OUT,
              workers = if (!is.na(workers) && workers > 1L) workers)
cat(sprintf("  %d runs in %.1f min\n", nrow(g), as.numeric(Sys.time() - t0, units = "mins")))
shipped <- vivax_params()$NU_HYPNOZOITE_BATCH
at_shipped <- g$set == paste0("nu", shipped)
ok("every run completes", all(g$ok))
ok(sprintf("no run at the shipped value, %g, dies out", shipped),
   !any(g$never_started[at_shipped]) && all(is.na(g$extinct_tick[at_shipped])))

per_run <- rbindlist(lapply(seq_len(nrow(g)), function(i) {
  all <- summarise_run(g$outdir[i], burn_frac = BURN_FRAC)
  unm <- summarise_run(g$outdir[i], burn_frac = BURN_FRAC, exclude_masked = TRUE)
  ## the genotypes each bite after the burn-in delivered
  bites <- load_events(g$outdir[i])[event_type == "vh_inoculation" & tick >= all$from_tick,
                                    .(k = uniqueN(lineage_id)), by = inoculation_id]
  data.table(nu = NU[match(g$set[i], names(sets))], seed = g$seed[i],
             extinct = !is.na(g$extinct_tick[i]),
             prevalence = all$prevalence_frac,
             relapses_per_py = all$relapses_per_py_at_risk,
             relapses_per_py_unmasked = unm$relapses_per_py_at_risk,
             molfob_per_py = all$molfob_clones_per_py,
             molfob_per_py_unmasked = unm$molfob_clones_per_py,
             moi_median = all$coi_median,
             relapse_frac = all$relapse_frac_acquisitions,
             relapse_frac_unmasked = unm$relapse_frac_acquisitions,
             masked_share = all$n_relapse_masked / all$n_relapse,
             eir_per_py = all$eir_per_py,
             multi_bites = if (nrow(bites)) mean(bites$k > 1) else NA_real_)
}))
fwrite(per_run, file.path(OUT, "calibrate_nu_runs.tsv"), sep = "\t")

m <- per_run[, c(list(extinct = sum(extinct)), lapply(.SD, mean, na.rm = TRUE)), by = nu,
             .SDcols = !c("seed", "extinct")]
inside <- function(x, target) x >= TARGETS[metric == target, lo] & x <= TARGETS[metric == target, hi]
m[, in_all := inside(relapses_per_py, "relapses_per_py") & inside(molfob_per_py, "molfob_per_py") &
              inside(moi_median, "moi_median") & inside(relapse_frac, "relapse_frac")]
m[, in_all_unmasked := inside(relapses_per_py_unmasked, "relapses_per_py") &
                       inside(molfob_per_py_unmasked, "molfob_per_py") & inside(moi_median, "moi_median") &
                       inside(relapse_frac_unmasked, "relapse_frac")]
fwrite(m, file.path(OUT, "calibrate_nu.tsv"), sep = "\t")
fwrite(TARGETS, file.path(OUT, "calibrate_nu_targets.tsv"), sep = "\t")

cat("\n  means over seeds", paste(SEEDS, collapse = ", "), "(targets:",
    paste(sprintf("%s %g-%g", TARGETS$label, TARGETS$lo, TARGETS$hi), collapse = ", "), ")\n\n")
show <- m[, .(nu, prevalence = round(prevalence, 2),
              `relapses/py` = round(relapses_per_py, 2), `relapses/py unmasked` = round(relapses_per_py_unmasked, 2),
              molFOB = round(molfob_per_py, 2), `molFOB unmasked` = round(molfob_per_py_unmasked, 2),
              `median MOI` = round(moi_median, 2),
              `relapse frac` = round(relapse_frac, 2), `relapse frac unmasked` = round(relapse_frac_unmasked, 2),
              `masked share` = round(masked_share, 2), EIR = round(eir_per_py, 2),
              `bites >1 genotype` = round(multi_bites, 3), `died out` = extinct,
              `in all 4` = in_all, `in all 4 unmasked` = in_all_unmasked)]
print(show, row.names = FALSE)
cat("\n  the four metrics over seeds (min - max):\n")
span <- function(x, fmt) sprintf(paste(fmt, "-", fmt), min(x, na.rm = TRUE), max(x, na.rm = TRUE))
rng <- per_run[, .(relapses_per_py = span(relapses_per_py, "%.2f"), molfob_per_py = span(molfob_per_py, "%.2f"),
                   moi_median = span(moi_median, "%g"), relapse_frac = span(relapse_frac, "%.3f")), by = nu]
print(rng, row.names = FALSE)

cat("\n")
ok(sprintf("the shipped value, %g, sits inside all four targets", shipped),
   isTRUE(m[nu == shipped, in_all]),
   sprintf("(counting every relapse; inside all four: nu = %s)", paste(m[in_all == TRUE, nu], collapse = ", ")))
cat(sprintf("\n  Tables: %s\n\n================ %d passed, %d failed ================\n", OUT, PASS, FAIL))
if (FAIL > 0) quit(status = 1)
