## check_parameters.R -- run the model across the parameter space and verify that
##   (1) every configuration completes and satisfies the structural invariants, and
##   (2) each parameter moves the model in the direction it is supposed to.
## Small, fast runs -- this is a behavioural smoke test, not a calibration.
##   Rscript -e 'source(system.file("validation", "check_parameters.R", package = "malariafwd"))'
suppressMessages({
  library(malariafwd)
  library(data.table)
})
OUT <- file.path(tempdir(), "mfw_params"); unlink(OUT, recursive = TRUE)
PASS <- 0; FAIL <- 0
ok <- function(lbl, cond, detail = "") {
  if (isTRUE(cond)) { PASS <<- PASS + 1; cat(sprintf("  PASS  %-52s %s\n", lbl, detail)) }
  else              { FAIL <<- FAIL + 1; cat(sprintf("  FAIL  %-52s %s\n", lbl, detail)) }
}

## every event type that ends an active blood-stage episode (public dictionary)
.CLEARERS <- event_types()$event_type[event_types()$ends_episode]

BASE <- list(N_HOSTS = 80L, N_MOSQ = 160L, BITING_RATE = 0.10, SIM_DAYS = 400L,
             N_FOUNDERS = 10L, N_FOUNDER_HOSTS = 10L, N_FOUNDER_LINEAGES = 10L)

## run one configuration and return its summary + invariant results
run_cfg <- function(tag, ...) {
  p <- do.call(vivax_params, utils::modifyList(BASE, list(...)))
  d <- file.path(OUT, tag)
  t0 <- Sys.time()
  run <- tryCatch(run_sim(p, seed = 1, outdir = d), error = function(e) e)
  if (inherits(run, "error")) { cat(sprintf("  FAIL  %-52s %s\n", tag, conditionMessage(run)))
    FAIL <<- FAIL + 1; return(NULL) }
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  ev <- as.data.table(load_events(run)); st <- as.data.table(load_state(run))
  m  <- summarise_run(run); cs <- classify_cases(ev)

  ## structural invariants that must hold for EVERY configuration
  est <- ev[grepl("^vh_inoculation", event_type) & is_dormant == 1,
            .(batch = pedigree_id, held = n_hypnozoites)]
  act <- ev[event_type == "hypnozoite_activation", .(batch = parent1_pedigree_id, n = n_hypnozoites)]
  gone <- ev[grepl("hypnozoite_clear", event_type), .(batch = pedigree_id, n = n_hypnozoites)]
  drew <- rbind(act, gone)[, .(d = sum(n)), by = batch]
  cons <- merge(est, drew, by = "batch", all.x = TRUE)[is.na(d), d := 0]
  inv <- c(over_draw = sum(cons$d > cons$held),
           dup_ids   = nrow(ev) - uniqueN(ev$event_id),
           dup_clear = sum(duplicated(ev[event_type %in% .CLEARERS, pedigree_id])))
  data.table(tag = tag, secs = round(secs, 1), events = nrow(ev),
             prev = m$prev, coi = m$COI, relapse_frac = m$relapse_frac,
             reservoir = m$mean_res, treated = m$n_bloodstage_treated,
             imported = cs[imported == TRUE, .N], cases = nrow(cs),
             recomb_pct = round(100 * mean(ev[grepl("^hv_uptake_(clonal|recombinant)$", event_type),
                                              event_type] == "hv_uptake_recombinant"), 1),
             lineages = uniqueN(cs$lineage_id),
             over_draw = inv[["over_draw"]], dup_ids = inv[["dup_ids"]], dup_clear = inv[["dup_clear"]])
}

cat("\n================ running configurations ================\n")
CFG <- list(
  baseline     = list(),
  biting_low   = list(BITING_RATE = 0.05),
  biting_high  = list(BITING_RATE = 0.18),
  founders_2   = list(N_FOUNDER_LINEAGES = 2L),
  founders_40  = list(N_FOUNDERS = 40L, N_FOUNDER_HOSTS = 40L, N_FOUNDER_LINEAGES = 40L),
  no_relapse   = list(NU_HYPNOZOITE_BATCH = 0),
  temperate    = list(K_LATENCY = 2L, DELTA_LATENCY = 0.02),
  imports      = list(IMPORT_START = 1L, IMPORT_END = 99999L, IMPORT_RATE_PER_TICK = 0.15),
  radcure      = list(TREATMENT_START = 100L, BLOODCURE_COVERAGE = 1,
                      RADCURE_COVERAGE = 0.9),
  bloodcure    = list(TREATMENT_START = 100L, BLOODCURE_COVERAGE = 0.9),
  intervention = list(INTERVENTION_START = 100L, INTERVENTION_END = 99999L,
                      INTERVENTION_EFFICACY = 0.9),
  mda          = list(MDA_TICK = 200L),
  seasonal     = list(SEASONAL_AMP = 0.8, SEASONAL_PHASE = 180L),
  resurgence   = list(RESURGENCE_TRIGGER_TICK = 200L, BITING_RATE_BOOST = 3),
  detection    = list(P_DETECT_PRIMARY = 0.5, P_DETECT_RELAPSE = 0.15),
  wsf_multi    = list(LAMBDA_OOCYSTS = 2.0, LAMBDA_SPOROZOITES = 2.0),
  wsf_nb       = list(LAMBDA_OOCYSTS = 2.0, K_OOCYSTS = 1.0, LAMBDA_SPOROZOITES = 2.0),
  genome_core  = list(GENOME_MODE = "core"),
  chr_subset   = list(CHROMOSOMES = "1-3"),
  gametocyte   = list(GAMETOCYTE_DELAY = 3L, P_GAMETOCYTE = 0.8)
)
R <- rbindlist(lapply(names(CFG), function(k) do.call(run_cfg, c(list(tag = k), CFG[[k]]))))

cat("\n================ summary ================\n")
print(R[, .(tag, secs, prev, coi, relapse_frac, reservoir, cases, imported, treated,
            recomb_pct, lineages)])

cat("\n================ structural invariants (must hold everywhere) ================\n")
ok("all configurations completed", nrow(R) == length(CFG), sprintf("%d/%d", nrow(R), length(CFG)))
ok("no hypnozoite batch over-draws in any config", all(R$over_draw == 0),
   sprintf("max %d", max(R$over_draw)))
ok("event ids unique in every config", all(R$dup_ids == 0), sprintf("max dup %d", max(R$dup_ids)))
ok("no double blood-stage clearance", all(R$dup_clear == 0), sprintf("max %d", max(R$dup_clear)))

cat("\n================ parameter RESPONSE (does each knob do its job?) ================\n")
g <- function(k, col) R[tag == k][[col]]
ok("higher biting -> higher prevalence", g("biting_high","prev") > g("biting_low","prev"),
   sprintf("%.2f vs %.2f", g("biting_high","prev"), g("biting_low","prev")))
ok("higher biting -> higher COI", g("biting_high","coi") >= g("biting_low","coi"),
   sprintf("%.2f vs %.2f", g("biting_high","coi"), g("biting_low","coi")))
ok("more founder lineages -> more circulating lineages",
   g("founders_40","lineages") > g("founders_2","lineages"),
   sprintf("%d vs %d", g("founders_40","lineages"), g("founders_2","lineages")))
ok("NU_HYPNOZOITE_BATCH = 0 -> no relapses at all", g("no_relapse","relapse_frac") == 0,
   sprintf("relapse_frac %.3f, reservoir %.2f", g("no_relapse","relapse_frac"), g("no_relapse","reservoir")))
ok("radical cure -> smaller reservoir than baseline",
   g("radcure","reservoir") < g("baseline","reservoir"),
   sprintf("%.2f vs %.2f", g("radcure","reservoir"), g("baseline","reservoir")))
ok("blood-stage treatment -> treated cases logged & lower prevalence",
   g("bloodcure","treated") > 0 && g("bloodcure","prev") < g("baseline","prev"),
   sprintf("treated %d, prev %.2f vs %.2f", g("bloodcure","treated"),
           g("bloodcure","prev"), g("baseline","prev")))
ok("vector-control intervention -> lower prevalence",
   g("intervention","prev") < g("baseline","prev"),
   sprintf("%.2f vs %.2f", g("intervention","prev"), g("baseline","prev")))
ok("importation -> imported cases appear", g("imports","imported") > 0,
   sprintf("%d imported (baseline %d)", g("imports","imported"), g("baseline","imported")))
ok("baseline is a closed system (no imports)", g("baseline","imported") == 0)
ok("multi-oocyst uptake -> more outcrossing than single",
   g("wsf_multi","recomb_pct") >= g("baseline","recomb_pct"),
   sprintf("%.1f%% vs %.1f%%", g("wsf_multi","recomb_pct"), g("baseline","recomb_pct")))
ok("resurgence trigger -> higher prevalence than baseline",
   g("resurgence","prev") > g("baseline","prev"),
   sprintf("%.2f vs %.2f", g("resurgence","prev"), g("baseline","prev")))

## temperate latency should DELAY the first relapse
ttr <- function(k) { r <- load_run(file.path(OUT, k))
  t <- time_to_relapse(load_events(r)); if (!nrow(t)) NA_real_ else
    stats::median(t[event == 1]$time, na.rm = TRUE) }
tb <- ttr("baseline"); tt <- ttr("temperate")
ok("temperate latency delays the first relapse", is.na(tt) || is.na(tb) || tt > tb,
   sprintf("median %s d (temperate) vs %s d (tropical)", round(tt), round(tb)))

## detection layer should thin the OBSERVED line list but not the truth
dl <- function(k) { ev <- as.data.table(load_events(load_run(file.path(OUT, k))))
  cs <- classify_cases(ev); c(true = nrow(cs), det = sum(cs$detected == 1)) }
d1 <- dl("baseline"); d2 <- dl("detection")
ok("detection layer thins the observed line list only",
   d2[["det"]] / d2[["true"]] < 0.9 && d1[["det"]] == d1[["true"]],
   sprintf("baseline %d/%d detected; P_DETECT run %d/%d", d1[["det"]], d1[["true"]],
           d2[["det"]], d2[["true"]]))

## genome options must change the simulated sequence length / chromosome count
cmb <- chrom_map(load_run(file.path(OUT, "baseline")))
cmc <- chrom_map(load_run(file.path(OUT, "genome_core")))
cms <- chrom_map(load_run(file.path(OUT, "chr_subset")))
ok("GENOME_MODE = core shortens the genome", max(cmc$end) < max(cmb$end),
   sprintf("%.1f vs %.1f Mb", max(cmc$end)/1e6, max(cmb$end)/1e6))
ok("CHROMOSOMES subset simulates fewer chromosomes", nrow(cms) == 3,
   sprintf("%d chromosomes", nrow(cms)))

## ==========================================================================================
## PART 2 -- reproducibility, biology knobs, treatment efficacy, and realistic COMBINATIONS
## ==========================================================================================
cat("\n================ reproducibility ================\n")
md5 <- function(k) tools::md5sum(file.path(OUT, k, "events.tsv"))[[1]]
invisible(run_cfg("seed_a1", BITING_RATE = 0.12))
p_same <- do.call(vivax_params, utils::modifyList(BASE, list(BITING_RATE = 0.12)))
run_sim(p_same, seed = 1, outdir = file.path(OUT, "seed_a2"))   # same seed, fresh dir
run_sim(p_same, seed = 2, outdir = file.path(OUT, "seed_b"))    # different seed
ok("same seed -> byte-identical event log", md5("seed_a1") == md5("seed_a2"))
ok("different seed -> different realisation", md5("seed_a1") != md5("seed_b"))

## A directional comparison on a noisy quantity needs more than one realisation: at these
## settings the standing reservoir swings roughly 0.7-1.9 between RNG streams, so a single-seed
## "<" decides nothing whenever the true effect is comparable to that spread.
mean_reservoir <- function(tag, ...) {
  mean(vapply(1:3, function(sd) {
    d <- file.path(OUT, paste0("res_", tag, "_", sd)); unlink(d, recursive = TRUE)
    p <- do.call(vivax_params, utils::modifyList(BASE, list(...)))
    summarise_run(run_sim(p, seed = sd, outdir = d))$mean_res
  }, numeric(1)))
}

cat("\n================ relapse biology + treatment efficacy ================\n")
R2 <- rbindlist(lapply(names(cfg2 <- list(
  alpha_fast   = list(ALPHA_ACTIVATE = 1/40),        # tropical-fast activation
  alpha_slow   = list(ALPHA_ACTIVATE = 1/400),       # slow activation
  mu_high      = list(MU_HYPNOZOITE = 1/60),         # hypnozoites die fast
  hypage_short = list(MAX_HYPNOZOITE_AGE = 120L),    # hard age cap bites
  radcure_weak = list(TREATMENT_START = 100L, BLOODCURE_COVERAGE = 1,
                      RADCURE_COVERAGE = 0.9, P_RAD_SURVIVAL = 0.8),
  bloodcure_weak = list(INTERVENTION_START = 100L, INTERVENTION_END = 99999L,
                      INTERVENTION_EFFICACY = 0, BLOODCURE_COVERAGE = 0.9, P_BLOODCURE_SURVIVAL = 0.8),
  boost_relapse = list(RESURGENCE_TRIGGER_TICK = 200L, ALPHA_ACTIVATE_BOOST = 5),
  boost_batch   = list(RESURGENCE_TRIGGER_TICK = 200L, NU_HYPNOZOITE_BATCH_BOOST = 3),
  boost_import  = list(IMPORT_START = 1L, IMPORT_END = 99999L, IMPORT_RATE_PER_TICK = 0.05,
                       RESURGENCE_TRIGGER_TICK = 200L, IMPORT_RATE_BOOST = 5),
  imports_fixed = list(IMPORT_START = 1L, IMPORT_END = 99999L, N_IMPORTS_PER_TICK = 1L),
  import_pulse  = list(IMPORT_TICK = 150L, N_IMPORTS = 20L),
  no_history    = list(RETAIN_FULL_HISTORY = FALSE),
  maxmoi_cap    = list(BITING_RATE = 0.25, MAX_MOI = 2L),
  ## the realistic composite: leaky border + control programme + full case management
  china_like    = list(IMPORT_START = 1L, IMPORT_END = 99999L, IMPORT_RATE_PER_TICK = 0.08,
                       INTERVENTION_START = 120L, INTERVENTION_END = 300L,
                       INTERVENTION_EFFICACY = 0.7, RADCURE_COVERAGE = 0.6,
                       BLOODCURE_COVERAGE = 0.8, P_DETECT_RELAPSE = 0.2,
                       RESURGENCE_TRIGGER_TICK = 300L, BITING_RATE_BOOST = 2)
)), function(k) do.call(run_cfg, c(list(tag = k), cfg2[[k]]))))
print(R2[, .(tag, secs, prev, coi, relapse_frac, reservoir, cases, imported, treated, lineages)])

ok("part-2 configurations all completed", nrow(R2) == length(cfg2),
   sprintf("%d/%d", nrow(R2), length(cfg2)))
ok("invariants hold in part 2", all(R2$over_draw == 0) && all(R2$dup_ids == 0) && all(R2$dup_clear == 0))
h <- function(k, col) R2[tag == k][[col]]
ok("faster activation -> higher relapse fraction",
   h("alpha_fast","relapse_frac") > h("alpha_slow","relapse_frac"),
   sprintf("%.2f vs %.2f", h("alpha_fast","relapse_frac"), h("alpha_slow","relapse_frac")))
ok("faster activation -> smaller standing reservoir",
   h("alpha_fast","reservoir") < h("alpha_slow","reservoir"),
   sprintf("%.2f vs %.2f", h("alpha_fast","reservoir"), h("alpha_slow","reservoir")))
mu_res <- mean_reservoir("muhigh", MU_HYPNOZOITE = 1/60)
bl_res <- mean_reservoir("base")
ok("high hypnozoite death -> smaller reservoir than baseline",
   mu_res < bl_res, sprintf("%.2f vs %.2f (3-seed means)", mu_res, bl_res))
ok("short MAX_HYPNOZOITE_AGE -> smaller reservoir",
   h("hypage_short","reservoir") < g("baseline","reservoir"),
   sprintf("%.2f vs %.2f", h("hypage_short","reservoir"), g("baseline","reservoir")))
## 3-seed, for the same reason mean_reservoir() exists: single-seed this comparison was a
## near coin flip (measured 5 of 6 seeds in the right direction on the pre-fence model), so a
## shifted RNG stream flipped it with no change in the underlying behaviour. Over 3 seeds the
## effect is unambiguous - weak 0.50 against strong 0.12, 6 of 6 seeds - and the check now
## measures radical-cure efficacy rather than which stream it drew.
rc_weak   <- mean_reservoir("rcw", TREATMENT_START = 100L, BLOODCURE_COVERAGE = 1,
                            RADCURE_COVERAGE = 0.9, P_RAD_SURVIVAL = 0.8)
rc_strong <- mean_reservoir("rcs", TREATMENT_START = 100L, BLOODCURE_COVERAGE = 1,
                            RADCURE_COVERAGE = 0.9)
ok("weak radical cure leaves a bigger reservoir than strong",
   rc_weak > rc_strong,
   sprintf("%.2f (P_RAD_SURVIVAL .8) vs %.2f (.05), 3-seed means", rc_weak, rc_strong))
ok("weak blood cure leaves higher prevalence than strong",
   h("bloodcure_weak","prev") > g("bloodcure","prev"),
   sprintf("%.2f vs %.2f", h("bloodcure_weak","prev"), g("bloodcure","prev")))
ok("relapse-activation boost raises the relapse fraction",
   h("boost_relapse","relapse_frac") >= g("baseline","relapse_frac"),
   sprintf("%.2f vs %.2f", h("boost_relapse","relapse_frac"), g("baseline","relapse_frac")))
ok("batch boost enlarges the reservoir", h("boost_batch","reservoir") > g("baseline","reservoir"),
   sprintf("%.2f vs %.2f", h("boost_batch","reservoir"), g("baseline","reservoir")))
ok("all three importation modes produce imports",
   h("boost_import","imported") > 0 && h("imports_fixed","imported") > 0 && h("import_pulse","imported") > 0,
   sprintf("rate+boost %d, fixed %d, pulse %d", h("boost_import","imported"),
           h("imports_fixed","imported"), h("import_pulse","imported")))
ok("RETAIN_FULL_HISTORY = FALSE still runs", h("no_history","cases") > 0,
   sprintf("%d cases (IBD-to-prior not guaranteed in this mode)", h("no_history","cases")))
ok("MAX_MOI caps COI", h("maxmoi_cap","coi") <= 2.0, sprintf("COI %.2f with MAX_MOI = 2", h("maxmoi_cap","coi")))
ok("composite China-like scenario runs and exercises every lever",
   h("china_like","imported") > 0 && h("china_like","treated") > 0 && h("china_like","cases") > 0,
   sprintf("%d cases, %d imported, %d treated", h("china_like","cases"),
           h("china_like","imported"), h("china_like","treated")))

## ==========================================================================================
## PART 3 -- run_grid(): the parallel sweep entry point most users will actually reach for.
## Worth its own check because a sweep failing would be a bad first experience, and because
## future/furrr workers load the INSTALLED package -- a failure here would not show up in any
## test run under devtools::load_all().
## ==========================================================================================
cat("\n================ run_grid (parallel sweep) ================\n")
gd <- file.path(OUT, "grid")
sets <- list(low  = do.call(vivax_params, utils::modifyList(BASE, list(BITING_RATE = 0.06))),
             high = do.call(vivax_params, utils::modifyList(BASE, list(BITING_RATE = 0.16))))
gr <- tryCatch(run_grid(sets, seeds = 1:2, outdir_base = gd, workers = 2),
               error = function(e) e)
if (inherits(gr, "error")) {
  ok("run_grid() completes", FALSE, conditionMessage(gr))
} else {
  gr <- data.table::as.data.table(gr)
  ok("run_grid() completes", TRUE, sprintf("%d runs", nrow(gr)))
  ok("run_grid returns one row per (set x seed)", nrow(gr) == length(sets) * 2,
     sprintf("%d rows", nrow(gr)))
  ok("every run_grid run reports ok", all(gr$ok %in% c(TRUE, 1)), "")
  outs <- list.dirs(gd, recursive = FALSE)
  ok("every run_grid run wrote its outputs",
     all(file.exists(file.path(outs, "events.tsv"))) &&
     all(file.exists(file.path(outs, "plasmo.trees"))), sprintf("%d dirs", length(outs)))
  ## the sweep must actually vary: different parameter sets give different prevalence,
  ## and different seeds give different realisations of the same set
  mm <- rbindlist(lapply(outs, function(o) {
    m <- summarise_run(o); data.table(dir = basename(o), prev = m$prev) }))
  lo <- mm[grepl("low", dir), mean(prev)]; hi <- mm[grepl("high", dir), mean(prev)]
  ok("sweep separates the parameter sets", hi > lo, sprintf("prev high %.2f vs low %.2f", hi, lo))
  ok("seeds within a set give different realisations",
     uniqueN(mm[grepl("high", dir), prev]) > 1 || uniqueN(mm[grepl("low", dir), prev]) > 1)
}

cat(sprintf("\n================ TOTAL %d passed, %d failed ================\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)
