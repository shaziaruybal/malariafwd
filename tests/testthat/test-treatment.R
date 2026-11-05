## Case management. These need SLiM (skip_if_no_slim() lives in helper-slim.R).

shared <- function(...) vivax_params(SIM_DAYS = 120L, N_HOSTS = 40L, N_MOSQ = 80L,
                                     BITING_RATE = 0.12, ...)

run_counts <- function(p, tag) {
  out <- file.path(tempdir(), paste0("mfw-tx-", tag))
  unlink(out, recursive = TRUE)
  r <- run_sim(p, seed = 2, outdir = out)
  ev <- data.table::as.data.table(load_events(r))
  list(bloodcure = ev[event_type == "case_bloodstage_clear", .N],
       radcure   = ev[event_type == "radcure_hypnozoite_clear", .N])
}

test_that("treatment is driven by coverage alone, not by the vector-control window", {
  skip_if_no_slim()
  ## D19. Every case-based treatment used to be gated on INTERVENTION_START/END - the same
  ## window that scales biting and importation - so setting a coverage without also opening
  ## that window produced zero treatments and no warning. That is how the issue was found.
  expect_equal(run_counts(shared(), "none")$bloodcure, 0)

  blood <- run_counts(shared(BLOODCURE_COVERAGE = 1), "blood")
  expect_gt(blood$bloodcure, 0)

  ## radical cure is NESTED inside blood-stage treatment, so it needs both coverages
  rad <- run_counts(shared(BLOODCURE_COVERAGE = 1, RADCURE_COVERAGE = 1), "rad")
  expect_gt(rad$radcure, 0)
})

test_that("the vector-control window does not switch treatment on", {
  skip_if_no_slim()
  ## ...and conversely: opening the intervention window must change transmission without
  ## producing a single treatment, or the two are still entangled.
  vc <- run_counts(shared(INTERVENTION_START = 1L, INTERVENTION_END = 120L,
                          INTERVENTION_EFFICACY = 0.9), "vc")
  expect_equal(vc$bloodcure, 0)
  expect_equal(vc$radcure, 0)
})

test_that("the case-management window bounds WHEN treatment happens", {
  skip_if_no_slim()
  ## Deliberately an assertion about TIMING, not about counts. Since treatment became
  ## host-level it suppresses transmission hard enough that an always-open window can
  ## produce FEWER courses than a late one - it prevents the cases it would have treated.
  out <- file.path(tempdir(), "mfw-tx-window")
  unlink(out, recursive = TRUE)
  r <- run_sim(shared(BLOODCURE_COVERAGE = 1, TREATMENT_START = 90L, TREATMENT_END = 120L),
               seed = 2, outdir = out)
  ev <- data.table::as.data.table(load_events(r))
  given <- ev[event_type == "bloodcure_given", tick]
  expect_gt(length(given), 0)
  ## Nothing before the window opens: no case can be DECIDED earlier.
  expect_true(all(given >= 90))
  ## Courses decided just before it shuts are delivered TREATMENT_DELAY days later, so a
  ## short tail past TREATMENT_END is correct - the window gates prescription, not delivery.
  expect_true(all(given < 120 + 5 * 7))
})

test_that("treatment and recrudescence delays are geometric with mean D days, at least one day", {
  skip_if_no_slim()
  ## One pulse of detected imports on day 5 is the only source of cases: a single mosquito
  ## that almost never bites, relapses and recrudescences go undetected, and founders are
  ## never treated. No course clears anything: every clone still patent when its course lands
  ## is suppressed, and each recrudescence is paired with its suppression through
  ## parent1_pedigree_id.
  delays <- function(d) {
    run <- run_sim(vivax_params(N_HOSTS = 1200, N_MOSQ = 1, BITING_RATE = 1e-6, SIM_DAYS = 45,
                                IMPORT_TICK = 5, N_IMPORTS = 1000,
                                P_DETECT_RELAPSE = 0, P_DETECT_RECRUD = 0,
                                BLOODCURE_COVERAGE = 1, P_BLOODCURE_SURVIVAL = 1,
                                P_RECRUD_GIVEN_FAIL = 1,
                                TREATMENT_DELAY = d, MEAN_RECRUD_DELAY = d),
                   seed = 3, outdir = file.path(tempdir(), paste0("mfw-delay-", d)))
    ev <- load_events(run)
    gap <- merge(ev[event_type == "recrudescence", .(holder = parent1_pedigree_id, t_rec = tick)],
                 ev[event_type == "bloodcure_subpatent", .(holder = pedigree_id, t_sub = tick)],
                 by = "holder")
    list(course = ev[event_type == "bloodcure_given", tick] - 5L, recrud = gap$t_rec - gap$t_sub)
  }
  one <- delays(1)
  expect_gt(length(one$course), 900)
  expect_gt(length(one$recrud), 800)
  ## D = 1: no delay comes later than the next day. A zero-day draw would also show as one day
  ## here (imports are seeded after the day's courses, and suppressions happen after the day's
  ## recrudescences), so the one-day share at D = 4 is what pins the floor.
  expect_true(all(one$course == 1))
  expect_true(all(one$recrud == 1))
  four <- delays(4)
  for (nm in names(four)) {
    x <- four[[nm]]
    ## Both at about 4 standard errors: the delay has mean 4 and variance 12, the one-day share
    ## is 1/4 with variance 3/16. max(1, rgeom(1/(1 + D))) gives a one-day share of 0.36, and a
    ## draw that can return 0 gives 0.44.
    expect_lt(abs(mean(x) - 4), 4 * sqrt(12 / length(x)), label = paste("mean", nm, "delay - 4"))
    expect_lt(abs(mean(x == 1) - 1/4), 4 * sqrt(3 / 16 / length(x)),
              label = paste("one-day share of", nm, "delays - 1/4"))
  }
})

test_that("a coverage against a closed treatment window is refused", {
  ## The failure mode D19 fixed was silence, so re-creating it in the NEW parameter has to
  ## be an error rather than a quietly empty result.
  expect_error(vivax_params(BLOODCURE_COVERAGE = 1, TREATMENT_START = 50L, TREATMENT_END = 50L),
               "case-management window is closed")
  ## suppressWarnings: this also trips the nested-radical-cure warning, which is correct
  ## but not what this expectation is about
  suppressWarnings(
    expect_error(vivax_params(RADCURE_COVERAGE = 1, SIM_DAYS = 100L, TREATMENT_START = 500L),
                 "case-management window is closed"))
  ## the defaults are "always on", so a coverage on its own is fine
  expect_no_error(vivax_params(BLOODCURE_COVERAGE = 1))
})

test_that("three-state recurrence: suppression, recrudescence, and a clonal mint", {
  skip_if_no_slim()
  base <- list(SIM_DAYS = 300L, N_HOSTS = 50L, N_MOSQ = 110L, BITING_RATE = 0.14,
               P_DETECT_PRIMARY = 0.9, P_DETECT_RELAPSE = 0.6)
  out <- file.path(tempdir(), "mfw-recrud"); unlink(out, recursive = TRUE)

  ## OFF is a strict no-op: no suppression, no recrudescence, no subpatent compartment
  off <- run_sim(do.call(vivax_params, c(base, list(BLOODCURE_COVERAGE = 1,
                 P_BLOODCURE_SURVIVAL = 0.5, P_RECRUD_GIVEN_FAIL = 0))),
                 seed = 11L, outdir = file.path(out, "off"))
  eo <- data.table::as.data.table(load_events(off))
  so <- data.table::as.data.table(load_state(off))
  expect_equal(nrow(eo[event_type == "recrudescence"]), 0L)
  expect_equal(nrow(eo[event_type == "bloodcure_subpatent"]), 0L)
  expect_equal(max(so$moi_subpatent), 0L)
  expect_gt(nrow(eo[event_type == "bloodcure_persist"]), 0)   # failures still happen

  ## ON: every failure suppresses, and suppressed clones re-expand
  on <- run_sim(do.call(vivax_params, c(base, list(BLOODCURE_COVERAGE = 1, RADCURE_COVERAGE = 1,
                P_BLOODCURE_SURVIVAL = 0.5, P_RECRUD_GIVEN_FAIL = 1))),
                seed = 11L, outdir = file.path(out, "on"))
  ev <- data.table::as.data.table(load_events(on))
  st <- data.table::as.data.table(load_state(on))
  sub <- ev[event_type == "bloodcure_subpatent"]
  rec <- ev[event_type == "recrudescence"]
  expect_gt(nrow(sub), 0); expect_gt(nrow(rec), 0)
  expect_equal(nrow(ev[event_type == "bloodcure_persist"]), 0L)  # all failures suppress
  expect_gt(max(st$moi_subpatent), 0L)

  ## D25: a recrudescence MINTS a new genome whose parent is the holder
  expect_true(all(rec$parent1_pedigree_id %in% sub$pedigree_id))
  expect_true(all(rec$pedigree_id != rec$parent1_pedigree_id))   # a NEW id, not the holder's
  ## a recrudescent clone is patent again, so a later course can fail on it and suppress it
  ## a second time - chains are expected. Any such suppression must postdate its own mint.
  cyc <- merge(rec[, .(pedigree_id, t_rec = tick)],
               sub[, .(pedigree_id, t_sub = tick)], by = "pedigree_id")
  if (nrow(cyc)) expect_true(all(cyc$t_sub > cyc$t_rec))
  expect_equal(nrow(rec), data.table::uniqueN(rec$pedigree_id))  # one episode, one id
  ## -1 for both, as a relapse does: not a V-H event, and inheriting them would place the
  ## recurrence inside another case's co-transmission group
  expect_true(all(rec$oocyst_id == -1)); expect_true(all(rec$inoculation_id == -1))
  ## it re-expands strictly after it was suppressed
  gap <- merge(rec[, .(holder = parent1_pedigree_id, t_rec = tick)],
               sub[, .(holder = pedigree_id, t_sub = tick)], by = "holder")
  expect_true(all(gap$t_rec > gap$t_sub))

  ## it is a CASE, with its own origin
  cs <- data.table::as.data.table(classify_cases(ev))
  expect_true("recrudescence" %in% cs$origin)
  expect_equal(nrow(cs[origin == "recrudescence"]), nrow(rec))

  ## ...but NOT an acquisition: it must not enter a force-of-infection rate
  fo <- data.table::as.data.table(force_of_infection(on))
  expect_true("is_acquisition" %in% names(fo))
  expect_false(any(fo[origin == "recrudescence"]$is_acquisition))
  rf <- relapse_fraction(ev)
  expect_gt(rf$relapse_fraction, rf$relapse_frac_episodes)   # bigger denominator, smaller value

  ## episode conservation: both sides of the vocabulary are derived, so nothing is orphaned
  onsets <- ev[event_type %in% .onset_events() & is_dormant == 0, pedigree_id]
  clears <- ev[event_type %in% .bloodstage_clear_events(), pedigree_id]
  expect_true(all(clears %in% onsets))
  expect_equal(data.table::uniqueN(clears), length(clears))   # nothing clears twice
})

test_that("a recrudescence is clonal to its holder and MOI-neutral", {
  skip_if_no_slim()
  skip_if_not(requireNamespace("reticulate", quietly = TRUE), "reticulate not available")
  out <- file.path(tempdir(), "mfw-recrud-ibd"); unlink(out, recursive = TRUE)
  run <- run_sim(vivax_params(SIM_DAYS = 300L, N_HOSTS = 50L, N_MOSQ = 110L,
                              BITING_RATE = 0.14, BLOODCURE_COVERAGE = 1,
                              P_BLOODCURE_SURVIVAL = 0.5, P_RECRUD_GIVEN_FAIL = 1),
                 seed = 11L, outdir = out)
  ev <- data.table::as.data.table(load_events(run))
  rec <- ev[event_type == "recrudescence"][1:min(.N, 12)]
  skip_if(nrow(rec) < 2, "too few recrudescences in this run")

  ## the payoff of minting: "IBD 1 to the prior episode" is a real measurement, not a
  ## self-pair that ibd_index() hard-codes to 1
  idx <- ibd_index(ibd_pairs(run, ids = unique(c(rec$pedigree_id, rec$parent1_pedigree_id))))
  expect_true(all(idx$getibd(rec$pedigree_id, rec$parent1_pedigree_id) > 0.999))
  tm <- data.table::as.data.table(arg_tmrca(run, rec$pedigree_id, rec$parent1_pedigree_id))
  expect_true(all(tm$coalesced > 0.99, na.rm = TRUE))
  expect_true(all(abs(tm$tmrca) < 1e-9, na.rm = TRUE))        # depth 0: parent -> child

  ## the mint inherits lineage_id, so retire-holder + mint is exactly MOI-neutral
  hold_lin <- ev[event_type == "bloodcure_subpatent", .(pedigree_id, lineage_id)]
  m <- merge(rec[, .(pedigree_id, lineage_id, holder = parent1_pedigree_id)],
             hold_lin, by.x = "holder", by.y = "pedigree_id", suffixes = c("", "_holder"))
  expect_equal(as.character(m$lineage_id), as.character(m$lineage_id_holder))
})
