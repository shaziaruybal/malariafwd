## A hand-made event log with every kind of transmission link, and a run directory holding it.

## An events table from `rows` (lists with at least tick and event_type), in load_events()'s columns
## and classes; every other column takes the value logEvent() writes when it is not given.
make_events <- function(rows = list()) {
  def <- list(inoculation_id = -1L, mosq_id = -1L, human_id = -1L, lineage_id = NA_integer_,
              parent1_lineage_id = NA_integer_, parent2_lineage_id = NA_integer_, is_recombinant = 0L,
              p_outcross = 0, recipient = "human", is_dormant = 0L, is_imported = 0L, pedigree_id = -1L,
              parent1_pedigree_id = -1L, parent2_pedigree_id = -1L, n_hypnozoites = 1L, oocyst_id = -1L,
              detected = -1L, into_patent_lineage = -1L)
  d <- if (length(rows)) data.table::rbindlist(rows, fill = TRUE)
       else data.table::data.table(tick = integer(), event_type = character())
  for (nm in names(def)) {
    if (!nm %in% names(d)) d[, (nm) := rep(def[[nm]], .N)]
    else if (!is.na(def[[nm]])) d[is.na(get(nm)), (nm) := def[[nm]]]
  }
  d <- d[order(tick)]
  d[, event_id := paste0(tick, "_", seq_len(.N)), by = tick]
  sch <- malariafwd:::.EVENTS_SCHEMA
  for (nm in names(sch))
    data.table::set(d, j = nm, value = switch(sch[[nm]], integer = as.integer, double = as.double,
                                              character = as.character)(d[[nm]]))
  data.table::setcolorder(d, names(sch))
  d[]
}

## Hosts 0-8, mosquitoes 100-106; numbers are pedigree ids, lineages in brackets.
##   0  founder in host 0 [0] (pedigree 0 is a real genome)
##   2  host 1, bite day 16 by mosquito 100 that fed on host 0 on day 5 (selfed oocyst 1) [0];
##   3  the same bite's batch; relapses as 4 (day 60, masked: clone 2 is patent) and 5 (day 120)
##   7  import into host 3, day 30 [7]; its batch 8 (day 31) relapses as 9 on day 90 (oocyst -1)
##  10  batch-only bite on host 4, day 50, by mosquito 102 that fed on host 3 on day 35 (oocyst 3) [7];
##      relapses as 11 on day 140
##  14  host 1 again, day 38, by mosquito 104 that fed on host 1 on day 25 (oocyst 5) [0]: a self link
##  15, 16  one bite on host 6, day 48, by mosquito 105 delivering two genotypes, from host 1
##      (oocyst 6, fed day 26) [0] and from host 3 (oocyst 7, fed day 36) [7]
##   6  host 2, day 85, by mosquito 101 that fed on host 6 on day 70: gametes of 15 and 16 outcrossed
##      (oocyst 2) [950]; 17 the same bite's batch, relapses as 18 on day 130 (masked: 12 is patent)
##  12  recrudescence of 6 in host 2 on day 105 (bloodcure_subpatent day 95) [950]
##  13  host 5, day 122, by mosquito 103 that fed on host 2 on day 110 (parent 12, oocyst 4) [950]
##  30  founder in host 7 [30]; 31 host 8, day 22, by mosquito 106 that fed on host 7 on day 10 (oocyst 8)
transmission_fixture <- function() {
  up <- function(tick, mosq, host, p1, ooc, p2 = -1L, l1, l2 = NA_integer_)
    list(tick = tick, event_type = if (p2 == -1L) "hv_uptake_clonal" else "hv_uptake_recombinant",
         mosq_id = mosq, human_id = host, recipient = "mosquito",
         lineage_id = if (p2 == -1L) l1 else NA_integer_, parent1_lineage_id = l1, parent2_lineage_id = l2,
         parent1_pedigree_id = p1, parent2_pedigree_id = p2, is_recombinant = as.integer(p2 != -1L),
         oocyst_id = ooc)
  inoc <- function(tick, mosq, host, ped, ooc, lin, inoc_id, dormant = 0L, n = 1L)
    list(tick = tick, event_type = "vh_inoculation", inoculation_id = inoc_id, mosq_id = mosq,
         human_id = host, lineage_id = lin, pedigree_id = ped, parent1_pedigree_id = 900L + ped,
         oocyst_id = ooc, is_dormant = dormant, n_hypnozoites = n, detected = if (dormant) -1L else 1L)
  act <- function(tick, host, ped, batch, ooc, lin, masked = 0L, imported = 0L)
    list(tick = tick, event_type = "hypnozoite_activation", human_id = host, lineage_id = lin,
         pedigree_id = ped, parent1_pedigree_id = batch, oocyst_id = ooc, n_hypnozoites = 1L,
         into_patent_lineage = masked, is_imported = imported, detected = 1L)
  make_events(list(
    list(tick = 1L, event_type = "vh_inoculation_seed", human_id = 0L, lineage_id = 0L, pedigree_id = 0L),
    list(tick = 1L, event_type = "vh_inoculation_seed", human_id = 7L, lineage_id = 30L, pedigree_id = 30L),
    up(5L, 100L, 0L, 0L, 1L, l1 = 0L),
    up(10L, 106L, 7L, 30L, 8L, l1 = 30L),
    inoc(16L, 100L, 1L, 2L, 1L, 0L, 1L), inoc(16L, 100L, 1L, 3L, 1L, 0L, 1L, dormant = 1L, n = 2L),
    inoc(22L, 106L, 8L, 31L, 8L, 30L, 2L),
    up(25L, 104L, 1L, 2L, 5L, l1 = 0L),
    up(26L, 105L, 1L, 2L, 6L, l1 = 0L),
    list(tick = 30L, event_type = "vh_inoculation_import", inoculation_id = 3L, human_id = 3L, lineage_id = 7L,
         pedigree_id = 7L, is_imported = 1L, detected = 1L),
    list(tick = 31L, event_type = "vh_inoculation_import", inoculation_id = 3L, human_id = 3L, lineage_id = 7L,
         pedigree_id = 8L, parent1_pedigree_id = 7L, is_imported = 1L, is_dormant = 1L, n_hypnozoites = 1L),
    up(35L, 102L, 3L, 7L, 3L, l1 = 7L),
    up(36L, 105L, 3L, 7L, 7L, l1 = 7L),
    inoc(38L, 104L, 1L, 14L, 5L, 0L, 4L),
    inoc(48L, 105L, 6L, 15L, 6L, 0L, 5L), inoc(48L, 105L, 6L, 16L, 7L, 7L, 5L),
    inoc(50L, 102L, 4L, 10L, 3L, 7L, 6L, dormant = 1L),
    act(60L, 1L, 4L, 3L, 1L, 0L, masked = 1L),
    up(70L, 101L, 6L, 15L, 2L, p2 = 16L, l1 = 0L, l2 = 7L),
    inoc(85L, 101L, 2L, 6L, 2L, 950L, 7L), inoc(85L, 101L, 2L, 17L, 2L, 950L, 7L, dormant = 1L),
    act(90L, 3L, 9L, 8L, -1L, 7L, imported = 1L),
    list(tick = 95L, event_type = "bloodcure_subpatent", human_id = 2L, lineage_id = 950L, pedigree_id = 6L),
    list(tick = 105L, event_type = "recrudescence", human_id = 2L, lineage_id = 950L, pedigree_id = 12L,
         parent1_pedigree_id = 6L, detected = 1L),
    up(110L, 103L, 2L, 12L, 4L, l1 = 950L),
    act(120L, 1L, 5L, 3L, 1L, 0L),
    inoc(122L, 103L, 5L, 13L, 4L, 950L, 8L),
    act(130L, 2L, 18L, 17L, 2L, 950L, masked = 1L),
    act(140L, 4L, 11L, 10L, 3L, 7L)
  ))
}

## The fixture as a run of 9 hosts and 150 days.
transmission_fixture_run <- function(events = transmission_fixture())
  load_run(fake_run(params = list(N_HOSTS = 9L, N_MOSQ = 7L, SIM_DAYS = 150L, N_FOUNDER_HOSTS = 2L),
                    events = events))

## A run with co-transmission, imports, recrudescence and batch-only bites, made once a session.
.links_runs <- new.env()
links_run <- function() {
  if (is.null(.links_runs$run))
    .links_runs$run <- suppressMessages(run_sim(vivax_params(
      N_HOSTS = 60, N_MOSQ = 120, BITING_RATE = 0.16, SIM_DAYS = 300, LAMBDA_SPOROZOITES = 2, IMPORT_START = 1L,
      IMPORT_END = 300L, IMPORT_RATE_PER_TICK = 0.1, BLOODCURE_COVERAGE = 0.5, P_BLOODCURE_SURVIVAL = 0.3,
      P_RECRUD_GIVEN_FAIL = 0.7, P_PRIMARY = 0.8, CHROMOSOMES = "1"), seed = 2, outdir = tempfile("mfwd_links_")))
  .links_runs$run
}
