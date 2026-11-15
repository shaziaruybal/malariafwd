## Runs shared by the window and episode tests, made once per session.

## W goes extinct early for seeds 1 to 6 (days 43 to 196), before its burn-in of 200 days ends.
w_params <- function() vivax_params(N_HOSTS = 60, N_MOSQ = 240, BITING_RATE = 0.12, SIM_DAYS = 400,
                                    N_FOUNDER_HOSTS = 12, NU_HYPNOZOITE_BATCH = 0,
                                    MEAN_BLOOD_STAGE_DURATION = 10, CHROMOSOMES = "1",
                                    LAMBDA_OOCYSTS = 0, LAMBDA_SPOROZOITES = 0)
.w_runs <- new.env()
w_run <- function(seed) {
  key <- as.character(seed)
  if (is.null(.w_runs[[key]]))
    .w_runs[[key]] <- suppressMessages(run_sim(w_params(), seed = seed, outdir = tempfile("mfwd_w_")))
  .w_runs[[key]]
}

## A hand-built run for value-level tests of episodes and relapses (no SLiM). Every row is one the
## model writes; columns not given take the model's defaults for that row kind. Host 1's relapse
## on day 20 (pedigree 104) is masked: its lineage's clone 101 is patent.
episodes_fixture <- function() {
  r <- function(tick, type, human, ped, lin = NA_integer_, dormant = 0L, imported = 0L, parent = -1L,
                n = 1L, detected = -1L, masked = -1L)
    data.table::data.table(tick = as.integer(tick), event_type = type, inoculation_id = -1L,
      mosq_id = -1L, human_id = as.integer(human), lineage_id = as.integer(lin),
      parent1_lineage_id = NA_integer_, parent2_lineage_id = NA_integer_, is_recombinant = 0L, p_outcross = 0,
      recipient = "human", is_dormant = as.integer(dormant), is_imported = as.integer(imported),
      pedigree_id = as.integer(ped), parent1_pedigree_id = as.integer(parent),
      parent2_pedigree_id = -1L, n_hypnozoites = as.integer(n), oocyst_id = -1L,
      detected = as.integer(detected), into_patent_lineage = as.integer(masked))
  ev <- data.table::rbindlist(list(
    r(1,  "vh_inoculation_seed", 2, 201, 9),                                   # founder host 2
    r(5,  "vh_inoculation", 3, 301, 4, detected = 1),                          # host 3: no batch
    r(7,  "vh_inoculation", 5, 501, 7, dormant = 1, n = 2),                    # host 5: batch only
    r(8,  "vh_inoculation", 4, 401, 5, detected = 1),
    r(8,  "vh_inoculation", 4, 402, 5, dormant = 1, n = 1),
    r(10, "vh_inoculation", 1, 101, 1, detected = 1),                          # host 1: two genotypes
    r(10, "vh_inoculation", 1, 102, 2, detected = 1),
    r(10, "vh_inoculation", 1, 103, 1, dormant = 1, n = 2),
    r(12, "bloodcure_given", 3, -1),
    r(12, "bloodcure_subpatent", 3, 301, 4),
    r(15, "human_clearance", 2, 201, 9),
    r(20, "hypnozoite_activation", 1, 104, 1, parent = 103, detected = 1, masked = 1),
    r(22, "recrudescence", 3, 302, 4, parent = 301, detected = 1),
    r(25, "vh_inoculation", 2, 202, 3, detected = 1),
    r(25, "vh_inoculation", 2, 203, 3, dormant = 1, n = 3),
    r(30, "hypnozoite_activation", 5, 502, 7, parent = 501, detected = 1, masked = 0),
    r(30, "human_clearance", 1, 101, 1),
    r(30, "human_clearance", 4, 401, 5),
    r(31, "hypnozoite_clearance", 5, 501, 7, dormant = 1, n = 1),
    r(33, "bloodcure_given", 2, -1),
    r(33, "case_bloodstage_clear", 2, 202, 3),
    r(33, "radcure_given", 2, -1),
    r(33, "radcure_hypnozoite_clear", 2, 203, 3, dormant = 1, n = 3),
    r(35, "human_clearance", 1, 104, 1),
    r(40, "human_clearance", 1, 102, 2),
    r(40, "human_clearance", 3, 302, 4),
    r(45, "hypnozoite_activation", 4, 403, 5, parent = 402, detected = 1, masked = 0),
    r(45, "vh_inoculation", 4, 404, 6, detected = 1),                          # same day: one episode
    r(50, "hypnozoite_activation", 1, 105, 1, parent = 103, detected = 0, masked = 0),
    r(60, "vh_inoculation_import", 5, 503, 8, imported = 1, detected = 1),
    r(61, "vh_inoculation_import", 5, 504, 8, dormant = 1, imported = 1, n = 1), # its batch, a day later
    r(70, "human_clearance", 4, 403, 5),
    r(80, "human_clearance", 5, 502, 7)))
  ev[, event_id := paste0(tick, "_", seq_len(.N)), by = tick]
  data.table::setcolorder(ev, c("tick", "event_id"))
  ev[]
}

## The fixture as a run of 6 hosts and `sim_days` days, with the state table `state` (header only
## by default).
episodes_fixture_run <- function(sim_days = 100L, events = episodes_fixture(), state = NULL)
  load_run(fake_run(params = list(N_HOSTS = 6L, N_MOSQ = 1L, SIM_DAYS = sim_days, N_FOUNDER_HOSTS = 1L),
                    events = events, state = state))

## Host rows for days 1 to `days` from clone spans [start, end): `blood` patent, `sub` subpatent,
## `batch` hypnozoite batches of `n`; `bites` the infectious bites (host, day). One mosquito row a
## day, which every host metric leaves out.
state_from_spans <- function(days, n_hosts, blood = spans(), sub = spans(), batch = spans(),
                             bites = data.table::data.table(entity_id = integer(), tick = integer())) {
  per_day <- function(s) {
    s <- s[start < end]
    if (!nrow(s)) return(data.table::data.table(entity_id = integer(), lineage = integer(), n = integer(),
                                                tick = integer()))
    s[, r := .I][, .(entity_id, lineage, n, tick = seq(start, end - 1L)), by = r][, r := NULL][]
  }
  b <- per_day(blood); u <- per_day(sub); z <- per_day(batch)
  g <- data.table::CJ(entity_id = seq_len(n_hosts), tick = seq_len(days))
  parts <- list(b[, .(moi_active = data.table::uniqueN(lineage)), by = .(entity_id, tick)],
                u[, .(moi_subpatent = data.table::uniqueN(lineage)), by = .(entity_id, tick)],
                z[, .(moi_dormant = data.table::uniqueN(lineage), hyp_reservoir_size = sum(n)), by = .(entity_id, tick)],
                data.table::rbindlist(list(b, u, z))[, .(lineage_ids = paste(unique(lineage), collapse = "|")),
                                                     by = .(entity_id, tick)],
                bites[, .(infectious_bites = .N), by = .(entity_id, tick)])
  for (x in parts) g <- merge(g, x, by = c("entity_id", "tick"), all.x = TRUE)
  for (cc in c("moi_active", "moi_subpatent", "moi_dormant", "hyp_reservoir_size", "infectious_bites"))
    data.table::set(g, which(is.na(g[[cc]])), cc, 0L)
  g[, `:=`(entity_type = "human", tag = as.integer(moi_active > 0 | moi_subpatent > 0),
           n_transmissible = NA_integer_, mosq_age_days = NA_integer_)]
  m <- data.table::data.table(tick = seq_len(days), entity_type = "mosquito", entity_id = 1L, tag = 1L,
                              moi_active = 1L, moi_subpatent = 0L, moi_dormant = 0L, hyp_reservoir_size = 0L,
                              lineage_ids = "9", n_transmissible = 1L, mosq_age_days = 3L,
                              infectious_bites = NA_integer_)
  out <- data.table::rbindlist(list(g[, st_cols, with = FALSE], m[, st_cols, with = FALSE]))
  out[order(tick, entity_type, entity_id)]
}
spans <- function(h = integer(), lin = integer(), start = integer(), end = integer(), n = 0L)
  data.table::data.table(entity_id = as.integer(h), lineage = as.integer(lin), start = as.integer(start),
                         end = as.integer(end), n = rep_len(as.integer(n), length(h)))

## The fixture's host rows, written by hand from its clones (each patent from its onset to the day
## before its end), the subpatent holder and the batches' counts, so that they depend on no function
## under test; and 9 infectious bites, two of them on host 6, whom no bite infected.
fixture_state <- function(sim_days = 100L) {
  e <- sim_days + 1L
  blood <- spans(c(1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5), c(1, 2, 1, 1, 9, 3, 4, 4, 5, 5, 6, 7, 8),
                 c(10, 10, 20, 50, 1, 25, 5, 22, 8, 45, 45, 30, 60), c(30, 40, 35, e, 15, 33, 12, 40, 30, 70, e, 80, e))
  batch <- spans(c(1, 1, 2, 4, 5, 5, 5), c(1, 1, 3, 5, 7, 7, 8), c(10, 20, 25, 8, 7, 30, 61),
                 c(20, 50, 33, 45, 30, 31, e), c(2, 1, 3, 1, 2, 1, 1))
  bites <- data.table::data.table(entity_id = c(3L, 5L, 4L, 1L, 2L, 4L, 6L, 1L, 6L),
                                  tick = c(5L, 7L, 8L, 10L, 25L, 45L, 50L, 60L, 90L))
  state_from_spans(sim_days, 6L, blood, spans(3, 4, 12, 22), batch, bites)
}

## Two hosts and 20 days: one local bite on host 1 on day 2, cleared on day 6, when the run went
## extinct; without `record`, the run stops there with no record of why.
extinct_fixture_run <- function(record = TRUE) {
  ev <- episodes_fixture()[pedigree_id == 301L & event_type == "vh_inoculation"]
  ev[, `:=`(tick = 2L, event_id = "2_1", human_id = 1L, pedigree_id = 11L, lineage_id = 1L)]
  ev <- rbind(ev, data.table::copy(ev)[, `:=`(tick = 6L, event_id = "6_1", event_type = "human_clearance",
                                              detected = -1L)])
  rec <- list(ACTUAL_SEED = 1L, SLIM_VERSION = "5.2", MALARIAFWD_VERSION = "0.1.0")
  if (record) rec$EXTINCT_TICK <- 6L
  load_run(fake_run(params = list(N_HOSTS = 2L, N_MOSQ = 1L, SIM_DAYS = 20L, N_FOUNDER_HOSTS = 1L),
                    record = rec, events = ev,
                    state = state_from_spans(6L, 2L, blood = spans(1, 1, 2, 6),
                                             bites = data.table::data.table(entity_id = 1L, tick = 2L))))
}

## Two hosts and 20 days with no founder and an import window that drew nothing: no event, ten days of
## empty hosts, and the record that the run never started.
never_fixture_run <- function()
  load_run(fake_run(params = list(N_HOSTS = 2L, N_MOSQ = 1L, SIM_DAYS = 20L, N_FOUNDER_HOSTS = 0L,
                                  IMPORT_START = 1L, IMPORT_END = 10L, IMPORT_RATE_PER_TICK = 1e-6),
                    record = list(ACTUAL_SEED = 1L, SLIM_VERSION = "5.2", MALARIAFWD_VERSION = "0.1.0",
                                  NEVER_STARTED = 10L),
                    state = state_from_spans(10L, 2L)))

## A treated run with imports, recrudescence and an MDA round; and one that never started.
.sum_runs <- new.env()
treated_run <- function() {
  if (is.null(.sum_runs$treated))
    .sum_runs$treated <- suppressMessages(run_sim(vivax_params(
      N_HOSTS = 60, N_MOSQ = 150, SIM_DAYS = 400, CHROMOSOMES = "1", BLOODCURE_COVERAGE = 0.5,
      P_BLOODCURE_SURVIVAL = 0.3, P_RECRUD_GIVEN_FAIL = 0.5, MDA_TICK = 200, MDA_COVERAGE = 0.5,
      IMPORT_START = 1, IMPORT_END = 400, IMPORT_RATE_PER_TICK = 0.05), seed = 1, outdir = tempfile("mfwd_treated_")))
  .sum_runs$treated
}
never_run <- function() {
  if (is.null(.sum_runs$never))
    .sum_runs$never <- suppressMessages(run_sim(vivax_params(
      N_HOSTS = 20, N_MOSQ = 20, SIM_DAYS = 60, N_FOUNDER_HOSTS = 0, IMPORT_START = 1, IMPORT_END = 20,
      IMPORT_RATE_PER_TICK = 1e-6, CHROMOSOMES = "1"), seed = 1, outdir = tempfile("mfwd_never_")))
  .sum_runs$never
}
