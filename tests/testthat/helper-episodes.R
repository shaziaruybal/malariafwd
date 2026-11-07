## Runs shared by the window and episode tests, made once per session.

## W goes extinct early for seeds 1 to 6 (days 43 to 196), before its burn-in of 200 days ends.
w_params <- function() vivax_params(N_HOSTS = 60, N_MOSQ = 240, BITING_RATE = 0.12, SIM_DAYS = 400,
                                    N_FOUNDER_HOSTS = 12, NU_HYPNOZOITE_BATCH = 0,
                                    MEAN_BLOOD_STAGE_DURATION = 10, CHROMOSOMES = "1")
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

## The fixture as a run of 6 hosts and `sim_days` days.
episodes_fixture_run <- function(sim_days = 100L, events = episodes_fixture())
  load_run(fake_run(params = list(N_HOSTS = 6L, N_MOSQ = 1L, SIM_DAYS = sim_days, N_FOUNDER_HOSTS = 1L),
                    events = events))
