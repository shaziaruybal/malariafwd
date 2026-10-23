## The four TRUE origins of a blood-stage episode, in a fixed order. Palettes, factor levels
## and grid pooling all key off this, so a run that happens to contain no recrudescences still
## produces the same shape as one that does.
.ORIGINS <- c("imported", "local reinfection", "relapse", "recrudescence")

## Origins that are an ACQUISITION of a new parasite clone. A recrudescence is a new clinical
## EPISODE but not a new clone - the genotype was already in the host - so it belongs in case
## counts and not in a force-of-infection rate.
.ACQUISITIONS <- c("imported", "local reinfection", "relapse")

## The origin vocabulary once relapses are split by the genotype's provenance. An imported
## relapse is a distinct truth category, so the split is worth having -- but it must not depend
## on whether a given run happens to contain imports, or pooling a grid puts the same quantity
## into `relapse`, `relapse (imported)` and `relapse (local)` as three separate columns.
.ORIGINS_SPLIT <- c("imported", "local reinfection",
                    "relapse (imported)", "relapse (local)", "recrudescence")

## Onset marks on a host timeline. Named so an unlisted `kind` is a visible missing entry
## rather than a silent shape = NA that draws nothing.
.KIND_SHAPE <- c("new infection" = 24, "relapse" = 21, "recrudescence" = 23)

## metrics.R -- epidemiological + genetic summaries of a run's event log / state table.
## Inputs are the data.tables from load_events()/load_state() (or pass a vivax_run to the
## run-level helpers). See DECISIONS.md for the identity model (pedigree/lineage/oocyst).

#' Distribution of complexity-of-infection (active MOI) across hosts
#' @param state a state data.table (from [load_state()]).
#' @param burn_in drop ticks before this (steady-state).
#' @return a data.table of `moi_active` counts and percentages.
#' @export
coi_distribution <- function(state, burn_in = 0) {
  h <- as.data.table(state)[entity_type == "human" & tick >= burn_in & moi_active > 0]
  h[, .N, by = moi_active][order(moi_active)][, pct := 100 * N / sum(N)][]
}

## Deterministic subsampling. Both metrics below cap the number of genomes they compare, and
## an unseeded sample() made repeat calls on the SAME run disagree -- new_infection_diversity()
## returned 0.0529 / 0.0544 / 0.0526 / 0.0523 across four identical calls. A metric must not do
## that. The draw stays random (so it is unbiased in time and lineage) but is reproducible, and
## the caller's RNG stream is restored so this cannot perturb a surrounding simulation.
.sample_fixed <- function(x, n, seed) {
  if (length(x) <= n) return(x)
  if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }
  set.seed(seed)
  sample(x, n)
}

#' Relapse fraction, in two units
#'
#' Relapses as a share of blood-stage episodes. Two readings are returned because a
#' recrudescence is a new clinical **episode** but not a newly **acquired** parasite -- the
#' genotype was already in the host (D25):
#'
#' * `relapse_fraction` -- relapses over **acquisitions** (imported + local reinfection +
#'   relapse). This is what the published 0.55-0.85 vivax band refers to, what `epi_summary()`
#'   in the analyses repo scores against, and what `summarise_run()$relapse_frac` reports.
#' * `relapse_frac_episodes` -- relapses over **every patent episode**, recrudescence included.
#'
#' With the shipped default (`P_RECRUD_GIVEN_FAIL = 0`) there are no recrudescences and the two
#' are identical; they separate only once treatment failure can suppress a clone.
#'
#' Both denominators include each host's first-ever episode, which is not a *recurrence*. For a
#' strict recurrences-only fraction, drop the first episode per host:
#' `cs[order(tick)][, .SD[-1], by = human_id][, mean(origin == "relapse")]`. That is a
#' materially different quantity, so never compare it against either reading here.
#'
#' @param events an events data.table (from [load_events()]).
#' @return a one-row data.table of `n_relapse`, `n_reinfection`, `n_recrudescence`,
#'   `relapse_fraction` (over acquisitions) and `relapse_frac_episodes` (over every patent
#'   episode). `summarise_run()$relapse_frac` reports the acquisition reading.
#' @export
relapse_fraction <- function(events) {
  cs <- classify_cases(events)
  n_rel <- cs[origin == "relapse", .N]
  n_new <- cs[origin %in% c("imported", "local reinfection"), .N]
  n_rec <- cs[origin == "recrudescence", .N]
  data.table(n_relapse = n_rel, n_reinfection = n_new, n_recrudescence = n_rec,
             ## the band-comparable number: relapses over ACQUISITIONS
             relapse_fraction = n_rel / max(n_rel + n_new, 1),
             ## every patent episode, recrudescence included
             relapse_frac_episodes = n_rel / max(n_rel + n_new + n_rec, 1))
}

#' Relapses per person-year among ever-infected hosts
#' @param events an events data.table.
#' @return a one-row data.table.
#' @export
relapses_per_person_year <- function(events) {
  ev <- as.data.table(events); endt <- max(ev$tick)
  inf <- ev[event_type %in% c("vh_inoculation","vh_inoculation_import","vh_inoculation_seed") &
              is_dormant == 0, .(t0 = min(tick)), by = human_id]
  rel <- ev[event_type == "hypnozoite_activation", .(nrel = .N), by = human_id]
  m <- merge(inf, rel, by = "human_id", all.x = TRUE)[is.na(nrel), nrel := 0]
  m[, py := (endt - t0) / 365]
  data.table(total_relapses = sum(m$nrel), person_years = round(sum(m$py), 1),
             relapses_per_py = round(sum(m$nrel) / sum(m$py), 2))
}

#' Mean hypnozoite reservoir per host over time
#' @param state a state data.table.
#' @return a data.table of `tick`, `mean_res`.
#' @export
reservoir_trajectory <- function(state) {
  as.data.table(state)[entity_type == "human", .(mean_res = mean(hyp_reservoir_size)), by = tick][order(tick)]
}

#' Hypnozoite dormancy durations (establishment -> activation, in days)
#' @param events an events data.table.
#' @return a data.table with a `duration` column.
#' @export
dormancy_durations <- function(events) {
  ev <- as.data.table(events)
  batch <- ev[event_type %in% c("vh_inoculation","vh_inoculation_import","vh_inoculation_seed") & is_dormant == 1,
              .(pedigree_id, est_tick = tick)]
  act <- ev[event_type == "hypnozoite_activation", .(parent1_pedigree_id, act_tick = tick)]
  m <- merge(act, batch, by.x = "parent1_pedigree_id", by.y = "pedigree_id")
  m[act_tick >= est_tick, .(duration = act_tick - est_tick)]
}

#' Classify each new blood-stage episode ("case") by TRUE origin
#'
#' `origin` gives the mechanism by which the episode arose, and takes **four** values:
#' `relapse` = a hypnozoite activated; `imported` = a new infection that arrived via
#' importation; `local reinfection` = a new local mosquito-borne infection; and
#' `recrudescence` = a clone that survived a treatment course, was suppressed below patency,
#' and re-expanded (D25). Note the fourth is **not** an acquisition -- no new parasite entered
#' the host -- so it is excluded from [force_of_infection()]'s rates and from
#' [new_infection_diversity()]. `imported` is a separate genotype-origin flag that is
#' `TRUE` whenever the parasite genotype is import-derived -- including **imported
#' relapses** (a relapse of an imported genotype), a key vivax truth category that `origin`
#' alone ("relapse") does not expose. So the four truth classes are: imported new infection
#' (`origin == "imported"`), local reinfection, imported relapse (`origin == "relapse"` &
#' `imported`), and local relapse (`origin == "relapse"` & `!imported`).
#' @param events an events data.table.
#' @return a data.table of `tick`, `human_id`, `lineage_id`, `oocyst_id`, `pedigree_id`,
#'   `origin`, `imported` and `detected` -- one row per blood-stage episode. Founder seeding
#'   (`vh_inoculation_seed`) is excluded: it is burn-in, not an acquired case.
#' @export
classify_cases <- function(events) {
  ev <- as.data.table(events)
  cases <- ev[(event_type %in% c("vh_inoculation","vh_inoculation_import") & is_dormant == 0) |
                event_type %in% c("hypnozoite_activation", "recrudescence")]
  cases[, origin := fcase(
    event_type == "recrudescence", "recrudescence",
    event_type == "hypnozoite_activation", "relapse",
    event_type == "vh_inoculation_import" | is_imported == 1, "imported",
    default = "local reinfection")]
  cases[, imported := is_imported == 1 | event_type == "vh_inoculation_import"]
  cases[, .(tick, human_id, lineage_id = as.character(lineage_id), oocyst_id, pedigree_id,
            origin, imported, detected)]
}

#' Binned incidence by origin
#' @param events an events data.table.
#' @param bin window width in days.
#' @param detected_only restrict to the reported line list.
#' @return a data.table of counts per (window, origin).
#' @export
incidence_by_origin <- function(events, bin = 30, detected_only = FALSE) {
  cs <- classify_cases(events)
  if (detected_only) cs <- cs[detected == 1]
  cs[, win := (tick %/% bin) * bin]
  cs[, .N, by = .(win, origin)]
}

#' Force of blood-stage infection, by origin (per host-year)
#'
#' New blood-stage clone acquisitions per host per year, split by origin.
#' `imported` + `local reinfection` give the mosquito-borne **force of infection (FOI)**;
#' adding `relapse` gives the **molecular force of blood-stage infection (molFOB)** -- a
#' compound measure that is elevated in *P. vivax* because relapses contribute new
#' blood-stage clones (Mueller et al. 2013).
#' Relapses are split into `relapse (imported)` and `relapse (local)` -- a relapse of an
#' imported genotype is a distinct truth category. Set `split_imported = FALSE` for a single
#' `relapse` class. Either way `origin` comes back as a factor with **all** levels present, so
#' the shape of the result depends on the argument and never on what the run happened to
#' contain: pooling a grid cannot split one quantity across three columns.
#'
#' A **recrudescence** is returned as its own origin with `is_acquisition = FALSE`. It is a new
#' patent episode but not a new clone, so including it would inflate molFOB -- which counts
#' clone acquisitions. Worth noting that a real molFOB, measured by comparing genotypes between
#' study visits, would often count a recrudescence as a new clone, because the clone was
#' undetectable at the previous visit; the gap between the two sums here estimates that bias.
#' @param run a run directory or `vivax_run`.
#' @param bin window width in days.
#' @param split_imported split relapses by the provenance of the genotype. `TRUE` (default)
#'   gives five origin levels, `FALSE` gives four.
#' @return a data.table of `win`, `origin`, `n`, `days` (the days the window covers, which is
#'   less than `bin` for a trailing window when the run length is not a multiple of `bin`),
#'   `rate` (per host-year, scaled by `days`) and `is_acquisition`. Sum only the
#'   `is_acquisition` rows for FOI or molFOB: a recrudescence is returned as its own row so
#'   incidence stays complete, but it is not a newly acquired clone.
#' @references Mueller I, et al. (2013) A high force of *Plasmodium vivax* blood-stage
#'   infection drives the rapid acquisition of immunity. PLoS Negl Trop Dis 7(9):e2403.
#' @export
force_of_infection <- function(run, bin = 30, split_imported = TRUE) {
  nh <- as.numeric(load_params(run)$N_HOSTS)
  cs <- classify_cases(load_events(run))
  ## Unconditional, and with fixed levels. The old "auto" relabelled relapses only when the run
  ## contained imports, so the returned vocabulary depended on run CONTENT: `baseline` came
  ## back with `relapse` and `imports` with `relapse (local)` for the same thing, and a pooled
  ## dcast over a grid made them separate columns without erroring.
  lv <- if (isTRUE(split_imported)) .ORIGINS_SPLIT else .ORIGINS
  if (isTRUE(split_imported))
    cs[origin == "relapse", origin := fifelse(imported, "relapse (imported)", "relapse (local)")]
  cs[, win := (tick %/% bin) * bin]
  ## Divide each window by the days it actually covers, not by `bin`. Unless the run length is
  ## a multiple of `bin` the trailing window is short - 10 days of a 30-day bin on a 400-day
  ## run - and dividing it by a full bin understated the final point 3x, so every series ended
  ## in a crash that was pure binning artefact, right where a reader looks to judge whether
  ## transmission stopped. The run end comes from the state table rather than the last case,
  ## so a run that ends with no transmission still reports its true trailing zero. Full windows
  ## are unchanged, because there days == bin.
  endt <- max(as.data.table(load_state(run))$tick)
  fo <- cs[, .(n = .N), by = .(win, origin)]
  fo[, days := pmin(bin, endt - win + 1)]
  fo[, rate := n / nh / (days / 365)]
  ## A recrudescence is a new EPISODE but not a new CLONE, so it must not enter a
  ## force-of-infection rate: molFOB (Mueller 2013) counts clone acquisitions. It is returned
  ## as its own row rather than dropped, so incidence stays complete and the caller sums the
  ## rows it means. Note a REAL molFOB, measured by comparing genotypes between visits, would
  ## often count a recrudescence as a new clone - the gap between the two sums here is an
  ## estimate of that bias.
  fo[, is_acquisition := sub(" \\(.*", "", as.character(origin)) %in% .ACQUISITIONS]
  fo[, origin := factor(origin, levels = lv)]
  stopifnot(!anyNA(fo$origin))   # an origin outside `lv` would silently become NA
  fo[order(win, origin)]
}

#' Mean pairwise IBD among cases active in each time window (relatedness trajectory)
#' @param events an events data.table.
#' @param idx an [ibd_index()] object.
#' @param bin window width (days).
#' @param cap max cases sampled per window.
#' @param min_pairs drop windows with fewer pairs.
#' @param seed seed for the per-window subsample when a window exceeds `cap`, so repeat calls
#'   on the same run agree. The caller's RNG stream is left untouched.
#' @return a data.table of `win`, `mean_ibd`, `n_pairs`.
#' @export
ibd_trajectory <- function(events, idx, bin = 90, cap = 40, min_pairs = 5, seed = 1) {
  cs <- classify_cases(events)[pedigree_id %in% idx$nm$pedigree_id]
  cs[, win := (tick %/% bin) * bin]
  out <- cs[, {
    peds <- .sample_fixed(pedigree_id, cap, seed + .GRP)
    if (length(peds) >= 2) {
      cb <- combn(length(peds), 2)
      ibd <- idx$getibd(peds[cb[1,]], peds[cb[2,]])
      .(mean_ibd = mean(ibd), n_pairs = length(ibd))
    } else .(mean_ibd = NA_real_, n_pairs = 0L)
  }, by = win][n_pairs >= min_pairs]
  out[order(win)]
}

#' Epidemiology-only (WHO-style) import classification of each case
#'
#' What a surveillance system would conclude from **travel history alone**, so that it can be
#' scored against the genetic truth in `origin` / `imported`. The two are *expected* to
#' disagree -- that disagreement is the quantity of interest, not an error.
#'
#' `window` is the dominant parameter, and the reason is specific to vivax: a hypnozoite can
#' activate long after the trip. On the audit run the gap from importation to an imported
#' relapse had a median of 43 days but an IQR of 18-95 and a maximum of 325, so a 60-day window
#' catches 59% of imported relapses, 90 days catches 73%, and 180 days 92% -- at the cost of
#' flagging steadily more locally-acquired cases as imported (520 rising to 935). Any travel
#' window is a trade between the two, and no window recovers an imported relapse that activates
#' a year later.
#'
#' @param events an events data.table.
#' @param window travel look-back window (days). See above: this choice dominates the result.
#' @return every column [classify_cases()] returns, plus `recent_travel` and a `surveillance`
#'   label to compare against `origin`. The schema does not depend on whether the run has
#'   imports, so results pool across a grid.
#' @export
classify_who <- function(events, window = 60) {
  cases <- classify_cases(events)
  ## Group by WHATEVER classify_cases() returns, so the output schema is the same table for
  ## every run. Naming four columns here meant an import-free run returned all of them plus
  ## `surveillance` while a run WITH imports returned only the four -- and `imported`, the
  ## truth flag, was present exactly where it is trivially FALSE and absent on the one run
  ## where it carries signal. Safe because pedigree_id is already unique per case row.
  keys <- names(cases)
  ## Any import establishment marks the trip, dormant or not. Requiring is_dormant == 0 meant a
  ## traveller who returned carrying only hypnozoites never registered as having travelled --
  ## 43 of 91 import-exposed hosts on the audit run, and only 59% of imported relapses were
  ## flagged. A surveillance system records the journey, not which compartment the parasite
  ## ended up in.
  imports <- as.data.table(events)[event_type == "vh_inoculation_import",
                                   .(human_id, imp_tick = tick)]
  if (nrow(imports) == 0)
    return(cases[, `:=`(recent_travel = FALSE, surveillance = "Indigenous")][])
  j <- imports[cases, on = "human_id", allow.cartesian = TRUE]
  j[, recent := !is.na(imp_tick) & imp_tick <= tick & imp_tick > tick - window]
  flag <- j[, .(recent_travel = any(recent, na.rm = TRUE)), by = keys]
  flag[, surveillance := fifelse(recent_travel, "Imported", "Indigenous")][]
}

#' Directed transmission links (source -> target) via shared oocyst
#' @param events an events data.table.
#' @return a data.table of `source`, `target` pedigree ids.
#' @export
transmission_links <- function(events) {
  ev <- as.data.table(events)
  up <- ev[event_type %in% c("hv_uptake_recombinant","hv_uptake_clonal") & oocyst_id != -1,
           .(oocyst_id, A1 = parent1_pedigree_id, A2 = parent2_pedigree_id)]
  up <- unique(up, by = "oocyst_id")
  inoc <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id != -1, .(oocyst_id, C = pedigree_id)]
  m <- merge(inoc, up, by = "oocyst_id")
  ## != -1, not > 0: -1 is the no-parent sentinel (malariafwd_vivax.slim ~1397) but pedigree
  ## 0 is a REAL genome -- the first seed founder -- whose onward edges were being dropped.
  ## Small but real: 8 of 1769 edges over the audit run set, and it cost `treat` an entire
  ## four-case transmission chain (7 components reported against 8).
  unique(rbindlist(list(m[A1 != -1, .(source = A1, target = C)],
                        m[A2 != -1, .(source = A2, target = C)])))
}

#' Effective reproduction number by true origin
#'
#' @details
#' `unit` chooses what a "source" is, and the two readings answer different questions.
#' `"episode"` credits onward transmission to the clinical episode it came from, so a
#' recrudescence is its own source -- the surveillance-targeting reading. `"acquisition"`
#' walks each transmitting genome back to the infection that originally brought it into the
#' host, so one acquisition gets one R -- the classic reading, and the one comparable to an
#' Rc estimated from routine case counts (Champagne et al. 2022). They differ substantially:
#' measured on a prototype run, local-reinfection Rc was 0.520 per acquisition against 0.116
#' per episode, with the difference re-attributed to a recrudescence class at 1.464.
#' @param events an events data.table.
#' @param unit `"episode"` (default) or `"acquisition"`; see Details.
#' @return a data.table of `origin`, `n_sources`, `Rc`, `pct_onward`.
#' @export
rc_by_origin <- function(events, unit = c("episode", "acquisition")) {
  unit <- match.arg(unit)
  ev <- as.data.table(events)
  links <- transmission_links(ev)
  cases <- classify_cases(ev)[, .(pedigree_id, origin)]
  if (unit == "acquisition") {
    ## Walk each transmitting genome back to the ACQUISITION it descends from, so a clone
    ## that was suppressed and recrudesced is credited to the infection that brought it in
    ## rather than to the recurrence. This is the unit a surveillance Rc estimated from
    ## routine case counts refers to (Champagne et al. 2022, Swiss TPH), so it is the one to
    ## compare against published values.
    rec <- ev[event_type == "recrudescence", .(ped = pedigree_id, holder = parent1_pedigree_id)]
    if (nrow(rec)) {
      m <- stats::setNames(rec$holder, as.character(rec$ped))
      root <- function(p) { k <- 0L
        while (!is.na(m[as.character(p)]) && k < 100L) { p <- unname(m[as.character(p)]); k <- k + 1L }
        p }
      links[, source := vapply(source, root, numeric(1))]
      cases <- cases[origin != "recrudescence"]
    }
  }
  outdeg <- links[, .(secondary = .N), by = source]
  rc <- merge(cases, outdeg, by.x = "pedigree_id", by.y = "source", all.x = TRUE)
  rc[is.na(secondary), secondary := 0]
  rc[, .(n_sources = .N, Rc = round(mean(secondary), 3),
         pct_onward = round(100 * mean(secondary > 0), 1)), by = origin]
}

#' Transmission-chain sizes (connected components of the source->target graph)
#' @param events an events data.table.
#' @return a data.table with a `chain_size` column.
#' @export
transmission_chains <- function(events) {
  links <- transmission_links(events)
  if (nrow(links) == 0) return(data.table(chain_size = integer(0)))
  nodes <- unique(c(links$source, links$target)); id <- setNames(seq_along(nodes), as.character(nodes))
  p <- seq_along(nodes)
  root <- function(x) { while (p[x] != x) x <- p[x]; x }
  s <- id[as.character(links$source)]; t <- id[as.character(links$target)]
  for (i in seq_along(s)) { ra <- root(s[i]); rb <- root(t[i]); if (ra != rb) p[ra] <- rb }
  comp <- vapply(seq_along(nodes), root, integer(1))
  data.table(node = nodes, comp = comp)[, .(chain_size = .N), by = comp][, .(chain_size)]
}

#' Host-level time to first relapse (for Kaplan-Meier), censored at end of sim
#' @param events an events data.table.
#' @return a data.table of `human_id`, `time`, `event`.
#' @export
time_to_relapse <- function(events) {
  ev <- as.data.table(events)
  inf <- ev[event_type %in% c("vh_inoculation","vh_inoculation_import","vh_inoculation_seed") &
              is_dormant == 0, .(t0 = min(tick)), by = human_id]
  rel <- ev[event_type == "hypnozoite_activation", .(tr = min(tick)), by = human_id]
  m <- merge(inf, rel, by = "human_id", all.x = TRUE)
  endt <- max(ev$tick)
  m[, `:=`(time = fifelse(!is.na(tr), tr - t0, endt - t0), event = as.integer(!is.na(tr)))]
  m[time >= 0]
}

#' Peak COI vs relapse count per host
#' @param state a state data.table.
#' @param events an events data.table.
#' @return a data.table of `human_id`, `peak_moi`, `n_relapse`.
#' @export
moi_vs_relapse <- function(state, events) {
  st <- as.data.table(state)[entity_type == "human"]
  peak <- st[, .(peak_moi = max(moi_active)), by = .(human_id = as.integer(entity_id))]
  rel <- as.data.table(events)[event_type == "hypnozoite_activation", .(n_relapse = .N), by = human_id]
  m <- merge(peak, rel, by = "human_id", all.x = TRUE)
  m[is.na(n_relapse), n_relapse := 0][peak_moi > 0]
}

#' Population host-state prevalence over time (the 4-way surveillance-blind-spot split)
#' @param state a state data.table.
#' @return a data.table of `tick`, `cls`, `N`, `frac`.
#' @export
host_state_prevalence <- function(state) {
  st <- as.data.table(state)[entity_type == "human"]
  st[, cls := fcase(
    moi_active > 0 & moi_dormant > 0, "active + reservoir",
    moi_active > 0 & moi_dormant == 0, "active only",
    moi_active == 0 & moi_dormant > 0, "reservoir only (invisible)",
    default = "uninfected")]
  st[, .N, by = .(tick, cls)][, frac := N / sum(N), by = tick][]
}

#' Pick a host with a rich relapse history (for a trajectory illustration)
#' @param events an events data.table.
#' @return a single host id.
#' @export
pick_host <- function(events) {
  as.data.table(events)[event_type == "hypnozoite_activation", .N, by = human_id][order(-N)][1]$human_id
}

#' One-call calibration metrics for a run (steady state = 2nd half)
#' @param run a `vivax_run` or run directory.
#' @param burn_frac steady-state fraction (drop the first `burn_frac` of the run).
#' @return a one-row data.table of prevalence, COI, polyclonality, reservoir, relapse metrics.
#' @export
summarise_run <- function(run, burn_frac = 0.5) {
  st <- load_state(run); ev <- load_events(run)
  endt <- max(st$tick); burn <- round(endt * burn_frac)
  h <- as.data.table(st)[entity_type == "human" & tick >= burn]
  cd <- coi_distribution(st, burn_in = burn)
  vi <- as.data.table(ev)[event_type %in% c("vh_inoculation","vh_inoculation_import") & is_dormant == 0]
  data.table(
    prev = round(mean(h$moi_active > 0), 3),
    COI = round(if (nrow(cd)) sum(cd$moi_active * cd$pct/100) else NA_real_, 2),
    poly = round(if (nrow(cd)) sum(cd[moi_active > 1]$pct) else NA_real_, 1),
    mean_res = round(mean(h$hyp_reservoir_size), 1),
    relapse_frac = round(relapse_fraction(ev)$relapse_fraction, 3),
    relapses_py = relapses_per_person_year(ev)$relapses_per_py,
    ## courses GIVEN, not clones cleared: one course clears every clone in the host, and a
    ## clone can survive it (bloodcure_persist), so counting clearances is wrong in both
    ## directions - and counts zero for MDA, which clears via mda_bloodstage_clear.
    n_bloodstage_treated = as.data.table(ev)[
      event_type %in% c("bloodcure_given", "mda_bloodcure_given"), .N],
    imported_frac = round(if (nrow(vi)) mean(vi$is_imported == 1) else 0, 3))
}

#' Relapse-vs-origin IBD: each relapse paired to a same-lineage prior primary
#' @param events an events data.table.
#' @param idx an [ibd_index()] object.
#' @return a data.table with an `ibd` column (should be ~1, clonal).
#' @export
ibd_relapse_origin <- function(events, idx) {
  ev <- as.data.table(events); nm <- idx$nm
  rel <- ev[event_type == "hypnozoite_activation",
            .(lineage_id = as.character(lineage_id), human_id, rp = pedigree_id, rt = tick)]
  ## every inoculation route, not just local bites: a relapse whose origin was an IMPORT
  ## was dropped by the merge below without warning, which silently excluded the pairs
  ## that are now the cleanest clonal check in the model.
  prm <- ev[event_type %in% c("vh_inoculation", "vh_inoculation_import",
                              "vh_inoculation_seed") & is_dormant == 0,
            .(lineage_id = as.character(lineage_id), human_id, pp = pedigree_id, pt = tick)]
  prm <- prm[pp %in% nm$pedigree_id][order(pt)][, .SD[1], by = .(lineage_id, human_id)]
  pr <- merge(rel[rp %in% nm$pedigree_id], prm, by = c("lineage_id","human_id"))[pt <= rt]
  pr[, ibd := idx$getibd(rp, pp)][]
}

#' New-infection genetic diversity (relatedness of the transmitting population)
#'
#' Pairwise IBD among non-relapse blood-stage infections in the steady-state window;
#' clonal fraction rises / lineage richness falls as transmission falls (the surveillance
#' signal). Requires the tree sequence (IBD).
#' @param run a `vivax_run` or run directory.
#' @param burn_frac steady-state fraction.
#' @param cap max cases sampled.
#' @param clonal IBD threshold for a "clonal" pair.
#' @param seed seed for the subsample when the pool exceeds `cap`, so repeat calls on the same
#'   run agree. The caller's RNG stream is left untouched.
#' @return a one-row data.table of `n`, `mean_ibd`, `clonal_frac`, `n_lineages`.
#' @export
new_infection_diversity <- function(run, burn_frac = 0.5, cap = 250, clonal = 0.5, seed = 1) {
  ev <- load_events(run)
  endt <- max(ev$tick); burn <- round(endt * burn_frac)
  ## Every non-relapse acquisition route, derived from classify_cases() rather than from a
  ## literal event type. Matching "vh_inoculation" alone dropped imported genomes, which
  ## overstated mean IBD and clonal fraction and understated lineage richness by 44% at a 17%
  ## import share and 4.3x at 61% -- all three biases pointing the same way as a genuine
  ## transmission decline, i.e. in the direction of the signal this metric is read for.
  ## acquisitions only. `origin != "relapse"` would now admit recrudescences, which are not
  ## newly acquired genomes - they were already in the host - and would bias this toward the
  ## clonality it exists to measure.
  peds <- classify_cases(ev)[origin %in% c("imported", "local reinfection") & tick >= burn,
                             unique(pedigree_id)]
  if (length(peds) < 6)
    return(data.table(n = length(peds), mean_ibd = NA_real_, clonal_frac = NA_real_, n_lineages = NA_integer_))
  peds <- .sample_fixed(peds, cap, seed)
  idx <- ibd_index(ibd_pairs(run, ids = peds))
  vp <- idx$nm$pedigree_id
  cb <- combn(length(vp), 2); di <- idx$getibd(vp[cb[1,]], vp[cb[2,]])
  data.table(n = length(vp), mean_ibd = mean(di), clonal_frac = 100 * mean(di >= clonal),
             n_lineages = uniqueN(ev[pedigree_id %in% vp, lineage_id]))
}
