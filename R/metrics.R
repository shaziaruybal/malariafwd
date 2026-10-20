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

#' Relapse fraction of recurrences (relapses / (relapses + local reinfections))
#' @param events an events data.table (from [load_events()]).
#' @return a one-row data.table.
#' @export
relapse_fraction <- function(events) {
  ev <- as.data.table(events)
  n_relapse  <- ev[event_type == "hypnozoite_activation", .N]
  n_reinfect <- ev[event_type == "vh_inoculation" & is_dormant == 0, .N]
  data.table(n_relapse = n_relapse, n_reinfection = n_reinfect,
             relapse_fraction = n_relapse / (n_relapse + n_reinfect))
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
#' `origin` gives the acquisition mechanism: relapse = hypnozoite activation;
#' imported = a new infection that arrived via importation; local reinfection = a new
#' local mosquito-borne infection. `imported` is a separate genotype-origin flag that is
#' `TRUE` whenever the parasite genotype is import-derived -- including **imported
#' relapses** (a relapse of an imported genotype), a key vivax truth category that `origin`
#' alone ("relapse") does not expose. So the four truth classes are: imported new infection
#' (`origin == "imported"`), local reinfection, imported relapse (`origin == "relapse"` &
#' `imported`), and local relapse (`origin == "relapse"` & `!imported`).
#' @param events an events data.table.
#' @return a data.table with `origin` and `imported` columns (plus case identifiers).
#' @export
classify_cases <- function(events) {
  ev <- as.data.table(events)
  cases <- ev[(event_type %in% c("vh_inoculation","vh_inoculation_import") & is_dormant == 0) |
                event_type == "hypnozoite_activation"]
  cases[, origin := fcase(
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
#' When the run contains imported cases, relapses are split into `relapse (imported)` and
#' `relapse (local)` -- a relapse of an imported genotype is a distinct truth category. Set
#' `split_imported = FALSE` to keep a single `relapse` class, or `TRUE` to force the split.
#' @param run a run directory or `vivax_run`.
#' @param bin window width in days.
#' @param split_imported `"auto"` (split only if the run has imports), `TRUE`, or `FALSE`.
#' @return a data.table of `win`, `origin`, `n`, and `rate` (new clones per host-year).
#' @references Mueller I, et al. (2013) A high force of *Plasmodium vivax* blood-stage
#'   infection drives the rapid acquisition of immunity. PLoS Negl Trop Dis 7(9):e2403.
#' @export
force_of_infection <- function(run, bin = 30, split_imported = "auto") {
  nh <- as.numeric(load_params(run)$N_HOSTS)
  cs <- classify_cases(load_events(run))
  if (isTRUE(split_imported) || (identical(split_imported, "auto") && any(cs$imported)))
    cs[origin == "relapse", origin := fifelse(imported, "relapse (imported)", "relapse (local)")]
  cs[, win := (tick %/% bin) * bin]
  fo <- cs[, .(n = .N), by = .(win, origin)]
  fo[, rate := n / nh / (bin / 365)][order(win, origin)]
}

#' Mean pairwise IBD among cases active in each time window (relatedness trajectory)
#' @param events an events data.table.
#' @param idx an [ibd_index()] object.
#' @param bin window width (days).
#' @param cap max cases sampled per window.
#' @param min_pairs drop windows with fewer pairs.
#' @return a data.table of `win`, `mean_ibd`, `n_pairs`.
#' @export
ibd_trajectory <- function(events, idx, bin = 90, cap = 40, min_pairs = 5) {
  cs <- classify_cases(events)[pedigree_id %in% idx$nm$pedigree_id]
  cs[, win := (tick %/% bin) * bin]
  out <- cs[, {
    peds <- if (.N > cap) sample(pedigree_id, cap) else pedigree_id
    if (length(peds) >= 2) {
      cb <- combn(length(peds), 2)
      ibd <- idx$getibd(peds[cb[1,]], peds[cb[2,]])
      .(mean_ibd = mean(ibd), n_pairs = length(ibd))
    } else .(mean_ibd = NA_real_, n_pairs = 0L)
  }, by = win][n_pairs >= min_pairs]
  out[order(win)]
}

#' Epidemiology-only (WHO-style) import classification of each case
#' @param events an events data.table.
#' @param window travel look-back window (days).
#' @return a data.table with a `surveillance` label to compare against `origin`.
#' @export
classify_who <- function(events, window = 60) {
  cases <- classify_cases(events)
  imports <- as.data.table(events)[event_type == "vh_inoculation_import" & is_dormant == 0,
                                   .(human_id, imp_tick = tick)]
  if (nrow(imports) == 0) return(cases[, surveillance := "Indigenous"][])
  j <- imports[cases, on = "human_id", allow.cartesian = TRUE]
  j[, recent := !is.na(imp_tick) & imp_tick <= tick & imp_tick > tick - window]
  flag <- j[, .(recent_travel = any(recent, na.rm = TRUE)), by = .(tick, human_id, pedigree_id, origin)]
  flag[, surveillance := fifelse(recent_travel, "Imported", "Indigenous")][]
}

#' Directed transmission links (source -> target) via shared oocyst
#' @param events an events data.table.
#' @return a data.table of `source`, `target` pedigree ids.
#' @export
transmission_links <- function(events) {
  ev <- as.data.table(events)
  up <- ev[event_type %in% c("hv_uptake_recombinant","hv_uptake_clonal") & oocyst_id > 0,
           .(oocyst_id, A1 = parent1_pedigree_id, A2 = parent2_pedigree_id)]
  up <- unique(up, by = "oocyst_id")
  inoc <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id > 0, .(oocyst_id, C = pedigree_id)]
  m <- merge(inoc, up, by = "oocyst_id")
  unique(rbindlist(list(m[A1 > 0, .(source = A1, target = C)], m[A2 > 0, .(source = A2, target = C)])))
}

#' Effective reproduction number by true origin
#' @param events an events data.table.
#' @return a data.table of `origin`, `n_sources`, `Rc`, `pct_onward`.
#' @export
rc_by_origin <- function(events) {
  links <- transmission_links(events)
  cases <- classify_cases(events)[, .(pedigree_id, origin)]
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
#' @return a one-row data.table of `n`, `mean_ibd`, `clonal_frac`, `n_lineages`.
#' @export
new_infection_diversity <- function(run, burn_frac = 0.5, cap = 250, clonal = 0.5) {
  ev <- load_events(run)
  endt <- max(ev$tick); burn <- round(endt * burn_frac)
  peds <- ev[event_type == "vh_inoculation" & is_dormant == 0 & tick >= burn, unique(pedigree_id)]
  if (length(peds) < 6)
    return(data.table(n = length(peds), mean_ibd = NA_real_, clonal_frac = NA_real_, n_lineages = NA_integer_))
  if (length(peds) > cap) peds <- sample(peds, cap)
  idx <- ibd_index(ibd_pairs(run, ids = peds))
  vp <- idx$nm$pedigree_id
  cb <- combn(length(vp), 2); di <- idx$getibd(vp[cb[1,]], vp[cb[2,]])
  data.table(n = length(vp), mean_ibd = mean(di), clonal_frac = 100 * mean(di >= clonal),
             n_lineages = uniqueN(ev[pedigree_id %in% vp, lineage_id]))
}
