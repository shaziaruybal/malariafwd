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

## Host states, in the order host_state_prevalence() gives its levels, each with its fill. Named
## for the same reason: an unlisted state would draw unfilled.
.HOST_STATE_FILL <- c("uninfected" = "grey85", "active only" = "goldenrod2", "active + reservoir" = "firebrick",
                      "subpatent" = "#8e6bbf", "reservoir only (invisible)" = "turquoise4")

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
## that. The draw stays random (so it is unbiased in time and lineage) but is reproducible.
.sample_fixed <- function(x, n, seed) if (length(x) <= n) x else .with_seed(seed, sample(x, n))

## `expr`, evaluated with the RNG seeded by `seed` under one fixed set of kinds, so that a seed gives
## the same draw whatever kinds the caller uses (a future or furrr worker runs L'Ecuyer-CMRG). The
## caller's kinds and stream, or its having none, are restored, so this cannot perturb a
## surrounding simulation.
.with_seed <- function(seed, expr) {
  kinds <- RNGkind()
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    suppressWarnings(RNGkind(kinds[1], kinds[2], kinds[3]))
    if (had) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  expr
}

#' Relapse fraction, by episodes and by acquisitions
#'
#' Relapses as a share of a run's blood-stage cases, in two readings, both counted per clone onset
#' (a row of [classify_cases()], so a bite that delivers two genotypes counts twice):
#'
#' * `relapse_frac_episodes` -- relapses over **every** case: imported, local reinfection, relapse
#'   and recrudescence.
#' * `relapse_frac_acquisitions` -- relapses over the **acquisitions**: imported, local reinfection
#'   and relapse. A recrudescence is a new clinical episode but not a newly acquired parasite (its
#'   genotype was already in the host), so it is left out. This is the unit of the published
#'   0.55-0.85 vivax band.
#'
#' The two are equal when a run has no recrudescence, as at the default `P_RECRUD_GIVEN_FAIL = 0`.
#' Both count each host's first case, which is not a recurrence; [recurrences()] gives the
#' recurrences alone.
#'
#' Masked relapses, into a patent clone of their own lineage (see [classify_cases()]), are counted
#' unless `exclude_masked = TRUE` leaves them out of the numerator and both denominators, as a study
#' that sees only the relapses that bring a new clone would. Neither is the band's own reading:
#' `NU_HYPNOZOITE_BATCH` was fitted with every relapse counted (0.66 at its shipped value of 3,
#' inside the band), and leaving the masked relapses out gives 0.39, below it (MODEL_DESCRIPTION.md,
#' section 13).
#'
#' @param events an events table from [load_events()], or a run (a `vivax_run` or its directory).
#' @param exclude_masked if `TRUE`, leave out the masked relapses.
#' @return a one-row data.table: `n_relapse`, `n_new_infection` (imported and local reinfection),
#'   `n_recrudescence`, `relapse_frac_episodes` and `relapse_frac_acquisitions`, each fraction `NA`
#'   when its denominator is 0 (a run with no case).
#' @seealso [classify_cases()], [summarise_run()]
#' @examples
#' \dontrun{relapse_fraction(load_events(run))}
#' @export
relapse_fraction <- function(events, exclude_masked = FALSE) {
  cs <- classify_cases(events, exclude_masked)
  n_rel <- cs[origin == "relapse", .N]
  n_new <- cs[origin %in% c("imported", "local reinfection"), .N]
  n_rec <- cs[origin == "recrudescence", .N]
  frac <- function(n, d) if (d > 0) n / d else NA_real_
  data.table(n_relapse = n_rel, n_new_infection = n_new, n_recrudescence = n_rec,
             relapse_frac_episodes = frac(n_rel, n_rel + n_new + n_rec),
             relapse_frac_acquisitions = frac(n_rel, n_rel + n_new))
}

## Relapses per person-year among the hosts ever infected, on days `from` to `to` (`NULL`: SIM_DAYS),
## which summarise_run() gives as relapses_per_py_at_risk (D66): every hypnozoite_activation row over
## the person-time from the day after each host's first establishment, a blood-stage primary or a
## hypnozoite batch from a bite or an import, to `to`, also after an extinction; founder seeding does
## not start it, and masked relapses are counted. One row: n_relapse, n_hosts (with person-time in
## the window), person_years and relapses_per_py (NA with no person-time).
relapses_per_person_year <- function(run, from = 0, to = NULL) {
  sim_days <- .sim_days(run, "relapses_per_person_year")
  from <- .check_day(from, "from")
  to <- if (is.null(to)) sim_days else .check_day(to, "to")
  if (to > sim_days)
    stop("`to` = ", to, " is after the run's last day, SIM_DAYS = ", sim_days, ".", call. = FALSE)
  if (from > to) stop("`from` = ", from, " is after `to` = ", to, ".", call. = FALSE)
  .relapses_per_person_year(load_events(run), sim_days, from, to)
}

## The core of relapses_per_person_year() (D66), for callers that hold the events. A host counts
## the days max(t0 + 1, from) to `to`, t0 its first establishment; masked relapses leave the
## numerator only.
.relapses_per_person_year <- function(events, sim_days, from = 0L, to = sim_days, exclude_masked = FALSE) {
  ev <- .as_events(events)
  est <- ev[event_type %in% setdiff(.inoculation_events(), "vh_inoculation_seed") & is_dormant %in% 0:1,
            .(human_id = as.integer(human_id), t0 = as.integer(tick))]
  est <- unique(data.table::setorder(est, human_id, t0), by = "human_id")   # each host's first
  days <- pmax(0, to - pmax(est$t0 + 1, from) + 1)
  n_rel <- classify_cases(ev, exclude_masked)[origin == "relapse" & tick >= from & tick <= to, .N]
  py <- sum(days) / 365
  data.table(n_relapse = n_rel, n_hosts = sum(days > 0), person_years = py,
             relapses_per_py = if (py > 0) n_rel / py else NA_real_)
}

#' Hypnozoite dormancy durations (establishment -> activation, in days)
#' @param events an events data.table.
#' @return a data.table with a `duration` column.
#' @seealso [time_to_relapse()]
#' @export
dormancy_durations <- function(events) {
  ev <- as.data.table(events)
  batch <- ev[event_type %in% .inoculation_events() & is_dormant == 1,
              .(pedigree_id, est_tick = tick)]
  act <- ev[event_type == "hypnozoite_activation", .(parent1_pedigree_id, act_tick = tick)]
  m <- merge(act, batch, by.x = "parent1_pedigree_id", by.y = "pedigree_id")
  m[act_tick >= est_tick, .(duration = act_tick - est_tick)]
}

#' Classify each blood-stage clone onset ("case") by its true origin
#'
#' One row per clone that became patent in a host, in the order of the event log: a bite or an
#' import that gave a blood-stage primary, a relapse, or a recrudescence. Clones that became patent
#' in one host on one day are one episode, sharing an `episode_id` (`"<human_id>:<tick>"`, the same
#' in every function and subset); a bite that delivered two genotypes gives two rows of one
#' episode. Founder seeding (`vh_inoculation_seed`) is not a case: the founders are where the run
#' starts. Detection is drawn per bite, activation, recrudescence or import, so the clones of one
#' episode can differ in `detected`.
#'
#' `origin` gives the mechanism by which the clone arose, and takes **four** values:
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
#'
#' A relapse into a patent clone of its own lineage, born on an earlier day, is *masked*: a study
#' that genotypes the host sees no new clone (MODEL_DESCRIPTION section 5). `into_patent_lineage`
#' marks it, and `exclude_masked = TRUE` leaves it out.
#' @param events an events table from [load_events()], or a run (a `vivax_run` or its directory).
#' @param exclude_masked if `TRUE`, leave out the masked relapses. It needs the flag on every
#'   relapse, which runs made before it was recorded lack.
#' @return a data.table, one row per clone onset: `tick` (the day it became patent), `human_id`,
#'   `episode_id`, `lineage_id`, `oocyst_id`, `pedigree_id`, `origin`, `imported`, `detected` (1
#'   or 0) and `into_patent_lineage` (`TRUE` for a masked relapse, `FALSE` for every other row,
#'   `NA` for a relapse whose run does not record the flag).
#' @seealso [episodes()] for every clone's span, [recurrences()] for each episode against the
#'   earlier ones.
#' @export
classify_cases <- function(events, exclude_masked = FALSE) {
  .check_flag(exclude_masked, "exclude_masked")
  ev <- .as_events(events)
  cases <- ev[event_type %in% setdiff(.onset_events(), "vh_inoculation_seed") & is_dormant == 0]
  cases[, origin := fcase(
    event_type == "recrudescence", "recrudescence",
    event_type == "hypnozoite_activation", "relapse",
    event_type == "vh_inoculation_import" | is_imported == 1, "imported",
    default = "local reinfection")]
  cases[, imported := is_imported == 1 | event_type == "vh_inoculation_import"]
  ## a relapse into a patent clone of its own lineage (D40); NA where the run does not record it
  flag <- if ("into_patent_lineage" %in% names(cases)) cases$into_patent_lineage else rep(NA_integer_, nrow(cases))
  cases[, into_patent_lineage := fifelse(event_type == "hypnozoite_activation", flag == 1L, FALSE)]
  if (exclude_masked) {
    gap <- cases[event_type == "hypnozoite_activation" & is.na(into_patent_lineage), .N]
    if (gap)
      stop("`exclude_masked = TRUE` needs into_patent_lineage on every hypnozoite_activation row, and ", gap,
           " of these events' activation rows lack it.", call. = FALSE)
    cases <- cases[into_patent_lineage == FALSE]
  }
  cases[, .(tick = as.integer(tick), human_id = as.integer(human_id),
            episode_id = paste0(human_id, ":", tick), lineage_id = as.character(lineage_id),
            oocyst_id = as.integer(oocyst_id), pedigree_id = as.integer(pedigree_id),
            origin = as.character(origin), imported = as.logical(imported),
            detected = as.integer(detected), into_patent_lineage = as.logical(into_patent_lineage))]
}

#' Incidence and force of blood-stage infection, by origin
#'
#' The cases of each window of the run, by origin, as counts (`n`, the incidence) and per host-year
#' (`rate`). Among the new clones, `imported` + `local reinfection` give the mosquito-borne **force
#' of infection (FOI)**; adding `relapse` gives the **molecular force of blood-stage infection
#' (molFOB)** -- a compound measure that is elevated in *P. vivax* because relapses contribute new
#' blood-stage clones (Koepfli et al. 2013). With `detected_only = TRUE` only the cases the
#' observation layer reported are counted: the incidence, FOI and molFOB of a line list.
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
#' @param exclude_masked if `TRUE`, leave out the masked relapses, into a patent clone of their
#'   own lineage (see [classify_cases()]): molFOB as a genotyping study measures it.
#' @param detected_only if `TRUE`, count only the reported cases (`detected == 1`). Detection is
#'   drawn once per bite, import, activation or recrudescence (`P_DETECT_PRIMARY`,
#'   `P_DETECT_RELAPSE`, `P_DETECT_RECRUD`).
#' @return a data.table with one row for every window of the run and every origin, by window:
#'   `win` (the window's first day, an integer: 0, `bin`, `2 * bin`, ...), `origin`, `n` (the
#'   cases, clone onsets as rows of [classify_cases()], 0 where there are none; [episodes()]
#'   groups them by host and day), `days` (the days the window covers within the run, which is
#'   `SIM_DAYS` days long whether or not it went extinct: `bin - 1` for the first window, which
#'   starts on day 1, and fewer for the last when the run ends inside it),
#'   `rate` (per host-year, scaled by `days`) and `is_acquisition`. Sum only the
#'   `is_acquisition` rows for FOI or molFOB: a recrudescence is returned as its own row so
#'   incidence stays complete, but it is not a newly acquired clone.
#' @references Koepfli C, et al. (2013) A high force of *Plasmodium vivax* blood-stage
#'   infection drives the rapid acquisition of immunity in Papua New Guinean children. PLoS Negl
#'   Trop Dis 7(9):e2403.
#' @seealso [plot_force_of_infection()]
#' @export
force_of_infection <- function(run, bin = 30, split_imported = TRUE, exclude_masked = FALSE,
                               detected_only = FALSE) {
  bin <- .check_bin(bin)
  .check_flag(exclude_masked, "exclude_masked")
  .check_flag(detected_only, "detected_only")
  nh <- as.numeric(.params_of(run)$N_HOSTS)
  cs <- classify_cases(load_events(run), exclude_masked)
  if (detected_only) cs <- cs[detected == 1L]
  ## Unconditional, and with fixed levels. The old "auto" relabelled relapses only when the run
  ## contained imports, so the returned vocabulary depended on run CONTENT: `baseline` came
  ## back with `relapse` and `imports` with `relapse (local)` for the same thing, and a pooled
  ## dcast over a grid made them separate columns without erroring.
  lv <- if (isTRUE(split_imported)) .ORIGINS_SPLIT else .ORIGINS
  if (isTRUE(split_imported))
    cs[origin == "relapse", origin := fifelse(imported, "relapse (imported)", "relapse (local)")]
  stopifnot(all(cs$origin %in% lv))   # an origin outside `lv` would be lost from the grid
  cs[, win := as.integer((tick %/% bin) * bin)]
  ## Every window of the run, for every origin, with n = 0 where there is none: a window with no
  ## case is a zero of the series, not a gap a line or an area would interpolate across. Each
  ## window is divided by the days it covers within the run, which ends at SIM_DAYS whether or not
  ## it went extinct: the first window holds days 1 to bin - 1 and the last ends at SIM_DAYS, so a
  ## series neither starts nor ends in a dip that is pure binning artefact. Full windows are
  ## unchanged, because there days == bin.
  grid <- .window_grid(.sim_days(run, "force_of_infection"), bin)
  fo <- merge(data.table::CJ(win = grid$win, origin = lv), cs[, .(n = .N), by = .(win, origin)],
              by = c("win", "origin"), all.x = TRUE)
  fo[is.na(n), n := 0L]
  fo <- merge(fo, grid, by = "win")
  fo[, rate := n / nh / (days / 365)]
  ## A recrudescence is a new EPISODE but not a new CLONE, so it must not enter a
  ## force-of-infection rate: molFOB (Koepfli et al. 2013) counts clone acquisitions. It is returned
  ## as its own row rather than dropped, so incidence stays complete and the caller sums the
  ## rows it means. Note a REAL molFOB, measured by comparing genotypes between visits, would
  ## often count a recrudescence as a new clone - the gap between the two sums here is an
  ## estimate of that bias.
  fo[, is_acquisition := sub(" \\(.*", "", origin) %in% .ACQUISITIONS]
  fo[, origin := factor(origin, levels = lv)]
  fo[order(win, origin), .(win, origin, n, days, rate, is_acquisition)]
}

## The windows, `bin` days wide, of a run of `sim_days` days that cover days `from` to `to`: each
## window's first day (0, bin, 2 bin, ...; day t falls in window (t %/% bin) * bin) and the days of
## [from, to] it covers. The run's first window holds days 1 to bin - 1.
.window_grid <- function(sim_days, bin, from = 1L, to = sim_days) {
  win <- seq((from %/% bin) * bin, (to %/% bin) * bin, by = bin)
  data.table::data.table(win = as.integer(win), days = as.integer(pmin(win + bin - 1, to) - pmax(win, from) + 1))
}

.check_bin <- function(bin) {
  if (!is.numeric(bin) || length(bin) != 1L || is.na(bin) || !is.finite(bin) || bin < 1 || bin != round(bin))
    stop("`bin` must be a whole number of days, at least 1, got ", paste(format(bin), collapse = ", "), ".",
         call. = FALSE)
  as.integer(bin)
}

#' Mean pairwise IBD among cases active in each time window (relatedness trajectory)
#' @param events an events data.table.
#' @param ibd a `vivax_ibd` object from [ibd_pairs()], of every pair among the run's cases.
#' @param bin window width (days).
#' @param cap max cases sampled per window.
#' @param min_pairs drop windows with fewer pairs.
#' @param seed seed for the per-window subsample when a window exceeds `cap`, so repeat calls
#'   on the same run agree. The caller's RNG stream is left untouched.
#' @return a data.table of `win`, `mean_ibd`, `n_pairs`.
#' @export
ibd_trajectory <- function(events, ibd, bin = 90, cap = 40, min_pairs = 5, seed = 1) {
  if (!inherits(ibd, "vivax_ibd") || !isTRUE(ibd$all_pairs))
    stop("`ibd` must be a vivax_ibd object of every pair among its genomes, from ibd_pairs().", call. = FALSE)
  x <- ibd
  cs <- classify_cases(events)[pedigree_id %in% x$samples$pedigree_id]
  cs[, win := (tick %/% bin) * bin]
  out <- cs[, {
    peds <- .sample_fixed(pedigree_id, cap, seed + .GRP)
    if (length(peds) >= 2) {
      cb <- combn(length(peds), 2)
      v <- ibd_lookup(x, peds[cb[1,]], peds[cb[2,]])
      .(mean_ibd = mean(v), n_pairs = length(v))
    } else .(mean_ibd = NA_real_, n_pairs = 0L)
  }, by = win][n_pairs >= min_pairs]
  out[order(win)]
}

#' Each case as a travel rule classes it: recent travel or not
#'
#' What a surveillance system that knows only travel history concludes: a case is
#' `"imported (recent travel)"` when its host brought an infection in during the `window` days up to
#' the case, and `"local (no recent travel)"` otherwise. Score it against the truth, the `imported`
#' column of [classify_cases()] or the `who_class` of [classify_chain()]. The two are expected to
#' disagree: the disagreement is what a travel rule costs.
#'
#' The trip is recorded when an import establishes in the host, as a blood-stage infection or as
#' hypnozoites only. The window decides most of the result, and for vivax no window is right: a
#' relapse of an imported batch can come long after the trip, so a short window misses the late
#' relapses, while a long one flags more of the cases the host acquired locally after the trip.
#'
#' @param events an events table from [load_events()], or a run (a `vivax_run` or its directory).
#' @param window the look-back, a whole number of days, 0 or more: a case on day `tick` is flagged
#'   when its host had an import on a day from `tick - window + 1` to `tick`. An import's own case
#'   is flagged for any `window` of 1 or more; at 0 no case is.
#' @return a data.table, one row per row of [classify_cases()] and in its order: its columns, then
#'   `recent_travel` (logical) and `surveillance`, `"imported (recent travel)"` or
#'   `"local (no recent travel)"`. A run without imports gives the same columns, every case local.
#' @seealso [classify_chain()], [classify_cases()]
#' @examples
#' \dontrun{
#' tr <- classify_travel(run, window = 90)
#' table(tr$surveillance, tr$imported)
#' }
#' @export
classify_travel <- function(events, window = 60) {
  ev <- .as_events(events)
  window <- .check_day(window, "window")
  cases <- classify_cases(ev)
  ## Group by WHATEVER classify_cases() returns, so the output schema is the same table for
  ## every run. Naming four columns here meant an import-free run returned all of them plus
  ## `surveillance` while a run WITH imports returned only the four -- and `imported`, the
  ## truth flag, was present exactly where it is trivially FALSE and absent on the one run
  ## where it carries signal. Safe because pedigree_id is already unique per case row.
  keys <- names(cases)
  ## Any import establishment marks the trip, dormant or not. Requiring is_dormant == 0 meant a
  ## traveller who returned carrying only hypnozoites never registered as having travelled --
  ## 43 of 91 import-exposed hosts on the audit run. A surveillance system records the journey,
  ## not which compartment the parasite ended up in.
  imports <- ev[event_type == "vh_inoculation_import", .(human_id, imp_tick = tick)]
  if (nrow(imports) == 0)
    return(cases[, `:=`(recent_travel = FALSE, surveillance = "local (no recent travel)")][])
  j <- imports[cases, on = "human_id", allow.cartesian = TRUE]
  j[, recent := !is.na(imp_tick) & imp_tick <= tick & imp_tick > tick - window]
  flag <- j[, .(recent_travel = any(recent, na.rm = TRUE)), by = keys]
  flag[, surveillance := fifelse(recent_travel, "imported (recent travel)", "local (no recent travel)")][]
}

#' Time from an infection to the first relapse of its hypnozoite batch
#'
#' One row per hypnozoite batch: the hypnozoites one genotype of a bite or an import laid down in a
#' host, so a bite that delivers two genotypes can lay down two. A batch's clock starts on the day it
#' was established, the day of the bite, or for an import that also brought a blood-stage infection,
#' the day after it. Founder hosts carry no batch until a bite or an import gives them one.
#'
#' The clock stops at the batch's first relapse (`event = 1`). Without one (`event = 0`) it stops
#' on the day the batch's last hypnozoite left it, by death or radical cure (`emptied =
#' "censor"`), or at the run's last day, `SIM_DAYS`, for a batch that still held hypnozoites then,
#' and with `emptied = "follow"` for one that emptied without relapsing too. Censoring at emptying
#' gives the time to relapse of the hypnozoites that live; following to the end gives the curve a
#' trial reports, where radical cure shows as relapses prevented.
#'
#' With `exclude_masked = TRUE` a masked relapse, into a patent clone of its own lineage (see
#' [classify_cases()]), does not stop the clock, although its hypnozoites leave the batch: the
#' first relapse is the first a study could see as a recurrence.
#'
#' @param run a `vivax_run` or its directory.
#' @param emptied `"censor"` (stop the clock on the day a batch empties without relapsing) or
#'   `"follow"` (follow it to `SIM_DAYS`).
#' @param exclude_masked if `TRUE`, the first relapse is the first that is not masked.
#' @return a data.table, one row per batch, by `batch_tick` and `batch_pedigree_id`:
#'   `batch_pedigree_id`, `human_id`, `lineage_id`, `imported` (an imported batch), `batch_tick`
#'   (the day it was established), `n_hypnozoites` (its size then), `time` (days from `batch_tick`
#'   to the first relapse or to censoring), `event` (1 relapse, 0 censored) and `end_reason`
#'   (`"relapse"`, `"emptied"` or `"run_end"`). A run with no batch gives no row.
#' @seealso [plot_relapse_km()], [dormancy_durations()]
#' @examples
#' \dontrun{plot_relapse_km(time_to_relapse(run))}
#' @export
time_to_relapse <- function(run, emptied = c("censor", "follow"), exclude_masked = FALSE) {
  emptied <- match.arg(emptied); .check_flag(exclude_masked, "exclude_masked")
  end <- .sim_days(run, "time_to_relapse"); ev <- load_events(run)
  b <- ev[event_type %in% .inoculation_events() & is_dormant == 1,
          .(batch_pedigree_id = as.integer(pedigree_id), human_id = as.integer(human_id),
            lineage_id = as.character(lineage_id), imported = is_imported == 1,
            batch_tick = as.integer(tick), n_hypnozoites = as.integer(n_hypnozoites))]
  act <- ev[event_type == "hypnozoite_activation",
            .(batch_pedigree_id = as.integer(parent1_pedigree_id), tick = as.integer(tick),
              n = as.integer(n_hypnozoites),
              masked = if ("into_patent_lineage" %in% names(ev)) into_patent_lineage == 1L else NA)]
  if (exclude_masked) {
    gap <- act[is.na(masked), .N]
    if (gap)
      stop("`exclude_masked = TRUE` needs into_patent_lineage on every hypnozoite_activation row, and ", gap,
           " of these events' activation rows lack it.", call. = FALSE)
  }
  ## every row that takes hypnozoites out of a batch, in day order; a masked relapse still does
  out <- rbindlist(list(
    act[, .(batch_pedigree_id, tick, n, relapse = !(exclude_masked & masked))],
    ev[event_type %in% .hypnozoite_exit_events(),
       .(batch_pedigree_id = as.integer(pedigree_id), tick = as.integer(tick), n = as.integer(n_hypnozoites),
         relapse = FALSE)]))
  out <- b[, .(batch_pedigree_id, n_hypnozoites)][out, on = "batch_pedigree_id", nomatch = NULL]
  setorder(out, batch_pedigree_id, tick)
  out[, left := n_hypnozoites - cumsum(n), by = batch_pedigree_id]
  first <- unique(out[relapse == TRUE, .(batch_pedigree_id, t_relapse = tick)], by = "batch_pedigree_id")
  empty <- unique(out[left <= 0L, .(batch_pedigree_id, t_empty = tick)], by = "batch_pedigree_id")
  m <- empty[first[b, on = "batch_pedigree_id"], on = "batch_pedigree_id"]
  m[!is.na(t_relapse) & !is.na(t_empty) & t_empty < t_relapse, t_relapse := NA_integer_]
  m[, end_reason := fcase(!is.na(t_relapse), "relapse", !is.na(t_empty), "emptied", default = "run_end")]
  m[, event := as.integer(end_reason == "relapse")]
  m[, time := fcase(end_reason == "relapse", t_relapse - batch_tick,
                    end_reason == "emptied" & emptied == "censor", t_empty - batch_tick,
                    default = NA_integer_)]
  m[is.na(time), time := end - batch_tick]
  m <- m[, .(batch_pedigree_id, human_id, lineage_id, imported, batch_tick, n_hypnozoites, time, event, end_reason)]
  m[order(batch_tick, batch_pedigree_id)]
}

#' Host-state prevalence and the hypnozoite reservoir over time
#'
#' The share of hosts in each of five states on each day of the state table: no parasite; patent
#' blood-stage clones, without or with hypnozoites; subpatent clones and no patent one, with or
#' without hypnozoites; and hypnozoites alone, a reservoir no blood survey can see. Each day also
#' carries its mean number of hypnozoites per host.
#' @param state a state table from [load_state()].
#' @return a data.table, one row per day of the state table and state, by day and state: `tick`,
#'   `cls` (a factor with levels `"uninfected"`, `"active only"`, `"active + reservoir"`,
#'   `"subpatent"` and `"reservoir only (invisible)"`), `n` (hosts, 0 when none), `frac` (`n`
#'   over the day's hosts) and `hypnozoites_per_host` (the day's mean `hyp_reservoir_size` over its
#'   hosts, the same on each of its rows, so `unique(x[, .(tick, hypnozoites_per_host)])` is the
#'   reservoir's trajectory). A state table without `moi_subpatent` counts no subpatent host.
#' @seealso [plot_host_state_prevalence()]
#' @export
host_state_prevalence <- function(state) {
  lv <- names(.HOST_STATE_FILL)
  st <- as.data.table(state)[entity_type == "human"]
  subp <- if ("moi_subpatent" %in% names(st)) fcoalesce(as.integer(st$moi_subpatent), 0L) else integer(nrow(st))
  st[, cls := fcase(moi_active > 0 & moi_dormant > 0, lv[3], moi_active > 0, lv[2],
                    subp > 0, lv[4], moi_dormant > 0, lv[5], default = lv[1])]
  cnt <- st[, .(n = .N), by = .(tick = as.integer(tick), cls)]
  out <- cnt[CJ(tick = unique(cnt$tick), cls = lv), on = .(tick, cls)][is.na(n), n := 0L]
  out[, frac := n / sum(n), by = tick]
  res <- st[, .(hypnozoites_per_host = mean(hyp_reservoir_size)), by = .(tick = as.integer(tick))]
  out <- res[out, on = "tick"]
  out[, cls := factor(cls, levels = lv)]
  out[order(tick, cls), .(tick, cls, n, frac, hypnozoites_per_host)]
}

#' A run's epidemiology over one window of days
#'
#' One row that summarises a run over a window of its days, by default the second half. Every
#' column covers the same days, so the rows of different runs, or of a grid, mean the same thing.
#'
#' **The window.** The first `burn_frac` of the run's `SIM_DAYS` days are burn-in, so the window is
#' days `round(SIM_DAYS * burn_frac) + 1` to `SIM_DAYS`; `from` and `to` set it directly. A run that
#' went extinct, or never started, is summarised over the same days, the days after its last state
#' row counting as days with no parasite. A run that stopped early without either record is refused
#' unless `to` is at most its last recorded day.
#'
#' **What is counted.** Cases are the clone onsets in the window, as [classify_cases()] gives them:
#' a bite that delivers two genotypes is two new clones, though one infectious bite, and founder
#' seeding is not an onset. A recrudescence is an episode but not an acquisition, so it is in
#' `molfob_episodes_per_py` and `relapse_frac_episodes` only. Masked relapses, into a patent clone
#' of their own lineage, are counted unless `exclude_masked = TRUE` leaves them out of every relapse
#' column; `n_relapse_masked` counts them either way. `NU_HYPNOZOITE_BATCH` was fitted with every
#' relapse counted (see [relapse_fraction()]).
#'
#' **Person-years.** The EIR, FOI and molFOB divide by the person-years of every host, `N_HOSTS`
#' times the window's days over 365. The relapse rate divides by the person-years at risk: each
#' host counts from the day after its first infection by a bite or an import, a blood-stage primary
#' or a hypnozoite batch, to the window's end; founder seeding does not put a host at risk.
#'
#' **Pooling.** To pool runs, sum the counts and the person-years and divide, or weight each rate by
#' its person-years: an average of rates weights a short or small run as much as a long one.
#'
#' @param run a `vivax_run` or its directory.
#' @param burn_frac the share of `SIM_DAYS` left out as burn-in, at least 0 and below 1: the window
#'   starts on day `round(SIM_DAYS * burn_frac) + 1`. Not with `from`.
#' @param from,to the window's first and last days, both included: `from = 0` is day 1, and
#'   `to = NULL` is `SIM_DAYS`.
#' @param exclude_masked if `TRUE`, leave masked relapses out of every relapse column.
#' @return a one-row data.table, nothing rounded:
#'   * the window and the run: `from_tick`, `to_tick` and `n_days` (integer), `exclude_masked`,
#'     `extinct` and `never_started` (logical), and `extinct_tick` (integer, `NA` for a run that
#'     did not go extinct);
#'   * host state, over the window's host-days: `prevalence_frac` (the share with a patent clone),
#'     `coi_mean` and `coi_median` (patent lineages per infected host-day) and `polyclonal_frac` (the
#'     share of infected host-days with two or more), each `NA` with no infected host-day, and
#'     `hypnozoites_per_host`;
#'   * rates per person-year: `eir_per_py` (infectious bites; `NA` for a state table without
#'     `infectious_bites`), `foi_per_py` (new clones from bites and imports),
#'     `molfob_clones_per_py` (and from relapses), `molfob_episodes_per_py` (and from
#'     recrudescences) and `relapses_per_py_at_risk` (`NA` with no time at risk);
#'   * fractions: `relapse_frac_episodes` and `relapse_frac_acquisitions`, as
#'     [relapse_fraction()] gives them, and `imported_frac` (imports among the new clones from bites
#'     and imports), each `NA` with nothing to divide;
#'   * counts (integer): `n_imported`, `n_local_reinfection`, `n_relapse`, `n_relapse_masked` (`NA`
#'     for an event log without the flag), `n_recrudescence` and `n_bloodstage_courses` (blood-stage
#'     courses, by case management or MDA);
#'   * denominators: `person_years` and `person_years_at_risk`.
#' @seealso [force_of_infection()], [relapse_fraction()], [classify_cases()], [plot_diagnostics()]
#' @examples
#' \dontrun{
#' summarise_run(run)                        # the second half of the run
#' summarise_run(run, from = 366, to = 730)  # its second year
#' }
#' @export
summarise_run <- function(run, burn_frac = 0.5, from = NULL, to = NULL, exclude_masked = FALSE) {
  .check_flag(exclude_masked, "exclude_masked")
  run <- load_run(run)
  w <- .run_window(.sim_days(run, "summarise_run"), burn_frac, from, to,
                   both = !missing(burn_frac) && !is.null(from))
  .summary_row(run, load_events(run), load_state(run), w, exclude_masked)
}

## The days a summary covers (D69), both included: from day round(SIM_DAYS * burn_frac) + 1, or
## `from` (0 is day 1), to `to`, SIM_DAYS by default.
.run_window <- function(sim_days, burn_frac = 0.5, from = NULL, to = NULL, both = FALSE) {
  if (both) stop("give `burn_frac` or `from`, not both.", call. = FALSE)
  if (is.null(from)) {
    if (!is.numeric(burn_frac) || length(burn_frac) != 1L || is.na(burn_frac) || burn_frac < 0 ||
        burn_frac >= 1)
      stop("`burn_frac` must be a number at least 0 and below 1: the share of the run's ", sim_days,
           " days (SIM_DAYS) left out as burn-in.", call. = FALSE)
    from <- round(sim_days * burn_frac) + 1
  } else from <- .check_day(from, "from")
  to <- if (is.null(to)) sim_days else .check_day(to, "to")
  if (to > sim_days)
    stop("`to` = ", to, " is after the run's last day, SIM_DAYS = ", sim_days, ".", call. = FALSE)
  from <- max(from, 1)
  if (from > to)
    stop("the window has no day: it starts on day ", from, " and ends on day ", to, ".", call. = FALSE)
  list(from = as.integer(from), to = as.integer(to), n_days = as.integer(to - from + 1))
}

## summarise_run()'s row for window `w`, from tables already read; plot_diagnostics() shares it.
.summary_row <- function(run, events, hosts, w, exclude_masked) {
  nh <- as.integer(.params_of(run)$N_HOSTS); sim_days <- .sim_days(run)
  ## the days after the last state row are days with no parasite only once the run has ended:
  ## extinct, or never started
  last <- if (nrow(hosts)) max(hosts$tick) else 0L
  ended <- if (!is.na(run$extinct_tick)) run$extinct_tick else if (isTRUE(run$never_started)) last else sim_days
  if (w$to > last && last < ended)
    stop(structure(class = c("malariafwd_unfinished_run", "error", "condition"), list(message = paste0(
      "state.tsv in ", .run_dir(run), " stops on day ", last, " of ", sim_days, ", and params.tsv ",
      "records no extinction and no start failure, so the run did not finish: summarise the days it ",
      "recorded with `to = ", last, "`."), call = NULL)))
  hw <- hosts[tick >= w$from & tick <= w$to]
  if (nrow(hw) != as.numeric(nh) * max(0L, min(w$to, last) - w$from + 1L))
    stop("state.tsv in ", .run_dir(run), " does not hold one row per host (N_HOSTS = ", nh,
         ") for each day from ", w$from, " to ", min(w$to, last), ".", call. = FALSE)
  host_days <- as.numeric(nh) * w$n_days; py <- host_days / 365
  inf <- hw$moi_active[hw$moi_active > 0]
  evw <- events[tick >= w$from & tick <= w$to]
  cs_all <- classify_cases(evw)
  cs <- if (exclude_masked) classify_cases(evw, exclude_masked = TRUE) else cs_all
  rf <- relapse_fraction(evw, exclude_masked)
  n_of <- function(o) cs[origin == o, .N]
  n_imp <- n_of("imported"); n_loc <- n_of("local reinfection")
  n_rel <- n_of("relapse"); n_rec <- n_of("recrudescence")
  flag <- cs_all[origin == "relapse"]$into_patent_lineage
  ## the whole log: a host's first establishment can lie in the burn-in
  risk <- .relapses_per_person_year(events, sim_days, w$from, w$to, exclude_masked)
  frac <- function(a, b) if (b > 0) a / b else NA_real_
  data.table(
    from_tick = w$from, to_tick = w$to, n_days = w$n_days, exclude_masked = exclude_masked,
    extinct = !is.na(run$extinct_tick), extinct_tick = as.integer(run$extinct_tick),
    never_started = isTRUE(run$never_started),
    prevalence_frac = sum(hw$moi_active > 0) / host_days,
    coi_mean = if (length(inf)) mean(inf) else NA_real_,
    coi_median = if (length(inf)) as.numeric(stats::median(inf)) else NA_real_,
    polyclonal_frac = if (length(inf)) mean(inf > 1) else NA_real_,
    hypnozoites_per_host = sum(as.numeric(hw$hyp_reservoir_size)) / host_days,
    eir_per_py = if (anyNA(hw$infectious_bites)) NA_real_ else sum(as.numeric(hw$infectious_bites)) / py,
    foi_per_py = (n_imp + n_loc) / py,
    molfob_clones_per_py = (n_imp + n_loc + n_rel) / py,
    molfob_episodes_per_py = (n_imp + n_loc + n_rel + n_rec) / py,
    relapses_per_py_at_risk = risk$relapses_per_py,
    relapse_frac_episodes = rf$relapse_frac_episodes,
    relapse_frac_acquisitions = rf$relapse_frac_acquisitions,
    imported_frac = frac(n_imp, n_imp + n_loc),
    n_imported = n_imp, n_local_reinfection = n_loc, n_relapse = n_rel,
    n_relapse_masked = if (anyNA(flag)) NA_integer_ else sum(flag),
    n_recrudescence = n_rec,
    ## courses given, not clones cleared: one course clears every clone in a host, and a clone can
    ## survive it (bloodcure_persist); an MDA round clears through mda_bloodstage_clear
    n_bloodstage_courses = evw[event_type %in% c("bloodcure_given", "mda_bloodcure_given"), .N],
    person_years = py, person_years_at_risk = risk$person_years)
}

#' Relapse-vs-origin IBD: each relapse paired to a same-lineage prior primary
#' @param events an events data.table.
#' @param ibd a `vivax_ibd` object from [ibd_pairs()] holding the relapses and their origins.
#' @return a data.table with an `ibd` column (should be ~1, clonal).
#' @export
ibd_relapse_origin <- function(events, ibd) {
  ev <- as.data.table(events); x <- ibd; nm <- x$samples
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
  pr[, ibd := ibd_lookup(x, rp, pp)][]
}

#' New-infection genetic diversity (relatedness of the transmitting population)
#'
#' Pairwise IBD among non-relapse blood-stage infections in the steady-state window;
#' clonal fraction rises / lineage richness falls as transmission falls (the surveillance
#' signal). Requires the tree sequence (IBD).
#' @param run a `vivax_run` or run directory.
#' @param burn_frac the fraction of `SIM_DAYS` discarded as burn-in.
#' @param cap max cases sampled.
#' @param clonal IBD threshold for a "clonal" pair.
#' @param seed seed for the subsample when the pool exceeds `cap`, so repeat calls on the same
#'   run agree. The caller's RNG stream is left untouched.
#' @return a one-row data.table of `n`, `mean_ibd`, `clonal_frac`, `n_lineages`.
#' @inheritSection tskit_available Python
#' @export
new_infection_diversity <- function(run, burn_frac = 0.5, cap = 250, clonal = 0.5, seed = 1) {
  ev <- load_events(run)
  ## the window summarise_run() gives the same burn_frac: from the day after the burn-in
  from <- .run_window(.sim_days(run, "new_infection_diversity"), burn_frac)$from
  ## Every non-relapse acquisition route, derived from classify_cases() rather than from a
  ## literal event type. Matching "vh_inoculation" alone dropped imported genomes, which
  ## overstated mean IBD and clonal fraction and understated lineage richness by 44% at a 17%
  ## import share and 4.3x at 61% -- all three biases pointing the same way as a genuine
  ## transmission decline, i.e. in the direction of the signal this metric is read for.
  ## acquisitions only. `origin != "relapse"` would now admit recrudescences, which are not
  ## newly acquired genomes - they were already in the host - and would bias this toward the
  ## clonality it exists to measure.
  peds <- classify_cases(ev)[origin %in% c("imported", "local reinfection") & tick >= from,
                             unique(pedigree_id)]
  if (length(peds) < 6)
    return(data.table(n = length(peds), mean_ibd = NA_real_, clonal_frac = NA_real_, n_lineages = NA_integer_))
  peds <- .sample_fixed(peds, cap, seed)
  .require_python("new_infection_diversity()", c("tskit", "numpy"))
  x <- suppressWarnings(.ibd_pairs(run, ids = peds, none_ok = TRUE))
  vp <- x$samples$pedigree_id
  if (length(vp) < 2)
    return(data.table(n = length(vp), mean_ibd = NA_real_, clonal_frac = NA_real_, n_lineages = NA_integer_))
  di <- x$pairs$ibd
  data.table(n = length(vp), mean_ibd = mean(di), clonal_frac = 100 * mean(di >= clonal),
             n_lineages = uniqueN(ev[pedigree_id %in% vp, lineage_id]))
}
