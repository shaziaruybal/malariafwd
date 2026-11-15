## episodes.R -- blood-stage clones and their spans, and each episode against the earlier ones of
## its host (D46). An episode is the clones that became patent in one host on one day, named by
## episode_id "<human_id>:<tick>".

#' Every blood-stage clone of a run, from onset to end
#'
#' One row per blood-stage clone, founder seeding included: when it became patent, how, and when
#' and how it ended. A clone counts in `state.tsv`'s `moi_active` from `onset_tick` to
#' `end_tick - 1`, so a clone that ends on its first day (an MDA round can clear one) is never in
#' `state.tsv`. A clone that ends below patency (`bloodcure_subpatent`) returns, if it does, as a
#' recrudescence: a row of its own.
#' @param events an events table from [load_events()], or a run (a `vivax_run` or its directory).
#' @return a data.table, by host, onset and pedigree id: `episode_id`, `human_id`, `pedigree_id`,
#'   `lineage_id`, `onset_tick`, `onset_event` (`vh_inoculation`, `vh_inoculation_import`,
#'   `vh_inoculation_seed`, `hypnozoite_activation` or `recrudescence`), `end_tick` (`NA` for a
#'   clone still patent when the run ended) and `end_reason` (`human_clearance`,
#'   `case_bloodstage_clear`, `mda_bloodstage_clear`, `bloodcure_subpatent`, or `"ongoing"`). Join
#'   [classify_cases()] on `pedigree_id` for each clone's origin; founders have none.
#' @seealso [recurrences()], [classify_cases()]
#' @examples
#' ep <- episodes(example_run())
#' head(ep)
#' @export
episodes <- function(events) {
  ev <- .as_events(events)
  on <- ev[event_type %in% .onset_events() & is_dormant == 0,
           .(human_id = as.integer(human_id), pedigree_id = as.integer(pedigree_id),
             lineage_id = as.character(lineage_id), onset_tick = as.integer(tick),
             onset_event = as.character(event_type))]
  en <- ev[event_type %in% .bloodstage_clear_events(),
           .(pedigree_id = as.integer(pedigree_id), end_tick = as.integer(tick), end_reason = as.character(event_type))]
  en <- en[order(pedigree_id, end_tick)][, .SD[1L], by = pedigree_id]   # a clone's first end
  out <- en[on, on = "pedigree_id"]
  out[is.na(end_reason), end_reason := "ongoing"]
  out[, episode_id := paste0(human_id, ":", onset_tick)]
  out <- out[, .(episode_id, human_id, pedigree_id, lineage_id, onset_tick, onset_event, end_tick, end_reason)]
  out[order(human_id, onset_tick, pedigree_id)]
}

## IBD of genome pairs from the run's vivax_ibd object, NA where either genome is not in it.
.ibd_values <- function(ibd, a, b) {
  if (!inherits(ibd, "vivax_ibd"))
    stop("`ibd` must be the run's IBD object, from ibd_pairs().", call. = FALSE)
  have <- ibd$samples$pedigree_id
  v <- rep(NA_real_, length(a)); k <- a %in% have & b %in% have
  if (any(k)) v[k] <- ibd_lookup(ibd, a[k], b[k])
  v
}

#' Each recurrence, against the earlier episodes of its host
#'
#' A recurrence is every episode of a host after its first: a relapse, a reinfection, an import
#' or a recrudescence. Each of its clones is paired with each clone of the earlier episodes (`prior
#' = "all"`) or of the latest one (`"previous"`); the clones of one episode are never paired, and
#' founder seeding is not an episode. Whether the earlier clone was still patent is not a condition:
#' join [episodes()] on `prior_pedigree_id` for that. Give the run's IBD object to read each pair's
#' IBD.
#'
#' With `prior = "all"` the table grows with the square of a host's episodes. To pair each
#' recurrence with the latest earlier episode on the observed line list, filter
#' `recurrences(ev)` to `prior_detected == 1` and keep each clone's row with the largest
#' `prior_tick`.
#' @param events an events table from [load_events()], or a run (a `vivax_run` or its directory).
#' @param prior `"all"` (every earlier episode) or `"previous"` (the latest one).
#' @param from keep the recurrences that begin on day `from` or later, paired with earlier
#'   episodes whenever they were.
#' @param ibd `NULL`, or the run's IBD object from [ibd_pairs()], for an `ibd` column.
#' @param exclude_masked if `TRUE`, masked relapses (into a patent clone of their own lineage, see
#'   [classify_cases()]) are neither recurrences nor earlier episodes.
#' @return a data.table, one row per (recurrence clone, earlier clone), by host, day and pedigree
#'   id: `human_id`, `episode_id`, `tick`, `pedigree_id`, `lineage_id`, `origin`, `imported`,
#'   `detected`, `into_patent_lineage`, the same of the earlier clone prefixed `prior_`
#'   (`prior_episode_id`, `prior_tick`, `prior_pedigree_id`, `prior_lineage_id`, `prior_origin`,
#'   `prior_imported`, `prior_detected`), `days_since` (1 or more), and with `ibd` given, `ibd`
#'   (`NA` where either genome is not in the object).
#' @seealso [episodes()], [classify_cases()], [ibd_pairs()]
#' @examples
#' head(recurrences(example_run()))
#' @examplesIf tskit_available()
#' run <- example_run()
#' x <- ibd_pairs(run, ids = classify_cases(run)$pedigree_id)
#' head(recurrences(run, ibd = x))   # with each pair's IBD
#' @export
recurrences <- function(events, prior = c("all", "previous"), from = 0, ibd = NULL, exclude_masked = FALSE) {
  prior <- match.arg(prior)
  from <- .check_day(from, "from")
  .check_flag(exclude_masked, "exclude_masked")
  if (!is.null(ibd) && !inherits(ibd, "vivax_ibd"))
    stop("`ibd` must be the run's IBD object, from ibd_pairs().", call. = FALSE)
  cs <- classify_cases(events, exclude_masked = exclude_masked)
  ep <- unique(cs[, .(human_id, tick)])[order(human_id, tick)]
  ep[, prev_tick := shift(tick), by = human_id]
  rec <- cs[ep[!is.na(prev_tick) & tick >= from], on = .(human_id, tick), nomatch = NULL]
  pr <- cs[, .(human_id, prior_episode_id = episode_id, prior_tick = tick, prior_pedigree_id = pedigree_id,
               prior_lineage_id = lineage_id, prior_origin = origin, prior_imported = imported,
               prior_detected = detected)]
  out <- if (prior == "previous")
    pr[rec, on = .(human_id, prior_tick = prev_tick), nomatch = NULL, allow.cartesian = TRUE]
  else
    pr[rec, on = .(human_id, prior_tick < tick), nomatch = NULL, allow.cartesian = TRUE,
       .(human_id, prior_episode_id, prior_tick = x.prior_tick, prior_pedigree_id, prior_lineage_id,
         prior_origin, prior_imported, prior_detected, tick = i.tick, episode_id, pedigree_id, lineage_id,
         origin, imported, detected, into_patent_lineage)]
  out <- out[, .(human_id, episode_id, tick, pedigree_id, lineage_id, origin, imported, detected,
                 into_patent_lineage, prior_episode_id, prior_tick, prior_pedigree_id, prior_lineage_id,
                 prior_origin, prior_imported, prior_detected, days_since = as.integer(tick - prior_tick))]
  data.table::setorder(out, human_id, tick, pedigree_id, prior_tick, prior_pedigree_id)
  if (!is.null(ibd)) out[, ibd := .ibd_values(ibd, pedigree_id, prior_pedigree_id)]
  out[]
}

## How a clone began, as the trajectory plots mark it.
.onset_kind <- function(onset_event)
  data.table::fcase(onset_event == "hypnozoite_activation", "relapse",
                    onset_event == "recrudescence", "recrudescence", default = "new infection")
