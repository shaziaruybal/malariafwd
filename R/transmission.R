## transmission.R -- who infected whom, read from the event log (D70): the links from each parent
## genome of an oocyst to each genome its products infected.

## The link table's columns and classes, which every path returns (D23).
.links0 <- function() data.table(
  source_pedigree_id = integer(), target_pedigree_id = integer(), kind = character(),
  source_human_id = integer(), target_human_id = integer(), mosq_id = integer(), oocyst_id = integer(),
  source_tick = integer(), uptake_tick = integer(), inoculation_tick = integer(), target_tick = integer(),
  weight = double())

#' Who infected whom: the ground-truth transmission links
#'
#' One row per link from a parent genome of an oocyst to a genome the oocyst's products infected. A
#' mosquito takes up gametes of one genome of the host it feeds on, or of two that it outcrosses,
#' as an oocyst, and a later bite delivers the oocyst's products to another host. The link is
#' **direct** when the bite starts a blood-stage infection, and **liver** when the bite leaves a
#' hypnozoite batch: each relapse of the batch is linked to the genomes the mosquito took up.
#'
#' Only local infections are targets. Imports, founders, the relapses of an import's batch and
#' recrudescences never are, though each can be a source. A host's own relapse is never linked to
#' its earlier infection, since no mosquito passes it on, but a mosquito that feeds on a host and
#' later bites the same host makes a link like any other. Each parent of an outcrossed oocyst is
#' linked with `weight` 1/2 and every other link has weight 1, so a target's weights sum to 1;
#' `unique(links, by = "target_pedigree_id")` gives one row per infected genome with the host that
#' infected it. Join [classify_cases()] on `target_pedigree_id = pedigree_id` for a target's
#' detection and whether a relapse was masked.
#'
#' Pass the whole event log, as [load_events()] gives it, and filter the result: a link needs the
#' rows of the feed, the bite and the source's onset, which can lie long before its target.
#'
#' @param events an events table from [load_events()].
#' @param kind `"all"`, `"direct"` or `"liver"`.
#' @return a data.table, one row per (source genome, target genome), by `target_tick`,
#'   `target_pedigree_id` and `source_pedigree_id`: `source_pedigree_id` (the infector genome),
#'   `target_pedigree_id` (the infected genome), `kind` (`"direct"` or `"liver"`),
#'   `source_human_id` (the host the mosquito fed on), `target_human_id` (the host it bit),
#'   `mosq_id`, `oocyst_id`, `source_tick` (the day the source's blood stage began), `uptake_tick`
#'   (the feed), `inoculation_tick` (the bite), `target_tick` (the bite, or for a liver link the
#'   relapse) and `weight` (1, or 1/2 for each parent of an outcrossed oocyst). Ids and days are
#'   integer.
#' @seealso [rc_by_origin()], [transmission_chains()]
#' @examples
#' \dontrun{transmission_links(load_events(run), kind = "liver")}
#' @export
transmission_links <- function(events, kind = c("all", "direct", "liver")) {
  kind <- match.arg(kind)
  if (!is.data.frame(events))
    stop("`events` must be an events table from load_events(); for a run, pass load_events(run).", call. = FALSE)
  need <- c("tick", "event_type", "mosq_id", "human_id", "pedigree_id", "parent1_pedigree_id",
            "parent2_pedigree_id", "is_dormant", "oocyst_id")
  miss <- setdiff(need, names(events))
  if (length(miss))
    stop("`events` lacks the column(s) ", paste(miss, collapse = ", "), ": pass the table load_events() gives.",
         call. = FALSE)
  ev <- as.data.table(events)
  ## one row per oocyst: the host the mosquito fed on, the mosquito, the day and the parent
  ## genome(s). -1 is logEvent()'s default, for a row with no oocyst or no second parent; oocyst
  ## ids start at 1 and pedigree 0 is a real genome, the first founder, so test the sentinel.
  up <- unique(ev[event_type %in% c("hv_uptake_clonal", "hv_uptake_recombinant") & oocyst_id != -1L,
                  .(oocyst_id, source_human_id = human_id, mosq_id, uptake_tick = tick,
                    p1 = parent1_pedigree_id, p2 = parent2_pedigree_id)], by = "oocyst_id")
  tg <- list()
  if (kind != "liver")
    tg$direct <- ev[event_type == "vh_inoculation" & is_dormant == 0L & oocyst_id != -1L,
                    .(oocyst_id, target_pedigree_id = pedigree_id, target_human_id = human_id,
                      inoculation_tick = tick, target_tick = tick, kind = "direct")]
  if (kind != "direct") {
    ## a relapse carries its batch's oocyst, and names the batch as its parent1_pedigree_id
    batch <- ev[event_type == "vh_inoculation" & is_dormant == 1L & oocyst_id != -1L,
                .(batch = pedigree_id, inoculation_tick = tick)]
    liver <- ev[event_type == "hypnozoite_activation" & oocyst_id != -1L,
                .(oocyst_id, target_pedigree_id = pedigree_id, target_human_id = human_id,
                  batch = parent1_pedigree_id, target_tick = tick, kind = "liver")]
    tg$liver <- batch[liver, on = "batch"][, batch := NULL]
  }
  tg <- rbindlist(tg, use.names = TRUE)
  no_feed <- if (nrow(tg)) tg[!oocyst_id %in% up$oocyst_id, uniqueN(target_pedigree_id)] else 0L
  m <- up[tg, on = "oocyst_id", nomatch = NULL]
  e <- rbindlist(list(m[, c(.SD, .(source_pedigree_id = p1))],
                      m[p2 != -1L, c(.SD, .(source_pedigree_id = p2))]), use.names = TRUE)
  no_day <- 0L
  if (nrow(e)) {
    e[, weight := 1 / .N, by = target_pedigree_id]
    onset <- unique(ev[event_type %in% .onset_events() & is_dormant == 0L,
                       .(source_pedigree_id = pedigree_id, source_tick = tick)], by = "source_pedigree_id")
    e <- onset[e, on = "source_pedigree_id"]
    no_day <- e[is.na(inoculation_tick) | is.na(source_tick), .N]
  }
  if (no_feed || no_day)
    warning("transmission_links(): this event log lacks the rows behind some links (", no_feed,
            " infected genomes have no feed, ", no_day, " links lack a bite or onset day): pass the whole ",
            "event log, as load_events() gives it, and filter the result.", call. = FALSE)
  if (!nrow(e)) return(.links0())
  setorder(e, target_tick, target_pedigree_id, source_pedigree_id)
  out <- e[, names(.links0()), with = FALSE]
  ## the template's types, so that pooled tables bind
  for (nm in setdiff(names(out), c("kind", "weight"))) set(out, j = nm, value = as.integer(out[[nm]]))
  out[]
}

## A genome as the clone it re-expanded from, walked back through repeated recrudescences: a
## recrudescence row's parent1_pedigree_id is the clone it came from.
.acquisition_root <- function(ped, events) {
  rec <- events[event_type == "recrudescence",
                .(pedigree_id = as.integer(pedigree_id), parent = as.integer(parent1_pedigree_id))]
  ped <- as.integer(ped)
  for (k in seq_len(nrow(rec))) {
    i <- match(ped, rec$pedigree_id)
    if (all(is.na(i))) break
    ped[!is.na(i)] <- rec$parent[i[!is.na(i)]]
  }
  ped
}

#' Case reproduction number (Rc) by true origin
#'
#' The mean number of secondary cases each case caused, by the case's true origin. A case is a
#' blood-stage clone onset, a row of [classify_cases()]; its secondary cases are the targets of its
#' links in [transmission_links()], infected by a bite's blood stage or by a relapse of the batch
#' the bite left.
#'
#' @details
#' **Cases.** With `unit = "episode"` every clone onset is a case, a recrudescence included. With
#' `unit = "acquisition"` a recrudescence is not a case: its onward transmission is credited to the
#' clone it re-expanded from, through repeated recrudescences, so each acquisition has one count. A
#' subpatent holder's transmission is credited to the case that began it. Founders infect, but are
#' not cases.
#'
#' **Secondary cases.** Each secondary case counts once: one from an outcrossed oocyst credits 1/2
#' to each parent. A bite that delivers two genotypes gives two secondary cases. Masked relapses
#' count like any other; join [classify_cases()] on the links to leave them out. `kind = "direct"`
#' and `kind = "liver"` add up to `kind = "all"`.
#'
#' **The end of the run.** A case near the end has had less time to transmit, above all through
#' relapses: a hypnozoite stays dormant `1 / (ALPHA_ACTIVATE + MU_HYPNOZOITE)` days on average,
#' about four months at the defaults, and more than a year for 1 in 20. `followup` counts only the
#' cases that began at least that many days before `SIM_DAYS`. A run that went extinct is not
#' censored.
#'
#' **The phase of the epidemic.** Rc reflects the phase the counted cases began in: above 1 while
#' an epidemic grows from its founders, near 1 at a self-sustaining steady state, and lower where
#' imports sustain transmission. To compare runs, choose `from` and `followup` inside a stretch that
#' [plot_diagnostics()] shows to be steady; `followup` alone keeps a growing run's early cases.
#'
#' @param run a `vivax_run` or its directory.
#' @param unit `"episode"` or `"acquisition"`: what a case is (see Details).
#' @param kind `"all"`, `"direct"` or `"liver"`: the links counted, as in [transmission_links()].
#' @param from count the cases that began on day `from` or later.
#' @param followup count the cases that began on day `SIM_DAYS - followup` or earlier.
#' @return a data.table, one row per origin, every origin present: `origin` (a factor of the
#'   origins [classify_cases()] gives, without `"recrudescence"` for `unit = "acquisition"`),
#'   `n_cases` (integer), `n_secondary` (the sum of the cases' link weights), `Rc`
#'   (`n_secondary / n_cases`) and `onward_frac` (the share of cases with a secondary case), the
#'   last two `NA` with no case. To pool runs, divide the summed `n_secondary` by the summed
#'   `n_cases`.
#' @references Fraser C (2007) Estimating individual and household reproduction numbers in an
#'   emerging epidemic. PLoS ONE 2(8): e758. \doi{10.1371/journal.pone.0000758}
#' @seealso [transmission_links()], [classify_cases()]
#' @examples
#' \dontrun{rc_by_origin(run, from = 365, followup = 365)}
#' @export
rc_by_origin <- function(run, unit = c("episode", "acquisition"), kind = c("all", "direct", "liver"),
                         from = 0, followup = 0) {
  unit <- match.arg(unit); kind <- match.arg(kind)
  from <- .check_day(from, "from"); followup <- .check_day(followup, "followup")
  ev <- load_events(run)
  last <- if (from > 0L || followup > 0L) .sim_days(run, "rc_by_origin") - followup else Inf
  if (max(from, 1L) > last)
    stop("no case can be counted: onsets from day ", max(from, 1L), " to day ", last,
         " (SIM_DAYS - followup) leave none.", call. = FALSE)
  links <- transmission_links(ev, kind)
  cases <- classify_cases(ev)[, .(pedigree_id, origin, tick)]
  lv <- if (unit == "episode") .ORIGINS else .ACQUISITIONS
  if (unit == "acquisition") {
    links[, source_pedigree_id := .acquisition_root(source_pedigree_id, ev)]
    cases <- cases[origin != "recrudescence"]
  }
  cases <- cases[tick >= from & tick <= last]
  sec <- links[, .(n_secondary = sum(weight)), by = .(pedigree_id = source_pedigree_id)]
  rc <- sec[cases, on = "pedigree_id"][is.na(n_secondary), n_secondary := 0]
  out <- rc[, .(n_cases = .N, n_secondary = sum(n_secondary), n_onward = sum(n_secondary > 0)), by = origin]
  out <- out[data.table(origin = lv), on = "origin"]
  out[is.na(n_cases), `:=`(n_cases = 0L, n_secondary = 0, n_onward = 0L)]
  out[, `:=`(Rc = fifelse(n_cases > 0L, n_secondary / n_cases, NA_real_),
             onward_frac = fifelse(n_cases > 0L, n_onward / n_cases, NA_real_))]
  out[, origin := factor(origin, levels = lv)]
  out[, .(origin, n_cases, n_secondary, Rc, onward_frac)][]
}
