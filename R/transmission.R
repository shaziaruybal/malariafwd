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
