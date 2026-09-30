## events.R -- the events.tsv event-type vocabulary + helpers.

# Event types that END a blood-stage episode (natural clearance + drug-induced).
# Use this wherever you reconstruct blood-stage episode spans, so treated / MDA'd
# cases are not silently missed -- they clear via these events, NOT `human_clearance`
# alone. (Reconstructing spans off `human_clearance` only over-counts durations when
# any blood-stage treatment or MDA is switched on.)
.bloodstage_clear_events <- c("human_clearance", "mda_bloodstage_clear", "case_bloodstage_clear")

#' The malariafwd event dictionary
#'
#' Every row in a run's `events.tsv` carries an `event_type`. This returns a tidy table
#' of the event types the model emits and what each means, so you can navigate an event
#' log. An inoculation that establishes neither a blood-stage infection nor a hypnozoite
#' batch is logged with a `_no_take` suffix (e.g. `vh_inoculation_no_take`).
#'
#' The three blood-stage clearance events (`human_clearance`, `mda_bloodstage_clear`,
#' `case_bloodstage_clear`) all end an active episode; treat them together when
#' reconstructing episode durations.
#'
#' @return a data.frame with columns `event_type`, `stage`, `ends_episode` (does this event
#'   end an active blood-stage episode?) and `description`.
#' @examples
#' event_types()
#' @export
event_types <- function() {
  data.frame(
    event_type = c(
      "vh_inoculation", "vh_inoculation_import", "vh_inoculation_seed",
      "hypnozoite_activation",
      "human_clearance", "case_bloodstage_clear", "mda_bloodstage_clear",
      "hypnozoite_clearance", "radcure_hypnozoite_clear",
      "mda_hypnozoite_clear", "mda_hypnozoite_survived",
      "hv_uptake_clonal", "hv_uptake_recombinant", "hv_uptake_seed",
      "mosq_death"),
    stage = c(
      "blood", "blood", "blood",
      "blood",
      "blood", "blood", "blood",
      "liver", "liver",
      "liver", "liver",
      "vector", "vector", "vector",
      "vector"),
    ## the three ways an active blood-stage episode can end; use these together when
    ## reconstructing episode durations, or treated / MDA'd cases are silently missed
    ends_episode = c(
      FALSE, FALSE, FALSE,
      FALSE,
      TRUE, TRUE, TRUE,
      FALSE, FALSE,
      FALSE, FALSE,
      FALSE, FALSE, FALSE,
      FALSE),
    description = c(
      "new local mosquito-borne blood-stage infection",
      "imported blood-stage infection",
      "founder blood-stage infection (burn-in seeding)",
      "hypnozoite activates -> a blood-stage relapse (clonal to its origin)",
      "blood-stage infection clears naturally",
      "blood-stage infection cleared by case-based treatment (BLOODCURE_COVERAGE)",
      "blood-stage infection cleared by mass drug administration (MDA_TICK)",
      "a hypnozoite dies / clears in the liver",
      "hypnozoite reservoir thinned by case-based radical cure (RADCURE_COVERAGE)",
      "hypnozoite reservoir thinned by MDA",
      "a hypnozoite survives MDA (logged for accounting)",
      "mosquito takes up a clonal (selfed) oocyst genotype",
      "mosquito takes up a recombinant (outcrossed) oocyst genotype",
      "mosquito takes up a founder (seed) genotype",
      "mosquito dies"),
    stringsAsFactors = FALSE)
}
