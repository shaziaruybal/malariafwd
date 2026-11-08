## events.R -- the events.tsv event-type vocabulary + helpers.

# Event types that END a blood-stage episode (natural clearance + drug-induced).
# Use this wherever you reconstruct blood-stage episode spans, so treated / MDA'd
# cases are not silently missed -- they clear via these events, NOT `human_clearance`
# alone. (Reconstructing spans off `human_clearance` only over-counts durations when
# any blood-stage treatment or MDA is switched on.)
## Derived, not restated: duplicating the dictionary here meant the two could disagree, and
## they did - bloodcure_subpatent ends a patent episode and would have been missed.
.bloodstage_clear_events <- function()
  event_types()$event_type[event_types()$ends_episode]

## ...and the same for the other side. Both sides of the episode vocabulary are now derived,
## so adding an event type to the dictionary is the only edit needed.
.onset_events <- function()
  event_types()$event_type[event_types()$is_onset]

## The onsets that bring a parasite into a host: a bite, an import or founder seeding.
.inoculation_events <- function()
  setdiff(.onset_events(), c("hypnozoite_activation", "recrudescence"))

#' The malariafwd event dictionary
#'
#' Every row in a run's `events.tsv` carries an `event_type`. This returns a tidy table
#' of the event types the model emits and what each means, so you can navigate an event
#' log. An inoculation that establishes neither a blood-stage infection nor a hypnozoite
#' batch is logged with a `_no_take` suffix (e.g. `vh_inoculation_no_take`), and so is an
#' import whose host already carries `MAX_MOI` blood-stage lineages
#' (`vh_inoculation_import_no_take`); a local bite refused at `MAX_MOI` leaves no row.
#'
#' **Derive filters from this table rather than restating event names.** `ends_episode`
#' marks every way a patent blood-stage episode can end and `is_onset` every way one can
#' begin, so episode reconstruction stays correct as the vocabulary grows.
#'
#' `mosq_death` is one row per parasite genome an infected mosquito carried when it died; the
#' deaths of uninfected mosquitoes are not logged, so its rows do not count deaths.
#'
#' @return a data.frame with columns `event_type`, `stage`, `ends_episode` (does this event
#'   end a patent blood-stage episode?), `is_onset` (does it begin one?) and `description`.
#' @examples
#' event_types()
#' subset(event_types(), is_onset)$event_type       # every way an episode can begin
#' @export
event_types <- function() {
  ## ONE RECORD PER ROW, deliberately. This used to be four positionally aligned vectors,
  ## and inserting seven rows into them at the wrong offset put another event's text on 9 of
  ## 23 rows of a public dictionary - caught by a regression sweep, not by the test suite.
  ## Here a mis-edit is a visibly short row, not a silent shift.
  r <- function(event_type, stage, ends_episode, is_onset, description)
    data.frame(event_type = event_type, stage = stage, ends_episode = ends_episode,
               is_onset = is_onset, description = description, stringsAsFactors = FALSE)
  rbind(
    r("vh_inoculation",                "blood",  FALSE, TRUE,
      "new local mosquito-borne blood-stage infection"),
    r("vh_inoculation_import",         "blood",  FALSE, TRUE,
      "imported blood-stage infection"),
    r("vh_inoculation_seed",           "blood",  FALSE, TRUE,
      "founder blood-stage infection (burn-in seeding)"),
    r("vh_inoculation_no_take",        "blood",  FALSE, FALSE,
      "a bite established neither a blood-stage infection nor hypnozoites (P_PRIMARY failed and the batch drew 0)"),
    r("vh_inoculation_import_no_take", "blood",  FALSE, FALSE,
      "an importation established neither compartment (P_PRIMARY failed and the batch drew 0), or its host was at MAX_MOI"),
    r("hypnozoite_activation",         "blood",  FALSE, TRUE,
      "hypnozoite activates -> a blood-stage relapse (clonal to its origin)"),
    r("recrudescence",                 "blood",  FALSE, TRUE,
      "a suppressed clone re-expands to patency after a failed course; mints a new genome clonal to the holder (P_RECRUD_GIVEN_FAIL)"),
    r("human_clearance",               "blood",  TRUE,  FALSE,
      "blood-stage infection clears naturally"),
    r("case_bloodstage_clear",         "blood",  TRUE,  FALSE,
      "blood-stage infection cleared by case-based treatment (BLOODCURE_COVERAGE)"),
    r("mda_bloodstage_clear",          "blood",  TRUE,  FALSE,
      "blood-stage infection cleared by mass drug administration (MDA_TICK), also one that became patent earlier that day (the day's imports arrive after the round)"),
    r("mda_subpatent_clear",           "blood",  FALSE, FALSE,
      "a suppressed clone cleared by MDA, so it never recrudesces; its patent episode had already ended at bloodcure_subpatent"),
    r("bloodcure_subpatent",           "blood",  TRUE,  FALSE,
      "a clone survived treatment but was suppressed BELOW PATENCY; the patent episode ends here and the clone recrudesces later (P_RECRUD_GIVEN_FAIL), unless an MDA round clears it first (mda_subpatent_clear)"),
    r("bloodcure_given",               "blood",  FALSE, FALSE,
      "a blood-stage course is given to a host (BLOODCURE_COVERAGE)"),
    r("radcure_given",                 "liver",  FALSE, FALSE,
      "radical cure is given to a host, nested inside a blood-stage course"),
    r("radcure_failure",               "liver",  FALSE, FALSE,
      "radical cure was given but the whole reservoir survived (RADCURE_MODE = person)"),
    r("bloodcure_persist",             "blood",  FALSE, FALSE,
      "a blood-stage clone survived treatment and continues at full density (P_BLOODCURE_SURVIVAL)"),
    r("mda_bloodcure_given",           "blood",  FALSE, FALSE,
      "a blood-stage course is given by MDA to a host not empty at the round (MDA_COVERAGE)"),
    r("mda_radcure_given",             "liver",  FALSE, FALSE,
      "radical cure is given by MDA"),
    r("mda_radcure_failure",           "liver",  FALSE, FALSE,
      "MDA radical cure given but the whole reservoir survived"),
    r("hypnozoite_clearance",          "liver",  FALSE, FALSE,
      "a hypnozoite dies / clears in the liver"),
    r("radcure_hypnozoite_clear",      "liver",  FALSE, FALSE,
      "hypnozoite reservoir thinned by case-based radical cure (RADCURE_COVERAGE)"),
    r("mda_hypnozoite_clear",          "liver",  FALSE, FALSE,
      "hypnozoite reservoir thinned by MDA"),
    r("hv_uptake_clonal",              "vector", FALSE, FALSE,
      "mosquito takes up a clonal (selfed) oocyst genotype"),
    r("hv_uptake_recombinant",         "vector", FALSE, FALSE,
      "mosquito takes up a recombinant (outcrossed) oocyst genotype"),
    r("mosq_death",                    "vector", FALSE, FALSE,
      "a parasite genome lost when an infected mosquito dies: one row per genome it carried"))
}
