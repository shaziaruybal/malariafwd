## param-dictionary.R -- one record per model parameter, from which the defaults, ?vivax_params,
## param_info() and the tables of MODEL_DESCRIPTION section 10 are made (D60). The records keep the
## order of the parameters, which is the order of names(vivax_params()) and of SLiM's arguments;
## the documentation reads them by group.
##
## A record: the typed default; its unit and meaning (ASCII); the single-value bounds, which
## validate_params() enforces and a test holds it to (min/max, open or closed; `values`, the
## allowed strings); `rule`, a joint or conditional rule in words; `shown`, the default as
## documented when it is not the plain value; and `source`, in prose.

.param <- function(default, unit, meaning, min = NA, max = NA, min_open = FALSE, max_open = FALSE,
                   values = NULL, rule = NA_character_, shown = NA_character_, source = NA_character_)
  list(default = default, unit = unit, meaning = meaning, min = min, max = max, min_open = min_open,
       max_open = max_open, values = values, rule = rule, shown = shown, source = source)

.param_dict <- list(
  SIM_ID = .param("default", "-", "A label for the run, written into params.tsv.",
    rule = "UTF-8, without control characters"),
  SIM_DAYS = .param(365L, "days",
    "Length of the run; one tick is one day. A run that goes extinct, or never starts, ends sooner.",
    min = 2, max = 1e9,
    rule = "a run that reaches day 99999 must move the days that default to 99999 past its end"),
  VERBOSE = .param(FALSE, "-", "Narrate every event in run.log. The outputs are the same either way."),
  RETAIN_FULL_HISTORY = .param(TRUE, "-",
    paste("Keep every genome in the tree sequence, which IBD to infections that have already cleared",
          "needs. With FALSE most earlier genomes leave the tree, and ibd_pairs() warns that the ids it",
          "was asked about are absent; they look up as NA.")),

  N_HOSTS = .param(40L, "hosts", "Humans; the number is constant.", min = 1),
  N_FOUNDER_HOSTS = .param(12L, "hosts",
    "Hosts infected at the start: establishment, whether the epidemic sustains itself.", min = 0,
    rule = "at most N_HOSTS; 0 needs an import path inside the run"),
  N_FOUNDER_LINEAGES = .param(12L, "genomes",
    paste("Distinct founding genomes among the founder hosts; the others carry clonal copies, seeded on",
          "day 2. The genealogical diversity knob: fewer lineages give higher background relatedness at",
          "the same establishment."),
    min = 0, shown = "N_FOUNDER_HOSTS",
    rule = "at most N_FOUNDER_HOSTS, and at least 1 when there are founder hosts"),
  N_MOSQ = .param(100L, "mosquitoes", "Mosquitoes; the number is constant, and every one starts uninfected.",
    min = 1),

  BITING_RATE = .param(0.1, "per mosquito per day",
    paste("Daily probability that a mosquito bites. Bites per human per day are BITING_RATE x N_MOSQ /",
          "N_HOSTS, not this value."),
    min = 0, min_open = TRUE,
    rule = "a mosquito bites at most once a day, so above 1 with its seasonal and resurgence multipliers it warns"),
  P_TRANS_VH = .param(0.5, "probability", "b: an infectious bite infects the human.", min = 0, max = 1),
  P_TRANS_HV = .param(0.5, "probability", "c: a mosquito that bites an infectious human takes up parasites.",
    min = 0, max = 1),
  MOSQ_LIFESPAN = .param(10, "days", "Mean mosquito lifespan: a mosquito dies each day with probability 1/MOSQ_LIFESPAN.",
    min = 1),
  EIP_DAYS = .param(10L, "days",
    "Extrinsic incubation period: days from uptake until a mosquito's products can be transmitted (fixed).",
    min = 0),
  MAX_MOI = .param(20L, "lineages",
    paste("A host already carrying this many blood-stage lineages, patent or subpatent, takes no new",
          "inoculation, counting the clones it gained earlier that day. Relapses bypass it, and a bite",
          "that delivers several genotypes is admitted whole. 1 means no superinfection; at 20 it is a",
          "guard that seldom binds."),
    min = 1),
  MAX_BROODS_PER_MOSQ = .param(20L, "lineages",
    paste("A mosquito carrying this many lineages takes up no more. One uptake can add up to four per",
          "oocyst, so the count can pass it."),
    min = 0),

  MEAN_BLOOD_STAGE_DURATION = .param(60, "days",
    "Mean duration of a blood-stage clone: it clears each day with probability 1/MEAN_BLOOD_STAGE_DURATION.",
    min = 1, source = "White et al. (2014)"),
  LAMBDA_OOCYSTS = .param(0, "oocysts",
    paste("Pre-truncation rate of the zero-truncated number of oocysts per successful uptake, each its",
          "own meiosis; 0 gives exactly one. See uptake_moments() and uptake_calibrate()."),
    min = 0, max = 1000, rule = "with K_OOCYSTS, the mean count per uptake is at most 1000",
    source = "membrane feeds of wild An. arabiensis in Ethiopia give a median of 7.2 oocysts per midgut (Ayele et al. 2023)"),
  LAMBDA_SPOROZOITES = .param(0, "products per bite",
    paste("Pre-truncation rate of the zero-truncated number of a mosquito's products an infectious bite",
          "samples, repeats of a lineage counted once. 0 gives one genotype per bite, so co-transmission",
          "is impossible at the default."),
    min = 0, max = 1000),
  K_OOCYSTS = .param(0, "-",
    paste("Negative-binomial dispersion of the oocyst count (variance lambda + lambda^2 / K); 0 is",
          "Poisson. It applies to oocysts only. Vivax oocyst counts are overdispersed."),
    min = 0, rule = "0, or at least 0.001"),
  GAMETOCYTE_DELAY = .param(0L, "days",
    paste("Days from blood-stage onset until a clone can infect mosquitoes (fixed). 0 acts as 1: a clone",
          "is first taken up the day after it appears. Vivax: about 2 to 3 days."),
    min = 0),
  P_GAMETOCYTE = .param(1, "probability",
    "An infection carries gametocytes, drawn once per infection. Vivax: about 0.9.", min = 0, max = 1),
  P_PRIMARY = .param(1, "probability",
    "A transmitted genotype causes a blood-stage primary; its hypnozoite batch is drawn either way.",
    min = 0, max = 1, source = "White et al. (2014), base model"),

  NU_HYPNOZOITE_BATCH = .param(3, "hypnozoites per bite",
    paste("Mean of the geometric hypnozoite batch an infectious bite seeds, split over the genotypes it",
          "delivers. Fitted to four target bands at moderate transmission (section 13 of the",
          "model description); re-fit it per setting."),
    min = 0, max = 1e6, rule = "times its boost, at most 1e6 once the trigger falls inside the run",
    source = "White et al. (2014), section 2"),
  ALPHA_ACTIVATE = .param(1/170, "per hypnozoite per day", "Daily probability that a hypnozoite activates: a relapse.",
    min = 0, max = 1, shown = "1/170",
    source = "an average carriage of about 170 days (Mehra et al. 2024b); White et al. (2014) estimate 1/332"),
  MU_HYPNOZOITE = .param(1/442, "per hypnozoite per day", "Daily probability that a hypnozoite dies, also while latent.",
    min = 0, max = 1, min_open = TRUE, shown = "1/442",
    source = "as in Mehra et al. (2024a), who credit White et al. (2014), whose table gives 1/425"),
  K_LATENCY = .param(0L, "stages",
    "Latency stages a batch passes before it can activate: 0 is tropical, 2 temperate.", min = 0,
    source = "White et al. (2014)"),
  DELTA_LATENCY = .param(1/100, "per day", "Daily probability that a latent batch moves on one stage.",
    min = 0, max = 1, min_open = TRUE, shown = "1/100", rule = "checked only when K_LATENCY > 0"),

  GENOME_MODE = .param("full", "-",
    paste("\"full\": the 14 PvP01 chromosomes at full length (24.2 Mb); \"core\": each at 90% of its",
          "length, a placeholder for the core genome."),
    values = c("full", "core")),
  CHROMOSOMES = .param("1-14", "-", "Which PvP01 chromosomes to simulate, such as \"1-14\" or \"1-3,5\".",
    rule = "chromosome numbers 1 to 14, each once, separated by commas, with low-high ranges"),
  RECOMB_RATE = .param(-1, "per bp",
    paste("Crossover rate within a chromosome; chromosomes assort independently. Run time grows steeply",
          "with it."),
    max = 0.5, rule = "0 or below means the default map rate, 7.4e-7", source = "Jiang et al. (2011)"),

  INTERVENTION_START = .param(99999L, "day", "First day of vector control.", shown = "99999 (never)"),
  INTERVENTION_END = .param(99999L, "day",
    "Vector control acts on days INTERVENTION_START to INTERVENTION_END - 1.", shown = "99999 (never)"),
  INTERVENTION_EFFICACY = .param(0, "fraction",
    paste("Biting is multiplied by 1 - INTERVENTION_EFFICACY inside the window; importation and",
          "treatment are not affected."),
    min = 0, max = 1),
  TREATMENT_START = .param(0L, "day",
    "First day on which a detected case can be prescribed a course; case management has its own window, separate from vector control."),
  TREATMENT_END = .param(99999L, "day",
    paste("Cases detected on days TREATMENT_START to TREATMENT_END - 1 are treated. The window gates the",
          "decision, not the delivery: a course prescribed before the end still arrives."),
    shown = "99999 (never)"),
  MDA_TICK = .param(99999L, "day",
    paste("Day of one MDA round, through the same course as case management: blood-stage drugs, and",
          "radical cure as RADCURE_COVERAGE gives it, for each host it reaches. It also reaches subpatent",
          "clones and clones that became patent that day, and it is not gated by the treatment window."),
    shown = "99999 (never)",
    rule = "not 1 when a founder host carries a copy of another founder's genome"),
  P_RAD_SURVIVAL = .param(0.29, "probability",
    paste("Radical cure fails: the patient's whole reservoir survives (RADCURE_MODE \"person\"), about",
          "71% cure; or each hypnozoite survives (\"hypnozoite\"), where the same number would mean 95%",
          "efficacy."),
    min = 0, max = 1, source = "White et al. (2018); malariasimulation"),
  RADCURE_MODE = .param("person", "-",
    paste("The unit of radical-cure failure: \"person\", or \"hypnozoite\", which thins the reservoir",
          "hypnozoite by hypnozoite and is the only mode that can clear part of it."),
    values = c("person", "hypnozoite")),
  RADCURE_COVERAGE = .param(0, "probability",
    paste("A treated case also gets radical cure: a probability given treatment, not a coverage of its",
          "own, so it does nothing unless BLOODCURE_COVERAGE > 0."),
    min = 0, max = 1),
  BLOODCURE_COVERAGE = .param(0, "probability",
    paste("A detected case is treated: a course clears each patent clone the host carries, unless the",
          "clone survives. 0 turns the whole treatment path off."),
    min = 0, max = 1),
  P_BLOODCURE_SURVIVAL = .param(0, "probability",
    "A clone survives a course: a patent one, or under MDA a subpatent one too.", min = 0, max = 1,
    source = "about 0.10 in the literature"),
  MDA_COVERAGE = .param(1, "probability", "The round reaches a host, drawn for each host.", min = 0, max = 1),
  P_RECRUD_GIVEN_FAIL = .param(0, "probability",
    paste("A clone that survives a course is suppressed below patency, to recrudesce later, rather than",
          "persisting. At 0 the third recurrence state is inert, which matches recrudescence being 0.2",
          "to 0.3% of recurrences after chloroquine or artesunate."),
    min = 0, max = 1, source = "Taylor et al. (2019)"),
  MEAN_RECRUD_DELAY = .param(10.5, "days",
    paste("Mean days from suppression to recrudescence, at least one; a literature assumption pending",
          "parameterisation. Geometric at the default: a median of 7 days, about 94% within 28 days."),
    min = 1, max = 1e6, source = "close to the prior of Taylor et al. (2019)"),
  RECRUD_DELAY_DIST = .param("geometric", "-",
    paste("Shape of that delay: \"geometric\", the discrete form of an exponential delay; \"poisson\",",
          "1 + Poisson(D - 1) with mean D, a tight-window sensitivity arm; or \"fixed\", round(D) days",
          "with halves rounded up (11 days at 10.5). Every shape is at least one day."),
    values = c("geometric", "poisson", "fixed")),
  P_DETECT_RECRUD = .param(0.5, "probability",
    "A recrudescence is detected; it should usually sit well below P_DETECT_PRIMARY.", min = 0, max = 1,
    source = "Ross et al. (2025) put recurrence detectability at 0.13 to 0.49 of a primary's"),
  SUBPATENT_INFECTIVITY = .param(0.15, "relative",
    paste("Infectiousness of a suppressed clone relative to a patent one. Neither 0 nor 1 is",
          "defensible."),
    min = 0, max = 1,
    source = "0.035 to 0.10 in the model of White et al. (2018); 0.10 to 0.20 measured (Koepfli et al. 2021)"),
  TREATMENT_DELAY = .param(7, "days",
    paste("Mean days from a case to its course: geometric, at least one day, so a course never lands on",
          "the day the case arises."),
    min = 1, max = 1e6),
  P_DETECT_PRIMARY = .param(1, "probability",
    "A primary infection or an import is detected, drawn once per bite (for all the genotypes it delivers) or per import.",
    min = 0, max = 1),
  P_DETECT_RELAPSE = .param(1, "probability",
    "A relapse is detected, drawn once per activation. Relapses are often subpatent, so a realistic value is lower.",
    min = 0, max = 1),
  SEASONAL_AMP = .param(0, "-",
    paste("Amplitude A of the seasonal factor 1 + A cos(2 pi (t - SEASONAL_PHASE) / 365), floored at 0,",
          "on biting and on the Poisson import rate; 0 means no seasonality."),
    min = 0),
  SEASONAL_PHASE = .param(0, "day of the year", "Day of the year of the seasonal peak."),
  RESURGENCE_TRIGGER_TICK = .param(99999L, "day", "Day from which the four boosts apply.", shown = "99999 (never)"),
  BITING_RATE_BOOST = .param(1, "multiplier", "Multiplies BITING_RATE from the trigger.", min = 0),
  ALPHA_ACTIVATE_BOOST = .param(1, "multiplier",
    paste("Multiplies ALPHA_ACTIVATE from the trigger, capped so that activation and death stay within",
          "one day's probability."),
    min = 0),
  NU_HYPNOZOITE_BATCH_BOOST = .param(1, "multiplier", "Multiplies NU_HYPNOZOITE_BATCH from the trigger.", min = 0),
  IMPORT_START = .param(99999L, "day", "First day of the import window.", shown = "99999 (never)"),
  IMPORT_END = .param(99999L, "day",
    "Imports in the window arrive on days max(1, IMPORT_START) to min(IMPORT_END - 1, SIM_DAYS).",
    shown = "99999 (never)"),
  N_IMPORTS_PER_TICK = .param(0L, "imports a day", "A fixed number of imports on each day of the window.",
    min = 0, rule = "warns when a day's fixed imports exceed N_HOSTS"),
  IMPORT_RATE_PER_TICK = .param(0, "imports a day",
    "Mean of a Poisson number of imports on each day of the window; it follows the season.",
    min = 0, max = 1e6, rule = "at its seasonal and boosted peak, at most 1e6 a day"),
  IMPORT_RATE_BOOST = .param(1, "multiplier", "Multiplies IMPORT_RATE_PER_TICK from RESURGENCE_TRIGGER_TICK.",
    min = 0),
  IMPORT_TICK = .param(99999L, "day", "Day of one pulse of N_IMPORTS imports.", shown = "99999 (never)"),
  N_IMPORTS = .param(0L, "imports", "Size of the pulse.", min = 0, rule = "warns above N_HOSTS")
)

## The documentation's groups, in reading order: a lead sentence, and the parameters in the order
## they are read. Every parameter is in exactly one.
.param_groups <- list(
  "Run" = list(lead = "Length, label and what is recorded.",
    params = c("SIM_DAYS", "SIM_ID", "VERBOSE", "RETAIN_FULL_HISTORY")),
  "Population" = list(lead = "Hosts, mosquitoes and the founding infections.",
    params = c("N_HOSTS", "N_MOSQ", "N_FOUNDER_HOSTS", "N_FOUNDER_LINEAGES")),
  "Transmission" = list(lead = "Bites, infection and the mosquito's incubation.",
    params = c("BITING_RATE", "P_TRANS_VH", "P_TRANS_HV", "MOSQ_LIFESPAN", "EIP_DAYS", "MAX_MOI",
               "MAX_BROODS_PER_MOSQ")),
  "Blood stage" = list(lead = "A clone's duration and its infectiousness to mosquitoes.",
    params = c("MEAN_BLOOD_STAGE_DURATION", "P_PRIMARY", "P_GAMETOCYTE", "GAMETOCYTE_DELAY")),
  "Mosquito uptake and co-transmission" = list(
    lead = paste("How many oocysts an uptake forms, and how many genotypes a bite delivers. Both lambdas",
                 "are pre-truncation rates, and 0 for both gives one oocyst per uptake and one genotype",
                 "per bite."),
    params = c("LAMBDA_OOCYSTS", "K_OOCYSTS", "LAMBDA_SPOROZOITES")),
  "Hypnozoite reservoir" = list(lead = "The liver-stage batches a bite seeds, and their activation and death.",
    params = c("NU_HYPNOZOITE_BATCH", "ALPHA_ACTIVATE", "MU_HYPNOZOITE", "K_LATENCY", "DELTA_LATENCY")),
  "Genome" = list(lead = "The neutral PvP01 genome the model simulates; the mutation rate is 0.",
    params = c("GENOME_MODE", "CHROMOSOMES", "RECOMB_RATE")),
  "Seasonality and vector control" = list(
    lead = "Factors on the biting rate; seasonality also acts on the Poisson import rate.",
    params = c("SEASONAL_AMP", "SEASONAL_PHASE", "INTERVENTION_START", "INTERVENTION_END",
               "INTERVENTION_EFFICACY")),
  "Detection and case management" = list(
    lead = paste("A case is detected, decided where it arises and treated later, on the host; radical",
                 "cure is nested in a blood-stage course. Off by default."),
    params = c("P_DETECT_PRIMARY", "P_DETECT_RELAPSE", "BLOODCURE_COVERAGE", "P_BLOODCURE_SURVIVAL",
               "RADCURE_COVERAGE", "P_RAD_SURVIVAL", "RADCURE_MODE", "TREATMENT_DELAY", "TREATMENT_START",
               "TREATMENT_END")),
  "Mass drug administration" = list(lead = "One round, through the same course as case management.",
    params = c("MDA_TICK", "MDA_COVERAGE")),
  "Recrudescence" = list(
    lead = "A clone that survives a course can be suppressed below patency and return. Off by default.",
    params = c("P_RECRUD_GIVEN_FAIL", "MEAN_RECRUD_DELAY", "RECRUD_DELAY_DIST", "P_DETECT_RECRUD",
               "SUBPATENT_INFECTIVITY")),
  "Resurgence" = list(lead = "Multipliers that switch on from one day, to start a resurgence.",
    params = c("RESURGENCE_TRIGGER_TICK", "BITING_RATE_BOOST", "ALPHA_ACTIVATE_BOOST",
               "NU_HYPNOZOITE_BATCH_BOOST")),
  "Importation" = list(
    lead = paste("Three routes, which add up: a pulse, and a fixed daily count and a Poisson rate in a",
                 "window. Each host takes at most one import a day, and a host at MAX_MOI takes none."),
    params = c("IMPORT_TICK", "N_IMPORTS", "IMPORT_START", "IMPORT_END", "N_IMPORTS_PER_TICK",
               "IMPORT_RATE_PER_TICK", "IMPORT_RATE_BOOST")))

## ---- the documentation, generated ----------------------------------------------------------
.param_num <- function(x) format(x, scientific = FALSE, trim = TRUE, digits = 15)
.param_cap <- function(x) if (x >= 1e6) sub("e\\+0?", "e", format(x, scientific = TRUE)) else .param_num(x)

.param_group_of <- function(name)
  names(.param_groups)[vapply(.param_groups, function(g) name %in% g$params, logical(1))]

.param_type <- function(r) {
  if (is.logical(r$default)) "TRUE or FALSE" else if (is.character(r$default)) "text"
  else if (is.integer(r$default)) "whole number" else "number"
}

## The default as documented.
.param_default_text <- function(r) {
  if (!is.na(r$shown)) return(r$shown)
  d <- r$default
  if (is.character(d)) paste0("\"", d, "\"") else if (is.logical(d)) as.character(d) else .param_num(d)
}

## What a value may be, in words: the type and the single-value bounds, then the rule.
.param_allowed <- function(r) {
  base <- if (length(r$values)) paste(paste0("\"", r$values, "\""), collapse = " or ")
  else if (!is.numeric(r$default)) .param_type(r)
  else {
    lo <- r$min; hi <- r$max
    rng <- if (!is.na(lo) && !is.na(hi))
      paste0("in ", if (r$min_open) "(" else "[", .param_cap(lo), ", ", .param_cap(hi), if (r$max_open) ")" else "]")
    else if (!is.na(lo)) paste(if (r$min_open) ">" else ">=", .param_cap(lo))
    else if (!is.na(hi)) paste(if (r$max_open) "<" else "<=", .param_cap(hi))
    else ""
    trimws(paste(.param_type(r), rng))
  }
  if (!is.na(r$rule)) paste0(base, "; ", r$rule) else base
}

.param_text <- function(name) {
  r <- .param_dict[[name]]
  paste0(r$meaning, if (r$unit != "-") paste0(" Unit: ", r$unit, "."), " Default ", .param_default_text(r),
         "; ", .param_allowed(r), ".", if (!is.na(r$source)) paste0(" Source: ", r$source, "."))
}

.rd_escape <- function(s) gsub("%", "\\\\%", gsub("([{}\\\\])", "\\\\\\1", s))

## The Parameters section of ?vivax_params: a subsection per group, a list of its parameters.
.param_rd <- function() {
  sub <- vapply(names(.param_groups), function(g) {
    items <- vapply(.param_groups[[g]]$params, function(nm)
      sprintf("\\item{\\code{%s}}{%s}", nm, .rd_escape(.param_text(nm))), "")
    sprintf("\\subsection{%s}{%s\n\\describe{\n%s\n}}", g, .rd_escape(.param_groups[[g]]$lead),
            paste(items, collapse = "\n"))
  }, "")
  sprintf("\\section{Parameters}{\nEvery parameter, by group. Section 10 of the \\href{https://shaziaruybal.github.io/malariafwd/MODEL_DESCRIPTION.html}{model description} holds the same tables.\n\n%s\n}",
          paste(sub, collapse = "\n\n"))
}

## The tables of MODEL_DESCRIPTION section 10: a heading, the lead and a table per group.
.param_md <- function() {
  out <- character()
  for (g in names(.param_groups)) {
    out <- c(out, paste0("### ", g), "", .param_groups[[g]]$lead, "",
             "| parameter | default | unit | allowed | meaning |", "|---|---|---|---|---|")
    for (nm in .param_groups[[g]]$params) {
      r <- .param_dict[[nm]]
      m <- if (is.na(r$source)) r$meaning else paste0(r$meaning, " Source: ", r$source, ".")
      out <- c(out, sprintf("| `%s` | %s | %s | %s | %s |", nm, .param_default_text(r), r$unit,
                            .param_allowed(r), gsub("|", "\\|", m, fixed = TRUE)))
    }
    out <- c(out, "")
  }
  paste(out, collapse = "\n")
}

#' What each parameter means, its unit and the values it may take
#'
#' The model's parameters as documented on [vivax_params()]: one row per parameter, by group.
#' @param params `NULL` for every parameter, or the names of the ones to give, in that order.
#' @return a `data.table` of character columns: `name`, `group`, `default` (as documented), `type`
#'   (`"whole number"`, `"number"`, `"TRUE or FALSE"` or `"text"`), `unit`, `allowed`, `meaning`
#'   and `source` (`NA` when none is given).
#' @examples
#' param_info(c("N_HOSTS", "BITING_RATE"))
#' @export
param_info <- function(params = NULL) {
  order <- unlist(lapply(.param_groups, `[[`, "params"), use.names = FALSE)
  if (!is.null(params)) {
    if (!is.character(params)) stop("`params` must be parameter names.", call. = FALSE)
    .refuse_unknown(params)
    order <- params
  }
  data.table::rbindlist(lapply(order, function(nm) {
    r <- .param_dict[[nm]]
    list(name = nm, group = .param_group_of(nm), default = .param_default_text(r), type = .param_type(r),
         unit = r$unit, allowed = .param_allowed(r), meaning = r$meaning, source = as.character(r$source))
  }))
}
