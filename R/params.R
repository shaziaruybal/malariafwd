## params.R -- the user-facing parameter surface. Defaults MIRROR the SLiM model's
## defineDefault() values (inst/slim/malariafwd_vivax.slim is the source of truth); the
## test in tests/testthat/ checks they stay in sync.

# All model defaults, grouped. Values match the SLiM defineDefault() calls exactly.
.malariafwd_defaults <- list(
  # run / meta
  SIM_ID = "default", SIM_DAYS = 365L, VERBOSE = FALSE, RETAIN_FULL_HISTORY = TRUE,
  # population
  N_HOSTS = 40L, N_FOUNDERS = 12L, N_FOUNDER_HOSTS = 12L, N_FOUNDER_LINEAGES = 12L, N_MOSQ = 100L,
  # transmission (Ross-Macdonald)
  BITING_RATE = 0.1, P_TRANS_VH = 0.5, P_TRANS_HV = 0.5, MOSQ_LIFESPAN = 10.0,
  EIP_DAYS = 10L, MAX_MOI = 20L, MAX_BROODS_PER_MOSQ = 20L,
  # parasite / infection
  MEAN_BLOOD_STAGE_DURATION = 60, N_OOCYSTS_PER_BITE = 1L,
  # mosquito uptake / inoculation sampling . lambda = 0 (default) -> exactly 1
  # oocyst per uptake and 1 genotype per bite. K_OOCYSTS > 0 makes the oocyst count
  # negative-binomial (var = lambda + lambda^2/K); 0 = Poisson.
  LAMBDA_OOCYSTS = 0.0, LAMBDA_SPOROZOITES = 0.0, K_OOCYSTS = 0.0,
  GAMETOCYTE_DELAY = 0L, P_GAMETOCYTE = 1.0, P_PRIMARY = 1.0,
  # hypnozoite reservoir
  NU_HYPNOZOITE_BATCH = 5.0, ALPHA_ACTIVATE = 1/170, MU_HYPNOZOITE = 1/442,
  K_LATENCY = 0L, DELTA_LATENCY = 1/100,
  # genome (see DECISIONS D2/D3/D13): pure coordinate switch; both modes rate 7.4e-7
  GENOME_MODE = "full", CHROMOSOMES = "1-14", RECOMB_RATE = -1.0,
  # vector control: biting rate and importation rate only, never case management
  INTERVENTION_START = 99999L, INTERVENTION_END = 99999L, INTERVENTION_EFFICACY = 0.0,
  ## case management has its own window (D19); INTERVENTION_* is vector control only
  TREATMENT_START = 0L, TREATMENT_END = 99999L,
  # case management: radical cure is NESTED inside blood-stage treatment, so
  # RADCURE_COVERAGE is P(radical cure | treated). P_RAD_SURVIVAL is per PERSON by
  # default - see RADCURE_MODE.
  MDA_TICK = 99999L, P_RAD_SURVIVAL = 0.29, RADCURE_MODE = "person",
  RADCURE_COVERAGE = 0.0, BLOODCURE_COVERAGE = 0.0, P_BLOODCURE_SURVIVAL = 0.0,
  MDA_COVERAGE = 1.0,
  P_RECRUD_GIVEN_FAIL = 0.0, MEAN_RECRUD_DELAY = 10.5, RECRUD_DELAY_DIST = "geometric",
  P_DETECT_RECRUD = 0.5,
  SUBPATENT_INFECTIVITY = 0.15,
  TREATMENT_DELAY = 7.0,
  # detection / observation
  P_DETECT_PRIMARY = 1.0, P_DETECT_RELAPSE = 1.0,
  # seasonality
  SEASONAL_AMP = 0.0, SEASONAL_PHASE = 0.0,
  # resurgence triggers
  RESURGENCE_TRIGGER_TICK = 99999L, BITING_RATE_BOOST = 1.0,
  ALPHA_ACTIVATE_BOOST = 1.0, NU_HYPNOZOITE_BATCH_BOOST = 1.0,
  # importation
  IMPORT_START = 99999L, IMPORT_END = 99999L, N_IMPORTS_PER_TICK = 0L,
  IMPORT_RATE_PER_TICK = 0.0, IMPORT_RATE_BOOST = 1.0, IMPORT_TICK = 99999L, N_IMPORTS = 0L
)

# Params that are Eidos STRINGS (must be quoted when passed to SLiM -d).
.malariafwd_string_params <- c("SIM_ID", "GENOME_MODE", "CHROMOSOMES", "RADCURE_MODE",
                               "RECRUD_DELAY_DIST", "OUTDIR")

#' Build a validated parameter set for the vivax model
#'
#' Returns a named list of all model parameters, starting from the defaults (which
#' mirror the SLiM model) and overriding any you name. `SEED` and `OUTDIR` are set by
#' [run_sim()], not here.
#'
#' Time is in days: one model tick is one day. A duration parameter `D` (such as
#' `MEAN_BLOOD_STAGE_DURATION` or `MOSQ_LIFESPAN`) is the mean of a geometric duration: the
#' event happens each day with probability `1/D`. The treatment delay (`TREATMENT_DELAY`) and,
#' under the default `RECRUD_DELAY_DIST = "geometric"`, the recrudescence delay
#' (`MEAN_RECRUD_DELAY`) follow the same rule from the day after the case or the suppression,
#' so they are at least one day; `EIP_DAYS` and `GAMETOCYTE_DELAY` are fixed delays. Rates
#' (`BITING_RATE`, `ALPHA_ACTIVATE`, `MU_HYPNOZOITE`, `DELTA_LATENCY`) are daily probabilities;
#' a mosquito bites at most once a day, so `BITING_RATE` times its seasonal and resurgence
#' multipliers is capped at 1.
#'
#' Values the model cannot run are refused here, before SLiM starts. Each value is a number,
#' `TRUE`/`FALSE` or a string, as its default is, and numbers are finite. Counts, `SIM_DAYS`,
#' the scheduling days (`MDA_TICK`, `IMPORT_TICK`, `RESURGENCE_TRIGGER_TICK` and the `*_START`
#' and `*_END` bounds) and the fixed delays `EIP_DAYS` and `GAMETOCYTE_DELAY` are whole numbers;
#' mean durations, rates and `SEASONAL_PHASE` need not be. Probabilities and coverages,
#' `ALPHA_ACTIVATE` and `INTERVENTION_EFFICACY` lie in \[0, 1\], `MU_HYPNOZOITE` in (0, 1\], and
#' `DELTA_LATENCY` in (0, 1\] when `K_LATENCY > 0`. `BITING_RATE` is positive (above 1 it only
#' warns). Boosts, import rates, `LAMBDA_*`, `SEASONAL_AMP` and counts are non-negative.
#' `RECOMB_RATE` is at most 0.5 per bp (any value <= 0 means the default map rate), and run time
#' grows steeply with it.
#'
#' A few sizes are capped far above anything realistic: a mean of 1000 oocysts per successful
#' uptake (see [uptake_moments()]), with `K_OOCYSTS` either 0 or at least 0.001;
#' `LAMBDA_SPOROZOITES` at 1000; the hypnozoite batch mean at 1e6, times its boost once the
#' trigger falls inside the run; the import rate at 1e6 a day; and `TREATMENT_DELAY` and
#' `MEAN_RECRUD_DELAY` at 1e6 days. `SIM_DAYS` is at most 1e9, SLiM's last tick, and a run that
#' reaches day 99999, the day the defaults use for "never", must move those days past its end. A
#' founder-free run whose imports never establish a parasite is not refused here; the model's
#' liveness check stops it.
#'
#' A whole number carrying floating-point noise (such as `1.1 * 100`), a number or factor given
#' for a string parameter, and spaces around the commas and dashes of `CHROMOSOMES` are stored
#' in the form SLiM is given.
#'
#' @param ... named parameter overrides, e.g. `N_HOSTS = 300, BITING_RATE = 0.09`.
#' @return An object of class `vivax_params` (a named list).
#' @examples
#' vivax_params(N_HOSTS = 300, BITING_RATE = 0.09, GENOME_MODE = "core")
#' @export
vivax_params <- function(...) {
  over <- list(...)
  if (length(over) && (is.null(names(over)) || !all(nzchar(names(over)))))
    stop("every argument must be named, e.g. vivax_params(N_HOSTS = 300).", call. = FALSE)
  unknown <- setdiff(names(over), names(.malariafwd_defaults))
  if (length(unknown))
    stop("unknown parameter(s): ", paste(unknown, collapse = ", "),
         "\nsee names(vivax_params()) for valid names.", call. = FALSE)
  twice <- unique(names(over)[duplicated(names(over))])
  if (length(twice))
    stop(paste(twice, collapse = ", "), " given more than once.", call. = FALSE)
  ## modifyList() would drop a NULL, leaving SLiM's own default in its place
  nulls <- names(over)[vapply(over, is.null, logical(1))]
  if (length(nulls))
    stop(paste(nulls, collapse = ", "), " cannot be NULL: leave a parameter out to use its default.",
         call. = FALSE)
  p <- utils::modifyList(.malariafwd_defaults, over)
  ## As in the model, N_FOUNDER_HOSTS defaults to N_FOUNDERS, and N_FOUNDER_LINEAGES to
  ## N_FOUNDER_HOSTS
  if (!"N_FOUNDER_HOSTS"    %in% names(over)) p$N_FOUNDER_HOSTS    <- p$N_FOUNDERS
  if (!"N_FOUNDER_LINEAGES" %in% names(over)) p$N_FOUNDER_LINEAGES <- p$N_FOUNDER_HOSTS
  p <- .normalise_params(p)
  validate_params(p)
  structure(p, class = "vivax_params")
}

## A number as as_slim_args() writes it: 15 significant digits, never in scientific notation.
.slim_num <- function(x) format(x, scientific = FALSE, trim = TRUE, digits = 15)

## Spellings that mean the same value are stored in the form SLiM is given, so the object,
## params.tsv and the run agree: a whole number with floating-point noise (1.1 * 100 is
## 110.00000000000001, written as 110), a number or factor for a string parameter (expand.grid()
## makes factors), and spaces around the commas and dashes of CHROMOSOMES. Spaces anywhere else
## are left for the check to refuse: "1 4" is not "14".
.normalise_params <- function(p) {
  d <- .malariafwd_defaults
  for (nm in intersect(names(p), names(d))) {
    v <- p[[nm]]
    if (length(v) != 1L || is.na(v)) next
    if (is.integer(d[[nm]]) && is.double(v) && is.finite(v) && grepl("^-?[0-9]+$", .slim_num(v)))
      p[[nm]] <- round(v)
    else if (is.character(d[[nm]]) && is.factor(v))
      p[[nm]] <- as.character(v)
    else if (is.character(d[[nm]]) && is.numeric(v))
      p[[nm]] <- .slim_num(v)
  }
  ch <- p$CHROMOSOMES
  if (is.character(ch) && length(ch) == 1L && !is.na(ch) &&
      grepl("^ *[0-9]+ *(- *[0-9]+ *)?(, *[0-9]+ *(- *[0-9]+ *)?)*$", ch))
    p$CHROMOSOMES <- gsub(" ", "", ch, fixed = TRUE)
  p
}

## The days that default to 99999, which the model reads as "never".
.never_params <- names(.malariafwd_defaults)[vapply(.malariafwd_defaults, identical, logical(1), 99999L)]

## The largest daily biting probability the parameters ask of a mosquito, worked out over the run
## the way the model applies them: seasonality only when SEASONAL_AMP > 0 (floored at 0), the
## intervention factor inside [INTERVENTION_START, INTERVENTION_END), the boost from
## RESURGENCE_TRIGGER_TICK. Biting starts on tick 2, once the populations exist. The factors change
## only on those days and the seasonal one repeats every 365 days, so a year from each change
## covers every value the run can take.
.peak_biting <- function(p) {
  last <- max(2, p$SIM_DAYS)
  starts <- c(2, p$INTERVENTION_START, p$INTERVENTION_END, p$RESURGENCE_TRIGGER_TICK)
  starts <- starts[starts >= 2 & starts <= last]
  t <- sort(unique(unlist(lapply(starts, function(a) seq(a, min(last, a + 364))))))
  r <- rep(p$BITING_RATE, length(t))
  if (p$SEASONAL_AMP > 0) {
    angle <- 2 * pi * (t - p$SEASONAL_PHASE) / 365
    if (!all(is.finite(angle))) return(c(rate = NaN, day = NA))
    r <- r * pmax(0, 1 + p$SEASONAL_AMP * cos(angle))
  }
  inside <- t >= p$INTERVENTION_START & t < p$INTERVENTION_END
  r[inside] <- r[inside] * (1 - p$INTERVENTION_EFFICACY)
  boosted <- t >= p$RESURGENCE_TRIGGER_TICK
  r[boosted] <- r[boosted] * p$BITING_RATE_BOOST
  c(rate = max(r), day = t[which.max(r)])
}

#' @keywords internal
validate_params <- function(p) {
  bad <- names(p)[!vapply(p, function(v) length(v) == 1L && !is.na(v), logical(1))]
  if (length(bad))
    stop(paste(bad, collapse = ", "), " must be a single non-missing value", call. = FALSE)
  .check_types(p)
  if (!p$GENOME_MODE %in% c("full", "core"))
    stop("GENOME_MODE must be 'full' or 'core', got '", p$GENOME_MODE, "'", call. = FALSE)
  .check_chromosomes(p$CHROMOSOMES)
  if (p$RECOMB_RATE > 0.5)
    stop("RECOMB_RATE must be at most 0.5 per bp, got ", p$RECOMB_RATE, " (a value <= 0 means the ",
         "default map rate, 7.4e-7).", call. = FALSE)
  for (pr in c("P_TRANS_VH", "P_TRANS_HV", "P_GAMETOCYTE", "P_PRIMARY",
               "P_DETECT_PRIMARY", "P_DETECT_RELAPSE", "RADCURE_COVERAGE", "P_RAD_SURVIVAL",
               "BLOODCURE_COVERAGE", "P_BLOODCURE_SURVIVAL",
               "P_RECRUD_GIVEN_FAIL", "P_DETECT_RECRUD", "SUBPATENT_INFECTIVITY",
               "MDA_COVERAGE", "ALPHA_ACTIVATE", "INTERVENTION_EFFICACY"))
    if (p[[pr]] < 0 || p[[pr]] > 1) stop(pr, " must be in [0, 1]", call. = FALSE)
  for (pr in c("SIM_DAYS", "N_HOSTS", "N_MOSQ", "BITING_RATE"))
    if (p[[pr]] <= 0) stop(pr, " must be > 0", call. = FALSE)
  for (pr in c("LAMBDA_OOCYSTS", "LAMBDA_SPOROZOITES", "K_OOCYSTS", "NU_HYPNOZOITE_BATCH",
               "ALPHA_ACTIVATE_BOOST", "NU_HYPNOZOITE_BATCH_BOOST", "BITING_RATE_BOOST",
               "IMPORT_RATE_PER_TICK", "IMPORT_RATE_BOOST", "N_IMPORTS", "N_IMPORTS_PER_TICK",
               "N_FOUNDERS", "N_FOUNDER_HOSTS", "N_FOUNDER_LINEAGES", "EIP_DAYS",
               "GAMETOCYTE_DELAY", "MAX_MOI", "MAX_BROODS_PER_MOSQ", "K_LATENCY"))
    if (p[[pr]] < 0) stop(pr, " must be >= 0", call. = FALSE)
  ## N_FOUNDERS itself is only the default of N_FOUNDER_HOSTS, so it is not checked
  if (p$N_FOUNDER_HOSTS > p$N_HOSTS)
    stop("N_FOUNDER_HOSTS (", p$N_FOUNDER_HOSTS, ") must be at most N_HOSTS (", p$N_HOSTS, "). ",
         "It defaults to N_FOUNDERS, so set either one.", call. = FALSE)
  if (p$N_FOUNDER_LINEAGES > p$N_FOUNDER_HOSTS)
    stop("N_FOUNDER_LINEAGES (", p$N_FOUNDER_LINEAGES, ") must be at most N_FOUNDER_HOSTS (",
         p$N_FOUNDER_HOSTS, "): each founder host carries one of the founding genomes.", call. = FALSE)
  if (p$SEASONAL_AMP < 0)
    stop("SEASONAL_AMP must be >= 0 (0 means no seasonality). SEASONAL_PHASE is the day of the ",
         "peak, so for a trough on day x use SEASONAL_PHASE = x + 182.5.", call. = FALSE)
  if (p$SIM_DAYS > 1e9)
    stop("SIM_DAYS must be at most 1e9, the last tick SLiM can schedule.", call. = FALSE)
  never <- .never_params[vapply(.never_params, function(nm) p[[nm]] == 99999, logical(1))]
  if (p$SIM_DAYS >= 99999 && length(never))
    stop("SIM_DAYS = ", .slim_num(p$SIM_DAYS), " reaches day 99999, which the defaults use for ",
         "'never', and ", paste(never, collapse = ", "), if (length(never) > 1) " are" else " is",
         " still 99999: set each to a day after SIM_DAYS, or to the day you mean.", call. = FALSE)
  if (!is.finite(p$MOSQ_LIFESPAN) || p$MOSQ_LIFESPAN < 1)
    stop("MOSQ_LIFESPAN must be >= 1 and finite (days): a mosquito dies each day with probability ",
         "1/MOSQ_LIFESPAN.", call. = FALSE)
  .check_sizes(p)
  if (p$MU_HYPNOZOITE <= 0)
    stop("MU_HYPNOZOITE must be > 0: hypnozoites that cannot die keep their batch alive forever.",
         call. = FALSE)
  if (p$MU_HYPNOZOITE > 1)
    stop("MU_HYPNOZOITE must be <= 1: it is a daily probability.", call. = FALSE)
  if (p$K_LATENCY > 0 && (p$DELTA_LATENCY <= 0 || p$DELTA_LATENCY > 1))
    stop("DELTA_LATENCY must be in (0, 1] when K_LATENCY > 0: a latent batch leaves each latency ",
         "stage with daily probability DELTA_LATENCY.", call. = FALSE)
  if (!is.finite(p$MEAN_BLOOD_STAGE_DURATION) || p$MEAN_BLOOD_STAGE_DURATION < 1)
    stop("MEAN_BLOOD_STAGE_DURATION must be >= 1 and finite (days): blood-stage infections clear each day ",
         "with probability 1/MEAN_BLOOD_STAGE_DURATION.", call. = FALSE)
  if (p$N_OOCYSTS_PER_BITE != 1)
    stop("N_OOCYSTS_PER_BITE is deprecated; use LAMBDA_SPOROZOITES instead.", call. = FALSE)
  if (!p$RADCURE_MODE %in% c("person", "hypnozoite"))
    stop("RADCURE_MODE must be 'person' or 'hypnozoite', got '", p$RADCURE_MODE, "'",
         call. = FALSE)
  if (!is.finite(p$TREATMENT_DELAY) || p$TREATMENT_DELAY < 1 || p$TREATMENT_DELAY > 1e6)
    stop("TREATMENT_DELAY must be >= 1 and finite (days, at most 1e6): a course arrives each day ",
         "after the case with probability 1/TREATMENT_DELAY, so 1 means the next day.", call. = FALSE)
  if (!is.finite(p$MEAN_RECRUD_DELAY) || p$MEAN_RECRUD_DELAY < 1 || p$MEAN_RECRUD_DELAY > 1e6)
    stop("MEAN_RECRUD_DELAY must be >= 1 and finite (days, at most 1e6): it is the mean delay from ",
         "suppression to recrudescence, which is at least one day.", call. = FALSE)
  if (!p$RECRUD_DELAY_DIST %in% c("geometric", "poisson", "fixed"))
    stop("RECRUD_DELAY_DIST must be 'geometric', 'poisson' or 'fixed', got '",
         p$RECRUD_DELAY_DIST, "'", call. = FALSE)

  ## ---- configurations that cannot produce a single parasite -------------------------------
  ## These used to run to completion with exit status 0 and an EMPTY event log, so a
  ## mis-specified cell in a sweep was indistinguishable from one that simply saw no
  ## transmission. Each is refused before SLiM is launched.
  if (p$SIM_DAYS < 2)
    stop("SIM_DAYS must be >= 2 (founder seeding completes on tick 2), got ", p$SIM_DAYS,
         call. = FALSE)
  if (p$N_FOUNDER_HOSTS > 0 && p$N_FOUNDER_LINEAGES < 1)
    stop("N_FOUNDER_LINEAGES must be >= 1 when N_FOUNDER_HOSTS > 0: ", p$N_FOUNDER_HOSTS,
         " founder host(s) cannot share 0 distinct genomes.", call. = FALSE)
  ## A naive population IS legitimate (D17: importation can now seed an empty population),
  ## but only if an import path is actually open. No founders AND no imports is a run that
  ## can never contain a parasite.
  if (p$N_FOUNDER_HOSTS < 1 && !has_import_path(p))
    stop("nothing can start this run: N_FOUNDER_HOSTS is 0 and no importation is configured.\n",
         "  Set N_FOUNDER_HOSTS >= 1, or open an import path (IMPORT_TICK + N_IMPORTS, or\n",
         "  IMPORT_START/IMPORT_END with N_IMPORTS_PER_TICK or IMPORT_RATE_PER_TICK).",
         call. = FALSE)

  ## A coverage set against a closed case-management window gives silently zero treatments.
  ## This is the shape of the bug D19 fixed, so it is worth refusing rather than repeating in
  ## a new parameter.
  if ((p$RADCURE_COVERAGE > 0 || p$BLOODCURE_COVERAGE > 0) &&
      (p$TREATMENT_END <= p$TREATMENT_START || p$TREATMENT_START > p$SIM_DAYS))
    stop("a treatment coverage is set but the case-management window is closed:\n",
         "  TREATMENT_START = ", p$TREATMENT_START, ", TREATMENT_END = ", p$TREATMENT_END,
         ", SIM_DAYS = ", p$SIM_DAYS, "\n",
         "  No case would ever be treated. Leave TREATMENT_START/END at their defaults for\n",
         "  treatment throughout the run, or set a window that falls inside it.", call. = FALSE)

  ## ---- warnings, once every refusal has had its turn ---------------------------------------
  peak <- .peak_biting(p)
  if (!is.finite(peak[["rate"]]))
    stop("BITING_RATE with its seasonal, intervention and resurgence multipliers cannot be ",
         "worked out (it overflows); check BITING_RATE, BITING_RATE_BOOST and SEASONAL_PHASE.",
         call. = FALSE)
  if (peak[["rate"]] > 1)
    warning("BITING_RATE with its seasonal, intervention and resurgence multipliers reaches ",
            signif(peak[["rate"]], 3), " on day ", peak[["day"]], ", but a mosquito bites at most ",
            "once a day, so the realised biting rate is capped below what the parameters imply.",
            call. = FALSE)
  if (p$P_RECRUD_GIVEN_FAIL > 0 && p$P_BLOODCURE_SURVIVAL <= 0)
    warning("P_RECRUD_GIVEN_FAIL > 0 has no effect while P_BLOODCURE_SURVIVAL is 0: ",
            "no course ever fails, so no clone is ever suppressed.", call. = FALSE)
  ## Radical cure is nested inside blood-stage treatment, so a coverage without one is
  ## silently inert - the same failure shape as the window bug in D19. But MDA delivers radical
  ## cure to everyone it reaches regardless of BLOODCURE_COVERAGE, so an MDA-only scenario is a
  ## legitimate configuration and must not be warned about: the fence was firing on exactly the
  ## scenario MDA exists for.
  mda_on <- p$MDA_TICK <= p$SIM_DAYS && p$MDA_COVERAGE > 0
  if (p$RADCURE_COVERAGE > 0 && p$BLOODCURE_COVERAGE == 0 && !mda_on)
    warning("RADCURE_COVERAGE > 0 has no effect while BLOODCURE_COVERAGE is 0 and no MDA round\n",
            "  is scheduled: case-based radical cure is given only to treated cases (it is never\n",
            "  given without a schizonticide). Set BLOODCURE_COVERAGE, or schedule MDA_TICK.",
            call. = FALSE)
  invisible(TRUE)
}

## Each value must have its default's type, since SLiM reads numbers, T/F and strings
## differently and a wrong type stops the run at its first use. Numbers must be finite and
## survive SLiM's command line, which cannot read a number between 0 and about 2.2e-308; the
## floor of 1e-300 leaves room for the 15 digits as_slim_args() writes. Where the
## default is an integer (a count, SIM_DAYS, a day the run is scheduled on, or a fixed delay) the
## value must be whole. A non-whole count stops the run or silently changes it (12.5 founder
## hosts seed 13), a non-whole SIM_DAYS stops SLiM at launch, a non-whole MDA_TICK or
## IMPORT_TICK never comes, and the other days and delays act as the next whole day.
.check_types <- function(p) {
  d <- .malariafwd_defaults
  day <- function(nm) nm %in% c("MDA_TICK", "IMPORT_TICK", "RESURGENCE_TRIGGER_TICK") ||
    grepl("_(START|END)$", nm)
  for (nm in intersect(names(p), names(d))) {
    v <- p[[nm]]; def <- d[[nm]]
    want <- if (is.numeric(def)) "a number" else if (is.logical(def)) "TRUE or FALSE" else "a string"
    ok <- if (is.numeric(def)) is.numeric(v) else if (is.logical(def)) is.logical(v) else is.character(v)
    if (!ok) stop(nm, " must be ", want, ", got ", deparse(v), call. = FALSE)
    ## byte-wise, so that accented and other non-ASCII labels pass in every locale
    if (is.character(v) && grepl("[\001-\037\177]", v, useBytes = TRUE))
      stop(nm, " must not contain control characters (newlines, tabs): SLiM and params.tsv ",
           "cannot carry them.", call. = FALSE)
    if (!is.numeric(v)) next
    if (!is.finite(v))
      stop(nm, " must be finite, got ", v,
           if (v > 0 && day(nm)) "; for 'never', use a day after SIM_DAYS"
           else if (v < 0 && grepl("_START$", nm)) "; for 'from the start', use 0"
           else if (nm == "K_OOCYSTS") "; use 0 for Poisson, the K -> infinity limit" else "",
           call. = FALSE)
    if (v != 0 && abs(v) < 1e-300)
      stop(nm, " = ", v, " is too small to pass to SLiM: values between 0 and 1e-300 are refused.",
           call. = FALSE)
    if (is.integer(def) && v != round(v))
      stop(nm, " must be a whole number, got ", .slim_num(v), "; use round(), floor() or ceiling().",
           call. = FALSE)
    if (is.integer(def) && abs(v) > 2^53)
      stop(nm, " must lie within +/- 2^53 (about 9e15), got ", v,
           if (v > 0 && day(nm)) "; for 'never', a day just after SIM_DAYS is enough" else "",
           call. = FALSE)
  }
}

## Sizes far beyond anything realistic, which would make a run crawl or fill memory, capped well
## below where SLiM's draws fail (binomial sizes are 32-bit; Poisson draws wrap above 2^32).
.check_sizes <- function(p) {
  if (p$LAMBDA_SPOROZOITES > 1000)
    stop("LAMBDA_SPOROZOITES must be <= 1000: it is the pre-truncation rate of oocyst products ",
         "sampled per infectious bite.", call. = FALSE)
  ## The gamma-Poisson sampler stops drawing the distribution it is asked for below this K:
  ## its gamma draws underflow to 0, and the run falls back to one oocyst per uptake.
  if (p$K_OOCYSTS > 0 && p$K_OOCYSTS < 1e-3)
    stop("K_OOCYSTS must be 0 (Poisson) or at least 0.001, got ", p$K_OOCYSTS, ": the model ",
         "cannot draw a more dispersed oocyst count.", call. = FALSE)
  ## A cap on the mean, not on LAMBDA_OOCYSTS alone: a small K gives a heavy tail whose mean, and
  ## so the work per uptake, can be many times lambda.
  m <- uptake_moments(p$LAMBDA_OOCYSTS, p$K_OOCYSTS)[["mean"]]
  if (m > 1000)
    stop("LAMBDA_OOCYSTS = ", p$LAMBDA_OOCYSTS, " with K_OOCYSTS = ", p$K_OOCYSTS, " gives a mean ",
         "of ", signif(m, 6), " oocysts per successful uptake (see uptake_moments()); the model ",
         "accepts at most 1000.", call. = FALSE)
  boost <- function(b) if (p$RESURGENCE_TRIGGER_TICK <= p$SIM_DAYS) max(1, b) else 1
  if (p$NU_HYPNOZOITE_BATCH * boost(p$NU_HYPNOZOITE_BATCH_BOOST) > 1e6)
    stop("NU_HYPNOZOITE_BATCH, times NU_HYPNOZOITE_BATCH_BOOST once the trigger falls inside the run, ",
         "must be at most 1e6 hypnozoites per bite, far below the 32-bit sizes of SLiM's binomial ",
         "draws.", call. = FALSE)
  peak <- p$IMPORT_RATE_PER_TICK * boost(p$IMPORT_RATE_BOOST) * (1 + p$SEASONAL_AMP)
  if (peak > 1e6)
    stop("IMPORT_RATE_PER_TICK x IMPORT_RATE_BOOST x (1 + SEASONAL_AMP) is ", signif(peak, 6),
         "; it must be at most 1e6 imports a day.", call. = FALSE)
}

## The model accepts "1-14", "1,3,5" or "5-8": comma-separated chromosome numbers or low-high
## ranges, each chromosome once. The ends are checked before a range is expanded.
.check_chromosomes <- function(x) {
  bad <- function(why) stop("CHROMOSOMES must be chromosome numbers 1-14, comma-separated, with ",
                            "low-high ranges, e.g. \"1-14\" or \"1-3,5\" (", why, "); got \"", x,
                            "\".", call. = FALSE)
  if (!grepl("^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$", x)) bad("unreadable")
  idx <- unlist(lapply(strsplit(strsplit(x, ",")[[1]], "-"), function(r) {
    r <- as.numeric(r)
    if (any(r < 1 | r > 14)) bad("out of range")
    if (length(r) == 2 && r[1] > r[2]) bad("a range runs high-low")
    seq(r[1], r[length(r)])
  }))
  if (anyDuplicated(idx)) bad("a chromosome is listed twice")
}

#' Is any importation path open in this parameter set?
#'
#' Either the one-shot pulse (`IMPORT_TICK` + `N_IMPORTS`) or the windowed rate
#' (`IMPORT_START`/`IMPORT_END` with `N_IMPORTS_PER_TICK` or `IMPORT_RATE_PER_TICK`);
#' in both cases the window has to fall inside the run.
#' @param p a parameter list.
#' @return `TRUE` if at least one importation can occur.
#' @keywords internal
has_import_path <- function(p) {
  pulse <- p$IMPORT_TICK <= p$SIM_DAYS && p$N_IMPORTS > 0
  windowed <- p$IMPORT_START <= p$SIM_DAYS && p$IMPORT_END > p$IMPORT_START &&
    (p$N_IMPORTS_PER_TICK > 0 || p$IMPORT_RATE_PER_TICK > 0)
  isTRUE(pulse) || isTRUE(windowed)
}

## An Eidos single-quoted string literal: backslashes first, then quotes, so a Windows path
## or an apostrophe in OUTDIR or SIM_ID reaches the model unchanged.
.eidos_quote <- function(s) paste0("'", gsub("'", "\\'", gsub("\\", "\\\\", s, fixed = TRUE), fixed = TRUE), "'")

#' Convert a parameter set into SLiM `-d` command-line arguments
#'
#' Logicals become `T`/`F`, strings are wrapped in Eidos single quotes (with quotes and
#' backslashes escaped, so Windows paths and apostrophes survive), and numbers are
#' written so that Eidos infers the type the model expects. Returned as a flat character
#' vector `c("-d","NAME=val", ...)` suitable for a shell-free process launch.
#'
#' @section Why float-typed parameters always carry a decimal point:
#' Eidos types `-d P=0` as an **integer** and does not promote it, so a probability set to
#' exactly 0 or 1 crashed the model mid-run in at least two different ways:
#' `rbinom()` rejecting an integer `prob` (reached by `P_RAD_SURVIVAL` with radical cure or
#' MDA on), and `min()` refusing mixed types (reached by `ALPHA_ACTIVATE = 0`). Both are
#' latent: the same parameter set completes if the code path
#' never fires, so a sweep fails in only some cells. "Perfect radical cure" and
#' "no relapses" are the two most natural settings in that space, and `validate_params()`
#' explicitly admits `[0, 1]`.
#'
#' The decision is keyed on the **default's** type, not the supplied value's, so
#' `N_HOSTS = 100` stays integral while `vivax_params(P_RAD_SURVIVAL = 0L)` — an integer
#' supplied by the user for a float parameter — is still emitted as `0.0`.
#'
#' @param params a `vivax_params` object (or named list).
#' @return character vector of `-d NAME=value` tokens.
#' @export
as_slim_args <- function(params) {
  ## 15 significant digits, not format()'s default 7: the old setting silently rounded
  ## ALPHA_ACTIVATE = 1/170 to 0.005882353, and params.tsv then recorded the rounded value.
  num <- .slim_num
  as_float <- function(x) {
    v <- num(x)
    if (grepl("^-?[0-9]+$", v)) paste0(v, ".0") else v
  }
  fmt <- function(name, val) {
    if (length(val) != 1L || is.na(val))
      stop(name, " must be a single non-missing value (got length ", length(val), ")", call. = FALSE)
    default <- .malariafwd_defaults[[name]]   # NULL for SEED/OUTDIR, added by run_sim()
    if (is.logical(val)) v <- if (isTRUE(val)) "T" else "F"
    else if (is.character(val) || name %in% .malariafwd_string_params) v <- .eidos_quote(val)
    else if (!is.null(default) && is.double(default)) v <- as_float(val)
    else v <- num(val)
    c("-d", paste0(name, "=", v))
  }
  unlist(Map(fmt, names(params), params), use.names = FALSE)
}

#' @export
print.vivax_params <- function(x, ...) {
  cat("<vivax_params>", length(x), "parameters\n")
  nd <- x[!vapply(names(x), function(n) identical(x[[n]], .malariafwd_defaults[[n]]), logical(1))]
  if (length(nd)) {
    cat("  non-default:\n")
    for (n in names(nd)) cat(sprintf("    %-24s = %s\n", n, format(nd[[n]])))
  } else cat("  (all defaults)\n")
  invisible(x)
}
