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
#' `TRUE`/`FALSE` or a string, as its default is; numbers are finite, and strings are UTF-8 text
#' without control characters. Counts, `SIM_DAYS`,
#' the scheduling days (`MDA_TICK`, `IMPORT_TICK`, `RESURGENCE_TRIGGER_TICK` and the `*_START`
#' and `*_END` bounds) and the fixed delays `EIP_DAYS` and `GAMETOCYTE_DELAY` are whole numbers;
#' mean durations, rates and `SEASONAL_PHASE` need not be. Probabilities and coverages,
#' `ALPHA_ACTIVATE` and `INTERVENTION_EFFICACY` lie in \[0, 1\], `MU_HYPNOZOITE` in (0, 1\], and
#' `DELTA_LATENCY` in (0, 1\] when `K_LATENCY > 0`. `BITING_RATE` is positive (above 1 it only
#' warns). Boosts, import rates, `LAMBDA_*`, `SEASONAL_AMP` and counts are non-negative, and
#' `MAX_MOI` is at least 1 (1 means no superinfection).
#' `RECOMB_RATE` is at most 0.5 per bp (any value <= 0 means the default map rate), and run time
#' grows steeply with it.
#'
#' A few sizes are capped far above anything realistic: a mean of 1000 oocysts per successful
#' uptake (see [uptake_moments()]), with `K_OOCYSTS` either 0 or at least 0.001;
#' `LAMBDA_SPOROZOITES` at 1000; the hypnozoite batch mean at 1e6, times its boost once the
#' trigger falls inside the run; the import rate at 1e6 a day; and `TREATMENT_DELAY` and
#' `MEAN_RECRUD_DELAY` at 1e6 days. `SIM_DAYS` is at most 1e9, SLiM's last tick, and a run that
#' reaches day 99999, the day the defaults use for "never", must move those days past its end and
#' cannot schedule anything on day 99999 itself. A founder-free run whose imports never establish
#' a parasite, by chance or because its settings let none take, is not refused here; it ends as
#' never started (see [load_run()]). Imports arrive on days 1 to `SIM_DAYS`: a pulse or an import
#' window set outside the run brings nothing, without a warning, which is how imports are switched
#' off in one arm of a grid. Each host takes at most one
#' import a day, so a day's imports beyond `N_HOSTS` never arrive. This warns when a day's fixed
#' imports inside the run exceed `N_HOSTS` (`N_IMPORTS` on `IMPORT_TICK`, `N_IMPORTS_PER_TICK` in
#' the window, or both on the same day), and when `IMPORT_RATE_PER_TICK` at its peak over the
#' window, with its seasonal, intervention and resurgence multipliers and with
#' `N_IMPORTS_PER_TICK`, draws more than `N_HOSTS` imports on a day with probability above 0.01.
#'
#' A named number (such as `coef(fit)["k"]`), a value wrapped in `I()`, a whole number carrying
#' floating-point noise (such as `1.1 * 100`, or `0.1 * 3 - 0.3` for 0), a number or factor given
#' for a string parameter, latin1 text, and spaces around the commas and dashes of `CHROMOSOMES` are
#' stored in the form SLiM is given. A number with any other class, such as a lubridate duration or
#' a `bit64::integer64`, is refused: give the number itself.
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
         if (any(c("SEED", "OUTDIR") %in% unknown))
           "\nSEED and OUTDIR are set by run_sim(seed =, outdir =), not in the parameter set.",
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
## params.tsv and the run agree: a named number (coef(fit)["k"]) or a 1 x 1 matrix, a whole number
## with floating-point noise (1.1 * 100 is 110.00000000000001, written as 110; 0.1 * 3 - 0.3 is
## 5.6e-17, written as 0), a number or factor for a string parameter (expand.grid() makes
## factors), text in another encoding (as UTF-8, which SLiM reads), and spaces around the commas
## and dashes of CHROMOSOMES. Spaces anywhere else are left for the check to refuse: "1 4" is not
## "14".
.normalise_params <- function(p) {
  d <- .malariafwd_defaults
  for (nm in intersect(names(p), names(d))) {
    v <- p[[nm]]
    ## I() only asks for a value to be taken as it is, so its class is dropped; oldClass() keeps
    ## I(factor(x)) a factor
    if (inherits(v, "AsIs")) v <- p[[nm]] <- `oldClass<-`(v, setdiff(oldClass(v), "AsIs"))
    if (is.atomic(v) && !is.object(v) && length(v)) v <- p[[nm]] <- as.vector(v)
    ## functions, symbols and calls are left for the type check to refuse
    if (!(is.atomic(v) || is.list(v)) || length(v) != 1L || is.na(v)) next
    ## a classed value (a Date, a duration, units) is left for the type check to refuse by name
    if (is.integer(d[[nm]]) && is.double(v) && !is.object(v) && is.finite(v) &&
        (abs(v) < 1e-9 || grepl("^-?[0-9]+$", .slim_num(v))))
      p[[nm]] <- round(v)
    else if (is.character(d[[nm]]) && is.factor(v))
      p[[nm]] <- as.character(v)
    else if (is.character(d[[nm]]) && is.numeric(v))
      p[[nm]] <- .slim_num(v)
    ## Text as UTF-8, after the factor step (expand.grid() of latin1 text gives a latin1 factor):
    ## latin1, which R reads as CP1252, as enc2utf8() converts it; UTF-8 marked as bytes; and text in
    ## a non-UTF-8 native encoding, through iconv(), which gives NA where it cannot convert (as in a
    ## C locale) and so leaves the text for the UTF-8 check to refuse.
    v <- p[[nm]]
    if (is.character(v) && length(v) == 1L && !is.na(v)) {
      enc <- Encoding(v)
      if (enc == "latin1") p[[nm]] <- enc2utf8(v)
      else if (enc == "bytes" && validUTF8(v)) p[[nm]] <- `Encoding<-`(v, "UTF-8")
      else if (enc == "unknown" && !validUTF8(v)) {
        u <- iconv(v, "", "UTF-8")
        if (!is.na(u)) p[[nm]] <- u
      }
    }
  }
  ch <- p$CHROMOSOMES
  if (is.character(ch) && length(ch) == 1L && !is.na(ch) &&
      grepl("^ *[0-9]+ *(- *[0-9]+ *)?(, *[0-9]+ *(- *[0-9]+ *)?)*$", ch, useBytes = TRUE))
    p$CHROMOSOMES <- gsub(" ", "", ch, fixed = TRUE)
  p
}

## The days that default to 99999, which the model reads as "never".
.never_params <- names(.malariafwd_defaults)[vapply(.malariafwd_defaults, identical, logical(1), 99999L)]

## The largest daily value of a rate over days first..last, worked out the way the model applies
## its factors to the biting and import rates: seasonality only when SEASONAL_AMP > 0 (floored at
## 0), the intervention factor inside [INTERVENTION_START, INTERVENTION_END), the boost from
## RESURGENCE_TRIGGER_TICK. The factors change only on those days and the seasonal one repeats every
## 365 days, so a year from each change covers every value the rate can take.
.peak_rate <- function(p, rate, boost, first, last) {
  starts <- c(first, p$INTERVENTION_START, p$INTERVENTION_END, p$RESURGENCE_TRIGGER_TICK)
  starts <- starts[starts >= first & starts <= last]
  t <- sort(unique(unlist(lapply(starts, function(a) seq(a, min(last, a + 364))))))
  r <- rep(rate, length(t))
  if (p$SEASONAL_AMP > 0) {
    angle <- 2 * pi * (t - p$SEASONAL_PHASE) / 365
    if (!all(is.finite(angle))) return(c(rate = NaN, day = NA))
    r <- r * pmax(0, 1 + p$SEASONAL_AMP * cos(angle))
  }
  inside <- t >= p$INTERVENTION_START & t < p$INTERVENTION_END
  r[inside] <- r[inside] * (1 - p$INTERVENTION_EFFICACY)
  boosted <- t >= p$RESURGENCE_TRIGGER_TICK
  r[boosted] <- r[boosted] * boost
  c(rate = max(r), day = t[which.max(r)])
}

## The largest daily biting probability the parameters ask of a mosquito. Biting starts on tick 2,
## once the populations exist.
.peak_biting <- function(p) .peak_rate(p, p$BITING_RATE, p$BITING_RATE_BOOST, 2, max(2, p$SIM_DAYS))

#' @keywords internal
validate_params <- function(p) {
  ## functions, symbols and calls are left for the type check, which names what was given
  bad <- names(p)[!vapply(p, function(v) !(is.atomic(v) || is.list(v)) ||
                            (length(v) == 1L && !is.na(v)), logical(1))]
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
               "GAMETOCYTE_DELAY", "MAX_BROODS_PER_MOSQ", "K_LATENCY"))
    if (p[[pr]] < 0) stop(pr, " must be >= 0", call. = FALSE)
  ## at 0 the gate would refuse every inoculation, an empty host's included (D38)
  if (p$MAX_MOI < 1)
    stop("MAX_MOI must be at least 1, got ", .slim_num(p$MAX_MOI), ": a host takes a new inoculation ",
         "only while it carries fewer than MAX_MOI blood-stage lineages, so 1 means no superinfection.",
         call. = FALSE)
  ## N_FOUNDERS itself is only the default of N_FOUNDER_HOSTS, so it is not checked
  if (p$N_FOUNDER_HOSTS > p$N_HOSTS)
    stop("N_FOUNDER_HOSTS (", .slim_num(p$N_FOUNDER_HOSTS), ") must be at most N_HOSTS (",
         .slim_num(p$N_HOSTS), "). It defaults to N_FOUNDERS, so set either one.", call. = FALSE)
  if (p$N_FOUNDER_LINEAGES > p$N_FOUNDER_HOSTS)
    stop("N_FOUNDER_LINEAGES (", .slim_num(p$N_FOUNDER_LINEAGES), ") must be at most N_FOUNDER_HOSTS (",
         .slim_num(p$N_FOUNDER_HOSTS), "): each founder host carries one of the founding genomes.",
         call. = FALSE)
  if (p$SEASONAL_AMP < 0)
    stop("SEASONAL_AMP must be >= 0 (0 means no seasonality). SEASONAL_PHASE is the day of the ",
         "peak, so for a trough on day x use SEASONAL_PHASE = x + 182.5.", call. = FALSE)
  if (p$SIM_DAYS > 1e9)
    stop("SIM_DAYS must be at most 1e9, the last tick SLiM can schedule.", call. = FALSE)
  never <- .never_params[vapply(.never_params, function(nm) p[[nm]] == 99999, logical(1))]
  if (p$SIM_DAYS >= 99999 && length(never))
    stop("SIM_DAYS = ", .slim_num(p$SIM_DAYS), " reaches day 99999, which the defaults use for ",
         "'never', and ", paste(never, collapse = ", "), if (length(never) > 1) " are" else " is",
         " still 99999: set each to a day after SIM_DAYS, or to the day you mean (in a run this ",
         "long, day 99999 itself cannot be used).", call. = FALSE)
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
    stop("nothing can start this run: N_FOUNDER_HOSTS is 0 and no importation is configured\n",
         "  inside the run. Set N_FOUNDER_HOSTS >= 1, or open an import path on a day from 1 to\n",
         "  SIM_DAYS (IMPORT_TICK + N_IMPORTS, or IMPORT_START/IMPORT_END with N_IMPORTS_PER_TICK\n",
         "  or IMPORT_RATE_PER_TICK).", call. = FALSE)

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
  ## Each host takes at most one import a day (D38), so a day's imports beyond N_HOSTS never
  ## arrive. The fixed counts of the run's days are checked: the pulse day (with the window's count
  ## when it falls in the window) and the window's other days, each named when it is over. A
  ## Poisson rate draws above N_HOSTS by chance, rarely while it is well below N_HOSTS but most days
  ## near or above it, so its peak over the window, with the window's count, is checked too.
  days <- .import_days(p); win <- days$window
  daily <- if (!is.null(win)) p$N_IMPORTS_PER_TICK else 0
  in_win <- !is.null(days$pulse) && !is.null(win) && days$pulse >= win[1] && days$pulse <= win[2]
  pulse <- if (!is.null(days$pulse)) p$N_IMPORTS + (if (in_win) daily else 0) else 0
  over <- c(if (pulse > p$N_HOSTS)
              paste0(if (pulse > p$N_IMPORTS) "N_IMPORTS and N_IMPORTS_PER_TICK bring " else "N_IMPORTS brings ",
                     .slim_num(pulse), " imports on day ", .slim_num(p$IMPORT_TICK)),
            if (daily > p$N_HOSTS && win[2] - win[1] + 1 > in_win)
              paste0("N_IMPORTS_PER_TICK brings ", .slim_num(daily), " imports a day"))
  if (length(over))
    warning(paste(over, collapse = ", and "), ", but each host takes at most one import a day, so ",
            "at most N_HOSTS = ", .slim_num(p$N_HOSTS), " of them arrive.", call. = FALSE)
  room <- p$N_HOSTS - daily
  if (!is.null(win) && p$IMPORT_RATE_PER_TICK > 0 && room >= 0) {
    peak <- .peak_rate(p, p$IMPORT_RATE_PER_TICK, p$IMPORT_RATE_BOOST, win[1], win[2])
    lost <- if (isTRUE(peak[["rate"]] > 0)) stats::ppois(room, peak[["rate"]], lower.tail = FALSE) else 0
    if (lost > 0.01)
      warning("IMPORT_RATE_PER_TICK with its seasonal, intervention and resurgence multipliers reaches ",
              signif(peak[["rate"]], 3), " imports a day on day ", .slim_num(peak[["day"]]),
              ": that day draws more than N_HOSTS = ", .slim_num(p$N_HOSTS), " imports",
              if (daily > 0) paste0(" together with N_IMPORTS_PER_TICK = ", .slim_num(daily)),
              " with probability ", signif(lost, 3), ", but each host takes at most one import a ",
              "day, so the excess never arrives.", call. = FALSE)
  }
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
## differently: a wrong type stops the run at its first use, or runs with another meaning (a string
## P_TRANS_VH is compared as text). Strings must be valid UTF-8, which SLiM and the tree-sequence
## provenance it writes require. Numbers must be finite and
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
    dur <- if (inherits(v, c("difftime", "Duration", "Period"))) "; for a duration give the number of days" else ""
    if (!ok) stop(nm, " must be ", want, ", got ",
                  if (is.object(v)) paste("an object of class", class(v)[1]) else deparse(v), dur, call. = FALSE)
    ## a number with a class need not mean its bare value (a lubridate duration counts seconds, and
    ## bit64's integer64 is not a double at all), so only plain numbers are taken; I() is dropped
    ## before this
    if (is.numeric(v) && is.object(v))
      stop(nm, " must be a plain number, got an object of class ", class(v)[1],
           if (nzchar(dur)) dur else "; give the number itself, for example with as.numeric()", ".",
           call. = FALSE)
    if (is.character(v) && !validUTF8(v))
      stop(nm, " must be valid UTF-8 text; for Windows or latin1 text, mark it with ",
           "Encoding(x) <- \"latin1\" (R reads it as Windows-1252) or convert it with ",
           "iconv(x, from = \"CP1252\", to = \"UTF-8\"), or read the file it came from with its ",
           "encoding declared.", call. = FALSE)
    ## byte-wise, so that accented and other non-ASCII labels pass in every locale: the C0
    ## controls and DEL, and the C1 controls U+0080-U+009F, which are C2 80 to C2 9F in UTF-8
    if (is.character(v) && grepl("[\001-\037\177]|\302[\200-\237]", v, useBytes = TRUE))
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
## below where SLiM's draws fail (binomial sizes are 32-bit; Poisson draws wrap above 2^32). A
## value over a cap only by rounding (a fitted or computed value a few ulps above it) is accepted,
## so a refused value prints visibly above the cap. NaN is over every cap.
.over_cap <- function(x, cap) !isTRUE(x <= cap * (1 + 1e-12))
.check_sizes <- function(p) {
  if (.over_cap(p$LAMBDA_SPOROZOITES, 1000))
    stop("LAMBDA_SPOROZOITES must be <= 1000, got ", format(p$LAMBDA_SPOROZOITES, digits = 15),
         ": it is the pre-truncation rate of oocyst products sampled per infectious bite.",
         call. = FALSE)
  ## K = 0.001 already comes within 0.5% of the largest variance the zero-truncated negative
  ## binomial reaches at its mean (the K -> 0 limit), for any mean up to 1000, while at a given
  ## mean the model's gamma-Poisson tries per uptake grow about as 1 / K below it.
  if (p$K_OOCYSTS > 0 && p$K_OOCYSTS < 1e-3 * (1 - 1e-12))
    stop("K_OOCYSTS must be 0 (Poisson) or at least 0.001, got ", format(p$K_OOCYSTS, digits = 15),
         ": at 0.001 the count already comes within 0.5% of the largest variance a zero-truncated ",
         "negative binomial reaches at its mean, and a smaller K makes each uptake take ",
         "proportionally more draws.", call. = FALSE)
  ## A cap on the mean, not on LAMBDA_OOCYSTS alone: a small K gives a heavy tail whose mean, and
  ## so the work per uptake, can be many times lambda.
  m <- uptake_moments(p$LAMBDA_OOCYSTS, p$K_OOCYSTS)[["mean"]]
  if (.over_cap(m, 1000))
    stop("LAMBDA_OOCYSTS = ", p$LAMBDA_OOCYSTS, " with K_OOCYSTS = ", p$K_OOCYSTS, " gives a mean ",
         "of ", format(m, digits = 15), " oocysts per successful uptake (see uptake_moments()); the ",
         "mean must be at most 1000.", call. = FALSE)
  ## a boost below 1 lowers the rate only after the trigger, so the peak is the unboosted rate
  boost <- function(b) if (p$RESURGENCE_TRIGGER_TICK <= p$SIM_DAYS) max(1, b) else 1
  nu <- p$NU_HYPNOZOITE_BATCH * boost(p$NU_HYPNOZOITE_BATCH_BOOST)
  if (.over_cap(nu, 1e6))
    stop("NU_HYPNOZOITE_BATCH, times NU_HYPNOZOITE_BATCH_BOOST (if above 1) once the trigger falls ",
         "inside the run, is ", format(nu, digits = 15), " hypnozoites per bite; it must be at most ",
         "1e6, far below the 32-bit sizes of SLiM's binomial draws.", call. = FALSE)
  peak <- p$IMPORT_RATE_PER_TICK * boost(p$IMPORT_RATE_BOOST) * (1 + p$SEASONAL_AMP)
  if (.over_cap(peak, 1e6))
    stop("IMPORT_RATE_PER_TICK, times IMPORT_RATE_BOOST (if above 1) once the trigger falls inside ",
         "the run and 1 + SEASONAL_AMP at the seasonal peak, reaches ", format(peak, digits = 15),
         " imports a day; it must be at most 1e6 imports a day.", call. = FALSE)
}

## The model accepts "1-14", "1,3,5" or "5-8": comma-separated chromosome numbers or low-high
## ranges, each chromosome once. The ends are checked before a range is expanded.
.check_chromosomes <- function(x) {
  bad <- function(why) stop("CHROMOSOMES must be chromosome numbers 1-14, comma-separated, with ",
                            "low-high ranges, e.g. \"1-14\" or \"1-3,5\" (", why, "); got \"", x,
                            "\".", call. = FALSE)
  if (!grepl("^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$", x, useBytes = TRUE)) bad("unreadable")
  idx <- unlist(lapply(strsplit(strsplit(x, ",")[[1]], "-"), function(r) {
    r <- as.numeric(r)
    if (any(r < 1 | r > 14)) bad("out of range")
    if (length(r) == 2 && r[1] > r[2]) bad("a range runs high-low")
    seq(r[1], r[length(r)])
  }))
  if (anyDuplicated(idx)) bad("a chromosome is listed twice")
}

## The days of the run on which the parameters bring imports, as the model reads them: the pulse
## day, when N_IMPORTS > 0 and IMPORT_TICK is a day from 1 to SIM_DAYS, and the first and last day
## of the window inside the run, max(1, IMPORT_START) to min(IMPORT_END - 1, SIM_DAYS), when it has
## any. Imports arrive in early(), which first runs on day 1. NULL for a pulse or window that brings
## nothing; whether the window brings anything also depends on its count and rate, which the callers
## check. The model's LAST_IMPORT_DAY uses the same days.
.import_days <- function(p) {
  first <- max(1, p$IMPORT_START); last <- min(p$IMPORT_END - 1, p$SIM_DAYS)
  list(pulse = if (p$N_IMPORTS > 0 && p$IMPORT_TICK >= 1 && p$IMPORT_TICK <= p$SIM_DAYS) p$IMPORT_TICK,
       window = if (first <= last) c(first, last))
}

#' Is any importation path open in this parameter set?
#'
#' Either the one-shot pulse (`IMPORT_TICK` + `N_IMPORTS`) or the windowed rate
#' (`IMPORT_START`/`IMPORT_END` with `N_IMPORTS_PER_TICK` or `IMPORT_RATE_PER_TICK`), on a day of
#' the run: imports arrive from day 1, so the pulse needs `IMPORT_TICK` from 1 to `SIM_DAYS`, and
#' the window needs a day from `max(1, IMPORT_START)` to `min(IMPORT_END - 1, SIM_DAYS)`.
#' @param p a parameter list.
#' @return `TRUE` if at least one importation can occur.
#' @keywords internal
has_import_path <- function(p) {
  d <- .import_days(p)
  !is.null(d$pulse) ||
    (!is.null(d$window) && (p$N_IMPORTS_PER_TICK > 0 || p$IMPORT_RATE_PER_TICK > 0))
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
