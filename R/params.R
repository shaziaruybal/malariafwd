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
  BITING_RATE = 0.1, P_TRANS_VH = 0.5, P_TRANS_HV = 0.5, MOSQ_LIFESPAN = 10L,
  EIP_DAYS = 10L, MAX_MOI = 20L, MAX_BROODS_PER_MOSQ = 20L,
  # parasite / infection
  MEAN_BLOOD_STAGE_DURATION = 60, MAX_BROOD_AGE = 180L, N_OOCYSTS_PER_BITE = 1L,
  # mosquito uptake / inoculation sampling (WS-F). lambda = 0 (default) -> exactly 1
  # oocyst per uptake and 1 genotype per bite. K_OOCYSTS > 0 makes the oocyst count
  # negative-binomial (var = lambda + lambda^2/K); 0 = Poisson.
  LAMBDA_OOCYSTS = 0.0, LAMBDA_SPOROZOITES = 0.0, K_OOCYSTS = 0.0,
  GAMETOCYTE_DELAY = 0L, P_GAMETOCYTE = 1.0, P_PRIMARY = 1.0,
  # hypnozoite reservoir
  NU_HYPNOZOITE_BATCH = 5.0, ALPHA_ACTIVATE = 1/170, MU_HYPNOZOITE = 1/442,
  K_LATENCY = 0L, DELTA_LATENCY = 1/100, MAX_HYPNOZOITE_AGE = 1500L,
  # genome (see DECISIONS D2/D3/D13): pure coordinate switch; both modes rate 7.4e-7
  GENOME_MODE = "full", CHROMOSOMES = "1-14", RECOMB_RATE = -1.0,
  # interventions (case-based radical cure = liver, blood-stage treatment = schizonticidal)
  INTERVENTION_START = 99999L, INTERVENTION_END = 99999L, INTERVENTION_EFFICACY = 0.0,
  MDA_TICK = 99999L, P_RAD_SURVIVAL = 0.05, RADCURE_COVERAGE = 0.0,
  BLOODCURE_COVERAGE = 0.0, P_BLOODCURE_SURVIVAL = 0.0,
  # detection / observation
  P_DETECT_PRIMARY = 1.0, P_DETECT_RELAPSE = 1.0,
  # seasonality
  SEASONAL_AMP = 0.0, SEASONAL_PHASE = 0L,
  # resurgence triggers
  RESURGENCE_TRIGGER_TICK = 99999L, BITING_RATE_BOOST = 1.0,
  ALPHA_ACTIVATE_BOOST = 1.0, NU_HYPNOZOITE_BATCH_BOOST = 1.0,
  # importation
  IMPORT_START = 99999L, IMPORT_END = 99999L, N_IMPORTS_PER_TICK = 0L,
  IMPORT_RATE_PER_TICK = 0.0, IMPORT_RATE_BOOST = 1.0, IMPORT_TICK = 99999L, N_IMPORTS = 0L
)

# Params that are Eidos STRINGS (must be quoted when passed to SLiM -d).
.malariafwd_string_params <- c("SIM_ID", "GENOME_MODE", "CHROMOSOMES", "OUTDIR")

#' Build a validated parameter set for the vivax model
#'
#' Returns a named list of all model parameters, starting from the defaults (which
#' mirror the SLiM model) and overriding any you name. `SEED` and `OUTDIR` are set by
#' [run_sim()], not here.
#'
#' @param ... named parameter overrides, e.g. `N_HOSTS = 300, BITING_RATE = 0.09`.
#' @return An object of class `vivax_params` (a named list).
#' @examples
#' vivax_params(N_HOSTS = 300, BITING_RATE = 0.09, GENOME_MODE = "core")
#' @export
vivax_params <- function(...) {
  over <- list(...)
  unknown <- setdiff(names(over), names(.malariafwd_defaults))
  if (length(unknown))
    stop("unknown parameter(s): ", paste(unknown, collapse = ", "),
         "\nsee names(vivax_params()) for valid names.", call. = FALSE)
  p <- utils::modifyList(.malariafwd_defaults, over)
  # v3 backward-compat: if the founder decoupling isn't explicitly set, mirror the legacy N_FOUNDERS
  if (!"N_FOUNDER_HOSTS"    %in% names(over)) p$N_FOUNDER_HOSTS    <- p$N_FOUNDERS
  if (!"N_FOUNDER_LINEAGES" %in% names(over)) p$N_FOUNDER_LINEAGES <- p$N_FOUNDER_HOSTS
  validate_params(p)
  structure(p, class = "vivax_params")
}

#' @keywords internal
validate_params <- function(p) {
  stopifnot(p$N_FOUNDERS <= p$N_HOSTS,
            p$N_FOUNDER_LINEAGES <= p$N_FOUNDER_HOSTS, p$N_FOUNDER_HOSTS <= p$N_HOSTS)
  if (!p$GENOME_MODE %in% c("full", "core"))
    stop("GENOME_MODE must be 'full' or 'core', got '", p$GENOME_MODE, "'", call. = FALSE)
  for (pr in c("P_TRANS_VH", "P_TRANS_HV", "P_GAMETOCYTE", "P_PRIMARY",
               "P_DETECT_PRIMARY", "P_DETECT_RELAPSE", "RADCURE_COVERAGE", "P_RAD_SURVIVAL",
               "BLOODCURE_COVERAGE", "P_BLOODCURE_SURVIVAL"))
    if (p[[pr]] < 0 || p[[pr]] > 1) stop(pr, " must be in [0, 1]", call. = FALSE)
  for (pr in c("SIM_DAYS", "N_HOSTS", "N_MOSQ", "BITING_RATE", "MOSQ_LIFESPAN"))
    if (p[[pr]] <= 0) stop(pr, " must be > 0", call. = FALSE)
  for (pr in c("LAMBDA_OOCYSTS", "LAMBDA_SPOROZOITES", "K_OOCYSTS"))
    if (p[[pr]] < 0) stop(pr, " must be >= 0", call. = FALSE)
  if (p$N_OOCYSTS_PER_BITE != 1)
    stop("N_OOCYSTS_PER_BITE is deprecated; use LAMBDA_SPOROZOITES instead.", call. = FALSE)

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
  invisible(TRUE)
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

#' Convert a parameter set into SLiM `-d` command-line arguments
#'
#' Numbers pass through as-is, logicals become `T`/`F`, and strings are wrapped in
#' Eidos single quotes. Returned as a flat character vector `c("-d","NAME=val", ...)`
#' suitable for a shell-free process launch (no quote-stripping).
#' @param params a `vivax_params` object (or named list).
#' @return character vector of `-d NAME=value` tokens.
#' @export
as_slim_args <- function(params) {
  fmt <- function(name, val) {
    if (is.logical(val)) v <- if (isTRUE(val)) "T" else "F"
    else if (is.character(val) || name %in% .malariafwd_string_params) v <- sprintf("'%s'", val)
    else v <- format(val, scientific = FALSE, trim = TRUE)
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
