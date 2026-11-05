test_that("vivax_params defaults stay in sync with the SLiM model", {
  model <- slim_model()
  skip_if(!nzchar(model) || !file.exists(model), "model template not found")
  src <- readLines(model, warn = FALSE)
  model_names <- unique(sub('.*defineDefault\\("([A-Z_]+)".*', "\\1",
                            grep('defineDefault\\("[A-Z_]+"', src, value = TRUE)))
  r_names <- names(malariafwd:::.malariafwd_defaults)
  # SEED and OUTDIR are run-level (set by run_sim), not part of vivax_params().
  expect_setequal(setdiff(model_names, c("SEED", "OUTDIR")), r_names)

  ## ...and the VALUES, not just the names. Comparing names alone let the two drift apart
  ## silently: the check passed with P_RAD_SURVIVAL perturbed 0.29 -> 0.05 on one side, which
  ## is exactly the divergence D20 was written to prevent (the same number means 77% efficacy
  ## per hypnozoite but 95% per person).
  dflt <- malariafwd:::.malariafwd_defaults
  num <- function(x) suppressWarnings(as.numeric(x))
  for (nm in setdiff(model_names, c("SEED", "OUTDIR"))) {
    line <- grep(paste0('defineDefault\\("', nm, '"'), src, value = TRUE)[1]
    raw  <- trimws(sub('.*defineDefault\\("[A-Z_]+",\\s*([^)]*)\\).*', "\\1", line))
    rval <- dflt[[nm]]
    sval <- if (raw %in% c("T", "F")) raw == "T"
            else if (grepl("^'.*'$|^\".*\"$", raw)) gsub("^['\"]|['\"]$", "", raw)
            else num(raw)
    if (is.na(sval[1]) && !is.logical(sval)) next          # expression, not a literal
    if (is.logical(rval) || is.character(rval))
      expect_identical(as.character(sval), as.character(rval),
                       info = paste(nm, "differs between the SLiM model and R"))
    else
      expect_equal(num(sval), num(rval), tolerance = 1e-9,
                   info = paste(nm, "differs between the SLiM model and R"))
  }
})

test_that("validation catches bad input", {
  expect_error(vivax_params(GENOME_MODE = "xyz"))
  expect_error(vivax_params(N_HOSTS = 5, N_FOUNDERS = 10))
  expect_error(vivax_params(FOO = 1))
  expect_error(vivax_params(P_TRANS_VH = 2))
  expect_error(vivax_params(SIM_DAYS = -1))
})

test_that("as_slim_args quotes strings, formats logicals, passes numbers", {
  a <- as_slim_args(vivax_params(GENOME_MODE = "core", VERBOSE = TRUE, N_HOSTS = 200))
  expect_true("GENOME_MODE='core'" %in% a)
  expect_true("VERBOSE=T" %in% a)
  expect_true("RETAIN_FULL_HISTORY=T" %in% a)
  expect_true("N_HOSTS=200" %in% a)
  expect_true(all(a[seq(1, length(a), 2)] == "-d"))  # every other token is the -d flag
})

test_that("print.vivax_params reports non-defaults", {
  expect_output(print(vivax_params(N_HOSTS = 300)), "N_HOSTS")
})

test_that("float parameters always carry a decimal point, integers never do", {
  ## Eidos types `-d P=0` as an integer and will not promote it, which crashed the model
  ## mid-run: rbinom() rejects an integer `prob`, and min() refuses mixed argument types.
  a <- as_slim_args(vivax_params(P_RAD_SURVIVAL = 0, BLOODCURE_COVERAGE = 1,
                                 ALPHA_ACTIVATE = 0, N_HOSTS = 100))
  expect_true("P_RAD_SURVIVAL=0.0" %in% a)
  expect_true("BLOODCURE_COVERAGE=1.0" %in% a)
  expect_true("ALPHA_ACTIVATE=0.0" %in% a)
  ## integer parameters must NOT gain one, or Eidos rejects them where it wants an integer
  expect_true("N_HOSTS=100" %in% a)
  expect_true("SIM_DAYS=365" %in% a)
  ## keyed on the DEFAULT's type, so an integer supplied for a float parameter still converts
  expect_true("P_RAD_SURVIVAL=0.0" %in% as_slim_args(vivax_params(P_RAD_SURVIVAL = 0L)))
})

test_that("as_slim_args does not round float parameters", {
  ## format()'s default of 7 significant digits silently truncated ALPHA_ACTIVATE = 1/170
  ## to 0.005882353, and params.tsv then recorded the rounded value as provenance.
  a <- as_slim_args(vivax_params(ALPHA_ACTIVATE = 1 / 170))
  v <- sub("^ALPHA_ACTIVATE=", "", grep("^ALPHA_ACTIVATE=", a, value = TRUE))
  expect_equal(as.numeric(v), 1 / 170, tolerance = 1e-12)
})

test_that("probabilities of exactly 0 and 1 run to completion", {
  skip_if_no_slim()
  ## Both crash paths. They are LATENT: the same parameter set completes if the code path
  ## never fires, so a sweep failed in only some cells. rbinom() needs radical cure or MDA
  ## to actually fire; min() is reached on any tick with a hypnozoite.
  base <- list(SIM_DAYS = 60L, N_HOSTS = 20L, N_MOSQ = 40L)
  cases <- list(
    rbinom_zero = c(base, list(P_RAD_SURVIVAL = 0, MDA_TICK = 30L)),
    rbinom_one  = c(base, list(P_RAD_SURVIVAL = 1, MDA_TICK = 30L)),
    min_alpha   = c(base, list(ALPHA_ACTIVATE = 0))
  )
  for (nm in names(cases)) {
    out <- file.path(tempdir(), paste0("mfw-zero-one-", nm))
    unlink(out, recursive = TRUE)
    expect_no_error(run_sim(do.call(vivax_params, cases[[nm]]), seed = 1, outdir = out))
  }
})

test_that("as_slim_args() escapes quotes and backslashes in string values", {
  a <- as_slim_args(list(SIM_ID = "Bob's sweep", OUTDIR = "C:\\temp\\new dir"))
  expect_equal(a, c("-d", "SIM_ID='Bob\\'s sweep'", "-d", "OUTDIR='C:\\\\temp\\\\new dir'"))
})

test_that("vector or missing parameter values are refused before SLiM is launched", {
  expect_error(vivax_params(EIP_DAYS = c(5, 10)), "EIP_DAYS must be a single non-missing value")
  expect_error(vivax_params(N_HOSTS = NA), "N_HOSTS must be a single non-missing value")
  expect_error(as_slim_args(list(EIP_DAYS = c(5, 10))), "EIP_DAYS must be a single non-missing value")
})

test_that("run_sim() refuses a seed that is not a single whole number", {
  expect_error(run_sim(seed = 1:3, outdir = tempfile()), "`seed` must be a single whole number")
  expect_error(run_sim(seed = 2.5, outdir = tempfile()), "`seed` must be a single whole number")
  expect_error(run_sim(seed = 2^31, outdir = tempfile()), "`seed` must be a single whole number")
})

test_that("a biting rate that would exceed one bite per mosquito per day warns", {
  expect_warning(vivax_params(BITING_RATE = 0.8, SEASONAL_AMP = 0.5), "bites at most once a day")
  expect_no_warning(vivax_params())
})

test_that("blood-stage infections have no age cap, and their mean duration is at least a day", {
  expect_error(vivax_params(MAX_BROOD_AGE = 180), "unknown parameter")
  expect_error(vivax_params(MEAN_BLOOD_STAGE_DURATION = 0.5), "MEAN_BLOOD_STAGE_DURATION must be >= 1")
})

test_that("hypnozoite batches have no age cap, and hypnozoites must be able to die", {
  expect_error(vivax_params(MAX_HYPNOZOITE_AGE = 1500), "unknown parameter")
  expect_error(vivax_params(MU_HYPNOZOITE = 0), "MU_HYPNOZOITE must be > 0")
})

test_that("treatment and recrudescence delays are at least a day, and finite", {
  expect_error(vivax_params(TREATMENT_DELAY = 0), "TREATMENT_DELAY must be >= 1")
  expect_error(vivax_params(TREATMENT_DELAY = 0.5), "TREATMENT_DELAY must be >= 1")
  expect_error(vivax_params(TREATMENT_DELAY = Inf), "TREATMENT_DELAY must be .*finite")
  expect_error(vivax_params(MEAN_RECRUD_DELAY = 0.5), "MEAN_RECRUD_DELAY must be >= 1")
  expect_error(vivax_params(MEAN_RECRUD_DELAY = Inf), "MEAN_RECRUD_DELAY must be .*finite")
  expect_no_error(vivax_params(TREATMENT_DELAY = 1, MEAN_RECRUD_DELAY = 1))
})

test_that("each value has its default's type, is finite, and is whole where the default is", {
  expect_error(vivax_params(ALPHA_ACTIVATE = FALSE), "ALPHA_ACTIVATE must be a number")
  expect_error(vivax_params(P_RAD_SURVIVAL = "0.29"), "P_RAD_SURVIVAL must be a number")
  expect_error(vivax_params(RETAIN_FULL_HISTORY = 1), "RETAIN_FULL_HISTORY must be TRUE or FALSE")
  expect_error(vivax_params(VERBOSE = "FALSE"), "VERBOSE must be TRUE or FALSE")
  expect_error(vivax_params(GENOME_MODE = TRUE), "GENOME_MODE must be a string")
  expect_error(vivax_params(MDA_TICK = Inf), "MDA_TICK must be finite.*use a day after SIM_DAYS")
  expect_error(vivax_params(TREATMENT_START = -Inf), "TREATMENT_START must be finite.*use 0")
  expect_error(vivax_params(K_OOCYSTS = Inf), "K_OOCYSTS must be finite.*use 0 for Poisson")
  expect_error(vivax_params(IMPORT_RATE_PER_TICK = Inf), "IMPORT_RATE_PER_TICK must be finite, got Inf$")
  expect_error(vivax_params(MU_HYPNOZOITE = 1e-310), "MU_HYPNOZOITE = .* is too small")
  ## SLiM cannot read the smallest normal double as as_slim_args() writes it
  expect_error(vivax_params(MU_HYPNOZOITE = .Machine$double.xmin), "is too small")
  expect_no_error(vivax_params(MU_HYPNOZOITE = 1e-300))
  expect_error(vivax_params(N_IMPORTS = 2.5), "N_IMPORTS must be a whole number, got 2.5; use round")
  expect_error(vivax_params(N_HOSTS = 100.5), "N_HOSTS must be a whole number")
  expect_error(vivax_params(N_HOSTS = 110.0001), "N_HOSTS must be a whole number, got 110.0001")
  expect_error(vivax_params(MDA_TICK = 150.5), "MDA_TICK must be a whole number")
  expect_error(vivax_params(MDA_TICK = 1e19), "MDA_TICK must lie within \\+/- 2\\^53")
  ## a mean duration and a day of the year need not be whole
  expect_no_error(vivax_params(MOSQ_LIFESPAN = 10.5, SEASONAL_PHASE = 40.5))
})

test_that("spellings of the same value are stored in the form SLiM is given", {
  ## a whole number with floating-point noise is whole as SLiM sees it
  p <- vivax_params(N_HOSTS = 100L, N_MOSQ = 1.1 * 100, N_FOUNDER_HOSTS = seq(0.1, 0.5, 0.1)[3] * 100,
                    SIM_DAYS = 365 * 2.2)
  expect_identical(c(p$N_MOSQ, p$N_FOUNDER_HOSTS, p$SIM_DAYS), c(110, 30, 803))
  expect_true("N_MOSQ=110" %in% as_slim_args(p))
  ## a number or a factor for a string parameter, as expand.grid() makes
  expect_identical(vivax_params(SIM_ID = 1)$SIM_ID, "1")
  expect_identical(vivax_params(SIM_ID = 100000)$SIM_ID, "100000")
  expect_identical(vivax_params(CHROMOSOMES = 14L)$CHROMOSOMES, "14")
  g <- expand.grid(RADCURE_MODE = c("person", "hypnozoite"), RADCURE_COVERAGE = 0.5,
                   BLOODCURE_COVERAGE = 1)
  modes <- vapply(seq_len(nrow(g)), function(i) do.call(vivax_params, as.list(g[i, ]))$RADCURE_MODE,
                  character(1))
  expect_identical(modes, c("person", "hypnozoite"))
  expect_error(vivax_params(GENOME_MODE = 1), "GENOME_MODE must be 'full' or 'core'")
  ## spaces around the commas and dashes of CHROMOSOMES, but not between digits
  expect_identical(vivax_params(CHROMOSOMES = " 1 - 3, 5 ")$CHROMOSOMES, "1-3,5")
  expect_error(vivax_params(CHROMOSOMES = "1 4"), "CHROMOSOMES.*unreadable")
})

test_that("rates and probabilities stay in range, and sizes stay within what SLiM can draw", {
  expect_error(vivax_params(ALPHA_ACTIVATE = -0.01), "ALPHA_ACTIVATE must be in \\[0, 1\\]")
  expect_error(vivax_params(ALPHA_ACTIVATE = 1.5), "ALPHA_ACTIVATE must be in \\[0, 1\\]")
  expect_error(vivax_params(ALPHA_ACTIVATE_BOOST = -1), "ALPHA_ACTIVATE_BOOST must be >= 0")
  expect_error(vivax_params(MU_HYPNOZOITE = 1.5), "MU_HYPNOZOITE must be <= 1")
  expect_error(vivax_params(INTERVENTION_EFFICACY = 1.2), "INTERVENTION_EFFICACY must be in \\[0, 1\\]")
  expect_error(vivax_params(IMPORT_RATE_BOOST = -1), "IMPORT_RATE_BOOST must be >= 0")
  expect_error(vivax_params(EIP_DAYS = -5L), "EIP_DAYS must be >= 0")
  expect_error(vivax_params(NU_HYPNOZOITE_BATCH = 1e9), "NU_HYPNOZOITE_BATCH")
  expect_error(vivax_params(NU_HYPNOZOITE_BATCH = 10, NU_HYPNOZOITE_BATCH_BOOST = 1e308,
                            RESURGENCE_TRIGGER_TICK = 5L), "NU_HYPNOZOITE_BATCH")
  expect_error(vivax_params(LAMBDA_SPOROZOITES = 1e19), "LAMBDA_SPOROZOITES must be <= 1000")
  expect_error(vivax_params(SEASONAL_AMP = -0.5), "SEASONAL_AMP must be >= 0")
  expect_error(vivax_params(RECRUD_DELAY_DIST = "fixed", MEAN_RECRUD_DELAY = 1e19),
               "MEAN_RECRUD_DELAY must be >= 1")
  expect_error(vivax_params(TREATMENT_DELAY = 1e10), "TREATMENT_DELAY must be >= 1")
  expect_error(vivax_params(MOSQ_LIFESPAN = 0.5), "MOSQ_LIFESPAN must be >= 1")
  expect_error(vivax_params(RECOMB_RATE = 0.6), "RECOMB_RATE must be at most 0.5")
  ## a boost counts only once its trigger falls inside the run
  expect_no_error(vivax_params(NU_HYPNOZOITE_BATCH = 10, NU_HYPNOZOITE_BATCH_BOOST = 1e9))
  ## large activation boosts stay legal: the model clamps them
  expect_no_error(vivax_params(ALPHA_ACTIVATE_BOOST = 1000, RESURGENCE_TRIGGER_TICK = 30L))
  ## a seasonal phase too large to work out the season with
  expect_error(vivax_params(SEASONAL_AMP = 0.5, SEASONAL_PHASE = 3e307), "cannot be worked out")
})

test_that("the oocyst count is capped by its mean, and K_OOCYSTS by what the sampler can draw", {
  expect_no_error(vivax_params(LAMBDA_OOCYSTS = 1000))
  expect_error(vivax_params(LAMBDA_OOCYSTS = 1000.5), "mean of 1000.5 oocysts per successful uptake")
  ## a small K gives a heavy tail: the mean is many times lambda
  expect_error(vivax_params(LAMBDA_OOCYSTS = 1000, K_OOCYSTS = 0.01), "oocysts per successful uptake")
  expect_error(vivax_params(LAMBDA_OOCYSTS = 1, K_OOCYSTS = 1e-4), "K_OOCYSTS must be 0 \\(Poisson\\) or at least 0.001")
  ## every value uptake_calibrate() returns is accepted
  for (target in list(c(3, 30), c(5, 100), c(10, 1000), c(2.5, 6))) {
    u <- suppressWarnings(uptake_calibrate(target[1], var = target[2]))  # some need a smaller k
    expect_no_error(vivax_params(LAMBDA_OOCYSTS = u$lambda, K_OOCYSTS = u$k))
  }
  expect_error(uptake_calibrate(1500), "at most 1000")
})

test_that("the import rate is capped where its draw would fail, not at the number of hosts", {
  ## one host, a window that closes before the boost starts: the model runs it
  expect_no_error(vivax_params(N_HOSTS = 1L, N_FOUNDERS = 1L, IMPORT_START = 1L, IMPORT_END = 50L,
                               IMPORT_RATE_PER_TICK = 0.6, RESURGENCE_TRIGGER_TICK = 100L,
                               IMPORT_RATE_BOOST = 2, SIM_DAYS = 200L))
  expect_error(vivax_params(IMPORT_START = 5L, IMPORT_END = 300L, IMPORT_RATE_PER_TICK = 1e19),
               "at most 1e6 imports a day")
})

test_that("a run of 99999 days or more must move the days that default to 'never'", {
  expect_error(vivax_params(SIM_DAYS = 99999L), "SIM_DAYS = 99999 reaches day 99999.*MDA_TICK")
  never <- malariafwd:::.never_params
  expect_length(never, 8L)
  later <- stats::setNames(as.list(rep(1e6, length(never))), never)
  expect_no_error(do.call(vivax_params, c(list(SIM_DAYS = 146000L), later)))
  expect_error(do.call(vivax_params, c(list(SIM_DAYS = 146000L), later[-1])), never[1])
  expect_error(do.call(vivax_params, c(list(SIM_DAYS = 1e9 + 1), later)), "at most 1e9")
})

test_that("strings that SLiM or params.tsv cannot carry are refused", {
  expect_error(vivax_params(SIM_ID = "a\nb"), "SIM_ID must not contain control characters")
  expect_error(vivax_params(SIM_ID = "a\tb"), "SIM_ID must not contain control characters")
  ## non-ASCII labels are not control characters, in any locale
  expect_no_error(vivax_params(SIM_ID = "Côte d’Ivoire 日本"))
  for (bad in c("", "1,", "a", "15", "0", "1:2", "1,1", "3-1", "1 4", "1-3 5", "1\t2",
                "99999999999", "1-99999999999"))
    expect_error(vivax_params(CHROMOSOMES = bad), "CHROMOSOMES must", label = shQuote(bad))
  expect_no_error(vivax_params(CHROMOSOMES = "1-3,5"))
  expect_no_error(vivax_params(CHROMOSOMES = "14"))
})

test_that("the biting warning follows the model's own seasonal, intervention and boost timeline", {
  expect_no_warning(vivax_params(BITING_RATE = 0.6, BITING_RATE_BOOST = 2))   # trigger after the run
  expect_warning(vivax_params(BITING_RATE = 0.6, BITING_RATE_BOOST = 2, RESURGENCE_TRIGGER_TICK = 100L),
                 "bites at most once a day")
  expect_no_warning(vivax_params(BITING_RATE = 0.07, BITING_RATE_BOOST = 16, RESURGENCE_TRIGGER_TICK = 100L,
                                 INTERVENTION_START = 50L, INTERVENTION_END = 400L,
                                 INTERVENTION_EFFICACY = 0.92))                # peak 0.07 x 0.08 x 16 = 0.09
  ## the boosted season peaks on day 100 + 8 x 365, not in the run's first year
  expect_warning(vivax_params(SIM_DAYS = 3650L, BITING_RATE = 0.6, SEASONAL_AMP = 0.5,
                              SEASONAL_PHASE = 100, BITING_RATE_BOOST = 1.2,
                              RESURGENCE_TRIGGER_TICK = 3000L), "reaches 1.08 on day 3020")
  ## a run far longer than any window is checked without building its whole timeline
  later <- stats::setNames(as.list(rep(1e9, 8)), malariafwd:::.never_params)
  expect_warning(do.call(vivax_params, c(list(SIM_DAYS = 1e8, BITING_RATE = 0.8, SEASONAL_AMP = 0.5),
                                         later)), "reaches 1.2 on day 365")
})

test_that("latent hypnozoites must be able to leave latency, and durations must be finite", {
  expect_error(vivax_params(K_LATENCY = 2L, DELTA_LATENCY = 0), "DELTA_LATENCY must be in \\(0, 1\\]")
  expect_no_error(vivax_params(K_LATENCY = 0L, DELTA_LATENCY = 0))
  expect_error(vivax_params(MEAN_BLOOD_STAGE_DURATION = Inf), "MEAN_BLOOD_STAGE_DURATION must be .*finite")
})
