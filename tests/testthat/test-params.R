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

test_that("the biting warning follows the model's own seasonal, intervention and boost timeline", {
  expect_no_warning(vivax_params(BITING_RATE = 0.8, SEASONAL_AMP = -0.5))     # the model ignores AMP < 0
  expect_no_warning(vivax_params(BITING_RATE = 0.6, BITING_RATE_BOOST = 2))   # trigger after the run
  expect_warning(vivax_params(BITING_RATE = 0.6, BITING_RATE_BOOST = 2, RESURGENCE_TRIGGER_TICK = 100L),
                 "bites at most once a day")
  expect_no_warning(vivax_params(BITING_RATE = 0.07, BITING_RATE_BOOST = 16, RESURGENCE_TRIGGER_TICK = 100L,
                                 INTERVENTION_START = 50L, INTERVENTION_END = 400L,
                                 INTERVENTION_EFFICACY = 0.92))                # peak 0.07 x 0.08 x 16 = 0.09
})

test_that("latent hypnozoites must be able to leave latency, and durations must be finite", {
  expect_error(vivax_params(K_LATENCY = 2L, DELTA_LATENCY = 0), "DELTA_LATENCY must be in \\(0, 1\\]")
  expect_no_error(vivax_params(K_LATENCY = 0L, DELTA_LATENCY = 0))
  expect_error(vivax_params(MEAN_BLOOD_STAGE_DURATION = Inf), "MEAN_BLOOD_STAGE_DURATION must be >= 1")
})
