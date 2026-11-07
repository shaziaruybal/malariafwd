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
  expect_error(vivax_params(N_HOSTS = 5, N_FOUNDER_HOSTS = 10))
  expect_error(vivax_params(FOO = 1))
  expect_error(vivax_params(P_TRANS_VH = 2))
  expect_error(vivax_params(SIM_DAYS = -1))
})

test_that("N_FOUNDERS is refused, naming the founder parameters to set instead", {
  expect_error(vivax_params(N_FOUNDERS = 12), "^N_FOUNDERS is not a parameter: set N_FOUNDER_HOSTS")
  e <- tryCatch(vivax_params(N_FOUNDERS = 12, FOO = 1), error = conditionMessage)
  expect_match(e, "N_FOUNDERS is not a parameter")
  expect_match(e, "unknown parameter\\(s\\): FOO")
  expect_false("N_FOUNDERS" %in% names(vivax_params()))
  expect_identical(vivax_params(N_FOUNDER_HOSTS = 5)$N_FOUNDER_LINEAGES, 5L)
  ## a set edited to hold it is refused when run_sim() rebuilds it
  p <- vivax_params(); p$N_FOUNDERS <- 5
  expect_error(run_sim(p, outdir = tempfile("mfwd_nf_")), "^N_FOUNDERS is not a parameter: set N_FOUNDER_HOSTS")
})

test_that("N_OOCYSTS_PER_BITE is refused, naming LAMBDA_SPOROZOITES", {
  expect_error(vivax_params(N_OOCYSTS_PER_BITE = 1), "^N_OOCYSTS_PER_BITE is not a parameter.*LAMBDA_SPOROZOITES")
  expect_false("N_OOCYSTS_PER_BITE" %in% names(vivax_params()))
})

test_that("every argument is named, given once, and not NULL", {
  expect_error(vivax_params(300), "every argument must be named")
  expect_error(vivax_params(N_HOSTS = 300, 0.5), "every argument must be named")
  expect_error(vivax_params(N_HOSTS = 300, N_HOSTS = 200), "N_HOSTS given more than once")
  expect_error(vivax_params(DELTA_LATENCY = NULL), "DELTA_LATENCY cannot be NULL")
})

test_that("founder counts are checked against each other, with messages that name them", {
  expect_no_error(vivax_params(N_HOSTS = 5L, N_FOUNDER_HOSTS = 2L))
  expect_error(vivax_params(N_HOSTS = 5L), "N_FOUNDER_HOSTS \\(12\\) must be at most N_HOSTS \\(5\\)")
  expect_error(vivax_params(N_HOSTS = 0L), "N_HOSTS must be > 0")
  expect_error(vivax_params(N_FOUNDER_HOSTS = 5L, N_FOUNDER_LINEAGES = 6L),
               "N_FOUNDER_LINEAGES \\(6\\) must be at most N_FOUNDER_HOSTS \\(5\\)")
  expect_error(vivax_params(N_HOSTS = 1e5, N_FOUNDER_HOSTS = 2e5),
               "N_FOUNDER_HOSTS \\(200000\\) must be at most N_HOSTS \\(100000\\)")
})

test_that("an import path is open only on a day of the run", {
  ## imports arrive from day 1, so a pulse on day 0 or a window that ends on day 1 never comes
  expect_error(vivax_params(N_FOUNDER_HOSTS = 0, IMPORT_TICK = 0, N_IMPORTS = 5), "nothing can start this run")
  expect_error(vivax_params(N_FOUNDER_HOSTS = 0, IMPORT_START = 0, IMPORT_END = 1, N_IMPORTS_PER_TICK = 5),
               "nothing can start this run")
  expect_error(vivax_params(N_FOUNDER_HOSTS = 0, SIM_DAYS = 100, IMPORT_TICK = 101, N_IMPORTS = 5),
               "nothing can start this run")
  expect_no_error(vivax_params(N_FOUNDER_HOSTS = 0, IMPORT_TICK = 1, N_IMPORTS = 5))
  expect_no_error(vivax_params(N_FOUNDER_HOSTS = 0, IMPORT_START = 0, IMPORT_END = 2, N_IMPORTS_PER_TICK = 5))
})

test_that("run_sim() checks a parameter set edited after it was built", {
  p <- vivax_params(SIM_DAYS = 10L)
  p$N_HOST <- 7L                                            # a misspelt name
  expect_error(run_sim(p, outdir = tempfile()), "unknown parameter\\(s\\): N_HOST")
  p <- vivax_params(SIM_DAYS = 10L)
  p$N_IMPORTS <- 2.5
  expect_error(run_sim(p, outdir = tempfile()), "N_IMPORTS must be a whole number")
  q <- utils::modifyList(vivax_params(), list(ALPHA_ACTIVATE = FALSE))   # keeps the class
  expect_s3_class(q, "vivax_params")
  expect_error(run_sim(q, outdir = tempfile()), "ALPHA_ACTIVATE must be a number")
  ## SEED and OUTDIR belong to run_sim()
  p <- vivax_params(SIM_DAYS = 10L)
  p$SEED <- 5L
  expect_error(run_sim(p, outdir = tempfile()), "SEED and OUTDIR are set by run_sim\\(seed =, outdir =\\)")
  ## a stored call or symbol is refused as it stands, not evaluated
  p <- vivax_params(SIM_DAYS = 10L)
  p$SIM_ID <- quote(paste("run", 1))
  expect_error(run_sim(p, outdir = tempfile()), "SIM_ID must be a string, got paste")
  p <- vivax_params(SIM_DAYS = 10L)
  p$N_HOSTS <- quote(x)
  expect_no_warning(expect_error(run_sim(p, outdir = tempfile()), "N_HOSTS must be a number, got x"))
  ## an unedited set is rebuilt unchanged
  p <- vivax_params(N_HOSTS = 30L, SEASONAL_AMP = 0.4)
  expect_identical(do.call(vivax_params, as.list(unclass(p)), quote = TRUE), p)
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
    expect_no_error(suppressMessages(run_sim(do.call(vivax_params, cases[[nm]]), seed = 1, outdir = out)))
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

test_that("a day's fixed imports above N_HOSTS warn, since a host takes one import a day", {
  expect_warning(vivax_params(N_HOSTS = 30, IMPORT_TICK = 3, N_IMPORTS = 50),
                 "N_IMPORTS brings 50 imports on day 3, but each host takes at most one import a day, so at most N_HOSTS = 30 of them arrive.",
                 fixed = TRUE)
  expect_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 10, N_IMPORTS_PER_TICK = 31),
                 "N_IMPORTS_PER_TICK brings 31 imports a day", fixed = TRUE)
  expect_warning(vivax_params(N_HOSTS = 30, IMPORT_TICK = 5, N_IMPORTS = 20, IMPORT_START = 1,
                              IMPORT_END = 10, N_IMPORTS_PER_TICK = 20),
                 "N_IMPORTS and N_IMPORTS_PER_TICK bring 40 imports on day 5", fixed = TRUE)
  ## a pulse and the window's count above N_HOSTS on different days are both named
  expect_warning(vivax_params(N_HOSTS = 30, SIM_DAYS = 10, IMPORT_TICK = 1, N_IMPORTS = 40, IMPORT_START = 2,
                              IMPORT_END = 11, N_IMPORTS_PER_TICK = 40),
                 paste("N_IMPORTS brings 40 imports on day 1, and N_IMPORTS_PER_TICK brings 40 imports a day,",
                       "but each host takes at most one import a day"), fixed = TRUE)
  ## exactly N_HOSTS, a pulse or a window after the run, a window before day 1, or a rate well below
  ## N_HOSTS (P(Poisson(5) > 30) is about 1e-13)
  expect_no_warning(vivax_params(N_HOSTS = 30, IMPORT_TICK = 3, N_IMPORTS = 30))
  expect_no_warning(vivax_params(N_HOSTS = 30, SIM_DAYS = 100, IMPORT_TICK = 200, N_IMPORTS = 50))
  expect_no_warning(vivax_params(N_HOSTS = 30, SIM_DAYS = 100, IMPORT_START = 150, IMPORT_END = 300,
                                 N_IMPORTS_PER_TICK = 50))
  expect_no_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 0, IMPORT_END = 1, N_IMPORTS_PER_TICK = 50))
  expect_no_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 10, IMPORT_RATE_PER_TICK = 5))
})

test_that("an import rate that often brings more than N_HOSTS imports a day warns", {
  ## P(Poisson(50) > 30) = 0.998, so about 20 imports are lost almost every day
  expect_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 10, IMPORT_RATE_PER_TICK = 50),
                 paste("IMPORT_RATE_PER_TICK with its seasonal and resurgence multipliers reaches 50",
                       "imports a day on day 1: that day draws more than N_HOSTS = 30 imports",
                       "with probability 0.998, but each host takes at most one import a day"), fixed = TRUE)
  ## the peak follows the model's timeline: the boost from its trigger, and the season
  expect_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 200, IMPORT_RATE_PER_TICK = 15,
                              IMPORT_RATE_BOOST = 3, RESURGENCE_TRIGGER_TICK = 100L),
                 "reaches 45 imports a day on day 100:", fixed = TRUE)
  expect_no_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 50, IMPORT_RATE_PER_TICK = 15,
                                 IMPORT_RATE_BOOST = 3, RESURGENCE_TRIGGER_TICK = 100L))
  expect_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 300, IMPORT_RATE_PER_TICK = 15,
                              SEASONAL_AMP = 1, SEASONAL_PHASE = 100), "reaches 30 imports a day on day 100:",
                 fixed = TRUE)
  ## vector control scales biting, not imports
  expect_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 100, IMPORT_RATE_PER_TICK = 40,
                              INTERVENTION_START = 1L, INTERVENTION_END = 100L, INTERVENTION_EFFICACY = 0.9),
                 "with its seasonal and resurgence multipliers reaches 40 imports a day on day 1:", fixed = TRUE)
  ## the window's fixed count leaves less room for the rate
  expect_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 10, IMPORT_RATE_PER_TICK = 10,
                              N_IMPORTS_PER_TICK = 20), "imports together with N_IMPORTS_PER_TICK = 20 with probability",
                 fixed = TRUE)
  expect_no_warning(vivax_params(N_HOSTS = 30, IMPORT_START = 1, IMPORT_END = 10, IMPORT_RATE_PER_TICK = 10))
})

test_that("MAX_MOI is at least 1", {
  expect_error(vivax_params(MAX_MOI = 0L), "MAX_MOI must be at least 1, got 0: .* so 1 means no superinfection")
  expect_error(vivax_params(MAX_MOI = -2L), "MAX_MOI must be at least 1")
  expect_no_error(vivax_params(MAX_MOI = 1L))
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
  ## a founder copy is seeded on day 2, so with copies a round comes on day 2 or later
  expect_error(vivax_params(MDA_TICK = 1L, N_FOUNDER_LINEAGES = 4L),
               "MDA_TICK must not be 1 when a founder host carries a copy")
  expect_no_error(vivax_params(MDA_TICK = 2L, N_FOUNDER_LINEAGES = 4L))
  expect_no_error(vivax_params(MDA_TICK = 1L))
  expect_no_error(vivax_params(MDA_TICK = 1L, N_FOUNDER_HOSTS = 0L, IMPORT_TICK = 5L, N_IMPORTS = 3L))
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
  ## a value that is not a vector at all is named in the refusal, without R's own warnings
  expect_no_warning(expect_error(vivax_params(N_HOSTS = function() 1), "N_HOSTS must be a number"))
  ## a Date is not a day number, and a classed number whose text is not a number cannot reach SLiM
  expect_error(vivax_params(MDA_TICK = as.Date("1970-01-03")), "MDA_TICK must be a number")
  expect_no_error(vivax_params(N_HOSTS = I(200 + 1e-13)))
  registerS3method("format", "mfwd_duration", function(x, ...) paste0(unclass(x) * 86400, "s"))
  expect_error(vivax_params(TREATMENT_DELAY = structure(7, class = "mfwd_duration")),
               "TREATMENT_DELAY must be a plain number, got an object of class mfwd_duration; give the number itself",
               fixed = TRUE)
  ## I() is dropped, whatever the value's type; a number with another class is refused
  expect_identical(vivax_params(ALPHA_ACTIVATE = I(1/170))$ALPHA_ACTIVATE, as.numeric("0.00588235294117647"))
  expect_identical(vivax_params(SIM_DAYS = I(365L))$SIM_DAYS, 365L)
  expect_identical(vivax_params(SIM_ID = I("run_a"))$SIM_ID, "run_a")
  expect_identical(vivax_params(SIM_ID = I(factor("run_a")))$SIM_ID, "run_a")
  expect_identical(vivax_params(RETAIN_FULL_HISTORY = I(TRUE))$RETAIN_FULL_HISTORY, TRUE)
  expect_error(vivax_params(N_HOSTS = ts(50)), "N_HOSTS must be a plain number, got an object of class ts")
  expect_error(vivax_params(SIM_DAYS = as.difftime(3, units = "days")),
               "SIM_DAYS must be a number, got an object of class difftime", fixed = TRUE)
})

test_that("spellings of the same value are stored in the form SLiM is given", {
  ## a whole number with floating-point noise is whole as SLiM sees it
  p <- vivax_params(N_HOSTS = 100L, N_MOSQ = 1.1 * 100, N_FOUNDER_HOSTS = seq(0.1, 0.5, 0.1)[3] * 100,
                    SIM_DAYS = 365 * 2.2)
  expect_identical(c(p$N_MOSQ, p$N_FOUNDER_HOSTS, p$SIM_DAYS), c(110L, 30L, 803L))
  expect_true("N_MOSQ=110" %in% as_slim_args(p))
  ## ...and so is 0 with noise, while a small fraction is still refused
  p <- vivax_params(TREATMENT_START = 0.1 * 3 - 0.3, N_IMPORTS = 10 * (1 - 0.9) - 1)
  expect_true(p$TREATMENT_START == 0 && p$N_IMPORTS == 0)
  expect_true("TREATMENT_START=0" %in% as_slim_args(p))
  expect_error(vivax_params(N_IMPORTS = 1e-5), "N_IMPORTS must be a whole number")
  ## a named number, as coef() or MASS::fitdistr() give it, and a 1 x 1 matrix
  est <- c(size = 0.8, mu = 2.3)
  p <- vivax_params(LAMBDA_OOCYSTS = est["mu"], K_OOCYSTS = est["size"], N_HOSTS = c(n = 50L),
                    INTERVENTION_START = matrix(10L), INTERVENTION_END = 20L,
                    INTERVENTION_EFFICACY = crossprod(0.5, 1))
  expect_identical(p[c("LAMBDA_OOCYSTS", "K_OOCYSTS", "N_HOSTS", "INTERVENTION_START", "INTERVENTION_EFFICACY")],
                   list(LAMBDA_OOCYSTS = 2.3, K_OOCYSTS = 0.8, N_HOSTS = 50L, INTERVENTION_START = 10L,
                        INTERVENTION_EFFICACY = 0.5))
  ## latin1 text is stored as UTF-8, which SLiM reads
  lat <- "Par\xe1"; Encoding(lat) <- "latin1"
  expect_identical(charToRaw(vivax_params(SIM_ID = lat)$SIM_ID), charToRaw("Pará"))
  ## also as a factor (expand.grid() of latin1 text), and with R's CP1252 reading of latin1
  expect_identical(charToRaw(vivax_params(SIM_ID = factor(lat))$SIM_ID), charToRaw("Pará"))
  dash <- "a\x96b"; Encoding(dash) <- "latin1"
  expect_identical(charToRaw(vivax_params(SIM_ID = dash)$SIM_ID), as.raw(c(0x61, 0xe2, 0x80, 0x93, 0x62)))
  ## UTF-8 marked as bytes is stored as UTF-8, and refused by the parameter's own check when wrong
  byt <- "C\xc3\xb4te"; Encoding(byt) <- "bytes"
  expect_identical(Encoding(vivax_params(SIM_ID = byt)$SIM_ID), "UTF-8")
  expect_error(vivax_params(GENOME_MODE = byt), "GENOME_MODE must be 'full' or 'core'")
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

test_that("the oocyst count is capped by its mean, and K_OOCYSTS floored at 0.001", {
  expect_no_error(vivax_params(LAMBDA_OOCYSTS = 1000))
  expect_error(vivax_params(LAMBDA_OOCYSTS = 1000.5), "mean of 1000.5 oocysts per successful uptake")
  ## a small K gives a heavy tail: the mean is many times lambda
  expect_error(vivax_params(LAMBDA_OOCYSTS = 1000, K_OOCYSTS = 0.01), "oocysts per successful uptake")
  expect_error(vivax_params(LAMBDA_OOCYSTS = 1, K_OOCYSTS = 1e-4), "K_OOCYSTS must be 0 \\(Poisson\\) or at least 0.001")
  expect_no_error(vivax_params(LAMBDA_OOCYSTS = 1, K_OOCYSTS = 1e-3 * (1 - 1e-13)))
  ## a mean too large to work out with lambda / K is still refused by name
  expect_error(vivax_params(LAMBDA_OOCYSTS = 1e306, K_OOCYSTS = 0.001), "LAMBDA_OOCYSTS = 1e\\+306 .* gives a mean of")
  ## a value over the cap only by rounding is accepted, and a refused one prints above the cap
  expect_no_error(vivax_params(LAMBDA_OOCYSTS = 1000 * (1 + 1e-13)))
  expect_error(vivax_params(LAMBDA_OOCYSTS = 1000.0004), "mean of 1000.0004 oocysts")
  ## uptake_calibrate() holds the target mean, so every value it returns is accepted
  for (target in list(c(2.5, 6), c(1000, 1500), c(1000, 3e5), c(1000, 2e6), c(1000, 5e6),
                      c(999.999, 5e6))) {
    expect_no_warning(u <- uptake_calibrate(target[1], var = target[2]))
    expect_equal(u$achieved[["mean"]], target[1], tolerance = 1e-12)
    expect_equal(u$achieved[["var"]], target[2], tolerance = 1e-6)
    expect_no_error(vivax_params(LAMBDA_OOCYSTS = u$lambda, K_OOCYSTS = u$k))
  }
  ## a variance more than the model's zero-truncated negative binomial reaches at that mean:
  ## k = 0.001, the mean held
  for (target in list(c(3, 30), c(5, 100), c(10, 1000), c(1000, 1e7), c(300, 1.7e7))) {
    expect_warning(u <- uptake_calibrate(target[1], var = target[2]),
                   "more than the model's zero-truncated negative binomial reaches at mean")
    expect_identical(u$k, 1e-3)
    expect_equal(u$achieved[["mean"]], target[1], tolerance = 1e-12)
    expect_no_error(vivax_params(LAMBDA_OOCYSTS = u$lambda, K_OOCYSTS = u$k))
  }
  ## near a mean of 1 the Poisson variance is tiny; a target a little above it is still out of reach
  expect_warning(u <- uptake_calibrate(1 + 1e-9, var = 1.9e-9), "zero-truncated negative binomial")
  expect_identical(u$k, 1e-3)
  expect_error(uptake_calibrate(1500), "at most 1000")
})

test_that("uptake_moments() keeps its precision for a tiny lambda or a huge K", {
  expect_equal(uptake_moments(2, k = 1), c(mean = 3, var = 6))
  expect_equal(uptake_moments(2, k = 1e17), uptake_moments(2), tolerance = 1e-12)
  expect_equal(uptake_moments(1e-17), c(mean = 1, var = 0))
  expect_equal(uptake_moments(1e-300, k = 0.5), c(mean = 1, var = 0))
  ## a zero-truncated count has a mean of at least 1 and a variance of at least 0
  grid <- expand.grid(lambda = 10^seq(-300, 3, by = 3), k = c(0, 10^seq(-3, 306, by = 3)))
  m <- mapply(uptake_moments, grid$lambda, grid$k)
  expect_true(all(is.finite(m) & m["mean", ] >= 1 & m["var", ] >= 0))
  ## so these are accepted rather than refused as a mean of Inf
  expect_no_error(vivax_params(LAMBDA_OOCYSTS = 2, K_OOCYSTS = 1e17))
  expect_no_error(vivax_params(LAMBDA_OOCYSTS = 1e-17))
  expect_no_error(vivax_params(LAMBDA_OOCYSTS = 1e-300, K_OOCYSTS = 1e300))
  expect_named(uptake_moments(c(lambda = 2), k = c(k = 1)), c("mean", "var"))
  ## where lambda / K overflows, and where the variance is a rounding error from 0
  expect_equal(uptake_moments(1e306, 1e-3)[["mean"]], 1.9642810272581538e306, tolerance = 1e-12)
  expect_gte(uptake_moments(1.6e-16)[["var"]], 0)
  expect_gte(uptake_moments(1.14838e-17, 0.0605668)[["var"]], 0)
  expect_equal(uptake_moments(1e-300, 1e-300)[["var"]], uptake_moments(1e-100, 1e-100)[["var"]], tolerance = 1e-9)
})

test_that("the import rate is capped where its draw would fail, not at the number of hosts", {
  ## one host, a window that closes before the boost starts: the model runs it, and since the one
  ## host takes one import a day, the warning gives the days that draw more
  expect_warning(vivax_params(N_HOSTS = 1L, N_FOUNDER_HOSTS = 1L, IMPORT_START = 1L, IMPORT_END = 50L,
                              IMPORT_RATE_PER_TICK = 0.6, RESURGENCE_TRIGGER_TICK = 100L,
                              IMPORT_RATE_BOOST = 2, SIM_DAYS = 200L),
                 "reaches 0.6 imports a day on day 1: that day draws more than N_HOSTS = 1 imports with probability 0.122",
                 fixed = TRUE)
  expect_error(vivax_params(IMPORT_START = 5L, IMPORT_END = 300L, IMPORT_RATE_PER_TICK = 1e19),
               "at most 1e6 imports a day")
  expect_error(vivax_params(IMPORT_RATE_PER_TICK = 1e6 + 0.4), "reaches 1000000.4 imports a day")
  ## a rate worked out to sit at the cap is accepted, though it lands an ulp above it
  amp <- 0.997; b <- 9.722
  expect_no_error(vivax_params(IMPORT_RATE_PER_TICK = 1e6 / (b * (1 + amp)), IMPORT_RATE_BOOST = b,
                               SEASONAL_AMP = amp, RESURGENCE_TRIGGER_TICK = 100L))
})

test_that("a run of 99999 days or more must move the days that default to 'never'", {
  expect_error(vivax_params(SIM_DAYS = 99999L), "SIM_DAYS = 99999 reaches day 99999.*MDA_TICK")
  never <- malariafwd:::.never_params
  expect_length(never, 8L)
  later <- stats::setNames(as.list(rep(1e6, length(never))), never)
  expect_no_error(do.call(vivax_params, c(list(SIM_DAYS = 146000L), later)))
  expect_error(do.call(vivax_params, c(list(SIM_DAYS = 146000L), later[-1])), never[1])
  expect_error(do.call(vivax_params, c(list(SIM_DAYS = 1e9 + 1), later)), "at most 1e9")
  ## in a run that long, day 99999 itself reads as 'never'
  later$MDA_TICK <- 99999L
  expect_error(do.call(vivax_params, c(list(SIM_DAYS = 120000L), later)),
               "MDA_TICK is still 99999.*day 99999 itself cannot be used")
})

test_that("strings that SLiM or params.tsv cannot carry are refused", {
  expect_error(vivax_params(SIM_ID = "a\nb"), "SIM_ID must not contain control characters")
  expect_error(vivax_params(SIM_ID = "a\tb"), "SIM_ID must not contain control characters")
  ## C1 controls too, such as the U+0096 that iconv(from = "latin1") makes of a Windows dash
  expect_error(vivax_params(SIM_ID = "a\u0096b"), "SIM_ID must not contain control characters")
  expect_identical(vivax_params(SIM_ID = "a\u00a0b\u00a3")$SIM_ID, "a\u00a0b\u00a3")
  ## non-ASCII labels are not control characters, in any locale
  expect_no_error(vivax_params(SIM_ID = "Côte d’Ivoire 日本"))
  ## bytes that are not UTF-8 (a latin1 file read without its encoding) cannot run
  expect_error(vivax_params(SIM_ID = "Par\xe1"), "SIM_ID must be valid UTF-8.*CP1252")
  expect_no_warning(expect_error(vivax_params(CHROMOSOMES = "\xe9"), "CHROMOSOMES must be valid UTF-8"))
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

test_that("a double is stored as SLiM is given it", {
  p <- vivax_params(ALPHA_ACTIVATE = 1/170, P_TRANS_VH = 1 + 2e-16, SEASONAL_AMP = 0.1 + 0.2,
                    MU_HYPNOZOITE = 1e-300)
  expect_identical(p$ALPHA_ACTIVATE, as.numeric("0.00588235294117647"))
  expect_identical(p$P_TRANS_VH, 1)
  expect_identical(p$SEASONAL_AMP, 0.3)
  expect_identical(p$MU_HYPNOZOITE, 1e-300)
  a <- as_slim_args(p)
  expect_true(all(c("ALPHA_ACTIVATE=0.00588235294117647", "P_TRANS_VH=1.0") %in% a))
  expect_error(vivax_params(P_TRANS_VH = 1 + 1e-14), "P_TRANS_VH must be in \\[0, 1\\]")
})

test_that("each value has its default's type", {
  nm <- c("N_HOSTS", "SIM_DAYS", "P_RAD_SURVIVAL", "MDA_TICK")
  p <- vivax_params(N_HOSTS = 40, SIM_DAYS = 100, P_RAD_SURVIVAL = 0L, MDA_TICK = 99999)
  expect_identical(vapply(p[nm], typeof, ""), vapply(.malariafwd_defaults[nm], typeof, ""))
  ## a whole number beyond the integer range stays a double, and is sent whole
  p <- vivax_params(INTERVENTION_START = 3e9, INTERVENTION_END = 4e9)
  expect_identical(p$INTERVENTION_START, 3e9)
  expect_true("INTERVENTION_START=3000000000" %in% as_slim_args(p))
})

test_that("the defaults are a parameter set", {
  expect_identical(vivax_params(), structure(.malariafwd_defaults, class = "vivax_params"))
  expect_output(print(vivax_params(N_HOSTS = 40)), "all defaults")
})

test_that("validation works in doubles", {
  for (n in list(2000000000L, 2e9))
    expect_warning(vivax_params(N_HOSTS = 30, IMPORT_TICK = 5L, N_IMPORTS = n, IMPORT_START = 1L,
                                IMPORT_END = 10L, N_IMPORTS_PER_TICK = n),
                   "each host takes at most one import a day")
})

test_that("a rebuilt set is the set", {
  for (p in list(vivax_params(ALPHA_ACTIVATE = 1/170, SEASONAL_AMP = 0.1 + 0.2),
                 vivax_params(INTERVENTION_START = 3e9, INTERVENTION_END = 4e9),
                 vivax_params(N_HOSTS = 40, P_RAD_SURVIVAL = 0L)))
    expect_identical(do.call(vivax_params, as.list(unclass(p)), quote = TRUE), p)
})
