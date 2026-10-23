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
                                 ALPHA_ACTIVATE = 0, MU_HYPNOZOITE = 0, N_HOSTS = 100))
  expect_true("P_RAD_SURVIVAL=0.0" %in% a)
  expect_true("BLOODCURE_COVERAGE=1.0" %in% a)
  expect_true("ALPHA_ACTIVATE=0.0" %in% a)
  expect_true("MU_HYPNOZOITE=0.0" %in% a)
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
    min_alpha   = c(base, list(ALPHA_ACTIVATE = 0)),
    min_mu      = c(base, list(MU_HYPNOZOITE = 0))
  )
  for (nm in names(cases)) {
    out <- file.path(tempdir(), paste0("mfw-zero-one-", nm))
    unlink(out, recursive = TRUE)
    expect_no_error(run_sim(do.call(vivax_params, cases[[nm]]), seed = 1, outdir = out))
  }
})
