## The parameter dictionary (R/param-dictionary.R): every parameter documented once, its documented
## bounds the ones vivax_params() enforces, and the generated documentation current.

dict <- malariafwd:::.param_dict

test_that("every parameter has a record and is in exactly one group", {
  expect_identical(names(dict), names(vivax_params()))
  gp <- unlist(lapply(malariafwd:::.param_groups, `[[`, "params"), use.names = FALSE)
  expect_setequal(gp, names(dict))
  expect_identical(anyDuplicated(gp), 0L)
  pi <- param_info()
  expect_s3_class(pi, "data.table")
  expect_identical(names(pi), c("name", "group", "default", "type", "unit", "allowed", "meaning",
                                "source", "setting"))
  expect_true(all(vapply(pi, is.character, logical(1))))
  expect_identical(pi$name, gp)
  expect_true(all(nzchar(pi$unit) & nzchar(pi$meaning) & nzchar(pi$allowed)))
  expect_identical(param_info("N_HOSTS")$unit, "hosts")
  expect_identical(param_info(c("BITING_RATE", "N_HOSTS"))$name, c("BITING_RATE", "N_HOSTS"))
  expect_error(param_info("FOO"), "unknown parameter\\(s\\): FOO")
})

test_that("the records' types are the model's", {
  src <- readLines(slim_model(), warn = FALSE)
  block <- function(name) {
    i <- grep(paste0("^\\s*", name, " = c\\("), src)[1]
    j <- i; while (!grepl("\\);\\s*$", src[j])) j <- j + 1L
    txt <- paste(src[i:j], collapse = " ")
    gsub('"', "", regmatches(txt, gregexpr('"[A-Z_]+"', txt))[[1]])
  }
  d <- vivax_params()
  expect_setequal(setdiff(block("whole_params"), "SEED"), names(Filter(is.integer, d)))
  expect_setequal(block("real_params"), names(Filter(is.double, d)))
  expect_setequal(block("flag_params"), names(Filter(is.logical, d)))
  expect_setequal(setdiff(block("text_params"), "OUTDIR"), names(Filter(is.character, d)))
  expect_setequal(block("_pn"), c(names(d), "SEED", "GENOME_LENGTH", "N_CHROMOSOMES", "RECOMB_WITHIN"))
})

test_that("the documented bounds are the ones vivax_params() enforces, and refusals name the parameter", {
  never <- stats::setNames(as.list(rep(2e9, length(malariafwd:::.never_params))), malariafwd:::.never_params)
  ## the other values a parameter's own bound needs to be reachable
  ctx <- list(N_HOSTS = list(N_FOUNDER_HOSTS = 1L, N_FOUNDER_LINEAGES = 1L),
              N_FOUNDER_HOSTS = list(IMPORT_TICK = 1L, N_IMPORTS = 1L, N_FOUNDER_LINEAGES = 0L),
              N_FOUNDER_LINEAGES = list(N_FOUNDER_HOSTS = 0L, IMPORT_TICK = 1L, N_IMPORTS = 1L),
              DELTA_LATENCY = list(K_LATENCY = 1L), SIM_DAYS = never)
  try1 <- function(nm, v) {
    r <- tryCatch(suppressWarnings(do.call(vivax_params, c(ctx[[nm]], stats::setNames(list(v), nm)))),
                  error = function(e) e)
    if (inherits(r, "error")) conditionMessage(r) else "ok"
  }
  for (nm in names(dict)) {
    r <- dict[[nm]]
    step <- function(x) if (is.integer(r$default)) 1 else max(1e-6, abs(x) * 1e-6)
    check <- function(v, ok) {
      got <- try1(nm, v)
      expect_identical(identical(got, "ok"), ok, label = paste(nm, "=", format(v), if (!ok) "refused" else "accepted"))
      if (!ok) expect_true(grepl(nm, got, fixed = TRUE), label = paste("the refusal of", nm, "=", format(v), "names it"))
    }
    if (length(r$values)) {
      for (v in r$values) check(v, TRUE)
      check("not-a-value", FALSE)
      next
    }
    if (!is.na(r$min)) {
      check(r$min, !r$min_open)
      if (r$min_open) check(r$min + step(r$min), TRUE)
      check(r$min - step(r$min), FALSE)
    }
    if (!is.na(r$max)) {
      check(r$max, TRUE)
      check(r$max + step(r$max), FALSE)
    }
  }
  ## K_OOCYSTS: 0, or at least 0.001
  expect_identical(try1("K_OOCYSTS", 0.0005) == "ok", FALSE)
  expect_identical(try1("K_OOCYSTS", 0.001), "ok")
})

test_that("the generated documentation is ASCII and lists each parameter once, by group", {
  for (txt in list(malariafwd:::.param_rd(), malariafwd:::.param_md())) {
    expect_false(grepl("[^\x01-\x7e]", txt, useBytes = TRUE))
    for (nm in names(dict))
      expect_identical(lengths(regmatches(txt, gregexpr(paste0("(`|\\\\code\\{)", nm, "(`|\\})"), txt))), 1L, label = nm)
  }
  md <- malariafwd:::.param_md()
  expect_identical(regmatches(md, gregexpr("(?m)^### .*$", md, perl = TRUE))[[1]],
                   paste("###", names(malariafwd:::.param_groups)))
})

test_that("?vivax_params and MODEL_DESCRIPTION section 10 hold the current tables", {
  rd <- test_path("..", "..", "man", "vivax_params.Rd")
  skip_if_not(file.exists(rd), "not in the source tree")
  expect_true(grepl(malariafwd:::.param_rd(), paste(readLines(rd, warn = FALSE), collapse = "\n"), fixed = TRUE),
              label = "man/vivax_params.Rd holds .param_rd() (run devtools::document())")
  md <- test_path("..", "..", "MODEL_DESCRIPTION.md")
  skip_if_not(file.exists(md), "not in the source tree")
  txt <- paste(readLines(md, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  block <- regmatches(txt, regexpr("(?s)<!-- BEGIN parameter tables.*?-->\n\n(.*?)\n<!-- END parameter tables -->", txt, perl = TRUE))
  expect_length(block, 1L)
  expect_identical(sub("(?s)^<!-- BEGIN parameter tables.*?-->\n\n(.*)\n<!-- END parameter tables -->$", "\\1", block, perl = TRUE),
                   malariafwd:::.param_md(), label = "MODEL_DESCRIPTION.md's tables (run dev/update-doc-tables.R)")
})

test_that("a parameter set prints what differs by value, and gives its table", {
  expect_output(print(vivax_params(N_HOSTS = 40)), "all defaults")
  p <- vivax_params(); p$N_HOSTS <- 40
  expect_output(print(p), "all defaults")
  d <- as.data.frame(vivax_params(N_HOSTS = 300))
  expect_identical(names(d), c("name", "value", "default", "changed", "group", "unit", "meaning"))
  expect_identical(d$name[d$changed], "N_HOSTS")
  expect_identical(d$value[d$name == "N_HOSTS"], "300")
  expect_identical(d$value[d$name == "ALPHA_ACTIVATE"], "0.00588235294117647")
})
