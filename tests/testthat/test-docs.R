## The documents: the help pages have examples and give reasons, not decision codes; the vignettes
## run their analysis on the bundled run; and the metadata files, the design decisions and the
## validation results agree with one another and with the repository.

## the source's man/ under load_all(), the installed help otherwise (R CMD check)
rd_db <- function() {
  man <- system.file("man", package = "malariafwd")
  if (nzchar(man) && dir.exists(man)) tools::Rd_db(dir = dirname(man)) else tools::Rd_db("malariafwd")
}
rd_flat <- function(x) if (is.list(x)) paste(vapply(x, rd_flat, ""), collapse = "") else paste(as.character(x), collapse = "")
rd_tagged <- function(rd, tag) Filter(function(e) identical(attr(e, "Rd_tag"), tag), rd)

test_that("the author's name, ORCID and roles agree across the metadata files", {
  root <- test_path("..", "..")
  skip_if_not(file.exists(file.path(root, "CITATION.cff")), "the metadata files are in the source tree only")
  read <- function(f) paste(readLines(file.path(root, f), encoding = "UTF-8", warn = FALSE), collapse = "\n")
  name <- "Ruybal-Pesántez"
  orcid <- "0000-0002-0495-179X"
  for (f in c("DESCRIPTION", "inst/CITATION", "LICENSE", "LICENSE.md", "CITATION.cff", ".zenodo.json", "README.md")) {
    txt <- read(f)
    expect_match(txt, name, fixed = TRUE, label = f)
    expect_no_match(txt, "Ruybal[- ]Pesantez|Ruybal Pes", label = f)
  }
  for (f in c("DESCRIPTION", "inst/CITATION", "CITATION.cff", ".zenodo.json"))
    expect_match(read(f), orcid, fixed = TRUE, label = f)
  au <- eval(parse(text = read.dcf(file.path(root, "DESCRIPTION"), fields = "Authors@R")[[1]]))
  expect_setequal(unlist(au$role), c("aut", "cre", "cph"))
  ## the funder is acknowledged where the release is described and archived
  for (f in c("README.md", ".zenodo.json"))
    expect_match(read(f), "LCF/BQ/PI25/12100035", fixed = TRUE, label = f)
})

## the numbers of the design decisions' entries, in file order
decision_numbers <- function() {
  f <- test_path("..", "..", "DECISIONS.md")
  skip_if_not(file.exists(f), "DECISIONS.md is in the source tree only")
  h <- grep("^### D[0-9]+ ", readLines(f, encoding = "UTF-8", warn = FALSE), value = TRUE)
  as.integer(sub("^### D([0-9]+) .*", "\\1", h))
}

test_that("every decision from D1 to the newest has one entry", {
  d <- decision_numbers()
  expect_setequal(d, seq_len(max(d)))
  expect_false(anyDuplicated(d) > 0)
})

test_that("every decision the code and the docs cite is an entry", {
  d <- decision_numbers()
  root <- test_path("..", "..")
  src <- c(file.path(root, c("MODEL_DESCRIPTION.md", "README.md", "NEWS.md")),
           list.files(file.path(root, c("R", "inst/slim", "inst/validation", "tests/testthat", "vignettes")),
                      "\\.(R|slim|Rmd)$", full.names = TRUE))
  stray <- character()
  for (f in src[file.exists(src)]) {
    txt <- readLines(f, encoding = "UTF-8", warn = FALSE)
    cited <- unlist(regmatches(txt, gregexpr("(?<![A-Za-z0-9_])D[0-9]{1,2}(?![0-9A-Za-z_])", txt, perl = TRUE)))
    n <- setdiff(as.integer(substring(cited, 2)), d)
    if (length(n)) stray <- c(stray, paste0(basename(f), ": D", n))
  }
  expect_identical(stray, character())
})

test_that("the help pages give reasons, not decision codes or private files", {
  db <- rd_db()
  txt <- vapply(db, function(rd) gsub("\\s+", " ", rd_flat(rd)), "")
  expect_identical(names(txt)[grepl("\\bD[0-9]{1,2}\\b", txt, perl = TRUE)], character())
  expect_identical(names(txt)[grepl("analys[ie]s repo|malariafwd-analyses|prototype run|audit run", txt)], character())
})

test_that("every export's help page has examples, and none is \\dontrun", {
  db <- rd_db()
  ex <- getNamespaceExports("malariafwd")
  tag <- function(e) { t <- attr(e, "Rd_tag"); if (is.null(t)) "" else t }
  has <- vapply(db, function(rd) {
    if (!any(vapply(rd_tagged(rd, "\\alias"), rd_flat, "") %in% ex)) return(NA_character_)
    exs <- rd_tagged(rd, "\\examples")
    if (!length(exs)) return("no examples")
    if ("\\dontrun" %in% vapply(exs[[1]], tag, "")) return("\\dontrun")
    "ok"
  }, "")
  has <- has[!is.na(has)]
  expect_gte(length(has), length(ex) - 1L)   # every export has a page (`malariafwd` is a package alias)
  expect_identical(names(has)[has != "ok"], character())
})

test_that("the vignettes run their analysis, on the bundled run", {
  dir <- test_path("..", "..", "vignettes")
  skip_if_not(dir.exists(dir), "the vignettes' sources are in the source tree only")
  for (f in list.files(dir, "\\.Rmd$", full.names = TRUE)) {
    txt <- readLines(f, encoding = "UTF-8", warn = FALSE)
    setup <- grep("opts_chunk\\$set", txt, value = TRUE)
    expect_false(any(grepl("eval\\s*=\\s*FALSE", setup)), label = basename(f))
    expect_true(any(grepl("example_run()", txt, fixed = TRUE)), label = basename(f))
  }
})

test_that("the validation results record a run of every checking script", {
  f <- system.file("validation", "RESULTS.md", package = "malariafwd")
  expect_true(nzchar(f))
  txt <- paste(readLines(f, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  expect_match(txt, "at commit [0-9a-f]{7}")
  scripts <- dir(dirname(f), "\\.R$")
  for (s in scripts) expect_match(txt, paste0("`", s, "`"), fixed = TRUE, label = s)
})

test_that("the design decisions name no file that is not in the repository", {
  root <- test_path("..", "..")
  f <- file.path(root, "DECISIONS.md")
  skip_if_not(file.exists(f), "DECISIONS.md is in the source tree only")
  ## the repository's files: from git in a clone (which leaves out the untracked .claude/), or the
  ## files of an exported copy, which holds only those
  ours <- if (dir.exists(file.path(root, ".git")))
    suppressWarnings(system2("git", c("-C", shQuote(root), "ls-files"), stdout = TRUE, stderr = FALSE))
  else list.files(root, recursive = TRUE, all.files = TRUE)
  skip_if(length(ours) == 0L, "the repository's files cannot be listed")
  txt <- paste(readLines(f, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  ticks <- gsub("^`|`$", "", regmatches(txt, gregexpr("`[^`\n]+`", txt))[[1]])
  files <- unique(sub(":[0-9, :]+$", "", grep("^[A-Za-z0-9_./-]+\\.(R|Rmd|md|slim|py|yaml|yml|json|cff)(:[0-9, :]+)?$",
                                               ticks, value = TRUE)))
  run_files <- "model.slim"   # the copy of the model in every run directory
  ok <- files %in% ours | basename(files) %in% basename(ours) | files %in% run_files
  expect_identical(files[!ok], character())
})
