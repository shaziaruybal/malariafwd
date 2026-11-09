## The help pages describe a run in one way, and the package page defines the terms they share.

## the source's man/ under load_all(), the installed help otherwise (R CMD check)
rd_db <- function() {
  man <- system.file("man", package = "malariafwd")
  if (nzchar(man) && dir.exists(man)) tools::Rd_db(dir = dirname(man)) else tools::Rd_db("malariafwd")
}
rd_flat <- function(x) if (is.list(x)) paste(vapply(x, rd_flat, ""), collapse = "") else paste(as.character(x), collapse = "")
rd_tagged <- function(rd, tag) Filter(function(e) identical(attr(e, "Rd_tag"), tag), rd)
## the text of a topic's \arguments item `run`, its white space squished; NULL without one
run_arg <- function(rd) {
  for (args in rd_tagged(rd, "\\arguments"))
    for (it in rd_tagged(args, "\\item"))
      if (identical(trimws(rd_flat(it[[1]])), "run")) return(gsub("\\s+", " ", trimws(rd_flat(it[[2]]))))
  NULL
}

test_that("every help page describes a run in the same words", {
  db <- rd_db()
  names(db) <- sub("\\.Rd$", "", basename(names(db)))
  txt <- unlist(lapply(db, run_arg))
  canon <- "a run: a vivax_run (from run_sim() or load_run()) or the path of a run directory."
  trees <- c("ibd_pairs", "arg_summary", "arg_tmrca")   # these also take a .trees file
  expect_identical(unique(unname(txt[!names(txt) %in% trees])), canon)
  expect_gte(sum(!names(txt) %in% trees), 15L)
  expect_true(all(startsWith(txt[trees], sub("\\.$", ",", canon))))
})

test_that("the package page is a page of concepts that names the shared terms", {
  db <- rd_db()
  pkg <- db[[grep("malariafwd-package(\\.Rd)?$", names(db))]]
  expect_false("internal" %in% vapply(rd_tagged(pkg, "\\keyword"), rd_flat, ""))
  text <- rd_flat(pkg)
  for (term in c("tick", "pedigree_id", "lineage_id", "COI", "MOI", "Rc", "vivax_run"))
    expect_match(text, term, fixed = TRUE, label = term)
})

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
