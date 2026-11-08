## Writes the generated tables of MODEL_DESCRIPTION.md: the parameter tables of section 10, from
## the parameter dictionary (R/param-dictionary.R), and the events.tsv and state.tsv tables of
## sections 16.1 and 16.2, from the column dictionary (R/outputs.R). Run from the package root after
## editing either:
##   Rscript dev/update-doc-tables.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
f <- "MODEL_DESCRIPTION.md"
txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
fill <- function(txt, key, body) {
  begin <- grep(paste0("^<!-- BEGIN ", key), txt, fixed = FALSE)
  end <- grep(paste0("^<!-- END ", key, " -->$"), txt)
  if (length(begin) != 1L || length(end) != 1L || end < begin)
    stop("MODEL_DESCRIPTION.md must hold one BEGIN and one END ", key, " marker", call. = FALSE)
  c(txt[seq_len(begin)], "", strsplit(body, "\n", fixed = TRUE)[[1]], "", txt[end:length(txt)])
}
txt <- fill(txt, "parameter tables", malariafwd:::.param_md())
txt <- fill(txt, "events.tsv table", malariafwd:::.dictionary_md("events"))
txt <- fill(txt, "state.tsv table", malariafwd:::.dictionary_md("state"))
writeLines(enc2utf8(txt), f, useBytes = TRUE)
cat("wrote the generated tables into", f, "\n")
