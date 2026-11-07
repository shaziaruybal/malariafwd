## Writes the parameter tables of MODEL_DESCRIPTION.md section 10 from the parameter dictionary
## (R/param-dictionary.R). Run from the package root after editing the dictionary:
##   Rscript dev/update-param-tables.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
f <- "MODEL_DESCRIPTION.md"
txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
begin <- grep("^<!-- BEGIN parameter tables", txt)
end <- grep("^<!-- END parameter tables -->$", txt)
if (length(begin) != 1L || length(end) != 1L || end < begin)
  stop("MODEL_DESCRIPTION.md must hold one BEGIN and one END parameter-tables marker", call. = FALSE)
out <- c(txt[seq_len(begin)], "", strsplit(malariafwd:::.param_md(), "\n", fixed = TRUE)[[1]], "",
         txt[end:length(txt)])
writeLines(enc2utf8(out), f, useBytes = TRUE)
cat("wrote the parameter tables into", f, "\n")
