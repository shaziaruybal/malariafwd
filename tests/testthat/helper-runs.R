## Run directories written by hand, as the model writes them, for tests that need no SLiM.

ev_cols <- c("tick", "event_id", "event_type", "inoculation_id", "mosq_id", "human_id", "lineage_id",
             "parent1_lineage_id", "parent2_lineage_id", "is_recombinant", "p_outcross", "recipient",
             "is_dormant", "is_imported", "pedigree_id", "parent1_pedigree_id", "parent2_pedigree_id",
             "n_hypnozoites", "oocyst_id", "detected", "into_patent_lineage")
st_cols <- c("tick", "entity_type", "entity_id", "tag", "moi_active", "moi_subpatent", "moi_dormant",
             "hyp_reservoir_size", "lineage_ids", "n_transmissible", "mosq_age_days", "infectious_bites")

## A value as the model's paramText() writes it: a float to 15 significant digits, with a decimal
## point when it is whole, and anything else as Eidos pastes it.
param_text <- function(v, float) {
  if (is.logical(v)) return(if (v) "T" else "F")
  if (!is.numeric(v)) return(v)
  if (!float) return(malariafwd:::.slim_num(v))
  s <- sprintf("%.15g", v)
  if (grepl("[.e]", s)) s else paste0(s, ".0")
}

## The lines of params.tsv for parameter set `p` and record rows `record`, without the rows named
## in `drop` and with the lines `extra` added.
params_tsv_lines <- function(p = vivax_params(), record = list(), drop = character(), extra = character()) {
  d <- malariafwd:::.malariafwd_defaults
  rec <- malariafwd:::.RECORD_ROWS
  rows <- c(vapply(names(p), function(k) paste0(k, "\t", param_text(p[[k]], is.double(d[[k]]))), ""),
            vapply(names(record), function(k) paste0(k, "\t", param_text(record[[k]], identical(rec[[k]], "double"))), ""))
  c("param\tvalue", rows[!sub("\t.*", "", rows) %in% drop], extra)
}

## A run directory with the given event and state rows (each a tab-separated line) under the given
## headers, and a params.tsv for `params` (overrides to vivax_params()) and `record`; with `done`,
## a DONE file that marks it a whole run.
fake_run <- function(ev_rows = character(), st_rows = character(), ev_head = ev_cols, st_head = st_cols,
                     params = list(), record = list(ACTUAL_SEED = 1L, SLIM_VERSION = "5.2",
                                                    MALARIAFWD_VERSION = "0.1.0"),
                     drop = character(), extra = character(), done = TRUE) {
  d <- tempfile("mfwd_fake_"); dir.create(d)
  writeLines(c(paste(ev_head, collapse = "\t"), ev_rows), file.path(d, "events.tsv"))
  writeLines(c(paste(st_head, collapse = "\t"), st_rows), file.path(d, "state.tsv"))
  writeLines(params_tsv_lines(do.call(vivax_params, params), record, drop, extra),
             file.path(d, "params.tsv"), useBytes = TRUE)
  if (done) refresh_done(d)
  d
}

## Writes a run directory's DONE again, for the files it now holds.
refresh_done <- function(d)
  malariafwd:::.write_done(d, seed = 1L, run_key = "fake", model_md5 = "fake", elapsed_s = 0,
                           slim_version = "5.2")
