## outputs.R -- the columns of the tables a run writes, and the reader that gives every run them.

## events.tsv and state.tsv, in the order the model writes them (logEvent() and
## writeStateSnapshot()), each column with its class.
.EVENTS_SCHEMA <- c(tick = "integer", event_id = "character", event_type = "character",
  inoculation_id = "integer", mosq_id = "integer", human_id = "integer", lineage_id = "integer",
  parent1_lineage_id = "integer", parent2_lineage_id = "integer", is_recombinant = "integer",
  p_outcross = "double", recipient = "character", is_dormant = "integer", is_imported = "integer",
  pedigree_id = "integer", parent1_pedigree_id = "integer", parent2_pedigree_id = "integer",
  n_hypnozoites = "integer", oocyst_id = "integer", detected = "integer",
  into_patent_lineage = "integer")

.STATE_SCHEMA <- c(tick = "integer", entity_type = "character", entity_id = "integer",
  tag = "integer", moi_active = "integer", moi_subpatent = "integer", moi_dormant = "integer",
  hyp_reservoir_size = "integer", lineage_ids = "character", n_transmissible = "integer",
  mosq_age_days = "integer", infectious_bites = "integer")

.na_of <- function(cls) switch(cls, integer = NA_integer_, double = NA_real_, character = NA_character_)

## One of a run's tables with the schema's columns, in its order and classes, whatever the run
## holds: fread() guesses a class from the values, so an empty file, or a column that happens to
## hold only NA, would otherwise come back logical. A column the file lacks is added as NA; one the
## schema lacks is kept after the schema's.
.read_output <- function(run, file, schema) {
  dir <- .run_dir(run)
  path <- file.path(dir, file)
  if (!file.exists(path))
    stop("no ", file, " in ", dir, ": it is not a malariafwd run directory, or the run did not finish.",
         call. = FALSE)
  if (file.size(path) == 0)
    stop(file, " in ", dir, " is empty: the run did not finish writing it.", call. = FALSE)
  hdr <- names(data.table::fread(path, sep = "\t", nrows = 0, showProgress = FALSE))
  want <- schema[intersect(hdr, names(schema))]
  cls <- stats::setNames(c(integer = "integer", double = "numeric", character = "character")[want],
                         names(want))
  w <- character()
  d <- withCallingHandlers(
    data.table::fread(path, sep = "\t", header = TRUE, quote = "", na.strings = "NA",
                      colClasses = split(names(cls), unname(cls)), showProgress = FALSE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  for (nm in names(want)) {
    got <- class(d[[nm]])[1]
    ok <- switch(want[[nm]], integer = got == "integer", double = got == "numeric",
                 character = got == "character")
    if (!ok) {
      bad <- if (is.character(d[[nm]]))
        which(is.na(suppressWarnings(as.numeric(d[[nm]]))) & d[[nm]] != "NA")[1] else NA
      stop(file, " in ", dir, ": column ", nm, " holds values that are not ",
           if (want[[nm]] == "double") "numbers" else paste0(want[[nm]], "s"),
           if (got %in% c("integer64", "numeric") && want[[nm]] == "integer")
             " (beyond the 32-bit integer range)"
           else if (!is.na(bad)) paste0(" (line ", bad + 1, ": \"", d[[nm]][bad], "\")"), call. = FALSE)
    }
  }
  for (nm in setdiff(names(schema), hdr)) data.table::set(d, j = nm, value = .na_of(schema[[nm]]))
  data.table::setcolorder(d, c(names(schema), setdiff(names(d), names(schema))))
  if (length(w)) warning(file, " in ", dir, ": ", paste(unique(w), collapse = "; "), call. = FALSE)
  d[]
}

## The rows of params.tsv that record the run rather than set it, each with its class: the seed SLiM
## ran with (beyond the 32-bit range it is kept as a double), the versions, the genome the model
## derived, and the outcome rows a run that went extinct or never started ends with.
.RECORD_ROWS <- c(ACTUAL_SEED = "integer", SEED = "integer", SLIM_VERSION = "character",
  MALARIAFWD_VERSION = "character", GENOME_LENGTH = "integer", N_CHROMOSOMES = "integer",
  RECOMB_WITHIN = "double", EXTINCT_TICK = "integer", NEVER_STARTED = "integer")

## A run's params.tsv, read as the parameter set it was made from: `params`, every parameter in
## the defaults' order and types (one the file lacks takes its default, and N_FOUNDER_LINEAGES
## follows N_FOUNDER_HOSTS as in vivax_params()); `record`, the record rows it has, typed; and
## `missing` and `dropped`, the parameters it lacks and the rows that are neither. Each line is
## split at its first tab, so a value keeps any tab or space it holds.
.read_params_tsv <- function(path) {
  dir <- dirname(path)
  where <- paste("params.tsv in", dir)
  if (!file.exists(path))
    stop("no params.tsv in ", dir, ": it is not a malariafwd run directory, or the run did not finish.",
         call. = FALSE)
  l <- sub("\r$", "", readLines(path, encoding = "UTF-8", warn = FALSE))
  if (!length(l) || l[1] != "param\tvalue")
    stop(where, " does not start with the line param<TAB>value.", call. = FALSE)
  l <- l[-1]
  tab <- regexpr("\t", l, fixed = TRUE)
  bad <- which(nzchar(l) & tab < 0)
  if (length(bad)) stop(where, ": line ", bad[1] + 1, " has no tab.", call. = FALSE)
  keep <- nzchar(l)
  l <- l[keep]; tab <- tab[keep]
  nm <- substr(l, 1, tab - 1)
  txt <- stats::setNames(substring(l, tab + 1), nm)
  if (anyDuplicated(nm)) stop(where, ": ", nm[duplicated(nm)][1], " is given twice.", call. = FALSE)
  d <- .malariafwd_defaults
  p <- utils::modifyList(d, .type_params(txt, where))
  if (!"N_FOUNDER_LINEAGES" %in% nm) p$N_FOUNDER_LINEAGES <- p$N_FOUNDER_HOSTS
  rec <- list()
  for (k in intersect(names(.RECORD_ROWS), nm)) {
    if (.RECORD_ROWS[[k]] == "character") { rec[[k]] <- txt[[k]]; next }
    v <- suppressWarnings(as.numeric(txt[[k]]))
    if (is.na(v)) stop(where, ": ", k, " is \"", txt[[k]], "\", not a number.", call. = FALSE)
    rec[[k]] <- if (.RECORD_ROWS[[k]] == "integer" && v == round(v) && abs(v) <= .Machine$integer.max)
      as.integer(v) else v
  }
  list(params = structure(p, class = "vivax_params"), record = rec,
       missing = setdiff(names(d), names(txt)),
       dropped = setdiff(nm, c(names(d), names(.RECORD_ROWS))))
}

## ---- the column dictionary (?malariafwd-outputs; MODEL_DESCRIPTION 16.1 and 16.2) ----

## What each column of events.tsv and state.tsv holds, in the schemas' order. One source for the
## help page and the manuscript tables, as the parameter dictionary is for the parameters.
.EVENTS_MEANING <- c(
  tick = "the simulation day",
  event_id = "`\"<tick>_<counter>\"`, unique",
  event_type = "the kind of event; `event_types()` lists them",
  inoculation_id = "the bite or import that brought the parasite, on bite rows (`vh_inoculation` and its no-take form) and import rows (an import's batch row, logged the next day, carries its primary's id); -1 elsewhere",
  mosq_id = "the mosquito, on uptake, bite and `mosq_death` rows; -1 elsewhere",
  human_id = "the host; -1 only on `mosq_death` rows",
  lineage_id = "the parasite's lineage, its clonal group; NA on recombinant uptakes, on course rows about the host (`bloodcure_given`, `radcure_given`, `radcure_failure` and their MDA forms) and on `vh_inoculation_import_no_take`",
  parent1_lineage_id = "on uptake rows, the lineage of the first gamete; NA elsewhere",
  parent2_lineage_id = "on recombinant uptake rows, the lineage of the second gamete; NA elsewhere",
  is_recombinant = "1 on an outcrossed uptake (`hv_uptake_recombinant`), 0 elsewhere",
  p_outcross = "on uptake rows, the Simpson diversity of the gametocyte pool the mosquito drew from, to 6 significant digits, a diagnostic; 0 elsewhere",
  recipient = "`\"mosquito\"` on uptake and `mosq_death` rows, `\"human\"` on the others",
  is_dormant = "1 on a hypnozoite batch's rows (its establishment, deaths and drug clearance), 0 on blood-stage and host rows, -1 on no-take rows",
  is_imported = "1 when the parasite was acquired by importation, carried within the host to its batches, relapses and recrudescences; a local infection by an import's offspring is 0",
  pedigree_id = "the genome, the key to the tree sequence (pedigree 0 is a real genome); -1 on uptake, course and no-take rows",
  parent1_pedigree_id = "the genome a row's genome came from: a gamete of an uptake, the mosquito product of a bite, the batch of a relapse, the suppressed clone of a recrudescence; -1 where none",
  parent2_pedigree_id = "the second gamete of a recombinant uptake; -1 elsewhere",
  n_hypnozoites = "the hypnozoites a row stands for: a batch's size when established, the number activating, dying or cleared, or the reservoir surviving a radical-cure failure; 1 on other rows",
  oocyst_id = "on an uptake row, the oocyst the mosquito made; on a local bite's blood-stage and batch rows, every relapse of that batch and the drug clearances of its clones, the oocyst of the infecting product; -1 for founders, imports and their batches and relapses, recrudescences, and rows that do not record it (`human_clearance`, `hypnozoite_clearance`, `mosq_death`, course rows)",
  detected = "on the blood-stage onset rows of bites, imports, relapses and recrudescences, 1 if the observation layer reported the case and 0 if not; -1 elsewhere, founder seeding included",
  into_patent_lineage = "on `hypnozoite_activation` rows, 1 if a patent clone of the same lineage, born on an earlier day, was still in the host (one that ends that day counts), 0 if not; -1 elsewhere")

.STATE_MEANING <- c(
  tick = "the simulation day",
  entity_type = "`\"human\"` or `\"mosquito\"`",
  entity_id = "the host's or the mosquito's id",
  tag = "hosts: 1 if the host carries a patent blood-stage clone, or a subpatent one when `SUBPATENT_INFECTIVITY > 0`, whether or not it carries gametocytes; 0 otherwise, hypnozoites alone included; mosquitoes: 1 if it carries a parasite product, inside or past the EIP",
  moi_active = "hosts: the distinct patent blood-stage lineages (the COI), suppressed clones excluded; mosquitoes: the distinct lineages it carries, inside the EIP included",
  moi_subpatent = "hosts: the distinct lineages suppressed below patency, present but undetectable; 0 unless `P_RECRUD_GIVEN_FAIL > 0`",
  moi_dormant = "hosts: the distinct genotypes in the hypnozoite reservoir",
  hyp_reservoir_size = "hosts: the hypnozoites, the sum of the batches' counts",
  lineage_ids = "the lineages present, `|`-separated and de-duplicated, even for one; NA when none",
  n_transmissible = "mosquitoes: the distinct lineages past the EIP; NA for hosts",
  mosq_age_days = "mosquitoes: the age in days; NA for hosts",
  infectious_bites = "hosts: the bites received that day from a mosquito carrying a product past the EIP, whatever the bite then did, the EIR's numerator; NA for mosquitoes")

## The dictionary of `which` ("events" or "state") as a markdown table (for MODEL_DESCRIPTION) or
## an Rd list (for the help, where a table does not wrap in the console).
.dictionary_md <- function(which) {
  sch <- if (which == "events") .EVENTS_SCHEMA else .STATE_SCHEMA
  mean <- if (which == "events") .EVENTS_MEANING else .STATE_MEANING
  stopifnot(identical(names(sch), names(mean)))
  ## a `|` in a meaning would end the table cell
  paste(c("| column | class | meaning |", "|---|---|---|",
          sprintf("| `%s` | %s | %s |", names(sch), unname(sch), gsub("|", "\\|", unname(mean), fixed = TRUE))),
        collapse = "\n")
}
.dictionary_rd <- function(which) {
  sch <- if (which == "events") .EVENTS_SCHEMA else .STATE_SCHEMA
  mean <- if (which == "events") .EVENTS_MEANING else .STATE_MEANING
  stopifnot(identical(names(sch), names(mean)))
  rd <- function(x) {
    x <- gsub("\\\\", "\\\\\\\\", x); x <- gsub("%", "\\\\%", x); x <- gsub("([{}])", "\\\\\\1", x)
    gsub("`([^`]*)`", "\\\\code{\\1}", x)
  }
  paste0("\\describe{\n", paste(sprintf("\\item{\\code{%s} (%s)}{%s}", names(sch), unname(sch), rd(unname(mean))),
                                collapse = "\n"), "\n}")
}
.outputs_rd <- function() {
  sec <- function(title, body) paste0("\\section{", title, "}{\n", body, "\n}")
  c(sec("events.tsv", paste0("One row per event, in the order the events were logged. -1 means not applicable and ",
                             "NA no lineage. A run with no event gives no row, and a column a file lacks reads as ",
                             "NA.\n\n", .dictionary_rd("events"))),
    sec("state.tsv", paste0("One row per host and per mosquito per day. ",
                            "\\code{load_state()} gives the host rows, and \\code{entity = \"mosquito\"} or ",
                            "\\code{\"all\"} the mosquitoes' or both.\n\n", .dictionary_rd("state"))))
}

#' The files a run writes, and their columns
#'
#' A run's directory holds:
#' * `events.tsv`, one row per event, read by [load_events()];
#' * `state.tsv`, one row per host and per mosquito per day, read by [load_state()];
#' * `params.tsv`, the parameter set and the run's record, read by [load_params()] and [load_run()];
#' * `plasmo.trees`, the tree sequence, found by [trees_path()] and read by [ibd_pairs()],
#'   [ibd_segments()], [arg_summary()] and [simulate_genotypes()];
#' * `model.slim`, `cmdline.txt`, `run.log` and `DONE`: the model run, the command, its log, and the
#'   record that the run finished, with the size of every file.
#'
#' `params.tsv` has a `param` and a `value` column: the record rows `ACTUAL_SEED`, `SLIM_VERSION` and
#' `MALARIAFWD_VERSION`, then every parameter with `SEED`, with the derived `GENOME_LENGTH`,
#' `N_CHROMOSOMES` and `RECOMB_WITHIN`, and last `EXTINCT_TICK` or `NEVER_STARTED` when the run went
#' extinct or never started. Whole numbers are written as such, other numbers to 15 significant
#' digits with `.0` on a whole one, flags as `T` or `F`, and text as given.
#'
#' `plasmo.trees` holds one node per genome, keyed by its pedigree id, with node times in ticks before
#' the end of the run; SLiM's second, vacant node of each genome is not a genome. A run with
#' `RETAIN_FULL_HISTORY = FALSE` keeps the genomes alive at its end, and most earlier ones leave it.
#' @evalRd .outputs_rd()
#' @name malariafwd-outputs
#' @aliases events.tsv state.tsv params.tsv plasmo.trees
#' @seealso [load_run()], [load_events()], [load_state()], [load_params()], [trees_path()],
#'   [event_types()]
NULL
