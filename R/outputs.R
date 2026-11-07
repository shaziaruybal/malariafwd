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

## older runs wrote lineage_ids under the name brood_ids
.OUTPUT_RENAMES <- list(events = character(), state = c(brood_ids = "lineage_ids"))

.na_of <- function(cls) switch(cls, integer = NA_integer_, double = NA_real_, character = NA_character_)

## One of a run's tables with the schema's columns, in its order and classes, whatever the run
## holds: fread() guesses a class from the values, so an empty file, or a column that happens to
## hold only NA, would otherwise come back logical. A column the file lacks is added as NA; one the
## schema lacks is kept after the schema's.
.read_output <- function(run, file, schema, renames = character()) {
  dir <- .run_dir(run)
  path <- file.path(dir, file)
  if (!file.exists(path))
    stop("no ", file, " in ", dir, ": it is not a malariafwd run directory, or the run did not finish.",
         call. = FALSE)
  if (file.size(path) == 0)
    stop(file, " in ", dir, " is empty: the run did not finish writing it.", call. = FALSE)
  hdr <- names(data.table::fread(path, sep = "\t", nrows = 0, showProgress = FALSE))
  now <- hdr
  for (old in names(renames))
    if (old %in% hdr && !renames[[old]] %in% hdr) now[hdr == old] <- renames[[old]]
  want <- schema[intersect(now, names(schema))]
  cls <- stats::setNames(c(integer = "integer", double = "numeric", character = "character")[want],
                         hdr[match(names(want), now)])
  w <- character()
  d <- withCallingHandlers(
    data.table::fread(path, sep = "\t", header = TRUE, quote = "", na.strings = "NA",
                      colClasses = split(names(cls), unname(cls)), showProgress = FALSE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  data.table::setnames(d, hdr, now)
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
  for (nm in setdiff(names(schema), now)) data.table::set(d, j = nm, value = .na_of(schema[[nm]]))
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
## split at its first tab, so a value keeps any tab or space it holds. The N_FOUNDERS row of a file
## made before it was retired gives its value to N_FOUNDER_HOSTS where the file has none.
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
  if ("N_FOUNDERS" %in% nm && !"N_FOUNDER_HOSTS" %in% nm) txt[["N_FOUNDER_HOSTS"]] <- txt[["N_FOUNDERS"]]
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
       missing = setdiff(names(d), c(names(txt), "N_FOUNDER_LINEAGES")),
       dropped = setdiff(nm, c(names(d), names(.RECORD_ROWS), "N_FOUNDERS")))
}
