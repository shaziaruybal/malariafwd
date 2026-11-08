## genotypes.R -- marker genotypes painted on the root genomes of a run's tree sequence, and their
## VCF (D75).

## Every root genome of the tree sequence (a node no edge has as its child: a day-1 founder template
## or an import), and the root genome each of the genomes `ids` descends from at each position.
## The tree is simplified to `ids` with keep_input_roots, so its roots are the root genomes, whatever
## genomes are kept. No random numbers here.
.genotypes_py <- '
def gt_roots(path, ids, positions):
    ts = tskit.load(path)
    ids = [int(p) for p in np.atleast_1d(ids)] if ids is not None else []
    pos = np.asarray([float(x) for x in np.atleast_1d(positions)] if positions is not None else [], float)
    out = dict(root_id=np.zeros(0, np.int64), root_tick=np.zeros(0, np.int64),
               roots=np.zeros((0, len(pos)), np.int32), missing=np.zeros(0, np.int64))
    if ts.num_nodes == 0:
        out["missing"] = np.array(ids, np.int64)
        return out
    md = ts.metadata
    end_tick = int(md["SLiM"]["tick"]) if isinstance(md, dict) and "SLiM" in md else 0
    p2n, ped = _ped_nodes(ts)
    nodes = ts.tables.nodes
    raw = tskit.unpack_bytes(nodes.metadata, nodes.metadata_offset)
    vac = np.array([m[8] & 1 for m in raw], dtype=bool)
    child = np.zeros(ts.num_nodes, bool); child[ts.edges_child] = True
    parent = np.zeros(ts.num_nodes, bool); parent[ts.edges_parent] = True
    sample = (nodes.flags & tskit.NODE_IS_SAMPLE) > 0
    rn = np.where(~child & ~vac & (parent | sample))[0]
    out["root_id"] = ped[rn].astype(np.int64)
    out["root_tick"] = np.rint(end_tick - ts.nodes_time[rn]).astype(np.int64)
    missing = [p for p in ids if p not in p2n]
    out["missing"] = np.array(missing, np.int64)
    if missing or not ids:
        out["roots"] = np.zeros((len(ids), len(pos)), np.int32)
        return out
    n = len(ids)
    sts = ts.simplify(samples=[p2n[p] for p in ids], keep_input_roots=True)
    del ts
    _, sped = _ped_nodes(sts)
    R = np.empty((n, len(pos)), np.int32)
    tree = sts.first(); last = -1; col = None
    for j in np.argsort(pos, kind="stable"):
        tree.seek(pos[j])
        if tree.index != last:
            col = np.full(n, -1, np.int32)
            for r in tree.roots:
                sub = tree.preorder(r)
                col[sub[sub < n]] = sped[r]
            if (col < 0).any():
                raise ValueError("a genome has no root genome at position %g" % pos[j])
            last = tree.index
        R[:, j] = col
    out["roots"] = R
    return out
'

## The run a genotype call reads: a vivax_run or its directory, whose params.tsv gives the
## chromosome layout and whose events.tsv the root genomes' origin.
.genotype_run <- function(run) {
  if (is.character(run) && length(run) == 1L && !dir.exists(run) && grepl("\\.trees$", run))
    stop("`run` must be a vivax_run or a run directory, not a .trees file: the genotypes need the run's ",
         "chromosome layout (params.tsv) and its events (events.tsv).", call. = FALSE)
  load_run(run)
}

## The loci table, checked against the run's chromosomes `cm`, with each locus's 0-based position on
## the concatenated genome added.
.genotype_loci <- function(cm, loci) {
  if (!is.data.frame(loci))
    stop("`loci` must be a data frame of chrom, pos and plaf, one row per biallelic locus.", call. = FALSE)
  lo <- data.table::copy(data.table::as.data.table(loci))
  if (!"plaf" %in% names(lo) && "PLAF" %in% names(lo))
    stop("`loci` needs `plaf`, the frequency of the ALT allele; for a table of the REF frequency (recombuddy's ",
         "PLAF), use plaf = 1 - PLAF.", call. = FALSE)
  miss <- setdiff(c("chrom", "pos", "plaf"), names(lo))
  if (length(miss)) stop("`loci` lacks the column(s) ", paste(miss, collapse = ", "), ".", call. = FALSE)
  row <- function(i, what) stop("loci row ", i, ": ", what, call. = FALSE)
  ch <- as.character(lo$chrom)
  pv <- regmatches(ch, regexec("^PvP01_([0-9]{2})_v[0-9]+$", ch))
  ch <- ifelse(lengths(pv) == 2L, vapply(pv, function(m) if (length(m) == 2L) m[2] else NA_character_, ""), ch)
  num <- suppressWarnings(as.numeric(ch))
  bad <- which(is.na(num) | num != round(num) | !num %in% 1:14)
  if (length(bad)) row(bad[1], paste0("chrom ", lo$chrom[bad[1]], " is not a PvP01 chromosome (1-14, or a name such as PvP01_03_v2)."))
  k <- match(num, cm$chrom)
  bad <- which(is.na(k))
  if (length(bad)) row(bad[1], paste0("chromosome ", num[bad[1]], " is not simulated in this run (CHROMOSOMES gives ",
                                      paste(cm$chrom, collapse = ", "), ")."))
  num_col <- function(x, name) {
    if (is.numeric(x) || all(is.na(x))) return(as.numeric(x))
    stop("`loci$", name, "` must be numbers.", call. = FALSE)
  }
  pos <- num_col(lo$pos, "pos")
  bad <- which(is.na(pos) | pos != round(pos) | pos < 1 | pos > cm$length[k])
  if (length(bad)) {
    i <- bad[1]
    row(i, if (is.na(pos[i]) || pos[i] != round(pos[i]) || pos[i] < 1)
             paste0("pos must be a whole number from 1 (1-based within the chromosome), not ", pos[i], ".")
           else paste0("pos ", .slim_num(pos[i]), " is beyond chromosome ", num[i], ", which is ",
                       .slim_num(cm$length[k[i]]), " bp in this run."))
  }
  pl <- num_col(lo$plaf, "plaf")
  bad <- which(is.na(pl) | !is.finite(pl) | pl < 0 | pl > 1)
  if (length(bad)) row(bad[1], paste0("plaf must be a frequency from 0 to 1, not ", pl[bad[1]], "."))
  dup <- which(duplicated(data.table::data.table(num, pos)))
  if (length(dup)) {
    j <- dup[1]; i <- which(num == num[j] & pos == pos[j])[1]
    stop("loci rows ", i, " and ", j, " are both chromosome ", num[j], ", pos ", .slim_num(pos[j]), ".", call. = FALSE)
  }
  ## set() rather than `[`: a column of `loci` named num, pos, pl, k or cm would take the place of
  ## these values inside `[`, and other columns are the user's to keep
  data.table::set(lo, j = c("chrom", "pos", "plaf"), value = list(as.integer(num), as.integer(pos), pl))
  data.table::set(lo, j = "locus", value = if ("locus" %in% names(lo)) as.character(lo$locus) else paste0(lo$chrom, ":", lo$pos))
  dup <- which(duplicated(lo$locus))
  if (length(dup)) stop("loci rows ", match(lo$locus[dup[1]], lo$locus), " and ", dup[1], " have the same locus id ",
                        lo$locus[dup[1]], ".", call. = FALSE)
  data.table::set(lo, j = "position", value = cm$start[k] + pos - 1)
  data.table::setcolorder(lo, c("locus", "chrom", "pos", "position", "plaf"))
  lo[]
}

## Host-sample calls from the clones' alleles: 0 when every clone of a sample carries REF, 1 when
## every one carries ALT, 0.5 when both are present. Samples in order of first appearance.
.pool_calls <- function(alleles, sample_id) {
  sids <- unique(sample_id); si <- match(sample_id, sids)
  if (!nrow(alleles)) return(matrix(numeric(), 0L, ncol(alleles), dimnames = list(character(), colnames(alleles))))
  s <- rowsum(alleles, si, reorder = FALSE); k <- tabulate(si, length(sids))
  out <- ifelse(s == 0, 0, ifelse(s == k, 1, 0.5))
  matrix(as.numeric(out), length(sids), ncol(alleles), dimnames = list(sids, colnames(alleles)))
}

## Observed calls: each non-missing call is wrong with probability `error` (a clone call flips; a host
## call becomes one of the other two), then missing with probability `missing`. A rate of 0 draws
## nothing.
.observe_calls <- function(true_calls, type, error, missing) {
  calls <- true_calls
  if (error > 0 && length(calls)) {
    e <- which(stats::runif(length(calls)) < error)
    if (type == "clone") calls[e] <- 1 - calls[e]
    else {
      ## one of the other two calls, at random
      t0 <- true_calls[e]; pick <- stats::runif(length(e)) < 0.5
      calls[e] <- ifelse(t0 == 0, ifelse(pick, 1, 0.5), ifelse(t0 == 1, ifelse(pick, 0.5, 0), ifelse(pick, 1, 0)))
    }
  }
  if (missing > 0 && length(calls)) calls[stats::runif(length(calls)) < missing] <- NA
  calls
}

## A vivax_genotypes object, its parts in order.
.new_vivax_genotypes <- function(calls, true_calls, alleles, roots, samples, clones, loci, root_genomes,
                                 root_alleles, chromosomes, type, error, missing, seed, run)
  structure(list(calls = calls, true_calls = true_calls, alleles = alleles, roots = roots, samples = samples,
                 clones = clones, loci = loci, root_genomes = root_genomes, root_alleles = root_alleles,
                 chromosomes = chromosomes, type = type, error = error, missing = missing, seed = seed, run = run),
            class = "vivax_genotypes")

## A single rate from 0 to 1, checked.
.check_rate <- function(x, name)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || x > 1)
    stop("`", name, "` must be a single number from 0 to 1.", call. = FALSE)

#' Simulate marker genotypes from a run's genealogy
#'
#' Biallelic genotypes at chosen loci for chosen genomes, painted on the genealogy the run recorded.
#' The model carries no mutation, so a run's genomes differ only by descent: at every position each
#' genome descends from exactly one **root genome**, a day-1 founder template or an imported genome,
#' and root genomes are unrelated. Each root genome draws its allele at each locus independently, ALT
#' with probability `plaf`, and every genome carries the allele of the root genome it descends from
#' there. `roots` records that root genome for every genome and locus, the truth behind the calls:
#' two genomes share a root genome at a locus exactly where [ibd_segments()] finds them sharing a
#' common ancestor.
#'
#' Each genome can be its own sample (a clone), or the genomes can be grouped into host samples by
#' `sample_id`, as a blood sample holds every clone in the host: a host call is 0 when every clone
#' carries REF, 1 when every one carries ALT, and 0.5 when both are present.
#'
#' @section Reproducibility:
#' The root genomes' alleles are drawn first, for every root genome of the run and every locus, so with
#' one seed, run and loci a genome's alleles do not depend on which other genomes are asked for, how
#' they are grouped, or the error and missingness rates; then the errors, then the missing calls.
#' Large sets of genomes can therefore be genotyped in batches with one seed, and several error rates
#' laid over one truth. A number in `seed` seeds this call alone and leaves the session's random
#' stream as it was; `NULL` draws from the session's stream.
#' @section What painting leaves out:
#' Root genomes are unrelated, so there is no background relatedness or linkage disequilibrium among
#' them; imports draw from the same frequencies as the founders; an allele never changes along the
#' genealogy; a host sample records which alleles are present, not their proportions; loci are
#' biallelic; and `plaf` is the frequency among the root genomes, which drifts among the genomes
#' sampled later.
#' @section Which genomes can be genotyped:
#' Any genome in the tree sequence: case clones, hypnozoite batches and mosquito products alike. A
#' run with `RETAIN_FULL_HISTORY = FALSE` records only the genomes alive at its end; a run that never
#' started holds none.
#'
#' @param run a `vivax_run` or its directory.
#' @param ids the pedigree ids of the genomes to genotype.
#' @param loci a data frame, one row per biallelic locus: `chrom` (the PvP01 chromosome, as a number or
#'   a name such as `PvP01_03_v2`), `pos` (1-based, within the chromosome as this run simulated it) and
#'   `plaf` (the frequency of the ALT allele among the root genomes); optionally `locus` (an id, by
#'   default `"<chrom>:<pos>"`), and `ref` and `alt` for [write_vcf()]. Other columns are kept.
#' @param sample_id `NULL` for clone samples, one per genome, or a vector as long as `ids` grouping the
#'   genomes into host samples; a genome can be in several samples.
#' @param error the probability that a call is wrong: a clone call flips, and a host call becomes one
#'   of the other two at random.
#' @param missing the probability that a call is missing, after the errors.
#' @param seed `NULL`, or a whole number that seeds this call alone.
#' @return a `vivax_genotypes` object, a list of:
#'   * `calls`: a sample x locus matrix of 0, 1, 0.5 (host samples only) and `NA`;
#'   * `true_calls`: the calls before error and missingness;
#'   * `alleles`: an integer genome x locus matrix of 0 (REF) and 1 (ALT), by pedigree id;
#'   * `roots`: an integer genome x locus matrix, the pedigree id of the root genome behind each allele;
#'   * `samples` (`sample_id`, `n_clones`) and `clones` (`sample_id`, `pedigree_id`, in input order);
#'   * `loci`: `locus`, `chrom`, `pos`, `position` (0-based bp along the concatenated genome, as in the
#'     IBD segments), `plaf`, then the other columns given;
#'   * `root_genomes` (`root_id`, `tick`: day 1 for a founder template, the arrival day for an import,
#'     and `imported`) and `root_alleles`, their alleles;
#'   * `chromosomes`, the run's [chrom_map()]; `type` (`"clone"` or `"host"`), `error`, `missing`,
#'     `seed` (`NA` when not set) and `run`.
#' @examples
#' \dontrun{
#' cm <- chrom_map(run)
#' loci <- data.frame(chrom = rep(cm$chrom, each = 5), pos = 1000 * rep(1:5, nrow(cm)), plaf = 0.3)
#' cs <- classify_cases(load_events(run))
#' g <- simulate_genotypes(run, cs$pedigree_id, loci, error = 0.01, missing = 0.05, seed = 1)
#' write_vcf(g, file.path(tempdir(), "cases.vcf"))
#' }
#' @seealso [write_vcf()], [ibd_segments()], [ibd_pairs()], [chrom_map()], [episodes()]
#' @inheritSection tskit_available Python
#' @export
simulate_genotypes <- function(run, ids, loci, sample_id = NULL, error = 0, missing = 0, seed = NULL) {
  r <- .genotype_run(run)
  ids <- .ibd_ids(ids, "ids")
  type <- if (is.null(sample_id)) "clone" else "host"
  if (type == "host") {
    if (!is.atomic(sample_id) || length(sample_id) != length(ids))
      stop("`sample_id` must be a vector as long as `ids`, one sample id per genome.", call. = FALSE)
    sample_id <- as.character(sample_id)
    if (anyNA(sample_id) || any(!nzchar(sample_id))) stop("`sample_id` must not be NA or empty.", call. = FALSE)
    if (anyDuplicated(paste(sample_id, ids, sep = "\r")))
      stop("a genome is in one host sample only once; each `ids` and `sample_id` pair must differ.", call. = FALSE)
  } else {
    if (anyDuplicated(ids))
      stop("each genome is genotyped once; give `sample_id` to put one genome in several host samples.", call. = FALSE)
    sample_id <- as.character(ids)
  }
  .check_rate(error, "error"); .check_rate(missing, "missing")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
                         seed != round(seed) || abs(seed) > .Machine$integer.max))
    stop("`seed` must be NULL or a single whole number.", call. = FALSE)
  cm <- chrom_map(r)
  lo <- .genotype_loci(cm, loci)
  if (isTRUE(r$never_started) && length(ids))
    stop("this run never started, so its tree sequence holds no genome to genotype (see ?load_run).", call. = FALSE)
  u <- unique(ids)
  if (isTRUE(r$never_started)) {
    res <- list(root_id = integer(), root_tick = integer(), roots = matrix(integer(), 0L, nrow(lo)), missing = integer())
  } else {
    .require_python("simulate_genotypes()", c("tskit", "numpy"))
    tf <- trees_path(r)
    if (!file.exists(tf)) stop("tree sequence not found: ", tf, call. = FALSE)
    m <- .py_module("_malariafwd_genotypes", .genotypes_py)
    res <- reticulate::py_to_r(m$gt_roots(tf, as.list(u), as.list(lo$position)))
  }
  absent <- as.integer(res$missing)
  if (length(absent)) {
    hint <- c(if (isFALSE(r$params$RETAIN_FULL_HISTORY))
                "a run with RETAIN_FULL_HISTORY = FALSE records only the genomes alive at its end",
              if (!is.na(r$extinct_tick)) paste0("this run went extinct on day ", r$extinct_tick, ", so none"))
    stop("pedigree id(s) ", paste(utils::head(absent, 5), collapse = ", "),
         if (length(absent) > 5) paste0(" and ", length(absent) - 5, " more"),
         " are not genomes in this run's tree sequence",
         if (length(hint)) paste0(" (", paste(hint, collapse = "; "), ")"), ".", call. = FALSE)
  }
  rg <- data.table::data.table(root_id = as.integer(res$root_id), tick = as.integer(res$root_tick))
  data.table::setorder(rg, root_id)
  ev <- load_events(r)
  rg[, imported := data.table::fifelse(root_id %in% ev[event_type == "vh_inoculation_import", pedigree_id], TRUE,
                       data.table::fifelse(root_id %in% ev[event_type == "vh_inoculation_seed", pedigree_id], FALSE, NA))]
  roots <- matrix(as.integer(res$roots), length(u), nrow(lo))
  if (length(roots) && !all(roots %in% rg$root_id))
    stop("internal: a genome's root is not a root genome of the run.", call. = FALSE)
  L <- nrow(lo); F <- nrow(rg)
  draw <- function() {
    ra <- matrix(stats::rbinom(F * L, 1L, rep(lo$plaf, each = F)), F, L)
    al <- matrix(ra[cbind(match(roots, rg$root_id), as.vector(col(roots)))], nrow(roots), L)
    storage.mode(al) <- "integer"
    tc <- if (type == "clone") { x <- al * 1; x } else .pool_calls(al[match(ids, u), , drop = FALSE], sample_id)
    list(ra = ra, al = al, tc = tc, calls = .observe_calls(tc, type, error, missing))
  }
  d <- if (is.null(seed)) draw() else .with_seed(seed, draw())
  sids <- unique(sample_id)
  dimnames(d$ra) <- list(as.character(rg$root_id), lo$locus)
  dimnames(d$al) <- dimnames(roots) <- list(as.character(u), lo$locus)
  if (type == "clone") d$tc <- d$tc[match(ids, u), , drop = FALSE]
  if (type == "clone") d$calls <- d$calls[match(ids, u), , drop = FALSE]
  dimnames(d$tc) <- dimnames(d$calls) <- list(sids, lo$locus)
  .new_vivax_genotypes(
    calls = d$calls, true_calls = d$tc, alleles = d$al, roots = roots,
    samples = data.table::data.table(sample_id = sids, n_clones = tabulate(match(sample_id, sids), length(sids))),
    clones = data.table::data.table(sample_id = sample_id, pedigree_id = ids), loci = lo, root_genomes = rg,
    root_alleles = d$ra, chromosomes = cm, type = type, error = error, missing = missing,
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed), run = r$outdir)
}

#' @rdname simulate_genotypes
#' @param x a `vivax_genotypes` object.
#' @param ... ignored.
#' @export
print.vivax_genotypes <- function(x, ...) {
  n <- nrow(x$calls)
  cat("<vivax_genotypes> ", n, if (x$type == "clone") " clone samples" else
        paste0(" host samples (", data.table::uniqueN(x$clones$pedigree_id), " genomes)"),
      " x ", ncol(x$calls), " loci on ", data.table::uniqueN(x$loci$chrom), " chromosomes\n", sep = "")
  cat("  painted on ", nrow(x$root_genomes), " root genomes (", sum(x$root_genomes$imported, na.rm = TRUE),
      " imported); error ", format(x$error), ", missing ", format(x$missing), ", ",
      if (is.na(x$seed)) "seed not set" else paste("seed", x$seed), "\n", sep = "")
  cat("  run: ", x$run, "\n", sep = "")
  invisible(x)
}

#' Write simulated genotypes as a VCF file
#'
#' Writes a [simulate_genotypes()] object as VCF 4.2: one contig per simulated PvP01 chromosome, with
#' the run's lengths; each locus at its 1-based position within its chromosome (so in `"core"` mode,
#' a position on the shortened chromosome), its id in `ID`, and the ALT frequency it was painted with
#' as `PLAF` in `INFO`, rather than `AF`, which tools read as the frequency among the samples. Records
#' are sorted by chromosome and position. Clone samples are written haploid (`0`, `1`, `.`) and host
#' samples unphased diploid, as VCFs called from mixed *Plasmodium* infections are: `0/0`, `1/1`, `0/1`
#' for both alleles present, `./.`. The model has no sequence, so without `ref` and `alt` columns in
#' the loci REF and ALT are the placeholders A and C, and the header says so.
#'
#' @param x a `vivax_genotypes` object, from [simulate_genotypes()].
#' @param file the path to write; a name ending in `.gz` is gzip-compressed (for tabix, write a plain
#'   file and compress it with bgzip).
#' @param calls `"observed"` (the calls with error and missingness) or `"true"`.
#' @param chrom_names `NULL`, or the 14 contig names by PvP01 number; by default `PvP01_01_v2` to
#'   `PvP01_14_v2`, and `as.character(1:14)` gives numbered contigs.
#' @return `file`, invisibly.
#' @examples
#' \dontrun{write_vcf(g, file.path(tempdir(), "genotypes.vcf.gz"))}
#' @seealso [simulate_genotypes()]
#' @export
write_vcf <- function(x, file, calls = c("observed", "true"), chrom_names = NULL) {
  if (!inherits(x, "vivax_genotypes"))
    stop("`x` must be a vivax_genotypes object, from simulate_genotypes(), not a ", class(x)[1], ".", call. = FALSE)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file))
    stop("`file` must be a single path.", call. = FALSE)
  calls <- match.arg(calls)
  if (is.null(chrom_names)) chrom_names <- sprintf("PvP01_%02d_v2", 1:14)
  if (!is.character(chrom_names) || length(chrom_names) != 14L || anyNA(chrom_names) || anyDuplicated(chrom_names) ||
      any(!nzchar(chrom_names)) || any(grepl("[[:space:]]", chrom_names)))
    stop("`chrom_names` must be 14 distinct contig names without spaces, one per PvP01 chromosome.", call. = FALSE)
  M <- if (calls == "observed") x$calls else x$true_calls
  sids <- rownames(M); if (is.null(sids)) sids <- character()
  if (any(grepl("[[:space:]]", sids))) stop("a sample id contains whitespace, which VCF does not allow: ",
                                          sids[grepl("[[:space:]]", sids)][1], ".", call. = FALSE)
  lo <- data.table::copy(x$loci)
  if (any(grepl("[[:space:];]", lo$locus)))
    stop("a locus id contains whitespace or ';', which VCF does not allow: ", lo$locus[grepl("[[:space:];]", lo$locus)][1],
         ".", call. = FALSE)
  given <- all(c("ref", "alt") %in% names(lo))
  if (given) {
    ok <- function(v) is.character(v) & grepl("^[ACGTN]$", v)
    if (!all(ok(lo$ref) & ok(lo$alt)) || any(lo$ref == lo$alt))
      stop("`ref` and `alt` must each be one of A, C, G, T or N, and differ.", call. = FALSE)
  }
  lo[, .row := .I]
  data.table::setorder(lo, chrom, pos)
  code <- M[, lo$.row, drop = FALSE]
  lev <- if (identical(x$type, "clone")) c("0", "1", "0/1", ".") else c("0/0", "1/1", "0/1", "./.")
  k <- ifelse(is.na(code), 4L, ifelse(code == 0, 1L, ifelse(code == 1, 2L, 3L)))
  gt <- matrix(lev[k], nrow(code))
  cm <- x$chromosomes
  hdr <- c("##fileformat=VCFv4.2",
           sprintf("##source=malariafwd %s simulate_genotypes()", utils::packageVersion("malariafwd")),
           paste0("##malariafwd_run=", x$run),
           paste0("##malariafwd_samples=", x$type),
           paste0("##malariafwd_calls=", if (calls == "true") "true" else
             sprintf("observed (error=%s, missing=%s)", format(x$error), format(x$missing))),
           paste0("##malariafwd_seed=", if (is.na(x$seed)) "NA" else x$seed),
           if (!given) "##malariafwd_alleles=REF and ALT are placeholders: the model has no sequence",
           sprintf("##contig=<ID=%s,length=%s>", chrom_names[cm$chrom], .slim_num(cm$length)),
           '##INFO=<ID=PLAF,Number=A,Type=Float,Description="Frequency of the ALT allele among the root genomes it was painted on">',
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
  body <- data.table::data.table(CHROM = chrom_names[lo$chrom], POS = lo$pos, ID = lo$locus,
                                 REF = if (given) lo$ref else rep("A", nrow(lo)),
                                 ALT = if (given) lo$alt else rep("C", nrow(lo)), QUAL = ".", FILTER = "PASS",
                                 INFO = paste0("PLAF=", formatC(lo$plaf, digits = 15, format = "g")))
  body[, INFO := gsub(" ", "", INFO, fixed = TRUE)]
  if (nrow(code)) {
    body[, FORMAT := "GT"]
    gd <- data.table::as.data.table(t(gt))
    data.table::setnames(gd, sids)
    body <- cbind(body, gd)
  }
  data.table::setnames(body, 1L, "#CHROM")
  gz <- grepl("\\.gz$", file)
  con <- if (gz) gzfile(file, "w") else file(file, "w")
  writeLines(hdr, con); close(con)
  data.table::fwrite(body, file, sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE,
                     compress = if (gz) "gzip" else "none", showProgress = FALSE)
  invisible(file)
}
