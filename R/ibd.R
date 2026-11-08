## ibd.R -- ground-truth identity-by-descent from the tree sequence (tskit via reticulate), and its
## segments (D74).
## reticulate is a Suggests, and Python with tskit and numpy a system requirement; each function
## checks them through .require_python() (R/python.R), which says what is missing.

.trees_of <- function(run) {
  if (inherits(run, "vivax_run")) return(trees_path(run))
  if (is.character(run) && length(run) == 1) {
    if (dir.exists(run)) return(trees_path(run))
    return(run)  # assume a .trees path
  }
  stop("`run` must be a vivax_run, a run directory, or a .trees path.", call. = FALSE)
}

## The tree sequence of `run`, which must exist.
.trees_file <- function(run) {
  tf <- .trees_of(run)
  if (!file.exists(tf)) stop("tree sequence not found: ", tf, call. = FALSE)
  tf
}

## The IBD module: the tree sequence simplified once to the genomes asked for, then one
## ibd_segments() call over every pair among them, or, for listed pairs, one call per genome with
## its partners, or, for groups of genomes, one call over the pairs within each group. min_span is 0 throughout: every shared stretch counts in a pair's fraction, and
## the segments are floored in R once merged and cut (D74). Returns numpy arrays, keyed by
## pedigree id.
.ibd_py_src <- '
def mfw_ibd(trees, ids, pa, pb, want_segments, groups=None):
    ts = tskit.load(trees)
    md = ts.metadata
    end_tick = int(md["SLiM"]["tick"]) if isinstance(md, dict) and "SLiM" in md else -1
    p2n, _ = _ped_nodes(ts)
    listed = pa is not None
    grouped = groups is not None
    if listed:
        pa = [int(x) for x in np.atleast_1d(pa)]; pb = [int(x) for x in np.atleast_1d(pb)]
        req = sorted(set(pa) | set(pb))
    elif grouped:
        groups = [[int(x) for x in np.atleast_1d(g)] for g in groups]
        req = sorted({x for g in groups for x in g})
    elif ids is None:
        req = sorted(p2n)
    else:
        req = sorted({int(x) for x in np.atleast_1d(ids)})
    found = [p for p in req if p in p2n]
    out = dict(end_tick=end_tick, sequence_length=float(ts.sequence_length),
               missing=np.array([p for p in req if p not in p2n], dtype=np.int64),
               ped=np.zeros(0, np.int64), node_time=np.zeros(0), pa=np.zeros(0, np.int64),
               pb=np.zeros(0, np.int64), span=np.zeros(0), sa=np.zeros(0, np.int64),
               sb=np.zeros(0, np.int64), left=np.zeros(0), right=np.zeros(0),
               mrca=np.zeros(0, np.int64), mrca_time=np.zeros(0))
    if not found:
        return out
    sts = ts.simplify(samples=sorted(p2n[p] for p in found))
    del ts
    _, sped = _ped_nodes(sts)
    samp = np.asarray(sts.samples(), dtype=np.int64)
    out["ped"] = sped[samp]
    out["node_time"] = sts.nodes_time[samp]
    A = []; B = []; S = []; SA = []; SB = []; LE = []; RI = []; MN = []
    def take(res):
        for (a, b), sl in res.items():
            A.append(sped[a]); B.append(sped[b]); S.append(sl.total_span)
            if want_segments:
                n = len(sl)
                SA.append(np.full(n, sped[a])); SB.append(np.full(n, sped[b]))
                LE.append(np.asarray(sl.left)); RI.append(np.asarray(sl.right)); MN.append(np.asarray(sl.node))
    # store_segments stays off unless they are asked for: total_span per pair needs no segments,
    # and storing every one only to sum it costs several times as long
    if listed:
        node = {int(sped[n]): int(n) for n in samp}
        grp = {}
        for x, y in zip(pa, pb):
            if x != y and x in node and y in node:
                grp.setdefault(min(x, y), set()).add(max(x, y))
        for u in sorted(grp):
            take(sts.ibd_segments(between=[[node[u]], [node[v] for v in sorted(grp[u])]], min_span=0,
                                  store_pairs=True, store_segments=bool(want_segments)))
    elif grouped:
        node = {int(sped[n]): int(n) for n in samp}
        for g in groups:
            within = sorted({node[x] for x in g if x in node})
            if len(within) >= 2:
                take(sts.ibd_segments(within=within, min_span=0, store_pairs=True,
                                      store_segments=bool(want_segments)))
    else:
        take(sts.ibd_segments(within=list(samp), min_span=0, store_pairs=True,
                              store_segments=bool(want_segments)))
    out["pa"] = np.array(A, dtype=np.int64); out["pb"] = np.array(B, dtype=np.int64)
    out["span"] = np.array(S, dtype=float)
    if want_segments and MN:
        mn = np.concatenate(MN).astype(np.int64)
        out.update(sa=np.concatenate(SA).astype(np.int64), sb=np.concatenate(SB).astype(np.int64),
                   left=np.concatenate(LE).astype(float), right=np.concatenate(RI).astype(float),
                   mrca=sped[mn], mrca_time=sts.nodes_time[mn])
    return out
'

## Pedigree ids given as `name`, checked.
.ibd_ids <- function(x, name) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x != round(x)) || any(abs(x) > .Machine$integer.max))
    stop("`", name, "` must be pedigree ids: whole numbers, none missing.", call. = FALSE)
  as.integer(x)
}

## The tree sequence's IBD for every pair of `ids`, for the listed pairs (pa[i], pb[i]), or for the
## pairs within each of the `groups` (a list of id vectors), with the segments when asked: the work
## ibd_pairs(), ibd_segments() and new_infection_diversity() share. The segments are the same
## whatever other genomes are passed. With `none_ok`, none of the genomes being in the tree sequence
## gives an object with no genome rather than an error. `got` holds only the pairs that share a
## segment.
.ibd_compute <- function(run, ids = NULL, pa = NULL, pb = NULL, min_span = 0, segments = FALSE,
                         absent = "they are listed in `missing` and look up as NA.", none_ok = FALSE,
                         groups = NULL) {
  tf <- .trees_of(run)
  if (!file.exists(tf)) stop("tree sequence not found: ", tf, call. = FALSE)
  m <- .py_module("_malariafwd_ibd", .ibd_py_src)
  arg <- function(x) if (is.null(x)) NULL else as.list(x)
  r <- reticulate::py_to_r(m$mfw_ibd(tf, arg(ids), arg(pa), arg(pb), segments,
                                     if (is.null(groups)) NULL else unname(lapply(groups, as.list))))
  missing <- as.integer(r$missing)
  req <- if (!is.null(groups)) unique(unlist(groups)) else if (is.null(pa)) ids else unique(c(pa, pb))
  if (!none_ok && !is.null(req) && length(req) && length(missing) == length(unique(req)))
    stop("none of the ", length(unique(req)), " requested ids is in the tree sequence",
         .absent_reason(tf), ".", call. = FALSE)
  if (length(missing))
    warning(length(missing), " of the requested ids are not in the tree sequence (such as ",
            paste(utils::head(missing, 3), collapse = ", "), "): ", absent, call. = FALSE)
  end_tick <- as.integer(r$end_tick); L <- r$sequence_length
  smp <- data.table::data.table(pedigree_id = as.integer(r$ped), node_time = as.numeric(r$node_time))
  smp[, tick := as.integer(end_tick - node_time)]
  data.table::setorder(smp, pedigree_id)
  pa_ <- as.integer(r$pa); pb_ <- as.integer(r$pb)
  got <- data.table::data.table(a = pmin(pa_, pb_), b = pmax(pa_, pb_), ibd = as.numeric(r$span) / L)
  sg <- NULL
  if (segments) {
    sa <- as.integer(r$sa); sb <- as.integer(r$sb)
    raw <- data.table::data.table(a = pmin(sa, sb), b = pmax(sa, sb), left = as.numeric(r$left),
                                  right = as.numeric(r$right), mrca_id = as.integer(r$mrca),
                                  mrca_time = as.numeric(r$mrca_time))
    sg <- .ibd_tracts(raw, .ibd_layout(run, tf, L), min_span)
    tt <- stats::setNames(smp$node_time, smp$pedigree_id)
    sg[, depth := mrca_time - pmax(tt[as.character(a)], tt[as.character(b)])]
  }
  list(got = got, samples = smp, missing = missing, segments = sg, sequence_length = L,
       end_tick = end_tick, trees = tf)
}

## Why a tree sequence can hold none of the genomes asked for, as the params.tsv beside it records.
.absent_reason <- function(tf) {
  p <- tryCatch(.read_params_tsv(file.path(dirname(tf), "params.tsv")), error = function(e) NULL)
  if (!is.null(p$record$NEVER_STARTED)) return("; the run never started, so it holds no genome")
  if (isFALSE(p$params$RETAIN_FULL_HISTORY))
    return("; with RETAIN_FULL_HISTORY = FALSE most earlier genomes leave the tree sequence")
  ""
}

## tskit's segment rows, which are pieces of the paths through the trees and so depend on the
## genomes simplified to, as segments (D74): maximal stretches over which a pair descends from one
## common ancestor genome, cut at each chromosome junction, with the pieces shorter than
## `min_span` left out after the cut. The ancestor is the node, not its time: two ancestors born on
## one tick must not merge.
.ibd_tracts <- function(raw, layout, min_span) {
  raw <- data.table::copy(data.table::as.data.table(raw))
  empty <- data.table::data.table(a = integer(), b = integer(), chrom = integer(), start = numeric(),
                                  end = numeric(), gstart = numeric(), gend = numeric(), length = numeric(),
                                  mrca_id = integer(), mrca_time = numeric())
  if (!nrow(raw)) return(empty)
  data.table::setorder(raw, a, b, left)
  n <- nrow(raw)
  new <- c(TRUE, raw$a[-1L] != raw$a[-n] | raw$b[-1L] != raw$b[-n] | raw$mrca_id[-1L] != raw$mrca_id[-n] |
                 raw$left[-1L] != raw$right[-n])
  raw[, tract := cumsum(new)]
  tr <- raw[, .(a = a[1L], b = b[1L], left = left[1L], right = right[.N], mrca_id = mrca_id[1L],
                mrca_time = mrca_time[1L]), by = tract]
  cs <- layout$start; ce <- layout$end
  cl <- findInterval(tr$left, cs)                       # the chromosome a tract starts in
  cr <- findInterval(tr$right, cs, left.open = TRUE)    # ... and the one its last base is in
  k <- cr - cl + 1L
  i <- rep(seq_len(nrow(tr)), k)
  ch <- cl[i] + sequence(k) - 1L
  gs <- pmax(tr$left[i], cs[ch]); ge <- pmin(tr$right[i], ce[ch])
  out <- data.table::data.table(a = tr$a[i], b = tr$b[i], chrom = as.integer(layout$chrom[ch]),
                                start = gs - cs[ch], end = ge - cs[ch], gstart = gs, gend = ge,
                                length = ge - gs, mrca_id = tr$mrca_id[i], mrca_time = tr$mrca_time[i])
  out[length >= min_span]
}

## The chromosome layout of the run behind tree sequence `tf`, of `L` bp: chrom_map() of its run
## directory, or of the directory a .trees file sits in.
.ibd_layout <- function(run, tf, L) {
  dir <- if (inherits(run, "vivax_run")) run$outdir else if (dir.exists(run)) run else dirname(tf)
  if (!file.exists(file.path(dir, "params.tsv")))
    stop("IBD segments are cut at the chromosome junctions, which need the run's params.tsv; there is ",
         "none beside ", tf, ".", call. = FALSE)
  cm <- chrom_map(dir)
  if (max(cm$end) != L)
    stop("the run's chromosomes end at ", .slim_num(max(cm$end)), " bp, but its tree sequence is ",
         .slim_num(L), " bp long.", call. = FALSE)
  cm
}

## A vivax_ibd object, its parts checked.
.new_vivax_ibd <- function(pairs, samples, missing, segments, sequence_length, min_span, all_pairs,
                           end_tick, trees) {
  stopifnot(identical(names(pairs), c("a", "b", "ibd")), is.integer(pairs$a), is.integer(pairs$b),
            is.double(pairs$ibd), identical(names(samples), c("pedigree_id", "node_time", "tick")),
            is.integer(samples$pedigree_id), is.double(samples$node_time), is.integer(samples$tick),
            is.integer(missing), is.null(segments) || data.table::is.data.table(segments))
  data.table::setkey(pairs, a, b)
  structure(list(pairs = pairs, samples = samples, missing = missing, segments = segments,
                 sequence_length = sequence_length, min_span = min_span, all_pairs = all_pairs,
                 end_tick = end_tick, trees = trees), class = "vivax_ibd")
}

#' Ground-truth IBD from a run's tree sequence
#'
#' Computes genome-wide identity-by-descent (shared-segment length / sequence length) among a set
#' of genomes, keyed by pedigree id, from the tskit tree sequence recorded by the model. The tree
#' is simplified to the genomes asked for, and SLiM's vacant nodes are left out.
#'
#' @section What is counted:
#' A pair's `ibd` counts **all** shared ancestry: every stretch over which the two genomes descend
#' from a common ancestor genome, of any age and any length, over the sequence length. `min_span`
#' does not change it. To impose a cutoff, ask for `segments` and filter them on their length or
#' their age (see *Segment ages*).
#' @section Segments:
#' A segment is a maximal stretch over which two genomes descend from one common ancestor genome,
#' cut at each chromosome junction: a change of ancestor starts a new segment even with no gap
#' between them. A pair gets the same segments whatever other genomes are passed, the same as
#' [ibd_segments()] gives, and at `min_span = 0` they add up to its `ibd`.
#'
#' @param run a run: a `vivax_run` (from [run_sim()] or [load_run()]) or the path of a run directory,
#'   or of a `.trees` file in its run directory, whose `params.tsv` gives the chromosome layout the
#'   segments need.
#' @param ids the pedigree ids of the genomes to compare, every pair among them; `NULL` for every
#'   genome the tree sequence records, which warns when that is more than 1e6 pairs.
#' @param pairs instead of `ids`, a data frame of pairs to compare, with the pedigree ids in columns
#'   `a` and `b`: only these pairs are computed, about 0.2 ms each, so use it for a few pairs among
#'   many genomes, such as each recurrence against its host's earlier episodes.
#' @param min_span the shortest segment returned, in bp, after the cut at the junctions: a segment
#'   of exactly `min_span` is kept. It does not change `ibd`.
#' @param segments if `TRUE`, also return every shared segment with its common ancestor and that
#'   ancestor's age.
#' @return a `vivax_ibd` object, a list of:
#'   * `pairs`, a data.table keyed by `a` and `b`: one row for every pair of distinct genomes
#'     found, `a < b` (integer pedigree ids), and `ibd` (double, in \[0, 1\]), 0 for a pair that
#'     shares nothing;
#'   * `samples`, one row per genome found, by pedigree id: `pedigree_id`, `node_time` (ticks
#'     before the end of the run) and `tick` (the day the genome arose, as in `events.tsv`);
#'   * `missing`, the requested ids the tree sequence does not hold;
#'   * `segments`, `NULL` or a data.table, by pair and position: `a`, `b`, `chrom` (the PvP01
#'     chromosome), `start` and `end` (bp within the chromosome, `end` exclusive), `gstart` and
#'     `gend` (bp along the concatenated genome), `length` (bp), `mrca_id` (the pedigree id of the
#'     common ancestor), `mrca_time` and `depth` (its age, see *Segment ages*);
#'   * `sequence_length`, `min_span`, `all_pairs` (`TRUE`: every pair of the genomes is in
#'     `pairs`; `FALSE` for listed `pairs`, which hold each listed pair of two distinct genomes
#'     found, once, with `a < b`), `end_tick` (the run's last tick) and `trees` (the file read).
#'
#'   Read a pair's value with [ibd_lookup()], or the table with `as.matrix()` or `as.data.frame()`.
#' @section Ids that are not in the tree sequence:
#' A run written with `RETAIN_FULL_HISTORY = FALSE` keeps only a fraction of genomes in the tree
#' sequence. A requested id it does not hold is listed in `missing`, gives one warning, and looks
#' up as `NA`, not as 0; asking only for such ids is an error.
#' @section Segment ages:
#' `mrca_time` is a tskit node time, ticks before the end of the run, not an age. A segment's
#' `depth` is how long before the older of the two genomes its common ancestor lived,
#' `mrca_time - max(node_time[a], node_time[b])`, in ticks. A recency filter is written against it, for
#' example `x$segments[depth <= 180]` for ancestry within about six months: filtering on
#' `mrca_time` alone scores a pair by when its ancestry last branched, not by how related the two
#' genomes are.
#' @section Cost:
#' One simplification and one pass over the pairs; `segments = TRUE` costs a few times as long and
#' gives one row per segment per pair, so restrict it with `ids` and `min_span`.
#' @examples
#' \dontrun{
#' x <- ibd_pairs(run, ids = classify_cases(load_events(run))$pedigree_id)
#' ibd_lookup(x, 12, 345)
#' }
#' @seealso [ibd_lookup()], [ibd_segments()]
#' @inheritSection tskit_available Python
#' @export
ibd_pairs <- function(run, ids = NULL, pairs = NULL, min_span = 1000, segments = FALSE)
  .ibd_pairs(run, ids, pairs, min_span, segments)

## ibd_pairs(); with `none_ok`, an object with no genome when none of the ids is in the tree
## sequence, for the readers that then give their NA row or empty table.
.ibd_pairs <- function(run, ids = NULL, pairs = NULL, min_span = 1000, segments = FALSE, none_ok = FALSE) {
  if (!is.null(ids) && !is.null(pairs)) stop("give `ids` or `pairs`, not both.", call. = FALSE)
  if (!is.null(ids)) ids <- sort(unique(.ibd_ids(ids, "ids")))
  if (!is.null(pairs)) {
    if (!is.data.frame(pairs) || !all(c("a", "b") %in% names(pairs)))
      stop("`pairs` must be a data frame of pedigree ids in columns `a` and `b`.", call. = FALSE)
    pa <- .ibd_ids(pairs$a, "pairs$a"); pb <- .ibd_ids(pairs$b, "pairs$b")
  }
  .ibd_check_min_span(min_span)
  .check_flag(segments, "segments")
  .require_python("ibd_pairs()", c("tskit", "numpy"))
  if (is.null(pairs)) {
    x <- .ibd_compute(run, ids = ids, min_span = min_span, segments = segments, none_ok = none_ok)
    p <- x$samples$pedigree_id; n <- length(p)
    if (is.null(ids) && n * (n - 1) / 2 > 1e6)
      warning("ibd_pairs() of every genome gave ", format(n * (n - 1) / 2, big.mark = ","), " pairs; pass ",
              "`ids` to compare only the genomes you need.", call. = FALSE)
    all <- if (n >= 2) data.table::data.table(a = rep(p[-n], (n - 1):1), b = p[sequence((n - 1):1, from = 2:n)])
           else data.table::data.table(a = integer(), b = integer())
  } else {
    x <- .ibd_compute(run, pa = pa, pb = pb, min_span = min_span, segments = segments, none_ok = none_ok)
    p <- x$samples$pedigree_id
    all <- unique(data.table::data.table(a = pmin(pa, pb), b = pmax(pa, pb)))[a != b & a %in% p & b %in% p]
    data.table::setorder(all, a, b)
  }
  pt <- merge(all, x$got, by = c("a", "b"), all.x = TRUE)
  pt[is.na(ibd), ibd := 0]
  .new_vivax_ibd(pt, x$samples, x$missing, x$segments, x$sequence_length, min_span, is.null(pairs), x$end_tick,
                 x$trees)
}

.ibd_check_min_span <- function(min_span)
  if (!is.numeric(min_span) || length(min_span) != 1L || is.na(min_span) || min_span < 0)
    stop("`min_span` must be a number of bp, 0 or more.", call. = FALSE)

#' The IBD of pairs of genomes, from ibd_pairs()
#'
#' `ibd_lookup()` gives the IBD of each pair (`a[i]`, `b[i]`): 1 for a genome with itself, the
#' pair's value otherwise, and `NA` when either genome is not in the object or, for an object of
#' listed pairs (`ibd_pairs(pairs =)`), when the pair was not listed. A length-1 argument is
#' recycled against the other. `as.matrix()` gives the pairs as a symmetric matrix named by
#' pedigree id, with a unit diagonal, and `as.data.frame()` the `pairs` table.
#' @param x a `vivax_ibd` object, from [ibd_pairs()].
#' @param a,b pedigree ids.
#' @param ids for `as.matrix()`, the genomes to give, in order; `NULL` for every one.
#' @param ... ignored.
#' @return `ibd_lookup()`: a numeric vector; `as.matrix()`: a matrix; `as.data.frame()`: a
#'   data.frame of `a`, `b` and `ibd`.
#' @examples
#' \dontrun{
#' x <- ibd_pairs(run, ids = c(12, 345, 678))
#' ibd_lookup(x, 12, c(345, 678))
#' as.matrix(x)
#' }
#' @seealso [ibd_pairs()]
#' @export
ibd_lookup <- function(x, a, b) {
  if (!inherits(x, "vivax_ibd")) stop("`x` must be a vivax_ibd object, from ibd_pairs().", call. = FALSE)
  a <- as.integer(a); b <- as.integer(b)
  if (length(a) != length(b)) {
    if (length(a) == 1L) a <- rep(a, length(b))
    else if (length(b) == 1L) b <- rep(b, length(a))
    else stop("`a` and `b` must have the same length, or one of them length 1.", call. = FALSE)
  }
  lo <- pmin(a, b); hi <- pmax(a, b)
  v <- x$pairs[data.table::data.table(a = lo, b = hi), on = c("a", "b")]$ibd
  known <- a %in% x$samples$pedigree_id & b %in% x$samples$pedigree_id
  v[known & a == b] <- 1
  v[!known] <- NA_real_
  unname(v)
}

#' @rdname ibd_lookup
#' @export
as.matrix.vivax_ibd <- function(x, ids = NULL, ...) {
  ids <- if (is.null(ids)) x$samples$pedigree_id else as.integer(ids)
  g <- expand.grid(i = seq_along(ids), j = seq_along(ids))
  matrix(ibd_lookup(x, ids[g$i], ids[g$j]), length(ids), length(ids),
         dimnames = list(as.character(ids), as.character(ids)))
}

#' @rdname ibd_lookup
#' @param row.names,optional ignored.
#' @export
as.data.frame.vivax_ibd <- function(x, row.names = NULL, optional = FALSE, ...)
  as.data.frame(x$pairs)

#' @export
print.vivax_ibd <- function(x, ...) {
  cat("<vivax_ibd> ", nrow(x$samples), " genomes, ", nrow(x$pairs), if (isTRUE(x$all_pairs)) " pairs" else " listed pairs",
      ", min_span ",
      .slim_num(x$min_span), " bp\n", sep = "")
  if (length(x$missing)) cat("  ", length(x$missing), " requested ids not in the tree sequence\n", sep = "")
  cat("  segments: ", if (is.null(x$segments)) "not computed" else paste(nrow(x$segments), "rows"), "\n", sep = "")
  cat("  from ", x$trees, ", ending on tick ", x$end_tick, "\n", sep = "")
  invisible(x)
}
