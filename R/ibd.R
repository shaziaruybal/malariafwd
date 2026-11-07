## ibd.R -- ground-truth identity-by-descent from the tree sequence (tskit via reticulate).
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

## The IBD module: the tree sequence simplified once to the requested genomes, and one
## ibd_segments() call over every pair among them. Returns numpy arrays, keyed by pedigree id.
.ibd_py_src <- '
def mfw_ibd(trees, ids, min_span, want_segments):
    ts = tskit.load(trees)
    md = ts.metadata
    end_tick = int(md["SLiM"]["tick"]) if isinstance(md, dict) and "SLiM" in md else -1
    p2n, _ = _ped_nodes(ts)
    req = sorted(p2n) if ids is None else sorted({int(x) for x in np.atleast_1d(ids)})
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
    # store_segments stays off unless they are asked for: total_span per pair needs no segments,
    # and storing every one only to sum it costs several times as long
    res = sts.ibd_segments(within=list(samp), min_span=min_span, store_pairs=True,
                           store_segments=bool(want_segments))
    A = []; B = []; S = []; SA = []; SB = []; LE = []; RI = []; MN = []
    for (a, b), sl in res.items():
        A.append(sped[a]); B.append(sped[b]); S.append(sl.total_span)
        if want_segments:
            n = len(sl)
            SA.append(np.full(n, sped[a])); SB.append(np.full(n, sped[b]))
            LE.append(np.asarray(sl.left)); RI.append(np.asarray(sl.right)); MN.append(np.asarray(sl.node))
    out["pa"] = np.array(A, dtype=np.int64); out["pb"] = np.array(B, dtype=np.int64)
    out["span"] = np.array(S, dtype=float)
    if want_segments and MN:
        mn = np.concatenate(MN).astype(np.int64)
        out.update(sa=np.concatenate(SA).astype(np.int64), sb=np.concatenate(SB).astype(np.int64),
                   left=np.concatenate(LE).astype(float), right=np.concatenate(RI).astype(float),
                   mrca=sped[mn], mrca_time=sts.nodes_time[mn])
    return out
'

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
#' is simplified to the genomes asked for, and SLiM's vacant nodes are left out. **All** shared
#' ancestry is counted, however far back it goes; to impose a cutoff, ask for `segments` and filter
#' them on their age (see *Segment ages*).
#'
#' @param run a `vivax_run`, a run directory, or a `.trees` file path.
#' @param ids the pedigree ids of the genomes to compare, every pair among them; `NULL` for every
#'   genome the tree sequence records, which warns when that is more than 1e6 pairs.
#' @param min_span minimum shared-segment length (bp) to count (default 1000).
#' @param segments if `TRUE`, also return every shared segment with its common ancestor and that
#'   ancestor's age.
#' @return a `vivax_ibd` object, a list of:
#'   * `pairs`, a data.table keyed by `a` and `b`: one row for every pair of distinct genomes
#'     found, `a < b` (integer pedigree ids), and `ibd` (double, in \[0, 1\]), 0 for a pair that
#'     shares nothing;
#'   * `samples`, one row per genome found, by pedigree id: `pedigree_id`, `node_time` (ticks
#'     before the end of the run) and `tick` (the day the genome arose, as in `events.tsv`);
#'   * `missing`, the requested ids the tree sequence does not hold;
#'   * `segments`, `NULL` or a data.table of `a`, `b`, `left`, `right` (bp), `mrca_id` (the pedigree
#'     id of the common ancestor), `mrca_time` and `depth` (its age, see *Segment ages*);
#'   * `sequence_length`, `min_span`, `all_pairs` (`TRUE`: every pair of the genomes is in
#'     `pairs`), `end_tick` (the run's last tick) and `trees` (the file read).
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
#' `mrca_time` alone scores a pair by when its lineage last branched, not by how related the two
#' genomes are.
#' @section Cost:
#' One simplification and one pass over the pairs; `segments = TRUE` costs several times as long
#' and gives one row per segment per pair, so restrict it with `ids` and `min_span`.
#' @examples
#' \dontrun{
#' x <- ibd_pairs(run, ids = classify_cases(load_events(run))$pedigree_id)
#' ibd_lookup(x, 12, 345)
#' }
#' @seealso [ibd_lookup()], [ibd_segments()]
#' @inheritSection tskit_available Python
#' @export
ibd_pairs <- function(run, ids = NULL, min_span = 1000, segments = FALSE) {
  if (!is.null(ids)) {
    if (!is.numeric(ids) || anyNA(ids) || any(!is.finite(ids)) || any(ids != round(ids)) ||
        any(abs(ids) > .Machine$integer.max))
      stop("`ids` must be pedigree ids: whole numbers, none missing.", call. = FALSE)
    ids <- sort(unique(as.integer(ids)))
  }
  if (!is.numeric(min_span) || length(min_span) != 1L || is.na(min_span) || min_span < 0)
    stop("`min_span` must be a number of bp, 0 or more.", call. = FALSE)
  .check_flag(segments, "segments")
  .require_python("ibd_pairs()", c("tskit", "numpy"))
  tf <- .trees_of(run)
  if (!file.exists(tf)) stop("tree sequence not found: ", tf, call. = FALSE)
  m <- .py_module("_malariafwd_ibd", .ibd_py_src)
  r <- reticulate::py_to_r(m$mfw_ibd(tf, if (is.null(ids)) NULL else as.list(ids), as.numeric(min_span), segments))
  missing <- as.integer(r$missing)
  if (!is.null(ids) && length(ids) && length(missing) == length(ids))
    stop("none of the ", length(ids), " requested ids is in the tree sequence; a run with ",
         "RETAIN_FULL_HISTORY = FALSE keeps only a fraction of genomes.", call. = FALSE)
  if (length(missing))
    warning(length(missing), " of the requested ids are not in the tree sequence (such as ",
            paste(utils::head(missing, 3), collapse = ", "), "): they are listed in `missing` and look ",
            "up as NA.", call. = FALSE)
  end_tick <- as.integer(r$end_tick)
  smp <- data.table::data.table(pedigree_id = as.integer(r$ped), node_time = as.numeric(r$node_time))
  smp[, tick := as.integer(end_tick - node_time)]
  data.table::setorder(smp, pedigree_id)
  p <- smp$pedigree_id; n <- length(p)
  if (is.null(ids) && n * (n - 1) / 2 > 1e6)
    warning("ibd_pairs() of every genome gave ", format(n * (n - 1) / 2, big.mark = ","), " pairs; pass ",
            "`ids` to compare only the genomes you need.", call. = FALSE)
  all <- if (n >= 2) data.table::data.table(a = rep(p[-n], (n - 1):1), b = p[sequence((n - 1):1, from = 2:n)])
         else data.table::data.table(a = integer(), b = integer())
  pa <- as.integer(r$pa); pb <- as.integer(r$pb)
  got <- data.table::data.table(a = pmin(pa, pb), b = pmax(pa, pb), ibd = as.numeric(r$span) / r$sequence_length)
  pairs <- merge(all, got, by = c("a", "b"), all.x = TRUE)
  pairs[is.na(ibd), ibd := 0]
  sg <- NULL
  if (segments) {
    sa <- as.integer(r$sa); sb <- as.integer(r$sb)
    sg <- data.table::data.table(a = pmin(sa, sb), b = pmax(sa, sb), left = as.numeric(r$left),
                                 right = as.numeric(r$right), mrca_id = as.integer(r$mrca),
                                 mrca_time = as.numeric(r$mrca_time))
    tt <- stats::setNames(smp$node_time, smp$pedigree_id)
    sg[, depth := mrca_time - pmax(tt[as.character(a)], tt[as.character(b)])]
    data.table::setorder(sg, a, b, left)
  }
  .new_vivax_ibd(pairs, smp, missing, sg, r$sequence_length, min_span, TRUE, end_tick, tf)
}

#' The IBD of pairs of genomes, from ibd_pairs()
#'
#' `ibd_lookup()` gives the IBD of each pair (`a[i]`, `b[i]`): 1 for a genome with itself, the
#' pair's value otherwise, and `NA` when either genome is not in the object. A length-1 argument is
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
  cat("<vivax_ibd> ", nrow(x$samples), " genomes, ", nrow(x$pairs), " pairs, min_span ",
      .slim_num(x$min_span), " bp\n", sep = "")
  if (length(x$missing)) cat("  ", length(x$missing), " requested ids not in the tree sequence\n", sep = "")
  cat("  segments: ", if (is.null(x$segments)) "not computed" else paste(nrow(x$segments), "rows"), "\n", sep = "")
  cat("  from ", x$trees, ", ending on tick ", x$end_tick, "\n", sep = "")
  invisible(x)
}
