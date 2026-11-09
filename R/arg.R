## arg.R -- ARG-native (tree) summaries of the recorded genealogy.
##
## The recorded genealogy is a FOREST, not a single tree: each founder lineage and each import is
## an independent root, so most pairs of genomes never coalesce within the simulated period. On a
## forest, tskit's branch statistics are not a relatedness measure: for a pair with no common
## ancestor they accumulate branch length up to each genome's own root, so pairs that coalesce and
## pairs that never do score alike. arg_tmrca() is the relatedness measure here: the time to the
## common ancestor where one exists, with `coalesced`, the fraction of the genome over which it
## does. arg_summary() shows the forest itself.

## The ARG module. Each function loads the tree sequence, simplifies it to the genomes asked for
## (every genome when `ids` is None), and says how many genomes the run recorded and how many it kept,
## so that R can give a run with no genome its empty result. Ids are passed as arguments.
.arg_py_src <- '
def _arg_load(trees, ids):
    ts = tskit.load(trees)
    p2n, _ = _ped_nodes(ts)
    keep = sorted(p2n.values()) if ids is None else sorted({p2n[int(p)] for p in np.atleast_1d(ids) if int(p) in p2n})
    sts = ts.simplify(samples=keep) if keep else None
    return ts, sts, len(p2n), len(keep)

def arg_summary(trees, ids):
    ts, sts, genomes, kept = _arg_load(trees, ids)
    out = dict(genomes=genomes, kept=kept, num_trees=int(ts.num_trees), num_nodes=int(ts.num_nodes),
               num_edges=int(ts.num_edges), num_samples=0, seq_length_mb=float(ts.sequence_length) / 1e6)
    if sts is None:
        return out
    nr = np.array([t.num_roots for t in sts.trees()])
    out.update(num_trees=int(sts.num_trees), num_nodes=int(sts.num_nodes), num_edges=int(sts.num_edges),
               num_samples=int(sts.num_samples), seq_length_mb=float(sts.sequence_length) / 1e6,
               roots_min=int(nr.min()), roots_median=float(np.median(nr)), roots_max=int(nr.max()),
               pct_trees_coalesced=float(100 * np.mean(nr == 1)), max_root_time=float(sts.max_root_time))
    return out

def arg_tmrca(trees, ids_a, ids_b):
    ids = sorted({int(x) for x in list(np.atleast_1d(ids_a)) + list(np.atleast_1d(ids_b))})
    ts, sts, genomes, kept = _arg_load(trees, ids)
    n = len(np.atleast_1d(ids_a))
    co = np.full(n, np.nan); tm = np.full(n, np.nan); nt = np.full(n, np.nan)
    if sts is not None:
        p2, _ = _ped_nodes(sts)
        for i, (x, y) in enumerate(zip(np.atleast_1d(ids_a), np.atleast_1d(ids_b))):
            a, b = p2.get(int(x)), p2.get(int(y))
            if a is None or b is None:
                continue
            if a == b:
                co[i] = 1.0; tm[i] = 0.0; nt[i] = float(sts.node(a).time); continue
            # t.tmrca() gives the ancestor node time, ticks before the end of the run; less the
            # older genome own node time, a depth, which is never negative
            base = max(sts.node(a).time, sts.node(b).time)
            tot = 0.0; c = 0.0; acc = 0.0
            for t in sts.trees():
                span = t.interval.right - t.interval.left
                tot += span
                try:
                    m = t.tmrca(a, b)
                except ValueError:
                    continue
                c += span; acc += span * m
            co[i] = c / tot if tot else np.nan
            if c > 0:
                tm[i] = acc / c - base; nt[i] = acc / c
    return dict(genomes=genomes, kept=kept, coalesced=co, tmrca=tm, tmrca_node_time=nt)
'

## A NaN from Python as NA.
.na_nan <- function(x) { x <- as.numeric(x); x[is.nan(x)] <- NA_real_; x }

## Stops when `ids` were given and none of them is a genome in a run that recorded some.
.arg_none_found <- function(r, ids, what)
  if (!is.null(ids) && length(ids) && r$genomes > 0 && r$kept == 0)
    stop(what, ": none of the ", length(unique(ids)), " ids given is a genome in this run's tree sequence (a run ",
         "with RETAIN_FULL_HISTORY = FALSE keeps only some genomes).", call. = FALSE)

#' Structure of the recorded ARG
#'
#' Whole-genealogy summary: how many local trees and edges, how deep it goes, and -- the key
#' diagnostic for a forward run -- how many ROOTS each local tree has. A fully coalesced
#' genealogy has one root; this model seeds independent founder lineages, so the expected
#' result is a forest with roughly one root per surviving founder lineage.
#'
#' @param run a run: a `vivax_run` (from [run_sim()] or [load_run()]) or the path of a run directory,
#'   or of a `.trees` file.
#' @param ids optional pedigree ids to restrict to (default: all recorded genomes); ids that are not
#'   genomes in the tree sequence are left out, and none being one is an error.
#' @return a one-row data.table: `num_trees`, `num_nodes`, `num_edges`, `num_samples`,
#'   `seq_length_mb`, `roots_min`/`roots_median`/`roots_max`, `pct_trees_coalesced`,
#'   `max_root_time` (ticks before the end of the run). A run that recorded no genome, such as
#'   one that never started, gives the tree sequence's counts, no samples, 0 roots, and `NA` for
#'   `pct_trees_coalesced` and `max_root_time`.
#' @examplesIf tskit_available()
#' arg_summary(example_run())
#' @inheritSection tskit_available Python
#' @export
arg_summary <- function(run, ids = NULL) {
  if (!is.null(ids)) ids <- .ibd_ids(ids, "ids")
  .require_python("arg_summary()", c("tskit", "numpy"))
  m <- .py_module("_malariafwd_arg", .arg_py_src)
  r <- reticulate::py_to_r(m$arg_summary(.trees_file(run), if (is.null(ids)) NULL else as.list(unique(ids))))
  .arg_none_found(r, ids, "arg_summary()")
  got <- !is.null(r$roots_min)
  data.table::data.table(num_trees = as.integer(r$num_trees), num_nodes = as.integer(r$num_nodes),
                         num_edges = as.integer(r$num_edges), num_samples = as.integer(r$num_samples),
                         seq_length_mb = as.numeric(r$seq_length_mb),
                         roots_min = if (got) as.integer(r$roots_min) else 0L,
                         roots_median = if (got) as.numeric(r$roots_median) else 0,
                         roots_max = if (got) as.integer(r$roots_max) else 0L,
                         pct_trees_coalesced = if (got) .na_nan(r$pct_trees_coalesced) else NA_real_,
                         max_root_time = if (got) .na_nan(r$max_root_time) else NA_real_)
}

#' Pairwise time to the most recent common ancestor (the ARG-native relatedness measure)
#'
#' For each pair, the span-weighted mean TMRCA **over the part of the genome where a common
#' ancestor exists**, plus `coalesced`, the fraction of the genome where it does. This is the
#' time-domain complement to IBD: IBD asks *how much* of the genome is shared,
#' `arg_tmrca()` asks *how far back* the sharing goes. Unlike tskit's branch-mode divergence it is
#' meaningful on a forest, because non-coalescing spans are excluded rather than scored as 0.
#'
#' Clonal pairs give `coalesced` = 1 with a very recent TMRCA; unrelated genomes give
#' `coalesced` = 0 and `tmrca` = `NA`.
#'
#' @param run a run: a `vivax_run` (from [run_sim()] or [load_run()]) or the path of a run directory,
#'   or of a `.trees` file.
#' @param ids_a,ids_b equal-length vectors of pedigree ids, compared pairwise.
#' @return a data.table, one row per pair: `a`, `b`, `coalesced` (0-1), `tmrca` and
#'   `tmrca_node_time` (both ticks; `NA` if the pair never coalesces or a genome is not in the tree
#'   sequence). A run that recorded no genome gives an `NA` row per pair, and none of the ids being a
#'   genome of a run that recorded some is an error.
#' @section Which TMRCA you get:
#' `tmrca` is a **depth**: ticks from the older of the two sampled genomes back to their
#' common ancestor, so an identical relationship returns an identical number wherever in the
#' run it occurred -- a clonal parent-child pair gives 0 whether it happened on tick 10 or
#' tick 390. `tmrca_node_time` is the raw tskit node time of that ancestor, i.e. ticks before
#' the **end of the run**; that is what `tree.tmrca()` returns natively and it is dominated by
#' when the pair was sampled rather than by how related it is. Use `tmrca` for relatedness.
#' @examplesIf tskit_available()
#' run <- example_run()
#' rel <- load_events(run)[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][1:3]
#' arg_tmrca(run, rel$pedigree_id, rel$parent1_pedigree_id)   # relapses vs their batches
#' @inheritSection tskit_available Python
#' @export
arg_tmrca <- function(run, ids_a, ids_b) {
  ids_a <- .ibd_ids(ids_a, "ids_a"); ids_b <- .ibd_ids(ids_b, "ids_b")
  if (length(ids_a) != length(ids_b))
    stop("`ids_a` and `ids_b` must have the same length, one pair per position; they have ", length(ids_a),
         " and ", length(ids_b), ".", call. = FALSE)
  .require_python("arg_tmrca()", c("tskit", "numpy"))
  n <- length(ids_a)
  out <- data.table::data.table(a = ids_a, b = ids_b, coalesced = rep(NA_real_, n), tmrca = rep(NA_real_, n),
                                tmrca_node_time = rep(NA_real_, n))
  if (!n) return(out)
  m <- .py_module("_malariafwd_arg", .arg_py_src)
  r <- reticulate::py_to_r(m$arg_tmrca(.trees_file(run), as.list(ids_a), as.list(ids_b)))
  .arg_none_found(r, c(ids_a, ids_b), "arg_tmrca()")
  out[, `:=`(coalesced = .na_nan(r$coalesced), tmrca = .na_nan(r$tmrca), tmrca_node_time = .na_nan(r$tmrca_node_time))]
  out[]
}
