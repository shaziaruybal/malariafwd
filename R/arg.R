## arg.R -- ARG-native (tree/branch) summaries of the recorded genealogy.
##
## The model is NEUTRAL (mutation rate 0), so every SITE-based statistic is identically zero.
## tskit's BRANCH mode is the right tool: it measures the statistic on the genealogy itself
## (in units of time), and equals the expected site value under any mutation rate. It is also
## the natural complement to IBD -- IBD asks "how much of the genome is shared", branch mode
## asks "how deep in time is the shared ancestry".
##
## IMPORTANT, and specific to a forward simulation like this one: the recorded genealogy is a
## FOREST, not a single tree. Each founder lineage (and each import) is an independent genome
## root, so most sample pairs never coalesce within the simulated period. tskit reports
## branch statistics that do NOT behave as relatedness measures. Measured on a real run:
## pairs that DO coalesce had mean branch divergence 29.1, pairs that never coalesce had 32.1
## -- indistinguishable, because for a disconnected pair the statistic accumulates branch
## length up to each sample's OWN root. So:
##   * do NOT read branch divergence as relatedness on a forest; use [arg_tmrca()], which
##     computes the time to the common ancestor explicitly and reports the fraction of the
##     genome over which one exists;
##   * every function here reports `coalesced` so the dilution is visible.

.arg_prelude <- '
import tskit, json, numpy as np
def _md(x):
    md = x.metadata
    if isinstance(md, bytes):
        try: md = json.loads(md.decode()) if md else {}
        except Exception: md = {}
    return md if isinstance(md, dict) else {}
def _vac(nd):
    v = _md(nd).get("is_vacant"); return bool(v[0]) if v else False
def _ped2node(ts):
    d = {}
    for n in ts.samples():
        nd = ts.node(n)
        if nd.individual == tskit.NULL or _vac(nd): continue
        m = _md(ts.individual(nd.individual))
        p = m.get("pedigree_id", m.get("pedigreeID"))
        if p is not None: d.setdefault(int(p), int(n))
    return d
def _aslist(x):
    # a length-1 R vector arrives as a python scalar, not a list
    if x is None: return None
    try: return [int(v) for v in x]
    except TypeError: return [int(x)]
def _coalesced_frac(sts, A, B=None):
    """Span-weighted fraction of PAIRS that share a common ancestor.

    The genealogy is a forest, so many pairs never coalesce; tskit scores those as 0
    divergence, which is why this fraction has to be reported next to any branch statistic.
    Within a local tree, two samples coalesce iff they sit under the same root, so the
    coalescing pairs are counted per root and span-weighted across trees.
    """
    A = set(A)
    within = B is None
    B = A if within else set(B)
    total = (len(A) * (len(A) - 1) / 2) if within else (len(A) * len(B))
    if total == 0: return float("nan")
    tot = 0.0; acc = 0.0
    for t in sts.trees():
        span = t.interval.right - t.interval.left
        tot += span
        c = 0
        for r in t.roots:
            ss = set(t.samples(r))
            a = len(ss & A)
            if within: c += a * (a - 1) / 2
            else:      c += a * len(ss & B)
        acc += span * (c / total)
    return acc / tot if tot else float("nan")
'

.arg_load <- function(run, ids = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("ARG summaries need the 'reticulate' package.", call. = FALSE)
  if ("py_require" %in% getNamespaceExports("reticulate"))
    try(reticulate::py_require(c("tskit", "numpy")), silent = TRUE)
  ok <- tryCatch({ reticulate::py_run_string("import tskit, numpy"); TRUE }, error = function(e) FALSE)
  if (!ok) stop("ARG summaries need Python with 'tskit' + 'numpy' ",
                "(e.g. reticulate::py_require(c('tskit','numpy'))).", call. = FALSE)
  reticulate::py_run_string(.arg_prelude)
  ## NB bind the python main module locally: `reticulate::py$x <- v` cannot assign into a
  ## namespace-qualified object.
  pymain <- reticulate::py
  pymain$`_trees` <- trees_path(run)
  pymain$`_ids` <- if (is.null(ids)) NULL else as.integer(unique(ids))
  reticulate::py_run_string('
ts = tskit.load(_trees)
p2n = _ped2node(ts)
if _ids is None:
    keep = sorted(p2n.values())
else:
    keep = sorted({p2n[int(p)] for p in _aslist(_ids) if int(p) in p2n})
sts = ts.simplify(samples=keep) if keep else None
n2p = {}
if sts is not None:
    q = _ped2node(sts)
    n2p = {v: k for k, v in q.items()}
')
  invisible(TRUE)
}

#' Structure of the recorded ARG
#'
#' Whole-genealogy summary: how many local trees and edges, how deep it goes, and -- the key
#' diagnostic for a forward run -- how many ROOTS each local tree has. A fully coalesced
#' genealogy has one root; this model seeds independent founder lineages, so the expected
#' result is a forest with roughly one root per surviving founder lineage.
#'
#' @param run a run directory or `vivax_run`.
#' @param ids optional pedigree ids to restrict to (default: all recorded genomes).
#' @return a one-row data.table: `num_trees`, `num_nodes`, `num_edges`, `num_samples`,
#'   `seq_length_mb`, `roots_min`/`roots_median`/`roots_max`, `pct_trees_coalesced`,
#'   `max_root_time` (ticks before the end of the run).
#' @examples
#' \dontrun{arg_summary(run)}
#' @export
arg_summary <- function(run, ids = NULL) {
  .arg_load(run, ids)
  reticulate::py_run_string('
nr = np.array([t.num_roots for t in sts.trees()])
_out = dict(num_trees=int(sts.num_trees), num_nodes=int(sts.num_nodes),
            num_edges=int(sts.num_edges), num_samples=int(sts.num_samples),
            seq_length_mb=float(sts.sequence_length)/1e6,
            roots_min=int(nr.min()), roots_median=float(np.median(nr)), roots_max=int(nr.max()),
            pct_trees_coalesced=float(100*np.mean(nr == 1)),
            max_root_time=float(sts.max_root_time))
')
  data.table::as.data.table(reticulate::py$`_out`)
}

#' Branch-mode genetic diversity of the genealogy
#'
#' Mutation-free diversity from the ARG: `diversity` is the mean pairwise branch length
#' separating two samples (= expected nucleotide diversity per unit mutation rate), and is
#' reported with `tajimas_d` and branch-mode `segregating_sites`. Because the genealogy is a
#' forest (see [arg_summary()]), these are diluted by pairs that never coalesce -- use
#' [arg_divergence()] if you need the coalescing fraction made explicit.
#'
#' @param run a run directory or `vivax_run`.
#' @param ids optional pedigree ids to restrict to.
#' @return a one-row data.table of `n_samples`, `diversity`, `segregating_sites`, `tajimas_d`,
#'   and `coalesced` (the fraction of sample pairs that share an ancestor at all).
#' @examples
#' \dontrun{arg_diversity(run, ids = classify_cases(load_events(run))$pedigree_id)}
#' @export
arg_diversity <- function(run, ids = NULL) {
  .arg_load(run, ids)
  reticulate::py_run_string('
S = list(sts.samples())
_out = dict(n_samples=int(sts.num_samples),
            diversity=float(sts.diversity(mode="branch")),
            segregating_sites=float(sts.segregating_sites(mode="branch")),
            tajimas_d=float(sts.Tajimas_D(mode="branch")),
            coalesced=float(_coalesced_frac(sts, S)))
')
  data.table::as.data.table(reticulate::py$`_out`)
}

#' Branch-mode divergence between two sets of genomes, with the coalescing fraction
#'
#' @section Warning:
#' On a non-coalesced genealogy (this model, see [arg_summary()]) branch divergence is **not**
#' a relatedness measure: for a pair with no common ancestor it accumulates branch length up
#' to each sample's own root, so coalescing and non-coalescing pairs score alike (measured:
#' 29.1 vs 32.1 on a real run). Use [arg_tmrca()] for relatedness in the time domain. This
#' function is kept for comparison with coalescent-based work, and always reports `coalesced`.
#'
#' @param run a run directory or `vivax_run`.
#' @param ids_a,ids_b pedigree ids of the two groups.
#' @return a one-row data.table: `n_a`, `n_b`, `divergence` (branch mode) and
#'   `coalesced` (0-1, the span-weighted fraction of A-B pairs sharing an ancestor).
#' @examples
#' \dontrun{arg_divergence(run, imported_ids, local_ids)}
#' @export
arg_divergence <- function(run, ids_a, ids_b) {
  .arg_load(run, unique(c(ids_a, ids_b)))
  pymain <- reticulate::py
  pymain$`_a` <- as.integer(unique(ids_a)); pymain$`_b` <- as.integer(unique(ids_b))
  reticulate::py_run_string('
p2 = _ped2node(sts)
A = sorted({p2[int(x)] for x in _aslist(_a) if int(x) in p2})
B = sorted({p2[int(x)] for x in _aslist(_b) if int(x) in p2})
if not A or not B:
    _out = dict(n_a=len(A), n_b=len(B), divergence=float("nan"), coalesced=float("nan"))
else:
    _out = dict(n_a=len(A), n_b=len(B),
                divergence=float(sts.divergence([A, B], mode="branch")),
                coalesced=float(_coalesced_frac(sts, A, B)))
')
  data.table::as.data.table(reticulate::py$`_out`)
}

#' Pairwise time to the most recent common ancestor (the ARG-native relatedness measure)
#'
#' For each pair, the span-weighted mean TMRCA **over the part of the genome where a common
#' ancestor exists**, plus `coalesced`, the fraction of the genome where it does. This is the
#' time-domain complement to IBD: IBD asks *how much* of the genome is shared,
#' `arg_tmrca()` asks *how far back* the sharing goes. Unlike [arg_divergence()] it is
#' meaningful on a forest, because non-coalescing spans are excluded rather than scored as 0.
#'
#' Clonal pairs give `coalesced` = 1 with a very recent TMRCA; unrelated genomes give
#' `coalesced` = 0 and `tmrca` = `NA`.
#'
#' @param run a run directory or `vivax_run`.
#' @param ids_a,ids_b equal-length vectors of pedigree ids, compared pairwise.
#' @return a data.table with `a`, `b`, `coalesced` (0-1), `tmrca` and `tmrca_node_time`
#'   (both ticks; `NA` if the pair never coalesces).
#' @section Which TMRCA you get:
#' `tmrca` is a **depth**: ticks from the older of the two sampled genomes back to their
#' common ancestor, so an identical relationship returns an identical number wherever in the
#' run it occurred -- a clonal parent-child pair gives 0 whether it happened on tick 10 or
#' tick 390. `tmrca_node_time` is the raw tskit node time of that ancestor, i.e. ticks before
#' the **end of the run**; that is what `tree.tmrca()` returns natively and it is dominated by
#' when the pair was sampled rather than by how related it is. Use `tmrca` for relatedness.
#' See D22 in `DECISIONS.md`.
#' @examples
#' \dontrun{
#' rel <- load_events(run)[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1]
#' arg_tmrca(run, rel$pedigree_id, rel$parent1_pedigree_id)   # relapse vs its origin
#' }
#' @export
arg_tmrca <- function(run, ids_a, ids_b) {
  stopifnot(length(ids_a) == length(ids_b))
  .arg_load(run, unique(c(ids_a, ids_b)))
  pymain <- reticulate::py
  pymain$`_a` <- as.integer(ids_a); pymain$`_b` <- as.integer(ids_b)
  reticulate::py_run_string('
p2 = _ped2node(sts)
rows = []
for x, y in zip(_aslist(_a), _aslist(_b)):
    a, b = p2.get(int(x)), p2.get(int(y))
    if a is None or b is None:
        rows.append(dict(a=int(x), b=int(y), coalesced=float("nan"),
                         tmrca=float("nan"), tmrca_node_time=float("nan"))); continue
    if a == b:
        rows.append(dict(a=int(x), b=int(y), coalesced=1.0,
                         tmrca=0.0, tmrca_node_time=float(sts.node(a).time))); continue
    # NOTE: no apostrophes in this block -- it is an R single-quoted string.
    # t.tmrca() returns the ancestor NODE TIME: ticks before the end of the run. Reported
    # raw, that makes an identical relationship look different depending on when it was
    # sampled, so subtract the older samples own node time to get a genealogical depth.
    # The ancestor is at least as old as both descendants, so this is never negative.
    base = max(sts.node(a).time, sts.node(b).time)
    tot = 0.0; co = 0.0; acc = 0.0
    for t in sts.trees():
        span = t.interval.right - t.interval.left
        tot += span
        try: tm = t.tmrca(a, b)          # raises if the pair has no common ancestor here
        except ValueError: continue
        co += span; acc += span * tm
    rows.append(dict(a=int(x), b=int(y),
                     coalesced=(co / tot) if tot else float("nan"),
                     tmrca=(acc / co - base) if co > 0 else float("nan"),
                     tmrca_node_time=(acc / co) if co > 0 else float("nan")))
_tm = rows
')
  data.table::rbindlist(lapply(reticulate::py$`_tm`, data.table::as.data.table))
}
