## arg.R -- ARG-native (tree) summaries of the recorded genealogy.
##
## The recorded genealogy is a FOREST, not a single tree: each founder lineage and each import is
## an independent root, so most pairs of genomes never coalesce within the simulated period. On a
## forest, tskit's branch statistics are not a relatedness measure: for a pair with no common
## ancestor they accumulate branch length up to each genome's own root, so pairs that coalesce and
## pairs that never do score alike. arg_tmrca() is the relatedness measure here: the time to the
## common ancestor where one exists, with `coalesced`, the fraction of the genome over which it
## does. arg_summary() shows the forest itself.

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
'

.arg_load <- function(run, ids = NULL, what) {
  .require_python(what, c("tskit", "numpy"))
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
#' @inheritSection tskit_available Python
#' @export
arg_summary <- function(run, ids = NULL) {
  .arg_load(run, ids, "arg_summary()")
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
#' @inheritSection tskit_available Python
#' @export
arg_tmrca <- function(run, ids_a, ids_b) {
  stopifnot(length(ids_a) == length(ids_b))
  .arg_load(run, unique(c(ids_a, ids_b)), "arg_tmrca()")
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
