## ibd.R -- ground-truth identity-by-descent from the tree sequence (tskit via reticulate).
## The tskit dependency is a Suggests: these functions error helpfully if it is missing.

.trees_of <- function(run) {
  if (inherits(run, "vivax_run")) return(trees_path(run))
  if (is.character(run) && length(run) == 1) {
    if (dir.exists(run)) return(trees_path(run))
    return(run)  # assume a .trees path
  }
  stop("`run` must be a vivax_run, a run directory, or a .trees path.", call. = FALSE)
}

.require_reticulate <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("IBD needs the 'reticulate' package: install.packages('reticulate').", call. = FALSE)
  if ("py_require" %in% getNamespaceExports("reticulate"))
    try(reticulate::py_require(c("tskit", "pandas")), silent = TRUE)
  ok <- tryCatch({ reticulate::py_run_string("import tskit, pandas"); TRUE }, error = function(e) FALSE)
  if (!ok) stop("IBD needs Python with 'tskit' + 'pandas' (e.g. reticulate::py_require(c('tskit','pandas'))).",
                call. = FALSE)
  invisible(TRUE)
}

#' Ground-truth IBD from a run's tree sequence
#'
#' Computes genome-wide identity-by-descent (shared-segment length / sequence length)
#' among a set of individuals, working on the tskit tree sequence recorded by the model.
#' If `ids` is given the tree is first simplified to just those individuals (memory-safe) and
#' SLiM's vacant null haplosome nodes are always excluded. **All** shared ancestry is counted,
#' however far back it goes -- there is no time cutoff. To impose one, ask for `segments` and
#' filter them yourself; see the note below, and D22 in `DECISIONS.md` for why a cutoff is no
#' longer offered as an argument.
#'
#' @param run a `vivax_run`, a run directory, or a `.trees` file path.
#' @param ids optional integer vector of `pedigree_id`s to restrict to (e.g. sampled cases).
#' @param min_span minimum shared-segment length (bp) to count (default 1000).
#' @param segments if `TRUE`, also return every shared segment with the time of the common
#'   ancestor it came from. Costs roughly 3-7x the default run time (the multiple grows with
#'   the number of samples) and returns one row per segment per pair, so restrict it with
#'   `ids` and `min_span`.
#' @return a list with `ibd` (data.table: `node_a`, `node_b`, `ibd_fraction`), `node_map`
#'   (data.table: `node_id`, `pedigree_id`) and `sequence_length`. With `segments = TRUE`,
#'   also `segments` (data.table: `node_a`, `node_b`, `left`, `right`, `mrca_node`,
#'   `mrca_time`). Pass to [ibd_index()] for lookups.
#' @section Segment times are not depths:
#' `mrca_time` is a tskit node time, i.e. **ticks before the end of the run**, not a
#' genealogical depth. A segment's depth relative to its own pair is
#' `mrca_time - max(t_a, t_b)`, where `t_a`/`t_b` are the two samples' own node times. A
#' recency filter must be written against that difference: filtering on `mrca_time` alone
#' scores a pair by when its lineage last branched, not by how related the two genomes are.
#' @export
ibd_pairs <- function(run, ids = NULL, min_span = 1000, segments = FALSE) {
  .require_reticulate()
  tf <- .trees_of(run)
  if (!file.exists(tf)) stop("tree sequence not found: ", tf, call. = FALSE)
  out <- tempfile(fileext = ".csv"); nmap <- tempfile(fileext = ".csv")
  segf <- tempfile(fileext = ".csv"); lenf <- tempfile(fileext = ".txt")
  want_expr <- if (is.null(ids)) "None" else {
    pf <- tempfile(); writeLines(as.character(unique(ids)), pf)
    paste0("set(int(x) for x in open(r'''", pf, "''').read().split())")
  }
  py <- paste0("
import tskit, pandas as pd, json, gc
ts = tskit.load(r'''", tf, "''')
want = ", want_expr, "
def _md(x):
    md = x.metadata
    if isinstance(md, bytes):
        try: md = json.loads(md.decode()) if md else {}
        except Exception: md = {}
    return md if isinstance(md, dict) else {}
def vac(nd):
    v = _md(nd).get('is_vacant'); return bool(v[0]) if v else False
def pedof(ind):
    md = _md(ind); return md.get('pedigree_id', md.get('pedigreeID'))
keep = []
for n in ts.samples():
    nd = ts.node(n)
    if nd.individual == tskit.NULL or vac(nd): continue
    p = pedof(ts.individual(nd.individual))
    if p is None: continue
    if want is not None and int(p) not in want: continue
    keep.append(int(n))
sts = ts.simplify(samples=keep); del ts; gc.collect()
n2p = {}
for n in sts.samples():
    nd = sts.node(n)
    if nd.individual != tskit.NULL and not vac(nd):
        p = pedof(sts.individual(nd.individual))
        if p is not None: n2p[int(n)] = int(p)
L = sts.sequence_length
open(r'''", lenf, "''', 'w').write(repr(L))
want_seg = ", if (isTRUE(segments)) "True" else "False", "
# store_segments stays OFF unless the caller asked for them: tskit exposes total_span per
# pair without it, and materialising every segment only to sum and discard it costs ~7x
# (38.1s vs 5.1s over 14.8M segments on a 5413-sample tree; 3.2x at 800 ids).
res = sts.ibd_segments(within=list(sts.samples()), min_span=", as.integer(min_span), ",
                       store_pairs=True, store_segments=want_seg)
keys = list(res.keys())
rows = [{'node_a': int(a), 'node_b': int(b), 'ibd_fraction': float(res[(a, b)].total_span) / L}
        for (a, b) in keys]
pd.DataFrame(rows, columns=['node_a','node_b','ibd_fraction']).to_csv(r'''", out, "''', index=False)
pd.DataFrame(list(n2p.items()), columns=['node_id','pedigree_id']).to_csv(r'''", nmap, "''', index=False)
if want_seg:
    nt = sts.nodes_time
    srows = [{'node_a': int(a), 'node_b': int(b), 'left': float(sg.left), 'right': float(sg.right),
              'mrca_node': int(sg.node), 'mrca_time': float(nt[int(sg.node)])}
             for (a, b) in keys for sg in res[(a, b)]]
    pd.DataFrame(srows, columns=['node_a','node_b','left','right','mrca_node','mrca_time']).to_csv(
        r'''", segf, "''', index=False)
")
  reticulate::py_run_string(py)
  ## Pin the column types. fread() on a header-only CSV -- which is what a run with no shared
  ## ancestry above min_span produces -- types every column as logical, so the table silently
  ## changes CLASS with run content and any rbind across runs mis-types or errors. Same
  ## reasoning as the coercion in segments.R.
  ibd <- data.table::fread(out, showProgress = FALSE)
  ibd[, `:=`(node_a = as.integer(node_a), node_b = as.integer(node_b),
             ibd_fraction = as.numeric(ibd_fraction))]
  nm <- data.table::fread(nmap, showProgress = FALSE)
  nm[, `:=`(node_id = as.integer(node_id), pedigree_id = as.integer(pedigree_id))]
  res <- list(ibd = ibd, node_map = nm,
              sequence_length = as.numeric(readLines(lenf, warn = FALSE)))
  if (isTRUE(segments)) {
    sg <- data.table::fread(segf, showProgress = FALSE)
    sg[, `:=`(node_a = as.integer(node_a), node_b = as.integer(node_b),
              left = as.numeric(left), right = as.numeric(right),
              mrca_node = as.integer(mrca_node), mrca_time = as.numeric(mrca_time))]
    res$segments <- sg
  }
  res
}

#' Build a fast pairwise-IBD lookup from [ibd_pairs()]
#'
#' @param ibd_obj the list returned by [ibd_pairs()].
#' @return a list with `nm` (the pedigree/node map) and `getibd(pa, pb)`, a vectorised
#'   function returning genome-wide IBD (0 for pairs sharing nothing).
#' @export
ibd_index <- function(ibd_obj) {
  nm <- data.table::as.data.table(ibd_obj$node_map)[!is.na(pedigree_id)]
  nm <- unique(nm, by = "pedigree_id")
  nmap <- stats::setNames(nm$node_id, as.character(nm$pedigree_id))
  ib <- data.table::as.data.table(ibd_obj$ibd)
  ib[, `:=`(node_a = as.double(node_a), node_b = as.double(node_b),
            ibd_fraction = as.double(ibd_fraction))]
  ib <- data.table::rbindlist(list(ib[, .(a = node_a, b = node_b, ibd = ibd_fraction)],
                                   ib[, .(a = node_b, b = node_a, ibd = ibd_fraction)]))
  data.table::setkey(ib, a, b)
  getibd <- function(pa, pb) {
    na <- nmap[as.character(pa)]; nb <- nmap[as.character(pb)]
    v <- ib[.(na, nb), ibd, on = .(a, b)]
    v <- data.table::fifelse(is.na(v), 0, v)   # 0-default = below the min_span detection floor
    ## tskit's ibd_segments(within=) only returns DISTINCT pairs, so a self-pair would fall
    ## through to the 0-default. A genome is identical to itself by definition, so return 1
    ## (otherwise a pairwise matrix built with the diagonal would be silently wrong).
    ## unname(): the node lookup carries pedigree-id names, which would otherwise ride along
    ## on the returned IBD values and trip up identity comparisons downstream.
    unname(data.table::fifelse(!is.na(na) & !is.na(nb) & na == nb, 1, v))
  }
  list(nm = nm, getibd = getibd)
}
