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
#' If `ids` is given the tree is first simplified to just those individuals (memory-safe);
#' SLiM's vacant null haplosome nodes are always excluded, and `max_time` restricts IBD to
#' recent shared ancestry (a transmission-specific, "recent" definition; see DECISIONS D10d).
#'
#' @param run a `vivax_run`, a run directory, or a `.trees` file path.
#' @param ids optional integer vector of `pedigree_id`s to restrict to (e.g. sampled cases).
#' @param min_span minimum shared-segment length (bp) to count (default 1000).
#' @param max_time if set, count only ancestry coalescing within this many time units
#'   (days) -- a recent-transmission IBD definition. `NULL` (default) = any coalescence time.
#' @return a list with `ibd` (data.table: `node_a`, `node_b`, `ibd_fraction`) and
#'   `node_map` (data.table: `node_id`, `pedigree_id`). Pass to [ibd_index()] for lookups.
#' @export
ibd_pairs <- function(run, ids = NULL, min_span = 1000, max_time = NULL) {
  .require_reticulate()
  tf <- .trees_of(run)
  if (!file.exists(tf)) stop("tree sequence not found: ", tf, call. = FALSE)
  out <- tempfile(fileext = ".csv"); nmap <- tempfile(fileext = ".csv")
  want_expr <- if (is.null(ids)) "None" else {
    pf <- tempfile(); writeLines(as.character(unique(ids)), pf)
    paste0("set(int(x) for x in open(r'''", pf, "''').read().split())")
  }
  mt <- if (is.null(max_time)) "None" else format(max_time, scientific = FALSE, trim = TRUE)
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
kw = dict(min_span=", as.integer(min_span), ", store_pairs=True, store_segments=True)
mt = ", mt, "
if mt is not None: kw['max_time'] = mt
rows = [{'node_a': int(a), 'node_b': int(b), 'ibd_fraction': float(sum(s.right-s.left for s in segs)/L)}
        for (a, b), segs in sts.ibd_segments(within=list(sts.samples()), **kw).items()]
pd.DataFrame(rows, columns=['node_a','node_b','ibd_fraction']).to_csv(r'''", out, "''', index=False)
pd.DataFrame(list(n2p.items()), columns=['node_id','pedigree_id']).to_csv(r'''", nmap, "''', index=False)
")
  reticulate::py_run_string(py)
  list(ibd = data.table::fread(out, showProgress = FALSE),
       node_map = data.table::fread(nmap, showProgress = FALSE))
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
