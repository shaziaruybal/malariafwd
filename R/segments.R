## segments.R -- pairwise IBD SEGMENTS (where along the genome two genomes share ancestry)
## and the two diagnostic plots built on them.
##
## NB distinct from the planned founder-mosaic `extract_segments()` (which paints each sample
## by the FOUNDER each block descends from). This is the pairwise view: for a given PAIR,
## which stretches are identical by descent.

## PvP01 chromosome lengths (bp), chr 1..14; GenBank LT635612.1..LT635625.1. Mirrors the
## PVP01_FULL constant in the SLiM model; "core" mode is 90% of each (see the model's GENOME
## block). The model concatenates these into one linkage group with r = 0.5 at each junction.
.PVP01_FULL <- c(1021664, 956327, 896704, 1012024, 1524814, 1042791, 1652210,
                 1761288, 2237066, 1548844, 2131221, 3182763, 2093556, 3153402)

#' Chromosome layout of a run's concatenated genome
#'
#' The model simulates the 14 PvP01 chromosomes as ONE concatenated linkage group. This
#' recovers the chromosome boundaries in those concatenated coordinates, from the genome
#' parameters the run itself recorded (`GENOME_MODE`, `CHROMOSOMES`).
#' @param run a run directory or `vivax_run`.
#' @return a data.table of `chrom` (PvP01 number), `start`, `end`, `length` in concatenated bp.
#' @examples
#' \dontrun{chrom_map(run)}
#' @export
chrom_map <- function(run) {
  p <- load_params(run)
  spec <- as.character(p[["CHROMOSOMES"]] %||% "1-14")
  idx <- unlist(lapply(strsplit(spec, ",")[[1]], function(tok) {
    q <- strsplit(trimws(tok), "-")[[1]]
    if (length(q) == 2) seq(as.integer(q[1]), as.integer(q[2])) else as.integer(q)
  }))
  len <- .PVP01_FULL[idx]
  if (identical(as.character(p[["GENOME_MODE"]] %||% "full"), "core")) len <- round(len * 0.90)
  data.table::data.table(chrom = idx, length = len,
                         start = cumsum(c(0, utils::head(len, -1))), end = cumsum(len))
}

#' Pairwise IBD segments along the genome
#'
#' The stretches of genome over which two parasites are identical by descent -- the detail
#' behind the single IBD number from [ibd_pairs()]. Segments are split at chromosome
#' boundaries and reported in both concatenated and within-chromosome coordinates.
#'
#' @param run a run directory or `vivax_run`.
#' @param ids_a,ids_b equal-length vectors of pedigree ids, compared pairwise.
#' @param min_span minimum segment length in bp (the detection floor; below this a shared
#'   stretch is treated as noise, mirroring real IBD inference).
#' @return a data.table of `pair`, `a`, `b`, `chrom`, `start`, `end` (within-chromosome bp),
#'   `gstart`, `gend` (concatenated bp) and `length`. A pair sharing nothing yields no rows.
#' @examples
#' \dontrun{
#' rel <- load_events(run)[event_type == "hypnozoite_activation" & parent1_pedigree_id > 0][1]
#' ibd_segments(run, rel$pedigree_id, rel$parent1_pedigree_id)
#' }
#' @export
ibd_segments <- function(run, ids_a, ids_b, min_span = 1000) {
  stopifnot(length(ids_a) == length(ids_b))
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("ibd_segments() needs the 'reticulate' package.", call. = FALSE)
  if ("py_require" %in% getNamespaceExports("reticulate"))
    try(reticulate::py_require(c("tskit", "pandas")), silent = TRUE)
  ok <- tryCatch({ reticulate::py_run_string("import tskit, pandas"); TRUE }, error = function(e) FALSE)
  if (!ok) stop("ibd_segments() needs Python with 'tskit' + 'pandas'.", call. = FALSE)
  pymain <- reticulate::py
  pymain$`_trees` <- trees_path(run)
  pymain$`_a` <- as.integer(ids_a); pymain$`_b` <- as.integer(ids_b)
  pymain$`_ms` <- as.integer(min_span)
  reticulate::py_run_string('
import tskit, pandas as pd, json
def _md(x):
    m = x.metadata
    if isinstance(m, bytes):
        try: m = json.loads(m.decode()) if m else {}
        except Exception: m = {}
    return m if isinstance(m, dict) else {}
def _vac(nd):
    v = _md(nd).get("is_vacant"); return bool(v[0]) if v else False
def _lst(x):
    if x is None: return []
    try: return [int(v) for v in x]
    except TypeError: return [int(x)]
ts = tskit.load(_trees)
# resolve pedigree_id -> node on the ORIGINAL tree (simplify would renumber nodes)
p2n = {}
for n in ts.samples():
    nd = ts.node(n)
    if nd.individual == tskit.NULL or _vac(nd): continue
    m = _md(ts.individual(nd.individual))
    p = m.get("pedigree_id", m.get("pedigreeID"))
    if p is not None: p2n.setdefault(int(p), int(n))
rows = []
for i, (x, y) in enumerate(zip(_lst(_a), _lst(_b))):
    na, nb = p2n.get(int(x)), p2n.get(int(y))
    if na is None or nb is None: continue
    if na == nb:                       # a genome is identical to itself end to end
        rows.append({"pair": i, "a": int(x), "b": int(y),
                     "gstart": 0.0, "gend": float(ts.sequence_length)}); continue
    for (_i, _j), segs in ts.ibd_segments(within=[na, nb], store_pairs=True,
                                          store_segments=True, min_span=_ms).items():
        for s in segs:
            rows.append({"pair": i, "a": int(x), "b": int(y),
                         "gstart": float(s.left), "gend": float(s.right)})
_seg = pd.DataFrame(rows, columns=["pair","a","b","gstart","gend"])
')
  seg <- data.table::as.data.table(reticulate::py$`_seg`)
  ## A pair sharing NOTHING yields an empty frame, and pandas types empty columns as object,
  ## which reticulate hands back as character. Every downstream join is on integer `pair`, so
  ## the types are pinned here -- otherwise an unrelated pair (a legitimate and common case)
  ## fails to merge in plot_ibd_pair() rather than drawing an empty karyotype.
  seg[, `:=`(pair = as.integer(pair), a = as.integer(a), b = as.integer(b),
             gstart = as.numeric(gstart), gend = as.numeric(gend))]
  ## Return the columns in the SAME ORDER as the non-empty path below. Appending them with
  ## `:=` left gstart/gend ahead of chrom/start/end, so a pair that shares nothing produced a
  ## different column ORDER from one that shares something -- and rbindlist's default
  ## use.names = "check" then binds by position, emitting only a message while pouring
  ## `gstart` values into `chrom` (1021664 instead of 1). Plausible-looking silent corruption
  ## in any loop over pairs or runs.
  if (!nrow(seg))
    return(seg[, .(pair, a, b, chrom = integer(), start = numeric(), end = numeric(),
                   gstart, gend, length = numeric())][])
  cm <- chrom_map(run)
  ## Split each segment at chromosome boundaries (a segment can span a junction when the same
  ## parent happens to continue across it). Done as a cross join + overlap filter, PER SEGMENT:
  ## comparing the 14-row chromosome table against a whole group's segment vector would
  ## silently recycle and mis-assign chromosomes for any pair with more than one segment.
  seg[, segid := .I]
  j <- seg[, cm[, .(chrom, cs = start, ce = end)], by = .(segid, pair, a, b, gstart, gend)]
  j <- j[ce > gstart & cs < gend]
  j[, `:=`(gs = pmax(gstart, cs), ge = pmin(gend, ce))]
  j[, .(pair, a, b, chrom, start = gs - cs, end = ge - cs,
        gstart = gs, gend = ge, length = ge - gs)][order(pair, gstart)]
}

#' Plot the IBD "karyotype" of one or more pairs
#'
#' Paints the genome for each pair and colours the stretches shared identical by descent, so
#' you can see *why* a pair's IBD value is what it is: a clonal pair (or a relapse against its
#' origin) is shared end to end, same-oocyst siblings show a recombination mosaic, and
#' unrelated genomes share nothing.
#'
#' `layout = "chromosome"` (the default) draws the 14 PvP01 chromosomes as separate bars with
#' each pair in its own facet -- the conventional karyotype view, and the one that shows
#' whether shared blocks respect chromosome boundaries. `layout = "concatenated"` draws the
#' single linkage group the model actually simulates, with chromosome dividers marked.
#'
#' @param run a run directory or `vivax_run`.
#' @param ids_a,ids_b equal-length vectors of pedigree ids, compared pairwise.
#' @param labels optional labels, one per pair (defaults to `a vs b`).
#' @param min_span detection floor in bp, passed to [ibd_segments()].
#' @param layout `"chromosome"` (one bar per chromosome, pairs as facets) or `"concatenated"`.
#' @param fill colour for the shared (IBD) blocks.
#' @return a ggplot object.
#' @examples
#' \dontrun{plot_ibd_pair(run, c(id1, id3), c(id2, id4), labels = c("relapse", "unrelated"))}
#' @export
plot_ibd_pair <- function(run, ids_a, ids_b, labels = NULL, min_span = 1000,
                          layout = c("chromosome", "concatenated"), fill = "#2166ac") {
  .need_ggplot()
  layout <- match.arg(layout)
  seg <- ibd_segments(run, ids_a, ids_b, min_span = min_span)
  cm <- chrom_map(run); L <- max(cm$end)
  if (is.null(labels)) labels <- sprintf("%s vs %s", ids_a, ids_b)
  lv <- unique(labels)
  key <- data.table::data.table(pair = seq_along(ids_a) - 1L, lab = factor(labels, lv),
                                x0 = 0, x1 = L / 1e6)
  ibd <- merge(key, seg[, .(ibd = sum(length) / L), by = pair], by = "pair", all.x = TRUE)
  ibd[is.na(ibd), ibd := 0]
  segl <- merge(seg, key[, .(pair, lab)], by = "pair")
  th <- theme_malariafwd(grid = FALSE)

  if (layout == "chromosome") {
    ## one bar per chromosome (karyotype), pairs as facet columns
    bg <- key[, cm[, .(chrom, len_mb = length / 1e6)], by = .(pair, lab)]
    bg[, ypos := factor(chrom, levels = rev(cm$chrom))]
    segl[, ypos := factor(chrom, levels = rev(cm$chrom))]
    return(
      ggplot2::ggplot() +
        ggplot2::geom_rect(data = bg, ggplot2::aes(xmin = 0, xmax = len_mb,
                           ymin = as.integer(ypos) - .38, ymax = as.integer(ypos) + .38),
                           fill = "grey92", colour = "grey15", linewidth = .3) +
        ggplot2::geom_rect(data = segl, ggplot2::aes(xmin = start / 1e6, xmax = end / 1e6,
                           ymin = as.integer(ypos) - .38, ymax = as.integer(ypos) + .38),
                           fill = fill, colour = "grey15", linewidth = .2) +
        ggplot2::geom_text(data = ibd, ggplot2::aes(x = Inf, y = 0.1,
                           label = sprintf("IBD = %.2f", ibd)), hjust = 1.05, size = 3.2,
                           colour = "grey25") +
        ggplot2::facet_wrap(~ lab, nrow = 1) +
        ggplot2::scale_y_continuous(NULL, breaks = seq_len(nrow(cm)),
                                    labels = rev(cm$chrom), expand = ggplot2::expansion(add = .8)) +
        ggplot2::scale_x_continuous("position within chromosome (Mb)", expand = c(0.01, 0)) +
        ggplot2::labs(caption = "rows = PvP01 chromosomes (1-14); grey = chromosome, colour = shared IBD") +
        th
    )
  }
  ## single concatenated linkage group, with chromosome dividers on top
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = key, ggplot2::aes(xmin = x0, xmax = x1, ymin = 0, ymax = 1),
                       fill = "grey92", colour = "grey15", linewidth = .3) +
    ggplot2::geom_rect(data = segl, ggplot2::aes(xmin = gstart / 1e6, xmax = gend / 1e6,
                       ymin = 0, ymax = 1), fill = fill) +
    ggplot2::geom_vline(data = cm[-nrow(cm)], ggplot2::aes(xintercept = end / 1e6),
                        colour = "grey15", linewidth = .3) +
    ggplot2::geom_text(data = ibd, ggplot2::aes(x = L / 1e6, y = .5,
                       label = sprintf("IBD = %.2f", ibd)), hjust = 1.05, size = 3.3,
                       colour = "grey25") +
    ggplot2::facet_wrap(~ lab, ncol = 1, strip.position = "left") +
    ggplot2::scale_x_continuous("position along the concatenated genome (Mb)", expand = c(0, 0)) +
    ggplot2::labs(y = NULL, caption = "dividers = PvP01 chromosome boundaries") +
    th + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1))
}

#' Pairs of genomes with a known true relationship, and their IBD
#'
#' Assembles pairs whose relationship is known from the event log -- a genome with itself,
#' clonal founder copies, a relapse against the hypnozoite batch it activated from, the same
#' lineage in two hosts, siblings co-transmitted from one oocyst, and distinct founder
#' lineages -- and attaches the ground-truth IBD of each. Feed it to
#' [plot_relatedness_ladder()] to check the genealogy against the textbook values.
#'
#' @param run a run directory or `vivax_run`.
#' @param max_per_class cap on pairs sampled per relationship class.
#' @return a data.table of `cls`, `a`, `b`, `ibd`.
#' @examples
#' \dontrun{plot_relatedness_ladder(ibd_truth_pairs(run))}
#' @export
ibd_truth_pairs <- function(run, max_per_class = 60) {
  ev <- data.table::as.data.table(load_events(run))
  cs <- classify_cases(ev)
  seeds <- ev[event_type == "vh_inoculation_seed"]
  P <- list(); add <- function(cls, a, b) if (length(a))
    P[[length(P) + 1]] <<- data.table::data.table(cls = cls, a = as.integer(a), b = as.integer(b))

  add(.LADDER[1], utils::head(cs$pedigree_id, 20), utils::head(cs$pedigree_id, 20))
  cl <- seeds[, .N, by = lineage_id][N > 1, lineage_id]
  for (l in utils::head(cl, max_per_class)) { pp <- seeds[lineage_id == l, pedigree_id]
    if (length(pp) >= 2) add(.LADDER[2], pp[1], pp[2]) }
  rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id > 0]
  add(.LADDER[3], utils::head(rel$pedigree_id, max_per_class),
                  utils::head(rel$parent1_pedigree_id, max_per_class))
  ct <- cs[origin != "relapse", .(ped = pedigree_id, human_id, lineage_id)]
  ctp <- ct[, if (data.table::uniqueN(human_id) > 1)
    .(a = ped[1], b = ped[which(human_id != human_id[1])[1]]), by = lineage_id][!is.na(b)]
  add(.LADDER[4], utils::head(ctp$a, max_per_class), utils::head(ctp$b, max_per_class))
  oo <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id > 0,
           .(ped = pedigree_id, oocyst_id, lineage_id)]
  sib <- oo[, if (.N > 1 && data.table::uniqueN(lineage_id) > 1) {
    i <- which(lineage_id != lineage_id[1])[1]; .(a = ped[1], b = ped[i]) },
    by = oocyst_id]
  ## A grouped expression that matches NO group returns a 0-row table carrying only the `by`
  ## column, so filtering on `b` errors rather than returning nothing. This happens whenever a
  ## run has no outcrossed oocyst that produced two distinct lineages - e.g. any run where
  ## treatment suppresses transmission hard enough.
  sib <- if (nrow(sib) && "b" %in% names(sib)) sib[!is.na(b)] else sib[0]
  add(.LADDER[5], utils::head(sib$a, max_per_class), utils::head(sib$b, max_per_class))
  fl <- unique(seeds, by = "lineage_id")$pedigree_id
  if (length(fl) >= 4) { cb <- utils::combn(length(fl), 2)
    add(.LADDER[6], fl[cb[1, ]], fl[cb[2, ]]) }

  out <- data.table::rbindlist(P)
  idx <- ibd_index(ibd_pairs(run, ids = unique(c(out$a, out$b))))
  out <- out[a %in% idx$nm$pedigree_id & b %in% idx$nm$pedigree_id]
  out[, ibd := idx$getibd(a, b)][, cls := factor(cls, .LADDER)][]
}

.LADDER <- c("self (same genome)", "clonal copies", "relapse vs its origin",
             "same lineage, 2 hosts", "same-oocyst sibs", "distinct founders")

#' Plot the relatedness ladder
#'
#' IBD by known true relationship, against the textbook rungs (0 unrelated, 1/4 half-sib,
#' 1/3 the mean over a meiotic tetrad, 1/2 full-sib, 1 clonal). The model-validation view:
#' clonal relationships should sit at 1, distinct founders at 0, and same-oocyst siblings
#' should split between ~1/2 and 0 (reciprocal recombinants share nothing).
#'
#' @param pairs output of [ibd_truth_pairs()].
#' @param style `"points"` plots every pair against the theoretical rungs; `"bars"` shows the
#'   composition of each class across relatedness bands, which reads better when classes are
#'   point masses (relapses are always clonal, distinct founders always unrelated). Both label
#'   the number of pairs per class, so a class resting on two observations is visible.
#' @return a ggplot object.
#' @examples
#' \dontrun{plot_relatedness_ladder(ibd_truth_pairs(run), style = "bars")}
#' @export
plot_relatedness_ladder <- function(pairs, style = c("points", "bars")) {
  .need_ggplot()
  style <- match.arg(style)
  pairs <- data.table::as.data.table(pairs)
  n_by <- pairs[, .(n = .N), by = cls]

  if (style == "bars") {
    ## Most relationship classes are POINT MASSES (relapses are always IBD 1, distinct founders
    ## always 0), so a box or a crossbar collapses to a line and implies a spread that is not
    ## there. Binning into relatedness bands and showing composition is honest about that, and
    ## shows the one genuinely interesting split -- same-oocyst sibs, which divide between high
    ## sharing and the reciprocal products that share nothing.
    brk <- c(-Inf, 1e-9, 0.25, 0.40, 0.60, 1 - 1e-9, Inf)
    lab <- c("0 (unrelated)", "(0, 1/4]", "(1/4, 0.4]", "(0.4, 0.6]", "(0.6, 1)", "1 (clonal)")
    d <- data.table::copy(pairs)[, band := cut(ibd, brk, labels = lab, right = TRUE)]
    d <- d[, .(n = .N), by = .(cls, band)][, frac := n / sum(n), by = cls]
    pal2 <- stats::setNames(
      c("#f0f0f0", "#d9e7f0", "#9ecae1", "#fdd49e", "#fc8d59", "#b2182b"), lab)
    return(
      ggplot2::ggplot(d, ggplot2::aes(cls, frac, fill = band)) +
        ggplot2::geom_col(width = .68, colour = "grey40", linewidth = .25) +
        ggplot2::geom_text(data = n_by, ggplot2::aes(cls, 1.04, label = paste0("n = ", n)),
                           inherit.aes = FALSE, size = 3, colour = "grey30") +
        ggplot2::scale_fill_manual(values = pal2, name = "IBD band", drop = FALSE) +
        ggplot2::scale_y_continuous("fraction of pairs",
                                    labels = function(x) paste0(round(100 * x), "%"),
                                    limits = c(0, 1.08), breaks = c(0, .25, .5, .75, 1)) +
        ggplot2::labs(x = NULL) +
        theme_malariafwd(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
    )
  }

  rungs <- data.table::data.table(
    y = c(0, 1/4, 1/3, 1/2, 1),
    lab = c("0 unrelated", "1/4 half-sib", "1/3 tetrad mean", "1/2 full-sib / parent", "1 clonal"))
  pal <- stats::setNames(c("#2166ac", "#4393c3", "#c8971f", "#7fbf7b", "#d6604d", "#999999"), .LADDER)
  ggplot2::ggplot(pairs, ggplot2::aes(cls, ibd, colour = cls)) +
    ggplot2::geom_text(data = n_by, ggplot2::aes(cls, 1.06, label = paste0("n = ", n)),
                       inherit.aes = FALSE, size = 3, colour = "grey30") +
    ggplot2::geom_hline(data = rungs, ggplot2::aes(yintercept = y), linetype = 2,
                        colour = "grey70", linewidth = .3) +
    ggplot2::geom_jitter(width = .18, height = 0, alpha = .35, size = 1.1) +
    ggplot2::stat_summary(fun = stats::median, geom = "crossbar", width = .5,
                          linewidth = .5, colour = "grey15") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous("IBD (fraction of genome shared from a recent common ancestor)",
                                breaks = rungs$y, labels = rungs$lab, limits = c(-0.02, 1.10)) +
    ggplot2::labs(x = NULL) +
    theme_malariafwd(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}
