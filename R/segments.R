## segments.R -- pairwise IBD SEGMENTS (where along the genome two genomes share ancestry, D74)
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
#' parameters the run itself recorded (`GENOME_MODE`, `CHROMOSOMES`), and stops if they do not
#' add up to the genome length the run recorded.
#' @param run a run directory or `vivax_run`.
#' @return a data.table of `chrom` (PvP01 number), `start`, `end`, `length` in concatenated bp.
#' @examples
#' \dontrun{chrom_map(run)}
#' @export
chrom_map <- function(run) {
  r <- .read_params_tsv(file.path(.run_dir(run), "params.tsv"))
  p <- r$params
  idx <- .check_chromosomes(p$CHROMOSOMES)
  len <- .PVP01_FULL[idx]
  if (p$GENOME_MODE == "core") len <- round(len * 0.90)
  gl <- r$record$GENOME_LENGTH
  if (!is.null(gl) && sum(len) != gl)
    stop("CHROMOSOMES \"", p$CHROMOSOMES, "\" in ", p$GENOME_MODE, " mode is ", .slim_num(sum(len)),
         " bp, but the run recorded GENOME_LENGTH ", .slim_num(gl), ".", call. = FALSE)
  data.table::data.table(chrom = idx, length = len,
                         start = cumsum(c(0, utils::head(len, -1))), end = cumsum(len))
}

#' Pairwise IBD segments along the genome
#'
#' The stretches of genome over which two parasites are identical by descent -- the detail
#' behind the single IBD number from [ibd_pairs()]. A segment is a maximal stretch over which the
#' two genomes descend from one common ancestor genome, cut at each chromosome junction, the same
#' rows as `ibd_pairs(segments = TRUE)` gives for the pair; each is in both concatenated and
#' within-chromosome coordinates.
#'
#' @param run a run directory or `vivax_run`.
#' @param ids_a,ids_b equal-length vectors of pedigree ids, compared position by position.
#' @param min_span the shortest segment returned, in bp, after the cut at the junctions (a segment
#'   of exactly `min_span` is kept): below it a shared stretch is treated as noise, mirroring the
#'   detection floor of IBD inference.
#' @return a data.table, by pair and position: `pair` (the position in `ids_a` and `ids_b`, from 1),
#'   `a` and `b` (as given), `chrom` (the PvP01 chromosome), `start` and `end` (bp within the
#'   chromosome, `end` exclusive), `gstart` and `gend` (bp along the concatenated genome), `length`
#'   (bp), `mrca_id` (the pedigree id of the common ancestor), `mrca_time` and `depth` (its age;
#'   see [ibd_pairs()]). A pair sharing nothing has no rows, and nor does a pair with a genome the
#'   tree sequence does not hold, which warns. A genome with itself is shared end to end: one row
#'   per chromosome, its own ancestor at depth 0. A pair given twice has its rows twice.
#' @examples
#' \dontrun{
#' rel <- load_events(run)[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][1]
#' ibd_segments(run, rel$pedigree_id, rel$parent1_pedigree_id)
#' }
#' @seealso [plot_ibd_pair()], [ibd_pairs()]
#' @inheritSection tskit_available Python
#' @export
ibd_segments <- function(run, ids_a, ids_b, min_span = 1000) {
  ids_a <- .ibd_ids(ids_a, "ids_a"); ids_b <- .ibd_ids(ids_b, "ids_b")
  if (length(ids_a) != length(ids_b))
    stop("`ids_a` and `ids_b` must have the same length, one pair per position; they have ",
         length(ids_a), " and ", length(ids_b), ".", call. = FALSE)
  .ibd_check_min_span(min_span)
  .require_python("ibd_segments()", c("tskit", "numpy"))
  x <- .ibd_compute(run, pa = ids_a, pb = ids_b, min_span = min_span, segments = TRUE,
                    absent = "their pairs have no rows.")
  key <- data.table::data.table(pair = seq_along(ids_a), a = ids_a, b = ids_b,
                                lo = pmin(ids_a, ids_b), hi = pmax(ids_a, ids_b))
  sg <- data.table::copy(x$segments)
  data.table::setnames(sg, c("a", "b"), c("lo", "hi"))
  two <- sg[key[a != b], on = c("lo", "hi"), nomatch = NULL, allow.cartesian = TRUE]
  ## a genome with itself: every chromosome, its own ancestor
  self <- key[a == b & a %in% x$samples$pedigree_id]
  cm <- if (nrow(self)) .ibd_layout(run, x$trees, x$sequence_length) else NULL
  one <- if (nrow(self))
    self[, .(chrom = as.integer(cm$chrom), start = 0, end = cm$length, gstart = cm$start, gend = cm$end,
             length = cm$length, mrca_id = a,
             mrca_time = x$samples$node_time[match(a, x$samples$pedigree_id)], depth = 0),
         by = .(pair, a, b)]
  cols <- c("pair", "a", "b", "chrom", "start", "end", "gstart", "gend", "length", "mrca_id", "mrca_time", "depth")
  out <- data.table::rbindlist(list(two[, cols, with = FALSE], if (!is.null(one)) one[, cols, with = FALSE]),
                               use.names = TRUE)
  data.table::setorder(out, pair, gstart)
  out[]
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
#' @inheritSection tskit_available Python
#' @export
plot_ibd_pair <- function(run, ids_a, ids_b, labels = NULL, min_span = 1000,
                          layout = c("chromosome", "concatenated"), fill = "#2166ac") {
  .need_ggplot()
  layout <- match.arg(layout)
  .require_python("plot_ibd_pair()", c("tskit", "numpy"))
  seg <- ibd_segments(run, ids_a, ids_b, min_span = min_span)
  cm <- chrom_map(run); L <- max(cm$end)
  if (is.null(labels)) labels <- sprintf("%s vs %s", ids_a, ids_b)
  lv <- unique(labels)
  key <- data.table::data.table(pair = seq_along(ids_a), lab = factor(labels, lv),
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
#' @return a data.table of `cls` (a factor of the six classes), `a` and `b` (integer pedigree ids)
#'   and `ibd` (double), one row per pair. A class with no pair in the run has no rows, and a run
#'   with none, such as one that never started, gives a table with no rows.
#' @examples
#' \dontrun{plot_relatedness_ladder(ibd_truth_pairs(run))}
#' @inheritSection tskit_available Python
#' @export
ibd_truth_pairs <- function(run, max_per_class = 60) {
  .require_python("ibd_truth_pairs()", c("tskit", "numpy"))
  ev <- data.table::as.data.table(load_events(run))
  cs <- classify_cases(ev)
  seeds <- ev[event_type == "vh_inoculation_seed"]
  P <- list(); add <- function(cls, a, b) if (length(a))
    P[[length(P) + 1]] <<- data.table::data.table(cls = cls, a = as.integer(a), b = as.integer(b))

  add(.LADDER[1], utils::head(cs$pedigree_id, 20), utils::head(cs$pedigree_id, 20))
  cl <- seeds[, .N, by = lineage_id][N > 1, lineage_id]
  for (l in utils::head(cl, max_per_class)) { pp <- seeds[lineage_id == l, pedigree_id]
    if (length(pp) >= 2) add(.LADDER[2], pp[1], pp[2]) }
  ## != -1, not > 0: -1 is the no-parent sentinel but pedigree 0 is a real genome (D23).
  rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1]
  add(.LADDER[3], utils::head(rel$pedigree_id, max_per_class),
                  utils::head(rel$parent1_pedigree_id, max_per_class))
  ## acquisitions only: `!= "relapse"` would admit recrudescences, which inherit their holder's
  ## lineage_id and are not transmissions.
  ct <- cs[origin %in% c("imported", "local reinfection"), .(ped = pedigree_id, human_id, lineage_id)]
  ## Each pair's second member is the first genome of another host (or lineage), NA when there is
  ## none. The expression is unconditional, so a and b are there, typed, also when no group has a
  ## pair: one that returned nothing for such a group left them out, and filtering on b failed.
  ctp <- ct[, .(a = ped[1], b = ped[which(human_id != human_id[1])[1]]), by = lineage_id][!is.na(b)]
  add(.LADDER[4], utils::head(ctp$a, max_per_class), utils::head(ctp$b, max_per_class))
  oo <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id != -1,
           .(ped = pedigree_id, oocyst_id, lineage_id)]
  sib <- oo[, .(a = ped[1], b = ped[which(lineage_id != lineage_id[1])[1]]), by = oocyst_id][!is.na(b)]
  add(.LADDER[5], utils::head(sib$a, max_per_class), utils::head(sib$b, max_per_class))
  fl <- unique(seeds, by = "lineage_id")$pedigree_id
  if (length(fl) >= 4) { cb <- utils::combn(length(fl), 2)
    add(.LADDER[6], fl[cb[1, ]], fl[cb[2, ]]) }

  ## a run with no pair, such as one that never started, gives the empty table at once: ids of NULL
  ## would ask ibd_pairs() for every genome
  if (!length(P))
    return(data.table::data.table(cls = factor(character(), .LADDER), a = integer(), b = integer(),
                                  ibd = numeric()))
  out <- data.table::rbindlist(P)
  x <- suppressWarnings(ibd_pairs(run, pairs = out[, .(a, b)]))
  out[, ibd := ibd_lookup(x, a, b)]
  out[!is.na(ibd)][, cls := factor(cls, .LADDER)][]
}

.LADDER <- c("self (same genome)", "clonal copies", "relapse vs its origin",
             "same lineage, 2 hosts", "same-oocyst sibs", "distinct founders")

#' Plot the relatedness ladder
#'
#' IBD by known true relationship, against the textbook rungs (0 unrelated, 1/4 half-sib,
#' 1/3 the mean over a meiotic tetrad, 1/2 full-sib, 1 clonal). The model-validation view:
#' clonal relationships should sit at 1, distinct founders at 0, and same-oocyst siblings (any
#' two products of one meiosis) in two groups: meiosis II sisters around 0.6 and the other pairs
#' around 0.2, a third on average for unrelated parents and more when the parents are related. On
#' the full map none shares nothing.
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
    ## shows where same-oocyst sibs fall: for unrelated parents, meiosis II sisters in (0.4, 0.8]
    ## and the other pairs below 0.4.
    brk <- c(-Inf, 1e-9, 0.25, 0.40, 0.80, 1 - 1e-9, Inf)
    lab <- c("0 (unrelated)", "(0, 1/4]", "(1/4, 0.4]", "(0.4, 0.8]", "(0.8, 1)", "1 (clonal)")
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
    ## the mean, which for same-oocyst sibs is the tetrad mean (their two groups put the median low)
    ggplot2::stat_summary(fun = mean, geom = "crossbar", width = .5,
                          linewidth = .5, colour = "grey15") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_y_continuous("IBD (fraction of genome shared identical by descent; all ancestry counted)",
                                breaks = rungs$y, labels = rungs$lab, limits = c(-0.02, 1.10)) +
    ggplot2::labs(x = NULL) +
    theme_malariafwd(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}
