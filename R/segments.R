## segments.R -- pairwise IBD SEGMENTS (where along the genome two genomes share ancestry, D74)
## and the two diagnostic plots built on them.
##
## This is the pairwise view: for a given PAIR, which stretches are identical by descent. The
## per-genome view, which root genome each genome descends from at each locus of a panel, is the
## `roots` matrix of simulate_genotypes() (genotypes.R).

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
#' @inheritParams load_events
#' @return a data.table of `chrom` (PvP01 number), `start`, `end`, `length` in concatenated bp.
#' @examples
#' chrom_map(example_run())
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
#' @inheritParams load_events
#' @param ids_a,ids_b equal-length vectors of pedigree ids, compared position by position.
#' @param min_span the shortest segment returned, in bp, after the cut at the junctions (a segment
#'   of exactly `min_span` is kept): below it a shared stretch is treated as noise, mirroring the
#'   detection floor of IBD inference.
#' @return a data.table, by pair and position: `pair` (the position in `ids_a` and `ids_b`, from 1),
#'   `a` and `b` (as given), `chrom` (the PvP01 chromosome), `start` and `end` (bp within the
#'   chromosome, `end` exclusive), `gstart` and `gend` (bp along the concatenated genome), `length`
#'   (bp), `mrca_id` (the pedigree id of the common ancestor), `mrca_time` and `depth` (its age;
#'   see [ibd_pairs()]). A pair sharing nothing has no rows, and nor does a pair with a genome the
#'   tree sequence does not hold, which warns; the table's `"missing"` attribute lists those ids. A
#'   genome with itself is shared end to end: one row per chromosome, its own ancestor at depth 0. A
#'   pair given twice has its rows twice.
#' @examplesIf tskit_available()
#' run <- example_run()
#' rel <- load_events(run)[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][1]
#' ibd_segments(run, rel$pedigree_id, rel$parent1_pedigree_id)   # all 14 chromosomes shared
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
  ## the requested ids the tree sequence lacks, whose pairs have no rows, for plot_ibd_pair()
  data.table::setattr(out, "missing", x$missing)
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
#' `group` sets kinds of pair side by side, such as relapses against their origin, siblings and
#' unrelated pairs: the pairs of a group are drawn in one panel under its name, one row per pair
#' from the top in the order given, and each group's shared blocks take a colour of their own. In
#' the chromosome layout each panel then has a column per chromosome, as wide as the chromosome is
#' long, and each row's label gives the pair's IBD.
#'
#' @inheritParams load_events
#' @param ids_a,ids_b equal-length vectors of pedigree ids, compared pairwise.
#' @param labels optional labels, one per pair (defaults to `a vs b`). A pair with a genome the
#'   tree sequence lacks is drawn empty, its IBD shown as `NA`.
#' @param group optional, one value per pair: pairs with the same value are drawn in one panel,
#'   the panels in the order the values first appear, or in the order of the levels of a factor.
#' @param min_span detection floor in bp, passed to [ibd_segments()].
#' @param layout `"chromosome"` (one bar per chromosome) or `"concatenated"`.
#' @param fill colour for the shared (IBD) blocks. With `group`, one colour for every group, or one
#'   for each group, in group order or named by group; left out, each group gets its own.
#' @return a ggplot object.
#' @examplesIf tskit_available() && requireNamespace("ggplot2", quietly = TRUE)
#' run <- example_run()
#' ev <- load_events(run)
#' rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1][1]
#' founders <- ev[event_type == "vh_inoculation_seed"]$pedigree_id[1:2]
#' plot_ibd_pair(run, c(rel$pedigree_id, founders[1]), c(rel$parent1_pedigree_id, founders[2]),
#'               labels = c("a relapse vs its batch", "two founders"))
#'
#' # three pairs of each kind, a panel and a colour each
#' tp <- ibd_truth_pairs(run)
#' tp <- tp[cls %in% c("relapse vs its origin", "same-oocyst sibs"), head(.SD, 3), by = cls]
#' plot_ibd_pair(run, tp$a, tp$b, group = tp$cls, layout = "concatenated")
#' @inheritSection tskit_available Python
#' @export
plot_ibd_pair <- function(run, ids_a, ids_b, labels = NULL, group = NULL, min_span = 1000,
                          layout = c("chromosome", "concatenated"), fill = "#2166ac") {
  .need_ggplot()
  layout <- match.arg(layout)
  if (!is.null(labels) && (length(labels) != length(ids_a) || anyNA(labels)))
    stop("`labels` must give one label per pair, with no NA: there are ", length(ids_a),
         " pairs, and it has ", length(labels), " values", if (anyNA(labels)) ", some NA", ".",
         call. = FALSE)
  if (!is.null(group)) {
    if (length(group) != length(ids_a) || anyNA(group))
      stop("`group` must give one group per pair, with no NA: there are ", length(ids_a),
           " pairs, and it has ", length(group), " values", if (anyNA(group)) ", some NA", ".",
           call. = FALSE)
    group <- if (is.factor(group)) droplevels(group) else factor(group, unique(group))
    cols <- .ibd_group_fill(fill, levels(group), missing(fill))
  }
  .require_python("plot_ibd_pair()", c("tskit", "numpy"))
  seg <- ibd_segments(run, ids_a, ids_b, min_span = min_span)
  ## a pair with a genome the tree sequence lacks has no IBD to show, rather than none shared
  absent <- ids_a %in% attr(seg, "missing") | ids_b %in% attr(seg, "missing")
  cm <- chrom_map(run); L <- max(cm$end)
  if (is.null(labels)) labels <- sprintf("%s vs %s", ids_a, ids_b)
  if (!is.null(group)) return(.plot_ibd_groups(seg, cm, L, labels, group, cols, layout, absent))
  lv <- unique(labels)
  key <- data.table::data.table(pair = seq_along(ids_a), lab = factor(labels, lv),
                                x0 = 0, x1 = L / 1e6)
  ibd <- merge(key, seg[, .(ibd = sum(length) / L), by = pair], by = "pair", all.x = TRUE)
  ibd[is.na(ibd), ibd := 0]
  ibd[pair %in% which(absent), ibd := NA_real_]
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
        ggplot2::scale_x_continuous("Position within chromosome (Mb)", expand = c(0.01, 0)) +
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
    ggplot2::scale_x_continuous("Position along the concatenated genome (Mb)", expand = c(0, 0)) +
    ggplot2::labs(y = NULL, caption = "dividers = PvP01 chromosome boundaries") +
    th + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1))
}

## plot_ibd_pair()'s colours for its groups, named by group: its own palette when `fill` is left
## out, else one colour for all, one per group in group order, or one named for each group.
.GROUP_FILL <- c("#2166ac", "#c8971f", "#2a9d8f", "#b4472e", "#8e6bbf", "#2e7d5b")
.ibd_group_fill <- function(fill, lv, default) {
  n <- length(lv)
  if (default) {
    fill <- if (n <= length(.GROUP_FILL)) .GROUP_FILL[seq_len(n)]
            else grDevices::colorRampPalette(.GROUP_FILL)(n)
  } else if (!is.null(names(fill))) {
    gone <- setdiff(lv, names(fill))
    if (length(gone))
      stop("`fill` names no colour for group", if (length(gone) > 1L) "s", " ",
           .and_list(paste0("\"", gone, "\"")), ".", call. = FALSE)
    fill <- fill[lv]
  } else if (length(fill) == 1L) {
    fill <- rep(fill, n)
  } else if (length(fill) != n) {
    stop("`fill` must be one colour, or one for each of the ", n, " groups (in group order, or ",
         "named by group); it has ", length(fill), ".", call. = FALSE)
  }
  stats::setNames(unname(fill), lv)
}

## plot_ibd_pair() with `group`: a panel per group, a row per pair (the first at the top; keyed by
## pair, so a label given twice keeps both rows), the shared blocks in the group's colour.
.plot_ibd_groups <- function(seg, cm, L, labels, group, cols, layout, absent) {
  n <- length(group)
  ibd <- numeric(n); v <- seg[, .(ibd = sum(length) / L), by = pair]; ibd[v$pair] <- v$ibd
  ibd[absent] <- NA_real_
  key <- data.table::data.table(pair = seq_len(n), grp = group, row = factor(seq_len(n), rev(seq_len(n))),
                                ibd = ibd)
  segl <- merge(seg, key[, .(pair, grp, row)], by = "pair")
  rowlab <- if (layout == "chromosome") sprintf("%s (IBD %.2f)", labels, ibd) else labels
  common <- list(
    ggplot2::scale_y_discrete(NULL, labels = function(b) rowlab[as.integer(b)]),
    ggplot2::scale_fill_manual(values = cols, limits = names(cols), guide = "none"),
    theme_malariafwd(grid = FALSE),
    ggplot2::theme(strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "bold")))

  if (layout == "chromosome") {
    ## a column per chromosome, as wide as it is long
    bg <- key[, cm[, .(chrom, len = length / 1e6)], by = .(pair, grp, row)]
    return(
      ggplot2::ggplot() +
        ggplot2::geom_tile(data = bg, ggplot2::aes(x = len / 2, width = len, y = row), height = .7,
                           fill = "grey92", colour = "grey15", linewidth = .25) +
        ggplot2::geom_tile(data = segl, ggplot2::aes(x = (start + end) / 2e6, width = (end - start) / 1e6,
                           y = row, fill = grp), height = .7) +
        ggplot2::facet_grid(grp ~ chrom, scales = "free", space = "free", switch = "y") +
        ggplot2::scale_x_continuous("PvP01 chromosome", breaks = NULL, expand = c(0.02, 0)) +
        ggplot2::labs(caption = "columns = PvP01 chromosomes (1-14), width by length; colour = shared IBD") +
        common + ggplot2::theme(panel.spacing.x = ggplot2::unit(1.5, "pt"))
    )
  }
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = key, ggplot2::aes(x = L / 2e6, width = L / 1e6, y = row), height = .7,
                       fill = "grey92", colour = "grey15", linewidth = .3) +
    ggplot2::geom_tile(data = segl, ggplot2::aes(x = (gstart + gend) / 2e6, width = (gend - gstart) / 1e6,
                       y = row, fill = grp), height = .7) +
    ggplot2::geom_vline(data = cm[-nrow(cm)], ggplot2::aes(xintercept = end / 1e6),
                        colour = "grey15", linewidth = .3) +
    ggplot2::geom_text(data = key, ggplot2::aes(x = L / 1e6, y = row, label = sprintf("IBD = %.2f", ibd)),
                       hjust = -0.1, size = 3.2, colour = "grey25") +
    ggplot2::facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_x_continuous("Position along the concatenated genome (Mb)",
                                expand = ggplot2::expansion(mult = c(0, .12))) +
    ggplot2::labs(caption = "dividers = PvP01 chromosome boundaries; colour = shared IBD") +
    common
}

#' Pairs of genomes with a known true relationship, and their IBD
#'
#' Assembles pairs whose relationship is known from the event log -- a genome with itself,
#' clonal founder copies, a relapse against the hypnozoite batch it activated from, the same
#' lineage in two hosts, siblings co-transmitted from one oocyst, and distinct founder
#' lineages -- and attaches the ground-truth IBD of each. Feed it to
#' [plot_relatedness_ladder()] to check the genealogy against the textbook values.
#'
#' @inheritParams load_events
#' @param max_per_class cap on pairs sampled per relationship class.
#' @return a data.table of `cls` (a factor of the six classes), `a` and `b` (integer pedigree ids)
#'   and `ibd` (double), one row per pair. A class with no pair in the run has no rows, and a run
#'   with none, such as one that never started, gives a table with no rows.
#' @examplesIf tskit_available()
#' tp <- ibd_truth_pairs(example_run())
#' tp[, .(pairs = .N, mean_ibd = mean(ibd)), by = cls]
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
  ## pair, and the filter on b always finds its column.
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
  x <- suppressWarnings(.ibd_pairs(run, pairs = out[, .(a, b)], none_ok = TRUE))
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
#' @param x a run, or the table [ibd_truth_pairs()] returns.
#' @param style `"points"` plots every pair against the theoretical rungs; `"bars"` shows the
#'   composition of each class across relatedness bands, which reads better when classes are
#'   point masses (relapses are always clonal, distinct founders always unrelated). Both label
#'   the number of pairs per class, so a class resting on two observations is visible.
#' @return a ggplot object.
#' @examplesIf tskit_available() && requireNamespace("ggplot2", quietly = TRUE)
#' plot_relatedness_ladder(ibd_truth_pairs(example_run()))
#' @export
plot_relatedness_ladder <- function(x, style = c("points", "bars")) {
  .need_ggplot()
  style <- match.arg(style)
  pairs <- .as_result(x, ibd_truth_pairs, c("cls", "ibd"))
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
        ## show.legend = TRUE: every band keeps its colour in the legend, also one no pair falls in
        ## (ggplot2 leaves the key of an unused level blank)
        ggplot2::geom_col(width = .68, colour = "grey40", linewidth = .25, show.legend = TRUE) +
        ggplot2::geom_text(data = n_by, ggplot2::aes(cls, 1.04, label = paste0("n = ", n)),
                           inherit.aes = FALSE, size = 3, colour = "grey30") +
        ggplot2::scale_fill_manual(values = pal2, name = "IBD band", drop = FALSE) +
        ggplot2::scale_y_continuous("Fraction of pairs",
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
