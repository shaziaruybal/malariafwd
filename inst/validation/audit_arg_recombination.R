## audit_arg_recombination.R -- does the recorded ARG recombine the way the map says?
## The genome is 14 PvP01 chromosomes CONCATENATED into one linkage group, with
## r = 0.5 at each of the 13 junctions (independent assortment) and 7.4e-7/bp within.
## So breakpoints in the tree sequence should show: (a) sharp spikes exactly at the 13
## junctions, (b) a ~uniform background within chromosomes at the map rate, and
## (c) a product of meiosis changing parent at rate r inside a chromosome, so that the
## segments it shares with one of its parents are exponential with mean 1/r, cut short at
## the chromosome ends, and it carries r x genome length (about 17.9) crossovers.
suppressMessages({
  library(malariafwd)
  library(data.table); library(ggplot2); library(patchwork); library(reticulate)
})
py_require(c("tskit", "pandas")); invisible(py_run_string("import tskit"))
FIGS <- Sys.getenv("MALARIAFWD_FIGS", file.path(tempdir(), "malariafwd-figs"))
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
OUT  <- file.path(tempdir(), "mfw_arg")

R_WITHIN <- 7.4e-7
PVP01 <- c(1021664, 956327, 896704, 1012024, 1524814, 1042791, 1652210,
           1761288, 2237066, 1548844, 2131221, 3182763, 2093556, 3153402)
bounds <- cumsum(PVP01)                       # junction after each chromosome (last = end)
junctions <- head(bounds, -1)                 # the 13 internal junctions

unlink(OUT, recursive = TRUE)
run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.14, SIM_DAYS = 500,
                            N_FOUNDER_HOSTS = 20, N_FOUNDER_LINEAGES = 12,
                            LAMBDA_OOCYSTS = 2.0, LAMBDA_SPOROZOITES = 2.0),
               seed = 7, outdir = OUT)

py$trees <- trees_path(run)
py_run_string("
import tskit, numpy as np
ts = tskit.load(trees)
bp = np.array(ts.breakpoints(as_array=True))
bp = bp[(bp > 0) & (bp < ts.sequence_length)]   # drop the two sequence ends
info = dict(num_trees=int(ts.num_trees), num_nodes=int(ts.num_nodes),
            num_edges=int(ts.num_edges), seqlen=float(ts.sequence_length),
            num_samples=int(ts.num_samples))
")
bp <- as.numeric(py$bp); info <- py$info
cat(sprintf("\nARG: %d local trees | %d nodes | %d edges | %d samples | %.2f Mb\n",
            info$num_trees, info$num_nodes, info$num_edges, info$num_samples, info$seqlen / 1e6))
cat(sprintf("distinct breakpoint positions: %d\n", length(bp)))

## ---- (a) how many breakpoints sit exactly on a chromosome junction? ----
on_j <- bp %in% junctions
cat(sprintf("\nbreakpoints exactly at a chromosome junction: %d / %d (%.1f%%)\n",
            sum(on_j), length(bp), 100 * mean(on_j)))
cat(sprintf("  all 13 junctions used? %s\n",
            if (length(unique(bp[on_j])) == 13) "YES" else
              sprintf("only %d of 13", length(unique(bp[on_j])))))
## expected share per product: 13 junctions x r=0.5 vs the within-chromosome map length
exp_j <- 13 * 0.5; exp_w <- sum(PVP01) * R_WITHIN
cat(sprintf("  expected junction share of crossovers = %.1f%% (%.1f junction vs %.1f within per product)\n",
            100 * exp_j / (exp_j + exp_w), exp_j, exp_w))

## ---- (b) is the WITHIN-chromosome background uniform? ----
w <- bp[!on_j]
chr_of <- findInterval(w, c(0, bounds), rightmost.closed = TRUE)
tab <- data.table(chr = chr_of)[, .N, by = chr][order(chr)]
tab[, len_mb := PVP01[chr] / 1e6][, per_mb := N / len_mb]
cat("\nwithin-chromosome breakpoint density (should be ~flat across chromosomes):\n")
print(tab[, .(chr, n = N, len_mb = round(len_mb, 2), per_mb = round(per_mb, 1))])
cat(sprintf("  CV of per-Mb density across chromosomes = %.3f (0 = perfectly uniform)\n",
            sd(tab$per_mb) / mean(tab$per_mb)))
## position WITHIN each chromosome should be uniform -> KS test on the pooled relative position
relpos <- (w - c(0, bounds)[chr_of]) / PVP01[chr_of]
ks <- suppressWarnings(stats::ks.test(relpos, "punif"))
cat(sprintf("  KS test of within-chromosome position vs Uniform(0,1): D = %.3f, p = %.3f%s\n",
            ks$statistic, ks$p.value, if (ks$p.value > 0.01) "  (consistent with uniform)" else "  <-- CHECK"))

## ---- (c) segments a meiotic product shares with its parent ~ Exponential(mean 1/r) ----
## Inside a chromosome a product of four-strand meiosis changes parent at rate r: chiasmata
## fall at rate 2r and each one involves the product's chromatid with probability 1/2. So a
## stretch it shares with one parent ends at rate r, and the stretches are exponential with
## mean 1/r, cut short at the chromosome ends. Take one product from each oocyst whose parents are two
## distinct founders, against its first parent: founders share nothing, so the product is IBD
## with that parent exactly where it carries that parent's genome. (Related parents would
## lengthen the shared stretches, and a stretch two products of one oocyst share ends at a
## crossover in either of them, rate 2r; the product-parent pair isolates the map.)
ev <- as.data.table(load_events(run))
founders <- unique(ev[event_type == "vh_inoculation_seed", lineage_id])
fp <- ev[event_type == "hv_uptake_recombinant" &
         parent1_lineage_id %in% founders & parent2_lineage_id %in% founders,
         .(oocyst_id, parent = parent1_pedigree_id)]
pp <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id %in% fp$oocyst_id,
         .(product = pedigree_id[1]), by = oocyst_id]
pp <- head(merge(pp, fp, by = "oocyst_id"), 100)
segs <- if (nrow(pp)) ibd_segments(run, pp$product, pp$parent) else data.table()
if (nrow(segs)) {
  ## expected lengths: the same process simulated on the map. On each chromosome a product
  ## starts on either parent with probability 1/2 and switches parent at rate r; keep the
  ## stretches it spends on one of them.
  sim_stretches <- function(L, n) {
    k  <- rpois(n, R_WITHIN * L)
    id <- c(rep(seq_len(n), k), seq_len(n))
    x  <- c(runif(sum(k), 0, L), rep(L, n))         # switch points, then the chromosome end
    o  <- order(id, x); id <- id[o]; x <- x[o]
    from <- c(0, head(x, -1)); from[!duplicated(id)] <- 0
    keep <- (rbinom(n, 1, 0.5)[id] + sequence(tabulate(id, n))) %% 2 == 0
    (x - from)[keep]
  }
  set.seed(1)
  expected_len <- unlist(lapply(PVP01, sim_stretches, n = 4000))
  ## each switch of parent inside a chromosome ends or starts a segment shared with this
  ## parent, so the segment endpoints that are not chromosome ends count the crossovers. (On
  ## the full map every product shares part of its genome with each parent; a product that
  ## is not in the tree sequence has no segments and is left out.)
  xo <- segs[, .(n = sum(start > 0) + sum(end < PVP01[chrom])), by = pair]
  xo_se <- sd(xo$n) / sqrt(nrow(xo))
  ks_c <- suppressWarnings(stats::ks.test(segs$length, expected_len))
  cat(sprintf(paste0("\nsegments a meiotic product shares with its parent (%d products, one per oocyst whose two\n",
                     "parents are distinct founders; %d segments):\n"), nrow(xo), nrow(segs)))
  cat(sprintf("  mean %.2f Mb   median %.2f Mb   | expected mean %.2f Mb   median %.2f Mb\n",
              mean(segs$length) / 1e6, median(segs$length) / 1e6,
              mean(expected_len) / 1e6, median(expected_len) / 1e6))
  cat(sprintf("  (expected = exponential with mean 1/r = %.2f Mb, cut short at the chromosome ends; simulated on the map)\n",
              1 / R_WITHIN / 1e6))
  cat(sprintf("  KS test against the expected lengths: D = %.3f, p = %.3f%s\n", ks_c$statistic, ks_c$p.value,
              if (ks_c$p.value > 0.01) "  (consistent with the map)" else "  <-- CHECK"))
  cat(sprintf("  crossovers per product: %.2f +/- %.2f | expected r x genome length = %.2f%s\n",
              mean(xo$n), xo_se, exp_w,
              if (abs(mean(xo$n) - exp_w) <= 3 * xo_se) "  (consistent with the map)" else "  <-- CHECK"))
}

## ---------------- figure ----------------
jd <- data.table(x = junctions / 1e6)
pA <- ggplot(data.table(pos = bp / 1e6), aes(pos)) +
  geom_vline(data = jd, aes(xintercept = x), colour = "#c0392b", linewidth = .35, alpha = .8) +
  geom_histogram(binwidth = 0.25, fill = "#2166ac", colour = NA) +
  scale_x_continuous("position along the concatenated genome (Mb)", expand = c(0.005, 0)) +
  labs(title = "A. Where the ARG recombines",
       subtitle = "red = the 13 chromosome junctions (r = 0.5); blue = breakpoints recorded in the tree sequence",
       y = "breakpoints / 0.25 Mb") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                        plot.title = element_text(face = "bold"))
pB <- ggplot(tab, aes(factor(chr), per_mb)) +
  geom_col(fill = "#4393c3") +
  geom_hline(yintercept = mean(tab$per_mb), linetype = 2, colour = "grey30") +
  labs(title = "B. Within-chromosome breakpoint density is flat",
       subtitle = sprintf("dashed = mean; CV = %.3f. A uniform map should give equal density per Mb on every chromosome.",
                          sd(tab$per_mb) / mean(tab$per_mb)),
       x = "PvP01 chromosome", y = "breakpoints / Mb") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                        plot.title = element_text(face = "bold"))
pl <- list(pA, pB)
if (nrow(segs)) {
  brk <- seq(0, ceiling(max(PVP01) / 1e5) / 10, by = 0.1)
  pC <- ggplot(segs, aes(length / 1e6)) +
    geom_histogram(aes(y = after_stat(density)), breaks = brk, fill = "#d6604d", colour = NA) +
    stat_bin(data = data.table(length = expected_len), aes(y = after_stat(density)), breaks = brk,
             geom = "step", direction = "mid", colour = "grey15", linewidth = .6, linetype = 2) +
    labs(title = "C. Segments a meiotic product shares with its parent follow the map",
         subtitle = sprintf(paste0("One product per oocyst against one of its two parents (unrelated founders). Dashed = expected: exponential\n",
                                   "with mean 1/r = %.2f Mb, cut short at the chromosome ends; the spikes are whole chromosomes from that parent."),
                            1 / R_WITHIN / 1e6),
         x = "IBD segment length (Mb)", y = "density") +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                          plot.title = element_text(face = "bold"))
  pl <- c(pl, list(pC))
}
ggsave(file.path(FIGS, "ARG_recombination.png"), wrap_plots(pl, ncol = 1),
       width = 11, height = 3.2 * length(pl), dpi = 140)
cat("\nwrote", file.path(FIGS, "ARG_recombination.png"), "\n")
