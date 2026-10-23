## audit_arg_recombination.R -- does the recorded ARG recombine the way the map says?
## The genome is 14 PvP01 chromosomes CONCATENATED into one linkage group, with
## r = 0.5 at each of the 13 junctions (independent assortment) and 7.4e-7/bp within.
## So breakpoints in the tree sequence should show: (a) sharp spikes exactly at the 13
## junctions, (b) a ~uniform background within chromosomes at the map rate, and
## (c) IBD segment lengths ~ exponential with mean 1/r.
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
                            N_FOUNDERS = 20, N_FOUNDER_HOSTS = 20, N_FOUNDER_LINEAGES = 12,
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
## expected share: 13 junctions x r=0.5 vs total map length
exp_j <- 13 * 0.5; exp_w <- sum(PVP01) * R_WITHIN
cat(sprintf("  expected junction share of crossovers = %.1f%% (%.1f junction vs %.1f within per meiosis)\n",
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

## ---- (c) IBD segment lengths from a single meiosis ~ Exponential(mean 1/r) ----
ev <- as.data.table(load_events(run))
oo <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id != -1,
         .(ped = pedigree_id, oocyst_id, lineage_id)]
sib <- oo[, if (.N > 1 && uniqueN(lineage_id) > 1) {
  i <- which(lineage_id != lineage_id[1])[1]; .(a = ped[1], b = ped[i]) }, by = oocyst_id][!is.na(b)]
segs <- data.table()
if (nrow(sib)) {
  py$pa <- as.integer(head(sib$a, 40)); py$pb <- as.integer(head(sib$b, 40))
  py_run_string("
import tskit, pandas as pd, json
ts = tskit.load(trees)
def _md(x):
    md = x.metadata
    if isinstance(md, bytes):
        try: md = json.loads(md.decode()) if md else {}
        except Exception: md = {}
    return md if isinstance(md, dict) else {}
def vac(nd):
    v = _md(nd).get('is_vacant'); return bool(v[0]) if v else False
p2n = {}
for n in ts.samples():
    nd = ts.node(n)
    if nd.individual == tskit.NULL or vac(nd): continue
    m = _md(ts.individual(nd.individual)); p = m.get('pedigree_id', m.get('pedigreeID'))
    if p is not None: p2n.setdefault(int(p), int(n))
rows=[]
for x,y in zip(pa,pb):
    a,b = p2n.get(int(x)), p2n.get(int(y))
    if a is None or b is None or a==b: continue
    for (_i,_j), ss in ts.ibd_segments(within=[a,b], store_pairs=True, store_segments=True, min_span=1000).items():
        for s in ss: rows.append({'len': float(s.right-s.left)})
seg_df = pd.DataFrame(rows, columns=['len'])
")
  segs <- as.data.table(py$seg_df)
}
if (nrow(segs)) {
  cat(sprintf("\nIBD segment lengths between same-oocyst sibs (n = %d segments):\n", nrow(segs)))
  cat(sprintf("  mean %.2f Mb   median %.2f Mb   | expected mean 1/r = %.2f Mb\n",
              mean(segs$len) / 1e6, median(segs$len) / 1e6, 1 / R_WITHIN / 1e6))
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
  pC <- ggplot(segs, aes(len / 1e6)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d6604d", colour = NA) +
    stat_function(fun = function(x) dexp(x, rate = 1 / (1 / R_WITHIN / 1e6)),
                  colour = "grey15", linewidth = .7, linetype = 2) +
    labs(title = "C. Shared-segment lengths follow the map",
         subtitle = sprintf("same-oocyst sibs (one meiosis). Dashed = Exponential with mean 1/r = %.2f Mb.",
                            1 / R_WITHIN / 1e6),
         x = "IBD segment length (Mb)", y = "density") +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                          plot.title = element_text(face = "bold"))
  pl <- c(pl, list(pC))
}
ggsave(file.path(FIGS, "ARG_recombination.png"), wrap_plots(pl, ncol = 1),
       width = 11, height = 3.2 * length(pl), dpi = 140)
cat("\nwrote", file.path(FIGS, "ARG_recombination.png"), "\n")
