## audit_arg_junctions.R -- two things the breakpoint histogram cannot show:
##  (1) INDEPENDENT ASSORTMENT at the 13 chromosome junctions (r = 0.5). Counting distinct
##      breakpoint positions cannot test this: a junction is ONE position, so it contributes
##      one breakpoint however many meioses crossed there. The real test is whether IBD state
##      DECORRELATES across a junction, while staying correlated over the same bp distance
##      within a chromosome.
##  (2) WHY same-oocyst sibs are sometimes IBD = 0 -- the tetrad structure.
suppressMessages({
  library(malariafwd)
  library(data.table); library(ggplot2); library(patchwork); library(reticulate)
})
py_require(c("tskit", "pandas")); invisible(py_run_string("import tskit"))
FIGS <- Sys.getenv("MALARIAFWD_FIGS", file.path(tempdir(), "malariafwd-figs"))
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
OUT  <- file.path(tempdir(), "mfw_junc")

PVP01 <- c(1021664, 956327, 896704, 1012024, 1524814, 1042791, 1652210,
           1761288, 2237066, 1548844, 2131221, 3182763, 2093556, 3153402)
bounds <- cumsum(PVP01); junctions <- head(bounds, -1)
D <- 100000   # flanking distance (bp) used on both tests

unlink(OUT, recursive = TRUE)
run <- run_sim(vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.14, SIM_DAYS = 500,
                            N_FOUNDERS = 20, N_FOUNDER_HOSTS = 20, N_FOUNDER_LINEAGES = 12,
                            LAMBDA_OOCYSTS = 2.0, LAMBDA_SPOROZOITES = 2.0),
               seed = 7, outdir = OUT)
ev <- as.data.table(load_events(run))

## ---------- identify tetrad products among co-transmitted infections ----------
## seedTetrad() makes 4 products per outcrossing meiosis: two RECIPROCAL recombinants
## (fresh lineage ids) and two PARENTAL CLONES (they inherit the parents' lineage ids).
rec <- ev[event_type == "hv_uptake_recombinant",
          .(oocyst_id, p1 = parent1_lineage_id, p2 = parent2_lineage_id)]
inoc <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id != -1,
           .(ped = pedigree_id, oocyst_id, lineage_id = as.character(lineage_id))]
tet <- merge(inoc, rec, by = "oocyst_id")
tet[, kind := fifelse(lineage_id == as.character(p1) | lineage_id == as.character(p2),
                      "parental clone", "recombinant")]
prs <- tet[, if (.N > 1) as.data.table(t(combn(.N, 2)))[, .(
  a = ped[V1], b = ped[V2], ka = kind[V1], kb = kind[V2],
  la = lineage_id[V1], lb = lineage_id[V2])], by = oocyst_id]
prs <- prs[la != lb]
prs[, pair_type := fcase(
  ka == "recombinant"    & kb == "recombinant",    "recombinant x recombinant\n(reciprocal halves)",
  ka == "parental clone" & kb == "parental clone", "clone x clone\n(the two parents)",
  default = "recombinant x parental clone")]

## ---------- IBD for those pairs ----------
idx <- ibd_index(ibd_pairs(run, ids = unique(c(prs$a, prs$b))))
inN <- function(x) x %in% idx$nm$pedigree_id
prs <- prs[inN(a) & inN(b)][, ibd := idx$getibd(a, b)]
cat("\n== same-oocyst sib pairs, split by which tetrad products they are ==\n")
print(prs[, .(n = .N, median_ibd = round(median(ibd), 3), mean_ibd = round(mean(ibd), 3),
              pct_zero = round(100 * mean(ibd < 0.01), 1)), by = pair_type][order(-n)])

## ---------- independent assortment across junctions ----------
## For each pair, is it IBD at (x - D) and at (x + D)? Junction positions vs control
## positions the same distance apart inside a chromosome.
ctrl <- sort(unique(pmin(pmax(round(c(bounds - PVP01) + PVP01 / 2), D + 1),
                         bounds - D - 1)))            # chromosome midpoints
## Use pairs related by EXACTLY ONE meiosis (a recombinant vs one of its parental clones).
## For those, an r = 0.5 junction must give P(same IBD state either side) = 0.5 exactly.
one <- prs[pair_type == "recombinant x parental clone"][1:min(.N, 120)]
py$pa2 <- as.integer(one$a); py$pb2 <- as.integer(one$b); py$trees <- trees_path(run)
py$jx <- as.numeric(junctions); py$cx <- as.numeric(ctrl); py$dd <- as.numeric(D)
py_run_string("
import tskit, pandas as pd, json, itertools
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
pairs = [(p2n[int(x)], p2n[int(y)]) for x, y in zip(pa2, pb2)
         if int(x) in p2n and int(y) in p2n and p2n[int(x)] != p2n[int(y)]]
rows = []
L = float(ts.sequence_length)
for a, b in pairs:
    segs = []
    for (_i,_j), ss in ts.ibd_segments(within=[a,b], store_pairs=True, store_segments=True, min_span=1000).items():
        segs += [(float(s.left), float(s.right)) for s in ss]
    if not segs: continue
    # CONDITION ON MOSAIC PAIRS. Clonal (IBD~1) and unrelated (IBD~0) pairs are in the SAME
    # state on both sides of every position by construction, which would swamp the linkage
    # signal we are testing. Only partially-IBD pairs carry information about recombination.
    frac = sum(r-l for l, r in segs) / L
    if not (0.05 < frac < 0.95): continue
    def is_ibd(x):
        return any(l <= x < r for l, r in segs)
    for kind, xs in (('junction', jx), ('within-chromosome control', cx)):
        for x in xs:
            rows.append({'kind': kind, 'x': float(x),
                         'left': bool(is_ibd(x-dd)), 'right': bool(is_ibd(x+dd))})
conc_df = pd.DataFrame(rows, columns=['kind','x','left','right'])
")
cc <- as.data.table(py$conc_df)
## concordance beyond chance: phi coefficient between IBD-left and IBD-right
res <- cc[, {
  ## as.numeric: these counts multiply to > .Machine$integer.max and would overflow
  a <- as.numeric(sum(left & right));  b <- as.numeric(sum(left & !right))
  c2 <- as.numeric(sum(!left & right)); d <- as.numeric(sum(!left & !right))
  den <- sqrt((a + b) * (c2 + d) * (a + c2) * (b + d))
  .(n = .N, pct_same = round(100 * mean(left == right), 1),
    phi = round(if (den > 0) (a * d - b * c2) / den else NA_real_, 3))
}, by = kind]
cat(sprintf("\n== IBD state %d kb either side of a position (n pairs x positions) ==\n", D / 1000))
print(res)
cat("  expected: junction r=0.5 -> phi ~ 0 (independent); within-chromosome r=7.4e-7*2e5=0.148 -> phi high\n")

## ---------------- figure ----------------
pal2 <- c("recombinant x parental clone" = "#4393c3",
          "recombinant x recombinant\n(reciprocal halves)" = "#c0392b",
          "clone x clone\n(the two parents)" = "#e08a1e")
pA <- ggplot(prs, aes(ibd, fill = pair_type)) +
  geom_histogram(binwidth = .05, colour = NA) +
  facet_wrap(~ pair_type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = pal2, guide = "none") +
  labs(title = "Why same-oocyst sibs are sometimes IBD = 0",
       subtitle = "One outcrossing meiosis makes 4 products: 2 reciprocal recombinants + 2 parental clones.\nOnly recombinant-vs-clone pairs share ancestry; the other two combinations share none.",
       x = "IBD", y = "pairs") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"),
                                        strip.text = element_text(size = 9))
pB <- ggplot(res, aes(kind, phi, fill = kind)) +
  geom_col(width = .55) + geom_hline(yintercept = 0, colour = "grey40") +
  geom_text(aes(label = sprintf("phi = %.3f\n%.0f%% same state", phi, pct_same)), vjust = -0.3, size = 3.4) +
  scale_fill_manual(values = c("junction" = "#c0392b", "within-chromosome control" = "#2166ac"),
                    guide = "none") +
  coord_cartesian(ylim = c(-0.1, 1.15)) +
  labs(title = sprintf("Independent assortment: IBD state %d kb either side", D / 1000),
       subtitle = "Pairs related by exactly ONE meiosis (a recombinant vs its parental clone).\nr = 0.5 at a junction should give 50% same state (phi ~ 0); 100 kb inside a chromosome should stay correlated.",
       x = NULL, y = "phi (correlation of IBD state across the position)") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(FIGS, "ARG_junctions_and_tetrads.png"), pA / pB + plot_layout(heights = c(1.3, 1)),
       width = 10, height = 9, dpi = 140)
cat("\nwrote", file.path(FIGS, "ARG_junctions_and_tetrads.png"), "\n")
