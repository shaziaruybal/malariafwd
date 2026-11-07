## explain_ibd.R -- what IBD is, and whether the model produces the textbook values.
## Fig 1: the actual shared genome SEGMENTS for one example pair of each relationship class.
## Fig 2: the relatedness LADDER -- IBD by known-truth class vs the theoretical expectation.
suppressMessages({
  library(malariafwd)
  library(data.table); library(ggplot2); library(patchwork); library(reticulate)
})
if (!malariafwd::tskit_available())
  stop("these checks need Python with tskit (>= 1.0), pandas and numpy; see malariafwd::check_setup()",
       call. = FALSE)
invisible(py_run_string("import tskit"))
FIGS <- Sys.getenv("MALARIAFWD_FIGS", file.path(tempdir(), "malariafwd-figs"))
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
OUT  <- file.path(tempdir(), "mfw_ibd_explain")

## ---- a run with enough outcrossing to populate the sib rungs ----
unlink(OUT, recursive = TRUE)
p <- vivax_params(N_HOSTS = 120, N_MOSQ = 240, BITING_RATE = 0.14, SIM_DAYS = 500,
                  N_FOUNDER_HOSTS = 20, N_FOUNDER_LINEAGES = 12,
                  LAMBDA_OOCYSTS = 2.0, LAMBDA_SPOROZOITES = 2.0)
run <- run_sim(p, seed = 7, outdir = OUT)
ev  <- as.data.table(load_events(run)); cs <- classify_cases(ev)
seeds <- ev[event_type == "vh_inoculation_seed"]

## ---------------- known-truth pair sets ----------------
pairs <- list()
add <- function(cls, a, b) if (length(a)) pairs[[length(pairs) + 1]] <<-
  data.table(cls = cls, a = as.integer(a), b = as.integer(b))

# 1. the same genome compared with itself
add("self (same genome)", cs$pedigree_id[1:20], cs$pedigree_id[1:20])

# 2. clonal founder copies: same lineage seeded into different hosts
cl <- seeds[, .N, by = lineage_id][N > 1, lineage_id]
for (l in head(cl, 8)) { pp <- seeds[lineage_id == l, pedigree_id]
  if (length(pp) >= 2) add("clonal copies", pp[1], pp[2]) }

# 3. a relapse vs the hypnozoite batch it activated from
rel <- ev[event_type == "hypnozoite_activation" & parent1_pedigree_id != -1,
          .(a = pedigree_id, b = parent1_pedigree_id)]
add("relapse vs its origin", head(rel$a, 60), head(rel$b, 60))

# 4. clonal transmission: same lineage, different hosts, via a bite
ct <- cs[origin != "relapse", .(ped = pedigree_id, human_id, lineage_id)]
ctp <- ct[, if (uniqueN(human_id) > 1) .(a = ped[1], b = ped[which(human_id != human_id[1])[1]]),
          by = lineage_id][!is.na(b)]
add("same lineage, 2 hosts", head(ctp$a, 60), head(ctp$b, 60))

# 5. same-oocyst sibs: two infections co-transmitted from ONE meiosis
oo <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id != -1,
         .(ped = pedigree_id, oocyst_id, lineage_id)]
sib <- oo[, if (.N > 1 && uniqueN(lineage_id) > 1) {
  i <- which(lineage_id != lineage_id[1])[1]; .(a = ped[1], b = ped[i]) }, by = oocyst_id][!is.na(b)]
add("same-oocyst sibs", head(sib$a, 60), head(sib$b, 60))

# 6. distinct founder lineages: independent genome roots
fl <- unique(seeds, by = "lineage_id")$pedigree_id
if (length(fl) >= 4) { cb <- combn(length(fl), 2)
  add("distinct founders", fl[cb[1, ]], fl[cb[2, ]]) }

P <- rbindlist(pairs)
idx <- ibd_index(ibd_pairs(run, ids = unique(c(P$a, P$b))))
inN <- function(x) x %in% idx$nm$pedigree_id
P <- P[inN(a) & inN(b)]
P[, ibd := idx$getibd(a, b)]

LV <- c("self (same genome)", "clonal copies", "relapse vs its origin", "same lineage, 2 hosts",
        "same-oocyst sibs", "distinct founders")
P[, cls := factor(cls, LV)]
cat("\n== IBD by known-truth relationship ==\n")
print(P[, .(n = .N, median = round(median(ibd), 3), mean = round(mean(ibd), 3)), by = cls][order(cls)])

## ---------------- Fig 2: the relatedness ladder ----------------
rungs <- data.table(y = c(0, 1/4, 1/3, 1/2, 1),
                    lab = c("0 unrelated", "1/4 half-sib", "1/3 mean over a tetrad",
                            "1/2 full-sib / parent", "1 clonal"))
pal <- c("self (same genome)" = "#2166ac", "clonal copies" = "#4393c3",
         "relapse vs its origin" = "#c8971f", "same lineage, 2 hosts" = "#7fbf7b",
         "same-oocyst sibs" = "#d6604d", "distinct founders" = "#999999")
p2 <- ggplot(P, aes(cls, ibd, colour = cls)) +
  geom_hline(data = rungs, aes(yintercept = y), linetype = 2, colour = "grey70", linewidth = .3) +
  geom_jitter(width = .18, height = 0, alpha = .35, size = 1.1) +
  stat_summary(fun = mean, geom = "crossbar", width = .5, linewidth = .5, colour = "grey15") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_y_continuous("IBD (fraction of genome shared from a recent common ancestor)",
                     breaks = rungs$y, labels = rungs$lab, limits = c(-0.02, 1.02)) +
  labs(title = "The relatedness ladder: IBD by known relationship",
       subtitle = "each dot is one pair of parasite genomes; bar = mean. Dashed lines are the textbook expectations.",
       x = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(FIGS, "IBD_ladder.png"), p2, width = 10, height = 6, dpi = 140)

## ---------------- Fig 1: the shared segments themselves ----------------
## one exemplar pair per class, with the actual IBD segments along the genome: the pair nearest
## the class median, which for same-oocyst sibs is a typical pair that is not a meiosis II sister
ex <- P[, .SD[which.min(abs(ibd - median(ibd)))], by = cls][order(cls)]
## NB ibd_pairs() simplifies the tree, which RENUMBERS nodes - its node_ids are not valid in a
## freshly loaded tree. So resolve pedigree_id -> node from the original tree's metadata here.
py$peds_a <- as.integer(ex$a); py$peds_b <- as.integer(ex$b)
py$trees  <- trees_path(run)
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
ped2node = {}
for n in ts.samples():
    nd = ts.node(n)
    if nd.individual == tskit.NULL or vac(nd): continue
    md = _md(ts.individual(nd.individual))
    p = md.get('pedigree_id', md.get('pedigreeID'))
    if p is not None: ped2node.setdefault(int(p), int(n))
rows=[]
for i,(pa,pb) in enumerate(zip(peds_a, peds_b)):
    a, b = ped2node.get(int(pa)), ped2node.get(int(pb))
    if a is None or b is None: continue
    if a == b:
        rows.append({'pair':i,'left':0.0,'right':float(ts.sequence_length)}); continue
    res = ts.ibd_segments(within=[a,b], store_pairs=True, store_segments=True, min_span=1000)
    for (_x,_y), segs in res.items():
        for s in segs: rows.append({'pair':i,'left':float(s.left),'right':float(s.right)})
segs_df = pd.DataFrame(rows, columns=['pair','left','right'])
seqlen = float(ts.sequence_length)
")
segs <- as.data.table(py$segs_df); L <- py$seqlen
ex[, pair := .I - 1L]
segs <- merge(segs, ex[, .(pair, cls, ibd)], by = "pair")
segs[, cls := factor(cls, LV)]
lab <- ex[, .(cls = factor(cls, LV), txt = sprintf("IBD = %.2f", ibd))]

bg <- data.table(cls = factor(LV, LV))   # one full-genome backdrop per class (a 0-IBD pair has no segments)
p1 <- ggplot(segs) +
  geom_rect(data = bg, aes(xmin = 0, xmax = L / 1e6, ymin = 0, ymax = 1), fill = "grey92") +
  geom_rect(aes(xmin = left / 1e6, xmax = right / 1e6, ymin = 0, ymax = 1, fill = cls)) +
  geom_text(data = lab, aes(x = L / 1e6, y = .5, label = txt), hjust = 1.05, size = 3.4, colour = "grey20") +
  facet_wrap(~ cls, ncol = 1, strip.position = "left") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous("position along the 14-chromosome genome (Mb)", expand = c(0, 0)) +
  labs(title = "What IBD is: the stretches of genome two parasites inherited from a shared ancestor",
       subtitle = "grey = the genome; coloured = segments shared identical-by-descent. IBD is just the coloured fraction.",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 10),
        plot.title = element_text(face = "bold", size = 12))
ggsave(file.path(FIGS, "IBD_segments.png"), p1, width = 11, height = 6, dpi = 140)
cat("\nwrote", file.path(FIGS, "IBD_segments.png"), "and IBD_ladder.png\n")
