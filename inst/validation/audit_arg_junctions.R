## audit_arg_junctions.R -- two things the breakpoint histogram cannot show:
##  (1) INDEPENDENT ASSORTMENT at the 13 chromosome junctions (r = 0.5). Counting distinct
##      breakpoint positions cannot test this: a junction is ONE position, so it contributes
##      one breakpoint however many meioses crossed there. The real test follows one product
##      of a meiosis against one of its parents: whether it carries that parent's genome must
##      DECORRELATE across a junction, while staying correlated over the same bp distance
##      within a chromosome, and the product must switch parent at half the junctions.
##  (2) how much the products of one meiosis share: a third of the genome on average when the
##      parents are unrelated, more when they are related, and on the full 14-chromosome map
##      at the default rate, never nothing.
suppressMessages({
  library(malariafwd)
  library(data.table); library(ggplot2); library(patchwork); library(reticulate)
})
py_require(c("tskit", "pandas")); invisible(py_run_string("import tskit"))
FIGS <- Sys.getenv("MALARIAFWD_FIGS", file.path(tempdir(), "malariafwd-figs"))
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)
OUT  <- file.path(tempdir(), "mfw_junc")

R_WITHIN <- 7.4e-7
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
## seedTetrad() makes 4 products per outcrossing meiosis by four-strand meiosis (D34): each is a
## mosaic of the two parents and a fresh lineage, and two of the four carry each parent's genome
## at every locus.
rec <- ev[event_type == "hv_uptake_recombinant",
          .(oocyst_id, p1 = as.character(parent1_lineage_id), p2 = as.character(parent2_lineage_id),
            par1 = parent1_pedigree_id, par2 = parent2_pedigree_id)]
founders <- unique(as.character(ev[event_type == "vh_inoculation_seed", lineage_id]))
## parents that are two distinct founders are unrelated; any other two parents may share part
## of their genome
rec[, unrelated_parents := p1 %in% founders & p2 %in% founders]
inoc <- ev[event_type == "vh_inoculation" & is_dormant == 0 & oocyst_id != -1,
           .(ped = pedigree_id, oocyst_id, lineage_id = as.character(lineage_id))]
tet <- merge(inoc, rec[, .(oocyst_id)], by = "oocyst_id")
prs <- tet[, if (.N > 1) as.data.table(t(combn(.N, 2)))[, .(
  a = ped[V1], b = ped[V2], la = lineage_id[V1], lb = lineage_id[V2])], by = oocyst_id]
prs <- merge(prs[la != lb], rec[, .(oocyst_id, par1, par2, unrelated_parents)], by = "oocyst_id")

## ---------- how much two products of one oocyst share ----------
## Where two products carry the same parent's genome they are IBD; elsewhere they are IBD only
## where their parents are. Two of the four products carry each parent's genome at every locus,
## so two products of unrelated parents share a third of the genome on average, and two
## products of parents that share a fraction f of theirs share 1/3 + 2f/3.
idx <- ibd_index(ibd_pairs(run, ids = unique(c(prs$a, prs$b, prs$par1, prs$par2))))
inN <- function(x) x %in% idx$nm$pedigree_id
prs <- prs[inN(a) & inN(b) & inN(par1) & inN(par2)]
prs[, `:=`(ibd = idx$getibd(a, b), f = idx$getibd(par1, par2))]
prs[, parents := factor(fifelse(unrelated_parents, "unrelated founders", "not both founders"),
                        c("unrelated founders", "not both founders"))]
sib <- prs[, .(n = .N, oocysts = uniqueN(oocyst_id), f = mean(f), mean_ibd = mean(ibd),
               expected = mean(1 / 3 + 2 * f / 3), se = sd(ibd - 2 * f / 3) / sqrt(.N),
               pct_zero = 100 * mean(ibd < 0.01)), by = parents][order(parents)]
sib[, check := fifelse(abs(mean_ibd - expected) <= 3 * se & pct_zero == 0, "ok", "<-- CHECK")]
cat("\n== same-oocyst sib pairs, by parentage ==\n")
print(sib[, .(parents, n, oocysts, f = round(f, 3), mean_ibd = round(mean_ibd, 3),
              expected = round(expected, 3), pct_zero = round(pct_zero, 1), check)], row.names = FALSE)
cat("  expected: mean_ibd = 1/3 when the parents are unrelated and 1/3 + 2f/3 when they share f\n",
    " (f = the parents' own IBD), within 3 standard errors; on the full 14-chromosome map at the\n",
    " default rate, pct_zero = 0\n")

## ---------- independent assortment across junctions ----------
## Take products of oocysts whose parents are two distinct founders, each against its first
## parent. Founders share nothing, so a product is IBD with that parent exactly where it carries
## that parent's genome. Each chromosome assorts independently, so the parent a product carries
## at the end of one chromosome says nothing about the start of the next: the state must be
## uncorrelated across a junction (phi ~ 0) and switch at half of the 13 junctions. Inside a
## chromosome the product changes parent only at a crossover (rate r), so 100 kb either side of
## a point the state stays correlated: phi = exp(-2r x 200 kb).
ctrl <- sort(unique(pmin(pmax(round(c(bounds - PVP01) + PVP01 / 2), D + 1),
                         bounds - D - 1)))            # chromosome midpoints
prod <- merge(unique(inoc, by = c("oocyst_id", "lineage_id")),
              rec[unrelated_parents == TRUE, .(oocyst_id, parent = par1)], by = "oocyst_id")
one <- prod[seq_len(min(.N, 120))]
seg <- ibd_segments(run, one$ped, one$parent)
## products are numbered from 0, as ibd_segments() numbers its pairs. On the full map every
## product shares part of its genome with each parent; a product that is not in the tree
## sequence has no segments and is left out.
used <- sort(unique(seg$pair))
cat(sprintf("\nproducts of founder-parent oocysts used for the junction test: %d (from %d oocysts)\n",
            length(used), uniqueN(one$oocyst_id[used + 1L])))
## is product i IBD with its parent at position x?
ibd_at <- function(i, x) {
  q <- data.table(pair = i, x = x, k = seq_along(x))
  q$k %in% seg[q, on = .(pair, gstart <= x, gend > x), k, nomatch = NULL]
}
## for each product, is it IBD with its parent at (x - D) and at (x + D)? Junction positions
## against control positions the same distance apart inside a chromosome.
pos <- data.table(kind = rep(c("junction", "within-chromosome control"),
                             c(length(junctions), length(ctrl))), x = c(junctions, ctrl))
cc <- CJ(pair = used, p = seq_len(nrow(pos)))
cc[, `:=`(kind = pos$kind[p], x = pos$x[p])]
cc[, `:=`(left = ibd_at(pair, x - D), right = ibd_at(pair, x + D))]
## concordance beyond chance: phi coefficient between IBD-left and IBD-right
phi_w <- exp(-2 * R_WITHIN * 2 * D)
phi_of <- function(l, r) {
  ## as.numeric: these counts multiply to > .Machine$integer.max and would overflow
  a <- as.numeric(sum(l & r));  b <- as.numeric(sum(l & !r))
  c2 <- as.numeric(sum(!l & r)); d <- as.numeric(sum(!l & !r))
  den <- sqrt((a + b) * (c2 + d) * (a + c2) * (b + d))
  if (den > 0) (a * d - b * c2) / den else NA_real_
}
## Products of one oocyst come from one meiosis (two sisters carry the same parent at every
## centromere), so they are not independent: standard errors come from resampling oocysts.
boot_se <- function(stat, rows, oocyst, B = 400) {
  by_o <- split(rows, oocyst)
  sd(replicate(B, stat(unlist(by_o[sample(names(by_o), replace = TRUE)], use.names = FALSE))), na.rm = TRUE)
}
set.seed(1)
cc[, oocyst := one$oocyst_id[pair + 1L]]
res <- cc[, .(n = .N, pct_same = round(100 * mean(left == right), 1), phi = round(phi_of(left, right), 3),
              se = round(boot_se(function(k) phi_of(left[k], right[k]), seq_len(.N), oocyst), 3)), by = kind]
res[, expected := round(fifelse(kind == "junction", 0, phi_w), 3)]
res[, check := fifelse(abs(phi - expected) <= 3 * se, "ok", "<-- CHECK")]
cat(sprintf("\n== a product's IBD state with its parent %d kb either side of a position (n products x positions) ==\n",
            D / 1000))
print(res)
cat(sprintf(paste0("  expected: phi ~ 0 across a junction (each chromosome assorts independently);\n",
                   "  phi = exp(-2r x %d kb) = %.3f within a chromosome (a product changes parent only at a crossover)\n"),
            2 * D / 1000, phi_w))
## a product switches parent at a junction when its state differs either side of it
sw <- CJ(pair = used, x = junctions)
sw[, switched := ibd_at(pair, x - 1) != ibd_at(pair, x)]
nsw <- sw[, .(n = sum(switched)), by = pair]
ns <- nsw$n
ns_se <- boot_se(function(k) mean(ns[k]), seq_along(ns), one$oocyst_id[nsw$pair + 1L])
cat(sprintf("parent switches per product at the 13 junctions: %.2f +/- %.2f | expected 13 x 1/2 = 6.5%s\n",
            mean(ns), ns_se, if (abs(mean(ns) - 6.5) <= 3 * ns_se) "  (consistent)" else "  <-- CHECK"))

## ---------------- figure ----------------
facet_lab <- sib[, setNames(fifelse(parents == "unrelated founders",
                                    sprintf("two founder parents (unrelated): %d pairs", n),
                                    sprintf("parents not both founders: %d pairs, mean parents' IBD f = %.2f", n, f)),
                            as.character(parents))]
pA <- ggplot(prs, aes(ibd)) +
  geom_histogram(binwidth = .05, boundary = 0, colour = NA, fill = "#4393c3") +
  geom_vline(data = sib, aes(xintercept = expected), colour = "grey30", linetype = 2) +
  facet_wrap(~ parents, nrow = 1, labeller = as_labeller(facet_lab)) +
  labs(title = "How much the products of one meiosis share",
       subtitle = paste0("Two of the four products carry each parent's genome at every locus, so two products share a third of the\n",
                         "genome on average when the parents are unrelated, and 1/3 + 2f/3 when the parents share f (dashed)."),
       x = "IBD between two products of one oocyst", y = "pairs") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"),
                                        strip.text = element_text(size = 9))
pB <- ggplot(res, aes(kind, phi, fill = kind)) +
  geom_col(width = .55) + geom_hline(yintercept = 0, colour = "grey40") +
  geom_errorbar(aes(ymin = expected, ymax = expected), width = .75, linewidth = .7, colour = "grey10") +
  geom_text(aes(y = pmax(phi, expected), label = sprintf("phi = %.3f\n%.0f%% same state", phi, pct_same)),
            vjust = -0.4, size = 3.4) +
  scale_fill_manual(values = c("junction" = "#c0392b", "within-chromosome control" = "#2166ac"),
                    guide = "none") +
  coord_cartesian(ylim = c(-0.1, 1.15)) +
  labs(title = sprintf("Independent assortment: a product's parent %d kb either side", D / 1000),
       subtitle = sprintf(paste0("Products of one meiosis against one of their parents (two unrelated founders). Chromosomes assort independently,\n",
                                 "so phi ~ 0 across a junction, and a product switches parent at %.2f of the 13 junctions (expected 6.5). Inside a\n",
                                 "chromosome it changes parent only at a crossover, so phi = exp(-2r x 200 kb) = %.2f. Black lines = expected."),
                          mean(ns), phi_w),
       x = NULL, y = "phi (correlation of IBD state across the position)") +
  theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(FIGS, "ARG_junctions_and_tetrads.png"), pA / pB + plot_layout(heights = c(1.3, 1)),
       width = 10, height = 9, dpi = 140)
cat("\nwrote", file.path(FIGS, "ARG_junctions_and_tetrads.png"), "\n")
