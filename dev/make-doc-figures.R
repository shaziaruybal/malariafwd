## Makes the figures that README.md and MODEL_DESCRIPTION.md embed, and the numbers the model
## description quotes beside them, from seeded runs. Run it from the package's root, with SLiM 5.2
## on the PATH and Python with tskit (see check_setup()):
##
##   Rscript dev/make-doc-figures.R
##
## The runs go one at a time (about 10 minutes): run_grid()'s workers would load the installed
## package rather than this source. Set MALARIAFWD_DOC_RUNS to a directory to keep the runs there:
## a second call then redraws from them, as run_grid() skips the cells that already hold their run
## (and refuses ones a different model made).
## The representative configuration of the model description's opening runs for seeds 1 to 10, and
## again with the observation layer of its detection example. Each quoted number is the mean over
## the ten seeds with its range, over the second half of each run, and is written to
## man/figures/representative-run.tsv; the figures draw seed 1, the karyotype the first seed that
## holds every kind of pair it shows. The schematic's PNG is rendered from its SVG with rsvg, or with
## qlmanage on macOS.

pkgload::load_all(".", quiet = TRUE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
SEEDS <- 1:10
FIG_SEED <- 1L
out <- Sys.getenv("MALARIAFWD_DOC_RUNS")
if (!nzchar(out)) {
  out <- file.path(tempdir(), "doc-figures")
  unlink(out, recursive = TRUE)
}
figs <- file.path("man", "figures")
save_png <- function(p, file, w, h) {
  ggsave(file.path(figs, file), p, width = w, height = h, dpi = 150, bg = "white",
         device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png")
  cat("wrote", file.path(figs, file), "\n")
}

## ---- the runs ------------------------------------------------------------------------------
base <- list(N_HOSTS = 200L, N_MOSQ = 400L, SIM_DAYS = 1095L, BITING_RATE = 0.08,
             N_FOUNDER_HOSTS = 40L, N_FOUNDER_LINEAGES = 10L, NU_HYPNOZOITE_BATCH = 3)
sets <- list(representative = do.call(vivax_params, base),
             detection = do.call(vivax_params, c(base, list(P_DETECT_PRIMARY = 0.6, P_DETECT_RELAPSE = 0.15))))
g <- run_grid(sets, seeds = SEEDS, outdir_base = out)
stopifnot(all(g$ok), all(is.na(g$extinct_tick)))
run_of <- function(set, seed) load_run(g$outdir[g$set == set & g$seed == seed])

## ---- the quoted numbers --------------------------------------------------------------------
second_half <- function(run) c(from = round(run$params$SIM_DAYS * 0.5) + 1, to = run$params$SIM_DAYS)
per_seed <- rbindlist(lapply(SEEDS, function(s) {
  r <- run_of("representative", s)
  m <- summarise_run(r)
  ev <- load_events(r)
  w0 <- second_half(r)
  up <- ev[tick >= w0[["from"]] & tick <= w0[["to"]] &
           event_type %in% c("hv_uptake_recombinant", "hv_uptake_clonal"), .N, by = event_type]
  n_up <- function(t) if (t %in% up$event_type) up[event_type == t, N] else 0L
  d <- run_of("detection", s)
  w <- second_half(d)
  cs <- classify_cases(load_events(d))[tick >= w[["from"]] & tick <= w[["to"]] &
                                       origin %in% c("relapse", "local reinfection")]
  data.table(seed = s,
             prevalence = m$prevalence_frac, coi_mean = m$coi_mean, polyclonal = m$polyclonal_frac,
             hypnozoites_per_host = m$hypnozoites_per_host,
             relapse_frac = m$relapse_frac_acquisitions, relapses_per_py = m$relapses_per_py_at_risk,
             molfob_per_py = m$molfob_clones_per_py, foi_per_py = m$foi_per_py,
             recombinant_uptakes = n_up("hv_uptake_recombinant"), clonal_uptakes = n_up("hv_uptake_clonal"),
             recombinant_frac = n_up("hv_uptake_recombinant") /
               (n_up("hv_uptake_recombinant") + n_up("hv_uptake_clonal")),
             relapse_share_true = cs[, mean(origin == "relapse")],
             relapse_share_detected = cs[detected == 1, mean(origin == "relapse")])
}))
long <- melt(per_seed[, lapply(.SD, as.numeric)], id.vars = "seed", variable.name = "metric")
tab <- long[, .(mean = mean(value), min = min(value), max = max(value),
                seed1 = value[seed == FIG_SEED]), by = metric]
fwrite(tab, file.path(figs, "representative-run.tsv"), sep = "\t")
cat("wrote", file.path(figs, "representative-run.tsv"), "\n")
print(tab[, .(metric, mean = signif(mean, 3), min = signif(min, 3), max = signif(max, 3))], row.names = FALSE)

## ---- the representative run ----------------------------------------------------------------
run <- run_of("representative", FIG_SEED)
ev <- load_events(run)
st <- load_state(run)
## A: the host with the most relapses (the first of them by id)
rel_by_host <- ev[event_type == "hypnozoite_activation", .N, by = human_id][order(-N, human_id)]
host <- rel_by_host$human_id[1]
pA <- plot_host_trajectory(run, host)
pB <- plot_host_state_prevalence(run)
## C: each host ever infected, its peak COI against the relapses it had
peak <- st[, .(peak_coi = max(moi_active)), by = .(human_id = entity_id)][peak_coi >= 1]
peak <- merge(peak, rel_by_host[, .(human_id, relapses = N)], by = "human_id", all.x = TRUE)
peak[is.na(relapses), relapses := 0L]
counts <- peak[, .N, by = peak_coi]
pC <- ggplot(peak, aes(factor(peak_coi), relapses)) +
  geom_boxplot(data = peak[peak_coi %in% counts[N >= 5, peak_coi]], outlier.shape = NA,
               fill = "grey88", colour = "grey30", width = 0.55) +
  ## seeded, so that a redraw puts every point where it was
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1),
             size = 0.8, alpha = 0.45) +
  geom_text(data = counts, aes(factor(peak_coi), Inf, label = paste("n =", N)), vjust = 1.4,
            size = 3, colour = "grey30", inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Peak COI (distinct lineages ever co-active in the host)", y = "Relapses in the host") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
pD <- plot_force_of_infection(run) + labs(subtitle = NULL)   # the caption says what the gap is
pE <- plot_relatedness_ladder(run, style = "bars")
rep_fig <- pA / (pB | pC) / (pD | pE) + plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1, 1.15))
save_png(rep_fig, "README-representative-run.png", 12, 13.5)

## ---- the IBD karyotype ---------------------------------------------------------------------
kinds <- c("relapse vs its origin", "same-oocyst sibs", "distinct founders")
karyo_seed <- NA_integer_
for (s in SEEDS) {
  tp <- ibd_truth_pairs(run_of("representative", s))
  if (all(kinds %in% as.character(tp$cls))) { karyo_seed <- s; break }
}
stopifnot(!is.na(karyo_seed))
krun <- run_of("representative", karyo_seed)
tp <- ibd_truth_pairs(krun)[as.character(cls) %in% kinds]
## for each kind, the pair whose IBD is nearest the kind's mean in this run (the first by ids)
pick <- tp[, .SD[order(abs(ibd - mean(ibd)), a, b)][1], by = cls]
pick <- pick[match(kinds, as.character(cls))]
labels <- sprintf("%s\nIBD %.2f", c("relapse vs its origin", "two products of one oocyst", "distinct founders"),
                  pick$ibd)
kfig <- plot_ibd_pair(krun, pick$a, pick$b, labels = labels)
save_png(kfig, "README-ibd-karyotype.png", 12, 7.5)

## ---- the observation layer -----------------------------------------------------------------
drun <- run_of("detection", FIG_SEED)
dcs <- classify_cases(load_events(drun))[origin %in% c("relapse", "local reinfection")]
## 30-day windows, each as a rate per 30 days: the last holds what is left of the run (15 days), and
## a window without an episode is a zero
dcs[, win := (tick - 1L) %/% 30L * 30L]
sim_days <- drun$params$SIM_DAYS
wins <- data.table(win = seq(0L, sim_days - 1L, by = 30L))
wins[, days := pmin(win + 30L, sim_days) - win]
per30 <- function(d, view) merge(wins, d[, .N, by = win], by = "win", all.x = TRUE)[
  is.na(N), N := 0L][, .(win, view = view, episodes = N / days * 30)]
series <- rbind(per30(dcs, "true"), per30(dcs[detected == 1], "detected"))
series[, view := factor(view, c("true", "detected"))]
pA2 <- ggplot(series, aes(win, episodes, colour = view)) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(true = "grey25", detected = "#b2182b"), name = NULL) +
  labs(x = "Day", y = "Blood-stage episodes per 30 d",
       subtitle = "true burden against what the surveillance system sees") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
w <- second_half(drun)
mix <- dcs[tick >= w[["from"]] & tick <= w[["to"]]]
mix <- rbind(mix[, .(view = "true", origin)], mix[detected == 1, .(view = "detected", origin)])
mix[, view := factor(view, c("true", "detected"))]
pB2 <- ggplot(mix, aes(view, fill = origin)) + geom_bar(position = "fill", width = 0.6, colour = "grey40") +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
  scale_fill_manual(values = c(relapse = "#C99722", `local reinfection` = "#4393C3"), name = NULL) +
  labs(x = NULL, y = "Share of episodes", subtitle = "detection shifts the observed mix away from relapses") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save_png(pA2 + pB2 + plot_annotation(tag_levels = "A"), "README-detection.png", 11, 3.8)
cat(sprintf("figures: seed %d; karyotype: seed %d (host %d in panel A)\n", FIG_SEED, karyo_seed, host))

## ---- the schematic -------------------------------------------------------------------------
svg <- file.path(figs, "model-structure.svg")
png <- file.path(figs, "model-structure.png")
if (requireNamespace("rsvg", quietly = TRUE)) {
  rsvg::rsvg_png(svg, png, width = 2400)
  cat("wrote", png, "\n")
} else if (Sys.info()[["sysname"]] == "Darwin" && nzchar(Sys.which("qlmanage"))) {
  source(file.path("dev", "make-readme-figures.R"))
} else {
  message("model-structure.png not rendered: install rsvg, or run on macOS")
}
