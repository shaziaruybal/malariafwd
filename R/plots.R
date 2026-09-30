## plots.R -- ggplot2 summary figures. ggplot2 is a Suggests dependency, so each helper
## checks for it and uses ggplot2:: explicitly (the package installs/runs without ggplot2).

.need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plotting needs the 'ggplot2' package: install.packages('ggplot2').", call. = FALSE)
}

#' Kaplan-Meier plot of time to first relapse
#'
#' No survival-package dependency: `S(t) = prod(1 - d_i/n_i)`.
#' @param ttr output of [time_to_relapse()].
#' @param title plot title.
#' @return a ggplot object.
#' @export
plot_relapse_km <- function(ttr, title = NULL) {
  .need_ggplot()
  d <- as.data.table(ttr)[order(time)]; N <- nrow(d)
  tab <- d[, .(d_i = sum(event), c_i = sum(1 - event)), by = time][order(time)]
  tab[, n_i := N - cumsum(shift(d_i + c_i, fill = 0))]
  tab[, surv := cumprod(1 - d_i / n_i)]
  ggplot2::ggplot(tab, ggplot2::aes(time, surv)) +
    ggplot2::geom_step(linewidth = .8, colour = "firebrick") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = title, x = "days since first infection", y = "fraction not yet relapsed")
}

#' Peak COI vs relapse count (boxplots)
#' @param mvr output of [moi_vs_relapse()].
#' @return a ggplot object.
#' @export
plot_moi_vs_relapse <- function(mvr) {
  .need_ggplot()
  ggplot2::ggplot(mvr[peak_moi <= 12], ggplot2::aes(factor(peak_moi), n_relapse)) +
    ggplot2::geom_boxplot(outlier.size = .5, fill = "grey85") +
    ggplot2::labs(x = "peak COI (distinct lineages ever co-active in host)",
                  y = "relapses in host")
}

#' Host-state prevalence over time (stacked area)
#' @param hsp output of [host_state_prevalence()].
#' @return a ggplot object.
#' @export
plot_host_state_prevalence <- function(hsp) {
  .need_ggplot()
  lv <- c("uninfected","active only","active + reservoir","reservoir only (invisible)")
  hsp <- copy(hsp)[, cls := factor(cls, levels = lv)]
  ggplot2::ggplot(hsp, ggplot2::aes(tick, frac, fill = cls)) + ggplot2::geom_area() +
    ggplot2::scale_fill_manual(values = c("grey85","goldenrod2","firebrick","turquoise4")) +
    ggplot2::labs(x = "day", y = "fraction of hosts", fill = NULL)
}

#' Within-host trajectory for one host
#'
#' `layout = "lineages"` (the default) stacks the individual blood-stage lineages, each in its
#' own colour, so you can see which genotypes are present at any moment, when a lineage clears
#' and when the same one returns; the hypnozoite reservoir is overlaid on a second axis.
#' `layout = "summary"` plots active COI and the reservoir as two lines instead.
#'
#' @param state a state data.table.
#' @param events an events data.table.
#' @param host_id host to plot.
#' @param layout `"lineages"` (stacked per-lineage areas) or `"summary"` (COI + reservoir lines).
#' @param palette colours cycled across lineages.
#' @return a ggplot object.
#' @export
plot_host_trajectory <- function(state, events, host_id,
                                 layout = c("lineages", "summary"),
                                 palette = c("#3B9AB2","#EBCC2A","#F21A00","#00A08A","#F98400",
                                             "#5BBCD6","#046C9A","#D69C4E","#FD6467","#5B1A18",
                                             "#E6A0C4","#7294D4","#02401B","#972D15","#9986A5")) {
  .need_ggplot()
  layout <- match.arg(layout)
  st <- as.data.table(state)[entity_type == "human" & as.integer(entity_id) == host_id]
  if (layout == "lineages") return(.host_lineage_plot(st, events, host_id, palette))
  ev <- as.data.table(events)[human_id == host_id &
        event_type %in% c("vh_inoculation","hypnozoite_activation")]
  ev <- ev[!(event_type == "vh_inoculation" & is_dormant == 1)]
  ev[, kind := fifelse(event_type == "hypnozoite_activation", "relapse", "new infection")]
  ggplot2::ggplot() +
    ggplot2::geom_line(data = st, ggplot2::aes(tick, hyp_reservoir_size, colour = "hypnozoite reservoir"), linewidth = .7) +
    ggplot2::geom_line(data = st, ggplot2::aes(tick, moi_active, colour = "active COI"), linewidth = .7) +
    ggplot2::geom_rug(data = ev, ggplot2::aes(x = tick, colour = kind), sides = "b",
                      length = grid::unit(.05, "npc")) +
    ggplot2::scale_colour_manual(values = c("active COI" = "goldenrod3", "hypnozoite reservoir" = "turquoise4",
                                            "relapse" = "firebrick", "new infection" = "grey40")) +
    ggplot2::labs(title = paste0("Within-host trajectory - host ", host_id),
                  x = "day", y = "count", colour = NULL)
}

.need_patchwork <- function() {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("plot_diagnostics() needs the 'patchwork' package: install.packages('patchwork').", call. = FALSE)
}

#' One-glance run-health dashboard
#'
#' A six-panel diagnostic for a single run: prevalence, incidence by origin, the force of
#' blood-stage infection by origin (FOI vs molFOB; see [force_of_infection()]), mean COI,
#' the hypnozoite reservoir, and circulating-lineage richness -- all over time, with the
#' burn-in marked. Surfaces the common failure modes at a glance: extinction (lines crash to
#' zero), saturation (prevalence pinned near 100%), and non-stationarity (metrics still
#' trending at the burn-in). The two "by origin" panels share a palette: imported /
#' local reinfection (= FOI) plus relapse (the extra molFOB contribution, Mueller et al. 2013).
#' @param run a run directory or `vivax_run`.
#' @param bin window width in days.
#' @param burn_frac fraction of the run marked as burn-in (a display + summary choice).
#' @return a patchwork object (combined ggplot).
#' @references Mueller I, et al. (2013) PLoS Negl Trop Dis 7(9):e2403 (molFOB).
#' @export
plot_diagnostics <- function(run, bin = 30, burn_frac = 0.5) {
  .need_ggplot(); .need_patchwork()
  st <- as.data.table(load_state(run)); ev <- as.data.table(load_events(run))
  rp <- load_params(run); nh <- as.numeric(rp$N_HOSTS); nm <- as.numeric(rp$N_MOSQ)
  h <- st[entity_type == "human"]; endt <- max(st$tick); burn <- round(endt * burn_frac)
  keep <- function(d) d[w <= endt - bin]
  h[, w := (tick %/% bin) * bin]
  prev <- keep(h[, .(y = 100 * mean(moi_active > 0)), by = w])
  coi  <- keep(h[moi_active > 0, .(y = mean(moi_active)), by = w])
  res  <- keep(h[, .(y = mean(hyp_reservoir_size)), by = w])
  olv  <- c("imported", "local reinfection", "relapse")
  ocol <- c("imported" = "#2a9d8f", "local reinfection" = "#737373", "relapse" = "#c8971f")
  cs   <- classify_cases(ev); cs[, w := (tick %/% bin) * bin]
  inc  <- keep(cs[, .(y = .N), by = .(w, origin)]); inc[, origin := factor(origin, olv)]
  fo   <- keep(cs[, .(y = .N / nh / (bin / 365)), by = .(w, origin)]); fo[, origin := factor(origin, olv)]
  lin  <- keep(cs[, .(y = uniqueN(lineage_id)), by = w])
  m    <- summarise_run(run)
  py   <- max((endt - burn) / 365, .001)
  foi  <- nrow(cs[origin != "relapse" & tick >= burn]) / nh / py
  mol  <- nrow(cs[tick >= burn]) / nh / py

  gg <- ggplot2::ggplot; ae <- ggplot2::aes
  vbl <- ggplot2::geom_vline(xintercept = burn, linetype = 2, colour = "grey65", linewidth = .3)
  fillsc <- ggplot2::scale_fill_manual(values = ocol, breaks = olv, limits = olv, drop = FALSE, name = NULL)
  th <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
  lp <- function(d, col, ttl, xlab, ylab, ymax = NA) {
    p <- gg(d, ae(w, y)) + vbl + ggplot2::geom_line(colour = col, linewidth = .9) +
      ggplot2::labs(title = ttl, x = xlab, y = ylab) + th
    if (!is.na(ymax)) p + ggplot2::coord_cartesian(ylim = c(0, ymax)) else p
  }
  p1 <- lp(prev, "#2166ac", "Prevalence", NULL, "% infected", 100)
  p2 <- gg(inc, ae(w, y, fill = origin)) + vbl + ggplot2::geom_area(alpha = .85) + fillsc +
    ggplot2::guides(fill = "none") +   # one shared origin legend (kept on the force panel)
    ggplot2::labs(title = "Incidence by origin (new cases)", x = NULL, y = paste0("cases / ", bin, " d")) + th
  p3 <- gg(fo, ae(w, y, fill = origin)) + vbl + ggplot2::geom_area(alpha = .9) + fillsc +
    ggplot2::labs(title = "Force of blood-stage infection by origin", x = NULL, y = "new clones / host-year") + th
  p4 <- lp(coi, "#8e6bbf", "Mean COI (infected hosts)", NULL, "MOI")
  p5 <- lp(res, "#cc5500", "Hypnozoite reservoir", "day", "mean batch/host")
  p6 <- lp(lin, "#2e7d5b", "Circulating lineages", "day", "distinct / window")

  patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = "malariafwd run diagnostics",
      subtitle = sprintf(
        paste0("%.0f hosts | %.0f mosquitoes | %d d (%.1f yr) | biting %.3g | founders %s | burn-in %d d (dashed)\n",
               "prev %.0f%% | COI %.2f | FOI %.1f/py (import+local) | molFOB %.1f/py (+relapse) | relapse frac %.2f"),
        nh, nm, endt, endt / 365, as.numeric(rp$BITING_RATE), rp$N_FOUNDER_LINEAGES, burn,
        100 * m$prev, m$COI, foi, mol, m$relapse_frac),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))) &
    ggplot2::theme(legend.position = "bottom")
}

## Stacked per-lineage blood-stage areas for one host, reservoir overlaid on a second axis.
## Spans come from each episode's onset event and its matching clearance (any of the three
## ways an episode can end); episodes still active at the end run to the last tick. The
## per-tick grid is completed with zeros so the areas drop to zero between episodes rather
## than being interpolated across the gap.
.host_lineage_plot <- function(st, events, host_id, palette) {
  ev <- as.data.table(events)
  last <- max(ev$tick)
  onset <- ev[human_id == host_id & is_dormant == 0 &
              event_type %in% c("vh_inoculation", "vh_inoculation_import",
                                "vh_inoculation_seed", "hypnozoite_activation"),
              .(pedigree_id, lineage_id = as.character(lineage_id), start = tick,
                kind = fifelse(event_type == "hypnozoite_activation", "relapse", "new infection"))]
  if (!nrow(onset)) stop("host ", host_id, " has no blood-stage episodes", call. = FALSE)
  clears <- event_types()$event_type[event_types()$ends_episode]
  ends <- ev[event_type %in% clears, .(pedigree_id, end = tick)]
  seg <- merge(onset, ends, by = "pedigree_id", all.x = TRUE)[is.na(end), end := last][end >= start]

  lv <- seg[order(start), unique(lineage_id)]
  grid <- seg[, .(tick = seq(start, end)), by = .(pedigree_id, lineage_id)]
  stk <- merge(data.table::CJ(tick = seq(min(seg$start), last), lineage_id = lv),
               grid[, .(n = .N), by = .(tick, lineage_id)],
               by = c("tick", "lineage_id"), all.x = TRUE)
  stk[is.na(n), n := 0L][, lineage := factor(lineage_id, lv)]

  pal <- stats::setNames(
    if (length(lv) <= length(palette)) palette[seq_along(lv)]
    else grDevices::colorRampPalette(palette)(length(lv)), lv)
  ymax <- max(1, stk[, sum(n), by = tick]$V1)
  res <- st[, .(tick, r = hyp_reservoir_size)]
  sf <- ymax / max(1, max(res$r))

  ggplot2::ggplot() +
    ggplot2::geom_area(data = stk, ggplot2::aes(tick, n, fill = lineage),
                       position = "stack", colour = "white", linewidth = .12) +
    ggplot2::geom_line(data = res, ggplot2::aes(tick, r * sf), linetype = 2, linewidth = .6,
                       colour = "grey20") +
    ggplot2::geom_point(data = seg, ggplot2::aes(start, -0.06 * ymax, shape = kind),
                        size = 1.7, colour = "grey20") +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_shape_manual(values = c("new infection" = 24, "relapse" = 21), name = NULL) +
    ggplot2::scale_y_continuous(
      "blood-stage lineages",
      sec.axis = ggplot2::sec_axis(~ . / sf, name = "hypnozoite reservoir")) +
    ggplot2::labs(x = "day") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   axis.title.y.right = ggplot2::element_text(colour = "grey30"))
}
