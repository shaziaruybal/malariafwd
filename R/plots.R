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
    ggplot2::labs(title = title, x = "days since first infection", y = "fraction not yet relapsed") +
    theme_malariafwd()
}

#' Peak COI vs relapse count
#'
#' Boxes are drawn only for COI levels with at least `min_n` hosts; sparser levels show their
#' raw points alone. Every level is labelled with its host count.
#'
#' @param mvr output of [moi_vs_relapse()].
#' @param min_n minimum hosts in a COI level before a box is drawn over it.
#' @return a ggplot object.
#' @export
plot_moi_vs_relapse <- function(mvr, min_n = 5) {
  .need_ggplot()
  d <- as.data.table(mvr)[peak_moi <= 12]
  n_by <- d[, .(n = .N), by = peak_moi]
  ## A box drawn over two hosts is not a distribution. Sparse COI levels (the tail is always
  ## thin) are shown as the raw points only, and every level carries its n, so the reader can
  ## see what the box is built on.
  d <- merge(d, n_by, by = "peak_moi")
  d[, grp := factor(peak_moi)]
  ggplot2::ggplot(d, ggplot2::aes(grp, n_relapse)) +
    ggplot2::geom_boxplot(data = d[n >= min_n], outlier.shape = NA, fill = "grey88",
                          colour = "grey35", width = .6, linewidth = .4) +
    ggplot2::geom_jitter(width = .14, height = .15, alpha = .35, size = .9, colour = "grey20") +
    ggplot2::geom_text(data = n_by, ggplot2::aes(factor(peak_moi), Inf, label = paste0("n = ", n)),
                       inherit.aes = FALSE, vjust = 1.4, size = 3, colour = "grey30") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.04, .14))) +
    ggplot2::labs(x = "peak COI (distinct lineages ever co-active in host)",
                  y = "relapses in host") +
    theme_malariafwd()
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
    ggplot2::labs(x = "day", y = "fraction of hosts", fill = NULL) +
    theme_malariafwd()
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
  ## DERIVED from the dictionary, not restated. A hard-coded list here has gone stale twice:
  ## once dropping imported onsets, and it would have dropped recrudescent ones too.
  ev <- as.data.table(events)[human_id == host_id & event_type %in% .onset_events()]
  ## is_dormant == 0 across ALL routes, not just vh_inoculation: a dormant row is a hypnozoite
  ## being seeded, which is not an episode onset, and keying the exclusion to one event type
  ## let dormant IMPORT inoculations through as onsets (5 spurious marks on one host).
  ev <- ev[is_dormant == 0]
  ev[, kind := data.table::fcase(event_type == "hypnozoite_activation", "relapse",
                                 event_type == "recrudescence", "recrudescence",
                                 default = "new infection")]
  ggplot2::ggplot() +
    ggplot2::geom_line(data = st, ggplot2::aes(tick, hyp_reservoir_size, colour = "hypnozoite reservoir"), linewidth = .7) +
    ggplot2::geom_line(data = st, ggplot2::aes(tick, moi_active, colour = "active COI"), linewidth = .7) +
    ggplot2::geom_rug(data = ev, ggplot2::aes(x = tick, colour = kind), sides = "b",
                      length = grid::unit(.05, "npc")) +
    ggplot2::scale_colour_manual(values = c("active COI" = "goldenrod3",
                                            "hypnozoite reservoir" = "turquoise4",
                                            "relapse" = "firebrick", "new infection" = "grey40",
                                            "recrudescence" = "#b4472e")) +
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
#' local reinfection (= FOI) plus relapse (the extra molFOB contribution, Koepfli et al. 2013),
#' plus recrudescence where a run produces any. Only the first three are clone ACQUISITIONS, so
#' the FOI panel sums those and not the fourth.
#' @param run a run directory or `vivax_run`.
#' @param bin window width in days.
#' The panels run to `SIM_DAYS`: a run that goes extinct shows lines that drop to zero and stay
#' there to the end of the run, and the subtitle says on which day it went extinct.
#' @param burn_frac the fraction of `SIM_DAYS` marked as burn-in (a display and summary choice).
#' @return a patchwork object (combined ggplot).
#' @references Koepfli C, et al. (2013) PLoS Negl Trop Dis 7(9):e2403 (molFOB).
#' @export
plot_diagnostics <- function(run, bin = 30, burn_frac = 0.5) {
  .need_ggplot(); .need_patchwork()
  st <- as.data.table(load_state(run)); ev <- as.data.table(load_events(run))
  rp <- .params_of(run); nh <- as.numeric(rp$N_HOSTS); nm <- as.numeric(rp$N_MOSQ)
  h <- st[entity_type == "human"]; endt <- .sim_days(run, "plot_diagnostics"); burn <- round(endt * burn_frac)
  ## every window of the run, to SIM_DAYS, with the days it covers; state.tsv stops at an extinction
  grid <- data.table::data.table(w = seq(0, (endt %/% bin) * bin, by = bin))
  grid[, days := pmin(w + bin - 1, endt) - pmax(w, 1) + 1]
  keep <- function(d) d[w <= endt - bin]
  h[, w := (tick %/% bin) * bin]
  ## prevalence and the reservoir over every host-day of a window, so that the days after an
  ## extinction count as days with no infected host; COI is NA where no host is infected
  hd <- merge(grid, h[, .(inf = sum(moi_active > 0), hyp = sum(hyp_reservoir_size)), by = w], by = "w", all.x = TRUE)
  hd[is.na(inf), `:=`(inf = 0, hyp = 0)]
  prev <- keep(hd[, .(w, y = 100 * inf / (nh * days))])
  res  <- keep(hd[, .(w, y = hyp / (nh * days))])
  coi  <- keep(merge(grid[, .(w)], h[moi_active > 0, .(y = mean(moi_active)), by = w], by = "w", all.x = TRUE))
  olv  <- .ORIGINS
  ocol <- c("imported" = "#2a9d8f", "local reinfection" = "#737373",
            "relapse" = "#c8971f", "recrudescence" = "#b4472e")
  cs   <- classify_cases(ev); cs[, w := (tick %/% bin) * bin]
  ## Zero-fill the (window x origin) grid before stacking. geom_area joins consecutive ROWS,
  ## so a window with no cases of one origin is not drawn as zero - it is interpolated
  ## straight across, inventing cases in exactly the empty windows you look at to ask whether
  ## transmission stopped. plot_force_of_infection() below carries zeros for the same reason.
  ## Every window of the run is on the grid, so a window with no cases at all is a zero too.
  fill0 <- function(d) {
    d <- merge(data.table::CJ(w = keep(grid)$w, origin = olv), d, by = c("w", "origin"),
               all.x = TRUE)[is.na(y), y := 0]
    d[, origin := factor(origin, olv)][order(w, origin)]
  }
  inc  <- fill0(keep(cs[, .(y = .N), by = .(w, origin)]))
  fo   <- fill0(keep(cs[, .(y = .N / nh / (bin / 365)), by = .(w, origin)]))
  lin  <- merge(keep(grid)[, .(w)], cs[, .(y = uniqueN(lineage_id)), by = w], by = "w", all.x = TRUE)[is.na(y), y := 0]
  m    <- summarise_run(run)
  r    <- if (inherits(run, "vivax_run")) run else load_run(run)
  outcome <- if (!is.na(r$extinct_tick)) paste0(", extinct on day ", r$extinct_tick)
             else if (isTRUE(r$never_started)) ", never started" else ""
  py   <- max((endt - burn) / 365, .001)
  ## a positive list, not a negation: `origin != "relapse"` silently admitted recrudescence,
  ## which is a new episode but not a new clone and so belongs in neither FOI nor molFOB.
  foi  <- nrow(cs[origin %in% c("imported", "local reinfection") & tick >= burn]) / nh / py
  mol  <- nrow(cs[origin %in% .ACQUISITIONS & tick >= burn]) / nh / py

  gg <- ggplot2::ggplot; ae <- ggplot2::aes
  vbl <- ggplot2::geom_vline(xintercept = burn, linetype = 2, colour = "grey65", linewidth = .3)
  fillsc <- ggplot2::scale_fill_manual(values = ocol, breaks = olv, limits = olv, drop = FALSE, name = NULL)
  th <- theme_malariafwd() +
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
        paste0("%.0f hosts | %.0f mosquitoes | %d d (%.1f yr)%s | biting %.3g | founders %s | burn-in %d d (dashed)\n",
               "prev %.0f%% | COI %.2f | FOI %.1f/py (import+local) | molFOB %.1f/py (+relapse) | relapse frac %.2f"),
        nh, nm, endt, endt / 365, outcome, as.numeric(rp$BITING_RATE), rp$N_FOUNDER_LINEAGES, burn,
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
  onset <- ev[human_id == host_id & is_dormant == 0 & event_type %in% .onset_events(),
              .(pedigree_id, lineage_id = as.character(lineage_id), start = tick,
                kind = data.table::fcase(event_type == "hypnozoite_activation", "relapse",
                                         event_type == "recrudescence", "recrudescence",
                                         default = "new infection"))]
  if (!nrow(onset)) stop("host ", host_id, " has no blood-stage episodes", call. = FALSE)
  stopifnot(all(onset$kind %in% names(.KIND_SHAPE)))   # an unlisted kind would draw as shape NA
  clears <- event_types()$event_type[event_types()$ends_episode]
  ends <- ev[event_type %in% clears, .(pedigree_id, end = tick)]
  ## end >= start: an MDA round can clear a clone on its first patent day (D44), a one-day
  ## episode that classify_cases() counts as a case too
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
    ## The legend shows only the kinds this host actually had - an untreated run has no
    ## recrudescences and should not advertise one. The guard against a kind with no shape is
    ## the stopifnot above, not drop = FALSE: ggplot2 maps an unlisted level to shape = NA and
    ## draws nothing, silently, so it has to be caught before it reaches the scale.
    ggplot2::scale_shape_manual(values = .KIND_SHAPE, name = NULL) +
    ggplot2::scale_y_continuous(
      "blood-stage lineages",
      sec.axis = ggplot2::sec_axis(~ . / sf, name = "hypnozoite reservoir")) +
    ggplot2::labs(x = "day") +
    theme_malariafwd() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   axis.title.y.right = ggplot2::element_text(colour = "grey30"))
}

#' Force of blood-stage infection: molFOB vs FOI
#'
#' Plots the compound force of blood-stage infection (**molFOB** -- every new blood-stage
#' episode, relapses included) against the mosquito-borne force of infection (**FOI** -- new
#' inoculations only). The distance between the lines is the relapse contribution, which in
#' *P. vivax* is typically the larger share and is invisible to any method that equates new
#' episodes with new transmission.
#'
#' Naming follows Koepfli et al. (2013): molFOB counts blood-stage infections however they arise,
#' FOI counts mosquito-borne infections only.
#'
#' @param foi output of [force_of_infection()].
#' @param smooth span for a loess trend through each series; `NA` draws the raw series.
#' @return a ggplot object.
#' @examples
#' \dontrun{plot_force_of_infection(force_of_infection(run))}
#' @export
plot_force_of_infection <- function(foi, smooth = 0.4) {
  .need_ggplot()
  d <- as.data.table(foi)
  ## POSITIVE lists, not negations. Two things have to be excluded for different reasons and
  ## a negation gets one of them wrong: force_of_infection(split_imported = "auto") relabels
  ## relapses as "relapse (local)" / "relapse (imported)" (so an exact setdiff removes
  ## nothing and FOI equals molFOB), while recrudescence is a new episode but not a new clone
  ## and belongs in neither series. `is_acquisition` is computed by force_of_infection()
  ## after any relabelling, so it survives both.
  acq  <- d[is_acquisition == TRUE]
  mosq <- grep("^relapse", unique(acq$origin), value = TRUE, invert = TRUE)
  ## named ONCE: these strings are both the series labels and the factor levels, and writing
  ## them twice is how the molFOB series ended up factored to NA - drawn grey, no legend entry.
  lv <- c("molFOB (all clone acquisitions)", "FOI (mosquito-borne only)")
  w <- rbind(
    acq[, .(rate = sum(rate), series = lv[1]), by = win],
    acq[origin %in% mosq, .(rate = sum(rate), series = lv[2]), by = win])
  ## carry zeros so a series does not break where a window has no episodes of one kind
  w <- merge(data.table::CJ(win = unique(d$win), series = unique(w$series)), w,
             by = c("win", "series"), all.x = TRUE)[is.na(rate), rate := 0]

  w[, series := factor(series, lv)]
  stopifnot(!anyNA(w$series))   # a label/level mismatch here is silent: NA draws grey, unlabelled
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = w, ggplot2::aes(win, rate, colour = series),
                       linewidth = .5, alpha = .45)
  if (!is.na(smooth))
    p <- p + ggplot2::geom_smooth(data = w, ggplot2::aes(win, rate, colour = series),
                                  method = "loess", span = smooth, se = FALSE,
                                  linewidth = 1, formula = y ~ x)
  p +
    ggplot2::scale_colour_manual(values = stats::setNames(c("#b2182b", "#2166ac"), lv), name = NULL) +
    ggplot2::labs(x = "day", y = "episodes per person-year",
                  subtitle = "the distance between the lines is the relapse contribution") +
    theme_malariafwd() +
    ggplot2::theme(legend.position = "top")
}
