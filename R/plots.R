## plots.R -- ggplot2 summary figures. ggplot2 is a Suggests dependency, so each helper
## checks for it and uses ggplot2:: explicitly (the package installs/runs without ggplot2).

.need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plotting needs the 'ggplot2' package: install.packages('ggplot2').", call. = FALSE)
}

#' Kaplan-Meier curve of the time to first relapse, per hypnozoite batch
#'
#' The fraction of hypnozoite batches not yet relapsed against the days since the infection that
#' laid each down, from [time_to_relapse()]. No survival-package dependency:
#' `S(t) = prod(1 - d_i/n_i)`.
#' @param x a run, or the table [time_to_relapse()] returns.
#' @param title plot title.
#' @return a ggplot object.
#' @seealso [time_to_relapse()]
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' plot_relapse_km(time_to_relapse(example_run()))
#' @export
plot_relapse_km <- function(x, title = NULL) {
  .need_ggplot()
  d <- .as_result(x, time_to_relapse, c("time", "event"))
  if (!nrow(d))
    stop("no hypnozoite batch was established in this run, so there is no time to relapse to plot.", call. = FALSE)
  N <- nrow(d)
  tab <- d[, .(d_i = sum(event), c_i = sum(1 - event)), by = time][order(time)]
  tab[, n_i := N - cumsum(shift(d_i + c_i, fill = 0))]
  tab[, surv := cumprod(1 - d_i / n_i)]
  tab <- rbind(data.table(time = 0L, d_i = 0, c_i = 0, n_i = N, surv = 1), tab[time > 0 | d_i > 0])
  ggplot2::ggplot(tab, ggplot2::aes(time, surv)) +
    ggplot2::geom_step(linewidth = .8, colour = "firebrick") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = title, x = "days since the infection that laid down the batch",
                  y = "fraction of batches not yet relapsed") +
    theme_malariafwd()
}

#' Host-state prevalence over time (stacked area)
#'
#' The share of hosts in each state, day by day, stacked; subpatent carriers and the hosts with
#' hypnozoites alone are bands of their own.
#' @param x a run, a state table from [load_state()], or the table [host_state_prevalence()] returns.
#' @return a ggplot object.
#' @seealso [host_state_prevalence()]
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' plot_host_state_prevalence(example_run())
#' @export
plot_host_state_prevalence <- function(x) {
  .need_ggplot()
  hsp <- if (is.data.frame(x) && "entity_type" %in% names(x)) host_state_prevalence(x)
         else copy(.as_result(x, host_state_prevalence, c("tick", "cls", "frac")))
  stopifnot(all(as.character(hsp$cls) %in% names(.HOST_STATE_FILL)))   # an unlisted state would draw unfilled
  hsp[, cls := factor(cls, levels = names(.HOST_STATE_FILL))]
  ggplot2::ggplot(hsp, ggplot2::aes(tick, frac, fill = cls)) + ggplot2::geom_area() +
    ggplot2::scale_fill_manual(values = .HOST_STATE_FILL) +
    ggplot2::labs(x = "day", y = "fraction of hosts", fill = NULL) +
    theme_malariafwd()
}

#' Within-host trajectory for one host
#'
#' `layout = "lineages"` (the default) stacks the individual blood-stage lineages, each in its
#' own colour, so you can see which genotypes are present at any moment, when a lineage clears
#' and when the same one returns; the hypnozoite reservoir is overlaid on a second axis.
#' `layout = "summary"` plots active COI and the reservoir as two lines instead. Each clone is
#' drawn from [episodes()]: from its onset to its end, and to the last day of the state table if it
#' was still patent when the run ended. The x axis runs to `SIM_DAYS`; on a run that went extinct the
#' lines end where `state.tsv` does.
#'
#' Each call reads the run's `state.tsv` and `events.tsv`.
#' @inheritParams load_events
#' @param host_id the host, one whole number (a run's hosts are 0 to `N_HOSTS - 1`).
#' @param layout `"lineages"` (stacked per-lineage areas) or `"summary"` (COI + reservoir lines).
#' @param palette colours cycled across lineages.
#' @return a ggplot object.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' run <- example_run()
#' n <- load_events(run)[event_type == "hypnozoite_activation", .N, by = human_id]
#' plot_host_trajectory(run, n$human_id[which.max(n$N)])   # the host with the most relapses
#' @export
plot_host_trajectory <- function(run, host_id,
                                 layout = c("lineages", "summary"),
                                 palette = c("#3B9AB2","#EBCC2A","#F21A00","#00A08A","#F98400",
                                             "#5BBCD6","#046C9A","#D69C4E","#FD6467","#5B1A18",
                                             "#E6A0C4","#7294D4","#02401B","#972D15","#9986A5")) {
  .need_ggplot()
  layout <- match.arg(layout)
  if (!is.numeric(host_id) || length(host_id) != 1L || is.na(host_id) || host_id != round(host_id) ||
      abs(host_id) > .Machine$integer.max)
    stop("`host_id` must be one whole number.", call. = FALSE)
  host_id <- as.integer(host_id)
  run <- load_run(run)
  sim_days <- .sim_days(run, "plot_host_trajectory")
  hosts <- load_state(run)
  st <- hosts[entity_id == host_id]
  if (!nrow(st))
    stop("host ", host_id, " is not a host of this run",
         if (nrow(hosts)) paste0(": host ids run from ", min(hosts$entity_id), " to ", max(hosts$entity_id)), ".",
         call. = FALSE)
  events <- load_events(run)
  if (layout == "lineages")
    return(.host_lineage_plot(st, events, host_id, palette) + ggplot2::expand_limits(x = sim_days))
  ## the host's clone onsets, from episodes(), which takes them from the event dictionary
  ev <- episodes(events)[human_id == host_id, .(tick = onset_tick, kind = .onset_kind(onset_event))]
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
                  x = "day", y = "count", colour = NULL) +
    ggplot2::expand_limits(x = sim_days) +
    theme_malariafwd()
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
#' trending at the burn-in). The two "by origin" panels share a palette. The force panel stacks
#' the new clones: imported and local reinfection (FOI), and relapse, which molFOB adds (Koepfli et
#' al. 2013). A recrudescence, a new episode but not a new clone, appears in the incidence panel
#' only. The subtitle gives [summarise_run()]'s numbers for the days after the burn-in.
#'
#' The panels run to `SIM_DAYS`: a run that goes extinct shows lines that drop to zero and stay
#' there to the end of the run, and the subtitle says on which day it went extinct. A window with no
#' case is drawn as zero.
#' @inheritParams load_events
#' @param bin window width in days.
#' @param burn_frac the share of the run's `SIM_DAYS` days marked as burn-in; the subtitle gives
#'   [summarise_run()]'s numbers for the rest.
#' @return a patchwork object (combined ggplot).
#' @references Koepfli C, et al. (2013) PLoS Negl Trop Dis 7(9):e2403 (molFOB).
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' plot_diagnostics(example_run())
#' @export
plot_diagnostics <- function(run, bin = 30, burn_frac = 0.5) {
  .need_ggplot(); .need_patchwork()
  run <- load_run(run); ev <- load_events(run); h <- load_state(run)
  rp <- .params_of(run); nh <- as.numeric(rp$N_HOSTS); nm <- as.numeric(rp$N_MOSQ)
  endt <- .sim_days(run, "plot_diagnostics")
  ## the subtitle's numbers are summarise_run()'s, over the days after the burn-in
  m <- tryCatch(.summary_row(run, ev, h, .run_window(endt, burn_frac), exclude_masked = FALSE),
                malariafwd_unfinished_run = function(e)
                  stop("plot_diagnostics() needs a finished run: ", conditionMessage(e), call. = FALSE))
  burn <- m$from_tick
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
  fill0 <- function(d, lv = olv) {
    d <- merge(data.table::CJ(w = keep(grid)$w, origin = lv), d, by = c("w", "origin"),
               all.x = TRUE)[is.na(y), y := 0]
    d[, origin := factor(origin, olv)][order(w, origin)]
  }
  inc  <- fill0(keep(cs[, .(y = .N), by = .(w, origin)]))
  ## new clones only: a recrudescence is a new episode but not a new clone
  ## per host-year of the days each window covers, as force_of_infection() gives it: the first holds
  ## days 1 to bin - 1
  fo   <- cs[origin %in% .ACQUISITIONS, .(n = .N), by = .(w, origin)][grid, on = "w", nomatch = NULL]
  fo   <- fill0(keep(fo[, .(w, origin, y = n / nh / (days / 365))]), .ACQUISITIONS)
  lin  <- merge(keep(grid)[, .(w)], cs[, .(y = uniqueN(lineage_id)), by = w], by = "w", all.x = TRUE)[is.na(y), y := 0]
  outcome <- if (!is.na(run$extinct_tick)) paste0(", extinct on day ", run$extinct_tick)
             else if (isTRUE(run$never_started)) ", never started" else ""

  gg <- ggplot2::ggplot; ae <- ggplot2::aes
  vbl <- ggplot2::geom_vline(xintercept = burn, linetype = 2, colour = "grey65", linewidth = .3)
  fillsc <- ggplot2::scale_fill_manual(values = ocol, breaks = olv, limits = olv, drop = FALSE, name = NULL)
  th <- theme_malariafwd() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom")
  lp <- function(d, col, ttl, xlab, ylab, ymax = NA) {
    ## a window with no infected host has no COI: the line breaks there
    p <- gg(d, ae(w, y)) + vbl + ggplot2::geom_line(colour = col, linewidth = .9, na.rm = TRUE) +
      ggplot2::labs(title = ttl, x = xlab, y = ylab) + th
    if (!is.na(ymax)) p + ggplot2::coord_cartesian(ylim = c(0, ymax)) else p
  }
  p1 <- lp(prev, "#2166ac", "Prevalence", NULL, "% infected", 100)
  p2 <- gg(inc, ae(w, y, fill = origin)) + vbl + ggplot2::geom_area(alpha = .85) + fillsc +
    ggplot2::guides(fill = "none") +   # one shared origin legend (kept on the force panel)
    ggplot2::labs(title = "Incidence by origin (new cases)", x = NULL, y = paste0("cases / ", bin, " d")) + th
  p3 <- gg(fo, ae(w, y, fill = origin)) + vbl + ggplot2::geom_area(alpha = .9) + fillsc +
    ggplot2::labs(title = "Force of blood-stage infection (new clones)", x = NULL, y = "new clones / host-year") + th
  p4 <- lp(coi, "#8e6bbf", "Mean COI (infected hosts)", NULL, "MOI")
  p5 <- lp(res, "#cc5500", "Hypnozoite reservoir", "day", "mean hypnozoites/host")
  p6 <- lp(lin, "#2e7d5b", "Circulating lineages", "day", "distinct / window")

  patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = "malariafwd run diagnostics",
      subtitle = sprintf(
        paste0("%.0f hosts | %.0f mosquitoes | %d d (%.1f yr)%s | biting %.3g | founders %s | burn-in %d d (dashed)\n",
               "days %d-%d: prevalence %.0f%% | COI %.2f | EIR %.1f/py | FOI %.1f/py | molFOB %.1f/py | relapse frac (episodes) %.2f"),
        nh, nm, endt, endt / 365, outcome, as.numeric(rp$BITING_RATE), rp$N_FOUNDER_LINEAGES, m$from_tick - 1L,
        m$from_tick, m$to_tick, 100 * m$prevalence_frac, m$coi_mean, m$eir_per_py, m$foi_per_py,
        m$molfob_clones_per_py, m$relapse_frac_episodes),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14))) &
    ggplot2::theme(legend.position = "bottom")
}

## Stacked per-lineage blood-stage areas for one host, reservoir overlaid on a second axis.
## Spans are the host's clones from episodes(), from onset to their first end (every way a patent
## episode can end, event_types()$ends_episode); a clone still patent when the run ended runs to
## the last day of state. An MDA round can end a clone on its first patent day (D44), a one-day
## span. The per-tick grid is completed with zeros so the areas drop to zero between episodes
## rather than being interpolated across the gap.
.host_lineage_plot <- function(st, events, host_id, palette) {
  ep <- episodes(events)[human_id == host_id]
  if (!nrow(ep)) stop("host ", host_id, " has no blood-stage episodes", call. = FALSE)
  last <- max(st$tick)
  seg <- ep[, .(pedigree_id, lineage_id, start = onset_tick, end = data.table::fifelse(is.na(end_tick), last, end_tick),
                kind = .onset_kind(onset_event))]
  stopifnot(all(seg$kind %in% names(.KIND_SHAPE)))   # an unlisted kind would draw as shape NA

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
#' Plots the compound force of blood-stage infection (**molFOB** -- every new blood-stage clone,
#' relapses included, recrudescences excluded) against the force of infection (**FOI** -- the new
#' clones from bites and imports). The distance between the lines is the relapse contribution,
#' which in *P. vivax* is typically the larger share and is invisible to any method that equates
#' new episodes with new transmission.
#'
#' Naming follows Koepfli et al. (2013): molFOB counts new blood-stage clones however they arise,
#' FOI those a bite or an import brought. A recrudescence is a new episode but not a new clone, so
#' it is in neither series. Both series are drawn in every window of the run, as zero where there
#' was no case.
#'
#' @param x a run, or the table [force_of_infection()] returns (for a run, with its defaults).
#' @param smooth span for a loess trend through each series; `NA` draws the raw series, as does a
#'   run with fewer than `5 / smooth` windows (13 at the default), too few for a trend.
#' @return a ggplot object.
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' run <- example_run()
#' plot_force_of_infection(run)
#' plot_force_of_infection(force_of_infection(run, bin = 15))
#' plot_force_of_infection(force_of_infection(run, bin = 60))   # 7 windows: no trend
#' @seealso [force_of_infection()]
#' @export
plot_force_of_infection <- function(x, smooth = 0.4) {
  .need_ggplot()
  d <- .as_result(x, force_of_infection, c("win", "origin", "rate", "is_acquisition"))
  ## POSITIVE lists, not negations. Two things have to be excluded for different reasons and
  ## a negation gets one of them wrong: force_of_infection() relabels relapses as
  ## "relapse (local)" / "relapse (imported)" by default (split_imported = TRUE), so an exact
  ## setdiff removes nothing and FOI equals molFOB, while recrudescence is a new episode but not
  ## a new clone and belongs in neither series. `is_acquisition` is computed by force_of_infection()
  ## after any relabelling, so it survives both.
  acq  <- d[is_acquisition == TRUE]
  mosq <- grep("^relapse", unique(acq$origin), value = TRUE, invert = TRUE)
  ## named ONCE: these strings are both the series labels and the factor levels, and writing
  ## them twice is how the molFOB series ended up factored to NA - drawn grey, no legend entry.
  lv <- c("molFOB (all clone acquisitions)", "FOI (bites and imports)")
  w <- rbind(
    acq[, .(rate = sum(rate), series = lv[1]), by = win],
    acq[origin %in% mosq, .(rate = sum(rate), series = lv[2]), by = win])
  ## carry zeros, for both series in every window, so a series neither breaks where a window has
  ## no episodes of one kind nor vanishes from a run with no acquisition
  w <- merge(data.table::CJ(win = unique(d$win), series = lv), w,
             by = c("win", "series"), all.x = TRUE)[is.na(rate), rate := 0]

  w[, series := factor(series, lv)]
  stopifnot(!anyNA(w$series))   # a label/level mismatch here is silent: NA draws grey, unlabelled
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = w, ggplot2::aes(win, rate, colour = series),
                       linewidth = .5, alpha = .45)
  ## a loess fit needs about five windows within its span: with fewer, loess warns and the curve
  ## says nothing the raw series does not
  if (!is.na(smooth) && data.table::uniqueN(w$win) * smooth >= 5)
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
