## plots.R -- ggplot2 summary figures. ggplot2 is a Suggests dependency, so each helper
## checks for it and uses ggplot2:: explicitly (the package installs/runs without ggplot2).

.need_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plotting needs the 'ggplot2' package: install.packages('ggplot2').", call. = FALSE)
}

#' Kaplan-Meier plot of time to first relapse
#'
#' No survival-package dependency: `S(t) = prod(1 - d_i/n_i)`.
#' @param ttr output of [compute_time_to_relapse()].
#' @param title plot title.
#' @return a ggplot object.
#' @export
plot_relapse_km <- function(ttr, title = "Time to first relapse (host-level Kaplan-Meier)") {
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
    ggplot2::labs(title = "Polyclonal infections predict more relapses",
                  x = "peak COI (distinct lineages ever co-active in host)", y = "relapses in host")
}

#' Host-state prevalence over time (stacked area; the reservoir-only class is invisible)
#' @param hsp output of [host_state_prevalence()].
#' @return a ggplot object.
#' @export
plot_host_state_prevalence <- function(hsp) {
  .need_ggplot()
  lv <- c("uninfected","active only","active + reservoir","reservoir only (invisible)")
  hsp <- copy(hsp)[, cls := factor(cls, levels = lv)]
  ggplot2::ggplot(hsp, ggplot2::aes(tick, frac, fill = cls)) + ggplot2::geom_area() +
    ggplot2::scale_fill_manual(values = c("grey85","goldenrod2","firebrick","turquoise4")) +
    ggplot2::labs(title = "Host states over time - the reservoir-only class is invisible to surveillance",
                  x = "day", y = "fraction of hosts", fill = NULL)
}

#' Within-host trajectory: active COI + hypnozoite reservoir, with event markers
#' @param state a state data.table.
#' @param events an events data.table.
#' @param host_id host to plot.
#' @return a ggplot object.
#' @export
plot_host_trajectory <- function(state, events, host_id) {
  .need_ggplot()
  st <- as.data.table(state)[entity_type == "human" & as.integer(entity_id) == host_id]
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
