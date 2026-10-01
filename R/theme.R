## theme.R -- one place that controls how every malariafwd plot looks.

#' malariafwd plot theme
#'
#' A `theme_bw()`-based theme applied by every `plot_*()` function in the package, so figures
#' from different functions compose cleanly in a panel. Set it once here rather than per plot.
#'
#' @param base_size base font size.
#' @param grid `TRUE` keeps light major gridlines; `FALSE` removes all gridlines (used by the
#'   karyotype and other plots where the grid competes with the data).
#' @return a ggplot2 theme object.
#' @examples
#' \dontrun{plot_relapse_km(time_to_relapse(ev)) + theme_malariafwd(base_size = 13)}
#' @export
theme_malariafwd <- function(base_size = 11, grid = TRUE) {
  .need_ggplot()
  th <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.border     = ggplot2::element_rect(colour = "grey60", fill = NA, linewidth = .4),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey60", linewidth = .4),
      strip.text       = ggplot2::element_text(size = base_size - 1),
      plot.title       = ggplot2::element_text(size = base_size, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1, colour = "grey30"),
      legend.key       = ggplot2::element_blank())
  if (!grid) th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  th
}
