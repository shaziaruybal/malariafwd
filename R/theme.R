## theme.R -- one place that controls how every malariafwd plot looks.

## The theme every plot_*() function applies, based on theme_bw(), so that figures from different
## functions compose in a panel. `base_size` is the base font size; `grid = FALSE` removes all
## gridlines, for the karyotype and other plots where the grid competes with the data.
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
