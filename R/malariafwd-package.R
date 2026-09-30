#' malariafwd: forward simulation of P. vivax transmission and genetics
#'
#' An R interface to a SLiM (Eidos) + tskit individual-based forward simulator of
#' Plasmodium vivax. Run the model from R ([run_sim()], [run_grid()]) and analyse the
#' outputs (loaders, epidemiological + genetic metrics, ground-truth IBD via tskit).
#'
#' The model itself is bundled at `inst/slim/malariafwd_vivax.slim` (see [slim_model()]);
#' design decisions are recorded in `DECISIONS.md`.
#'
#' @keywords internal
#' @import data.table
#' @importFrom stats setNames
#' @importFrom utils combn
"_PACKAGE"
