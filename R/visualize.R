#' Visualize a Matrix with Group Boundaries
#'
#' @description
#' Visualize a matrix (e.g., a precision, covariance, or adjacency matrix) as a
#' heatmap with optional group-based reordering and dashed boundary lines
#' separating group blocks.
#'
#' @param mat A p-by-p matrix.
#'
#' @param membership An integer vector specifying the group membership.
#' The length of \code{membership} must be consistent with the dimension p.
#'
#' @param reorder_by_group A boolean (default = \code{FALSE}) specifying whether
#' to reorder both rows and columns of \code{mat} according to \code{membership},
#' such that variables from the same group are contiguous.
#'
#' @param colors A vector of colors specifying an n-color gradient scale for
#' the fill aesthetics.
#'
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom scales rescale
#'
#' @return
#' A \code{ggplot2} heatmap showing the matrix entries. Dashed lines indicate
#' group boundaries. The plot title also reports matrix dimension and sparsity.
#'
#' @examples
#' ## reproducibility for everything
#' set.seed(1234)
#'
#' ## user-defined sampler
#' my_gamma <- function(n) {
#'   rgamma(n, shape = 20, scale = 5)
#' }
#'
#' ## block-structured precision matrix based on SBM
#' sim <- gen_prec_sbm(p = 100, K = 5,
#'                     within.prob = 0.5, between.prob = 0.05,
#'                     weight.dists = list(my_gamma, "unif"),
#'                     weight.paras = list(NULL, c(min = 0, max = 1)),
#'                     min.eig = 0.1)
#'
#' ## visualization
#' visualize(sim$Omega, sim$membership)
#'
#' @export

visualize <- function(mat, membership, reorder_by_group = FALSE,
                      colors = colorRampPalette(
                        c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                          "yellow", "#FF7F00", "red", "#7F0000"))(256)) {

  p <- ncol(mat)
  if (length(membership) != p) {
    stop(sprintf("Length of 'membership' (%d) must equal the matrix dimension p (%d).",
                 length(membership), p))
  }

  ## optionally reorder rows/cols by group membership
  if (reorder_by_group) {
    o <- order(membership)
    mat <- mat[o, o]
    membership <- membership[o]
  }

  ## compute group sizes and boundary positions for dashed lines
  grp_sizes <- table(membership)
  cuts <- cumsum(grp_sizes)
  bnds <- cuts[-length(cuts)] + 0.5 ## boundaries between groups
  y_bnds <- p - bnds + 1  ## flipped y, note: scale_y_discrete(limits = rev)

  ## declare
  Col <- Row <- value <- NULL

  ## plot data
  labs <- paste0("V", seq_len(p))
  plotData <- data.frame(
    Row = factor(rep(labs, times = p),  levels = labs),
    Col = factor(rep(labs, each = p), levels = labs),
    value = as.vector(mat),
    check.names = FALSE
  )
  # plotData <- as.data.frame(mat) %>%
  #   mutate(Row = factor(paste0("V", rownames(.)),
  #                       levels = paste0("V", rownames(.)))) %>%
  #   pivot_longer(-Row, names_to = "Col", values_to = "value") %>%
  #   mutate(Col = factor(Col, levels = levels(Row)))

  ## sparsity (proportion of zero entries)
  sparsity <- round(sum(mat == 0) / length(mat), 4)

  ## zero -> NA for better white rendering
  plotData$value[plotData$value == 0] <- NA

  ## color scaling range
  vmin <- min(plotData$value, na.rm = TRUE)
  vmax <- max(plotData$value, na.rm = TRUE)

  ggplot(plotData, aes(x = Col, y = Row, fill = value)) +
    coord_fixed() +
    geom_tile() +
    guides(fill = guide_colourbar(title = NULL, barwidth = 0.5, barheight = 5)) +
    scale_y_discrete(limits = rev) +
    scale_fill_gradientn(colours = colors,
                         values = rescale(c(vmin, 0, vmax)),
                         limits = c(vmin, vmax),
                         na.value = "white") +
    geom_vline(xintercept = bnds, linetype = "dashed") +
    geom_hline(yintercept = y_bnds, linetype = "dashed") +
    labs(x = NULL, y = NULL,
         title = sprintf("p = %d, Sparsity = %s", ncol(mat), sparsity)) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.ticks = element_blank(),
          # plot.margin = margin(1, 1, 1, 1, "mm"),
          plot.title = element_text(hjust = .5))
}
