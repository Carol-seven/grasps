#' Plot Function for S3 Class "penderiv"
#'
#' @description
#' Generate a visualization of penalty functions produced by
#' \code{\link[grasps]{compute_penalty}}, or penalty derivatives produced by
#' \code{\link[grasps]{compute_derivative}}.
#' The plot automatically summarizes multiple configurations of penalty type,
#' \eqn{\lambda}, and \eqn{\gamma}. Optional zooming is supported through
#' \code{\link[ggforce]{facet_zoom}}.
#'
#' @param x An object inheriting from S3 class \code{"penderiv"}, typically
#' returned by \code{\link[grasps]{compute_penalty}}, or
#' \code{\link[grasps]{compute_derivative}}.
#'
#' @param ... Optional arguments passed to \code{\link[ggforce]{facet_zoom}}
#' to zoom in on a subset of the data, while keeping the view of the full
#' dataset as a separate panel.
#'
#' @return
#' An object of class \code{ggplot}.
#'
#' @example
#' inst/example/ex-plot.penderiv.R
#'
#' @import ggplot2
#' @import ggforce
#'
#' @export

plot.penderiv <- function(x, ...) {

  ## identify which configuration columns vary across rows
  target <- c("penalty", "lambda", "gamma")
  cols_keep <- target[sapply(x[target], function(col) {
    col_na_rm <- col[!is.na(col)]
    length(col_na_rm) && (anyNA(col) || any(col_na_rm != col_na_rm[1]))
  })]

  ## clean data
  df <- x[, (names(x) %in% c("omega", "value", cols_keep)), drop = FALSE]
  key <- paste(cols_keep, collapse = "_")

  ## build facet_zoom() layer if zooming arguments are provided
  fz <- if (length(list(...)) > 0L) {
    do.call(facet_zoom, list(...))
  } else {
    NULL
  }

  ## declare
  group <- omega <- value <- NULL

  ## all configurations identical (only one curve)
  if (key == "") {
    ggplot(df, aes(x = omega, y = value)) +
      geom_line() +
      fz +
      labs(x = expression(italic(omega)),
           y = if (inherits(x, "penalty")) {
             expression("Penalty Function" ~ italic(lambda) * italic(p) ~ "(" * italic(omega) * ")")
           } else if (inherits(x, "derivative")) {
             expression("Derivative Function" ~ italic(lambda) * italic(p) ~ "'(" * italic(omega) * ")")
           } else {
             NULL
           }) +
      theme_bw() +
      theme(legend.position = "bottom")

  ## multiple configurations, use color grouping and legend
  } else {

    ## nicely formatted labels for lambda and gamma
    lambda_fmt <- paste0("\u03BB = ", df$lambda)
    gamma_fmt <- ifelse(is.na(df$gamma), "", paste0("\u03B3 = ", df$gamma))

    ## group labels based on which parameter(s) vary
    df$group <- switch(
      key,
      "penalty" = df$penalty,
      "lambda" = lambda_fmt,
      "gamma" = gamma_fmt,
      "penalty_lambda" = paste0(df$penalty, " (", lambda_fmt, ")"),
      "penalty_gamma" = paste0(df$penalty, " (", gamma_fmt, ")"),
      "lambda_gamma" = paste0(lambda_fmt, ", ", gamma_fmt),
      "penalty_lambda_gamma" = paste0(df$penalty, " (", lambda_fmt, ", ", gamma_fmt, ")")
    )

    ## legend label based on which parameter(s) vary
    legend_name <- switch(
      key,
      lambda = "\u03BB",
      gamma  = "\u03B3",
      "Penalty Type"
    )

    ggplot(df, aes(x = omega, y = value, color = group)) +
      geom_line() +
      fz +
      labs(x = expression(italic(omega)),
           y = if (inherits(x, "penalty")) {
             expression("Penalty Function" ~ italic(lambda) * italic(p) ~ "(" * italic(omega) * ")")
           } else if (inherits(x, "derivative")) {
             expression("Derivative Function" ~ italic(lambda) * italic(p) ~ "'(" * italic(omega) * ")")
           } else {
             NULL
           },
           color = legend_name) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
}
