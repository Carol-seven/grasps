library(grasps)
library(ggplot2)

pen_df <- compute_penalty(
  omega = seq(-4, 4, by = 0.01),
  penalty = c("lasso", "atan", "capped", "exp",
              "lq", "lsp", "mcp", "scad",
              "selo", "sica"),
  lambda = 1
)

plot(pen_df, xlim = c(-1, 1), ylim = c(0, 1), zoom.size = 1) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))
