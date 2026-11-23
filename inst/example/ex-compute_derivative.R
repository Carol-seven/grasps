library(grasps)
library(ggplot2)

deriv_df <- compute_derivative(
  omega = seq(0, 4, by = 0.01),
  penalty = c("lasso", "atan", "capped", "exp",
              "lq", "lsp", "mcp", "scad",
              "selo", "sica"),
  lambda = 1
)

plot(deriv_df) +
  scale_y_continuous(limits = c(0, 1.5)) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))
