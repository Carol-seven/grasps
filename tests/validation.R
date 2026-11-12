# Standalone regression tests for input validation and failed fits.
# The broader published-reference suite remains on test/numerical-validation.

library(grasps)

## test helpers ----------------------------------------------------------------

# Check that an expression raises an error containing the expected text.
expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(pattern, message, fixed = TRUE))
}

# Evaluate an expression while capturing and muffling its warning messages.
with_warnings <- function(expr) {
  messages <- character()
  result <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(result = result, messages = messages)
}

# Fit a small default model, replacing arguments supplied through `...`.
fit <- function(...) {
  args <- list(X = diag(2), n = 20, membership = 1:2,
               penalty = "lasso", lambda = 0.1, alpha = 0.5)
  changes <- list(...)
  args[names(changes)] <- changes
  do.call(grasps, args)
}

## reject malformed data before conversion to Armadillo matrices ---------------

expect_error(fit(X = data.frame(x = 1:3)), "`X`")
expect_error(fit(X = matrix(NA_real_, 2, 2)), "`X`")
expect_error(fit(X = matrix(numeric(), 0, 2)), "`X`")
expect_error(fit(X = diag(2) + 1i), "`X`")
expect_error(fit(X = diag(c(0, 1))), "positive diagonal")
expect_error(fit(X = matrix(c(1, 2, 2, 1), 2)), "positive semidefinite")
expect_error(fit(membership = c(1, NA)), "`membership`")
expect_error(fit(membership = c(1, 1.5)), "`membership`")
expect_error(fit(membership = c("a", "b")), "`membership`")
expect_error(fit(penalty = c("lasso", "scad")), "`penalty`")
expect_error(fit(crit = NA_character_), "`crit`")
expect_error(fit(diag.ind = 1), "`diag.ind`")

## validate scalar and active controls -----------------------------------------

bad.controls <- list(rho = Inf, tau.incr = NA_real_, tau.decr = 1,
                     nu = c(2, 3), tol.abs = NULL, tol.rel = "small",
                     maxiter = 0)
for (name in names(bad.controls)) {
  expect_error(do.call(fit, bad.controls[name]), paste0("`", name, "`"))
}
expect_error(fit(maxiter = 1.5), "`maxiter`")
expect_error(fit(maxiter = 2^31), "`maxiter`")
expect_error(fit(lambda = 1i), "`lambda`")
expect_error(fit(alpha = NaN), "`alpha`")
expect_error(fit(crit = "EBIC", ebic.tuning = 1.1), "`ebic.tuning`")
expect_error(fit(lambda = NULL, nlambda = "ten"), "`nlambda`")
expect_error(fit(lambda = NULL, tol.lambda = 0), "`tol.lambda`")
expect_error(fit(lambda = NULL, maxiter.lambda = 0), "`maxiter.lambda`")
expect_error(fit(lambda = NULL, growiter.lambda = -1), "`growiter.lambda`")
expect_error(fit(lambda = NULL, lambda.min.ratio = 0), "`lambda.min.ratio`")
bad.gamma <- list(adapt = 0, atan = Inf, exp = -1, lq = 1, lsp = 0,
                  mcp = 1, scad = 2)
for (penalty in names(bad.gamma)) {
  expect_error(fit(penalty = penalty, gamma = bad.gamma[[penalty]]), "`gamma`")
}
expect_error(fit(penalty = "scad", gamma = c(3, 4)), "`gamma`")

## validate initial estimators -------------------------------------------------

expect_error(fit(penalty = "scad", initial = "invalid"), "`initial`")
expect_error(fit(penalty = "scad", initial = diag(3)), "`initial`")
expect_error(fit(penalty = "scad",
                 initial = matrix(c(1, 0, 1, 1), 2)), "`initial`")
expect_error(fit(penalty = "scad",
                 initial = matrix(c(1, Inf, Inf, 1), 2)), "`initial`")
expect_error(fit(penalty = "scad", initial = diag(2) + 1i), "`initial`")
expect_error(fit(X = diag(2), n = 2, penalty = "scad", initial = "invS"),
             "requires n > p")
expect_error(fit(X = matrix(1, 2, 2), penalty = "scad", initial = "invS"),
             "not invertible")

lasso.initial <- fit(penalty = "scad", initial = "lasso")
inverse.initial <- fit(penalty = "scad", initial = "invS")
matrix.initial <- fit(penalty = "scad", initial = diag(2))
stopifnot(isTRUE(lasso.initial$initial.converged),
          is.na(inverse.initial$initial.converged),
          is.na(matrix.initial$initial.converged),
          isTRUE(all.equal(inverse.initial$initial, diag(2))),
          isTRUE(all.equal(matrix.initial$initial, diag(2))))

# `initial` is inactive for the lasso penalty and must therefore be ignored.
lasso.initial.inactive <- fit(initial = "invalid")
stopifnot(lasso.initial.inactive$converged)

inverse.grid <- fit(penalty = "scad", lambda = c(0.1, 0.2),
                    initial = "invS")
stopifnot(length(inverse.grid$converged.grid) == 2L,
          is.na(inverse.grid$initial.converged))

## validate iterative LLA controls only when they are active -------------------

expect_error(fit(penalty = "scad", lla = 1), "`lla`")
expect_error(fit(penalty = "scad", lla = NA), "`lla`")
expect_error(fit(penalty = "scad", lla = c(TRUE, FALSE)), "`lla`")
expect_error(fit(penalty = "scad", lla = TRUE, tol.lla = 0), "`tol.lla`")
expect_error(fit(penalty = "scad", lla = TRUE, tol.lla = Inf), "`tol.lla`")
expect_error(fit(penalty = "scad", lla = TRUE, maxiter.lla = 0), "`maxiter.lla`")
expect_error(fit(penalty = "scad", lla = TRUE, maxiter.lla = 1.5), "`maxiter.lla`")

## do not range-check LLA controls when iterative LLA is inactive --------------

lla.inactive <- fit(penalty = "scad", lla = FALSE, tol.lla = 0, maxiter.lla = 0)
stopifnot(lla.inactive$lla.iterations == 1L, is.na(lla.inactive$lla.converged))
lasso.lla <- fit(penalty = "lasso", lla = TRUE, tol.lla = 0, maxiter.lla = 0)
stopifnot(is.na(lasso.lla$lla.iterations), is.na(lasso.lla$lla.converged))

## adaptive-lasso weights remain fixed regardless of `lla` ---------------------

adapt.one.step <- fit(penalty = "adapt", lla = FALSE)
adapt.lla <- fit(penalty = "adapt", lla = TRUE, tol.lla = 0, maxiter.lla = 0)
stopifnot(isTRUE(all.equal(adapt.lla$hatOmega, adapt.one.step$hatOmega)),
          is.na(adapt.lla$lla.iterations), is.na(adapt.lla$lla.converged))

## validate sample-size and cross-validation controls --------------------------

expect_error(grasps(diag(2), membership = 1:2, penalty = "lasso",
                    lambda = 0.1, alpha = 0.5), "`n` must")
expect_error(fit(n = 1), "`n`")
X <- cbind(sin(1:12), cos((1:12) / 2))
set.seed(2026)
expect_error(fit(X = X, n = 11), "must equal `nrow(X)`")
expect_error(fit(X = cbind(1:12, rep(1, 12)), n = 12), "positive diagonal")
expect_error(fit(crit = "CV"), "CV requires")
expect_error(fit(X = X, n = 12, crit = "CV", kfold = 2.5), "`kfold`")
expect_error(fit(X = X, n = 12, crit = "CV", kfold = 12),
             "at least 2 observations")
expect_error(fit(X = X, n = 12, penalty = "scad", crit = "CV",
                 initial = diag(2)), "cannot be used")
expect_error(fit(X = X[1:4, ], n = 4, penalty = "scad", crit = "CV",
                 kfold = 2, initial = "invS"), "training observations")

set.seed(2026)
inverse.cv <- fit(X = X, n = 12, penalty = "scad", crit = "CV", kfold = 3,
                  lambda = c(0.1, 0.2), initial = "invS")
stopifnot(identical(dim(inverse.cv$converged.grid), c(2L, 3L)),
          is.na(inverse.cv$initial.converged))

## retain valid legacy inputs --------------------------------------------------

# Singular covariances, nonconsecutive group labels, and inactive controls used
# to work and must not be rejected by stricter validation.

valid <- with_warnings(fit(X = matrix(1, 2, 2), membership = c(-2L, 7L),
                           gamma = "ignored", nlambda = NA,
                           lambda.min.ratio = -1, tol.lambda = 0))
stopifnot(valid$result$converged, length(valid$messages) == 0L)
raw <- grasps(X, membership = 1:2, penalty = "lasso", lambda = 0.1, alpha = 0.5)
covariance <- fit(X = crossprod(scale(X, scale = FALSE)) / nrow(X), n = nrow(X))
stopifnot(isTRUE(all.equal(raw$hatOmega, covariance$hatOmega)))

## distinguish one-step from iterative LLA for scalar SCAD ---------------------

scad.one.step <- grasps(matrix(1), n = 20, membership = 1L, penalty = "scad",
                        lambda = 0.3, alpha = 1, lla = FALSE)
scad.iterative <- grasps(matrix(1), n = 20, membership = 1L, penalty = "scad",
                         lambda = 0.3, alpha = 1, lla = TRUE, maxiter.lla = 100)
stopifnot(isTRUE(all.equal(scad.one.step$hatOmega[1, 1], 0.8879,
                           tolerance = 1e-4)),
          scad.one.step$lla.iterations == 1L,
          is.na(scad.one.step$lla.converged),
          isTRUE(all.equal(scad.iterative$hatOmega[1, 1], 0.9411,
                           tolerance = 1e-4)),
          scad.iterative$lla.iterations > 1L,
          isTRUE(scad.iterative$lla.converged))

## report failure to satisfy the iterative LLA stopping criterion --------------

scad.limited <- with_warnings(grasps(
  matrix(1), n = 20, membership = 1L, penalty = "scad",
  lambda = 0.3, alpha = 1, lla = TRUE,
  tol.lla = 1e-15, maxiter.lla = 1
))
stopifnot(scad.limited$result$lla.iterations == 1L,
          isFALSE(scad.limited$result$lla.converged),
          length(scad.limited$messages) == 1L,
          grepl("iterative LLA convergence criterion",
                scad.limited$messages, fixed = TRUE),
          isTRUE(all.equal(scad.limited$result$hatOmega,
                           scad.one.step$hatOmega, tolerance = 1e-8)))

## retain iterative LLA convergence status for every IC candidate --------------

scad.grid <- with_warnings(fit(penalty = "scad",
                               lambda = c(0.1, 0.2), alpha = 1,
                               lla = TRUE, tol.lla = 1e-15,
                               maxiter.lla = 1))
stopifnot(length(scad.grid$result$lla.converged.grid) == 2L,
          all(scad.grid$result$lla.converged.grid == FALSE),
          isFALSE(scad.grid$result$lla.converged),
          any(grepl("iterative LLA convergence criterion",
                    scad.grid$messages, fixed = TRUE)))

## reject invalid model-selection candidates -----------------------------------

expect_error(fit(lambda = c(1e10, 2e10), maxiter = 1),
             "No candidate has a finite information-criterion score")
expect_error(fit(X = X, n = 12, crit = "CV", kfold = 3,
                  lambda = c(1e10, 2e10), maxiter = 1),
             "No candidate has a finite cross-validation score")
expect_error(fit(lambda = 1e10, maxiter = 1),
             "positive-definite precision matrix")
selected <- with_warnings(fit(lambda = c(1e10, 0.1), maxiter = 1))
stopifnot(selected$result$lambda == 0.1,
          is.infinite(selected$result$IC.score[1]),
          is.finite(selected$result$IC.score[2]),
          !selected$result$converged,
          length(selected$messages) == 1L,
          grepl("did not satisfy the ADMM convergence", selected$messages,
                fixed = TRUE))
single <- with_warnings(fit(maxiter = 1))
stopifnot(!single$result$converged, length(single$messages) == 1L)
set.seed(2026)
cv <- with_warnings(fit(X = X, n = 12, crit = "CV", kfold = 3,
                         lambda = c(1e10, 0.1), maxiter = 1))
stopifnot(cv$result$lambda == 0.1,
          all(is.infinite(cv$result$CV.loss[1, ])),
          all(is.finite(cv$result$CV.loss[2, ])),
          length(cv$messages) == 1L)

## CV records one iterative LLA status per candidate and fold ------------------

set.seed(2026)
scad.cv <- with_warnings(fit(X = X, n = 12, penalty = "scad",
                             crit = "CV", kfold = 3,
                             lambda = c(0.1, 0.2), alpha = 1,
                             lla = TRUE, tol.lla = 1e-15,
                             maxiter.lla = 1))
stopifnot(identical(dim(scad.cv$result$lla.converged.grid), c(2L, 3L)),
          all(scad.cv$result$lla.converged.grid == FALSE),
          any(grepl("iterative LLA convergence criterion",
                    scad.cv$messages, fixed = TRUE)))

## reject indefinite estimates with a positive determinant ---------------------

# A positive determinant does not establish positive definiteness.
# Construct a covariance whose first ADMM precision update has two identical
# SPD blocks.
# Thresholding each block by 0.15 gives eigenvalues of opposite signs,
# so the full thresholded matrix has two negative eigenvalues but a positive
# determinant.

first.update <- matrix(0, 4, 4)
first.update[1:2, 1:2] <- first.update[3:4, 3:4] <-
  matrix(c(0.1, 0.2, 0.2, 0.5), 2, 2)
S <- solve(first.update) - 2 * first.update
invalid <- pmax(first.update - 0.15, 0)
stopifnot(det(invalid) > 0,
          sum(eigen(invalid, symmetric = TRUE,
                    only.values = TRUE)$values < 0) == 2L)
indefinite <- with_warnings(fit(X = S, membership = 1:4, alpha = 1,
                                 lambda = c(0.3, 0.001), maxiter = 1))
stopifnot(is.infinite(indefinite$result$IC.score[1]),
          is.finite(indefinite$result$IC.score[2]),
          indefinite$result$lambda == 0.001,
          length(indefinite$messages) == 1L)
