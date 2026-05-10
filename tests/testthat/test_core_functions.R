context("emglasso core functions")
library(emglasso)

test_that("rmmnorm generates correct dimensions", {
  N <- 100
  D <- 5
  Tau <- c(0.5, 0.5)
  Mu <- list(rep(0, D), rep(1, D))
  Sigma <- list(diag(D), diag(D))
  
  X <- rmmnorm(N, D, Tau, Mu, Sigma)
  
  expect_equal(dim(X), c(N, D))
  expect_true(all(is.finite(X)))
})

test_that("rspdmatrix generates positive definite matrix", {
  D <- 5
  lambda <- 0.3
  
  M <- rspdmatrix(D, lambda)
  
  expect_equal(dim(M), c(D, D))
  expect_true(all(is.finite(M)))
  
  # Check symmetry
  expect_equal(M, t(M), tolerance = 1e-10)
  
  # Check positive definiteness by eigenvalues
  eigenvalues <- eigen(M)$values
  expect_true(all(eigenvalues > 0))
})

test_that("estep returns proper dimensions", {
  set.seed(42)
  N <- 50
  D <- 3
  K <- 2
  
  X <- matrix(rnorm(N * D), N, D)
  Tau <- c(0.5, 0.5)
  Mu <- list(rep(0, D), rep(1, D))
  Sigma <- list(diag(D), diag(D))
  Theta <- list(Tau = Tau, Mu = Mu, Sigma = Sigma)
  
  result <- estep(X, Theta, Sigma)
  
  expect_equal(dim(result$LogLike), c(N, K))
  expect_equal(dim(result$Tnk), c(N, K))
  
  # Check that responsibilities sum to 1
  row_sums <- rowSums(result$Tnk)
  expect_equal(row_sums, rep(1, N), tolerance = 1e-10)
  
  # Check that values are in [0, 1]
  expect_true(all(result$Tnk >= 0 & result$Tnk <= 1))
})

test_that("estep handles degenerate densities without NaN", {
  X <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  Theta <- list(
    Tau = c(0.5, 0.5),
    Mu = list(c(1000, 1000), c(-1000, -1000)),
    Sigma = list(diag(2), diag(2))
  )

  result <- estep(X, Theta, Theta$Sigma)

  expect_true(all(is.finite(result$Tnk)))
  expect_equal(rowSums(result$Tnk), rep(1, nrow(X)), tolerance = 1e-10)
})

test_that("mstep returns parameter list with correct structure", {
  set.seed(42)
  N <- 50
  K <- 2
  D <- 3
  
  # Create posterior responsibilities
  Tnk <- matrix(0, N, K)
  Tnk[1:(N/2), 1] <- 1
  Tnk[(N/2+1):N, 2] <- 1
  
  X <- matrix(rnorm(N * D), N, D)
  
  result <- mstep(X, Tnk)
  
  expect_true("Tau" %in% names(result))
  expect_true("Mu" %in% names(result))
  expect_true("Sigma" %in% names(result))
  
  expect_length(result$Tau, K)
  expect_length(result$Mu, K)
  expect_length(result$Sigma, K)
  
  # Check that tau sums to 1
  expect_equal(sum(result$Tau), 1, tolerance = 1e-10)
  
  # Check dimensions of means and covariances
  for (k in 1:K) {
    expect_equal(length(result$Mu[[k]]), D)
    expect_equal(dim(result$Sigma[[k]]), c(D, D))
  }
})

test_that("mstep validates malformed Tnk", {
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)

  expect_error(
    mstep(X, matrix(c(1, 2, 3), nrow = 3)),
    "nrow(Tnk) == nrow(x)",
    fixed = TRUE
  )
  expect_error(mstep(X, matrix(-1, nrow = 10, ncol = 2)), "non-negative")

  bad <- matrix(0, nrow = 10, ncol = 2)
  bad[, 1] <- 1 / 10
  expect_error(mstep(X, bad), "positive sum")
})

test_that("emglasso runs without error on small synthetic data", {
  set.seed(42)
  N <- 100
  D <- 3
  Tau0 <- c(0.4, 0.6)
  K0 <- length(Tau0)
  
  Mu0 <- list(rep(0, D), rep(2, D))
  Sigma0 <- list(diag(D), diag(D))
  
  X <- rmmnorm(N, D, Tau0, Mu0, Sigma0)
  
  # Run emglasso
  result <- emglasso(X, iTau = Tau0, Tol = 1e-2, MaxIt = 10)
  
  expect_true("Theta" %in% names(result))
  expect_true("BIC" %in% names(result))
  
  expect_true("Tau" %in% names(result$Theta))
  expect_true("Mu" %in% names(result$Theta))
  expect_true("Sigma" %in% names(result$Theta))
  
  expect_length(result$Theta$Tau, K0)
  expect_equal(sum(result$Theta$Tau), 1, tolerance = 1e-10)
})

test_that("emglasso handles higher dimension correctly", {
  set.seed(123)
  N <- 200
  D <- 10
  Tau0 <- c(0.5, 0.5)
  K0 <- length(Tau0)
  
  Mu0 <- list(rep(0, D), rep(3, D))
  Sigma0 <- list(diag(D), diag(D))
  
  X <- rmmnorm(N, D, Tau0, Mu0, Sigma0)
  
  result <- emglasso(X, iTau = Tau0, Tol = 1e-2, MaxIt = 10)
  
  expect_equal(length(result$Theta$Sigma), K0)
  for (k in 1:K0) {
    expect_equal(dim(result$Theta$Sigma[[k]]), c(D, D))
  }
})

test_that("rmmnorm normalizes Tau and returns exactly N samples", {
  set.seed(99)
  N <- 120
  D <- 4
  Tau <- c(2, 3, 5)
  Mu <- list(rep(0, D), rep(1, D), rep(2, D))
  Sigma <- list(diag(D), diag(D), diag(D))

  X <- rmmnorm(N, D, Tau, Mu, Sigma)

  expect_equal(nrow(X), N)
  expect_equal(ncol(X), D)
  expect_true(all(is.finite(X)))
})

test_that("rmmnorm rejects malformed inputs", {
  D <- 3
  Mu <- list(rep(0, D), rep(1, D))
  Sigma <- list(diag(D), diag(D))

  expect_error(rmmnorm(0, D, c(0.5, 0.5), Mu, Sigma), "N must")
  expect_error(rmmnorm(10, 0, c(0.5, 0.5), Mu, Sigma), "D must")
  expect_error(rmmnorm(10, D, c(-1, 2), Mu, Sigma), "Tau must")
  expect_error(rmmnorm(10, D, c(0, 0), Mu, Sigma), "positive sum")
  expect_error(rmmnorm(10, D, c(0.5, 0.5), list(rep(0, D)), Sigma), "Mu and Sigma")
})

test_that("rspdmatrix validates arguments", {
  expect_error(rspdmatrix(0, 0.5), "D must")
  expect_error(rspdmatrix(5, -0.1), "lambda must")
  expect_error(rspdmatrix(5, 1.1), "lambda must")
  expect_error(rspdmatrix(5, 0.3, epsilon = 0), "epsilon must")
})

test_that("emglasso normalizes iTau and validates inputs", {
  set.seed(123)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)

  fit <- emglasso(X, iTau = c(2, 3), Tol = 1e-2, MaxIt = 10)
  expect_equal(sum(fit$Theta$Tau), 1, tolerance = 1e-10)

  expect_error(emglasso(X, iTau = c(0, 0)), "positive sum")
  expect_error(emglasso(X, iTau = c(-1, 2)), "non-negative")
  expect_error(emglasso(X, iTau = c(0.5, 0.5), Tol = -1), "Tol must")
  expect_error(emglasso(X, iTau = c(0.5, 0.5), MaxIt = 0), "MaxIt must")
})
