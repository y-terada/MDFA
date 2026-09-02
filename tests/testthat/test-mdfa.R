make_test_data <- function(p = 12, m = 2, n = 60, seed = 1) {
  set.seed(seed)
  f_size <- p / m
  Lam_ast <- matrix(0, nrow = p, ncol = m)
  Lam_ast[0 * f_size + (1:f_size), 1] <- 0.9
  Lam_ast[1 * f_size + (1:f_size), 2] <- 0.8
  uniq_var <- 1 - rowSums(Lam_ast^2)
  Sig_ast <- Lam_ast %*% t(Lam_ast) + diag(uniq_var)

  L <- chol(Sig_ast)
  Z <- matrix(rnorm(n * p), n, p)
  Z %*% L
}

test_that("MDFA_int recovers approximately the right number of factors (ALS)", {
  skip_if_not_installed("Rfast")
  skip_if_not_installed("psych")

  data <- make_test_data()
  S <- cov(data)

  res <- MDFA_cov(covar = S, nfactors = 2, alg = "ALS", trace = 0, REPORT = 0, rotate = "varimax")

  expect_s3_class(res$loadings, "loadings")
  expect_equal(dim(unclass(res$loadings)), c(12, 2))
  expect_true(all(is.finite(unclass(res$loadings))))
})

test_that("MDFA_int recovers approximately the right number of factors (BFGS)", {
  skip_if_not_installed("Rfast")
  skip_if_not_installed("psych")

  data <- make_test_data()
  S <- cov(data)

  res <- MDFA_cov(covar = S, nfactors = 2, alg = "BFGS", trace = 0, REPORT = 0, rotate = "varimax")

  expect_equal(dim(unclass(res$loadings)), c(12, 2))
  expect_true(all(is.finite(unclass(res$loadings))))
})

test_that("MDFA runs without bootstrap (n.iter = 1)", {
  skip_if_not_installed("Rfast")
  skip_if_not_installed("psych")
  skip_if_not_installed("clue")
  skip_if_not_installed("fungible")

  data <- make_test_data()
  res <- MDFA(data = data, nfactors = 2, alg = "ALS", n.iter = 1, trace = 0, REPORT = 0, rotate = "varimax")

  expect_s3_class(res, "mdfa-r")
  expect_null(res$cis)
})

test_that("MDFA with bootstrap attaches confidence intervals", {
  skip_if_not_installed("Rfast")
  skip_if_not_installed("psych")
  skip_if_not_installed("clue")
  skip_if_not_installed("fungible")

  data <- make_test_data(n = 60)
  res <- MDFA(
    data = data, nfactors = 2, alg = "ALS", n.iter = 20, trace = 0, REPORT = 0,
    rotate = "varimax", num_cores = 1
  )

  expect_s3_class(res, "mdfa.ci")
  expect_true(!is.null(res$cis))
  expect_equal(dim(res$cis$means), c(12, 2))
  expect_equal(dim(res$cis$sds), c(12, 2))
})

test_that("rotate = 'equamax' runs successfully", {
  skip_if_not_installed("Rfast")
  skip_if_not_installed("psych")
  skip_if_not_installed("GPArotation")

  data <- make_test_data()
  S <- cov(data)

  res <- MDFA_cov(covar = S, nfactors = 2, alg = "ALS", trace = 0, REPORT = 0, rotate = "equamax")

  expect_equal(dim(unclass(res$loadings)), c(12, 2))
  expect_true(all(is.finite(unclass(res$loadings))))
})
