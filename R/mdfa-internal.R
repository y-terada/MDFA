# Internal helper functions for MDFA
# These are not exported; they are used internally by MDFA_int() and MDFA().

#' Initial loadings via squared multiple correlations
#'
#' @param S A covariance matrix.
#' @param nfactors Number of factors.
#' @param tol Numerical tolerance (currently unused, kept for API stability).
#'
#' @return A list with `loadings` and `uniquenesses`.
#' @noRd
smc_loading_init <- function(S, nfactors, tol = 1e-12) {
  res_smc <- psych::smc(R = S, covar = TRUE)
  psi2 <- diag(S) - res_smc
  Sc <- S - diag(psi2)

  eg <- Rfast::eigen.sym(Sc, k = nfactors, vectors = TRUE)
  Lam <- sweep(eg$vectors, 2, sqrt(pmax(eg$values, 0)), "*")

  list(loadings = Lam, uniquenesses = psi2)
}

#' MDFA loss function
#'
#' @param par Parameter vector (loadings followed by sqrt-uniquenesses).
#' @param Sig Covariance matrix.
#' @param num_factor Number of factors.
#' @param L Optional lower-triangular Cholesky factor of `Sig` (`t(chol(Sig))`).
#'   Computed internally if not supplied.
#'
#' @return A single numeric loss value.
#' @noRd
loss_mdfa <- function(par, Sig, num_factor, L = NULL) {
  p <- nrow(Sig)
  m <- num_factor
  num_lam <- m * p

  lam_vec <- par[1:num_lam]
  psi_vec <- par[-(1:num_lam)]

  Lam <- matrix(lam_vec, p, m)

  SigPhi <- tcrossprod(Lam) + diag(psi_vec^2)

  if (is.null(L)) {
    L <- t(chol(Sig))
  }

  K <- crossprod(L, SigPhi %*% L)
  eg <- eigen((K + t(K)) / 2, symmetric = TRUE)

  sum(diag(Sig)) + sum(diag(SigPhi)) -
    2 * sum(sqrt(pmax(eg$values, 0)))
}

#' Inverse square root of a symmetric positive (semi-)definite matrix
#'
#' @param K A symmetric matrix.
#' @param eps Numerical floor applied to eigenvalues before inversion.
#'
#' @return The matrix \eqn{K^{-1/2}}.
#' @noRd
spd_inv_half <- function(K, eps = 1e-12) {
  K <- (K + t(K)) / 2
  eg <- eigen(K, symmetric = TRUE)
  Q <- eg$vectors
  d <- pmax(eg$values, eps)

  t(t(Q) / sqrt(d)) %*% t(Q)
}

#' Gradient of the MDFA loss function
#'
#' @inheritParams loss_mdfa
#' @param eps Numerical floor used inside [spd_inv_half()].
#'
#' @return A numeric gradient vector, in the same order as `par`.
#' @noRd
grad_mdfa <- function(par, Sig, num_factor, eps = 1e-12, L = NULL) {
  p <- nrow(Sig)
  m <- num_factor
  num_lam <- m * p

  lam_vec <- par[1:num_lam]
  psi_vec <- par[-(1:num_lam)]

  Lam <- matrix(lam_vec, p, m)

  if (is.null(L)) {
    L <- t(chol(Sig))
  }

  SigPhi <- tcrossprod(Lam)
  diag(SigPhi) <- diag(SigPhi) + psi_vec^2

  K <- crossprod(L, SigPhi %*% L)

  K_inv_half <- spd_inv_half(K, eps = eps)

  G <- L %*% K_inv_half %*% t(L)
  G <- (G + t(G)) / 2

  grad_Lam <- 2 * (Lam - G %*% Lam)
  grad_psi <- 2 * psi_vec * (1 - diag(G))

  c(
    grad_Lam,
    grad_psi
  )
}