
#' Matrix Decomposition Factor Analysis for a covariance matrix
#'
#' Performs Matrix Decomposition Factor Analysis (MDFA) for a given
#' covariance matrix. Unlike [MDFA()], this function does not perform
#' bootstrap replications and therefore does not provide standard errors
#' for the estimated loadings. Loadings and unique
#' variances are estimated from a data matrix (or its sample covariance
#' matrix) by minimizing a discrepancy between the data and the model,
#' using either an Alternating Least Squares (ALS) algorithm (see Adachi
#' and Trendafilov (2018)) or a quasi-Newton (BFGS)
#' algorithm. Here, we note that the ALS algorithm proposed by Adachi (2012) is employed.
#' Thus, the common and unique factor scores are not estimated.
#'
#' @param covar A covariance matrix.
#' @param nfactors Number of factors.
#' @param tol The absolute convergence tolerance.
#' @param maxit The maximum number of iterations.
#' @param SMC If `TRUE`, use squared multiple correlations as the initial
#'   communality estimate. If `FALSE`, initial unique variances are drawn
#'   randomly from `uniq_range`.
#' @param trace Non-negative integer. If positive, tracing information on
#'   the progress of the optimization is produced.
#' @param REPORT Reporting frequency: progress is printed every `REPORT`
#'   iterations when `trace > 0`.
#' @param rotate Rotation method applied to the loadings. One of `"none"`,
#'   `"varimax"`, `"Varimax"`, `"quartimax"`, `"bentlerT"`, `"targetT"`,
#'   `"TargetT"`, or `"equamax"`. Rotations other than `"varimax"` and
#'   `"none"` require the \pkg{GPArotation} package.
#' @param n.rotations If greater than 1, the rotation is instead performed
#'   `n.rotations` times from random orthogonal starting configurations via
#'   \code{\link[psych]{faRotations}} (borrowed from the \pkg{psych}
#'   package), and the best solution (by hyperplane count) is kept. This can
#'   help avoid poor local optima in the rotation step. Requires the
#'   \pkg{psych} package.
#' @param hyper Hyperplane count threshold passed to
#'   \code{\link[psych]{faRotations}} when `n.rotations > 1`; unused
#'   otherwise.
#' @param uniq_range Interval from which random initial unique variances are
#'   drawn when `SMC = FALSE`.
#' @param alg Estimation algorithm: `"ALS"` for Alternating Least Squares, or
#'   `"BFGS"` for a quasi-Newton gradient-based optimizer (via [stats::optim()]).
#' @param ... Additional arguments passed on to the rotation function.
#'
#' @return A list with elements:
#' \describe{
#'   \item{loadings}{A matrix of loadings, one column per factor.}
#'   \item{uniquenesses}{`diag(covar - model)`; note this differs from the
#'     unique variances estimated directly by MDFA, which are not returned
#'     by this function.}
#'   \item{communality, communalities}{Communalities implied by the model.}
#'   \item{e.values}{Eigenvalues of `covar`.}
#'   \item{model}{The model-implied covariance matrix, `loadings \%*\% t(loadings)`.}
#'   \item{rotation}{The rotation method used.}
#'   \item{rot.mat}{The rotation matrix, if a rotation was applied.}
#'   \item{fm}{Always `"mdfa"`; identifies the factoring method.}
#' }
#'
#' @references
#' Adachi, K. (2012). Some contributions to data-fitting factor analysis
#' with empirical comparisons to covariance-fitting factor analysis.
#' Journal of the Japanese Society of Computational Statistics, 25, 25-38.
#' \doi{10.5183/jjscs.1106001_197}
#'
#' Adachi, K. and Trendafilov, N. T. (2018). Some mathematical properties of
#' the matrix decomposition solution in factor analysis. Psychometrika, 83,
#' 407-424. \doi{10.1007/s11336-017-9600-y}
#'
#' @export
MDFA_cov <- function(covar, nfactors = 3, tol = 1e-5, maxit = 1000, SMC = TRUE, trace = 1, REPORT = 1,
                      rotate = "varimax", n.rotations = 1, hyper = 0.15, uniq_range = c(0.2, 1), alg = "ALS", ...) {
  p <- ncol(covar)
  S <- covar
  num_lam <- nfactors * p

  # initial value: smc & one step ols
  # ----------------------------------------------------------------
  if (SMC == TRUE) {
    ini_para <- smc_loading_init(S = S, nfactors = nfactors, tol = tol)
  } else {
    psi2 <- stats::runif(p, min = uniq_range[1], max = uniq_range[2])
    Sc <- S - diag(psi2)
    res_eigen <- Rfast::eigen.sym(Sc, k = nfactors, vectors = TRUE)
    d_tmp <- c(pmax(res_eigen$values, 0))
    Lam <- t(t(res_eigen$vectors) * sqrt(d_tmp))
    ini_para <- list(loadings = Lam, uniquenesses = psi2)
  }
  # ----------------------------------------------------------------

  Lam <- ini_para$loadings
  psi2 <- ini_para$uniquenesses

  if (alg[1] == "ALS") {
    loss <- Inf
    loss_vec <- rep(NA, maxit)
    cnt <- 0
    B <- cbind(Lam, diag(sqrt(psi2)))
    tmp_B <- B
    S_diag <- diag(S)

    for (t in 1:maxit) {
      cnt <- cnt + 1

      BtSB <- crossprod(B, S %*% B)
      eig_BtSB <- Rfast::eigen.sym(BtSB, k = p)
      V <- eig_BtSB$vectors
      a <- as.numeric(sqrt(pmax(eig_BtSB$values, 0)))
      BBt <- tcrossprod(B)
      Rchol <- chol(BBt)
      BV <- B %*% V
      A <- backsolve(Rchol, forwardsolve(t(Rchol), BV))
      S_XZ <- (t(t(A) * a)) %*% t(V)

      tmp_Lam <- S_XZ[, 1:nfactors]
      tmp_sig <- diag(S_XZ[, -(1:nfactors)])
      tmp_B[, 1:nfactors] <- tmp_Lam
      diag(tmp_B[, -(1:nfactors)]) <- tmp_sig
      tmp_cov <- tmp_Lam %*% t(tmp_Lam)
      tmp_diag <- diag(tmp_cov) + tmp_sig^2
      loss_vec[t] <- tmp_loss <- sum(S_diag - tmp_diag)
      if (loss - tmp_loss < tol * (abs(loss) + tol)) {
        loss <- tmp_loss
        B <- tmp_B
        break
      }

      loss <- tmp_loss
      B <- tmp_B
      if (trace > 0) {
        if (t %% REPORT == 0) {
          print(paste("t =", t, "loss =", tmp_loss))
        }
      }
    } # END: for t

    Lam <- B[, 1:nfactors]
    uniq <- diag(B[, -(1:nfactors)])^2
    fval <- loss
  } else if (alg[1] == "BFGS") {
    ini_par <- c(Lam, sqrt(psi2))
    L <- tryCatch(
      t(chol(S)),
      error = function(e) {
        eigS <- eigen((S + t(S)) / 2, symmetric = TRUE)
        eigS$vectors %*%
          (sqrt(pmax(eigS$values, 1e-12)) * t(eigS$vectors))
      }
    )
    fit_mdfa <- stats::optim(
      par = ini_par, fn = loss_mdfa, gr = grad_mdfa,
      Sig = S, num_factor = nfactors, L = L,
      method = "BFGS", control = list(trace = trace, maxit = maxit, reltol = tol, REPORT = max(1L, REPORT))
    )

    Lam[1:num_lam] <- fit_mdfa$par[1:num_lam]
    uniq <- fit_mdfa$par[-(1:num_lam)]^2
    fval <- fit_mdfa$value
  } else {
    stop('`alg` must be one of "ALS" or "BFGS".')
  }

  if (REPORT > 0) {
    print(paste("Final loss =", fval))
  }

  # fac like processing
  # ----------------------------------
  e.values <- eigen(covar, symmetric = TRUE)$values
  loadings <- Lam
  uniquenesses <- uniq
  model <- loadings %*% t(loadings)

  if (!is.double(loadings)) {
    warning("the matrix has produced imaginary results -- proceed with caution")
    loadings <- matrix(as.double(loadings), ncol = nfactors)
  }
  if (nfactors > 1) {
    sign.tot <- vector(mode = "numeric", length = nfactors)
    sign.tot <- sign(colSums(loadings))
    sign.tot[sign.tot == 0] <- 1
    loadings <- loadings %*% diag(sign.tot)
  } else {
    if (sum(loadings) < 0) {
      loadings <- -as.matrix(loadings)
    } else {
      loadings <- as.matrix(loadings)
    }
    colnames(loadings) <- "MR1"
  }

  colnames(loadings) <- paste("MDFA", 1:nfactors, sep = "")
  rownames(loadings) <- rownames(covar)

  model <- loadings %*% t(loadings)
  rot.mat <- NULL
  rotated <- NULL
  if (rotate != "none") {
    if (nfactors > 1) {
      if (n.rotations > 1) {
        if (!requireNamespace("psych", quietly = TRUE)) {
          stop("I am sorry, n.rotations > 1 requires the psych package to be installed")
        }
        rotated <- psych::faRotations(loadings, r = covar, n.rotations = n.rotations, rotate = rotate, hyper = hyper, ...)
        loadings <- rotated$loadings
        rot.mat <- rotated$rot.mat
      } else {
        rotated <- NULL
        if (rotate %in% c("varimax", "Varimax", "quartimax", "bentlerT", "targetT", "TargetT", "equamax")) {
          Phi <- NULL
          switch(rotate,
            varimax = {
              rotated <- stats::varimax(loadings)
              loadings <- rotated$loadings
              rot.mat <- rotated$rotmat
            },
            Varimax = {
              if (!requireNamespace("GPArotation", quietly = TRUE)) {
                stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
              }
              rotated <- GPArotation::Varimax(loadings, ...)
              loadings <- rotated$loadings
              rot.mat <- t(solve(rotated$Th))
            },
            quartimax = {
              if (!requireNamespace("GPArotation", quietly = TRUE)) {
                stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
              }
              rotated <- GPArotation::quartimax(loadings, ...)
              loadings <- rotated$loadings
              rot.mat <- t(solve(rotated$Th))
            },
            bentlerT = {
              if (!requireNamespace("GPArotation", quietly = TRUE)) {
                stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
              }
              rotated <- GPArotation::bentlerT(loadings, ...)
              loadings <- rotated$loadings
              rot.mat <- t(solve(rotated$Th))
            },
            targetT = {
              if (!requireNamespace("GPArotation", quietly = TRUE)) {
                stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
              }
              rotated <- GPArotation::targetT(loadings, Tmat = diag(ncol(loadings)), ...)
              loadings <- rotated$loadings
              rot.mat <- t(solve(rotated$Th))
            },
            TargetT = {
              if (!requireNamespace("GPArotation", quietly = TRUE)) {
                stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
              }
              rot <- GPArotation::targetT(loadings, Tmat = diag(ncol(loadings)), ...)
              loadings <- rot$loadings
              rot.mat <- t(solve(rot$Th))
            },
            equamax = {
              if (!requireNamespace("GPArotation", quietly = TRUE)) {
                stop("I am sorry, to do this rotation requires the GPArotation package to be installed")
              }
              rotated <- GPArotation::equamax(loadings, ...)
              loadings <- rotated$loadings
              rot.mat <- t(solve(rotated$Th))
            }
          )
        }
      }
    }
  } else {
    rotated <- NULL
  }
  # ----------------------------------

  signed <- sign(colSums(loadings))
  signed[signed == 0] <- 1
  loadings <- loadings %*% diag(signed)

  colnames(loadings) <- paste("MDFA", 1:nfactors, sep = "")
  rownames(loadings) <- rownames(covar)

  if (nfactors > 1) {
    ev.rotated <- diag(t(loadings) %*% loadings)
    ev.order <- order(ev.rotated, decreasing = TRUE)
    loadings <- loadings[, ev.order]
  }
  rownames(loadings) <- colnames(covar)
  class(loadings) <- "loadings"

  result <- list()
  result$rotation <- rotate
  result$communality <- diag(model)
  result$communalities <- diag(model)
  result$uniquenesses <- diag(covar - model) # Note: differs from the unique variances estimated by MDFA itself.
  result$e.values <- e.values
  result$loadings <- loadings
  result$model <- model
  result$fm <- "mdfa"
  result$rot.mat <- rot.mat

  return(result)
}

