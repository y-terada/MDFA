#' Matrix Decomposition Factor Analysis
#'
#' Performs exploratory factor analysis using Matrix Decomposition Factor
#' Analysis (MDFA), a matrix factorization-based estimator for the classical
#' factor model \eqn{x = \Lambda f + \Psi \epsilon}. Loadings and unique
#' variances are estimated from a data matrix by minimizing a discrepancy between the data and the model,
#' using either an Alternating Least Squares (ALS) algorithm (see Adachi
#' and Trendafilov (2018)) or a quasi-Newton (BFGS)
#' algorithm. Here, we note that the ALS algorithm proposed by Adachi (2012) is employed.
#' Thus, the common and unique factor scores are not estimated in MDFA. 
#' Non-parametric bootstrap replication of the estimation procedure can
#' optionally be used to obtain standard errors and confidence intervals
#' for the estimated loadings. This bootstrap inference is justified by the
#' consistency and asymptotic normality of the MDFA estimator established
#' in Terada (2025).
#'
#' @details
#' Matrix Decomposition Factor Analysis (MDFA) is a matrix
#' factorization-based approach to factor analysis introduced in the early
#' 2000s. Terada (2025) shows that the MDFA estimator is a minimum
#' discrepancy function estimator based on the squared
#' Bures-Wasserstein distance between the sample and modeled covariance
#' matrices.
#'
#' This package uses the unbiased sample covariance matrix (i.e.,
#' \code{\link[stats]{cov}}, with divisor `n - 1`) throughout.
#'
#' Bootstrap replication uses \code{\link[parallel]{mclapply}}, which relies
#' on forking to duplicate the current R process across `num_cores` workers.
#' This is fast, simple, and requires no extra setup on Unix-like systems
#' (Linux, macOS). Windows does not support forking, so on Windows
#' \code{mclapply()} automatically falls back to running serially (i.e.,
#' \code{num_cores} is effectively ignored); results are unaffected, but the
#' computation will be slower than on Linux/macOS with the same
#' \code{num_cores}.
#'
#'
#' @param data A numeric data matrix or data frame (rows = observations,
#'   columns = variables).
#' @param nfactors Number of factors.
#' @param tol The absolute convergence tolerance.
#' @param n.iter Number of bootstrap replications. If `n.iter <= 1`, no
#'   bootstrap is performed and only the estimate on the original data is
#'   returned.
#' @param maxit The maximum number of iterations for the ALS/BFGS algorithm.
#' @param SMC If `TRUE`, use squared multiple correlations as the initial
#'   communality estimate.
#' @param trace Non-negative integer controlling optimizer verbosity.
#' @param REPORT Reporting frequency; used when `trace > 0`.
#' @param uniq_range Interval from which random initial unique variances are
#'   drawn when `SMC = FALSE`.
#' @param rotate Rotation method: one of `"none"`, `"varimax"`, `"Varimax"`,
#'   `"quartimax"`, `"bentlerT"`, `"targetT"`, `"TargetT"`, or `"equamax"`.
#'   Rotations other than `"varimax"` and `"none"` require the
#'   \pkg{GPArotation} package.
#' @param n.rotations If greater than 1, the rotation is instead performed
#'   `n.rotations` times from random orthogonal starting configurations via
#'   \code{\link[psych]{faRotations}} (borrowed from the \pkg{psych}
#'   package), and the best solution (by hyperplane count) is kept. This can
#'   help avoid poor local optima in the rotation step. Requires the
#'   \pkg{psych} package.
#' @param hyper Hyperplane count threshold passed to
#'   \code{\link[psych]{faRotations}} when `n.rotations > 1`; unused
#'   otherwise.
#' @param alg Estimation algorithm, `"BFGS"` or `"ALS"`.
#' @param num_cores Number of cores used for bootstrap replications via
#'   [parallel::mclapply()]. Defaults to half of the cores detected by
#'   [parallel::detectCores()]. Note that `mclapply()` forks processes and
#'   therefore runs serially (with a warning) on Windows.
#' @param p Significance level used for the bootstrap confidence intervals
#'   (e.g. `p = 0.05` gives 95% confidence intervals).
#' @param ... Additional arguments passed on to the internal fitting
#'   routine and, from there, to the rotation function.
#'
#' @return An object of class `c("mdfa-r", "mdfa")` (no bootstrap) or
#'   `c("mdfa-r", "mdfa.ci")` (with bootstrap). It is a list with elements:
#' \describe{
#'   \item{loadings}{A matrix of loadings (class `"loadings"`), one column
#'     per factor, fitted on the original data.}
#'   \item{uniquenesses}{`diag(covar - model)` for the original data; note
#'     this differs from the unique variances estimated directly by MDFA,
#'     which are not returned by this function.}
#'   \item{communality, communalities}{Communalities implied by the model
#'     (identical values under both names).}
#'   \item{e.values}{Eigenvalues of `cov(data)`.}
#'   \item{model}{The model-implied covariance matrix, `loadings \%*\% t(loadings)`.}
#'   \item{rotation}{The rotation method used (the `rotate` argument).}
#'   \item{rot.mat}{The rotation matrix, or `NULL` if `rotate = "none"`.}
#'   \item{fm}{Always `"mdfa"`; identifies the factoring method.}
#'   \item{cis}{Only present when `n.iter > 1`. A list with:
#'   \describe{
#'     \item{cis$means}{Bootstrap mean loadings.}
#'     \item{cis$sds}{Bootstrap standard deviations of the loadings.}
#'     \item{cis$ci}{Bootstrap confidence intervals (`lower`, `upper`) for the loadings.}
#'     \item{cis$p}{Two-sided bootstrap p-values for the loadings.}
#'     \item{cis$replicates}{The raw bootstrap replicate loadings.}
#'   }}
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
#' Terada, Y. (2025). Statistical properties of matrix decomposition factor
#' analysis. \doi{10.48550/arXiv.2403.06968}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' p <- 12
#' m <- 2
#' n <- 80
#'
#' f_size <- p / m
#' Lam_ast <- matrix(0, nrow = p, ncol = m)
#' Lam_ast[0 * f_size + (1:f_size), 1] <- 0.9
#' Lam_ast[1 * f_size + (1:f_size), 2] <- 0.8
#' uniq_var <- 1 - rowSums(Lam_ast^2)
#' Sig_ast <- Lam_ast %*% t(Lam_ast) + diag(uniq_var)
#'
#' data <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sig_ast)
#' colnames(data) <- paste0("V", seq_len(p))
#'
#' res <- MDFA(
#'   data = data, nfactors = m, SMC = TRUE, alg = "ALS",
#'   n.iter = 10^3, trace = 0, REPORT = 0, num_cores = 1
#' )
#' res$loadings
#'
#' # Visualize the loadings and their bootstrap standard errors side by side
#' plot_heatmap <- function(M, pal, zlim, main) {
#'   nv <- nrow(M); nf <- ncol(M)
#'   image(
#'     x = seq_len(nf), y = seq_len(nv), z = t(M[nv:1, ]),
#'     zlim = zlim, axes = FALSE, xlab = "", ylab = "", col = pal
#'   )
#'   axis(1, at = seq_len(nf), labels = colnames(M), tick = FALSE)
#'   axis(2, at = seq_len(nv), labels = rev(rownames(M)), las = 2, tick = FALSE)
#'   title(main)
#' }
#' plot_colorbar <- function(pal, zlim) {
#'   legend_seq <- seq(zlim[1], zlim[2], length.out = 100)
#'   image(
#'     x = 1, y = legend_seq, z = matrix(legend_seq, nrow = 1),
#'     col = pal, axes = FALSE, xlab = "", ylab = ""
#'   )
#'   axis(4, las = 1)
#' }
#'
#' Lam <- unclass(res$loadings)
#' S <- res$cis$sds
#'
#' pal_lam <- colorRampPalette(c("blue", "white", "red"))(30)
#' maxabs <- max(abs(Lam))
#' zlim_lam <- c(-maxabs, maxabs) # symmetric around 0 (blue-white-red)
#'
#' pal_sd <- colorRampPalette(c("white", "red"))(30)
#' zlim_sd <- c(0, max(S) * 1.1) # a little headroom above the max
#'
#' layout(matrix(1:4, nrow = 1), widths = c(3, 1, 3, 1))
#'
#' par(mar = c(4, 5, 3, 1))
#' plot_heatmap(Lam, pal_lam, zlim_lam, "Loadings")
#' par(mar = c(4, 1, 3, 3))
#' plot_colorbar(pal_lam, zlim_lam)
#'
#' par(mar = c(4, 3, 3, 1))
#' plot_heatmap(S, pal_sd, zlim_sd, "Bootstrap SE")
#' par(mar = c(4, 1, 3, 3))
#' plot_colorbar(pal_sd, zlim_sd)
#'
#' layout(1)
#' }
#'
#'
#'
#' @export
MDFA <- function(data, nfactors = 3, tol = 1e-5, n.iter = 10^3, maxit = 1000, SMC = TRUE, trace = 1, REPORT = 1,
                  uniq_range = c(0.2, 1), rotate = "quartimax", n.rotations = 1, hyper = 0.15, alg = c("BFGS", "ALS"),
                  num_cores = parallel::detectCores() / 2, p = 0.05, ...) {
  nvar <- ncol(data)
  n.obs <- nrow(data)

  # MDFA for original data
  # ------------------------------------------------------------------------
  ecov <- stats::cov(data)
  res_org <- MDFA_cov(
    covar = ecov, nfactors = nfactors, tol = tol, maxit = maxit, SMC = SMC,
    trace = trace, REPORT = REPORT, rotate = rotate, n.rotations = n.rotations, hyper = hyper, uniq_range = uniq_range, alg = alg, ... = ...
  )
  fl <- res_org$loadings
  # ------------------------------------------------------------------------

  if (n.iter > 1) {
    replicateslist <- parallel::mclapply(
      1:n.iter,
      function(x) {
        data_bp <- data[sample(n.obs, n.obs, replace = TRUE), ]
        ecov_bp <- stats::cov(data_bp)
        fs <- MDFA_cov(
          covar = ecov_bp, nfactors = nfactors, tol = tol, maxit = maxit, SMC = SMC,
          trace = 0, REPORT = 0, uniq_range = uniq_range, alg = alg, rotate = rotate, n.rotations = n.rotations, hyper = hyper
        )
        # permutation and sign alignment
        # ------------------------------------
        Lbp <- as.matrix(fs$loadings)
        Lref <- as.matrix(fl)

        C <- crossprod(Lbp, Lref)

        perm <- max.col(abs(C))

        if (length(unique(perm)) != nfactors) {
          perm <- as.integer(clue::solve_LSAP(-abs(C)))
        }

        Lbp <- Lbp[, perm, drop = FALSE]

        sgn <- sign(diag(C[, perm, drop = FALSE]))
        sgn[sgn == 0] <- 1

        Lbp <- t(t(Lbp) * sgn)
        # ------------------------------------

        t.rot <- psych::target.rot(Lbp, Lref)
        list(loadings = t.rot$loadings)
      }, # END: function
      mc.cores = num_cores
    )

    replicates <- matrix(unlist(replicateslist), nrow = n.iter, byrow = TRUE)
    means <- colMeans(replicates, na.rm = TRUE)
    sds <- apply(replicates, 2, stats::sd, na.rm = TRUE)
    means <- matrix(means[1:(nvar * nfactors)], ncol = nfactors)
    sds <- matrix(sds[1:(nvar * nfactors)], ncol = nfactors)
    tci <- abs(means) / sds
    ptci <- 1 - stats::pnorm(tci)
    ci.lower <- means + stats::qnorm(p / 2) * sds
    ci.upper <- means + stats::qnorm(1 - p / 2) * sds
    ci <- data.frame(lower = ci.lower, upper = ci.upper)
    class(means) <- "loadings"
    colnames(means) <- colnames(sds) <- colnames(fl)
    rownames(means) <- rownames(sds) <- rownames(fl)
    res_org$cis <- list(means = means, sds = sds, ci = ci, p = 2 * ptci, replicates = replicates)
    results <- res_org
    class(results) <- c("mdfa-r", "mdfa.ci")
  } else {
    results <- res_org
    class(results) <- c("mdfa-r", "mdfa")
  }

  return(results)
}
