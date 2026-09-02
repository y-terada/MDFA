# MDFA

Matrix Decomposition Factor Analysis (MDFA) for R.

## Model

In factor analysis, we consider the following model for a $p$-dimensional
observation $x$:

$$
x = \mu + \Lambda f + \epsilon,
$$

where $\mu \in \mathbb{R}^p$ is a mean vector, $m$ is the number of
factors ($m < p$), $\Lambda \in \mathbb{R}^{p\times m}$ is a factor
loading matrix, $f$ is an $m$-dimensional centered random vector with
the identity covariance, and $\epsilon$ is a $p$-dimensional
uncorrelated centered random vector, independent of $f$, with diagonal
covariance matrix $\mathrm{Var}(\epsilon) = \Psi^2 =
\mathrm{diag}(\sigma_1^2,\dots,\sigma_p^2)$. Each component of $f$ and
$\epsilon$ is called a common factor and a unique factor, respectively.

MDFA is a matrix factorization-based approach to factor analysis,
introduced in the early 2000s and actively developed in computational
statistics. In MDFA, the estimator is obtained by minimizing the
following principal-component-analysis-like loss function:

$$
\mathcal{L}_n(\Phi, Z) = \frac{1}{n}\|X_n - (F\Lambda^\top + E\Psi)\|_F^2
= \frac{1}{n}\sum_{i=1}^n \|x_i - (\Lambda f_i + \Psi e_i)\|^2,
$$

where $\Phi = [\Lambda, \Psi]$ collects the loading matrix $\Lambda$ and
the (diagonal) unique-variance matrix $\Psi$, $X_n$ is the (centered)
data matrix with rows $x_i$, and $Z = [F, E]$ collects the common factor
scores $F \in \mathbb{R}^{n\times m}$ and normalized unique factor scores
$E \in \mathbb{R}^{n\times p}$, constrained by

$$
\mathbf{1}_n^\top F = \mathbf{0}_m^\top,\quad
\mathbf{1}_n^\top E = \mathbf{0}_p^\top,\quad
\frac{1}{n}F^\top F = I_m,\quad
\frac{1}{n}E^\top E = I_p,\quad \text{and}\quad
F^\top E = O_{m\times p}.
$$

The minimization problem for $\mathcal{L}_n$ can be solved by a simple
alternating minimization algorithm. For estimating only $\Phi = [\Lambda,
\Psi]$, Adachi (2012) shows that the algorithm can be performed using
only the sample covariance matrix $\widehat{S}_n$. 
This package uses that algorithm and does not estimate $Z$.

Terada (2025) shows that the MDFA estimator can be formulated as a
minimum discrepancy function (MDF) estimator, using the squared
Bures-Wasserstein distance between the sample and modeled covariance
matrices. More precisely, 

$$
\begin{aligned}
\mathcal{L}_n(\Phi) &= \min_{Z} \mathcal{L}_n(\Phi, Z)
= \operatorname{tr}(\widehat{S}_n) + \operatorname{tr}(\Phi\Phi^\top)
- 2\operatorname{tr}\left\{(\Phi^\top \widehat{S}_n \Phi)^{1/2}\right\} \\
&= \operatorname{tr}(\widehat{S}_n) + \operatorname{tr}\{\Sigma(\Phi)\}
- 2\operatorname{tr}\left\{\left(\widehat{S}_n^{1/2}\Sigma(\Phi)
\widehat{S}_n^{1/2}\right)^{1/2}\right\}
= d_{\mathrm{BW}}^2\left(\widehat{S}_n, \Sigma(\Phi)\right),
\end{aligned}
$$
where $\Sigma(\Phi) := \Phi\Phi^\top = \Lambda\Lambda^\top + \Psi^2$ and $d_{\mathrm{BW}}(A, B)$ is the
Bures-Wasserstein distance between positive semidefinite matrices $A$
and $B$; see Bhatia et al. (2019).

This package estimates factor loadings and unique variances for a data
matrix (or its sample covariance matrix), using either the Alternating
Least Squares (ALS) algorithm above (Adachi, 2012; Adachi and
Trendafilov, 2018) or a quasi-Newton (BFGS) algorithm. Non-parametric
bootstrap replication of the estimation procedure can optionally be used
to obtain standard errors and confidence intervals for the loadings;
this bootstrap inference is justified by the consistency and asymptotic
normality of the MDFA estimator established in Terada (2025).

## Installation

This package is not yet on CRAN. Install the development version from
GitHub:

```r
# install.packages("devtools")
devtools::install_github("y-terada/MDFA")
```

Or, from a local checkout:

```r
devtools::install_local("path/to/MDFA")
```

### Dependencies

`MDFA` imports `Rfast`, `psych`, `clue`, and `parallel`. Install
them first if you don't already have them:

```r
install.packages(c("Rfast", "psych", "clue"))
```

`GPArotation` and `mvtnorm` are used for some rotation methods and in the
examples/vignette; install them too if you plan to use those:

```r
install.packages(c("GPArotation", "mvtnorm"))
```

## Usage

```r
library(MDFA)

set.seed(1)
p <- 20; m <- 3; n <- 100

f_size <- p / m
Lam_ast <- matrix(0, nrow = p, ncol = m)
Lam_ast[0 * f_size + (1:f_size), 1] <- 0.9
Lam_ast[1 * f_size + (1:f_size), 2] <- 0.8
Lam_ast[2 * f_size + (1:f_size), 3] <- 0.6

uniq_var <- 1 - rowSums(Lam_ast^2)
Sig_ast <- Lam_ast %*% t(Lam_ast) + diag(uniq_var)

data <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sig_ast)

# Point estimate only
res <- MDFA(data = data, nfactors = m, alg = "ALS", n.iter = 1, rotate = "varimax", trace = 0, REPORT = 0)
res$loadings

# With bootstrap confidence intervals
res_ci <- MDFA(
  data = data, nfactors = m, alg = "ALS", n.iter = 10^3,
  rotate = "varimax", trace = 0, REPORT = 0, num_cores = 1
)
res_ci$cis$means
res_ci$cis$sds
```

See `vignette("intro", package = "MDFA")` for a longer walkthrough, and
`?MDFA` / `?MDFA_cov` for full argument documentation. `MDFA()` fits from raw
data and supports bootstrap inference; `MDFA_cov()` fits directly from a
covariance (or correlation) matrix without bootstrapping.

## Known limitations

- Bootstrap replication uses `parallel::mclapply()`, which forks processes
  and therefore runs serially (with a warning) on Windows.

## References

* Bhatia, R., Jain, T. and Lim, Y. (2019). On the Bures-Wasserstein
  distance between positive definite matrices. *Expositiones
  Mathematicae*, 37(2), 165-191.
  [doi:10.1016/j.exmath.2018.01.002](https://doi.org/10.1016/j.exmath.2018.01.002)
* Terada, Y. (2025). Statistical properties of matrix decomposition
  factor analysis. [doi:10.48550/arXiv.2403.06968](https://doi.org/10.48550/arXiv.2403.06968)
* Adachi, K. (2012). Some contributions to data-fitting factor analysis
  with empirical comparisons to covariance-fitting factor analysis.
  *Journal of the Japanese Society of Computational Statistics*, 25,
  25-38. [doi:10.5183/jjscs.1106001_197](https://doi.org/10.5183/jjscs.1106001_197)
* Adachi, K. and Trendafilov, N. T. (2018). Some mathematical properties
  of the matrix decomposition solution in factor analysis.
  *Psychometrika*, 83, 407-424. [doi:10.1007/s11336-017-9600-y](https://doi.org/10.1007/s11336-017-9600-y)

## License

MIT © Yoshikazu Terada