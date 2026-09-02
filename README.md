# MDFA

Matrix Decomposition Factor Analysis (MDFA) for R.

MDFA estimates factor loadings and uniquenesses for a covariance matrix by
minimizing a matrix-decomposition based loss function, using either an
Alternating Least Squares (ALS) algorithm or a quasi-Newton (BFGS) algorithm.
Non-parametric bootstrap replications can be used to obtain standard errors
and confidence intervals for the loadings.

## Installation

This package is not yet on CRAN. Install the development version from a
local source directory or from GitHub once published:

```r
# From a local checkout:
devtools::install_local("path/to/MDFA")

# Or, after this becomes a GitHub repo:
# devtools::install_github("yourusername/MDFA")
```

### Dependencies

`MDFA` imports `Rfast`, `psych`, `clue`, `fungible`, and `parallel`. Install
them first if you don't already have them:

```r
install.packages(c("Rfast", "psych", "clue", "fungible"))
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
res_ci <- MDFA(data = data, nfactors = m, alg = "ALS", n.iter = 10^3, rotate = "varimax", trace = 0, REPORT = 0)
res_ci$cis$means
res_ci$cis$sds
```

See `vignette("intro", package = "MDFA")` for a longer walkthrough, and
`?MDFA` for full argument documentation. (`MDFA()` is the only exported
function; the covariance-matrix-only fitting routine it calls internally,
`MDFA_int()`, is not part of the public interface, but is reachable via
`MDFA:::MDFA_int()` if needed.)

## Known limitations

- Bootstrap replication uses `parallel::mclapply()`, which forks processes
  and therefore runs serially (with a warning) on Windows.

## License

MIT © Your Name
