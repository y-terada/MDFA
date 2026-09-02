# MDFA 0.1.0

* Initial CRAN release.
* `MDFA()`: fits Matrix Decomposition Factor Analysis (MDFA) to a data
  matrix using either an Alternating Least Squares (ALS) algorithm
  (Adachi, 2012; Adachi and Trendafilov, 2018) or a quasi-Newton (BFGS)
  algorithm, with optional non-parametric bootstrap replication for
  standard errors and confidence intervals on the loadings.
* `MDFA_cov()`: fits MDFA directly from a covariance (or correlation)
  matrix, without bootstrapping.