
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication: Lee and Ng (2022, ICML)

Replication files for all numerical results in Lee and Ng (2022, ICML).

## Installation

Before running the replication code, it is necessary to install the R
package **sketching**.

You can install the released version of sketching from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sketching")
```

Alternatively, you can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools") if devtools is not installed
devtools::install_github("https://github.com/sokbae/sketching")
```

## Replication Code

The latest version of the paper is included in this Github repository
(see paper.pdf).

-   Table 1 can be replicated by running monte_carlo/monte_carlo_ols.R.

-   Table 2 can be replicated by running
    monte_carlo/monte_carlo_2sls_F\_test.R.

-   Table 3 can be replicated by running monte_carlo/monte_carlo_2sls.R.

-   Tables 4 and 5 can be replicated by running
    empirical_example/empirical_example.R.

## Reference

-   Lee, S. and Ng, S. (2022). “Least Squares Estimation Using Sketched
    Data with Heteroskedastic Errors,”
    [arXiv:2007.07781](https://arxiv.org/abs/2007.07781), accepted for
    presentation at the Thirty-ninth International Conference on Machine
    Learning (ICML 2022).
