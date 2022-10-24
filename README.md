
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dunnR - Dunn test for multiple comparisons using rank sums

<!-- badges: start -->

<!-- badges: end -->

Non-parametric testing procedure for identifying differences between
multiple groups in a population. Implementation of Dunn JO (1964).
Multiple Comparisons Using Rank Sums, Technometrics, 6:3, 241-252. doi:
[10.1080/00401706.1964.10490181](https://doi.org/10.1080/00401706.1964.10490181)

## Installation

dunnR is not yet available on CRAN, but you can install it using
devtools:

``` r
install.packages("devtools")
devtools::install_github("bnanes/dunnR")
```

## Example

``` r
library(dunnR)
dunn.test(mpg ~ cyl, data = mtcars)
#> 
#>  Dunn test for multiple comparisons using rank sums
#> 
#> data:  mpg by cyl
#> Kruskal-Wallis chi-squared = 25.746, df = 2, p-value = 2.566e-06
#> 
#> -----------------------------------------
#> Post-hoc comparisons between groups
#> P value adjustment method: classic Dunn test (Bonferroni)
#> -----------------------------------------
#>   Group.A Group.B      Z p.unadjust   p.value Sig.
#> 1       4       6 2.1016    0.03559   0.10676     
#> 2       4       8 5.0654  4.075e-07 1.223e-06  ***
#> 3       6       8 2.2138    0.02684   0.08053     
#> -----------------------------------------
#> *, P<0.05; **, P<0.01; ***, P<0.001
```
