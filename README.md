# tsacor

Trial Sequential Analysis (TSA) for meta-analyses of Pearson product-moment
correlations (r) and Spearman rank correlations (rho), pooled on Fisher's z
scale.

tsacor is the sister package of
[tsahr](https://github.com/tarak-dhaouadi/tsahr) (TSA for meta-analyses of
hazard ratios): it reuses the same RTSA-derived alpha-/beta-spending boundary
engine, adapted here to the information scale of correlations (inverse
variance of Fisher's z) and to cumulative *participants* rather than
cumulative *events*.

## Installation

```r
# install.packages("remotes")
remotes::install_github("tarak-dhaouadi/tsacor")
```

## Usage

```r
library(tsacor)

path <- tsacor_example_data()      # bundled example: r_meta.xlsx, 20 studies
res  <- tsa_cor(path, target_r = 0.10, order_by = "Year")

summary(res)
plot(res)
```

`data` needs, at minimum, the columns `Study`, `r` (a Pearson correlation or
Spearman rho) and `n_subjects`; by default (`se_source = "ci"`) it also needs
`lbound`/`ubound`, the 95% confidence limits of each study's correlation.
See `?tsa_cor` for the full set of options (`cor_type`, `se_source`,
`spearman_variance`, `method`, `re_inference`, `boundary_route`, ...).

## Methodology

* Wetterslev J, Thorlund K, Brok J, Gluud C. "Estimating required information
  size by quantifying diversity in random-effects model meta-analyses." BMC
  Med Res Methodol. 2009;9:86.
* Fieller EC, Hartley HO, Pearson ES. "Tests for rank correlation
  coefficients. I." Biometrika. 1957;44:470-481.
* Bonett DG, Wright TA. "Sample size requirements for estimating Pearson,
  Kendall and Spearman correlations." Psychometrika. 2000;65:23-28.
* The boundary engine is a compiled port of
  [RTSA](https://CRAN.R-project.org/package=RTSA) (Soerensen, Olsen, Lange
  and Gluud), the R implementation of the Copenhagen Trial Unit's Trial
  Sequential Analysis software. See `inst/COPYRIGHTS`.

## License

GPL (>= 2). See `inst/COPYRIGHTS` for the origin and copyright of the parts
ported from RTSA.
