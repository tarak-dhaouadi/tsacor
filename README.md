# tsacor

Trial Sequential Analysis (TSA) for meta-analyses of correlations, in R:
Pearson product-moment correlations (r) and Spearman rank correlations
(rho), pooled on Fisher's z scale.

**Standalone Java application**
A standalone Java implementation of `tsacor` is available as [tsacor-java](https://github.com/tarak-dhaouadi/tsacor-java), providing the same TSA methodology without requiring R.

Adapts the classical Wetterslev/Thorlund/Copenhagen Trial Unit TSA
framework to correlation coefficients. It's worth distinguishing what's
established methodology versus what this package specifically contributes:

**Established components** (from the cited literature):
- the **Fisher z transformation** (`atanh(r)`) of Fisher (1921), on which
  the study-level correlations are pooled and back-transformed with
  `tanh()`,
- the single-study sample-size formula for correlations,
  `n = c (z_(1-alpha/2) + z_power)^2 / atanh(r0)^2 + 3`, here used as the
  required-information-size analogue of the Schoenfeld formula for hazard
  ratios,
- the **variance factor `c`** of Fisher's z for Spearman's rho
  (1.06, Fieller-Hartley-Pearson 1957; or `1 + rho^2/2`, Bonett-Wright
  2000),
- the **Diversity (D²)** heterogeneity adjustment of Wetterslev et al.
  (2009),
- the general trial-sequential-monitoring (alpha/beta-spending) framework.

**This package's specific implementation choices:**
- a correlation-specific adaptation combining the above into one workflow,
- **inverse-variance information** based on each study's own standard error
  of Fisher's z (`sum(1/se_z^2)`), taken either from the study's reported
  confidence interval (`se_source = "ci"`, the default) or from its sample
  size (`se_source = "n"`), used as the accrued-information measure,
- a Z-curve drawn against the **cumulative number of participants**
  (`n_subjects`) rather than cumulative events, with the DARIS translated
  to a single-study-equivalent number of participants via the variance
  factor `c`,
- a compiled (C++) recursive numerical integration engine, ported from
  RTSA, for the O'Brien-Fleming-type alpha- and beta-spending boundaries
  (no external group-sequential-design package, and no fixed
  software-imposed limit on the number of looks, subject to available
  computational resources) -- see `?tsa_cor`, `src/rtsa_core.h` and
  `inst/REVERSE_ENGINEERING_RTSA.md`. An R-only fallback engine
  (`R/obf_boundaries.R`) is kept and is ON by default
  (`legacy_fallback = TRUE`); it only runs if the compiled engine fails,
  and is clearly flagged in the result whenever it is used; set
  `legacy_fallback = FALSE` for confirmatory or RTSA-parity work,
- applying this monitoring framework to a **cumulative random-effects**
  meta-analysis Z-curve -- see the Caveats section below, this is an
  approximation shared with the official Copenhagen Trial Unit TSA
  software, not an exact result.

## Installation

```r
# install.packages("remotes")
remotes::install_github("tarak-dhaouadi/tsacor")
```

You'll also need its dependencies if you don't already have them:

```r
install.packages(c("metafor", "readxl", "ggplot2"))
```

## Usage

```r
library(tsacor)

# Try it on the bundled example dataset: 20 studies (2005-2024)
path <- tsacor_example_data()

res <- tsa_cor(
  data            = path,
  target_r        = 0.10,     # pre-specified anticipated correlation (recommended)
  alpha_two_sided = 0.05,
  power           = 0.80,
  order_by        = "Year"    # TSA is order-dependent: sort chronologically
)

summary(res)   # full results table
plot(res)      # the TSA chart
```

The subtitle of the chart shows the model/design summary and, on a second
line, the pooled random-effects correlation with its 95% CI, the p-value,
tau² (Fisher z scale) and I². With `boundary_route = "analysis"`, the
position and size of the "Analysis-route endpoint ... reached" label can be
set with `endpoint_label_x`, `endpoint_label_y` and `endpoint_label_size`.

By default the four DARIS-related labels (theoretical DARIS, historical-rate
projection, "DARIS information reached" and analysis-route endpoint) are drawn
in the upper part of the plot when the final Z-score is negative and in the
lower part when it is positive, with the "Participants accrued" label on the
opposite side (over the curve). Use `daris_label_y`, `historical_label_y`,
`info_threshold_label_y`, `endpoint_label_y` and `participants_label_y` to
override these defaults.

### Pearson or Spearman

`cor_type = "pearson"` (default; shorthand `"r"`) or `"spearman"`
(shorthand `"rho"`) states which coefficient the `r` column contains; all
studies must report the same kind. For Spearman's rho, `spearman_variance`
selects the variance model of Fisher's z, used to translate information
into participants and, with `se_source = "n"`, to compute the standard
errors:

```r
res <- tsa_cor(path, target_r = 0.10, cor_type = "spearman",
               spearman_variance = "bonett_wright", order_by = "Year")
```

### Source of the standard errors

By default (`se_source = "ci"`) the standard error of each study's Fisher z
is derived from its reported confidence interval,
`(atanh(ubound) - atanh(lbound)) / (2 q)`, which also works for Spearman
coefficients whose intervals were computed by other methods. With
`se_source = "n"` it is derived from the sample size only,
`sqrt(c / (n - 3))`. Both versions are returned in `res$data` (`se_z_ci`
and `se_z_n`), so the choice can be checked as a sensitivity analysis; if
your intervals are not 95% intervals, set `ci_level` accordingly.

### Random-effects inference

By default (`re_inference = "standard"`) the random-effects model uses the
usual normal-theory inference. `re_inference = "hksj"` (alias `"knha"`)
switches to the Hartung-Knapp-Sidik-Jonkman adjustment, and
`re_inference = "hksj_adhoc"` (alias `"knha_adhoc"`) to its ad hoc variant
in which the variance multiplier is never below 1:

```r
res <- tsa_cor(path, target_r = 0.10, order_by = "Year", re_inference = "hksj")
plot(res)   # the caption names the inference option
```

The tau² estimator (`method`), DARIS and the boundaries are unaffected; the
pooled CI/p-value and the cumulative Z-curve (shown as the normal-equivalent
of the HKSJ t statistic) change. The HKSJ statistic is undefined at the
first look (`cumulative$Z` is `NA` there, and no point is drawn), and early
looks generally are unstable (very few degrees of freedom) -- `tsa_cor()`
warns about this whenever a non-standard option is used; see `?tsa_cor`.

### Using your own data

`data` can be a data.frame or a path to an `.xlsx` file with one row per
study and (at least) these columns:

| Column       | Meaning                                                                                  |
|--------------|-------------------------------------------------------------------------------------------|
| `Study`      | Unique study label                                                                        |
| `r`          | Pearson correlation or Spearman rho of that study (strictly between -1 and 1)             |
| `n_subjects` | number of subjects the correlation is based on (whole number, greater than 3)             |
| `lbound`     | lower limit of the study's confidence interval, on the correlation scale (needed when `se_source = "ci"`, the default) |
| `ubound`     | upper limit of the study's confidence interval, on the correlation scale (needed when `se_source = "ci"`, the default) |

Any other columns (the bundled example also has `Year`, `Ethnicity` and
`Age`) are kept in the returned `data` and can be used with `order_by`.
Column names have spaces replaced with underscores on load. Each row is
treated as one independent study; at least two studies are required, and
fewer than 10 raise a warning because the heterogeneity and Diversity
estimates, and therefore DARIS and the boundaries, can be very unstable.

Rows should be in **chronological (publication) order** — TSA results
are order-dependent. Either pre-sort your data yourself, or pass
`order_by = "<column name>"` (e.g. a publication-year column) to have
`tsa_cor()` sort it for you explicitly:

```r
res <- tsa_cor(path, target_r = 0.10, order_by = "Year")
```

### Structure of the example dataset

`tsacor_example_data()` returns the path to the bundled `r_meta.xlsx`
(20 studies, 2005-2024):

| Column       | Type      | Meaning                                                 | Used by `tsa_cor()`            |
|--------------|-----------|----------------------------------------------------------|--------------------------------|
| `Study`      | character | Unique study label (`Study_01` ... `Study_20`)           | yes (required)                 |
| `Year`       | numeric   | Publication year                                         | via `order_by = "Year"`        |
| `Ethnicity`  | character | Ethnic group of the study population                     | no (kept in `res$data`)        |
| `Age`        | numeric   | Age of the study population                              | no (kept in `res$data`)        |
| `r`          | numeric   | Correlation coefficient of the study                     | yes (required)                 |
| `lbound`     | numeric   | Lower limit of the 95% CI of `r`                         | yes (`se_source = "ci"`)       |
| `ubound`     | numeric   | Upper limit of the 95% CI of `r`                         | yes (`se_source = "ci"`)       |
| `n_subjects` | numeric   | Number of subjects the correlation is based on           | yes (required)                 |

For example, the first study is `Study_01` (2005, European, mean age 42):
`r = 0.534`, 95% CI 0.384 to 0.654, `n_subjects = 120`.

### Saving the plot / results

`tsa_cor()` does not write files itself; use standard R tools:

```r
p <- plot(res)
ggplot2::ggsave("tsa_plot.png", p, width = 11, height = 7.5, dpi = 300)

write.csv(res$cumulative,    "tsa_cumulative_results.csv", row.names = FALSE)
write.csv(res$summary_table, "tsa_summary.csv",             row.names = FALSE)
```

## Important caveats

- **Circularity of `target_r = NA`:** if you don't specify `target_r`,
  the required information size is calculated from the *observed* pooled
  correlation, which is circular and can make the TSA boundary collapse to
  the conventional boundary almost immediately. Set `target_r` to a
  pre-specified, scientifically meaningful correlation for a standard,
  publication-quality TSA. Because the required information grows as
  `1 / atanh(r0)^2`, targets close to 0 (|r0| < 0.10 raises a warning)
  demand very large amounts of information. See `?tsa_cor` for details.
- **Random-effects approximation:** the plotted Z-curve is a cumulative
  random-effects meta-analysis, whose between-study variance is
  re-estimated at every step. The Lan-DeMets/O'Brien-Fleming monitoring
  boundaries strictly assume a fixed-information canonical process, so
  applying them to a random-effects Z-curve is a standard approximation
  (shared with the official Copenhagen Trial Unit TSA software), not an
  exact result.
- **Participants are a single-study equivalent:** the participant-based
  DARIS (`c * DARIS + 3`) is a single-study-equivalent reference, not an
  exact count of the participants needed across studies, because a
  meta-analysis of k studies loses 3 units of information per study whereas
  the equivalent loses them only once. All "reached" verdicts use the
  observed information (`sum(1/se_z^2)`) against DARIS instead.
- **`method` changes more than just the pooled effect estimate:** the
  heterogeneity-variance (`tau^2`) estimator selected via `method` (default
  `"DL"`) does **not** change the mathematical alpha-spending function or
  the boundary-calculation algorithm -- those are a fixed part of the
  group-sequential design chosen up front (`alpha_two_sided`, `power`).
  However, because `method` changes `tau^2`, it also changes the pooled
  SE, the random-effects cumulative Z-curve, D-squared, DARIS, and the
  cumulative information schedule -- and therefore *which study
  corresponds to which information fraction*. So while the spending
  function itself is unaffected, switching, e.g., `method = "DL"` to
  `method = "REML"` can still change the boundary values attached to
  the observed looks indirectly, by changing the information schedule
  those looks land on, and therefore the practical timing of a boundary
  crossing.
- **Projections are indicative:** the estimated additional studies and
  participants needed to reach DARIS are linear extrapolations at the
  observed historical rate of information; they do not model changes in
  tau² or in the size and precision of future studies.

## References

Wetterslev J, Thorlund K, Brok J, Gluud C. "Estimating required
information size by quantifying diversity in random-effects model
meta-analyses." *BMC Med Res Methodol.* 2009;9:86.

Fisher RA. "On the 'probable error' of a coefficient of correlation
deduced from a small sample." *Metron.* 1921;1:3-32.

Fieller EC, Hartley HO, Pearson ES. "Tests for rank correlation
coefficients. I." *Biometrika.* 1957;44:470-481.

Bonett DG, Wright TA. "Sample size requirements for estimating Pearson,
Kendall and Spearman correlations." *Psychometrika.* 2000;65:23-28.

Hartung J, Knapp G. "On tests of the overall treatment effect in
meta-analysis with normally distributed responses." *Stat Med.*
2001;20:1771-1782.

Sidik K, Jonkman JN. "A simple confidence interval for meta-analysis."
*Stat Med.* 2002;21:3153-3164.

## Attribution and license

tsacor contains code ported from the R package
[RTSA](https://cran.r-project.org/package=RTSA) (Anne Lyngholm Soerensen,
Markus Harboe Olsen, Theis Lange and Christian Gluud), licensed GPL (>= 2):
the C++ boundary engine (`src/rtsa_core.h`, `src/rtsa_engine.cpp`), its R
orchestration (`R/rtsa_engine.R`) and the R-only reconstruction of it
(`R/obf_boundaries.R`). RTSA is the R version of Trial Sequential Analysis
(TSA), originally developed as a stand-alone Java program by the Copenhagen
Trial Unit; the RTSA manual is heavily inspired by the user manual for TSA by
Kristian Thorlund, Janus Engstrøm, Jørn Wetterslev, Jesper Brok, Georgina
Imberger and Christian Gluud. The original TSA software is available at
<https://ctu.dk/tools>:

> Copenhagen Trial Unit, Centre for Clinical Intervention Research,
> Department 3344, Rigshospitalet, DK-2100 Copenhagen Ø, Denmark.
> Tel. +45 3545 7171, Fax +45 3545 7101, E-mail: tsa@ctu.dk

Because of that, **tsacor is licensed GPL (>= 2)**. See `inst/COPYRIGHTS`
for the file-by-file provenance. If you use tsacor for boundary computations
please also cite RTSA and the TSA software.
