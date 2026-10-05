# tsacor 0.1.1

* `plot.tsa_cor()`: the default vertical placement of the labels now depends on
  the sign of the final cumulative Z-score. When the Z-curve is positive, the
  four DARIS-related labels (theoretical DARIS participant-equivalent,
  historical-rate projection, "DARIS information reached" and analysis-route
  endpoint) are placed in the lower part of the plot and the "Participants
  accrued" label over the curve in the upper part; when the Z-curve is negative
  (or zero) the placement is unchanged (DARIS labels in the upper part,
  "Participants accrued" in the lower part). These are defaults only: the
  `daris_label_y`, `info_threshold_label_y`, `endpoint_label_y`,
  `historical_label_y` and `participants_label_y` arguments still override them.

# tsacor 0.1.0

* First release of tsacor: Trial Sequential Analysis for meta-analyses of
  Pearson product-moment correlations and Spearman rank correlations, pooled
  on Fisher's z scale.
* `tsa_cor()`: random-effects meta-analysis of Fisher's z (DerSimonian-Laird by
  default, with `method` offering the other `metafor::rma()` estimators that
  need no extra arguments; optional Hartung-Knapp-Sidik-Jonkman inference via
  `re_inference`); heterogeneity (Q, I^2, tau^2) and the Wetterslev et al.
  (2009) Diversity (D^2) adjustment; required information size and
  Diversity-Adjusted Required Information Size (DARIS) from a pre-specified
  or observed target correlation, translated to a single-study-equivalent
  number of participants via the Fisher-z variance factor `c` (1 for Pearson
  r; 1.06, Fieller-Hartley-Pearson 1957, or `1 + rho^2/2`, Bonett-Wright 2000,
  for Spearman rho); cumulative (sequential) analysis against the cumulative
  number of participants; O'Brien-Fleming-type alpha- and beta-spending
  trial sequential monitoring boundaries, computed with the RTSA-derived
  recursive-integration engine shared with the sister package tsahr (same
  boundary computations tsahr uses for hazard ratios; see
  `inst/COPYRIGHTS` and `inst/REVERSE_ENGINEERING_RTSA.md`), including both
  RTSA's "design" and "analysis" retrospective boundary routes and a
  legacy R-only fallback engine.
* `print.tsa_cor()`, `summary.tsa_cor()` and `plot.tsa_cor()` methods; the plot
  shows the cumulative Z-curve, the monitoring boundaries, the conventional
  significance boundary and both the theoretical and observed-information
  DARIS participant-equivalent reference lines.
* `tsacor_example_data()` and the bundled `r_meta` dataset (20 studies,
  2005-2024) for examples and unit tests.
* `cor_type` also accepts the shorthands `"r"` (for `"pearson"`) and `"rho"`
  (for `"spearman"`), matched case- and whitespace-insensitively.
* Documentation: the `tsa_cor()` help page now documents every argument in the
  detail of its sister function `tsa_hr()`, in particular the `method` (tau^2
  estimator) argument (what it changes, what it does not, aliases,
  unsupported strings, convergence), the `re_inference` options, the data
  validation rules and warnings, the retrospective projection and the RTSA
  quirks kept in the boundary engine. `?plot.tsa_cor` now lists all plot
  arguments.
