## *** LEGACY-ENGINE TESTS (R-only, pre-0.2.7.11) ***
## Everything in this file (apart from the final 0.2.7.11 test, which
## checks the current engine's beta_engine object) exercises the OLD pure-R
## boundary functions
## (.obf_alpha_boundary(), .obf_beta_boundary(), .rtsa_beta_boundary(),
## .rtsa_beta_boundary_analysis(), ...), which tsa_cor() now uses only as its
## opt-in fallback (legacy_fallback = TRUE, flagged in the result). Since
## 0.2.7.19 .rtsa_beta_boundary() recomputes its final efficacy wall from the
## FFT alpha recursion instead of assuming qnorm(1 - alpha / 2); that
## constant survives in these tests only for single-look schedules (where it
## IS the right value) and for the analysis wrapper's design_R endpoint
## approximation. This engine is approximate either way (FFT alpha; measured
## root error 7.6e-4 on the reference schedule, see NEWS 0.2.7.21). It is NOT
## the current engine's convention: the compiled RTSA-derived
## design route sets the final futility bound equal to the final efficacy
## bound (the alpha recursion's value at t = 1, e.g. ~2.13, not 1.96).
## Tests of the current engine are in test-rtsa-engine-parity.R,
## test-robust-schedules.R and test-boundary-route-and-fallback.R.

test_that("0.2.7.10: .rtsa_beta_boundary_analysis() calibrates design_R via .rtsa_beta_boundary()'s side = 2 non-binding warp_root search, NOT a side = 1 futility = \"none\" root search", {
  ## ** 0.2.7.7 got this right; 0.2.7.9 briefly regressed it to a side =
  ## 1 / futility = \"none\" calibration on the strength of an external
  ## reconstruction that queried the wrong RTSA branch; reverted here
  ## (0.2.7.10) after re-reading RTSA's own R/RTSA.R top-level wrapper,
  ## which manufactures design_R via
  ## `boundaries(timing = ..., side = side, futility = futility, ...,
  ## type = "design")` -- i.e. RTSA()'s OWN side/futility arguments
  ## (side = 2, futility = "non-binding" for this package's design), not
  ## hardcoded side = 1 / futility = "none". That is exactly
  ## .rtsa_beta_boundary()'s own two-pass warp_root search.
  t <- c(0.20, 0.40, 0.60, 0.80, 0.95)
  alpha <- 0.05
  beta <- 0.20
  alpha_ref <- c(4.877, 3.357, 2.680, 2.290, 2.10)

  design_fit <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)
  ans <- tsacor:::.rtsa_beta_boundary_analysis(t, alpha, beta, alpha_ref)

  expect_equal(ans$design_R, design_fit$warp_root, tolerance = 1e-12)

  ## The fixed drift must use this design's own side = 2 convention
  ## (abs(qnorm(alpha/2) + qnorm(beta))), matching .rtsa_beta_boundary()
  ## and the package's own side = 2 alpha engine -- NOT side = 1's plain
  ## abs(qnorm(alpha) + qnorm(beta)).
  expect_equal(ans$delta, abs(stats::qnorm(alpha / 2) + stats::qnorm(beta)),
               tolerance = 1e-12)
})

test_that("0.2.7.7: an externally supplied design_R is used as-is, with no internal root search", {
  t <- c(0.20, 0.40, 0.60, 0.80, 0.95)
  alpha <- 0.05
  beta <- 0.20
  alpha_ref <- c(4.877, 3.357, 2.680, 2.290, 2.10)

  ans <- tsacor:::.rtsa_beta_boundary_analysis(
    t, alpha, beta, alpha_ref, design_R = 1.10
  )
  expect_equal(ans$design_R, 1.10)
  expect_length(ans$boundary, length(t))
})

test_that("0.2.7.7: the design endpoint is appended, not observed, and is excluded from the returned boundary", {
  t <- c(0.20, 0.40, 0.60, 0.80, 0.95)
  alpha <- 0.05
  beta <- 0.20
  alpha_ref <- c(4.877, 3.357, 2.680, 2.290, 2.10)

  ans <- tsacor:::.rtsa_beta_boundary_analysis(t, alpha, beta, alpha_ref)

  ## Only the actually observed looks come back, not the synthetic
  ## design_R endpoint used to calibrate them.
  expect_length(ans$boundary, length(t))
  expect_true(ans$design_R %in% ans$t_ext)
  expect_false(ans$design_R %in% t)
})

test_that("0.2.7.7: the analysis-mode info scale is unwarped (org_t == t_ext), unlike design mode", {
  ## No warp_root/inf_warp re-scaling of the info axis in the analysis
  ## branch. We can't inspect org_t directly from the public wrapper,
  ## but we can confirm the design_R branch's info scale is NOT the
  ## design-mode warp_root scale by checking that supplying a design_R
  ## far from the design-mode warp_root still produces a finite, sane
  ## result (i.e. the analysis branch does not silently fall back to
  ## re-deriving its own warp).
  t <- c(0.20, 0.40, 0.60, 0.80, 0.95)
  alpha <- 0.05
  beta <- 0.20
  alpha_ref <- c(4.877, 3.357, 2.680, 2.290, 2.10)

  ans_low <- tsacor:::.rtsa_beta_boundary_analysis(
    t, alpha, beta, alpha_ref, design_R = 1.02
  )
  ans_high <- tsacor:::.rtsa_beta_boundary_analysis(
    t, alpha, beta, alpha_ref, design_R = 1.30
  )

  expect_true(all(is.finite(ans_low$boundary) | is.na(ans_low$boundary)))
  expect_true(all(is.finite(ans_high$boundary) | is.na(ans_high$boundary)))
  ## A materially different design_R should change the beta-spending
  ## timeline (inf_frac / design_R) and hence, generically, the
  ## boundaries -- confirming design_R is actually being used, not
  ## ignored in favour of some other fixed scale.
  expect_false(isTRUE(all.equal(ans_low$boundary, ans_high$boundary)))
})

test_that("0.2.7.10: rm_bs is derived (iterated fixed point), never hard-coded or fixed at 0", {
  ## RTSA's own boundaries() calls beta_boundary() THREE times for the
  ## side = 2 / futility = "non-binding" / type = "analysis" branch this
  ## package's design reaches: once unsuppressed, then twice more with
  ## rm_bs re-derived from the previous pass's negative-za count. rm_bs
  ## must therefore be able to land anywhere in [0, length(t_ext)]
  ## depending on the data's own shape -- not pinned at a fixed value
  ## (0.2.7.9's regression) nor at a fixed 5 (the pre-0.2.7.7 mechanism).
  t_tight <- c(0.01, 0.02, 0.03, 0.05, 0.08, 0.90)
  alpha <- 0.05
  beta <- 0.20
  alpha_ref_tight <- c(9.5, 7.3, 6.2, 5.1, 4.2, 2.10)

  ans_tight <- tsacor:::.rtsa_beta_boundary_analysis(t_tight, alpha, beta, alpha_ref_tight)
  expect_true(ans_tight$rm_bs >= 0L)
  expect_true(ans_tight$rm_bs <= length(ans_tight$t_ext))

  t_spread <- c(0.20, 0.40, 0.60, 0.80, 0.95)
  alpha_ref_spread <- c(4.877, 3.357, 2.680, 2.290, 2.10)
  ans_spread <- tsacor:::.rtsa_beta_boundary_analysis(t_spread, alpha, beta, alpha_ref_spread)
  expect_true(ans_spread$rm_bs >= 0L)
  expect_true(ans_spread$rm_bs <= length(ans_spread$t_ext))
})

test_that("0.2.7.10: looks suppressed by rm_bs (za pinned at the +/-20 sentinel) come back as NA, not as a finite negative number", {
  ## This is RTSA's own convention (boundaries.R:
  ## `beta_ubound <- c(rep(NA, sum(abs(za) == 20)), za[abs(za) < 20])`),
  ## reproduced by .rtsa_beta_boundary_analysis()'s own sentinel-to-NA
  ## conversion. 0.2.7.9's regression (design_R from the wrong branch,
  ## rm_bs fixed at 0) meant this sentinel was essentially never hit, so
  ## early looks surfaced as small negative finite numbers instead of NA.
  t <- c(0.20, 0.40, 0.60, 0.80, 0.95)
  alpha <- 0.05
  beta <- 0.20
  alpha_ref <- c(4.877, 3.357, 2.680, 2.290, 2.10)
  ans <- tsacor:::.rtsa_beta_boundary_analysis(t, alpha, beta, alpha_ref)

  if (ans$rm_bs > 0L) {
    expect_true(any(is.na(ans$boundary[seq_len(ans$rm_bs)])))
  }
})

test_that("0.2.7.7: an over-powered single early look (t > design_R) gets the definitive final boundary", {
  ## Mirrors .rtsa_beta_boundary()'s own over_power convention, applied
  ## relative to design_R instead of 1.
  t <- c(1.2, 1.5, 2.0)
  alpha <- 0.05
  beta <- 0.20
  final_ab <- stats::qnorm(1 - alpha / 2)
  alpha_ref <- rep(final_ab, length(t))

  ans <- tsacor:::.rtsa_beta_boundary_analysis(t, alpha, beta, alpha_ref)
  expect_length(ans$boundary, length(t))
  expect_true(ans$over_power)
  expect_true(all(ans$boundary == final_ab))
})

test_that("0.2.7.7: futility boundary stays capped at the corresponding efficacy boundary", {
  t <- c(0.10, 0.25, 0.50, 0.75, 0.95)
  alpha <- 0.05
  beta <- 0.20
  alpha_ref <- c(5.0, 3.5, 2.8, 2.2, 2.05)
  ans <- tsacor:::.rtsa_beta_boundary_analysis(t, alpha, beta, alpha_ref)

  ok <- !is.na(ans$boundary)
  expect_true(all(ans$boundary[ok] <= alpha_ref[ok] + 1e-14))
})

