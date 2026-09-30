## *** LEGACY-ENGINE TESTS ***  These check the spending FORMULA and the OLD
## FFT alpha engine (.obf_alpha_boundary(), now only tsa_cor()'s fallback). The
## current engine's alpha recursion is checked in test-rtsa-engine-parity.R
## (RTSA's own alpha reference 4.877 3.357 2.680 2.290 2.031 to 1e-3).

## Regression test for the alpha-spending formula correction: verifies
## the two-sided (side = 2) O'Brien-Fleming-type spending function
## against a concrete, exact, approximation-free comparison point.
##
## The FIRST look of any group-sequential design has no prior boundary
## to condition on, so its boundary is an exact closed-form function of
## the spending function value at that look -- c_1 = qnorm(1 - A(t_1)/2)
## -- with no recursive-integration approximation involved. This makes
## it a clean, engine-independent check of the spending FORMULA itself
## (as opposed to the recursive engine's numerics), which is exactly
## what needed checking here.
##
## Reference value is live output from RTSA::boundaries() (Copenhagen
## Trial Unit's R package -- the reference TSA implementation this
## package models itself on):
##   RTSA::boundaries(timing = c(0.2,0.4,0.6,0.8,1), alpha = 0.05,
##                     side = 2, es_alpha = "esOF")
## reported a first-look upper boundary of 4.877 at an (internally
## event-count-rounded) timing of t = 0.205.

test_that("alpha-spending formula matches RTSA::boundaries() at the first look", {
  t1 <- 0.205
  alpha <- 0.05

  A1 <- tsacor:::.alpha_spend_OF(t1, alpha)
  c1 <- stats::qnorm(1 - A1 / 2)

  rtsa_reported <- 4.877

  ## Loose-ish tolerance: RTSA's reported timing (0.205) is itself
  ## rounded from an underlying integer event count, so an exact bit-for-
  ## bit match isn't expected -- but the corrected formula should land
  ## within a few hundredths, not the ~0.5 discrepancy the pre-correction
  ## formula had.
  expect_equal(c1, rtsa_reported, tolerance = 0.1)
})

test_that("alpha-spending formula reaches exactly alpha at t = 1", {
  alpha <- 0.05
  expect_equal(tsacor:::.alpha_spend_OF(1, alpha), alpha, tolerance = 1e-12)
})

test_that("alpha-spending formula is the side=2 (two-sided) OF-type form", {
  t <- c(0.1, 0.3, 0.5, 0.7, 1)
  alpha <- 0.05
  got <- tsacor:::.alpha_spend_OF(t, alpha)
  expected <- 4 * stats::pnorm(
    stats::qnorm(1 - alpha / 4) / sqrt(t),
    lower.tail = FALSE
  )
  expect_equal(got, expected, tolerance = 1e-14)
})

## The checks above pin the spending FORMULA, but a coarse default grid
## in the recursive engine could in principle still reproduce a
## formula-correct-but-numerically-off set of boundaries. This test pins
## the full recursive ENGINE output against all 5 RTSA reference values
## directly (using the nominal timing actually passed to
## RTSA::boundaries(), not its rounded reported SMA_Timing column), so a
## regression to a coarser grid or any other engine-level numerical
## issue would be caught here even though the formula-only checks above
## would still pass. NOTE: an earlier draft of this file's comments (and
## of the VALIDATION note in R/obf_boundaries.R) claimed the
## then-default n_grid=2000 produced ~0.03 error at the final look,
## specifically. That specific figure could not be reproduced by an
## independent from-scratch Python port of this exact algorithm (which
## instead got ~0.003-0.006 max error at n_grid=2000 using the correct
## nominal timing) and was never confirmed by an actual run of this R
## code -- treat it as unverified rather than as an established
## regression case. n_grid's default is kept at 16000 regardless, as a
## conservative, low-cost accuracy margin rather than a fix for a
## confirmed problem at 2000.
test_that("full recursive alpha boundary engine matches RTSA::boundaries() at all 5 looks", {
  t <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  alpha <- 0.05

  got <- tsacor:::.obf_alpha_boundary(t, alpha)
  rtsa_reported <- c(4.877, 3.357, 2.680, 2.290, 2.031)

  ## Tolerance reflects RTSA's own 3-decimal reporting precision, not an
  ## exact bit-for-bit target. Based on an independent Python port of
  ## this algorithm this is expected to pass even at the old
  ## n_grid=2000 default, not just at the current 16000 -- see the note
  ## above and the VALIDATION note in R/obf_boundaries.R.
  expect_equal(got, rtsa_reported, tolerance = 0.01)
})
