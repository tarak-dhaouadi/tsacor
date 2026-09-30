## 0.2.7.11 -- RTSA-derived compiled engine (src/rtsa_core.h, R/rtsa_engine.R).
##
## Reference values below were produced by a Python port of
## RTSA 0.2.2's R sources (tools/rtsa_port.py) (alpha_boundary(), beta_boundary(), boundaries()),
## which itself reproduces (i) RTSA::boundaries()'s alpha reference used
## elsewhere in this test suite, and (ii) the printed outputs of RTSA's
## "futility" vignette (SMA timing 0.541/0.812/1.083, FutLower
## 0.332/1.292/2.014).  The compiled C++ core agreed with that port to
## ~1e-14 in a standalone (no R) build.  The R orchestration was NOT run
## when these values were written (no R interpreter available) -- if a test
## here fails, compare against a live RTSA::boundaries() call first
## (see the last block).

t8 <- c(0.461092, 0.527034, 0.610745, 0.653952,
        0.754067, 0.825699, 0.903255, 0.961096)

test_that("alpha engine reproduces RTSA's alpha reference (Simpson recursion)", {
  got <- tsacor:::.rtsa_alpha_cpp(c(0.2, 0.4, 0.6, 0.8, 1), side = 2L, alpha = 0.05)
  expect_equal(got$alpha_ubound, c(4.877, 3.357, 2.680, 2.290, 2.031),
               tolerance = 1e-3)
  ## 4 equally spaced looks, one-sided 0.025: 2.963 2.359 2.014 (RTSA vignette)
  got1 <- tsacor:::.rtsa_alpha_cpp(c(0.5, 0.75, 1), side = 1L, alpha = 0.025)
  expect_equal(got1$alpha_ubound, c(2.963, 2.359, 2.014), tolerance = 1e-3)
})

test_that("side = 1 non-binding design reproduces RTSA's published vignette output", {
  alpha <- 0.025; beta <- 0.1; t <- c(0.5, 0.75, 1)
  ub <- tsacor:::.rtsa_alpha_cpp(t, 1L, alpha)$alpha_ubound
  delta <- abs(stats::qnorm(alpha) + stats::qnorm(beta))
  f <- function(x) {
    ub[3] - tsacor:::.rtsa_beta_cpp(t, ub, beta, delta, warp_root = x)$za[3]
  }
  root <- stats::uniroot(f, lower = 0.9, upper = 1.5, tol = 1e-9)$root
  za <- tsacor:::.rtsa_beta_cpp(t, ub, beta, delta, warp_root = root)$za
  expect_equal(round(t * root, 3), c(0.541, 0.812, 1.083))
  expect_equal(round(za, 3), c(0.332, 1.292, 2.014))
})

test_that("design route: final efficacy wall is the alpha recursion's value, not qnorm(1 - alpha/2)", {
  des <- tsacor:::.rtsa_design_bounds(t8, alpha = 0.05, beta = 0.20)
  expect_equal(des$timing, c(t8, 1))
  expect_equal(des$root, 1.210717, tolerance = 1e-5)
  expect_equal(des$rm_bs, 0L)
  expect_equal(des$alpha_ubound,
               c(3.101131, 2.923183, 2.699739, 2.634717, 2.412119,
                 2.316074, 2.209634, 2.153588, 2.127034), tolerance = 2e-4)
  ## the beta bounds reported by RTSA for these fractions
  expect_equal(des$beta_ubound[1:8],
               c(0.531029, 0.697635, 0.980490, 1.076945,
                 1.400023, 1.567680, 1.768009, 1.929516), tolerance = 2e-4)
  ## by construction the futility bound meets the efficacy wall at t = 1
  expect_equal(des$beta_ubound[9], des$alpha_ubound[9], tolerance = 1e-7)
  expect_gt(des$alpha_ubound[9], stats::qnorm(0.975) + 0.1)
})

test_that("design route: early looks with negative first-pass bounds are suppressed to NA", {
  t7 <- c(0.05, 0.12, 0.20, 0.35, 0.50, 0.70, 0.90)
  des <- tsacor:::.rtsa_design_bounds(t7, alpha = 0.05, beta = 0.20)
  expect_equal(des$rm_bs, 4L)
  expect_true(all(is.na(des$beta_ubound[1:4])))
  expect_equal(des$beta_ubound[5:8],
               c(0.671083, 1.292845, 1.799751, 2.070360), tolerance = 2e-4)
  expect_equal(des$root, 1.175054, tolerance = 1e-5)
})

test_that("analysis route (RTSA() retrospective): design_R, t / design_R spending, unwarped information", {
  ret <- tsacor:::.rtsa_retrospective(t8, alpha = 0.05, beta = 0.20)
  expect_equal(ret$design_R, 1.210717, tolerance = 1e-5)
  ana <- ret$analysis
  expect_equal(ana$timing, c(t8, ret$design_R), tolerance = 1e-12)
  ## alpha bounds are recomputed on t / design_R -- NOT the design-type ones
  expect_equal(ana$alpha_ubound,
               c(3.449109, 3.240996, 2.993992, 2.917153, 2.676108,
                 2.567387, 2.449628, 2.385720, 2.035924), tolerance = 2e-4)
  expect_equal(ana$beta_ubound[1:8],
               c(0.125985, 0.303716, 0.587548, 0.683736,
                 1.001699, 1.163209, 1.344194, 1.455898), tolerance = 2e-4)
  ## last look clipped to the efficacy wall, as RTSA's boundaries() does
  expect_equal(ana$beta_ubound[9], ana$alpha_ubound[9], tolerance = 1e-9)
  ## the two routes are genuinely different computations
  expect_gt(max(abs(ret$design$beta_ubound[1:8] - ana$beta_ubound[1:8])), 0.3)
})

test_that("compiled init_int/recur_int/prob agree with a direct R evaluation of first.cpp's formulas", {
  set.seed(1)
  zj <- seq(-3, 3, length.out = 21); wj <- rep(0.3, 21)
  stdv <- cbind(sqrt(c(0.4, 0.3, 0.3)), sqrt(c(0.4, 0.7, 1.0)))
  delta <- 2.8
  last <- tsacor:::rtsa_init_int_cpp(wj, zj, delta, stdv[, 1])
  expect_equal(last, wj * stats::dnorm(zj, delta * stdv[1, 1], 1), tolerance = 1e-14)

  zj_up <- seq(-2.5, 2.5, length.out = 15); wj_up <- rep(0.2, 15)
  up <- tsacor:::rtsa_recur_int_cpp(2L, stdv, zj, last, zj_up, wj_up, delta, FALSE)
  ref <- vapply(seq_along(zj_up), function(i) {
    sum(last * stdv[2, 2] / stdv[2, 1] *
          stats::dnorm((zj_up[i] * stdv[2, 2] - zj * stdv[1, 2]) / stdv[2, 1],
                       delta * stdv[2, 1], 1)) * wj_up[i]
  }, numeric(1))
  expect_equal(up, ref, tolerance = 1e-13)

  xq <- 0.4
  p <- tsacor:::rtsa_prob_cpp(xq, last, zj, 2L, stdv, TRUE, delta)
  p_ref <- sum(last * stats::pnorm((xq - zj * stdv[1, 2]) / stdv[2, 1],
                                   delta * stdv[2, 1], 1))
  expect_equal(p, p_ref, tolerance = 1e-13)
})

test_that("compiled functions are identical to RTSA's own first.cpp exports (when RTSA is installed)", {
  skip_if_not_installed("RTSA")
  zj <- seq(-3, 3, length.out = 21); wj <- rep(0.3, 21)
  stdv <- cbind(sqrt(c(0.4, 0.3, 0.3)), sqrt(c(0.4, 0.7, 1.0)))
  delta <- 2.8
  last <- RTSA:::init_int(wj, zj, delta, stdv[, 1])
  expect_equal(tsacor:::rtsa_init_int_cpp(wj, zj, delta, stdv[, 1]), last,
               tolerance = 1e-15)
  zj_up <- seq(-2.5, 2.5, length.out = 15); wj_up <- rep(0.2, 15)
  expect_equal(
    tsacor:::rtsa_recur_int_cpp(2L, stdv, zj, last, zj_up, wj_up, delta, FALSE),
    RTSA:::recur_int(2L, stdv, zj, last, zj_up, wj_up, delta, FALSE),
    tolerance = 1e-15)
  expect_equal(
    tsacor:::rtsa_prob_cpp(0.4, last, zj, 2L, stdv, TRUE, delta),
    RTSA:::prob(0.4, last, zj, 2L, stdv, TRUE, delta),
    tolerance = 1e-15)
})

test_that("design and analysis routes match a live RTSA::boundaries() call (when RTSA is installed)", {
  skip_if_not_installed("RTSA")
  b <- RTSA::boundaries(timing = t8, alpha = 0.05, beta = 0.20, side = 2,
                        futility = "non-binding", es_alpha = "esOF",
                        es_beta = "esOF", type = "design")
  des <- tsacor:::.rtsa_design_bounds(t8, 0.05, 0.20)
  expect_equal(des$root, b$root, tolerance = 1e-7)
  expect_equal(des$alpha_ubound, b$alpha_ubound, tolerance = 1e-7)
  expect_equal(des$beta_ubound, b$beta_ubound, tolerance = 1e-6)

  ret <- tsacor:::.rtsa_retrospective(t8, 0.05, 0.20)
  a <- RTSA::boundaries(timing = ret$analysis$timing, alpha = 0.05, beta = 0.20,
                        side = 2, futility = "non-binding", es_alpha = "esOF",
                        es_beta = "esOF", type = "analysis",
                        design_R = ret$design_R)
  expect_equal(ret$analysis$beta_ubound, a$beta_ubound, tolerance = 1e-6)
})



