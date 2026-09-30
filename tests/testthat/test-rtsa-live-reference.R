## Frozen LIVE-RTSA reference (inst/extdata/rtsa_0.2.2_reference.R): numbers
## printed by RTSA 0.2.2 itself, so these tests do NOT depend on RTSA being
## installed and do NOT go through the Python reconstruction. Until 0.2.7.17
## the package's numerical parity was established against a Python port written
## from RTSA's sources; this file compares the compiled engine with RTSA's own
## output. (Alpha bounds are deterministic recursions, compared tightly. Design
## beta bounds sit on a uniroot() result with tol = 1e-9 and RTSA's own
## uniroot() iterates differ from ours at rounding level, so they and the root
## use 1e-7. The analysis route has no root search, so its beta bounds are
## compared at 1e-9 -- observed agreement is ~1e-15.)

ref <- dget(system.file("extdata", "rtsa_0.2.2_reference.R", package = "tsacor"))

test_that("the frozen reference is what it claims to be", {
  expect_identical(ref$rtsa_version, "0.2.2")
  expect_equal(ref$design$timing, c(0.25, 0.50, 0.75, 1.00))
  expect_length(ref$design$beta_ubound, 4L)
  ## the analysis call's timing stopped short of design_R -> 5 beta elements
  expect_length(ref$analysis$beta_ubound, 5L)
  expect_equal(ref$analysis$design_R, ref$design$root)
})

test_that("design route reproduces live RTSA 0.2.2 (root, alpha, beta, spending)", {
  des <- suppressWarnings(
    tsacor:::.rtsa_design_bounds(ref$design$timing, ref$alpha, ref$beta))
  expect_equal(des$root, ref$design$root, tolerance = 1e-7)
  expect_equal(des$alpha_ubound, ref$design$alpha_ubound, tolerance = 1e-10)
  expect_equal(des$beta_ubound, ref$design$beta_ubound, tolerance = 1e-7)
  expect_equal(des$rm_bs, 1L)            # RTSA suppressed exactly the first look
  expect_equal(des$beta_spent, ref$design$bs_cum, tolerance = 1e-9)
  expect_equal(des$beta_spent_delta, ref$design$bs_incr, tolerance = 1e-9)
  ## and the calibration's defining property, as RTSA prints it
  n <- length(des$alpha_ubound)
  expect_equal(des$beta_ubound[n], des$alpha_ubound[n], tolerance = 1e-7)
})

test_that("alpha engine reproduces live RTSA 0.2.2 exactly (design and analysis spending)", {
  a_des <- tsacor:::.rtsa_alpha_cpp(ref$design$timing, side = 2L, alpha = ref$alpha)
  expect_equal(a_des$alpha_ubound, ref$design$alpha_ubound, tolerance = 1e-12)
  ## analysis-route alpha differs from design-route alpha at the ~1e-8 level
  ## (information scale t * design_R changes the grid), and RTSA reports that
  a_ana <- tsacor:::.rtsa_alpha_cpp(ref$analysis$timing, side = 2L, alpha = ref$alpha,
                                   design_R = ref$analysis$design_R)
  expect_equal(a_ana$alpha_ubound, ref$analysis$alpha_ubound, tolerance = 1e-12)
  expect_gt(max(abs(a_ana$alpha_ubound - a_des$alpha_ubound)), 1e-9)
})

test_that("analysis route reproduces live RTSA 0.2.2 (RTSA's own call shape)", {
  ana <- suppressWarnings(tsacor:::.rtsa_analysis_bounds(
    ref$analysis$timing, ref$analysis$design_R, ref$alpha, ref$beta))
  expect_equal(ana$alpha_ubound, ref$analysis$alpha_ubound, tolerance = 1e-12)
  expect_length(ana$beta_ubound, 5L)
  expect_equal(ana$beta_ubound, ref$analysis$beta_ubound, tolerance = 1e-9)
  expect_equal(ana$beta_spent, ref$analysis$bs_cum, tolerance = 1e-9)
  expect_equal(ana$beta_spent_delta, ref$analysis$bs_incr, tolerance = 1e-9)
  ## RTSA's final clamp compares beta_ubound[length(alpha_ubound)], i.e. look 4
  ## here, and look 4 (1.72) is below its wall (2.01): no clamp; the appended
  ## look's 2.1228 (above 2.01) is deliberately left unclamped, as RTSA prints it
  expect_gt(ana$beta_ubound[5], ana$alpha_ubound[4])
})

test_that("a length mismatch other than RTSA's appended-look shape is still rejected", {
  ub <- tsacor:::.rtsa_alpha_cpp(c(0.25, 0.5, 0.75, 1), side = 2L, alpha = 0.05)$alpha_ubound
  delta <- abs(stats::qnorm(0.025) + stats::qnorm(0.20))
  ## design mode (warp_root): alpha bound must match the timeline exactly
  expect_error(
    tsacor:::.rtsa_beta_cpp(c(0.25, 0.5, 0.75, 1, 1.1), ub, 0.2, delta,
                           warp_root = 1.1, warn = FALSE),
    "equal length")
  ## analysis mode: one short is RTSA's shape; two short is not
  expect_error(
    tsacor:::.rtsa_beta_cpp(c(0.25, 0.5, 0.75, 1), ub[1:2], 0.2, delta,
                           design_R = 1.13, warn = FALSE),
    "equal length")
})

test_that("LEGACY fallback (0.2.7.19 fix) is close to live RTSA and no longer sits on the qnorm(1 - alpha/2) wall", {
  ## The legacy R-only engine uses the approximate FFT alpha recursion, so it
  ## cannot reproduce RTSA to more than ~1e-3; what this test protects is the
  ## 0.2.7.19 fix. Measured in R (0.2.7.21) on this schedule: FFT alpha
  ## 4.332634 2.963388 2.358980 2.012955 (final-look error -1.14e-3); legacy
  ## warp root 1.132483 vs live RTSA 1.133242 (error 7.6e-4). Before the fix the
  ## root was 1.0972 (error 3.6e-2) and the final wall 1.95996, so the
  ## tolerances below (2e-3 / 3e-3) separate fixed from unfixed.
  t4 <- ref$design$timing
  a_fft <- tsacor:::.obf_alpha_boundary(t4, ref$alpha)
  expect_lt(max(abs(a_fft - ref$design$alpha_ubound)), 3e-3)

  leg <- tsacor:::.rtsa_beta_boundary(t4, ref$alpha, ref$beta, a_fft)
  expect_lt(abs(leg$warp_root - ref$design$root), 2e-3)

  ok <- !is.na(ref$design$beta_ubound)
  expect_false(anyNA(leg$boundary[ok]))
  expect_lt(max(abs(leg$boundary[ok] - ref$design$beta_ubound[ok])), 3e-3)
  ## the final wall is the recomputed one (~2.013), not qnorm(1 - alpha/2) (1.960)
  expect_gt(utils::tail(leg$boundary, 1) - stats::qnorm(1 - ref$alpha / 2), 0.03)
})

