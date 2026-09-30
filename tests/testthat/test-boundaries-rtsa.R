## *** LEGACY-ENGINE TESTS (R-only, pre-0.2.7.11) ***
## Everything in this file exercises the OLD pure-R boundary functions
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

test_that("RTSA beta spending formula is used", {
  t <- c(0.1, 0.25, 0.5, 0.75, 1)
  beta <- 0.20
  got <- tsacor:::.rtsa_beta_spend_OF(t, beta)$betaValuesCumulated
  expected <- 2 * stats::pnorm(
    stats::qnorm(1 - beta / 2) / sqrt(t),
    lower.tail = FALSE
  )
  expect_equal(got, expected, tolerance = 1e-14)
  expect_equal(got[length(got)], beta, tolerance = 1e-14)
})

test_that("RTSA retrospective inner-wedge engine has definitive alpha boundary", {
  ## Uses the live RTSA::boundaries() reference alpha boundaries (a
  ## realistic decreasing-then-flattening schedule) rather than a flat
  ## wall: with the 0.2.7.4 root-search now driving the recursion, a
  ## flat wall at every look is an unrealistic, harder-to-bracket case
  ## worth avoiding in a test whose point is the definitive-final-look
  ## property, not root-search robustness (that gets its own tests
  ## below).
  ##
  ## ** 0.2.7.19: ** the final entry of `alpha_ref` below is now the TRUE
  ## final alpha value (computed the same way .rtsa_beta_boundary()
  ## computes it internally), not qnorm(1-alpha/2) -- that was the
  ## 0.2.7.18-and-earlier bug (see the .rtsa_beta_boundary() block
  ## comment), and using the (now strictly smaller) approximation here
  ## would trigger the function's own defensive
  ## `pmin(boundary, c_vec_alpha)` safeguard, masking exactly the
  ## behaviour this test checks. .rtsa_beta_boundary() in fact ignores
  ## this final entry entirely now and recomputes it fresh; this test
  ## checks its result against that same direct recomputation.
  t <- c(0.20, 0.40, 0.60, 0.80, 1)
  alpha <- 0.05
  beta <- 0.20

  true_final_alpha <- utils::tail(tsacor:::.obf_alpha_boundary(t, alpha = alpha), 1L)

  alpha_ref <- c(4.877, 3.357, 2.680, 2.290, true_final_alpha)
  ans <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)

  expect_length(ans$boundary, length(t))
  expect_equal(ans$boundary[length(t)], true_final_alpha, tolerance = 1e-10)
  expect_equal(ans$beta_spent[length(t)], beta, tolerance = 1e-12)
  ## The whole point of the fix: the true value differs materially from
  ## the old qnorm(1-alpha/2) approximation for this schedule.
  expect_true(abs(true_final_alpha - stats::qnorm(1 - alpha / 2)) > 0.01)
})

test_that("RTSA over-powered analysis uses the true (recursion-based) definitive boundary", {
  ## ** 0.2.7.19: ** renamed from "...uses conventional definitive
  ## boundary" -- that name described the pre-fix qnorm(1-alpha/2)
  ## approximation, which is no longer what this function returns.
  t <- c(0.20, 0.40, 0.60, 0.90, 1.10, 1.50)
  alpha <- 0.05
  beta <- 0.20

  true_final_alpha <- utils::tail(
    tsacor:::.obf_alpha_boundary(sort(unique(c(t[t < 1], 1))), alpha = alpha), 1L
  )

  ## Realistic decreasing alpha wall for the four pre-DARIS looks
  ## (approximate OF-type shape); the two over-powered looks (t >= 1) get
  ## the TRUE final alpha value here, matching how tsa_cor.R itself builds
  ## c_vec_alpha for post-DARIS entries in practice (`alpha_bounds_design[
  ## match(pmin(t, 1), boundary_timing)]`, where alpha_bounds_design comes
  ## from the same .obf_alpha_boundary() recursion) -- NOT
  ## qnorm(1-alpha/2): that would now be a strictly SMALLER value than
  ## the true final boundary .rtsa_beta_boundary() computes internally,
  ## and would trigger the function's own defensive
  ## `pmin(boundary, c_vec_alpha)` safeguard, masking exactly the
  ## behaviour this test checks.
  alpha_ref <- c(4.20, 3.10, 2.55, 2.10, true_final_alpha, true_final_alpha)
  ans <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)

  expect_true(ans$over_power)
  expect_equal(ans$fakeIFY, true_final_alpha, tolerance = 1e-10)
  expect_true(all(ans$boundary[t >= 1] == ans$fakeIFY))
  ## Confirms the fix actually changed behaviour for this schedule too
  ## (a conservative, not the maximum measured, threshold: this specific
  ## schedule's true value has not been independently cross-checked
  ## against a live RTSA run the way the c(0.2,0.4,0.6,0.8,1) schedule
  ## in the previous test has been -- only that it is NOT the old
  ## qnorm(1-alpha/2) approximation).
  expect_true(abs(true_final_alpha - stats::qnorm(1 - alpha / 2)) > 1e-6)
})

test_that("RTSA beta boundary returns NA only for negligible early spend, not for every non-positive value", {
  ## ** Updated for the 0.2.7 reconstruction. ** The previous ("inner
  ## wedge") engine forced every non-positive futility value to NA
  ## (`b[b <= 0] <- NA_real_`); real RTSA does not do this -- it returns
  ## NA only where the incremental beta-spend at a look is negligible
  ## enough to be clipped to exactly 0 by the spending-function
  ## tolerance (the "zninf" sentinel case; see boundaries.R's
  ## `abs(lb$za) == 20` checks), OR (0.2.7.4) where the two-pass
  ## rm_bs suppression zeroes an early look's spend because it came
  ## back negative on the unconstrained first pass. A finite, even
  ## negative, non-binding futility boundary is otherwise a legitimate
  ## result and must NOT be suppressed just for being non-positive.
  ##
  ## An extremely early first look (t = 1e-6) is a case we can reason
  ## about in closed form regardless of the recursive integration's
  ## numerics: the incremental spend there underflows to exactly 0 in
  ## .rtsa_beta_spend_OF() well before any recursion or root search
  ## runs, so the NA rule must fire at that look specifically under
  ## either mechanism.
  ## ** FIXED in 0.2.7.20 (not an R CMD check failure -- masked by the
  ## function's own defensive pmin(boundary, c_vec_alpha) clip, but
  ## testing the wrong thing for the wrong reason). ** alpha_ref's last
  ## entry used to be the stale qnorm(1-alpha/2) constant; since that is
  ## SMALLER than the true recomputed final_alpha_bound (~2.031, see the
  ## "RTSA retrospective inner-wedge engine" test above), the function's
  ## own pmin() safeguard clipped the true value back down to 1.96 before
  ## it reached this test's assertion -- so `expect_equal(...,
  ## qnorm(1-alpha/2))` was passing by coincidentally re-deriving the OLD
  ## wrong constant via the clip, not by actually exercising the fixed
  ## recomputation. Using the true value here (as the two tests above
  ## already do) means alpha_ref no longer clips anything, and this
  ## assertion tests what it claims to.
  t <- c(1e-6, 0.20, 0.40, 0.60, 0.80, 1.0)
  alpha <- 0.05
  beta <- 0.20
  true_final_alpha <- utils::tail(tsacor:::.obf_alpha_boundary(t, alpha = alpha), 1L)
  alpha_ref <- c(6.0, 4.877, 3.357, 2.680, 2.290, true_final_alpha)
  spend <- tsacor:::.rtsa_beta_spend_OF(t, beta)

  ans <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)

  expect_equal(spend$betaValuesDelta[1], 0)
  expect_true(is.na(ans$boundary[1]))
  expect_equal(tail(ans$boundary, 1), true_final_alpha, tolerance = 1e-10)

  ## Deliberately NOT asserted: that every finite boundary value is
  ## positive. Forcing that was the bug.
})

test_that("futility boundary is capped at the corresponding alpha boundary", {
  ## This remains a hard guarantee: not because RTSA's own recursion
  ## enforces it at every step (it doesn't -- see the block comment in
  ## R/obf_boundaries.R), but because tsacor applies an explicit,
  ## disclosed pmin() safeguard against the supplied alpha boundary on
  ## top of the reconstructed RTSA recursion.
  t <- c(0.10, 0.25, 0.50, 0.75, 1.00)
  alpha <- 0.05
  beta <- 0.20
  alpha_ref <- c(5.0, 3.5, 2.8, 2.2, qnorm(1 - alpha / 2))
  ans <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)

  ok <- !is.na(ans$boundary)
  expect_true(all(ans$boundary[ok] <= alpha_ref[ok] + 1e-14))
})

test_that("RTSA beta engine uses the fixed theoretical drift, not an empirical one", {
  ## ** Core regression test for the 0.2.7 reconstruction. ** The
  ## previous engine derived its drift empirically from the shape of a
  ## symmetric, null-referenced futility wedge ("testDrift"). RTSA's own
  ## beta_boundary() uses a FIXED, closed-form drift computed directly
  ## from alpha, beta, and side -- delta = |qnorm(alpha/side) +
  ## qnorm(beta)| -- independent of the data or of the wedge shape.
  alpha <- 0.05
  beta <- 0.20
  t <- c(0.20, 0.40, 0.60, 0.80, 1.00)
  ## The live RTSA::boundaries() reference alpha boundaries for this
  ## schedule (see test-alpha-spend-rtsa.R / NEWS.md 0.2.6).
  alpha_ref <- c(4.877, 3.357, 2.680, 2.290, 2.031)
  ans <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)

  expect_equal(ans$delta, abs(qnorm(alpha / 2) + qnorm(beta)), tolerance = 1e-14)
})

test_that("0.2.7.4: the information-scale root search actually engages and hits its target", {
  ## Core regression test for the 0.2.7.4 fix: .rtsa2_find_warp_root()
  ## must find an x such that .rtsa2_inf_warp(x, ...) is (numerically)
  ## zero, i.e. the final look's futility boundary, computed under that
  ## warp, lands on the fixed efficacy boundary there.
  ##
  ## ** FIXED in 0.2.7.20. ** alpha_ref's final entry used to be hardcoded
  ## to stats::qnorm(1 - alpha/2) -- exactly the wrong "continuous-limit"
  ## constant .rtsa_beta_boundary() itself was fixed (0.2.7.19) to stop
  ## assuming. Once that fix landed, `ans <- .rtsa_beta_boundary(...)`
  ## started recomputing its OWN, correct final_alpha_bound internally
  ## (via .obf_alpha_boundary(), see obf_boundaries.R) rather than trusting
  ## this test's alpha_ref -- so ans$warp_root solves a DIFFERENT final
  ## look than the one gap_ans below was checking it against, and the two
  ## stopped agreeing (off by 0.072). Fixed by deriving alpha_ref from the
  ## same .obf_alpha_boundary() recursion .rtsa_beta_boundary() itself now
  ## uses, so both sides of this check are testing the same equation.
  alpha <- 0.05
  beta <- 0.20
  t <- c(0.20, 0.40, 0.60, 0.80, 1.00)
  alpha_ref <- tsacor:::.obf_alpha_boundary(t, alpha = alpha)
  delta <- abs(qnorm(alpha / 2) + qnorm(beta))

  root <- tsacor:::.rtsa2_find_warp_root(
    t = t, beta = beta, delta = delta, alpha_ubound = alpha_ref,
    rm_bs = 0L, start = 0.95, step = 0.02, max_iter = 50L
  )
  gap <- tsacor:::.rtsa2_inf_warp(root, t = t, beta = beta, delta = delta,
                                  alpha_ubound = alpha_ref, rm_bs = 0L)

  expect_true(is.finite(root))
  expect_equal(gap, 0, tolerance = 1e-6)

  ## And the full public wrapper's own root search converges to a root
  ## of the SAME objective function AT ITS OWN rm_bs (not necessarily
  ## rm_bs = 0: PASS 2 of .rtsa_beta_boundary()'s pipeline suppresses
  ## whichever early looks came back negative on PASS 1, which changes
  ## the objective function itself, so ans$warp_root and the rm_bs = 0
  ## `root` above are, correctly, roots of two DIFFERENT equations and
  ## are not expected to coincide -- asserting near-equality between
  ## them here was the actual test bug (found by R CMD check on
  ## 0.2.7.5), not a package bug: ans$warp_root and `root` differed by
  ## about 4.6e-4, larger than this test's original 1e-4 tolerance,
  ## simply because this specific design's PASS 1 suppressed one early
  ## look (ans$rm_bs > 0) while the isolated `root` above was computed
  ## at rm_bs = 0 throughout. The correct invariant -- and the one this
  ## test now checks -- is that ans$warp_root actually solves the
  ## pipeline's own (rm_bs-adjusted) equation, not that it matches a
  ## different equation's root.
  ##
  ## `alpha_ref` is passed here only as .rtsa_beta_boundary()'s c_vec_alpha
  ## argument (used for the function's own defensive pmin() clip, NOT for
  ## its internal final_alpha_bound -- see obf_boundaries.R); gap_ans below
  ## is checked against the SAME .obf_alpha_boundary()-derived alpha_ref
  ## the function recomputes for itself internally, which is why they
  ## still have to agree regardless of what c_vec_alpha the caller passed.
  ans <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)
  gap_ans <- tsacor:::.rtsa2_inf_warp(
    ans$warp_root, t = t, beta = beta, delta = delta,
    alpha_ubound = alpha_ref, rm_bs = ans$rm_bs
  )
  expect_equal(gap_ans, 0, tolerance = 1e-6)
})

test_that("0.2.7.4: the naive warp_root = 1 case is measurably different from the root-found solution", {
  ## Confirms the root search is not a no-op: evaluating the SAME
  ## design at warp_root = 1 (what 0.2.7-0.2.7.3 implicitly always did)
  ## gives a materially different final-look gap than the root-found
  ## value does at (numerically) zero -- i.e. there is a real
  ## discrepancy for the root search to close, matching the reported
  ## "beta-bounds are far from RTSA" symptom.
  alpha <- 0.05
  beta <- 0.20
  t <- c(0.20, 0.40, 0.60, 0.80, 1.00)
  alpha_ref <- c(4.877, 3.357, 2.680, 2.290, stats::qnorm(1 - alpha / 2))
  delta <- abs(qnorm(alpha / 2) + qnorm(beta))

  gap_naive <- tsacor:::.rtsa2_inf_warp(1, t = t, beta = beta, delta = delta,
                                        alpha_ubound = alpha_ref, rm_bs = 0L)
  expect_gt(abs(gap_naive), 1e-3)
})

test_that("0.2.7.5: .rtsa_beta_spend_OF() accepts t = 0 (needed by rm_bs suppression)", {
  ## Regression test for the near-universal R CMD check failure fixed
  ## in 0.2.7.5: PASS 2 of .rtsa_beta_boundary()'s root search (any
  ## rm_bs > 0) zeroes the first `rm_bs` entries of its timing vector
  ## by design, and .rtsa_beta_spend_OF() previously rejected any t
  ## <= 0 outright, which made every such call fail with this
  ## function's own error -- not a genuine root-search non-convergence,
  ## even though it surfaced as "root search did not converge" one
  ## level up. The underlying formula was always well-defined at t = 0
  ## (qnorm(1-beta/2)/sqrt(0) = +Inf, pnorm(Inf, lower.tail=FALSE) = 0,
  ## so cumulative spend is exactly 0, the correct "no spend yet"
  ## value) -- only the validation guard was wrong.
  beta <- 0.20
  out <- tsacor:::.rtsa_beta_spend_OF(c(0, 0, 0.3, 0.6, 1), beta)
  expect_equal(out$betaValuesCumulated[1:2], c(0, 0))
  expect_equal(out$betaValuesDelta[1:2], c(0, 0))
  expect_true(all(is.finite(out$betaValuesCumulated)))

  ## A genuinely negative information fraction must still be rejected.
  expect_error(tsacor:::.rtsa_beta_spend_OF(c(-0.1, 0.5, 1), beta))
})

test_that("0.2.7.5: the two-pass root search converges on a realistic multi-look design", {
  ## End-to-end regression test for the same bug, exercised through the
  ## full .rtsa_beta_boundary() pipeline (PASS 1 -> rm_bs -> PASS 2)
  ## rather than the isolated spend function, on a schedule chosen to
  ## produce rm_bs > 0 on the first pass (early, tightly-spaced looks),
  ## which is exactly the condition that triggered the bug.
  alpha <- 0.05
  beta <- 0.20
  t <- c(0.01, 0.02, 0.03, 0.05, 0.08, 0.90, 1.00)
  ## ** FIXED in 0.2.7.20 (not an R CMD check failure -- masked by pmin(),
  ## same as the fix immediately above). ** alpha_ref's last entry is now
  ## the true recomputed final alpha value rather than the stale
  ## qnorm(1-alpha/2) constant, for the same reason: the old constant is
  ## smaller than the truth, so the function's own defensive clip was
  ## silently reproducing it regardless of whether the recomputation
  ## fix was working.
  true_final_alpha <- utils::tail(tsacor:::.obf_alpha_boundary(t, alpha = alpha), 1L)
  ## A plausible decreasing-then-flattening alpha wall for the interior
  ## looks (not RTSA-verified numbers, just realistic shape/magnitude for
  ## this test) with the TRUE final value at the end.
  alpha_ref <- c(9.5, 7.3, 6.2, 5.1, 4.2, 2.10, true_final_alpha)

  ans <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)
  expect_length(ans$boundary, length(t))
  expect_true(ans$rm_bs > 0)
  expect_equal(ans$boundary[length(t)], true_final_alpha, tolerance = 1e-10)
  ## Suppressed early looks must come back NA (routed through the
  ## zninf sentinel by the rm_bs mechanism), not error out.
  expect_true(any(is.na(ans$boundary[seq_len(5)])))
})

test_that("0.2.7.5: a single-look schedule (DARIS already reached at the first study) does not error", {
  ## Regression test for the second bug found alongside the t = 0 fix:
  ## when every observed information fraction is already >= 1 (no
  ## pre-DARIS interim look at all), timing_beta collapses to the
  ## single synthetic point t = 1, where the incremental beta-spend is
  ## exactly `beta` by construction regardless of warp_root -- so
  ## .rtsa2_inf_warp() is a constant function of the warp factor and
  ## can never bracket a root. .rtsa_beta_boundary() now detects this
  ## and skips the root search entirely rather than exhausting all 50
  ## widening attempts and erroring.
  alpha <- 0.05
  beta <- 0.20
  t <- c(1.2, 1.5, 2.0)  ## every study individually exceeds DARIS
  alpha_ref <- rep(stats::qnorm(1 - alpha / 2), length(t))

  ans <- tsacor:::.rtsa_beta_boundary(t, alpha, beta, alpha_ref)
  expect_length(ans$boundary, length(t))
  expect_true(all(ans$boundary == stats::qnorm(1 - alpha / 2)))
  expect_true(ans$over_power)
})

test_that("legacy (non-RTSA-matching) beta engine internals have been removed", {
  ## The 0.2.4-0.2.6.9 "inner wedge" engine (and its own cited sources,
  ## none of which exist in actual RTSA 0.2.2) is fully superseded by
  ## the 0.2.7 reconstruction; its now-dead helper functions should not
  ## still be present.
  removed <- c(".rtsa_gfunc", ".rtsa_trap", ".rtsa_fcab", ".rtsa_qpos",
               ".rtsa_first", ".rtsa_other", ".rtsa_searchfunc_old",
               ".rtsa_get_inner_wedge")
  ns <- asNamespace("tsacor")
  present <- removed[vapply(removed, exists, logical(1), envir = ns,
                             inherits = FALSE)]
  expect_length(present, 0)
})


test_that("0.2.7.1: .rtsa2_z_n_w() tolerates collapsed za/zb without erroring", {
  ## Direct unit test of the safety net in .rtsa2_z_n_w(): even if two
  ## boundaries are passed in already equal or reversed (bypassing the
  ## .rtsa2_beta_boundary_core() clamp entirely, e.g. a future direct
  ## caller), the function must return a valid grid rather than erroring.
  info <- list(sd_incr = c(1, 1), sd_proc = c(1, 1))

  ## Equal za/zb at look 2 (would previously collapse the trimmed grid).
  out_equal <- tsacor:::.rtsa2_z_n_w(
    r = 18, info = info, za = c(0, 1.5), zb = c(3, 1.5), i = 2L, delta = 0
  )
  expect_true(all(is.finite(out_equal$zj)))
  expect_true(all(is.finite(out_equal$wj)))
  expect_true(length(out_equal$zj) >= 3L)

  ## Reversed za/zb at look 2 (za > zb).
  out_reversed <- tsacor:::.rtsa2_z_n_w(
    r = 18, info = info, za = c(0, 1.6), zb = c(3, 1.5), i = 2L, delta = 0
  )
  expect_true(all(is.finite(out_reversed$zj)))
  expect_true(all(is.finite(out_reversed$wj)))
  expect_true(length(out_reversed$zj) >= 3L)
})

test_that("0.2.7.1: the za/zb clamp is a no-op for well-separated boundaries", {
  ## The preventive clamp must not perturb any normal, well-separated
  ## boundary pair -- it should only ever bite within its tiny margin.
  expect_equal(tsacor:::.rtsa2_clamp_za_below_zb(0, 3), 0)
  expect_equal(tsacor:::.rtsa2_clamp_za_below_zb(-20, 2.5), -20)
  ## But it clamps when the two are within (or past) the margin.
  expect_lt(tsacor:::.rtsa2_clamp_za_below_zb(1.9999995, 2), 2)
  expect_lt(tsacor:::.rtsa2_clamp_za_below_zb(2.5, 2), 2)
})

test_that("0.2.7.2: .rtsa2_seq_by2() returns empty rather than erroring when to < from", {
  ## Direct unit test for the second seq()-by-argument crash, fixed in
  ## 0.2.7.2: raw seq(from, to, 2) throws "wrong sign in 'by' argument"
  ## when to < from, instead of returning an empty index range.
  expect_identical(tsacor:::.rtsa2_seq_by2(3, 1), integer(0))
  expect_identical(tsacor:::.rtsa2_seq_by2(2, 0), integer(0))
  ## Unaffected (non-empty-range) cases are untouched.
  expect_identical(tsacor:::.rtsa2_seq_by2(3, 3), 3)
  expect_identical(tsacor:::.rtsa2_seq_by2(2, 4), c(2, 4))
})

test_that("0.2.7.2: a near-collapsed [za, zb] window produces a valid, non-erroring Simpson grid", {
  ## Regression test for the exact crash in the R CMD check log: a very
  ## narrow [za[i], zb[i]] window (as the 0.2.7.1 clamp, or a naturally
  ## tight boundary pair, can produce) previously errored inside the
  ## Simpson-weight loop at seq(3, m - 2, 2) = seq(3, 1, 2), once the
  ## 0.2.7.1 fix made a 2-node `xi` (m = 3) a reachable case.
  ##
  ## This deliberately does NOT hard-code the exact number of surviving
  ## grid nodes: an earlier version of this test used lo = 1.5, which
  ## (for r = 18, delta = 0, sd_incr = 1) happens to land EXACTLY on one
  ## of .rtsa2_z_n_w()'s own log-spaced grid nodes (node j = 72 evaluates
  ## to precisely -3 + (72 - 18)/12 = 1.5), so one extra node coincided
  ## with `lo` and survived trimming, giving m = 5 instead of the m = 3
  ## the test assumed -- a correct result from the code, just an
  ## incorrect prediction in the test. `lo` below is chosen off any such
  ## grid-node alignment, but more importantly this test now checks the
  ## invariants that must hold for ANY valid Simpson grid built on
  ## [lo, hi], regardless of exactly how many interior nodes (zero, one,
  ## or by rarer coincidence more) survive trimming: no error, an odd
  ## node count >= 3, endpoints exactly at lo/hi, a non-decreasing
  ## sequence, and -- the property Simpson's rule actually guarantees,
  ## and what the downstream recursion actually relies on -- total
  ## weight equal to the window width.
  info <- list(sd_incr = c(1, 1), sd_proc = c(1, 1))
  lo <- 1.234567  ## deliberately off the r = 18 grid's node spacing
  hi <- lo + 1e-6  ## mimics the 0.2.7.1 clamp margin: no interior node survives
  out <- tsacor:::.rtsa2_z_n_w(
    r = 18, info = info, za = c(0, lo), zb = c(3, hi), i = 2L, delta = 0
  )
  m <- length(out$zj)

  expect_true(m >= 3L && m %% 2L == 1L)
  expect_length(out$wj, m)
  expect_true(all(is.finite(out$zj)))
  expect_true(all(is.finite(out$wj)))
  expect_equal(out$zj[1], lo, tolerance = 1e-12)
  expect_equal(out$zj[m], hi, tolerance = 1e-12)
  expect_true(all(diff(out$zj) >= -1e-12))  ## non-decreasing
  expect_equal(sum(out$wj), hi - lo, tolerance = 1e-9)
})
