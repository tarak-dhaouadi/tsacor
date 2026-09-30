## Copyright (C) the RTSA authors (Anne Lyngholm Soerensen, Markus Harboe Olsen,
## Theis Lange, Christian Gluud) for the algorithms and code this file is derived
## from (RTSA 0.2.2, GPL (>= 2)); copyright (C) Tarak Dhaouadi for the
## adaptation. This file is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by the Free
## Software Foundation; either version 2 of the License, or (at your option) any
## later version. See DESCRIPTION and inst/COPYRIGHTS.
## -------------------------------------------------------------------------
## NOTE (tsacor 0.1.0): this file is inherited from the tsahr package (same
## maintainer), which in turn ports RTSA. Version numbers such as "0.2.7.13" in
## the comments below refer to the development history of tsahr, NOT of tsacor,
## and comments that name tsa_cor() / "tsacor" describe what were originally
## tsa_hr() / tsahr statements in that history (the boundary engine itself is
## effect-measure agnostic: it works on information fractions only). The only
## functional changes made when tsacor was created are the package/registration
## names and the rename of the events-based helpers to their participants-based
## tsacor counterparts (see .tsacor_participants_at_fraction()).
## -------------------------------------------------------------------------
## Internal O'Brien-Fleming-type alpha-spending and non-binding
## beta-spending group sequential boundary calculator.
##
## This replaces an earlier dependency on 'rpact' for computing trial
## sequential monitoring boundaries. rpact's getDesignGroupSequential()
## imposes a hard, software-specific limit of kMax <= 20 looks (with only
## an advisory warning above kMax=10) -- a package validation choice, not
## a statistical one. A meta-analysis with many included studies (each
## study being one "look" in TSA) routinely exceeds this. The functions
## below implement the same underlying, textbook methodology directly, so
## the number of looks is limited only by computation time, not an
## artificial cap.
##
## Methodology:
##   - Alpha-spending function (Lan & DeMets 1983 O'Brien-Fleming-type
##     approximation, RTSA/CTU parameterisation, TWO-SIDED, side = 2 --
##     see VALIDATION below):
##       alpha*(t) = 4*(1 - Phi(z_{alpha/4}/sqrt(t)))
##     ** CORRECTED in 0.2.6. ** Versions 0.2.0-0.2.5.1 used
##     alpha*(t) = 2*(1 - Phi(z_{alpha/2}/sqrt(t))) -- the more
##     commonly-seen textbook two-argument form, quoted directly in
##     gsDesign's documentation and in several published trial
##     statistical analysis plans for "two-sided alpha spending". That
##     form is internally valid as *an* alpha-spending function (it does
##     reach exactly `alpha` at t=1), and would be the right choice in
##     those other contexts, but it does NOT match the specific TSA
##     methodology (Copenhagen Trial Unit / RTSA / Thorlund et al.) that
##     this package models itself on and that its own documentation
##     cites (Miladinovic et al. 2013, Wetterslev et al. 2009). See
##     VALIDATION below for the live RTSA::boundaries() comparison that
##     confirms the corrected formula and quantifies the error in the
##     old one.
##   - Beta-spending function (non-binding futility, O'Brien-Fleming-type,
##     targeting a central inner wedge rather than two outer tails, and
##     matching RTSA's own getInnerWedge()/beta-spending construction --
##     see inst/REVERSE_ENGINEERING_RTSA.md for the reverse-engineering
##     write-up, and the "RTSA / original CTU-TSA non-binding futility
##     engine" section below for the implementation):
##       beta*(t) = 2*(1 - Phi(z_{beta/2}/sqrt(t)))
##     implemented in .rtsa_beta_spend_OF() below. (An earlier,
##     un-doubled candidate formula, beta*(t) = 1 - Phi(z_beta/sqrt(t)),
##     was implemented in this file as .beta_spend_OF() but never used by
##     the actual boundary-solving code path above -- it was dead code
##     and has been removed.)
##   - Sequential boundaries at each information fraction are obtained via
##     the standard recursive numerical integration (Armitage, McPherson &
##     Rowe 1969; generalised for spending functions by Reboussin, DeMets,
##     Kim & Lan 2000), carried out here in "B-value" (Brownian motion)
##     space, where increments between looks are simple additive Gaussians
##     -- computed efficiently via FFT-based convolution (stats::convolve)
##     rather than an O(K*N^2) dense matrix, so this comfortably scales to
##     many looks (tested to K > 100 in well under a second per call).
##   - Non-binding futility (beta) boundaries are computed holding the
##     already-fixed alpha boundary, under the alternative hypothesis with
##     standardised drift theta = z_{alpha/2} + z_beta (matching the
##     single-look power calculation used elsewhere in this package), as
##     is standard for non-binding futility monitoring.
##
## VALIDATION:
##   - Alpha-spending boundary engine, FORMULA correctness (0.2.6 fix):
##     confirmed against a live call to RTSA::boundaries() -- the
##     reference TSA implementation this package documents itself as
##     following -- across the full 5-look schedule, not just one point:
##       > RTSA::boundaries(timing=c(0.2,0.4,0.6,0.8,1), alpha=0.05,
##                           side=2, es_alpha="esOF")
##       Upper: 4.877  3.357  2.680  2.290  2.031
##     The corrected formula alpha*(t) = 4*(1-Phi(z_{alpha/4}/sqrt(t)))
##     is clearly distinguishable from the pre-correction
##     2*(1-Phi(z_{alpha/2}/sqrt(t))) form, which gives 4.383, 3.099,
##     2.554, 2.254, 2.063 for the same schedule -- measurably
##     different, especially at the early looks, where it matters most.
##     Both forms independently satisfy alpha*(1) = alpha, which is why
##     this bug was not caught by the pre-0.2.6 Monte Carlo/closed-form
##     checks below -- those confirm the *total* spend is correct, not
##     the *shape* across interim looks. A first-look regression test
##     (`test-alpha-spend-rtsa.R`) pins the corrected formula against
##     this RTSA reference in exact closed form (the first look has no
##     prior boundary to condition on, so needs no recursive-engine
##     approximation).
##   - Alpha-spending boundary engine, full 5-LOOK RECURSIVE-ENGINE match
##     against the same RTSA reference: the full schedule above was
##     independently re-implemented in Python (a faithful line-by-line
##     port of `.convolve_step()`/`.obf_alpha_boundary()`'s math -- R's
##     `convolve(h, rev(kernel), type="open")` is, by construction of
##     that call, a standard convolution sum, which the Python port
##     computes directly) and run end-to-end (all 5 looks through the
##     recursion), not just checked in closed form at the first look.
##     IMPORTANT -- use the nominal timing actually passed to
##     RTSA::boundaries() (0.2, 0.4, 0.6, 0.8, 1.0), not the rounded
##     `SMA_Timing` column RTSA *reports back* (0.205, 0.409, ...): the
##     two were confused at one point during this package's development
##     and produced a spurious ~0.06 "discrepancy" at the first look
##     that was really just a units mismatch, not a real formula or
##     engine problem -- using the correct nominal timing, the first
##     look matches RTSA to within 0.0001, as it must (it's an exact
##     closed-form quantity, independent of the recursive engine or its
##     grid resolution entirely).
##     Using the correct nominal timing, the Python port's result at
##     n_grid=2000 (this engine's long-standing default) was already
##     close to RTSA across all 5 looks (max absolute error ~0.003-0.006
##     in the runs checked), NOT the ~0.03 final-look error that was
##     reported in an earlier draft of this note. That earlier ~0.03
##     figure, and the accompanying "error shrinks from 0.0325 at
##     n_grid=2000 down to 0.0020 at n_grid=32000" progression, could
##     NOT be reproduced by the independent Python port and are not
##     currently believed to be accurate -- they were not confirmed by
##     an actual run of this package's R code (no R interpreter was
##     available in the environment that wrote that draft) and should be
##     treated as unverified until someone with R access runs the
##     snippet below and reports the real numbers:
##       t <- c(0.2, 0.4, 0.6, 0.8, 1.0); alpha <- 0.05
##       tsacor:::.obf_alpha_boundary(t, alpha, n_grid = 2000)
##       tsacor:::.obf_alpha_boundary(t, alpha, n_grid = 16000)
##     Pending that confirmation, `n_grid`'s default is nonetheless kept
##     at 16000 (raised from 2000) purely as a conservative, low-cost
##     safety margin: FFT-based convolution makes the accuracy gain
##     (finer discretisation of the recursive integral) essentially free
##     in practice for this package's actual usage pattern (a handful of
##     interim looks per `tsa_cor()` call, not a hot simulation loop), so
##     there is no real reason to prefer the smaller grid even without a
##     confirmed problem at 2000. A regression test
##     (`test-alpha-spend-rtsa.R`) pins the full recursive-engine output
##     against all 5 RTSA reference values at a 0.01 tolerance using the
##     correct nominal timing; this is expected to pass at both
##     n_grid=2000 and n_grid=16000 based on the Python check above, but
##     has not been confirmed against the actual R implementation.
##     ** Measured in R (0.2.7.21, at the default n_grid = 16000): **
##       tsacor:::.obf_alpha_boundary(c(0.25, 0.5, 0.75, 1), 0.05)
##         = 4.332634 2.963388 2.358980 2.012955
##       live RTSA 0.2.2 (inst/extdata/rtsa_0.2.2_reference.R)
##         = 4.332634 2.963131 2.359044 2.014090
##       error (FFT - RTSA): 4e-7, 2.6e-4, -6.4e-5, -1.14e-3 -- largest at the
##       final look. One schedule only: not a general accuracy bound. The
##       5-look snippet above (t = 0.2 ... 1) has still not been run.
##   - The recursive integration engine itself (the FFT-based recursion
##     below) was independently reproduced line-by-line in Python (scipy)
##     and checked two ways. These checks validate the recursion
##     MECHANICS -- i.e. that the engine correctly turns a given spending
##     function into boundaries that spend exactly that much alpha
##     overall -- independent of which spending-function formula is
##     plugged in; they do NOT by themselves validate which formula
##     should be plugged in, since (as the fix above found) more than one
##     formula can each self-consistently reach the target alpha at t=1
##     while implying materially different boundaries at earlier looks.
##     That formula-choice question was instead checked directly against
##     live RTSA::boundaries() output (above), not by Monte Carlo:
##       1. The K=5 equally-spaced two-sided alpha=0.05 case, run under
##          the PRE-0.2.6 spending function, matches the classical,
##          widely-published O'Brien-Fleming boundary constant (~2.040
##          final-look value) to within grid-discretisation error. The
##          formula-choice question itself (old vs. RTSA-matched) was
##          checked separately and directly against live
##          RTSA::boundaries() output, both in closed form at the first
##          look and, as of 0.2.6.1, end-to-end through the full
##          recursive engine at all 5 looks (see the VALIDATION entries
##          above for both).
##       2. A direct (non-simulation) 2D numerical integration for K=2
##          confirmed the constructed boundaries yield the intended
##          overall two-sided alpha almost exactly (0.0506 vs a 0.05
##          target).
##       3. A 100,000-replicate Monte Carlo across five configurations
##          (K=2, K=3 unequally spaced, K=5, K=10, K=39 unequally spaced)
##          gave empirical type-I error of 5.02%-5.22% against a 0.05
##          nominal target (Monte Carlo SE ~=0.07%) -- i.e. within 2-3
##          simulation SEs in every configuration checked, and with no
##          degradation at larger K. This was run under the pre-0.2.6
##          formula. The recursion machinery itself was NOT changed by
##          this fix (only the four-line .alpha_spend_OF() body was), so
##          this remains valid evidence for the recursion -- and a
##          smaller-scale re-check (20,000 replicates, K=2/3/5/10) was
##          run directly against the corrected, RTSA-matched formula
##          specifically to confirm the fix didn't break exact-alpha
##          control: empirical type-I error came back at 4.35%-4.93%
##          against the 5% nominal target (MC SE ~=0.15%), i.e. within
##          about 0.5-4 simulation SEs, consistent with correct
##          behaviour.
##       This engine-level validation was carried out without R (no
##       rpact comparison was possible in that environment); it is
##       independent of this package's own R test suite, but is not a
##       comparison against rpact specifically. Users who need a
##       like-for-like comparison against rpact should still do so
##       themselves for their exact design if that specific software
##       agreement matters to them.
##   - Non-binding beta/futility boundary engine: this is a materially
##     harder quantity to calibrate than the alpha boundaries, because
##     the realised futility-stopping probability under H1 depends on
##     conditioning on non-crossing of the (already-fixed) alpha
##     boundary at every prior look, not just on the beta-spending
##     function's target value in isolation. A Monte Carlo check under
##     H1 (K=5, theta matching the design's power calculation) found
##     that the realised probability of stopping for non-binding
##     futility, and the realised power (probability of crossing the
##     efficacy boundary), both deviate from their nominal targets by a
##     non-trivial amount (single-digit percentage points in the
##     configuration checked). Because the futility boundary is
##     explicitly NON-BINDING in this package (crossing it is reported
##     as advisory only -- see `entered_futility_region` -- and never
##     triggers an automatic stopping decision, nor affects the alpha
##     boundary or the "DARIS reached" verdict), this imprecision does
##     not compromise the package's type-I error control, which is
##     governed solely by the independently-validated alpha engine
##     above. It does mean the futility band drawn on the plot should be
##     read as an approximate, illustrative guide to when the trend is
##     unpromising, not as a boundary with a precisely calibrated
##     operating characteristic. See `?tsa_cor` for further discussion.

## Two-sided (side = 2) O'Brien-Fleming-type Lan-DeMets alpha-spending
## function, matching RTSA's own parameterization (verified directly
## against a live RTSA::boundaries() call -- see VALIDATION above):
## the per-side alpha (alpha/2) is run through the one-sided OF-type
## spending shape -- i.e. divided by 2 again inside the normal
## quantile, giving alpha/4 -- and the result is doubled to combine the
## two symmetric (upper/lower) sides. This is algebraically RTSA's
## commented `alphas()` helper evaluated at side = 2:
##   2*(1 - pnorm(qnorm(1 - alpha/side/2)/sqrt(t))) * side
## which collapses to the 4*(1 - Phi(z_{alpha/4}/sqrt(t))) used here.
## ** CORRECTED in 0.2.6. ** Versions 0.2.0-0.2.5.1 used the more
## commonly-seen textbook two-argument form
## 2*(1-Phi(z_{alpha/2}/sqrt(t))) instead, which reaches alpha*(1) =
## alpha just as exactly but spends alpha meaningfully faster at early
## looks than RTSA's own boundaries do for the same schedule -- a
## genuine correctness bug relative to this package's own documented
## reference methodology, not a stylistic or parameterisation
## preference. See VALIDATION above for the numeric comparison.
.alpha_spend_OF <- function(t, alpha) {
  z <- stats::qnorm(1 - alpha / 4)
  4 * (1 - stats::pnorm(z / sqrt(t)))
}

## FFT-based Gaussian convolution step: given a (sub-)density h defined on
## an equally-spaced grid, returns the density after one Brownian-motion
## increment of variance dt (and optional mean shift), evaluated back on
## the SAME grid.
.convolve_step <- function(h, bgrid, db, dt, shift = 0) {
  n <- length(bgrid)
  lags <- seq(-(n - 1), (n - 1)) * db
  kernel <- stats::dnorm(lags - shift, sd = sqrt(dt))
  conv_result <- stats::convolve(h, rev(kernel), type = "open") * db
  conv_result[n:(2 * n - 1)]
}

## Alpha (efficacy) boundary at each information fraction in t (t[K]
## should be 1, i.e. the design's final/reference information size).
## Returns critical values in Z-space; Inf where the information fraction
## is too small for any finite boundary to be reached under the spending
## function (this is the mathematically correct answer -- no amount of
## evidence that early would justify stopping -- not a failure).
.obf_alpha_boundary <- function(t, alpha, bmax = 15, n_grid = 16000) {
  K <- length(t)
  bgrid <- seq(-bmax, bmax, length.out = n_grid)
  db <- bgrid[2] - bgrid[1]

  A <- .alpha_spend_OF(t, alpha)
  a_vec <- numeric(K)

  sd1 <- sqrt(t[1])
  a_vec[1] <- sd1 * stats::qnorm(1 - A[1] / 2)
  if (!is.finite(a_vec[1])) a_vec[1] <- Inf
  h <- stats::dnorm(bgrid, mean = 0, sd = sd1)
  h[abs(bgrid) >= a_vec[1]] <- 0

  if (K >= 2) {
    for (k in 2:K) {
      dt <- t[k] - t[k - 1]
      h_pre <- .convolve_step(h, bgrid, db, dt)

      target <- A[k] - A[k - 1]
      tail_prob <- function(a) 2 * sum(h_pre[bgrid >= a]) * db
      a_vec[k] <- if (tail_prob(bmax - 0.01) > target) {
        Inf
      } else {
        tryCatch(
          stats::uniroot(function(a) tail_prob(a) - target,
                          lower = 0, upper = bmax - 0.01)$root,
          error = function(e) Inf)
      }
      h <- h_pre
      h[abs(bgrid) >= a_vec[k]] <- 0
    }
  }
  a_vec / sqrt(t)
}

## Beta (non-binding futility) boundary at each information fraction,
## holding the already-computed alpha boundary (c_vec_alpha, in Z-space)
## fixed. Returns Z-space values; NA where no futility bound is computed
## (too early for a formal call, or -- by construction -- at the final
## look, since there is never a futility DECISION at a definitive final
## analysis).

## -------------------------------------------------------------------------
## -------------------------------------------------------------------------
## RTSA's actual non-binding beta-spending / futility-boundary engine
##
## ** RECONSTRUCTED in 0.2.7. ** Checked directly against the real RTSA
## R package source (RTSA 0.2.2, https://github.com/AnneLyng/RTSA:
## R/RTSA_helperfunctions.R, R/boundaries.R, src/first.cpp), which is NOT
## what the previous engine (0.2.4-0.2.6.9) was built from. That engine's
## own comments (later softened in 0.2.6.6-0.2.6.9 to "adapted from",
## but never actually rebuilt) described it as "a literal
## R implementation of the retrospective/analysis-mode 'inner wedge'
## algorithm supplied from RTSA's old TSA functions ('translated from
## java')" -- and indeed, none of the function names it cited as its
## source (betas_Obf(), sdfunc(), trap_old(), gfunc(), fcab_old(),
## qpos_old(), first_old(), other_old(), searchfunc_old(), getInnerWedge(),
## fakeIFY, tsa_beta_bound) appear anywhere in the actual RTSA 0.2.2
## source. It was a faithful implementation of a DIFFERENT, older tool,
## not of this package's own documented reference. See
## inst/REVERSE_ENGINEERING_RTSA.md for the full write-up of what that
## meant in practice and what changed here.
##
## This is a line-by-line port of RTSA's own beta_boundary(), z_n_w(),
## and searchfunc() (R/RTSA_helperfunctions.R), and of init_int(),
## recur_int(), and prob() (src/first.cpp -- ported to pure R here,
## since tsacor has no compiled-code dependency), specialised to the one
## configuration this package actually needs:
##   * es_beta = "esOF" (O'Brien-Fleming-type beta spending) -- the only
##     spending family tsacor exposes, matching .rtsa_beta_spend_OF()
##     below (unchanged from previous versions; its formula already
##     matched RTSA's real esOF(), just not the recursion around it);
##   * side = 1 passed internally to RTSA's own beta_boundary(), even
##     though the overall design is two-sided (side = 2). This is
##     RTSA's OWN convention, not a tsacor simplification: see
##     boundaries(), the side == 2 / futility == "non-binding" / type
##     == "design" branch, where beta_boundary() is always called with
##     side = 1 and an EXPLICIT delta override -- the futility "inner
##     wedge" is a one-sided construction against a fixed drift, not a
##     beta budget split across two sides;
##   * delta = |qnorm(alpha/2) + qnorm(beta)|, the same fixed
##     standardised drift RTSA computes explicitly in that branch
##     (`delta <- abs(qnorm(alpha/side)+qnorm(beta))` with the design's
##     actual side = 2) and passes into beta_boundary() -- a FIXED
##     theoretical quantity, not something derived from the shape of
##     the futility wedge itself (which is what the previous
##     "testDrift" engine did).
##
## The single most consequential difference from the previous engine:
## RTSA solves for each futility boundary DIRECTLY under this fixed
## alternative-hypothesis drift, using the ACTUAL alpha efficacy
## boundary (c_vec_alpha -- generally a decreasing-then-flattening
## sequence, not a constant) as the upper wall of the recursive
## integration at EVERY step, not a symmetric null-referenced
## construction that only gets reconciled with the alpha boundary,
## approximately, after the fact.
##
## RTSA's real algorithm also does NOT hide every non-positive futility
## value: it returns NA only where the pinned "zninf" sentinel (-20) was
## used because the incremental beta spend at that look was negligible
## (see boundaries.R's `abs(lb$za) == 20` checks). A negative non-binding
## futility boundary is a legitimate result -- it says the trial would
## only be flagged for advisory futility if the cumulative Z-statistic
## had already crossed to the "wrong" side of the null by that point --
## and previous versions of this package incorrectly suppressed it
## (`b[b <= 0] <- NA_real_`).
##
## One deliberate departure from a literal port: RTSA's own `boundaries()`
## achieves an exact meeting of the alpha and beta boundaries at the
## final, definitive look by root-finding (`uniroot(inf_warp, ...)`,
## rescaling its own "design" timing via `warp_root` until the two meet
## exactly) -- machinery that exists to solve a DESIGN problem (find the
## sample-size inflation that hits a target power for a not-yet-observed
## trial). tsacor has no such design phase: it operates directly on
## observed information fractions, with t = 1 defined as the DARIS/HARIS
## point itself. Rather than port RTSA's root-finding wholesale (a
## substantially different, prospective-design undertaking, not part of
## reconstructing the beta-spending recursion itself), this package
## keeps its previous, explicit convention: at the single definitive
## final look, there is no distinct "non-binding early stop for
## futility" apart from the main efficacy decision, so that look's
## futility boundary is set to equal the (exact, closed-form) final
## efficacy boundary qnorm(1-alpha/2) directly, exactly as previous
## versions of this package already did.
##
## VALIDATION STATUS (honest -- please read before trusting this for a
## real analysis): the alpha engine above was checked against a LIVE
## RTSA::boundaries() call (see VALIDATION notes higher in this file).
## This beta/futility reconstruction has NOT been -- no R interpreter is
## available in the environment that wrote this port, so none of the
## R code below has actually been executed, let alone compared against
## RTSA's real output. It is a careful, line-by-line reading of RTSA's
## published source, not a numerically confirmed match. Before relying
## on this for anything but an approximate, illustrative futility band
## (which is how this package has always described its futility output
## -- see the "Non-binding beta/futility boundary engine" VALIDATION
## note above), please run a direct comparison, e.g.:
##   design <- RTSA::boundaries(timing = c(0.2,0.4,0.6,0.8,1), alpha = 0.05,
##                               beta = 0.2, side = 2, futility = "non-binding",
##                               es_alpha = "esOF", es_beta = "esOF")
##   design$beta_ubound
##   tsacor:::.rtsa_beta_boundary(c(0.2,0.4,0.6,0.8,1), alpha = 0.05, beta = 0.2,
##                                c_vec_alpha = design$alpha_ubound)$boundary
## and compare the first K-1 entries directly (the final entry is
## tsacor's own definitive-look convention described above, not a
## quantity RTSA's un-rooted beta_boundary() output should be expected
## to match without its own root-finding applied).
## -------------------------------------------------------------------------

## Small defensive helper (0.2.7.2, NOT part of RTSA's own code): a
## seq(from, to, 2) that returns an empty integer vector instead of
## erroring when `to < from`, i.e. when the requested index range is
## empty. Needed because R's seq() with an explicit `by` argument
## throws "wrong sign in 'by' argument" in that case rather than
## returning integer(0) the way e.g. seq_len(0) does -- see
## .rtsa2_z_n_w()'s own safety-net comment for the general issue this
## belongs to. Used inside .rtsa2_z_n_w() in place of the raw
## seq(3, m - 2, 2) / seq(2, m - 1, 2) index-membership checks in its
## Simpson-weight loop, which are legitimately EMPTY (not erroring)
## ranges whenever the integration grid has few points -- m == 3 in
## particular (i.e. exactly 2 nodes in `xi`), which both the 0.2.7.1
## za/zb clamp and safety-net fix, and, independently, a naturally very
## narrow [za[i], zb[i]] interval near boundary convergence, can
## legitimately produce. Not exported.
.rtsa2_seq_by2 <- function(from, to) {
  if (to < from) integer(0) else seq(from, to, 2)
}

## Ported from RTSA's z_n_w() (R/RTSA_helperfunctions.R): builds a
## Simpson's-rule integration grid on [za[i], zb[i]] (Z-space at look
## i), with a log-spaced tail extension so the grid still reaches a
## possibly-distant boundary without wasting points on the bulk of the
## density.
.rtsa2_z_n_w <- function(r, info, za, zb, i, delta) {
  j <- 1:(6 * r - 1)
  xi <- delta * info$sd_incr[i] +
    (j < r) * (-3 - 4 * log(r / j)) +
    (r <= j & j <= 5 * r) * (-3 + 3 * (j - r) / (2 * r)) +
    (5 * r < j) * (3 + 4 * log(r / (6 * r - j)))

  if (sum(xi < za[i]) > 0) {
    indi <- max(which(xi < za[i]))
    xi <- xi[indi:length(xi)]
    xi[1] <- za[i]
  }
  if (sum(xi > zb[i]) > 0) {
    indi <- min(which(xi > zb[i]))
    xi <- xi[1:indi]
    xi[indi] <- zb[i]
  }

  ## Safety net for grid collapse: the two trims above are only
  ## guaranteed to leave >= 2 points in `xi` when za[i] < zb[i] with
  ## enough of a gap that at least one of the log-spaced nodes `j`
  ## falls strictly between them. If za[i] and zb[i] end up equal, or
  ## reversed (za[i] >= zb[i]), or close enough together that every
  ## node lands outside [za[i], zb[i]], `xi` can collapse to a single
  ## point (or, in principle, zero). That is a real edge case in this
  ## port (observed with the package's own example data under
  ## target_r = NA / the circular-target scenario), not merely a
  ## translation slip: RTSA's own R source has the identical unguarded
  ## seq(1, length(xi) - 1, 1) pattern below, which throws "wrong sign
  ## in 'by' argument" for length(xi) <= 1 rather than the empty
  ## sequence a length-based idiom like seq_len(length(xi) - 1) would
  ## give (NB: 1:0 is NOT such an idiom -- it evaluates to c(1, 0), a
  ## length-2 *descending* sequence, not an empty one; seq_len(0) is
  ## the correct empty-sequence idiom and is what the fix below uses).
  ## As of 0.2.7, this collapse is also actively prevented from
  ## occurring during normal operation by a clamp in
  ## .rtsa2_beta_boundary_core() that keeps za[i] from ever getting
  ## closer than a small margin to zb[i] in the first place -- this
  ## block is the last-resort fallback if it still happens (e.g. from
  ## a future caller invoking this function directly with unclamped
  ## za/zb). Rather than erroring, or silently returning a physically
  ## meaningless single-point "grid" (which would feed NA/garbage
  ## Simpson weights into the recursion via out-of-range zj[3]/zj[m-2]
  ## lookups below), we widen back out to the full two-point interval
  ## [za[i], zb[i]] -- inserting a minimal positive width if the two
  ## are equal or reversed -- so the Simpson's-rule construction below
  ## always has a valid (if minimal, 3-node) grid to build.
  if (length(xi) < 2L) {
    lo <- za[i]
    hi <- zb[i]
    if (!(hi > lo)) hi <- lo + 1e-8
    xi <- c(lo, hi)
  }

  m <- length(xi) * 2 - 1
  zj <- numeric(m)
  zj[seq(1, m, 2)] <- xi
  zj[seq(2, m - 1, 2)] <- (xi[seq_len(length(xi) - 1)] +
                             xi[seq(2, length(xi), 1)]) / 2

  ij <- 1:m
  wj <- numeric(m)
  for (k in ij) {
    if (k == 1) {
      wj[k] <- (1 / 6) * (zj[3] - zj[1])
    } else if (k %in% .rtsa2_seq_by2(3, m - 2)) {
      wj[k] <- (1 / 6) * (zj[k + 2] - zj[k - 2])
    } else if (k %in% .rtsa2_seq_by2(2, m - 1)) {
      wj[k] <- (4 / 6) * (zj[k + 1] - zj[k - 1])
    } else {
      wj[k] <- (1 / 6) * (zj[m] - zj[m - 2])
    }
  }
  list(zj = zj, wj = wj)
}

## Ported from src/first.cpp's init_int(): seeds the running density at
## the transition into look 2, evaluated on look 1's z_n_w() grid, as a
## Normal(delta * sd_incr[1], 1) density (i.e. under the fixed
## alternative-hypothesis drift), weighted by the Simpson weights.
## `stdv1` is the scalar info$sd_incr[1] -- the C++ original receives
## the full sd_incr vector but only ever reads its first element.
.rtsa2_init_int <- function(wj, zj, delta, stdv1) {
  wj * stats::dnorm(zj, mean = delta * stdv1, sd = 1)
}

## Ported from src/first.cpp's recur_int(): propagates the running
## density from look (k-1)'s grid (zj, on look (k-1)'s Z-scale) onto
## look k's grid (zj_up, on look k's Z-scale), using the independent-
## increment Gaussian transition density between the two looks'
## Y-scale (information-scale) values, under the fixed drift delta.
## `stdv` is a matrix with column 1 = sd_incr, column 2 = sd_proc (one
## row per look), matching RTSA's own `matrix(unlist(info), ncol = 2)`.
## RTSA's own beta_boundary() always calls this with bs = FALSE.
.rtsa2_recur_int <- function(k, stdv, zj, last, zj_up, wj_up, delta, bs) {
  y_prev <- zj * stdv[k - 1L, 2L]
  y_curr <- zj_up * stdv[k, 2L]
  sd_k <- stdv[k, 1L]
  ratio <- stdv[k, 2L] / sd_k

  last_up <- vapply(seq_along(zj_up), function(idx) {
    arg <- if (isTRUE(bs)) {
      (y_prev - y_curr[idx]) / sd_k
    } else {
      (y_curr[idx] - y_prev) / sd_k
    }
    sum(last * ratio * stats::dnorm(arg, mean = delta * sd_k, sd = 1))
  }, numeric(1))

  last_up * wj_up
}

## Ported from src/first.cpp's prob(): the cumulative crossing
## probability of a candidate boundary xq at look k, given the running
## density on look (k-1)'s grid (zj) and the fixed drift delta.
## searchfunc() always calls this with bs = TRUE for the beta/futility
## search.
.rtsa2_prob <- function(xq, last, zj, k, stdv, bs, delta) {
  y_prev <- zj * stdv[k - 1L, 2L]
  sd_k <- stdv[k, 1L]
  p <- if (isTRUE(bs)) {
    stats::pnorm((xq - y_prev) / sd_k, mean = delta * sd_k, sd = 1)
  } else if (delta != 0) {
    stats::pnorm((y_prev - xq) / sd_k, mean = -delta * sd_k, sd = 1)
  } else {
    stats::pnorm((y_prev - xq) / sd_k, mean = delta * sd_k, sd = 1)
  }
  sum(last * p)
}

## Ported from RTSA's searchfunc(): finds the boundary value (on the
## Y/information scale) at look i whose cumulative crossing probability
## (.rtsa2_prob()) equals the target incremental spend `as`, by an
## expanding-then-refining search -- ported literally, including its
## specific step-halving search pattern, rather than replaced with a
## smooth root-finder, so the numerics track RTSA's rather than merely
## its target. The `iter > 100000` safety cutoff below is NOT in RTSA's
## own code (its while(cond) loop has no iteration cap); it is a
## disclosed tsacor-side robustness addition so a pathological input
## cannot hang the R session, and only ever fires with a warning on
## non-convergence -- it does not alter the result in any case that
## would have converged under RTSA's own code.
.rtsa2_searchfunc <- function(last, zj, i, as, stdv, za, zb, tol, bs, delta) {
  maxnn <- 50
  upper <- zb[i - 1L] * stdv[i, 2L]
  if (isTRUE(bs)) upper <- za[i - 1L] * stdv[i, 2L]
  del <- 10
  qout <- .rtsa2_prob(xq = upper, last = last, zj = zj, k = i, stdv = stdv,
                       bs = bs, delta = delta)

  cond <- TRUE
  iter <- 0L
  while (cond) {
    iter <- iter + 1L
    if (abs(qout - as) <= tol) {
      cond <- FALSE
      break
    }
    if (qout > as + tol) {
      del <- del / 10
      for (k in 1:maxnn) {
        if (isTRUE(bs)) upper <- upper - 2 * del
        upper <- upper + del
        qout <- .rtsa2_prob(xq = upper, last = last, zj = zj, k = i,
                             stdv = stdv, bs = bs, delta = delta)
        if (qout <= as + tol) break
      }
    }
    if (qout < as - tol) {
      del <- del / 10
      for (k in 1:maxnn) {
        if (isTRUE(bs)) upper <- upper + 2 * del
        upper <- upper - del
        qout <- .rtsa2_prob(xq = upper, last = last, zj = zj, k = i,
                             stdv = stdv, bs = bs, delta = delta)
        if (qout >= as - tol) break
      }
    }
    if (iter > 100000L) {
      warning(sprintf(
        "RTSA-ported futility search did not converge at look %d (residual %.3g); using the closest value found.",
        i, abs(qout - as)
      ), call. = FALSE)
      break
    }
  }
  upper / stdv[i, 2L]
}

## RTSA's beta-spending increments (unchanged from previous versions --
## this formula already matched RTSA's real esOF(beta/side, timing) at
## the side = 1 convention beta_boundary() is actually called with; see
## the module-level VALIDATION note near the top of this file).
##
## RE-VERIFIED in 0.2.7.4 against a specific claim that this should be
## `qnorm(1 - beta)` (i.e. no `/2`), on the grounds that RTSA's
## `esOF(alpha, timing)` is itself `2*(1-pnorm(qnorm(1-alpha)/sqrt(t)))`.
## Checked directly against the actual RTSA 0.2.2 source
## (R/RTSA_helperfunctions.R): `esOF()` is
##   as_cum[i] <- 2*(1 - pnorm(qnorm(1-alpha/2)/sqrt(timing[i])))
## i.e. the `/2` IS present in the real function -- the claim was
## incorrect, and applying it would have reintroduced a real error
## (dividing the effective beta-spend rate in half, disproportionately
## understating early-look futility boundaries -- consistent with the
## "beta-bounds are far from RTSA, and further off at earlier looks"
## symptom that prompted the claim in the first place, but that
## symptom's actual cause was the missing information-scale root search
## below, not this formula). Do not change this back to `qnorm(1-beta)`
## without a fresh, direct read of RTSA's actual esOF() source.
.rtsa_beta_spend_OF <- function(t, beta, tol = 1e-13) {
  t <- as.numeric(t)
  if (any(!is.finite(t)) || any(t < 0))
    stop("information fractions must be finite and >= 0")
  if (!is.finite(beta) || beta <= 0 || beta >= 1)
    stop("beta must be strictly between 0 and 1")

  ## t == 0 is deliberately ALLOWED (0.2.7.4 fix): the 0.2.7.4 rm_bs
  ## suppression mechanism in .rtsa2_beta_boundary_core() (see below)
  ## zeroes the first `rm_bs` entries of the timing vector it passes in
  ## here, by design -- that is how RTSA's own two-pass root search
  ## suppresses an early look's spend after the first pass found a
  ## negative futility boundary there. This function's formula already
  ## handles t = 0 correctly with no special-casing needed:
  ## qnorm(1-beta/2) is always > 0 for beta in (0,1), so dividing by
  ## sqrt(0) = 0 gives +Inf (not NaN, and R does not warn on this),
  ## and pnorm(Inf, lower.tail = FALSE) = 0 exactly -- i.e. t = 0
  ## naturally gives zero cumulative spend, exactly the "no information,
  ## no spend yet" semantics that entry is meant to have. A stricter
  ## `t <= 0` guard here (as in 0.2.7-0.2.7.3) rejected these entries
  ## outright, which made every rm_bs > 0 call -- i.e. the entire second
  ## pass of the 0.2.7.4 root search -- fail unconditionally with this
  ## function's own error. That was the actual, near-universal cause of
  ## the "root search did not converge" failures reported against the
  ## first 0.2.7.4 build, not a genuine non-convergence of the root
  ## search itself (confirmed by direct reproduction: the same design,
  ## re-run in an independent Python port of this exact algorithm,
  ## converges immediately once this guard allows t = 0 through).
  cum <- 2 * stats::pnorm(
    stats::qnorm(1 - beta / 2) / sqrt(t),
    lower.tail = FALSE
  )
  delta <- c(cum[1L], diff(cum))
  delta[delta < tol] <- 0
  list(betaValuesCumulated = cum, betaValuesDelta = delta)
}

## Preventive clamp (0.2.7): keeps the non-binding futility boundary
## za[i] from ever landing closer than `gap` to the fixed efficacy wall
## zb[i] at the same look, which is what actually caused the grid
## collapse .rtsa2_z_n_w() now also guards against defensively (see its
## own comment) -- this is the primary fix; that one is the safety net.
## za[i] is set from one of three branches in the core recursion above
## (a fixed sentinel, a fixed 0, or an unconstrained root found by
## .rtsa2_searchfunc()); none of those three is otherwise guaranteed to
## stay below zb[i], and the sentinel/0 cases were not observed to
## collide with any realistic zb[i] in testing, but are clamped too for
## uniformity and future-proofing rather than relying on that holding.
## The observed failure mode (package example data, target_r = NA) was
## the .rtsa2_searchfunc() branch returning a value that landed at or
## past zb[i] at a late look, where the non-binding futility boundary is
## expected to approach the efficacy boundary closely by design (they
## are constructed to meet exactly at the final, t = 1 look) -- `gap` is
## deliberately tiny (Z-scale, not probability-scale) so this clamp only
## ever bites in that near-convergence regime and does not perturb any
## well-separated boundary value.
##
## IMPORTANT (0.2.7.4): this clamp is applied to every look EXCEPT the
## final one. The information-scale root search added in 0.2.7.4
## (.rtsa2_find_warp_root(), .rtsa2_inf_warp()) needs za[nn] to be able
## to reach -- and, for `uniroot()` to bracket a root at all, briefly
## cross -- zb[nn] exactly as the candidate warp factor is varied; a
## fixed-gap ceiling on za[nn] would make that impossible (the search
## objective, za[nn] - zb[nn], could then never change sign), so the
## final look is deliberately left unclamped here and instead protected
## purely by .rtsa2_z_n_w()'s own degenerate-grid fallback. Not
## exported.
.rtsa2_clamp_za_below_zb <- function(za_i, zb_i, gap = 1e-6) {
  min(za_i, zb_i - gap)
}

## Ported from RTSA's beta_boundary() (R/RTSA_helperfunctions.R),
## specialised to es_beta = "esOF" and to RTSA's own side = 1 convention
## for the futility construction (see the block comment above).
## `alpha_ubound` is the ALREADY-COMPUTED efficacy (alpha) boundary at
## every entry of `t`, used as the fixed upper wall throughout the
## recursion -- this (not a symmetric substitute) is what RTSA's own
## `zb <- alpha_bound` uses.
## `org_t` and `beta_timing`, when supplied directly, override the
## t*warp_root / rm_bs-zeroing derivation below -- this is what lets
## .rtsa_beta_boundary_analysis() (the type = "analysis"/design_R
## branch, added in 0.2.7.7) reuse this exact recursion with its own,
## differently-derived info-scale (the observed information fractions
## themselves, unwarped) and beta-spending timeline (observed fraction
## of the DESIGN's total information, inf_frac / design_R) instead of
## the design-mode t*warp_root pairing. When both are NULL (the
## original, design-mode call signature), behaviour is unchanged.
.rtsa2_beta_boundary_core <- function(t, beta, delta, alpha_ubound,
                                       warp_root = 1, rm_bs = 0L,
                                       org_t = NULL, beta_timing = NULL,
                                       zninf = -20, tol = 1e-15, r = 18) {
  t <- as.numeric(t)
  nn <- length(t)
  if (nn < 1L) stop("at least one information fraction is required")
  if (length(alpha_ubound) != nn)
    stop("alpha_ubound must have the same length as t")

  ## RTSA's own beta_boundary() keeps two DIFFERENT fraction scales
  ## alive at once: `beta_timing` (== `inf_frac`, possibly with its
  ## first `rm_bs` entries zeroed out) drives which fraction of beta is
  ## considered spent at each look, while `org_inf_frac` (==
  ## `inf_frac * warp_root`) drives the actual information/standard-
  ## deviation scale (info$sd_incr, info$sd_proc) the recursion runs
  ## on. These are NOT the same vector whenever warp_root != 1 -- see
  ## the block comment on .rtsa_beta_boundary() below for why a
  ## warp_root != 1 is needed at all.
  if (is.null(beta_timing)) {
    beta_timing <- t
    if (rm_bs > 0L) {
      beta_timing <- c(rep(0, rm_bs), beta_timing[-seq_len(rm_bs)])
    }
  }
  if (length(beta_timing) != nn)
    stop("beta_timing must have the same length as t")
  outbeta <- .rtsa_beta_spend_OF(beta_timing, beta, tol = tol)

  if (is.null(org_t)) org_t <- t * warp_root
  if (length(org_t) != nn)
    stop("org_t must have the same length as t")
  info <- list(sd_incr = sqrt(c(org_t[1L], diff(org_t))), sd_proc = sqrt(org_t))
  stdv <- cbind(info$sd_incr, info$sd_proc)  ## col 1 = sd_incr, col 2 = sd_proc

  za <- numeric(nn)
  zb <- alpha_ubound
  ya <- numeric(nn)
  yb <- zb * info$sd_proc

  d1 <- outbeta$betaValuesDelta[1L]
  d1 <- min(beta, max(0, d1))
  if (d1 == 0) {
    za[1L] <- zninf
  } else if (d1 == beta) {
    za[1L] <- 0
  } else {
    za[1L] <- stats::qnorm(d1, mean = info$sd_proc[1L] * delta, sd = 1)
  }
  ## Preventive clamp (0.2.7), skipped on the final look (0.2.7.4) --
  ## see .rtsa2_clamp_za_below_zb() above.
  if (nn > 1L) za[1L] <- .rtsa2_clamp_za_below_zb(za[1L], zb[1L])
  ya[1L] <- za[1L] * info$sd_incr[1L]

  zj_wj <- .rtsa2_z_n_w(r = r, info = info, za = za, zb = zb, i = 1L, delta = delta)
  zj <- zj_wj$zj
  wj <- zj_wj$wj
  last <- NULL

  if (nn >= 2L) {
    for (i in 2:nn) {
      if (i == 2L) {
        last <- .rtsa2_init_int(wj = wj, zj = zj, delta = delta,
                                 stdv1 = info$sd_incr[1L])
      }

      di <- outbeta$betaValuesDelta[i]
      if (di <= 0 || di >= 1) {
        di <- min(1, di)
        di <- max(0, di)
      }

      if (di < tol) {
        za[i] <- zninf
        if (i < nn) za[i] <- .rtsa2_clamp_za_below_zb(za[i], zb[i])
        ya[i] <- za[i] * info$sd_incr[i]
      } else if (di == beta) {
        za[i] <- 0
        if (i < nn) za[i] <- .rtsa2_clamp_za_below_zb(za[i], zb[i])
        ya[i] <- za[i] * info$sd_incr[i]
      } else {
        za[i] <- .rtsa2_searchfunc(
          last = last, zj = zj, i = i, as = di, stdv = stdv,
          za = za, zb = zb, tol = tol, bs = TRUE, delta = delta
        )
        ## Preventive clamp (0.2.7), skipped on the final look (0.2.7.4):
        ## this is the branch that actually triggered the grid-collapse
        ## crash for looks before the final one -- see
        ## .rtsa2_clamp_za_below_zb() above for why.
        if (i < nn) za[i] <- .rtsa2_clamp_za_below_zb(za[i], zb[i])
        ya[i] <- za[i] * info$sd_proc[i]
      }

      if (i != nn) {
        zj_wj_up <- .rtsa2_z_n_w(r = r, info = info, za = za, zb = zb,
                                  i = i, delta = delta)
        last <- .rtsa2_recur_int(
          k = i, stdv = stdv, zj = zj, last = last,
          zj_up = zj_wj_up$zj, wj_up = zj_wj_up$wj, delta = delta, bs = FALSE
        )
        zj <- zj_wj_up$zj
        wj <- zj_wj_up$wj
      }
    }
  }

  list(
    inf_frac = t,
    za = za,
    ya = ya,
    zb = zb,
    yb = yb,
    delta = delta,
    betaValuesCumulated = outbeta$betaValuesCumulated,
    betaValuesDelta = outbeta$betaValuesDelta,
    sdincr = info$sd_incr,
    sdproc = info$sd_proc
  )
}

## Ported from RTSA's inf_warp() (R/RTSA_helperfunctions.R): the
## objective function RTSA's own uniroot() call drives to zero. Returns
## the gap between the final look's futility boundary (under a
## candidate information-scale warp factor `x`) and the final look's
## fixed efficacy boundary -- zero exactly when the non-binding
## futility construction reaches the efficacy boundary precisely at the
## definitive final look, which is what a properly power-matched
## two-sided non-binding design requires.
.rtsa2_inf_warp <- function(x, t, beta, delta, alpha_ubound, rm_bs = 0L) {
  ans <- .rtsa2_beta_boundary_core(
    t = t, beta = beta, delta = delta, alpha_ubound = alpha_ubound,
    warp_root = x, rm_bs = rm_bs
  )
  ans$za[length(t)] - alpha_ubound[length(t)]
}

## Ported from the root-finding loop inside RTSA's boundaries()
## (side == 2, futility == "non-binding", type == "design" branch):
## `uniroot()` needs a bracket whose endpoints have opposite signs, and
## RTSA's own code does not know that bracket in advance, so it slides
## a narrow window of width `step` upward (starting at
## [start - step, start]) until `uniroot()` succeeds, for up to
## `max_iter` attempts. This is what actually determines the
## information-scale inflation ("warp_root", RTSA's `root`) needed for
## the non-binding futility construction to reach the fixed efficacy
## boundary exactly at the definitive final look -- see the block
## comment on .rtsa_beta_boundary() below for why this step exists at
## all and what happens without it.
.rtsa2_find_warp_root <- function(t, beta, delta, alpha_ubound, rm_bs = 0L,
                                   start = 0.95, step = 0.02, max_iter = 50L,
                                   tol_root = 1e-9) {
  f <- function(x) {
    .rtsa2_inf_warp(x, t = t, beta = beta, delta = delta,
                     alpha_ubound = alpha_ubound, rm_bs = rm_bs)
  }
  upper <- start
  for (n_itr in seq_len(max_iter)) {
    lower <- upper - step
    root <- tryCatch(
      stats::uniroot(f, lower = lower, upper = upper, tol = tol_root)$root,
      error = function(e) NULL
    )
    if (!is.null(root)) return(root)
    upper <- upper + step
  }
  stop(
    "Non-binding futility boundaries could not be computed (the RTSA-style ",
    "information-scale root search did not converge). Consider setting a ",
    "different target_r/power, or treat the futility band as unavailable ",
    "for this design."
  )
}

## Public internal entry point used by tsa_cor().
##
## Constructs RTSA's non-binding futility boundary sequence for the
## observed information-fraction timeline `t`, using the corresponding
## already-computed alpha (efficacy) boundary `c_vec_alpha` (same
## length/order as `t`) as the fixed upper wall throughout the
## recursion, and RTSA's own fixed standardised drift
## delta = |qnorm(alpha/2) + qnorm(beta)| (side = 2, i.e. this
## package's two-sided efficacy design).
##
## ** RTSA-style information-scale warp, added in 0.2.7.4. ** A naive
## reading of RTSA's beta_boundary() recursion (what 0.2.7-0.2.7.3
## implemented) uses the SAME information fractions both for looking up
## how much beta has been spent at each look and for the standard-
## deviation scale the recursion runs on -- i.e. implicitly always
## "warp_root = 1". RTSA's own boundaries() does NOT do this: it always
## calls beta_boundary() through a root-finding step
## (`uniroot(inf_warp, ...)`, ported above as .rtsa2_find_warp_root())
## that searches for a scalar inflation of the information SCALE alone
## (org_inf_frac = inf_frac * root; the beta-spending fractions
## themselves stay unwarped) such that the resulting futility boundary
## lands EXACTLY on the fixed efficacy boundary at the definitive final
## look. Without this, the futility boundaries a naive port produces
## are systematically too conservative relative to a live
## RTSA::boundaries() call across the whole schedule, and the final
## look generally does not meet the efficacy boundary at all -- exactly
## the two discrepancies reported against 0.2.7.3. RTSA's own code then
## repeats this root search a second time, now suppressing (zeroing the
## beta-spending fraction of) whichever early looks came back with a
## negative futility boundary on the first pass
## (`rm_bs = sum(lb$za < 0)`), and re-solves for the root again with
## that adjustment -- this is also why RTSA's real output does not show
## a small negative number at an early look where essentially no
## futility band exists yet: that look's beta-spend gets zeroed on the
## second pass, which routes it through the zninf sentinel (see
## .rtsa2_beta_boundary_core() above) instead. Both passes are ported
## faithfully below, including running the suppression step exactly
## once (RTSA's own code does not iterate this to convergence, so
## neither does this port).
.rtsa_beta_boundary <- function(t, alpha, beta, c_vec_alpha) {
  t <- as.numeric(t)
  if (length(t) != length(c_vec_alpha))
    stop("t and c_vec_alpha must have the same length")

  ## RTSA's own fixed drift for the non-binding futility construction
  ## (boundaries(), side == 2, futility == "non-binding", type ==
  ## "design" branch: `delta <- abs(qnorm(alpha/side)+qnorm(beta))`
  ## with side = 2), NOT derived empirically from the futility wedge
  ## itself.
  delta <- abs(stats::qnorm(alpha / 2) + stats::qnorm(beta))

  over_power <- any(t > 1)

  ## ** FIXED in 0.2.7.19. ** This used to assume "by construction of any
  ## properly normalised two-sided alpha-spending function, the efficacy
  ## boundary at t = 1 is exactly qnorm(1-alpha/2)" and substituted that
  ## constant here. That assumption is WRONG for RTSA's actual discretised
  ## O'Brien-Fleming-type spending recursion: the equality holds only for a
  ## SINGLE look (no earlier spending); with more looks the final wall is
  ## LARGER, and it grows with the number of looks (e.g. 2.014 for 4 even
  ## looks, 2.185 for 100 even looks at alpha = 0.05) -- it depends on the
  ## schedule, not on a continuous-monitoring limit. Confirmed against a
  ## live RTSA reconstruction: for one real
  ## schedule, RTSA's own final alpha boundary was 2.014090377368289, not
  ## qnorm(1-alpha/2) = 1.959963984540054 -- a difference large enough to
  ## materially shift the warp_root this function solves for (1.133242 vs.
  ## the wrong 1.097193) and, via the tightening step below, the reported
  ## final futility boundary itself (2.014090 vs. 1.959964).
  ##
  ## Fixed by recomputing the TRUE, schedule-dependent alpha boundary via
  ## .obf_alpha_boundary() -- the same side = 2 FFT recursion used
  ## everywhere else in this legacy engine, which is an APPROXIMATION of
  ## RTSA's own Simpson recursion, not a validated reproduction of it
  ## (final-look error 1.1e-3 on the 4-look schedule 0.25/0.5/0.75/1,
  ## measured 0.2.7.21) -- run on timing_beta itself (the
  ## observed pre-1 looks plus the definitive t = 1 point), rather than
  ## reusing whichever `c_vec_alpha` the caller happened to already have.
  ## This is correct and self-contained whether or not the caller's `t`
  ## already contains an exact t = 1 entry (compare .tsacor_legacy_
  ## boundaries(), which sometimes supplies one and sometimes does not),
  ## and it reproduces `c_vec_alpha`'s own pre-1 values exactly (same
  ## recursion, same input points) wherever a direct comparison is
  ## possible -- it is not a second, independent alpha engine.
  timing_beta <- sort(unique(c(t[t < 1], 1)))
  alpha_ubound_beta <- .obf_alpha_boundary(timing_beta, alpha = alpha)
  final_alpha_bound <- utils::tail(alpha_ubound_beta, 1L)

  if (length(timing_beta) == 1L) {
    ## Degenerate case (0.2.7.4): no pre-DARIS interim look at all --
    ## every observed information fraction already reached/exceeded 1
    ## before this call (e.g. a single very informative early study).
    ## The root search below has nothing to calibrate against: with
    ## beta_timing == 1 exactly, the incremental beta-spend at the one
    ## and only look is exactly `beta` by construction of any complete
    ## spending function, which always hits the `d1 == beta` special
    ## case (za = 0) in .rtsa2_beta_boundary_core() REGARDLESS of
    ## warp_root -- so .rtsa2_inf_warp() is a constant function of the
    ## warp factor and can never cross zero for any bracket, which is a
    ## genuine mathematical degeneracy, not a search failure. Skip the
    ## root search entirely and go straight to the definitive-look
    ## convention used below in every other case too.
    root2 <- 1
    rm_bs <- 0L
    ans <- list(
      za = final_alpha_bound,
      betaValuesCumulated = beta,
      betaValuesDelta = beta,
      sdincr = sqrt(timing_beta),
      sdproc = sqrt(timing_beta),
      ya = final_alpha_bound * sqrt(timing_beta)
    )
  } else {
    ## PASS 1: find the information-scale warp with no suppression yet
    ## (rm_bs = 0), matching RTSA's first uniroot()/beta_boundary() pair.
    root1 <- .rtsa2_find_warp_root(
      t = timing_beta, beta = beta, delta = delta,
      alpha_ubound = alpha_ubound_beta, rm_bs = 0L,
      start = 0.95, step = 0.02, max_iter = 50L
    )
    lb1 <- .rtsa2_beta_boundary_core(
      t = timing_beta, beta = beta, delta = delta,
      alpha_ubound = alpha_ubound_beta, warp_root = root1, rm_bs = 0L
    )

    ## PASS 2: re-find the warp with the early negative-za looks from
    ## PASS 1 suppressed (their beta-spend fraction zeroed), matching
    ## RTSA's second uniroot()/beta_boundary() pair. RTSA widens its
    ## search window by 0.05 (not 0.02) for this second pass.
    rm_bs <- sum(lb1$za < 0)
    root2 <- .rtsa2_find_warp_root(
      t = timing_beta, beta = beta, delta = delta,
      alpha_ubound = alpha_ubound_beta, rm_bs = rm_bs,
      start = 0.95, step = 0.05, max_iter = 50L
    )
    ans <- .rtsa2_beta_boundary_core(
      t = timing_beta, beta = beta, delta = delta,
      alpha_ubound = alpha_ubound_beta, warp_root = root2, rm_bs = rm_bs
    )
  }

  ## Definitive-final-look convention (kept from previous versions): no
  ## distinct non-binding futility zone separate from the efficacy
  ## decision at t = 1. With the warp/root search above, ans$za[last]
  ## should already be extremely close to final_alpha_bound BY
  ## CONSTRUCTION (that is exactly what the root search solves for);
  ## this line only tightens that to exact equality rather than leaving
  ## it at the root search's own numerical tolerance.
  za_seq <- ans$za
  za_seq[length(za_seq)] <- final_alpha_bound

  match_pre <- match(t, timing_beta)
  za_final <- utils::tail(za_seq, 1L)
  boundary <- ifelse(t < 1, za_seq[match_pre], za_final)

  ## Defensive tsacor-side safeguard, kept from previous versions and NOT
  ## part of RTSA's own algorithm: a non-binding futility boundary should
  ## never lie above the corresponding efficacy boundary. RTSA's own
  ## recursion does not hard-enforce this at every step -- only the grid
  ## construction two steps back clips indirectly -- so a pathological or
  ## internally inconsistent alpha/beta/delta combination could in
  ## principle let the search push a futility value above the efficacy
  ## wall; this pmin() makes that impossible regardless. For any
  ## self-consistent design (alpha, power, and delta mutually agreeing,
  ## as tsacor's own DARIS/HARIS construction intends) this should be a
  ## no-op in practice.
  boundary <- pmin(boundary, c_vec_alpha)

  ## RTSA returns NA only where the incremental beta spend at that look
  ## was negligible (za pinned to the zninf = -20 sentinel) -- NOT for
  ## every non-positive value (see the block comment above). With the
  ## rm_bs suppression from PASS 2 now wired in, early looks that would
  ## otherwise have shown a small negative number are routed through
  ## this same sentinel and hidden here, matching RTSA's real behaviour;
  ## any OTHER finite negative value RTSA's own (single-pass, non-
  ## iterated) suppression does not happen to catch is, correctly, left
  ## visible rather than hidden.
  zninf <- -20
  boundary[abs(boundary - zninf) < 1e-8] <- NA_real_

  beta_cum <- ans$betaValuesCumulated
  beta_delta <- ans$betaValuesDelta
  beta_spent <- ifelse(t < 1, beta_cum[match_pre], beta)
  beta_spent_delta <- ifelse(
    t < 1, beta_delta[match_pre],
    beta - utils::tail(beta_cum, 2L)[1L]
  )
  auxiliary <- ifelse(t < 1, za_seq[match_pre], NA_real_)

  list(
    boundary = boundary,
    beta_spent = beta_spent,
    beta_spent_raw = beta_spent,
    beta_spent_delta = beta_spent_delta,
    auxiliary_futility = auxiliary,
    delta = delta,
    fakeIFY = final_alpha_bound,
    over_power = over_power,
    timing_beta = timing_beta,
    warp_root = root2,
    rm_bs = rm_bs,
    sdincr = ans$sdincr,
    sdproc = ans$sdproc,
    ya = ans$ya,
    yb = ans$yb
  )
}

## Backward-compatible helper: returns only the boundary vector.
.obf_beta_boundary <- function(t, alpha, beta, c_vec_alpha) {
  .rtsa_beta_boundary(
    t = t, alpha = alpha, beta = beta, c_vec_alpha = c_vec_alpha
  )$boundary
}

## -------------------------------------------------------------------------
## RTSA's type = "analysis" (design_R) beta/futility branch
##
## ** Added in 0.2.7.7. 0.2.7.9 introduced a REGRESSION here (reverted in
## 0.2.7.10) by switching design_R/delta/rm_bs to a side = 1,
## futility = "none", right_power()-based calibration, on the strength of
## an external "live RTSA reconstruction" that turned out to have queried
## the wrong branch of RTSA's own source. That reconstruction is WRONG,
## and has been directly falsified by re-reading RTSA's actual R/RTSA.R
## top-level wrapper (not just boundaries.R in isolation), specifically
## the branch that manufactures design_R when no prior design object is
## supplied for a retrospective analysis:
##
##   bounds <- boundaries(timing = timing, alpha = alpha, beta = beta,
##                         side = side, futility = futility,
##                         es_alpha = es_alpha, es_beta = es_beta,
##                         type = "design")
##   design_R <- bounds$root
##
## `side = side` and `futility = futility` here are RTSA()'s OWN top-level
## arguments -- i.e. whatever the caller passed to RTSA() itself (side = 2,
## futility = "non-binding" for the two-sided, non-binding-futility TSA
## this package performs), NOT hardcoded to side = 1 / futility = "none".
## So the design_R calibration call this package needs to reproduce is
## boundaries(..., side = 2, futility = "non-binding", type = "design")
## -- exactly RTSA's non-binding-futility, two-pass warp_root search,
## i.e. exactly .rtsa_beta_boundary() above, unchanged, called on the
## observed timing. That is what 0.2.7.7 did. 0.2.7.9's side = 1 /
## futility = "none" / right_power() path (boundaries.R lines ~69-90) is
## a genuinely different RTSA code path -- the plain power/sample-size
## root search used for a design that has NO futility boundaries at all
## -- and is simply the wrong branch for this package's non-binding,
## two-sided design, regardless of what any external reconstruction
## computed against it.
##
## .rtsa_beta_boundary() above (unchanged since 0.2.7.6) is a faithful
## port of RTSA's boundaries(..., side = 2, futility = "non-binding",
## type = "design") branch: given only a timing vector, it root-finds its
## own information-scale inflation (warp_root) from scratch so the
## non-binding futility construction meets the efficacy boundary exactly
## at the final look. That IS RTSA's own design_R for this design, per
## the RTSA.R trace above -- not a mismatched substitute for it.
##
## RTSA's actual type = "analysis" branch this package's own design
## reaches is boundaries.R's side == 2 / futility == "non-binding" /
## `else` (i.e. type == "analysis") branch (lines ~436-492): it does NOT
## run its own root search. It takes the design_R already solved above
## and calls beta_boundary() THREE times against it -- once unsuppressed,
## then twice more with rm_bs re-derived from the previous pass's
## negative-za count -- see RTSA_helperfunctions.R's beta_boundary(), the
## `if(!is.null(design_R))` block: the OBSERVED information fractions
## (`inf_frac`) drive the actual info/sd scale UNWARPED
## (`org_inf_frac <- inf_frac`, only appending `design_R` itself as a
## final point if the observed data has not yet reached it), while the
## beta-SPENDING budget consumed at each look is looked up against
## `inf_frac / design_R` -- i.e. against how much of the eventually-
## planned total information has actually accrued, not against the raw
## observed fraction directly. The fixed drift `delta` this branch uses
## is `abs(qnorm(alpha/side)+qnorm(beta))` evaluated at side = 2 (the
## design's own overall sidedness, matching this package's alpha engine
## and .rtsa_beta_boundary() above) -- NOT side = 1.
##
## tsacor has no separate prospective "design" phase of its own -- each
## tsa_cor() call operates retrospectively on whatever information
## fractions the included studies actually produced. This function
## reproduces RTSA's own two-step pipeline for that situation: when
## design_R is not supplied, it is calibrated by calling
## .rtsa_beta_boundary() (the side = 2, non-binding-futility, two-pass
## warp_root search) on the same observed timing `t`, then the design_R
## branch below is run against that root.
##
## VALIDATION STATUS: a careful, line-by-line reading of RTSA's published
## alpha_boundary()/beta_boundary()/boundaries()/RTSA() source, re-checked
## against R/RTSA.R's own top-level wrapper (not just boundaries.R read in
## isolation) after 0.2.7.9's regression. No R interpreter has been
## available in any environment that has worked on this port, so this
## remains a source trace, not a numerically executed match. Before
## relying on this for anything but an approximate, illustrative futility
## band, run a direct comparison against a live two-step RTSA call on
## your own data, e.g.:
##   design <- RTSA::boundaries(timing = c(<observed fractions>, 1),
##                               alpha = 0.05, beta = 0.2, side = 2,
##                               futility = "non-binding", es_alpha = "esOF",
##                               es_beta = "esOF", type = "design")
##   analysis <- RTSA::boundaries(timing = c(<observed fractions capped
##                                  at design$root>), alpha = 0.05,
##                                  beta = 0.2, side = 2,
##                                  futility = "non-binding",
##                                  es_alpha = "esOF", es_beta = "esOF",
##                                  type = "analysis", design_R = design$root)
##   analysis$beta_ubound; analysis$root
##   tsacor:::.rtsa_beta_boundary_analysis(<same observed fractions>,
##                                          alpha = 0.05, beta = 0.2,
##                                          c_vec_alpha = <matching
##                                          alpha_ubound>)
## If you can run this comparison and it disagrees, please report the
## EXACT boundaries() call you used (side, futility, es_alpha, es_beta)
## alongside the numbers -- the 0.2.7.9 regression happened precisely
## because a side/futility mismatch was not stated explicitly.
## -------------------------------------------------------------------------

.rtsa_beta_boundary_analysis <- function(t, alpha, beta, c_vec_alpha,
                                          design_R = NULL) {
  t <- as.numeric(t)
  if (length(t) != length(c_vec_alpha))
    stop("t and c_vec_alpha must have the same length")

  ## Same fixed theoretical drift as the design-mode engine, at this
  ## design's own side = 2 -- see the block comment above for why side = 1
  ## (introduced in 0.2.7.9, reverted here) was wrong.
  delta <- abs(stats::qnorm(alpha / 2) + stats::qnorm(beta))
  final_alpha_bound <- stats::qnorm(1 - alpha / 2, lower.tail = TRUE)

  ## STEP 1 (only when design_R is not externally supplied): manufacture
  ## a design_R the same way RTSA's own RTSA() wrapper does when no prior
  ## design object exists -- run the side = 2, futility = "non-binding",
  ## type = "design" root search (.rtsa_beta_boundary(), PASS 1 + PASS 2)
  ## on the observed timing itself, and take its warp_root as design_R.
  if (is.null(design_R)) {
    design_fit <- .rtsa_beta_boundary(
      t = t, alpha = alpha, beta = beta, c_vec_alpha = c_vec_alpha
    )
    design_R <- design_fit$warp_root
  }
  if (!is.finite(design_R) || design_R <= 0)
    stop("design_R must be a finite, strictly positive scalar")

  ## STEP 2: extend/trim the observed timing onto the design endpoint,
  ## exactly as RTSA's own RTSA() wrapper does before its type =
  ## "analysis" call (R/RTSA.R): append design_R if the observed data has
  ## not yet reached it; otherwise keep only the sub-design_R looks plus
  ## design_R itself as the definitive point.
  if (max(t) < design_R) {
    t_ext <- c(t, design_R)
  } else if (max(t) > design_R) {
    t_ext <- c(t[t < design_R], design_R)
  } else {
    t_ext <- t
  }

  ## Alpha (efficacy) boundary at each point of t_ext: reuse the
  ## already-computed c_vec_alpha wherever t_ext coincides with an
  ## original observed look, and tsacor's own established definitive-
  ## endpoint convention (final_alpha_bound = qnorm(1-alpha/2)) for the
  ## design_R endpoint itself -- the same convention .rtsa_beta_boundary()
  ## already applies at its own t = 1 endpoint above.
  match_ext <- match(t_ext, t)
  alpha_ubound_ext <- ifelse(
    is.na(match_ext), final_alpha_bound, c_vec_alpha[match_ext]
  )

  if (length(t_ext) == 1L) {
    ## Degenerate case, mirroring .rtsa_beta_boundary()'s own single-look
    ## fallback above: nothing to calibrate a beta-spending timeline
    ## against. Route straight to the definitive-look convention.
    za_ext <- final_alpha_bound
    beta_cum_ext <- beta
    beta_delta_ext <- beta
    rm_bs <- 0L
  } else {
    ## RTSA-style beta-spending timeline: fraction of the DESIGN's total
    ## planned information consumed so far (inf_frac / design_R), not
    ## the raw observed fraction -- the analysis-mode counterpart of the
    ## design engine's t / warp_root pairing. Trimmed to strictly-less-
    ## than-1 entries plus one synthetic final point at exactly 1
    ## (design fully spent), matching RTSA's own beta_boundary() design_R
    ## branch.
    beta_timing_raw <- t_ext / design_R
    beta_timing_raw <- c(beta_timing_raw[beta_timing_raw < 1], 1)
    if (length(beta_timing_raw) != length(t_ext))
      stop("internal error: beta timing and info-scale lengths diverged")

    ## RTSA's own boundaries() calls beta_boundary() THREE times for this
    ## (side = 2, futility = "non-binding", type = "analysis") branch:
    ## once unsuppressed, then twice more with rm_bs re-derived from the
    ## previous pass's negative-za count -- converging the suppression
    ## fixed point WITHOUT any further root search (design_R itself stays
    ## fixed throughout; only which early looks get suppressed can
    ## change). No warp_root/inf_warp search happens at this level at
    ## all.
    rm_bs <- 0L
    lb <- NULL
    for (pass in 1:3) {
      beta_timing <- beta_timing_raw
      if (rm_bs > 0L) {
        beta_timing <- c(rep(0, rm_bs), beta_timing[-seq_len(rm_bs)])
      }
      lb <- .rtsa2_beta_boundary_core(
        t = t_ext, beta = beta, delta = delta, alpha_ubound = alpha_ubound_ext,
        org_t = t_ext, beta_timing = beta_timing
      )
      rm_bs <- sum(lb$za < 0)
    }

    za_ext <- lb$za
    beta_cum_ext <- lb$betaValuesCumulated
    beta_delta_ext <- lb$betaValuesDelta

    ## Same non-binding-vs-efficacy safeguard RTSA's own analysis branch
    ## applies at its last look (boundaries.R: if the futility bound at
    ## the final look exceeds the efficacy bound there, clip it down to
    ## the efficacy bound).
    nnn <- length(za_ext)
    if (za_ext[nnn] > alpha_ubound_ext[nnn]) za_ext[nnn] <- alpha_ubound_ext[nnn]
  }

  ## STEP 3: return only the boundaries at the ACTUALLY OBSERVED looks
  ## (t), dropping the synthetic design endpoint whenever it was
  ## appended rather than genuinely observed.
  match_pre <- match(t, t_ext)
  boundary <- za_ext[match_pre]

  ## Over-powered observed looks (t > design_R, e.g. a single very
  ## informative early study that alone exceeds the design's planned
  ## total information) were trimmed out of t_ext entirely above and so
  ## have no match there; mirror .rtsa_beta_boundary()'s own convention
  ## for this case (its `over_power`/t > 1 branch) by giving them the
  ## same definitive final boundary as the design endpoint itself.
  over_power_look <- t > design_R
  if (any(over_power_look)) boundary[over_power_look] <- final_alpha_bound

  ## RTSA's own boundaries() converts every look pinned exactly at the
  ## +/-20 sentinel (i.e. the recursion never bothered to solve a real
  ## boundary there because rm_bs suppression zeroed its target spend) to
  ## NA before returning: `beta_ubound <- c(rep(NA, sum(abs(za)==20)),
  ## za[abs(za)<20])`. This is what makes RTSA's own printed early-look
  ## futility boundaries show as blank/NA rather than as some finite
  ## value close to -20 -- it is not a display-layer nicety tsacor adds on
  ## top; it is RTSA's own numeric convention, reproduced here.
  zninf <- -20
  boundary[abs(boundary - zninf) < 1e-8] <- NA_real_

  ## Same defensive tsacor-side safeguard as the design-mode engine: a
  ## non-binding futility boundary should never lie above the
  ## corresponding efficacy boundary.
  boundary <- pmin(boundary, c_vec_alpha)

  ## Same "fully spent" convention .rtsa_beta_boundary() uses for its own
  ## t >= 1 over-powered entries, here relative to design_R instead of 1.
  beta_spent <- ifelse(over_power_look, beta, beta_cum_ext[match_pre])
  beta_spent_delta <- ifelse(
    over_power_look, beta - utils::tail(beta_cum_ext, 2L)[1L],
    beta_delta_ext[match_pre]
  )

  list(
    boundary = boundary,
    delta = delta,
    design_R = design_R,
    over_power = any(over_power_look),
    t_ext = t_ext,
    alpha_ubound_ext = alpha_ubound_ext,
    rm_bs = rm_bs,
    beta_spent = beta_spent,
    beta_spent_delta = beta_spent_delta
  )
}

## Backward-compatible helper: returns only the boundary vector.
.obf_beta_boundary_analysis <- function(t, alpha, beta, c_vec_alpha,
                                         design_R = NULL) {
  .rtsa_beta_boundary_analysis(
    t = t, alpha = alpha, beta = beta, c_vec_alpha = c_vec_alpha,
    design_R = design_R
  )$boundary
}

## -------------------------------------------------------------------------
## Pre-0.2.7.11 boundary pipeline, kept ONLY as a fallback for tsa_cor() when
## the RTSA-exact compiled engine (R/rtsa_engine.R) cannot produce a result.
## Approximate: its alpha engine is an FFT convolution (not RTSA's Simpson
## recursion). Until 0.2.7.19 its futility route also substituted
## qnorm(1 - alpha/2) for the final efficacy wall; .rtsa_beta_boundary() now
## recomputes that wall from the FFT alpha recursion (the analysis wrapper's
## design_R endpoint still uses the constant). Measured on the reference
## schedule 0.25/0.5/0.75/1 (0.2.7.21): warp root 1.132483 vs live RTSA
## 1.133242 (error 7.6e-4; 3.6e-2 before the 0.2.7.19 fix), and the whole
## residual comes from the FFT alpha input (final-look alpha error 1.1e-3).
## See NEWS.md, 0.2.7.19 and 0.2.7.21.
## -------------------------------------------------------------------------
.tsacor_legacy_boundaries <- function(info_fracs, boundary_timing, alpha, beta) {
  alpha_bounds_design <- .obf_alpha_boundary(boundary_timing, alpha = alpha)
  beta_unique_fracs <- sort(unique(info_fracs))
  beta_alpha_ref_observed <- alpha_bounds_design[
    match(pmin(beta_unique_fracs, 1), boundary_timing)
  ]
  beta_engine <- .rtsa_beta_boundary_analysis(
    beta_unique_fracs, alpha = alpha, beta = beta,
    c_vec_alpha = beta_alpha_ref_observed
  )
  beta_pre_daris <- beta_engine$boundary[
    match(boundary_timing[boundary_timing < 1], beta_unique_fracs)
  ]
  list(alpha_bounds_design = alpha_bounds_design,
       beta_pre_daris = beta_pre_daris,
       beta_engine = beta_engine)
}
