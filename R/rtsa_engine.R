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
## -------------------------------------------------------------------------
## RTSA-derived boundary engine (added in 0.2.7.11)
##
## Thin R layer over src/rtsa_core.h / src/rtsa_engine.cpp, which port RTSA
## 0.2.2's recursive-integration engine (alpha_boundary(), beta_boundary(),
## z_n_w(), searchfunc(), init_int(), recur_int(), prob()) end-to-end, and
## over R's own stats::uniroot(), exactly as RTSA's boundaries() uses it.
## This file only reproduces RTSA's *orchestration* (which root searches are
## run, in which order, on which information/spending scales).
##
## THE TWO ROUTES RTSA HAS.  RTSA::boundaries(side = 2, futility =
## "non-binding") computes beta (futility) bounds in two different ways,
## both through the same beta_boundary() called with side = 1:
##
##   type = "design"   -> .rtsa_design_bounds()
##       spending timeline = the looks `t` themselves; information scale =
##       t * warp_root, where warp_root ("root") is found by uniroot so that
##       the futility bound meets the efficacy bound at the final look
##       (t = 1).  Two root searches (second with the negative early
##       looks' spend suppressed, rm_bs).
##   type = "analysis" -> .rtsa_analysis_bounds()
##       needs a design_R (the design's root); spending timeline = t /
##       design_R, information scale = t (unwarped), delta fixed, rm_bs
##       iterated three times, no root search.  The alpha bounds this route
##       is run against are recomputed on the timeline t / design_R (they
##       are NOT the design-type alpha bounds).
##
## ** 0.2.7.13: ** tsa_cor() gained a `boundary_route = c("design", "analysis")`
## argument. "design" (the default, unchanged) uses .rtsa_design_bounds()
## only (DARIS = t = 1). "analysis" uses .rtsa_retrospective() (design pass
## for design_R, then the analysis pass against it) and moves the formal
## endpoint to design_R * DARIS, exactly as RTSA's own RTSA(type =
## "analysis", design = NULL) does -- see tsa_cor()'s own documentation and
## NEWS.md (0.2.7.13) for what that changes downstream (the "DARIS reached"
## verdict, the plotted endpoint, and the decision layer, not just the
## futility numbers).
##
## .rtsa_retrospective() chains them the way RTSA::RTSA(type = "analysis",
## design = NULL) does: a design pass on the observed looks gives design_R,
## then the analysis pass runs against it.
##
## IMPORTANT (the gap this file closes).  In both routes the efficacy wall
## at the FINAL look is the value the alpha recursion returns for that look
## (e.g. 2.127 for a 9-look schedule whose looks are dense near 1), not
## qnorm(1 - alpha/2) = 1.96.  tsahr <= 0.2.7.10 (the ancestor of tsacor) substituted 1.96; that
## changes design_R and therefore every beta bound.
## -------------------------------------------------------------------------

.rtsa_null_to_na <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)

## ** 0.2.7.13: ** every call into the compiled engine below can now come
## back carrying `grid_collapses`/`slow_searches` diagnostic counters (see
## src/rtsa_core.h). Both are departures from RTSA's own (uncapped / never-
## widening) behaviour, so both are surfaced immediately as R warnings --
## never silently absorbed -- right where the compiled result is received,
## so every caller (design pass, analysis pass, either boundary) gets them
## without having to remember to check.
.rtsa_warn_diagnostics <- function(res, where) {
  ## ** 0.2.7.14: ** a REVERSED integration interval (lower wall above the
  ## upper wall) is reported on its own, more alarming warning, distinct
  ## from a merely degenerate (zero-width) grid: widening a reversed interval
  ## makes an inverted boundary state look like an ordinary computation.
  ## It warns rather than errors because the same state can occur
  ## transiently at a bad candidate information scale inside a root search
  ## (those candidate evaluations are not warned about -- see `warn` in
  ## .rtsa_beta_cpp() -- only the converged passes are).
  gr_n <- res$grid_reversed
  if (!is.null(gr_n) && isTRUE(gr_n > 0L)) {
    warning(sprintf(
      paste0("*** RTSA-ported integration interval was REVERSED (lower wall ",
             "above the upper wall) %d time(s) while computing %s. The ",
             "interval was widened to a minimal non-zero window so the ",
             "calculation could continue, but a reversed configuration is ",
             "NOT a valid RTSA computation (RTSA's own R code would have ",
             "errored): the boundaries from the look where this happened ",
             "onward must not be trusted or reported as RTSA-equivalent. ",
             "Check the design (information fractions, alpha, power). ***"),
      gr_n, where
    ), call. = FALSE, immediate. = TRUE)
  }
  gc_n <- res$grid_collapses
  if (!is.null(gc_n) && isTRUE(gc_n > 0L)) {
    warning(sprintf(
      paste0("RTSA-ported integration grid collapsed to a degenerate ",
             "(fewer-than-two-node or zero-width) interval %d time(s) while ",
             "computing %s; the calculation continued (an interval with no ",
             "positive width is widened to a minimal non-zero window). RTSA's ",
             "own R code would have errored (\"wrong sign in 'by' argument\") ",
             "at the first one instead -- treat this boundary sequence with ",
             "extra caution at the look(s) where it happened."),
      gc_n, where
    ), call. = FALSE)
  }
  ss_n <- res$slow_searches
  if (!is.null(ss_n) && isTRUE(ss_n > 0L)) {
    warning(sprintf(
      paste0("RTSA-ported boundary search needed the slow (iteration-capped) ",
             "path %d time(s) while computing %s; the result still converged ",
             "within a loose tolerance and was accepted, but this is slower ",
             "and less certain than RTSA's own uncapped search."),
      ss_n, where
    ), call. = FALSE)
  }
  res
}

## RTSA alpha_boundary(side, esOF) -- compiled.
.rtsa_alpha_cpp <- function(inf_frac, side, alpha, design_R = NULL,
                            tol = 1e-9, r = 18L, where = "the alpha (efficacy) boundary") {
  inf_frac <- as.numeric(inf_frac)
  ## 0.2.7.16: !is.finite() (not anyNA()) so Inf / -Inf / NaN are rejected too;
  ## anyNA() lets Inf through and the compiled recursion then returns garbage.
  if (!length(inf_frac) || any(!is.finite(inf_frac)) || any(inf_frac <= 0))
    stop("information fractions must be finite and strictly positive")
  res <- rtsa_alpha_boundary_cpp(inf_frac, as.integer(side), as.numeric(alpha),
                                 .rtsa_null_to_na(design_R), as.numeric(tol),
                                 as.integer(r))
  .rtsa_warn_diagnostics(res, where)
}

## RTSA beta_boundary(side = 1, esOF) -- compiled.  `warp_root` (design
## route) and `design_R` (analysis route) are NULL when unused, as in RTSA.
##
## `warn = FALSE` (0.2.7.14) skips the diagnostics warnings.  It is used ONLY
## for the candidate evaluations inside the root searches, where a transient
## bad candidate information scale can legitimately produce a degenerate or
## reversed interior grid without affecting the converged answer; the
## converged passes always run with warn = TRUE, so a persistent problem at
## the accepted root is still reported.
##
## ** 0.2.7.16 -- READ BEFORE USING `$za`. ** When `$unreachable_look > 0` the
## entries of `$za` from that look onward are NOT boundaries: they are the
## deliberately non-physical sentinel `zb + 1` ("beyond the efficacy wall"),
## present only so that the sign of the final gap `zb[n] - za[n]` -- what the
## information-scale root searches use -- is negative. They must never be
## reported, plotted or interpreted as futility bounds. (Reachability is
## decided in C++ directly from the analytic limit -- target spend > sum of the
## surviving mass, up to the numerical tolerance.)
.rtsa_beta_cpp <- function(inf_frac, alpha_bound, beta, delta, rm_bs = 0L,
                           design_R = NULL, warp_root = NULL, side = 1L,
                           zninf = -20, tol = 1e-15, r = 18L,
                           where = "the beta (futility) boundary",
                           warn = TRUE) {
  if (!length(inf_frac) || any(!is.finite(inf_frac)))
    stop("information fractions must be finite")
  res <- rtsa_beta_boundary_cpp(as.numeric(inf_frac), as.numeric(alpha_bound),
                                as.numeric(beta), as.integer(side), as.numeric(delta),
                                as.integer(rm_bs), .rtsa_null_to_na(design_R),
                                .rtsa_null_to_na(warp_root), as.numeric(zninf),
                                as.numeric(tol), as.integer(r))
  if (isTRUE(warn)) .rtsa_warn_diagnostics(res, where) else res
}

## RTSA's "slide a narrow bracket upward until uniroot() succeeds" loop.
.rtsa_slide_root <- function(f, start, step, max_iter = 50L, tol = 1e-9) {
  upper <- start
  last_err <- NULL   # 0.2.7.15: last failure that was NOT just "no sign change"
  for (n_itr in seq_len(max_iter)) {
    root <- tryCatch(
      stats::uniroot(f, lower = upper - step, upper = upper, tol = tol)$root,
      error = function(e) {
        if (!grepl("opposite sign", conditionMessage(e), fixed = TRUE))
          last_err <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(root)) return(root)
    upper <- upper + step
  }
  stop("RTSA-style information-scale root search did not converge ",
       "(no sign change found in [", start - step, ", ",
       start + step * (max_iter - 1L), "])",
       if (!is.null(last_err)) paste0("; last engine error: ", last_err) else "",
       ".")
}

## RTSA::boundaries(timing = t, alpha, beta, side = 2,
##                  futility = "non-binding", es_alpha = "esOF",
##                  es_beta = "esOF", type = "design")
## A final look at t = 1 is appended when max(t) < 1, as RTSA does.
.rtsa_design_bounds <- function(t, alpha, beta) {
  t <- as.numeric(t)
  if (!length(t) || any(!is.finite(t)) || any(t <= 0))
    stop("information fractions must be finite and strictly positive")
  if (is.unsorted(t, strictly = TRUE))
    stop("information fractions must be strictly increasing")
  if (max(t) < 1) t <- c(t, 1)
  nt <- length(t)
  .rtsa_check_look_spacing(t, "the design-route look schedule")

  ab <- .rtsa_alpha_cpp(t, side = 2L, alpha = alpha,
                       where = "the design-route alpha (efficacy) boundary")
  ub <- ab$alpha_ubound
  delta <- abs(stats::qnorm(alpha / 2) + stats::qnorm(beta))

  if (nt == 1L) {
    ## RTSA cannot run its beta recursion on a single look (`for (i in
    ## 2:nn)` with nn = 1); there is nothing to calibrate.
    return(list(timing = t, alpha_ubound = ub, za = ub, beta_ubound = ub,
                root = 1, rm_bs = 0L, delta = delta, beta_spent = beta,
                beta_spent_delta = beta, sentinel = FALSE))
  }

  gap <- function(x, rm) {
    lb <- .rtsa_beta_cpp(t, ub, beta, delta, rm_bs = rm, warp_root = x,
                        where = "the design-route beta (futility) boundary (root search)",
                        warn = FALSE)
    ub[nt] - lb$za[nt]
  }
  ## PASS 1: RTSA, nn_max = 50 windows of width 0.02 starting at [.93, .95]
  root1 <- .rtsa_slide_root(function(x) gap(x, 0L), start = 0.95, step = 0.02)
  lb1 <- .rtsa_beta_cpp(t, ub, beta, delta, rm_bs = 0L, warp_root = root1,
                       where = "the design-route beta (futility) boundary (pass 1)")
  .rtsa_check_converged_pass(lb1, ub[nt], "the design-route calibration (pass 1)")
  rm_bs <- sum(lb1$za < 0)
  if (rm_bs >= nt)
    stop("every look has a negative futility bound on the first pass; ",
         "the non-binding design cannot be computed.")
  ## PASS 2: negative early looks get zero beta spend; windows of width 0.05
  root2 <- .rtsa_slide_root(function(x) gap(x, as.integer(rm_bs)),
                            start = 0.95, step = 0.05)
  lb <- .rtsa_beta_cpp(t, ub, beta, delta, rm_bs = as.integer(rm_bs),
                       warp_root = root2,
                       where = "the design-route beta (futility) boundary (pass 2)")
  .rtsa_check_converged_pass(lb, ub[nt], "the design-route calibration (pass 2)")

  za <- lb$za
  sentinel <- abs(za) == 20
  beta_ubound <- za
  beta_ubound[sentinel] <- NA_real_
  list(timing = t, alpha_ubound = ub, za = za, beta_ubound = beta_ubound,
       root = root2, rm_bs = as.integer(rm_bs), delta = delta,
       beta_spent = lb$as_cum, beta_spent_delta = lb$as_incr,
       sentinel = sentinel)
}

## RTSA::boundaries(timing = t_ext, ..., side = 2, futility = "non-binding",
##                  type = "analysis", design_R = design_R)
## `t_ext` must already end at design_R (RTSA() appends it before calling).
.rtsa_analysis_bounds <- function(t_ext, design_R, alpha, beta) {
  t_ext <- as.numeric(t_ext)
  if (!length(t_ext) || any(!is.finite(t_ext)) || any(t_ext <= 0))
    stop("information fractions must be finite and strictly positive")
  if (!is.finite(design_R) || design_R <= 0)
    stop("design_R must be a finite, strictly positive scalar")
  ## `t_ext` normally ends at design_R (as RTSA() builds it). RTSA's own
  ## boundaries(type = "analysis") also accepts a timing that stops short of
  ## design_R and appends it itself; that call shape is supported (0.2.7.17)
  ## so the frozen live-RTSA reference, produced that way, can be reproduced.
  if (length(t_ext) < 2L && !(max(t_ext) < design_R))
    stop("the analysis route needs at least one observed look plus design_R")
  .rtsa_check_look_spacing(t_ext, "the analysis-route look schedule")
  ab <- .rtsa_alpha_cpp(t_ext, side = 2L, alpha = alpha, design_R = design_R,
                       where = "the analysis-route alpha (efficacy) boundary")
  ub <- ab$alpha_ubound
  delta <- abs(stats::qnorm(alpha / 2) + stats::qnorm(beta))

  rm_bs <- 0L
  lb <- NULL
  for (pass in 1:3) {        # RTSA: three beta_boundary() calls
    lb <- .rtsa_beta_cpp(t_ext, ub, beta, delta, rm_bs = rm_bs,
                         design_R = design_R,
                         where = sprintf("the analysis-route beta (futility) boundary (pass %d)", pass))
    ## 0.2.7.15: an unreachable beta target at an INTERIOR look means the
    ## futility bound reaches the efficacy wall before the final look -- a
    ## genuine failure of this analysis pass. At the FINAL look it means the
    ## bound lies beyond the wall, which RTSA's own final clamp below
    ## (beta_ubound[last] > alpha_ubound[last] -> alpha_ubound[last]) already
    ## handles, so it is accepted.
    unr <- if (is.null(lb$unreachable_look)) 0L else lb$unreachable_look
    n_looks <- length(lb$za)   # = length(t_ext), or +1 when RTSA appends design_R
    if (unr > 0L && unr < n_looks)
      stop(sprintf(paste0("the futility bound reaches the efficacy wall at look ",
                          "%d of %d in the analysis-route beta calculation (pass %d); ",
                          "the analysis route cannot be computed for this schedule."),
                   unr, n_looks, pass), call. = FALSE)
    rm_bs <- as.integer(sum(lb$za < 0))
  }
  final_beyond_wall <- unr > 0L
  za <- lb$za
  sentinel <- abs(za) == 20
  beta_ubound <- za
  beta_ubound[sentinel] <- NA_real_
  nl <- length(ub)
  if (!is.na(beta_ubound[nl]) && beta_ubound[nl] > ub[nl])
    beta_ubound[nl] <- ub[nl]
  ## NB `za` is the RAW recursion output. When `final_beyond_wall` is TRUE its
  ## last element is the non-physical sentinel zb + 1 (see .rtsa_beta_cpp());
  ## only `beta_ubound` (RTSA's final clamp applied, sentinels -> NA) is
  ## meaningful as a boundary.
  list(timing = t_ext, alpha_ubound = ub, za = za, beta_ubound = beta_ubound,
       root = design_R, rm_bs = rm_bs, delta = delta,
       beta_spent = lb$as_cum, beta_spent_delta = lb$as_incr,
       sentinel = sentinel, final_beyond_wall = final_beyond_wall)
}

## RTSA::RTSA(type = "analysis", design = NULL, power_adj = TRUE) for a
## two-sided, non-binding-futility design: design pass -> design_R ->
## analysis pass.  `t_obs` = observed information / required information.
.rtsa_retrospective <- function(t_obs, alpha, beta) {
  t_obs <- as.numeric(t_obs)
  t_design <- t_obs[t_obs <= 1]
  if (!length(t_design))
    stop("the required information size is already reached at the first look")
  des <- .rtsa_design_bounds(t_design, alpha, beta)
  R <- des$root
  t_ext <- if (max(t_obs) < R) {
    c(t_obs, R)
  } else if (max(t_obs) > R) {
    c(t_obs[t_obs < R], R)
  } else {
    t_obs
  }
  ana <- .rtsa_analysis_bounds(t_ext, R, alpha, beta)
  list(design = des, analysis = ana, design_R = R)
}

## Estimated cumulative number of PARTICIPANTS at which the observed
## information first reaches `target` (an information fraction), by linear
## interpolation between the two bracketing looks; NA if never reached.
## (tsacor: the participants-based counterpart of tsahr's events-based helper of
## the same purpose; `cum_n` is the cumulative sum of n_subjects.) DARIS --
## target 1 -- and the route endpoint -- target route_endpoint -- are
## interpolated separately.
.tsacor_participants_at_fraction <- function(cumul_df, target) {
  if (!any(cumul_df$info_fraction >= target)) return(NA_real_)
  reach_idx <- which(cumul_df$info_fraction >= target)[1]
  if (reach_idx == 1) return(cumul_df$cum_n[1])
  f0 <- cumul_df$info_fraction[reach_idx - 1]
  f1 <- cumul_df$info_fraction[reach_idx]
  n0 <- cumul_df$cum_n[reach_idx - 1]
  n1 <- cumul_df$cum_n[reach_idx]
  w  <- if (f1 > f0) (target - f0) / (f1 - f0) else 0
  n0 + w * (n1 - n0)
}

## Definitive-look decision fields (0.2.7.14). `z` is the full cumulative Z
## vector; `boundary_upper` / `futility_upper` are the boundaries the formal
## decision compares against, one per look up to the definitive look
## (i.e. indexable by `final_tsa_look`). Everything refers to the DEFINITIVE
## look ONLY -- unlike crossed_tsa / entered_futility_region, which are
## "at ANY formal look". All three fields are NA when the route endpoint has
## not been reached (there is then no definitive look to report on).
.tsacor_definitive_look <- function(z, final_reached, final_tsa_look,
                                   boundary_upper, futility_upper) {
  na_out <- list(final_crossed_efficacy = NA,
                 final_entered_futility_region = NA,
                 final_non_efficacy = NA)
  if (!isTRUE(final_reached)) return(na_out)
  z_final <- abs(z[final_tsa_look])
  b_final <- boundary_upper[final_tsa_look]
  f_final <- futility_upper[final_tsa_look]
  crossed <- if (is.na(b_final) || is.na(z_final)) NA else (z_final >= b_final)
  entered <- if (is.na(f_final) || is.na(z_final)) NA else (z_final <= f_final)
  list(final_crossed_efficacy = crossed,
       final_entered_futility_region = entered,
       final_non_efficacy = if (is.na(crossed)) NA else !crossed)
}

## ** 0.2.7.15 ** A CONVERGED calibration pass (the design route's pass 1 /
## pass 2 at the accepted root) must never be "unreachable" -- the root search
## brackets a sign change of the final gap, so an unreachable final target at
## the accepted root would mean the search converged onto the edge of the
## infeasible region rather than onto a root. Also re-verifies the residual
## final gap directly (uniroot's tolerance is on the root, not on the gap).
.rtsa_check_converged_pass <- function(lb, final_wall, where, tol_gap = 1e-6) {
  if (isTRUE(lb$unreachable_look > 0L))
    stop(sprintf(paste0("%s reached an unreachable futility target at look %d: ",
                        "the information-scale root search converged onto the ",
                        "infeasible region, not onto a root."),
                 where, lb$unreachable_look), call. = FALSE)
  res_gap <- final_wall - lb$za[length(lb$za)]
  if (!is.finite(res_gap) || abs(res_gap) > tol_gap)
    stop(sprintf(paste0("%s: the final futility bound does not meet the final ",
                        "efficacy bound at the accepted root (residual gap %.3g)."),
                 where, res_gap), call. = FALSE)
  invisible(lb)
}

## ** 0.2.7.15 ** Look-spacing diagnostic. RTSA's own RTSA(type = "design")
## STOPS when any look adds less than 1% of the required information, and
## RTSA() drops such looks before computing boundaries; tsacor keeps every
## study. The 0.25% threshold below is NOT an RTSA rule: it is an empirical
## warning level measured for THIS implementation's default integration grid
## (r = 18 vs r = 72). That grid resolves increments down to ~0.5% of the
## required information to within ~1e-3 in the boundary and down to ~0.25% to
## within ~1e-5 for a single close pair; below that the recursion can become
## numerically unreliable (in the cases measured, grossly so: e.g. a final
## efficacy wall of 3-12 for a schedule whose true wall is about 2.1). The
## evidence is a set of measurements on even schedules and on a single close
## pair, not a systematic benchmark of the failure region. Warn -- do not alter
## the schedule.
.rtsa_check_look_spacing <- function(t, where, threshold = 0.0025) {
  t <- as.numeric(t)
  if (length(t) < 2L) return(invisible(t))
  ## Increments BETWEEN looks: the first look's own size is irrelevant to the
  ## grid accuracy (its density is evaluated directly, not through a kernel).
  dt <- diff(t)
  small <- which(dt < threshold)
  if (length(small)) {
    warning(sprintf(
      paste0("%s has %d look(s) that add less than %.2f%% of the required ",
             "information (smallest increment %.4g of the required ",
             "information, at look %d of %d). Below this spacing the ",
             "recursive integration at the default grid resolution can become ",
             "numerically unreliable (in the cases measured, sometimes grossly ",
             "so), so the boundaries should not be trusted. (The %.2f%% level is an ",
             "empirical warning threshold for this implementation's default ",
             "integration grid, not an RTSA rule; RTSA itself refuses such ",
             "schedules (design) or drops looks adding < 1%% (RTSA()).) ",
             "Consider merging near-simultaneous studies."),
      where, length(small), threshold * 100, min(dt), which.min(dt) + 1L, length(t),
      threshold * 100
    ), call. = FALSE)
  }
  invisible(t)
}

## Formatters for tsa_cor()'s summary_table `Value` column (0.2.7.22). Element by
## element on purpose: format() on a vector applies ONE common number of
## decimals to all of it (0.509 -> "0.509" but 30478 -> "30478.000").
.tsacor_format_numeric <- function(x) {
  x <- as.numeric(x)
  vapply(x, function(v) {
    if (is.na(v)) NA_character_ else format(v, scientific = FALSE, trim = TRUE, digits = 15)
  }, character(1))
}
.tsacor_format_logical <- function(x) as.character(as.logical(x))

