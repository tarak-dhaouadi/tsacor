> **Inherited document.** This file is carried over unchanged from the tsahr package (same maintainer) because tsacor reuses the same RTSA-derived boundary engine. Version numbers such as "0.2.7.13" and every mention of `tsahr`/`tsa_hr()` refer to the development history of tsahr; the engine itself works on information fractions only and is effect-measure agnostic.

# RTSA-matched boundary engines in tsahr

## Current implementation (0.2.7.11 - 0.2.7.18)

This section describes the code as shipped now. Everything under
"HISTORICAL" further down documents the earlier pure-R reconstruction
(0.2.4-0.2.7.10) and is kept only as a debugging record; where it contradicts
this section, this section is right.

**Architecture.** `tsa_hr()` -> `R/rtsa_engine.R` (RTSA's orchestration: which
root searches run, in which order, on which timeline; `stats::uniroot` as RTSA
uses it) -> `src/rtsa_engine.cpp` (Rcpp glue) -> `src/rtsa_core.h` (C++ port of
RTSA 0.2.2's `alpha_boundary()`, `beta_boundary()`, `z_n_w()`, `searchfunc()`,
`esOF()`, `sd_inf()` and the three used `first.cpp` functions `init_int()`,
`recur_int()`, `prob()`). The package needs a C++ toolchain
(`NeedsCompilation: yes`).

**What the engine does**

* *Efficacy (alpha):* RTSA's Simpson-grid recursion (`r = 18`), spending
  `esOF(alpha/side)`. The FFT engine in `R/obf_boundaries.R` is legacy.
* *Futility (beta), `boundary_route = "design"` (default):* RTSA's
  `boundaries(type = "design")` -- alpha bounds on the (t < 1, 1) timeline, then
  the two-pass information-scale root search (`rm_bs` suppresses the spend at
  early looks whose first-pass bound is negative). The wall the search targets
  is the alpha recursion's OWN value at t = 1 (2.127 for the 9-look schedule
  below), never `qnorm(1 - alpha/2)`. From 0.2.7.12 the final futility bound
  equals the final efficacy bound, as in RTSA. `qnorm(1 - alpha/2)` survives
  only in the legacy fallback engine.
* *`boundary_route = "analysis"`:* RTSA's `RTSA(type = "analysis", design =
  NULL)` chain (design pass -> `design_R` -> analysis pass on `t / design_R`
  with alpha respent on that timeline); the formal endpoint is `design_R *
  DARIS`, reported separately from DARIS itself (0.2.7.14).
* *Decision layer:* `crossed_tsa` / `entered_futility_region` are "at any
  formal look"; `final_crossed_efficacy`, `final_non_efficacy` and
  `final_entered_futility_region` are the definitive look only (0.2.7.14).

**Disclosed departures from RTSA, and safeguards** (none changes a converged,
ordinary result)

* `searchfunc()` has an iteration cap (RTSA's loop has none; RTSA-type runs
  converge in <= 9 rounds). Beyond the cap it throws unless the residual is
  within a loose tolerance (1e-6), in which case it is counted
  (`slow_searches`).
* **Unreachable beta target (0.2.7.15, analytic test 0.2.7.16).** When the spend
  a look asks for exceeds the probability mass still alive -- decided directly:
  target > `qmax = sum(last)`, the limit of the beta probability as the boundary
  goes to +Inf, tested directly up to the numerical tolerance -- the futility
  bound would sit beyond the efficacy wall.
  `beta_boundary()` reports this (`unreachable_look`) and returns `zb + 1` from
  that look on. That `zb + 1` is a non-physical SENTINEL, not a boundary: it only
  makes the final gap the root search works on negative -- the correct sign --
  instead of the candidate aborting. A reachable target that the search fails to
  find is a hard error. A
  *converged* calibration pass must never be unreachable, and its residual gap
  is re-verified (`.rtsa_check_converged_pass()`). In the analysis route an
  unreachable FINAL look is handled by RTSA's own final clamp; an interior one
  is an error. (0.2.7.13/14 threw here instead: the bracketing root search then
  failed for many schedules, including a real 37-look schedule, and `tsa_hr()`
  fell back to the legacy engine.)
* **`spend == beta` shortcut not ported (0.2.7.22).** RTSA sets `za = 0` when a
  look's spend equals `beta` bit for bit. With all interim looks suppressed that
  made the design calibration's gap constant ("no root bracket") for
  `beta = 1 - power` at power 0.80 / 0.95 (RTSA's literal `beta = 0.2` does not
  trigger it). Dropped at every look; unreachable in the frozen live reference.
* `z_n_w()`: a degenerate grid (zero width) is counted as `grid_collapses`; a
  REVERSED interval (lower wall above the upper wall) is counted separately as
  `grid_reversed` and warned about loudly. Neither throws, because the root
  searches evaluate `z_n_w()` at transient candidate scales; those candidate
  evaluations do not warn, the converged passes do.
* **Look spacing (0.2.7.15).** RTSA refuses design timings that add < 1% of the
  required information and drops such looks in `RTSA()`; tsahr keeps every study
  but warns (`.rtsa_check_look_spacing()`) when the increment between looks is
  < 0.25% of the required information, because the default grid then loses
  reliability (see the measurements below). The 0.25% level is an
  EMPIRICAL warning threshold measured for this implementation's default grid
  (r = 18 vs r = 72); it is not an RTSA rule.

**Robustness measurements (0.2.7.15, compiled core through the design-route
orchestration).** 300 random schedules of 3-45 looks: 186 calibrated with the
0.2.7.13/14 behaviour (measured on that older build, not shipped), 300 with the
current one (`tools/cpp_vs_py.py`, section 9, standalone -- no R). Schedules that failed before and
calibrate now include the real 37-look (40-study) schedule (root 1.2278, 14
suppressed looks, final wall 2.1703), 100 evenly spaced looks (root 1.2337, 39
suppressed), and a dense-at-end schedule. Grid accuracy of the final efficacy
wall, `r = 18` vs `r = 72`: 10, 40 and 100 even looks agree to <= 5e-5; 200
even looks (increment 0.005) differ by 6e-4; 500 even looks (0.002) give 11.7
vs 2.216 -- numerically unreliable, and with no other warning from the engine --
and a single extra look 0.001 or 0.0005 after another differs by 3e-4 / 6e-3.
Hence the 0.25% warning threshold. (These are measurements on even schedules
and on one close-pair family, not a systematic benchmark of the failure
region; the threshold is empirical.)

**Validation status.** Three independent-of-each-other lines of evidence:

1. *Live RTSA 0.2.2 (0.2.7.17).* `inst/extdata/rtsa_0.2.2_reference.R` freezes
   numbers printed by RTSA itself (`boundaries()`, timing 0.25/0.5/0.75/1,
   alpha 0.05, beta 0.2, non-binding, esOF; design route and RTSA's analysis
   call). The compiled core reproduces them: design root to 2e-16, alpha bounds
   to 2e-15, design beta bounds to 2e-16 (root search included); analysis alpha
   bounds to 2e-15 and analysis beta bounds exactly (difference 0.0).
   `tests/testthat/test-rtsa-live-reference.R` checks this without RTSA being
   installed; `tools/cpp_vs_py.py` section 10 checks it without R.
2. *Published vignette output.* RTSA's futility vignette (SMA timing
   0.541/0.812/1.083, futility 0.332/1.292/2.014) and alpha reference
   (4.877 3.357 2.680 2.290 2.031) are reproduced.
3. *Python port.* `tools/rtsa_port.py` (written from RTSA's R sources) agrees
   with the core to ~1e-14 and with the live reference above to the same
   precision -- so the reading of RTSA that both share is now confirmed on this
   case, not merely assumed.

What the frozen reference does NOT cover: many-look schedules, schedules with
`rm_bs` > 1, and the extended-timing analysis call that `RTSA()` itself uses
internally (timing ending at `design_R`). Those rest on evidence 2 and 3 (and on
the real-R runs of the maintainer's own 40-study data). RTSA's analysis call
shape frozen here (timing stopping short of `design_R`) makes `beta_boundary()`
append `design_R` itself, giving five beta bounds against four alpha bounds;
the C++ port accepts exactly that shape (alpha bound one element short, design-R
mode only, never read at the appended look).

### Evidence (0.2.7.11): what RTSA really computes, and where tsahr diverged

Evidence (a Python port of RTSA 0.2.2's R sources, written for this package
and checked against RTSA's published futility vignette and alpha reference;
the C++ core agrees with the port to ~1e-14, which shows faithful
translation of the port, not independent confirmation of the reading of RTSA), alpha = 0.05, beta = 0.20, looks
0.461092, 0.527034, 0.610745, 0.653952, 0.754067, 0.825699, 0.903255, 0.961096:

| route | first look ... last pre-1 look |
|---|---|
| RTSA `type = "design"` (root 1.2107, final wall 2.127) | 0.531 0.698 0.981 1.077 1.400 1.568 1.768 1.930 |
| RTSA `type = "analysis"` (design_R = 1.2107) | 0.126 0.304 0.588 0.684 1.002 1.163 1.344 1.456 |
| tsahr design pass with final wall forced to 1.96 (root 1.1517) | 0.479 0.642 0.921 1.015 1.334 1.498 1.690 1.832 |

* Both RTSA routes call the same `beta_boundary(side = 1)`; they differ in
  spending timeline (`t` vs `t / design_R`), information scale (`t * root`
  vs `t`) and in the alpha bounds the futility recursion is run against.
* The wall at the final look is the alpha recursion's value there, never
  `qnorm(1 - alpha/2)`.
* `first.cpp`'s `first`, `other`, `fcab`, `qpos`, `trap` are dead code in RTSA.
* Reproduce: `tools/rtsa_port.py`, `tools/cpp_vs_py.py`.

### 0.2.7.12 note
`tsa_hr()` uses the design route only (design pass applied to the observed
information fractions, DARIS = t = 1). RTSA's retrospective chain (design pass
-> design_R -> analysis pass at design_R, with alpha recomputed on t/design_R)
is implemented internally (`.rtsa_retrospective()`) but not used by `tsa_hr()`.
From 0.2.7.12 the final futility bound equals the final efficacy bound, as in
RTSA's design pass.

### 0.2.7.14 note
The compiled engine now distinguishes three numerical states in `z_n_w()`:
ordinary interval; degenerate (zero-width, za == zb) grid, counted as
`grid_collapses`; and REVERSED interval (za > zb), counted as `grid_reversed`
and reported by its own warning. None of them throws (the root searches
evaluate `z_n_w()` at transient candidate information scales), but only the
first is an ordinary RTSA computation; `tools/cpp_vs_py.py` checks the
classification. `searchfunc()` still throws on non-convergence beyond the loose
tolerance, except for the unreachable-beta-target case, which 0.2.7.15 turns into
a "beyond the wall" candidate (see "Current implementation").

## HISTORICAL: pure-R reconstruction of the boundary engines (0.2.4 - 0.2.7.10)

> Everything below documents how the earlier, R-only engines were built and
> debugged. It is NOT a description of the current implementation (see
> "Current implementation" above); some statements are explicitly marked
> superseded where they would otherwise mislead.

## Alpha-spending fix (0.2.6)

`tsahr` versions 0.2.0-0.2.5.1 computed the two-sided alpha-spending
boundaries with the more commonly-seen textbook two-argument Lan-DeMets
O'Brien-Fleming form, `2*(1-pnorm(qnorm(1-alpha/2)/sqrt(t)))`. That form
is internally valid as *an* alpha-spending function (it does reach
exactly `alpha` at t=1), and is the form quoted directly in some
general group-sequential-design references (e.g. gsDesign's
documentation) -- but it does not match the specific TSA methodology
(Copenhagen Trial Unit / RTSA, Thorlund et al.) that this package
documents itself as following (Miladinovic et al. 2013, Wetterslev
et al. 2009), and does not match RTSA's own `side`-parameterised form
evaluated at `side=2`: `2*(1-pnorm(qnorm(1-alpha/side/2)/sqrt(t)))*side`,
which at `side=2` reduces to `4*(1-pnorm(qnorm(1-alpha/4)/sqrt(t)))`.
0.2.6 corrects `.alpha_spend_OF()` to use this RTSA-matched form. This
was confirmed directly against a live call, across the full 5-look
schedule:

```
> RTSA::boundaries(timing=c(0.2,0.4,0.6,0.8,1), alpha=0.05, side=2,
                    es_alpha="esOF")
Upper: 4.877  3.357  2.680  2.290  2.031
```

All 5 boundaries match the corrected formula (run through tsahr's
existing recursive integration engine) essentially exactly, and are
clearly distinguishable from the old formula's 4.383, 3.099, 2.554,
2.254, 2.063 for the same schedule. A first-look regression test
(`test-alpha-spend-rtsa.R`) additionally pins the corrected formula
against this RTSA reference in closed form, since the first look has
no prior boundary to condition on and so needs no recursive-engine
approximation.

Both forms independently reach exactly `alpha` at t=1; they differ in
how alpha is allocated across interim looks -- the old form spent
roughly 9-11x more alpha at the first look of that schedule than RTSA's
boundaries call for, which works directly against the reason people
choose an O'Brien-Fleming design (strong protection against declaring
"significance" from sparse early evidence). Because both forms
satisfy `alpha*(1) = alpha` by construction, this discrepancy was not
caught by the existing total-alpha Monte Carlo/closed-form checks in
`R/obf_boundaries.R`, which check the *total* spend rather than the
*shape* across interim looks. This is a genuine correctness bug
relative to this package's own documented reference methodology, not a
stylistic or reference-preference choice -- anyone who used `tsa_hr()`
from a version between 0.2.0 and 0.2.5.1 for TSA boundaries at an
interim look should re-run their analysis with 0.2.6 or later. Only
`.alpha_spend_OF()` changed; the surrounding recursive integration
engine (the actual boundary-solving machinery) is untouched, as is the
beta/futility engine described below.

## Default grid resolution (0.2.6.1)

The "all 5 boundaries match essentially exactly" claim above was
re-checked with an independent from-scratch Python port of
`.obf_alpha_boundary()`, run against the nominal timing actually passed
to `RTSA::boundaries()` (0.2, 0.4, 0.6, 0.8, 1.0) rather than the
rounded `SMA_Timing` column it reports back (0.205, 0.409, ...) --
mixing those two up produces a spurious ~0.06 discrepancy at the first
look that has nothing to do with the formula or engine, since the
first look is an exact closed-form quantity determined entirely by the
timing value fed in. Using the correct nominal timing, the Python port
matched RTSA's boundaries to within ~0.003-0.006 across all 5 looks
even at this engine's long-standing default of `n_grid=2000` -- i.e.
"essentially exactly" holds up, and does not require the larger grid
introduced below to be true.

A previous draft of this note additionally claimed a specific,
more severe problem at n_grid=2000 -- boundaries of 4.877, 3.358,
2.703, 2.307, 2.064 against RTSA's 4.877, 3.357, 2.680, 2.290, 2.031
(~0.03 error at the final look), with a cited progression of
0.0325/0.0133/0.0074/0.0044/0.0020 as n_grid rose from 2000 to 32000.
That specific set of numbers could not be reproduced by the Python
port above and was never confirmed by an actual run of this package's
R code -- it is retracted here as unverified rather than repeated.
`.obf_alpha_boundary()`'s default `n_grid` is nonetheless still raised
from 2000 to 16000 in this release, on narrower grounds: FFT-based
convolution makes the extra grid resolution essentially free for this
package's actual usage pattern (a handful of interim looks per
`tsa_hr()` call), so there's no real cost to the more conservative
choice even without a confirmed problem at 2000. Anyone with a working
R installation is encouraged to run the snippet in the VALIDATION note
of `R/obf_boundaries.R` and confirm (or correct) the numbers above.

A 20,000-replicate Monte Carlo re-check (K=2, K=3 unequally spaced,
K=5, K=10) run directly against the corrected formula gave empirical
type-I error of 4.35%-4.93% against the 5% nominal target (Monte Carlo
SE ~=0.15%), consistent with correct behaviour. See the `VALIDATION`
note at the top of `R/obf_boundaries.R` for the fuller numeric account.
The full 5-look, engine-level comparison against live RTSA output has
been checked against the corrected formula specifically (see "Default
grid resolution" above); the one comparison that has only been run
under the pre-0.2.6 formula and not repeated for the corrected one is
the classical published ~2.040 O'Brien-Fleming constant check, which is
a different, narrower reference point (a single textbook number, not
RTSA's own multi-look output) -- that narrower gap does not need to
hold anyone up, since the direct RTSA comparison above is the stronger
and more relevant check for this package's specific goal of matching
RTSA rather than the general OF literature.

## Beta/futility boundary engine (added 0.2.4, RECONSTRUCTED in 0.2.7)

### What 0.2.4-0.2.6.9 actually shipped

`tsahr` 0.2.4 claimed to implement the retrospective/analysis-mode
`tsa_beta_bound` algorithm from RTSA, citing `getInnerWedge()` and the
helper functions it calls (`betas_Obf`, `sdfunc`, `first_old`,
`other_old`, `searchfunc_old`, `qpos_old`, `fcab_old`, `trap_old`, and
`gfunc`) as its source, described as RTSA's "old TSA functions
(translated from java)".

Having now obtained and read the actual RTSA package source (RTSA
0.2.2, https://github.com/AnneLyng/RTSA: `R/RTSA_helperfunctions.R`,
`R/boundaries.R`, `src/first.cpp`) directly: **none of those function
names exist anywhere in it.** `tsa_beta_bound`, `getInnerWedge`,
`betas_Obf`, `sdfunc`, `first_old`, `other_old`, `searchfunc_old`,
`qpos_old`, `fcab_old`, `gfunc`, `trap_old`, and `fakeIFY` (as a
standalone RTSA identifier) are simply not present. RTSA's real
non-binding-futility machinery is `beta_boundary()`, `z_n_w()`,
`searchfunc()` (R/RTSA_helperfunctions.R), and `init_int()`,
`recur_int()`, `prob()` (src/first.cpp) -- entirely different function
names, and, in several material respects described below, a different
algorithm from what 0.2.4-0.2.6.9 described and implemented. Whatever
the original source for the 0.2.4 write-up was, it was not this
package's own stated reference. 0.2.6.6-0.2.6.9 softened the file-level
comment from "literal R implementation" to "adapted from", on an
external audit's advice, but did not re-examine or rebuild the engine
itself -- the algorithm below is what 0.2.4-0.2.6.9 actually computed.

The old algorithm, for the record:
1. Compute cumulative O'Brien-Fleming beta spending
   (`2 * upper_tail(qnorm(1-beta/2) / sqrt(t))` -- this part was, and
   remains, correct: see below).
2. Convert it to incremental beta spending.
3. Build a **symmetric, null-referenced** inner wedge recursively (i.e.
   the upper and lower edges of the wedge were `+-za`, not the actual
   asymmetric efficacy boundary).
4. Propagate the surviving density through successive information
   increments using a custom trapezoidal numerical integration.
5. Iteratively solve each new wedge edge so its incremental probability
   equals the incremental beta spend.
6. Calculate an **empirical** drift, `testDrift = fakeIFY + abs(ya[last])`,
   from the shape of the wedge itself.
7. Return `za + sqrt(t) * testDrift` as the futility boundary -- i.e.
   shift the null-referenced wedge by this empirical drift after the
   fact.
8. If observed information exceeds 1, use `qnorm(1-alpha/2)` as
   `fakeIFY`, compute the wedge only for `t < 1`, and append the
   definitive conventional boundary; otherwise replace the final wedge
   value by `qnorm(1-alpha/2)` directly.
9. Suppress (set to `NA`) every non-positive resulting boundary.

### What real RTSA does, and what 0.2.7 reconstructs

RTSA's actual `beta_boundary()` (called from `boundaries()` with
`side = 1` and an explicit `delta` override even for a two-sided
design -- RTSA's own convention, not a simplification) differs in three
consequential ways:

1. **Fixed, not empirical, drift.** `delta <- abs(qnorm(alpha/side) +
   qnorm(beta))` is computed directly from alpha, beta, and side before
   any recursion runs. It does not depend on the shape of the futility
   construction at all.
2. **The true efficacy boundary as the wall throughout.** `zb <-
   alpha_boundaries$alpha_ubound` -- the already-computed, generally
   decreasing-then-flattening alpha (efficacy) boundary -- bounds the
   recursive integration grid (`z_n_w()`) at *every* look, not just
   after the fact. Each futility boundary is solved directly under the
   fixed alternative-hypothesis drift from (1), conditioning on not yet
   having crossed this real wall, via `searchfunc()` calling
   `recur_int()`/`prob()`/`init_int()` (`src/first.cpp`) to propagate
   and query the running density on a Simpson's-rule grid.
3. **NA means "negligible spend," not "non-positive."** RTSA pins `za`
   to a sentinel (`zninf = -20`) only when the incremental beta spend at
   a look underflows to (numerically) zero, and reports `NA` only for
   that sentinel (`abs(lb$za) == 20` in `boundaries.R`). A finite,
   possibly negative, futility boundary is a legitimate, informative
   result -- it says the trial would only be flagged for advisory
   futility if the cumulative Z-statistic had already crossed to the
   "wrong" side of the null by that point -- and is not suppressed.

0.2.7 replaces the entire wedge/testDrift construction with a
line-by-line port of RTSA's real `beta_boundary()`, `z_n_w()`, and
`searchfunc()` (ported near-verbatim, since they are already pure R),
and of `init_int()`, `recur_int()`, and `prob()` (ported from
`src/first.cpp` to pure R, since tsahr has no compiled-code
dependency), specialised to `es_beta = "esOF"` -- the only spending
family tsahr exposes, and the formula this package's beta-spending
function already correctly implemented (point 1 in the old algorithm
above was never the bug; everything built on top of it was).

> **[SUPERSEDED in 0.2.7.11]** tsahr now DOES have a compiled-code
> dependency: `src/rtsa_core.h` / `src/rtsa_engine.cpp` (Rcpp,
> `NeedsCompilation: yes`) port these functions -- and `alpha_boundary()`,
> `beta_boundary()`, `z_n_w()`, `searchfunc()`, `esOF()`, `sd_inf()` -- to
> C++. The pure-R port described here survives only as the legacy fallback.

Two deliberate, disclosed departures from a literal port, both kept
from previous versions of this package:

> **[SUPERSEDED in 0.2.7.11/0.2.7.12 -- this is no longer true of the
> current engine.]** The current design engine DOES port RTSA's
> information-scale root search (`.rtsa_design_bounds()`), and the final
> futility bound is the final value of the alpha recursion (equal to the final
> efficacy bound), not `qnorm(1-alpha/2)`. The convention described in this
> bullet survives only in the legacy R-only fallback (`.obf_beta_boundary()`
> and friends). See "Current implementation" above.
>
> * **Definitive final look.** RTSA's own `boundaries()` achieves an
  exact meeting of the alpha and beta boundaries at the final look via
  root-finding (`uniroot(inf_warp, ...)`, rescaling its own "design"
  timing until the two meet) -- machinery for a *design* problem (find
  the sample-size inflation hitting a target power for a not-yet-run
  trial) that has no counterpart in tsahr's retrospective,
  observed-data architecture. Rather than port that root-finding
  wholesale, 0.2.7 keeps this package's previous convention: the single
  definitive final look (t = 1) has no distinct non-binding futility
  zone apart from the main efficacy decision, so its boundary is set to
  the closed-form `qnorm(1-alpha/2)` directly.
* **Defensive `pmin()` against the efficacy boundary.** RTSA's own
  recursion does not hard-enforce `za[i] <= zb[i]` at every step (only
  the grid construction two steps back clips indirectly), so a
  pathological or internally-inconsistent alpha/beta/delta combination
  could in principle let the search push a futility value above the
  efficacy wall. tsahr applies `pmin(boundary, c_vec_alpha)` on top of
  the reconstructed recursion to make that impossible regardless --
  this was already present in 0.2.4-0.2.6.9 and is kept.

`t` values past 1 (observed studies accrued after DARIS/HARIS) still
have no counterpart in RTSA's own output at all (its timing vector is
capped at its design/root); as before, the single definitive t = 1
boundary is repeated for every t >= 1, purely for display continuity.

### The information-scale root search (0.2.7.4)

0.2.7-0.2.7.3 ported the per-look recursion above faithfully but missed
something upstream of it: RTSA's `boundaries()` never runs
`beta_boundary()` directly on the observed information fractions. It
always wraps it in a root search first.

Concretely, `beta_boundary()` keeps two different fraction scales
alive: `beta_timing` (== the raw information fractions, used to look up
how much beta has been spent at each look) and `org_inf_frac` (==
`inf_frac * warp_root`, used for the actual standard-deviation scale --
`info$sd_incr`, `info$sd_proc` -- the recursion runs on). Whenever
`warp_root != 1`, these are genuinely different vectors. RTSA's
`boundaries()` (side = 2, futility = "non-binding", type = "design")
finds `warp_root` via `uniroot(inf_warp, ...)`, where `inf_warp(x)` runs
the full recursion at `warp_root = x` and returns the gap between the
resulting final-look futility boundary and the fixed efficacy boundary
there -- i.e. it searches for the information-scale inflation under
which the two boundaries meet EXACTLY at the definitive final look.
`uniroot()` needs a bracket, and RTSA doesn't know one in advance, so it
slides a narrow window (`[start - step, start]`, `step = 0.02`, `start`
initially 0.95) upward on failure, up to 50 attempts.

It then repeats this once more: after the first root-found pass, any
look whose futility boundary came back negative has its beta-spend
fraction zeroed (`rm_bs = sum(lb$za < 0)`, `beta_timing <- c(rep(0,
rm_bs), beta_timing[-(1:rm_bs)])`) -- which routes it through the
"negligible spend" sentinel instead -- and the root is re-found (window
width 0.05 this time) under that adjustment. This second pass is *why*
RTSA's real output never shows a small negative early futility value:
it doesn't detect and hide one after the fact, it prevents one from
surviving in the first place.

0.2.7-0.2.7.3 implicitly always used `warp_root = 1` and never ran this
suppression pass, which is why the beta boundaries it produced were
systematically far from a live `RTSA::boundaries()` call across the
whole schedule (wrong information scale throughout), did not meet the
efficacy boundary at the final look (nothing was solving for that), and
occasionally showed a small negative value at an early look (nothing
was suppressing it). 0.2.7.4 ports both passes as `.rtsa2_inf_warp()`
and `.rtsa2_find_warp_root()`, and threads `warp_root`/`rm_bs` through
`.rtsa2_beta_boundary_core()`.

One interaction worth flagging: the 0.2.7.1 preventive clamp that keeps
`za[i]` a small margin below `zb[i]` (added to stop a Simpson-grid
collapse crash -- see below) is no longer applied to the *final* look,
because the root search needs `za[last]` to be able to reach, and
briefly cross, `zb[last]` as the candidate warp factor is varied, or
`uniroot()` can never bracket a sign change there at all. The clamp is
kept for every earlier look, where the original crash was actually
triggered.

### A formula claim that was checked and rejected (0.2.7.4)

A bug report proposed changing `.rtsa_beta_spend_OF()`'s
`qnorm(1 - beta / 2)` to `qnorm(1 - beta)`, on the grounds that RTSA's
own `esOF(alpha, timing)` is `2*(1-pnorm(qnorm(1-alpha)/sqrt(t)))` (no
`/2`). Checked directly against the actual RTSA 0.2.2 source
(`R/RTSA_helperfunctions.R`):

```r
esOF <- function(alpha, timing){
  ...
  as_cum[i] <- 2*(1 - pnorm(qnorm(1-alpha/2)/sqrt(timing[i]), mean = 0, sd = 1))
  ...
}
```

The `/2` is present in the real function. `.rtsa_beta_boundary()` calls
`beta_boundary()` with `side = 1` (RTSA's own convention for the
non-binding futility construction, even in a two-sided design -- see
above), so `beta_spend <- esOF(beta/side, beta_timing)` evaluates to
`esOF(beta, beta_timing)`, which is exactly
`.rtsa_beta_spend_OF()`'s existing `2*(1-pnorm(qnorm(1-beta/2)/sqrt(t)))`.
This was already correct before 0.2.7.4 and is unchanged by it. Applying
the proposed fix would have halved the effective beta-spend rate and
made the actual bug (the missing root search above) worse, not better --
a smaller effective beta-spend pushes every futility boundary even
further from RTSA's real output, which plausibly is what made the claim
seem plausible from the symptom alone.

### 0.2.7.5: two bugs found by actually executing the algorithm

0.2.7.4's root search shipped with a bug severe enough that it failed
on almost any realistic dataset (`R CMD check` showed the second root
search pass failing on nearly every test that called `tsa_hr()`). This
time, rather than reasoning through another fix from source alone, the
exact algorithm was re-implemented from scratch in Python
(numpy/scipy) and actually run in the environment available here, on
the specific schedule from the failing check log and several other
realistic designs. That surfaced two real bugs, both fixed in 0.2.7.5:

1. `.rtsa_beta_spend_OF()` rejected `t = 0` outright
   (`if (any(t <= 0)) stop(...)`), a guard predating the `rm_bs`
   mechanism, which deliberately zeroes early timing entries -- so
   every PASS-2 call with `rm_bs > 0` (i.e. most realistic designs)
   failed immediately on this function's own error, which then
   surfaced one level up, misleadingly, as "root search did not
   converge." The spending formula itself was always correct at
   `t = 0` (`qnorm(1-beta/2)/sqrt(0) = +Inf`, `pnorm(Inf,
   lower.tail = FALSE) = 0` -- exactly the right "no spend yet"
   answer); only the guard was wrong.
2. A single-observed-look schedule (every information fraction
   already `>= 1`) is a genuine mathematical degeneracy for the root
   search, not a bug to search harder for: at the one synthetic point
   `t = 1`, the incremental spend is always exactly `beta`, which
   fixes `za = 0` regardless of the warp factor, so the search
   objective is constant and can never bracket a root. Detected and
   handled as a special case now.

See NEWS.md's 0.2.7.5 entry for the full account, including which
schedules were checked in the Python re-implementation and confirmed
to converge post-fix (and fail pre-fix, confirming the diagnosis
rather than just the cure).

### Validation status (honest)

The alpha engine (above) was checked against a **live** call to
`RTSA::boundaries()` across a full 5-look schedule and matched
essentially exactly. This beta/futility reconstruction has not been --
no R interpreter has been available at any point while writing any of
this, so none of the actual R package code has been run, let alone
compared against RTSA's real numeric output. 0.2.7.5 raises the bar
somewhat: the underlying algorithm was independently re-implemented
and executed (in Python, not R) against several realistic schedules,
which is how the two bugs above were actually found rather than
guessed at -- but that is still not the same as running this R package
itself, and it is not a comparison against RTSA's real output. Before
relying on this for anything beyond an approximate, illustrative
futility band, run a direct comparison, e.g.:

```r
design <- RTSA::boundaries(timing = c(0.2,0.4,0.6,0.8,1), alpha = 0.05,
                            beta = 0.2, side = 2, futility = "non-binding",
                            es_alpha = "esOF", es_beta = "esOF")
design$beta_ubound
design$root
tsahr:::.rtsa_beta_boundary(c(0.2,0.4,0.6,0.8,1), alpha = 0.05, beta = 0.2,
                             c_vec_alpha = design$alpha_ubound)$boundary
```

As of 0.2.7.4, every entry should be directly comparable (not just the
first K-1 as in 0.2.7-0.2.7.3): the final entry is now expected to
match `design$beta_ubound`'s own final value too, since both are
solving the same root-matching condition, and `warp_root`/`root` should
agree as well.

### Functions removed in 0.2.7

Superseded by the reconstruction and deleted as dead code:
`.rtsa_gfunc`, `.rtsa_trap`, `.rtsa_fcab`, `.rtsa_qpos`, `.rtsa_first`,
`.rtsa_other`, `.rtsa_searchfunc_old`, `.rtsa_get_inner_wedge`.
`.rtsa_beta_spend_OF()` is unchanged (its formula was already correct).

The package deliberately retains its existing alpha/effect-size/DARIS
machinery; only the beta/futility engine has been replaced.
### What this reconstruction actually matches: `type = "design"`, used for an `analysis`-shaped purpose (fixed in 0.2.7.7, regressed in 0.2.7.9, corrected in 0.2.7.10)

Everything above in this section describes `.rtsa_beta_boundary()`,
which is a faithful port of RTSA's `boundaries(..., type = "design")`
branch -- the code path that solves a non-binding futility construction
from scratch given only `alpha`, `beta`, and a timing vector, by
root-finding its own information-scale inflation (`warp_root`). That
is RTSA's *prospective-design* math.

`tsa_hr()` itself, however, is always retrospective: it takes whatever
information fractions the included studies actually produced and
treats that sequence directly as the timing, with no prior design
phase of its own. Purpose-wise, that is RTSA's `type = "analysis"` use
case (retrospective monitoring of real accruing data). But
`type = "analysis"` has a hard prerequisite in RTSA's own source: a
`design_R` (a sample-size-inflation root) that was already solved by
an earlier, separate call -- see `boundaries()`'s `side == 2` /
`futility == "non-binding"` branch, the `else` arm under
`if(type == "design")` (which consumes `design_R`, it does not solve
it), and `beta_boundary()`'s own `if(!is.null(design_R))` block.
Before 0.2.7.7, `tsa_hr()` had no such prior design call to get a
`design_R` from, so its beta engine ran `type = "design"`'s math on
the observed timing every single call -- mechanically the wrong branch
for what is conceptually an `analysis` call, even though it was the
only branch available with the information `tsa_hr()` had.

0.2.7.7 introduced `.rtsa_beta_boundary_analysis()` to close this gap,
calling `.rtsa_beta_boundary()` (the `side = 2`, `futility =
"non-binding"`, two-pass `warp_root` design engine described above) to
manufacture `design_R`, then running a 3-pass `rm_bs`-re-derivation
loop around the beta calculation.

**0.2.7.9 regressed this**, on the strength of an external "live RTSA
0.2.2 reconstruction" claiming `delta`, `design_R`, and `rm_bs` were
all wrong (`2.801585` vs. `2.486475`; `1.151571` vs. `1.057434`; `5`
vs. `0`), and switched the calibration to a `side = 1`,
`futility = "none"`, `uniroot(right_power, ...)` power-only root
search (`boundaries.R` lines ~69-90) with `rm_bs` fixed at `0`.

**That external reconstruction queried the wrong branch of RTSA's
source**, and this was caught by going back to RTSA's own top-level
`RTSA()` wrapper (`R/RTSA.R`) rather than reading `boundaries()` in
isolation. The exact block that manufactures `design_R` when no
design object is supplied for a retrospective analysis is:

```r
bounds <-
  boundaries(
    timing = timing,
    alpha = alpha,
    beta = beta,
    side = side,
    futility = futility,
    es_alpha = es_alpha,
    es_beta = es_beta,
    type = "design"
  )
design_R <- bounds$root
```

`side = side` and `futility = futility` here are `RTSA()`'s *own*
top-level arguments -- i.e. whatever the caller passed to `RTSA()`
itself. For this package's design (a two-sided efficacy test with a
non-binding futility boundary), that is `side = 2`,
`futility = "non-binding"` -- **not** `side = 1`,
`futility = "none"`. The calibration call this package needs to
reproduce is therefore `boundaries(..., side = 2,
futility = "non-binding", type = "design")`, which is exactly
`.rtsa_beta_boundary()`'s own two-pass `warp_root` search -- the thing
0.2.7.7 already called, and 0.2.7.9 replaced with a calculation from a
genuinely different RTSA code path (the plain power/sample-size root
search for a design with *no* futility boundaries at all, which is not
what this package has).

As of 0.2.7.10, `.rtsa_beta_boundary_analysis()` reverts to 0.2.7.7's
`design_R` source (`.rtsa_beta_boundary()`, side = 2, unchanged) and
`delta` (`abs(qnorm(alpha/2) + qnorm(beta))`, side = 2, matching the
package's own alpha engine and the outer design's own sidedness), and
restores the 3-pass `rm_bs` fixed-point iteration around the beta
calculation -- RTSA's real `side == 2` / `futility == "non-binding"` /
`type == "analysis"` branch (`boundaries.R` lines ~436-492) calls
`beta_boundary()` three times, re-deriving `rm_bs` from the previous
pass's negative-`za` count each time, not once with `rm_bs` fixed at
`0` (that no-suppression, single-call pattern belongs to the
`side == 1` branch RTSA reaches for a *different* design than this
package's). The `side = 1`/`futility = "none"`/`right_power()`
functions 0.2.7.9 added (`.rtsa_design_R()`,
`.rtsa2_alpha_boundary_design_side1()`, `.rtsa2_esOF()`,
`.rtsa2_ma_power_upper()`, `.rtsa2_right_power()`) have been removed
rather than left in as dead, misleading code.

One concrete, checkable symptom of the 0.2.7.9 regression: RTSA's own
`boundaries()` converts every look whose `za` is pinned exactly at the
`+/-20` sentinel (i.e. suppressed by `rm_bs`) to `NA` before returning
(`beta_ubound <- c(rep(NA, sum(abs(za) == 20)), za[abs(za) < 20])`).
With `rm_bs` fixed at `0` (0.2.7.9), that sentinel was essentially
never hit, so early-look futility boundaries that should have been
blank/`NA` (matching RTSA's own printed output) instead surfaced as
small, genuine finite negative numbers in `tsa_hr()`'s printed table
and plot. This was never a separate display-layer bug to patch --
`tsa_hr()` has no NA-hiding logic of its own beyond what
`.rtsa_beta_boundary_analysis()` returns -- it was a direct consequence
of the wrong calibration, and disappears once `design_R`/`rm_bs` are
computed correctly.

A related question raised alongside the 0.2.7.9 regression was whether
tsahr's own R port of the recursive numerical integration
(`.rtsa2_init_int()`, `.rtsa2_recur_int()`, `.rtsa2_prob()`) might
itself be numerically diverging from RTSA's compiled
`src/first.cpp`/Rcpp equivalents (`init_int()`, `recur_int()`,
`prob()`), and whether tsahr needed its own compiled C++ to match
exactly. A line-by-line comparison against `src/first.cpp` (also
checking `RcppExports.cpp`'s registration table) shows this is *not*
the case: RTSA's own R code (`RTSA_helperfunctions.R`) only ever calls
three of `first.cpp`'s eight exported functions --
`init_int()`, `recur_int()`, and `prob()` -- the rest (`first()`,
`trap()`, `fcab()`, `other()`, `qpos()`) are unused/orphaned exports,
never called from any RTSA R function. All three functions that *are*
used call `R::dnorm()`/`R::pnorm()` internally, which is the exact
same underlying Rmath C library that R's own `stats::dnorm()`/
`stats::pnorm()` call -- there is no floating-point difference between
computing a normal density/CDF from R versus from C++ via Rcpp, since
both paths reach the identical compiled routine. `.rtsa2_init_int()`,
`.rtsa2_recur_int()`, and `.rtsa2_prob()` were checked term-by-term
against `init_int()`/`recur_int()`/`prob()` (accounting for the
0/1-indexing offset between C++ and R, and RTSA's `stdv` matrix column
convention -- column 1 = `sd_incr`, column 2 = `sd_proc`) and match
exactly. The "beta bounds too small" symptom that prompted the C++
question was fully explained by the `design_R`/`delta`/`rm_bs`
regression above; no compiled code was needed or added.

> **[SUPERSEDED in 0.2.7.11]** The conclusion that "no compiled code was
> needed" concerned the numerical *kernels* (which do match term for term).
> The gap that mattered was orchestration -- above all the final efficacy
> wall -- and once that was closed the whole recursion was moved to C++ for
> speed and to run RTSA's exact alpha recursion, not the FFT approximation.

Source reference:
https://github.com/AnneLyng/RTSA (RTSA 0.2.2) --
`R/RTSA_helperfunctions.R`, `R/boundaries.R`, `R/RTSA.R`, `src/first.cpp`
