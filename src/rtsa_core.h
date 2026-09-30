// Copyright (C) the RTSA authors (Anne Lyngholm Soerensen, Markus Harboe Olsen,
// Theis Lange, Christian Gluud) for the algorithms and code this file is derived
// from (RTSA 0.2.2, GPL (>= 2)); copyright (C) Tarak Dhaouadi for the
// adaptation. This file is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option) any
// later version. See DESCRIPTION and inst/COPYRIGHTS.
// -----------------------------------------------------------------------------
// NOTE (tsacor 0.1.0): this file is inherited unchanged from the tsahr package
// (same maintainer), which in turn ports RTSA; only package/registration names
// were changed. Names such as "tsa_cor()" / "tsacor" in comments stand for what were
// originally tsa_hr() / tsahr in that package's history. The engine works on
// information fractions only and is therefore effect-measure agnostic.
// -----------------------------------------------------------------------------
// rtsa_core.h -- C++ port of RTSA's recursive-integration boundary engine.
//
// Source of truth: RTSA 0.2.2 (Soerensen, Olsen, Lange, Gluud),
//   src/first.cpp            init_int(), recur_int(), prob()
//   R/RTSA_helperfunctions.R z_n_w(), searchfunc(), alpha_boundary(),
//                            beta_boundary(), esOF(), sd_inf()
// The three first.cpp functions are reproduced term-for-term (same operand
// order, same Rmath calls). first.cpp's other five exports (first, other, fcab,
// qpos, trap) are NEVER called from any R file of RTSA and are therefore not
// ported here (see NEWS.md, 0.2.7.11).
//
// The two boundary drivers (alpha_boundary(), beta_boundary()) are ported
// end-to-end -- spending function, Simpson grid, information scale -- so that
// the whole numerical pathway is executed in compiled code, exactly the way
// RTSA executes it in R + C++.  Root finding (uniroot) stays in R.
//
// Deliberate, disclosed differences from RTSA (none changes a converged
// result):
//   * searchfunc(): RTSA's `while(cond)` has no iteration cap; here it stops
//     after `max_outer` rounds (default 400; RTSA's own runs converge in <= 9).
//     ** 0.2.7.13: ** hitting that cap is no longer silently accepted. The
//     residual |qout - as| is checked against a loose tolerance
//     (kLooseSearchTol, 1e-6 -- far coarser than the 1e-9/1e-15 convergence
//     tolerances RTSA's own recursion is run at, so genuine floating-point
//     noise on an otherwise-converged search still passes): within it, the
//     value is accepted and counted as a "slow" search (Diagnostics::
//     slow_searches); beyond it, this throws std::runtime_error (an R error
//     via Rcpp) rather than returning a value RTSA would never have produced
//     (RTSA's own uncapped loop would simply never terminate there).
//   * z_n_w(): RTSA errors ("wrong sign in 'by'") when the grid collapses to
//     fewer than 2 nodes; here every such event is counted
//     (Diagnostics::grid_collapses, or grid_reversed when the interval is
//     reversed, hi < lo) rather than erroring, and an interval with no positive
//     width (!(hi > lo)) is widened to the minimal non-zero window
//     [lo, lo+1e-8]; a zero-width two-node grid is only counted, and a
//     collapsed grid whose interval already has positive width is used as is.
//     So a single
//     degenerate look does not abort an otherwise-computable boundary
//     sequence. The count is surfaced to the caller as an R warning (see
//     R/rtsa_engine.R) precisely because this is a real departure from RTSA.
//     ** 0.2.7.15: ** an UNREACHABLE beta target -- the futility bound would
//     have to sit beyond the efficacy wall because the spend asked for at a
//     look exceeds the probability mass still alive (0.2.7.16: decided
//     directly from the analytic limit, by comparing the target with
//     qmax = sum(last), the limit of the beta probability as the boundary
//     goes to +Inf, up to the numerical tolerance) -- is no longer
//     reported as a generic non-convergence error. A reachable target that
//     the search merely fails to find remains a hard error. beta_boundary()
//     catches it
//     (SearchUnreachable), sets BetaOut::unreachable_look and fills the
//     remaining futility bounds with zb + 1 ("beyond the wall"), so the final
//     gap zb[nn] - za[nn] the information-scale root searches work on comes
//     out NEGATIVE, which is exactly the right sign, instead of the whole
//     candidate evaluation aborting. (0.2.7.13/14 threw here, which made the
//     bracketing root search fail for many schedules -- e.g. 37-look and
//     100-look schedules -- and pushed tsa_cor() onto the legacy fallback;
//     0.2.7.12's silent cap-and-accept returned the right sign only by
//     accident and returned junk elsewhere.) A converged calibration pass must
//     never be unreachable; R checks that.
//     ** 0.2.7.14: ** a REVERSED interval (lower wall above the upper wall,
//     za > zb) is a different, more serious state than a merely degenerate
//     one: widening it to [za, za + 1e-8] turns an inverted boundary
//     configuration into something that looks like an ordinary computation.
//     It is therefore counted SEPARATELY (Diagnostics::grid_reversed) and
//     surfaced by its own, more alarming R warning. It deliberately does NOT
//     throw here: z_n_w() is also evaluated at every candidate information
//     scale tried by the root searches, where an interior look can be
//     transiently reversed at a bad candidate without the search (which only
//     targets the final look) being unable to continue past it and converge.
//   * Array LENGTHS are validated on entry (inf_frac, alpha_bound and the beta
//     timeline must agree, apart from RTSA's documented appended-look shape;
//     rm_bs must not exceed the number of looks), and reachability/convergence
//     failures are reported. The inner numerical loops still use unchecked
//     operator[] on vectors whose sizes were established by those checks;
//     RTSA's R code indexes past vector ends silently, which this port does
//     not rely on.
//   * ** 0.2.7.22: ** the exact-float shortcuts `spend == beta -> za = 0` that
//     RTSA has at the first look and at later looks are NOT ported (see
//     beta_boundary()).
//
// Compiled inside R (Rcpp) it calls R::dnorm/R::pnorm/R::qnorm -- the very
// same Rmath routines stats::dnorm/pnorm/qnorm and RTSA's C++ use.  With
// -DRTSA_STANDALONE it uses a self-contained fallback (tools/ only) so the
// numerics can be regression-tested without an R installation.
// -----------------------------------------------------------------------------
#ifndef TSACOR_RTSA_CORE_H
#define TSACOR_RTSA_CORE_H

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <string>

#ifndef RTSA_STANDALONE
#include <Rcpp.h>
#endif

namespace rtsa {

typedef std::vector<double> vec;

// ---------------------------------------------------------------- Rmath layer
#ifdef RTSA_STANDALONE
inline double zdens(double x, double mu, double sd) {
  const double z = (x - mu) / sd;
  return std::exp(-0.5 * z * z) / (sd * 2.5066282746310002);  // sqrt(2*pi)
}
inline double zcdf(double x, double mu, double sd) {
  return 0.5 * std::erfc(-(x - mu) / sd * 0.70710678118654752440);
}
// bisection on the (accurate-in-both-tails) erfc-based pnorm; test use only
inline double zquant(double p, double mu, double sd, bool lower) {
  if (!lower) p = 1.0 - p;   // adequate for the standalone tests (p not tiny)
  double lo = -40.0, hi = 40.0;
  for (int it = 0; it < 300; ++it) {
    double mid = 0.5 * (lo + hi);
    if (zcdf(mid, 0.0, 1.0) < p) lo = mid; else hi = mid;
  }
  return mu + sd * 0.5 * (lo + hi);
}
inline double zquant_upper(double p) {          // z with P(Z > z) = p
  double lo = -40.0, hi = 40.0;
  for (int it = 0; it < 300; ++it) {
    double mid = 0.5 * (lo + hi);
    if (0.5 * std::erfc(mid * 0.70710678118654752440) > p) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}
#else
inline double zdens(double x, double mu, double sd) { return R::dnorm(x, mu, sd, 0); }
inline double zcdf(double x, double mu, double sd) { return R::pnorm(x, mu, sd, 1, 0); }
inline double zquant(double p, double mu, double sd, bool lower) {
  return R::qnorm(p, mu, sd, lower ? 1 : 0, 0);
}
inline double zquant_upper(double p) { return R::qnorm(p, 0.0, 1.0, 0, 0); }
#endif

const double NaN = std::numeric_limits<double>::quiet_NaN();
inline bool is_nan(double x) { return x != x; }

// ------------------------------------------------------------- Diagnostics
// Counters threaded (by pointer, optional) through z_n_w()/searchfunc() and
// accumulated into AlphaOut/BetaOut, so the R layer can warn on anything
// that departs from RTSA's own (uncapped / non-widening) behaviour. Passing
// nullptr disables counting (used by the standalone/parity test harness,
// which wants the raw numerics without diagnostic bookkeeping).
// Thrown by searchfunc() (beta searches only) when the search cannot converge
// because the target spend exceeds the probability mass that can still be
// reached at that look (qout stays below `as`).
struct SearchUnreachable : public std::runtime_error {
  explicit SearchUnreachable(const std::string& m) : std::runtime_error(m) {}
};

struct Diagnostics {
  int grid_collapses = 0;   // degenerate (zero-width / not-reversed) grids
  int grid_reversed = 0;    // 0.2.7.14: REVERSED grids (lower wall above upper)
  int slow_searches = 0;
};

// How far a non-converged searchfunc() residual may be from the target
// spend before it is treated as a genuine failure rather than harmless
// floating-point noise on an otherwise-converged search. Deliberately far
// coarser than the 1e-9 (alpha) / 1e-15 (beta) convergence tolerances the
// recursion itself is run at.
const double kLooseSearchTol = 1e-6;

// ------------------------------------------------------------ spending / info
// RTSA esOF(): Lan-DeMets O'Brien-Fleming-type spending.  NB the `1 - zcdf()`
// form (not lower.tail = FALSE) is RTSA's own and is kept for exact parity.
inline void esOF(double a, const vec& timing, vec& as_cum, vec& as_incr) {
  const size_t n = timing.size();
  as_cum.assign(n, 0.0);
  as_incr.assign(n, 0.0);
  const double q = zquant(1.0 - a / 2.0, 0.0, 1.0, true);
  for (size_t i = 0; i < n; ++i) {
    as_cum[i] = 2.0 * (1.0 - zcdf(q / std::sqrt(timing[i]), 0.0, 1.0));
    as_incr[i] = (i == 0) ? as_cum[i] : as_cum[i] - as_cum[i - 1];
  }
}

// RTSA sd_inf(): sd of the increments and of the process
inline void sd_inf(const vec& timing, vec& sd_incr, vec& sd_proc) {
  const size_t n = timing.size();
  sd_incr.assign(n, 0.0);
  sd_proc.assign(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    const double prev = (i == 0) ? 0.0 : timing[i - 1];
    sd_incr[i] = std::sqrt(timing[i] - prev);
    sd_proc[i] = std::sqrt(timing[i]);
  }
}

// -------------------------------------------------------------------- z_n_w()
struct Grid { vec zj, wj; };

// `i` is the 1-based look index used throughout RTSA's R code.
inline Grid z_n_w(int r, const vec& sd_incr, const vec& za, const vec& zb,
                  int i, double delta, Diagnostics* diag = nullptr) {
  const int n0 = 6 * r - 1;
  vec xi(n0);
  for (int j = 1; j <= n0; ++j) {
    double v = delta * sd_incr[i - 1];
    if (j < r)                 v += (-3.0 - 4.0 * std::log((double)r / (double)j));
    if (r <= j && j <= 5 * r)  v += (-3.0 + 3.0 * (double)(j - r) / (2.0 * (double)r));
    if (5 * r < j)             v += (3.0 + 4.0 * std::log((double)r / (double)(6 * r - j)));
    xi[j - 1] = v;
  }
  // lower trim
  int last_below = -1;
  for (int k = 0; k < n0; ++k) if (xi[k] < za[i - 1]) last_below = k;
  if (last_below >= 0) {
    xi.erase(xi.begin(), xi.begin() + last_below);
    xi[0] = za[i - 1];
  }
  // upper trim
  int first_above = -1;
  for (size_t k = 0; k < xi.size(); ++k) if (xi[k] > zb[i - 1]) { first_above = (int)k; break; }
  if (first_above >= 0) {
    xi.resize(first_above + 1);
    xi[first_above] = zb[i - 1];
  }
  // grid-collapse safety net (RTSA would throw here)
  if (xi.size() < 2) {
    double lo = za[i - 1], hi = zb[i - 1];
    if (hi < lo) {
      // 0.2.7.14: reversed interval (za > zb) -- counted on its own, NOT
      // folded into grid_collapses; see the header note above.
      if (diag) diag->grid_reversed += 1;
    } else {
      if (diag) diag->grid_collapses += 1;
    }
    if (!(hi > lo)) hi = lo + 1e-8;
    xi.assign(2, 0.0);
    xi[0] = lo; xi[1] = hi;
  } else if (!(xi.back() > xi.front())) {
    // zero-width two-node grid (za == zb): kept as-is (RTSA arithmetic
    // unchanged) but now visible in the diagnostics.
    if (diag) diag->grid_collapses += 1;
  }
  const int n = (int)xi.size();
  const int m = 2 * n - 1;
  Grid g;
  g.zj.assign(m, 0.0);
  g.wj.assign(m, 0.0);
  for (int k = 0; k < n; ++k) g.zj[2 * k] = xi[k];
  for (int k = 0; k < n - 1; ++k) g.zj[2 * k + 1] = (xi[k] + xi[k + 1]) / 2.0;
  const vec& zj = g.zj;
  for (int k = 1; k <= m; ++k) {            // 1-based k as in R
    if (k == 1) {
      g.wj[k - 1] = (1.0 / 6.0) * (zj[2] - zj[0]);
    } else if ((k % 2 == 1) && k >= 3 && k <= m - 2) {
      g.wj[k - 1] = (1.0 / 6.0) * (zj[k + 1] - zj[k - 3]);
    } else if ((k % 2 == 0) && k >= 2 && k <= m - 1) {
      g.wj[k - 1] = (4.0 / 6.0) * (zj[k] - zj[k - 2]);
    } else {
      g.wj[k - 1] = (1.0 / 6.0) * (zj[m - 1] - zj[m - 3]);
    }
  }
  return g;
}

// ----------------------------------------------- first.cpp: init_int/recur_int/prob
inline vec init_int(const vec& wj, const vec& zj, double delta, double sd1) {
  const size_t n = zj.size();
  vec last(n);
  for (size_t i = 0; i < n; ++i) last[i] = wj[i] * (zdens(zj[i], delta * sd1, 1.0));
  return last;
}

// k is the 1-based index of the NEW look (RTSA's `k`); stdv(k-1,1) == sd_proc[k-1] etc.
inline vec recur_int(int k, const vec& sd_incr, const vec& sd_proc, const vec& zj,
                     const vec& last, const vec& zj_up, const vec& wj_up,
                     double delta, bool bs) {
  const size_t nu = zj_up.size(), nl = last.size();
  vec last_up(nu, 0.0);
  const double sdk = sd_incr[k - 1];       // stdv(k-1,0)
  const double spk = sd_proc[k - 1];       // stdv(k-1,1)
  const double spp = sd_proc[k - 2];       // stdv(k-2,1)
  for (size_t i = 0; i < nu; ++i) {
    for (size_t j = 0; j < nl; ++j) {
      if (bs) {
        last_up[i] += last[j] * spk / sdk *
          zdens((zj[j] * spp - zj_up[i] * spk) / sdk, delta * sdk, 1.0);
      } else {
        last_up[i] += last[j] * spk / sdk *
          zdens((zj_up[i] * spk - zj[j] * spp) / sdk, delta * sdk, 1.0);
      }
    }
    last_up[i] *= wj_up[i];
  }
  return last_up;
}

inline double prob(double xq, const vec& last, const vec& zj, int k,
                   const vec& sd_incr, const vec& sd_proc, bool bs, double delta) {
  double p_out = 0.0;
  const double sdk = sd_incr[k - 1];
  const double spp = sd_proc[k - 2];
  for (size_t i = 0; i < zj.size(); ++i) {
    if (!bs && delta != 0.0) {
      p_out += last[i] * (zcdf((zj[i] * spp - xq) / sdk, -delta * sdk, 1.0));
    } else if (bs) {
      p_out += last[i] * (zcdf((xq - zj[i] * spp) / sdk, delta * sdk, 1.0));
    } else {
      p_out += last[i] * (zcdf((zj[i] * spp - xq) / sdk, delta * sdk, 1.0));
    }
  }
  return p_out;
}

// ---------------------------------------------------------------- searchfunc()
inline double searchfunc(const vec& last, const vec& zj, int i, double as,
                         const vec& sd_incr, const vec& sd_proc,
                         const vec& za, const vec& zb, double tol, bool bs,
                         double delta, Diagnostics* diag = nullptr,
                         int max_outer = 400) {
  const int maxnn = 50;
  // ** 0.2.7.16: analytic reachability test for a BETA search. **
  // As the boundary xq -> +Inf every Phi(...) in prob() tends to 1, so the
  // largest spend any boundary can produce at this look is
  // qmax = sum(last), the probability mass still alive. If the target spend
  // exceeds it (beyond the summation tolerance), NO boundary value can reach
  // the target: the futility bound would lie beyond the efficacy wall. This is
  // decided here, from that analytic limit and up to the numerical tolerance
  // (`qmax < as - tol`, tol = 1e-15 for beta searches, plus floating-point
  // summation error), before any iterating -- not inferred afterwards
  // from where an unconverged search happened to stop. (0.2.7.15 inferred it
  // from `qout < as` at the iteration cap, which is not a proof: a reachable
  // target the search merely failed to find would also end below the target.)
  if (bs) {
    double qmax = 0.0;
    for (size_t q = 0; q < last.size(); ++q) qmax += last[q];
    if (qmax < as - tol) {
      throw SearchUnreachable(
        "RTSA-ported beta search: the target spend (" + std::to_string(as) +
        ") exceeds the probability mass still alive (" + std::to_string(qmax) +
        ") at look " + std::to_string(i) +
        " (the futility bound would lie beyond the efficacy wall)");
    }
  }
  double upper = zb[i - 2] * sd_proc[i - 1];
  if (bs) upper = za[i - 2] * sd_proc[i - 1];
  double del = 10.0;
  double qout = prob(upper, last, zj, i, sd_incr, sd_proc, bs, delta);
  int outer = 0;
  bool converged = false;
  for (;;) {
    ++outer;
    if (std::fabs(qout - as) <= tol) { converged = true; break; }
    if (qout > as + tol) {
      del = del / 10.0;
      for (int k = 0; k < maxnn; ++k) {
        if (bs) upper = upper - 2.0 * del;
        upper = upper + del;
        qout = prob(upper, last, zj, i, sd_incr, sd_proc, bs, delta);
        if (qout <= as + tol) break;
      }
    }
    if (qout < as - tol) {
      del = del / 10.0;
      for (int k = 0; k < maxnn; ++k) {
        if (bs) upper = upper + 2.0 * del;
        upper = upper - del;
        qout = prob(upper, last, zj, i, sd_incr, sd_proc, bs, delta);
        if (qout >= as - tol) break;
      }
    }
    if (outer > max_outer) break;           // disclosed tsacor-side safety cap
  }
  // ** 0.2.7.13: ** the cap above is RTSA's uncapped while(cond) loop made
  // finite; what it means to STOP there without having converged is new
  // territory RTSA itself never visits. Accept silently only within a loose
  // tolerance (harmless floating-point residue on what is, for all
  // practical purposes, a converged search); beyond it, this is a genuine
  // failure to find the boundary value RTSA's own algorithm defines, and is
  // reported as such rather than returned as if it were a real answer.
  if (!converged) {
    const double resid = std::fabs(qout - as);
    if (resid > kLooseSearchTol) {
      // Reaching this point means the target IS reachable (for beta searches
      // qmax >= as - tol was verified above; alpha searches have no such
      // ceiling) but the search failed to find it: a genuine numerical
      // failure, never to be passed off as "beyond the wall".
      throw std::runtime_error(
        "RTSA-ported boundary search did not converge at look " +
        std::to_string(i) + " after " + std::to_string(max_outer) +
        " rounds (residual " + std::to_string(resid) +
        " exceeds the loose tolerance " + std::to_string(kLooseSearchTol) +
        "). This boundary sequence cannot be trusted as computed; consider a "
        "different design (information fractions, alpha, power) or report "
        "the problem.");
    }
    if (diag) diag->slow_searches += 1;
  }
  return upper / sd_proc[i - 1];
}

// ------------------------------------------------------------- alpha_boundary()
struct AlphaOut { vec zb, as_incr, as_cum; int grid_collapses = 0; int grid_reversed = 0; int slow_searches = 0; };

// design_R = NaN  -> type "design";   design_R finite -> type "analysis"
// (alpha bounds are scale-invariant in theory; the analysis info scale
//  inf_frac * design_R is kept anyway so the arithmetic is RTSA's).
inline AlphaOut alpha_boundary(const vec& inf_frac, int side, double alpha,
                               double design_R, double tol = 1e-9, int r = 18,
                               Diagnostics* diag = nullptr) {
  const double zninf = -20.0;
  const size_t nn = inf_frac.size();
  if (nn < 1) throw std::invalid_argument("inf_frac must not be empty");
  vec alpha_timing = inf_frac;
  const double mx = *std::max_element(inf_frac.begin(), inf_frac.end());
  if (mx > 1.0) for (size_t i = 0; i < nn; ++i) alpha_timing[i] = inf_frac[i] / mx;

  AlphaOut out;
  esOF(alpha / side, alpha_timing, out.as_cum, out.as_incr);

  vec scale = inf_frac;
  if (!is_nan(design_R)) for (size_t i = 0; i < nn; ++i) scale[i] = inf_frac[i] * design_R;
  vec sd_incr, sd_proc;
  sd_inf(scale, sd_incr, sd_proc);

  out.as_incr[0] = std::min(alpha, out.as_incr[0]);
  out.as_incr[0] = std::max(0.0, out.as_incr[0]);

  vec za(nn, 0.0), zb(nn, 0.0);
  zb[0] = (out.as_incr[0] < tol) ? -zninf : zquant_upper(out.as_incr[0]);
  if (side == 1) za[0] = zninf; else za[0] = -zb[0];

  Grid g = z_n_w(r, sd_incr, za, zb, 1, 0.0, diag);
  vec last;
  for (size_t ii = 2; ii <= nn; ++ii) {
    const int i = (int)ii;
    if (i == 2) last = init_int(g.wj, g.zj, 0.0, sd_incr[0]);
    double a = out.as_incr[i - 1];
    if (a <= 0.0 || a >= 1.0) {
      a = std::max(0.0, std::min(1.0, a));
      out.as_incr[i - 1] = a;
    }
    if (a < tol) {
      zb[i - 1] = -zninf;
    } else if (a == 1.0) {
      zb[i - 1] = 0.0;
    } else {
      zb[i - 1] = searchfunc(last, g.zj, i, a, sd_incr, sd_proc, za, zb, tol,
                             false, 0.0, diag);
    }
    za[i - 1] = (side == 1) ? zninf : -zb[i - 1];
    if (i != (int)nn) {
      Grid up = z_n_w(r, sd_incr, za, zb, i, 0.0, diag);
      last = recur_int(i, sd_incr, sd_proc, g.zj, last, up.zj, up.wj, 0.0, false);
      g = up;
    }
  }
  out.zb = zb;
  if (diag) { out.grid_collapses = diag->grid_collapses; out.grid_reversed = diag->grid_reversed; out.slow_searches = diag->slow_searches; }
  return out;
}

// -------------------------------------------------------------- beta_boundary()
struct BetaOut {
  vec za, zb, as_incr, as_cum, beta_timing, info_frac_used, sd_incr, sd_proc;
  int grid_collapses = 0; int grid_reversed = 0; int slow_searches = 0;
  int unreachable_look = 0;   // 0.2.7.15: 1-based look where the beta target was unreachable (0 = none)
  // ** za is NOT a valid boundary from `unreachable_look` onward. ** Those
  // entries are the deliberately non-physical sentinel zb + 1: they exist
  // only to encode "beyond the efficacy wall" so that the sign of the final
  // gap zb[nn] - za[nn] (what the information-scale root searches use) comes
  // out negative. They must never be reported, plotted or interpreted as
  // futility bounds. R discards every such candidate; a converged calibration
  // pass with unreachable_look > 0 is an error (R/rtsa_engine.R).
};

// side: the `side` argument RTSA's boundaries() hands to beta_boundary()
//       (always 1 in the non-binding branches); delta must be supplied.
// warp_root / design_R: NaN when RTSA passes NULL.
inline BetaOut beta_boundary(const vec& inf_frac, const vec& alpha_bound,
                             double beta, int side, double delta, int rm_bs,
                             double design_R, double warp_root,
                             double zninf = -20.0, double tol = 1e-15, int r = 18,
                             Diagnostics* diag = nullptr) {
  size_t nn = inf_frac.size();
  if (nn < 1) throw std::invalid_argument("inf_frac must not be empty");

  vec org = inf_frac;
  if (!is_nan(warp_root)) for (size_t i = 0; i < nn; ++i) org[i] = inf_frac[i] * warp_root;

  vec beta_timing = inf_frac;
  if (!is_nan(design_R)) {
    beta_timing.assign(nn, 0.0);
    for (size_t i = 0; i < nn; ++i) beta_timing[i] = inf_frac[i] / design_R;
    org = inf_frac;
    const double mo = *std::max_element(org.begin(), org.end());
    if (mo < design_R) org.push_back(design_R);
    vec bt;
    for (size_t i = 0; i < beta_timing.size(); ++i) if (beta_timing[i] < 1.0) bt.push_back(beta_timing[i]);
    bt.push_back(1.0);
    beta_timing = bt;
    nn = beta_timing.size();
  }
  bool any_gt1 = false;
  for (size_t i = 0; i < beta_timing.size(); ++i) if (beta_timing[i] > 1.0) any_gt1 = true;
  if (any_gt1) {
    vec bt;
    for (size_t i = 0; i < beta_timing.size(); ++i) if (beta_timing[i] < 1.0) bt.push_back(beta_timing[i]);
    bt.push_back(1.0);
    beta_timing = bt;
    nn = beta_timing.size();
  }
  if (rm_bs != 0) {
    if ((size_t)rm_bs > beta_timing.size()) throw std::invalid_argument("rm_bs exceeds number of looks");
    for (int i = 0; i < rm_bs; ++i) beta_timing[i] = 0.0;
  }
  // ** 0.2.7.17 -- RTSA's unextended analysis call. ** RTSA::boundaries(type =
  // "analysis", design_R) may be given a timing that does NOT end at design_R;
  // beta_boundary() then appends design_R to the information scale itself
  // (so nn = looks + 1) while the alpha bounds it is handed cover only the
  // looks (RTSA itself recycles them with an R warning that is harmless because
  // zb is never read at the appended look). That exact shape -- alpha_bound one
  // element shorter than the beta timeline, design-R mode only -- is accepted
  // and padded with NaN at the appended look, which is never read. Any other
  // length mismatch is still an error.
  const bool zb_short = !is_nan(design_R) && alpha_bound.size() + 1 == nn;
  if (org.size() != nn || (alpha_bound.size() != nn && !zb_short))
    throw std::invalid_argument("inf_frac, alpha_bound and the beta timeline must have equal length");

  BetaOut out;
  esOF(beta / side, beta_timing, out.as_cum, out.as_incr);
  vec sd_incr, sd_proc;
  sd_inf(org, sd_incr, sd_proc);

  if (out.as_incr[0] <= 0.0 || out.as_incr[0] >= beta) {
    out.as_incr[0] = std::min(beta, out.as_incr[0]);
    out.as_incr[0] = std::max(0.0, out.as_incr[0]);
  }

  vec za(nn, 0.0);
  vec zb = alpha_bound;
  if (zb_short) zb.push_back(NaN);   // appended design_R look: never read
  // ** 0.2.7.22 -- deliberate departure from RTSA: no `spend == beta` shortcut.
  // RTSA has `else if (as_incr[i] == beta) za[i] <- 0` at every look. When all
  // interim looks are suppressed (rm_bs = nt - 1: the low-information,
  // "evidence still insufficient" case) the final look carries the whole beta
  // spend, and whether that spend equals `beta` bit for bit decides whether the
  // shortcut fires. tsacor computes beta = 1 - power (0.19999999999999996 for
  // power 0.8), which equals RTSA's own spend arithmetic in the last bit,
  // whereas RTSA is normally called with a literal beta = 0.2 that does not. When
  // it fires, za is 0 whatever the information scale, so the calibration gap
  // zb[nn] - za[nn] the root search works on is constant and can never change
  // sign ("no root bracket"), and a design that calibrates at power 0.8 + 1e-12
  // fails at 0.8. za = 0 has no statistical meaning here; the searched value
  // (qnorm at the first look, searchfunc() later) is the correct one, so the
  // shortcut is dropped. It is unreachable in every case the frozen RTSA
  // reference covers, so parity there is unchanged.
  if (out.as_incr[0] == 0.0) {
    za[0] = zninf;
  } else {
    za[0] = zquant(out.as_incr[0], sd_proc[0] * delta, 1.0, true);
  }

  Grid g = z_n_w(r, sd_incr, za, zb, 1, delta, diag);
  vec last;
  for (size_t ii = 2; ii <= nn; ++ii) {
    const int i = (int)ii;
    if (i == 2) last = init_int(g.wj, g.zj, delta, sd_incr[0]);
    double a = out.as_incr[i - 1];
    if (a <= 0.0 || a >= 1.0) {
      a = std::max(0.0, std::min(1.0, a));
      out.as_incr[i - 1] = a;
    }
    if (a < tol) {
      za[i - 1] = zninf;
    } else {           // 0.2.7.22: RTSA's `a == beta -> za = 0` shortcut dropped (see above)
      try {
        za[i - 1] = searchfunc(last, g.zj, i, a, sd_incr, sd_proc, za, zb, tol,
                               true, delta, diag);
      } catch (const SearchUnreachable&) {
        // 0.2.7.15: candidate is "beyond the wall" from this look on. The
        // zb + 1 written here is a SENTINEL, not a boundary (see BetaOut).
        out.unreachable_look = i;
        for (size_t j = (size_t)(i - 1); j < nn; ++j) {
          // the padded (appended-look) wall is NaN: use the last real wall
          const double wall = is_nan(zb[j]) ? zb[nn - 2] : zb[j];
          za[j] = wall + 1.0;
        }
        break;
      }
    }
    if (i != (int)nn) {
      Grid up = z_n_w(r, sd_incr, za, zb, i, delta, diag);
      last = recur_int(i, sd_incr, sd_proc, g.zj, last, up.zj, up.wj, delta, false);
      g = up;
    }
  }
  out.za = za;
  out.zb = zb;
  out.beta_timing = beta_timing;
  out.info_frac_used = org;
  out.sd_incr = sd_incr;
  out.sd_proc = sd_proc;
  if (diag) { out.grid_collapses = diag->grid_collapses; out.grid_reversed = diag->grid_reversed; out.slow_searches = diag->slow_searches; }
  return out;
}

}  // namespace rtsa
#endif  // TSACOR_RTSA_CORE_H
