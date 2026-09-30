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
// rtsa_engine.cpp -- Rcpp glue for the RTSA boundary engine (rtsa_core.h).
//
// Exports:
//   rtsa_alpha_boundary_cpp()  RTSA alpha_boundary()  (esOF spending)
//   rtsa_beta_boundary_cpp()   RTSA beta_boundary()   (esOF spending; design
//                              mode via warp_root, analysis mode via design_R)
//   rtsa_init_int_cpp(), rtsa_recur_int_cpp(), rtsa_prob_cpp()
//                              signature-compatible copies of the three
//                              RTSA:::init_int / recur_int / prob functions
//                              (parity tests call both side by side)
// -----------------------------------------------------------------------------
#include <Rcpp.h>
#include "rtsa_core.h"
using namespace Rcpp;

static rtsa::vec as_vec(const NumericVector& x) {
  return rtsa::vec(x.begin(), x.end());
}

// [[Rcpp::export]]
List rtsa_alpha_boundary_cpp(NumericVector inf_frac, int side, double alpha,
                             double design_R, double tol, int r) {
  rtsa::Diagnostics diag;
  rtsa::AlphaOut a = rtsa::alpha_boundary(as_vec(inf_frac), side, alpha,
                                          design_R, tol, r, &diag);
  return List::create(_["alpha_ubound"]   = wrap(a.zb),
                      _["as_incr"]        = wrap(a.as_incr),
                      _["as_cum"]         = wrap(a.as_cum),
                      _["grid_collapses"] = a.grid_collapses,
                      _["grid_reversed"]  = a.grid_reversed,
                      _["slow_searches"]  = a.slow_searches);
}

// [[Rcpp::export]]
List rtsa_beta_boundary_cpp(NumericVector inf_frac, NumericVector alpha_bound,
                            double beta, int side, double delta, int rm_bs,
                            double design_R, double warp_root, double zninf,
                            double tol, int r) {
  rtsa::Diagnostics diag;
  rtsa::BetaOut b = rtsa::beta_boundary(as_vec(inf_frac), as_vec(alpha_bound),
                                        beta, side, delta, rm_bs, design_R,
                                        warp_root, zninf, tol, r, &diag);
  return List::create(_["za"]             = wrap(b.za),
                      _["zb"]             = wrap(b.zb),
                      _["as_incr"]        = wrap(b.as_incr),
                      _["as_cum"]         = wrap(b.as_cum),
                      _["beta_timing"]    = wrap(b.beta_timing),
                      _["info_frac"]      = wrap(b.info_frac_used),
                      _["sd_incr"]        = wrap(b.sd_incr),
                      _["sd_proc"]        = wrap(b.sd_proc),
                      _["grid_collapses"] = b.grid_collapses,
                      _["grid_reversed"]  = b.grid_reversed,
                      _["slow_searches"]  = b.slow_searches,
                      _["unreachable_look"] = b.unreachable_look);
}

// ---- signature-compatible copies of RTSA's first.cpp used functions ---------
// [[Rcpp::export]]
NumericVector rtsa_init_int_cpp(NumericVector wj, NumericVector zj,
                                double delta, NumericVector stdv) {
  return wrap(rtsa::init_int(as_vec(wj), as_vec(zj), delta, stdv[0]));
}

// [[Rcpp::export]]
NumericVector rtsa_recur_int_cpp(int k, NumericMatrix stdv, NumericVector zj,
                                 NumericVector last, NumericVector zj_up,
                                 NumericVector wj_up, double delta, bool bs) {
  rtsa::vec sdi(stdv.nrow()), sdp(stdv.nrow());
  for (int i = 0; i < stdv.nrow(); ++i) { sdi[i] = stdv(i, 0); sdp[i] = stdv(i, 1); }
  return wrap(rtsa::recur_int(k, sdi, sdp, as_vec(zj), as_vec(last),
                              as_vec(zj_up), as_vec(wj_up), delta, bs));
}

// [[Rcpp::export]]
double rtsa_prob_cpp(double xq, NumericVector last, NumericVector zj, int k,
                     NumericMatrix stdv, bool bs, double delta) {
  rtsa::vec sdi(stdv.nrow()), sdp(stdv.nrow());
  for (int i = 0; i < stdv.nrow(); ++i) { sdi[i] = stdv(i, 0); sdp[i] = stdv(i, 1); }
  return rtsa::prob(xq, as_vec(last), as_vec(zj), k, sdi, sdp, bs, delta);
}
