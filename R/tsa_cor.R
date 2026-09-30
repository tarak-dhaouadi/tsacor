## Copyright (C) the RTSA authors (Anne Lyngholm Soerensen, Markus Harboe Olsen,
## Theis Lange, Christian Gluud) for the algorithms and code the boundary
## computations called from this file are derived from (RTSA 0.2.2, GPL
## (>= 2)); copyright (C) Tarak Dhaouadi for the adaptation. See DESCRIPTION
## and inst/COPYRIGHTS.
##
## tsa_cor(): Trial Sequential Analysis for meta-analyses of Pearson r / Spearman
## rho correlations pooled on Fisher's z scale. The structure, the random-
## effects / HKSJ machinery, the Diversity adjustment, the cumulative analysis,
## the RTSA-derived boundary orchestration and the decision logic are inherited
## from tsa_hr() of the tsahr package (same maintainer); the effect-measure
## specific parts (data handling, Fisher z, information size, participants
## instead of events, projection) are new in tsacor 0.1.0.

## Internal helper: map a metafor `method` code to a human-readable label
## for use in printed/plotted output (e.g. "DL" -> "DerSimonian-Laird").
## Falls back to the bare code itself (quoted) for any method not in the
## table, so an unrecognised-but-valid metafor method code never errors
## here -- it just prints less prettily. Not exported.
.tsacor_method_label <- function(method) {
  labels <- c(
    DL    = "DerSimonian-Laird",
    HE    = "Hedges",
    HS    = "Hunter-Schmidt",
    HSk   = "Hunter-Schmidt (with k correction)",
    SJ    = "Sidik-Jonkman",
    ML    = "maximum likelihood",
    REML  = "restricted maximum likelihood",
    EB    = "empirical Bayes",
    PM    = "Paule-Mandel",
    GENQ  = "generalized Q-statistic",
    PMM   = "Paule-Mandel (median-unbiased)",
    GENQM = "generalized Q-statistic (median-unbiased)",
    ## Aliases for the Hedges estimator; tsa_cor() normalises these to
    ## "HE" before they reach here, but keep them mapped so the helper is
    ## correct if called directly with an un-normalised code.
    CO    = "Hedges (Cochran alias)",
    VC    = "Hedges (variance-component alias)"
  )
  if (method %in% names(labels)) labels[[method]] else paste0("'", method, "'")
}

## Internal helpers for the random-effects inference option (0.2.8.11).
## Not exported.

## Human-readable label for a (canonical) re_inference value.
.tsacor_re_inference_label <- function(re_inference) {
  labels <- c(
    standard   = "standard (Wald-type z test)",
    hksj       = "Hartung-Knapp-Sidik-Jonkman (HKSJ)",
    hksj_adhoc = "HKSJ with ad hoc variance correction"
  )
  if (re_inference %in% names(labels)) labels[[re_inference]] else paste0("'", re_inference, "'")
}

## Validate/normalise `re_inference`. Case-insensitive; "knha" is an alias
## for "hksj" and "knha_adhoc" for "hksj_adhoc". Returns the canonical value
## (one of "standard", "hksj", "hksj_adhoc").
.tsacor_normalise_re_inference <- function(re_inference) {
  map <- c(standard   = "standard",
           hksj       = "hksj",
           knha       = "hksj",
           hksj_adhoc = "hksj_adhoc",
           knha_adhoc = "hksj_adhoc")
  msg <- paste0("re_inference must be one of: \"standard\", \"hksj\" (alias \"knha\"), ",
                "\"hksj_adhoc\" (alias \"knha_adhoc\"); matching is case-insensitive.")
  if (!is.character(re_inference) || length(re_inference) != 1L ||
      is.na(re_inference)) {
    stop(msg, call. = FALSE)
  }
  key <- tolower(trimws(re_inference))
  if (!(key %in% names(map))) stop(msg, call. = FALSE)
  unname(map[[key]])
}

## Hartung-Knapp-Sidik-Jonkman scale factor q for a random-effects fit with
## given tau^2: q = sum(w_i (y_i - mu)^2) / (k - 1), w_i = 1/(sei_i^2 + tau2)
## and mu the corresponding weighted mean. This is the multiplier that the
## HKSJ adjustment applies to the usual variance of the pooled estimate
## (identical to what metafor::rma(test = "knha") uses). NA for k < 2.
.tsacor_hksj_q <- function(yi, sei, tau2) {
  k <- length(yi)
  if (k < 2L || !is.finite(tau2)) return(NA_real_)
  w  <- 1 / (sei^2 + tau2)
  mu <- sum(w * yi) / sum(w)
  sum(w * (yi - mu)^2) / (k - 1)
}

## Re-do the inference of the cumulative (sequential) random-effects table
## under HKSJ (0.2.8.11). `cumul_df` is the data frame from
## metafor::cumul() of the STANDARD random-effects fit, so the cumulative
## point estimates and tau^2 values are untouched; only the standard error,
## test statistic, p-value and 95% CI are replaced, look by look, by their
## HKSJ counterparts (t distribution with k - 1 df, variance multiplied by q
## for "hksj" or by max(1, q) for "hksj_adhoc"). The Z column is the
## NORMAL-EQUIVALENT of the resulting t statistic -- sign(estimate) times the
## normal quantile with the same (two-sided) p-value -- so that the Z-curve
## stays on the scale the monitoring boundaries and the conventional
## boundary (qnorm(1 - alpha/2)) are defined on. Looks at which HKSJ is not
## defined (the first look, k = 1; or a non-positive/non-finite scale
## factor) keep the standard z-based values (re_df = Inf, re_scale = 1).
.tsacor_cumulative_re_inference <- function(cumul_df, yi, sei, re_inference, method) {
  k_tot   <- length(yi)
  z_eq    <- cumul_df$estimate / cumul_df$se
  scale_v <- rep(1, k_tot)
  df_v    <- rep(Inf, k_tot)
  se_v    <- cumul_df$se
  stat_v  <- z_eq
  p_v     <- cumul_df$pval
  lb_v    <- cumul_df$ci.lb
  ub_v    <- cumul_df$ci.ub

  for (i in seq_len(k_tot)) {
    if (i < 2L) next
    idx <- seq_len(i)
    tau2_i <- if ("tau2" %in% names(cumul_df)) cumul_df$tau2[i] else NA_real_
    if (!is.finite(tau2_i)) {
      tau2_i <- tryCatch(
        metafor::rma(yi = yi[idx], sei = sei[idx], method = method)$tau2,
        error = function(e) NA_real_)
    }
    q <- .tsacor_hksj_q(yi[idx], sei[idx], tau2_i)
    if (!is.finite(q) || q <= 0) next
    scale_i <- if (identical(re_inference, "hksj_adhoc")) max(1, q) else q
    w    <- 1 / (sei[idx]^2 + tau2_i)
    est  <- cumul_df$estimate[i]
    se_i <- sqrt(scale_i / sum(w))
    df_i <- i - 1
    t_i  <- est / se_i
    ## log-scale upper-tail probability keeps the normal-equivalent accurate
    ## even for very small p-values
    lp1  <- stats::pt(abs(t_i), df = df_i, lower.tail = FALSE, log.p = TRUE)
    crit <- stats::qt(0.975, df = df_i)

    scale_v[i] <- scale_i
    df_v[i]    <- df_i
    se_v[i]    <- se_i
    stat_v[i]  <- t_i
    p_v[i]     <- min(1, 2 * exp(lp1))
    lb_v[i]    <- est - crit * se_i
    ub_v[i]    <- est + crit * se_i
    z_eq[i]    <- sign(est) * (-stats::qnorm(lp1, log.p = TRUE))
  }

  ## 0.2.8.12: the HKSJ statistic is undefined at k = 1 (a t distribution
  ## with 0 df has no defined quantile), so the first look's Z is reported
  ## as NA rather than silently falling back to the z-based value -- both
  ## in the printed cumulative tables (which show this column directly)
  ## and in the plotted Z-curve (a single NA point/segment is simply not
  ## drawn; see plot.tsa_cor()). This is deliberate under every non-standard
  ## re_inference option, so it is set unconditionally here (this function
  ## is only ever called for "hksj"/"hksj_adhoc"; see tsa_cor()).
  z_eq[1] <- NA_real_

  cumul_df$se       <- se_v
  cumul_df$zval     <- stat_v
  cumul_df$pval     <- p_v
  cumul_df$ci.lb    <- lb_v
  cumul_df$ci.ub    <- ub_v
  cumul_df$re_scale <- scale_v
  cumul_df$re_df    <- df_v
  cumul_df$Z        <- z_eq
  cumul_df
}

## Internal helpers specific to correlations (tsacor 0.1.0). Not exported.

## Variance-inflation factor c of the Fisher-z transformed correlation,
## Var(z) = c / (n - 3):
##   Pearson r ............................................ c = 1
##   Spearman rho, Fieller-Hartley-Pearson (1957) ......... c = 1.06
##   Spearman rho, Bonett-Wright (2000) ................... c = 1 + rho^2 / 2
## `r_ref` is the (anticipated or observed) correlation the Bonett-Wright
## factor is evaluated at; it is ignored by the other two.
.tsacor_variance_factor <- function(cor_type, spearman_variance, r_ref = 0) {
  if (identical(cor_type, "pearson")) return(1)
  if (identical(spearman_variance, "fieller")) return(1.06)
  1 + r_ref^2 / 2
}

## Standard error of Fisher's z from the sample size (study by study).
.tsacor_se_from_n <- function(r, n, cor_type, spearman_variance) {
  if (identical(cor_type, "pearson")) return(1 / sqrt(n - 3))
  if (identical(spearman_variance, "fieller")) return(sqrt(1.06 / (n - 3)))
  sqrt((1 + r^2 / 2) / (n - 3))
}

## Standard error of Fisher's z implied by a reported confidence interval for
## the correlation: the CI is transformed with atanh() and its width divided by
## 2 * qnorm(1 - (1 - ci_level) / 2).
.tsacor_se_from_ci <- function(lbound, ubound, ci_level = 0.95) {
  (atanh(ubound) - atanh(lbound)) / (2 * stats::qnorm(1 - (1 - ci_level) / 2))
}

## Human-readable label for the correlation type.
.tsacor_cor_label <- function(cor_type) {
  if (identical(cor_type, "spearman")) "Spearman rho" else "Pearson r"
}

#' Trial Sequential Analysis for a meta-analysis of correlations
#'
#' Performs a Trial Sequential Analysis (TSA) for a meta-analysis of Pearson
#' product-moment correlations (\eqn{r}) or Spearman rank correlations
#' (\eqn{\rho}) pooled on Fisher's \eqn{z} scale. Adapts the classical
#' Wetterslev/Thorlund/CTU TSA framework to correlations: study-level
#' correlations are transformed with \eqn{z = \mathrm{atanh}(r)}, pooled with a
#' random-effects model, and monitored against O'Brien-Fleming-type alpha- and
#' beta-spending boundaries computed on the inverse-variance information scale
#' with the recursive numerical integration engine ported from the R package
#' 'RTSA' (Soerensen, Olsen, Lange and Gluud). The required information size is
#' corrected for heterogeneity with the Diversity (D-squared) adjustment of
#' Wetterslev et al. (2009).
#'
#' @param data A data.frame, or a path to an .xlsx file (read with
#'   \code{readxl::read_excel()}), containing one row per study with (at least)
#'   the columns \code{Study}, \code{r} (the Pearson correlation or Spearman
#'   rho) and \code{n_subjects} (the number of subjects the correlation is based
#'   on). When \code{se_source = "ci"} (the default) the columns \code{lbound}
#'   and \code{ubound} (the lower and upper limits of the confidence interval of
#'   the correlation, on the \emph{correlation} scale, not on the Fisher
#'   \eqn{z} scale) are also required. Any other columns (e.g. \code{Year},
#'   \code{Ethnicity}, \code{Age}) are kept in the returned \code{data} and can
#'   be used with \code{order_by}. Each row is treated as one independent study,
#'   and all studies must report the same kind of coefficient (see
#'   \code{cor_type}). Rows are treated as being in chronological (publication)
#'   order; reorder your data, or use \code{order_by}, if the row order in your
#'   file is not chronological.
#'
#'   The data are validated before any model is fitted. \code{Study} must be
#'   non-missing, non-empty and unique; at least two studies are required (an
#'   error), and fewer than 10 studies raise a warning because the heterogeneity
#'   and Diversity estimates, and therefore DARIS and the monitoring
#'   boundaries, can be very unstable with so few studies. \code{r},
#'   \code{n_subjects} (and \code{lbound}/\code{ubound} when they are used) must
#'   be numeric and finite; \code{r} must lie strictly between -1 and 1 (Fisher's
#'   \eqn{z} is infinite at \eqn{|r| = 1}); \code{n_subjects} must be whole
#'   numbers greater than 3 (the variance of Fisher's \eqn{z} is
#'   \eqn{1/(n - 3)}); and, when used, \code{lbound} and \code{ubound} must lie
#'   strictly between -1 and 1 with \code{lbound < ubound}. Column names in
#'   \code{data} have spaces replaced with underscores on load (so an Excel
#'   header \dQuote{n subjects} becomes \code{n_subjects}). If two distinct
#'   headers would collide once spaces become underscores, \code{tsa_cor()}
#'   stops rather than silently using whichever column came first.
#' @param alpha_two_sided Overall two-sided type I error for the TSA
#'   monitoring boundaries; it also fixes the critical value
#'   \eqn{z_{1-\alpha/2}} in the required information size. A single finite
#'   value strictly between 0 and 1. Default \code{0.05}.
#' @param power Desired power (\eqn{1 - \beta}) for the required information
#'   size calculation; \eqn{\beta = 1 - } \code{power} is the error that the
#'   beta-spending (futility) boundaries spend. A single finite value strictly
#'   between 0 and 1. Default \code{0.80}.
#' @param target_r Anticipated (target) correlation used for the required
#'   information size calculation, on the correlation scale (not on the
#'   Fisher \eqn{z} scale). Must be strictly between -1 and 1 and non-zero
#'   (\eqn{\mathrm{atanh}(0) = 0} makes the required information infinite, so
#'   \code{target_r = 0} is an error); only its absolute value matters, because
#'   the boundaries are two-sided and the required information depends on
#'   \eqn{\mathrm{atanh}(r_0)^2}. A value with \eqn{|r_0| < 0.10} is accepted
#'   but raises a warning, since the required information size grows as
#'   \eqn{1 / \mathrm{atanh}(r_0)^2} and therefore increases rapidly as the
#'   target approaches 0. Default \code{NA}, which uses the \emph{observed}
#'   pooled correlation from the random-effects meta-analysis -- see Details
#'   for an important caution about this default (a warning is always raised,
#'   and an observed pooled correlation of exactly 0 is an error). Set to a
#'   pre-specified, clinically or scientifically meaningful value (e.g.
#'   \code{0.20}) for a standard, non-circular, protocol-driven TSA. Note that
#'   the target drives the whole design: a target that is small relative to the
#'   observed pooled correlation demands far more information than the studies
#'   supply (DARIS is then not reached and the retrospective projection is
#'   reported), whereas a target that is large relative to it makes DARIS small
#'   and easily reached.
#' @param cor_type Character string, \code{"pearson"} (default) or
#'   \code{"spearman"}: the kind of correlation coefficient in the
#'   \code{r} column. The shorthands \code{"r"} (for \code{"pearson"}) and
#'   \code{"rho"} (for \code{"spearman"}) are also accepted, and matching is
#'   case- and whitespace-insensitive (e.g. \code{"R"}, \code{" Rho "},
#'   \code{"PEARSON"} all work). Any other value (e.g. Kendall's tau, which is
#'   not supported) is an error. It determines the variance model of
#'   Fisher's \eqn{z} (see \code{se_source} and \code{spearman_variance}), the
#'   variance factor \eqn{c} used to translate information into participants
#'   (\eqn{c = 1} for Pearson), and the wording of the printed and plotted
#'   output. All studies must report the same kind of coefficient.
#' @param se_source Character string, \code{"ci"} (default) or \code{"n"},
#'   giving where the standard error of each study's Fisher \eqn{z} comes
#'   from. \code{"ci"}: from the reported confidence interval,
#'   \eqn{SE(z) = (\mathrm{atanh}(ubound) - \mathrm{atanh}(lbound)) /
#'   (2 q)} with \eqn{q} the normal quantile of \code{ci_level}; this
#'   uses whatever variance model produced the published interval and
#'   therefore also works for Spearman coefficients whose intervals were
#'   computed by other methods. \code{"n"}: from the sample size only,
#'   \eqn{SE(z) = \sqrt{c / (n - 3)}} with \eqn{c = 1} for Pearson r and
#'   \eqn{c} chosen by \code{spearman_variance} for Spearman rho. Both
#'   versions are always returned in \code{data} (\code{se_z_ci} and
#'   \code{se_z_n}) when the data allow it (\code{se_z_ci} is \code{NA} without
#'   usable interval columns), so the choice can be checked as a sensitivity
#'   analysis; the column \code{se_z} holds the version actually used. With
#'   \code{"ci"}, two consistency checks raise warnings: a correlation lying
#'   outside its own interval, and an interval whose midpoint on the Fisher
#'   \eqn{z} scale differs from \eqn{\mathrm{atanh}(r)} by more than 0.25
#'   standard errors (it then does not look like a Fisher-\eqn{z} interval, for
#'   example a bootstrap interval, and the width-based standard error is only
#'   approximate; \code{se_source = "n"} is then worth considering). When
#'   \code{verbose = TRUE} the median and range of the ratio of the CI-based to
#'   the n-based standard error are printed.
#' @param spearman_variance Character string, \code{"fieller"} (default) or
#'   \code{"bonett_wright"}. Only used when \code{cor_type = "spearman"}: the
#'   variance model of the Fisher-transformed Spearman coefficient,
#'   \eqn{Var(z) = 1.06 / (n - 3)} (Fieller, Hartley and Pearson, 1957) or
#'   \eqn{Var(z) = (1 + \rho^2 / 2) / (n - 3)} (Bonett and Wright, 2000). It
#'   is used (i) to compute the standard errors when \code{se_source = "n"}
#'   (the Bonett-Wright version uses each study's own \eqn{\rho}) and
#'   (ii) always, to translate information units into an equivalent number of
#'   participants (see Details), including when the standard errors come from
#'   the confidence intervals; for that translation the Bonett-Wright factor is
#'   evaluated at the anticipated correlation. Ignored for Pearson correlations
#'   (\eqn{c = 1}).
#' @param ci_level Confidence level of the reported intervals in the
#'   \code{lbound}/\code{ubound} columns, used only when
#'   \code{se_source = "ci"} (and for the consistency check between the two
#'   standard-error versions). A single finite value strictly between 0 and 1.
#'   Default \code{0.95}. It must match the level at which the intervals in
#'   \code{data} were actually computed: the standard errors are scaled by the
#'   corresponding normal quantile, so a wrong level rescales every standard
#'   error by the same factor.
#'
#' @param method Character string specifying the heterogeneity-variance
#'   (tau^2) estimator used for the \emph{random-effects} meta-analysis and
#'   cumulative (sequential) TSA model, passed to \code{method} in
#'   \code{metafor::rma()} after validation and, for the aliases
#'   \code{"CO"}/\code{"VC"}, normalisation to \code{"HE"}. One of \code{"DL"}
#'   (DerSimonian-Laird, the default -- kept as the default here for consistency
#'   with the sister package tsahr, whose analyses always used DL; note this
#'   differs from \code{metafor::rma()}'s own default of \code{"REML"}),
#'   \code{"HE"}, \code{"HS"}, \code{"HSk"}, \code{"SJ"}, \code{"ML"},
#'   \code{"REML"}, \code{"EB"}, \code{"PM"}, or \code{"PMM"}. See
#'   \code{?metafor::rma} for the definition of each estimator. The string is
#'   matched exactly (case-sensitively, e.g. \code{"REML"}, not \code{"reml"});
#'   anything else is an error that lists the supported values.
#'
#'   \emph{What \code{method} changes.} The estimator determines tau^2 (and
#'   I^2), hence the random-effects weights, the pooled Fisher \eqn{z}
#'   and its standard error (\code{res_re}), and every row of the
#'   cumulative analysis (estimate, tau^2, standard error, \code{Z}). Because
#'   the Diversity \eqn{D^2 = (Var_{RE} - Var_{FE}) / Var_{RE}} is computed from
#'   the variance of the random-effects pooled estimate obtained with the
#'   chosen estimator (always from the standard, z-based fit, whatever
#'   \code{re_inference} says), \code{method} also changes \eqn{D^2}, the
#'   adjustment factor \eqn{1 / (1 - D^2)}, DARIS, the information fractions
#'   and, through them, the position of the alpha-/beta-spending boundaries
#'   (which are computed on the information-fraction timeline), the
#'   \dQuote{DARIS reached} verdict and the retrospective projection. The
#'   boundary engine itself does not depend on \code{method}; only the
#'   information fractions it is run on do.
#'
#'   \emph{What \code{method} does not change.} Cochran's Q, the study-level
#'   standard errors, the accrued information (the sum of
#'   \code{1 / se_z^2}, a fixed-effect quantity), the unadjusted required
#'   information size and its participant equivalent (unless
#'   \code{target_r = NA}, in which case the anticipated correlation is the
#'   pooled estimate, which does depend on \code{method}), and the equal-effects
#'   model used as the comparator for \eqn{D^2}, which is always fitted with
#'   \code{method = "FE"}. \code{"FE"} is deliberately not offered for the
#'   random-effects model: the two variances would coincide and \eqn{D^2} would
#'   be 0 by construction.
#'
#'   \emph{Choosing an estimator.} \code{"DL"} is a simple, non-iterative
#'   moment estimator that can underestimate tau^2, notably with few studies
#'   or substantial heterogeneity; iterative estimators such as \code{"REML"}
#'   or \code{"PM"} are frequently preferred in the methodological literature.
#'   Because \eqn{D^2} and DARIS depend on the estimator, repeating the
#'   analysis with more than one \code{method} is a worthwhile sensitivity
#'   analysis. The iterative estimators (\code{"ML"}, \code{"REML"},
#'   \code{"EB"}, \code{"PM"}, \code{"PMM"}) can fail to converge, especially
#'   with few studies; because the cumulative analysis re-fits the model at
#'   every look, the estimator must also be usable on the first few studies,
#'   and an error from \code{metafor} is passed on unchanged.
#'
#'   \emph{Aliases.} \code{"CO"} and \code{"VC"} are accepted as aliases for
#'   \code{"HE"} and are normalised to \code{"HE"} before being passed to
#'   \code{metafor::rma()}. The \code{metafor} documentation notes that the
#'   Hedges estimator is also known as the Cochran (\code{"CO"}) or
#'   variance-component (\code{"VC"}) estimator, and that those strings may be
#'   used to select it -- but that alias is not accepted by every
#'   \code{metafor} version (older releases reject a bare \code{"CO"} with
#'   \dQuote{Unknown 'method' specified}). Normalising here makes
#'   \code{tsa_cor()} behave identically across \code{metafor} versions rather
#'   than inheriting that version skew, and avoids pinning a minimum
#'   \code{metafor} version purely for an alias. All three strings denote the
#'   same estimator, so this has no numerical consequence. The returned object
#'   records both \code{parameters$method} (the normalised string actually
#'   used, i.e. \code{"HE"}) and \code{parameters$method_requested} (what the
#'   caller passed).
#'
#'   \emph{Unsupported strings.} Two method strings accepted by
#'   \code{metafor::rma()} are deliberately NOT supported here:
#'   \code{"GENQ"} and \code{"GENQM"} require the caller to also supply a
#'   \code{weights} argument to \code{metafor::rma()}, which \code{tsa_cor()}
#'   does not currently collect or pass through, so passing them here raises an
#'   explicit error explaining why rather than silently forwarding to
#'   \code{metafor::rma()} and surfacing its own unrelated error.
#' @param order_by Optional name of a column in \code{data} to sort by
#'   (ascending) before the cumulative analysis, e.g. \code{"Year"}. TSA is
#'   order-dependent, so getting the chronological order right matters.
#'   Default \code{NULL}, which uses the row order already present in
#'   \code{data} and assumes it is chronological (with no way for the package
#'   to verify this). Column names in \code{data} have spaces replaced with
#'   underscores on load; \code{order_by} is normalised the same way, so either
#'   \code{"Publication Year"} or \code{"Publication_Year"} will match that
#'   column. If two distinct headers would collide once spaces become
#'   underscores, \code{tsa_cor()} stops rather than silently using whichever
#'   column came first. A column that does not exist is an error. Warnings are
#'   raised when the column contains \code{NA} (those rows sort to one end under
#'   R's default handling), is neither numeric nor a Date/POSIXct (R's default
#'   ordering for its type, e.g. lexical for character, may then not be
#'   chronological), or contains tied values (tied studies keep their original
#'   relative order, which may affect the cumulative analysis; a finer-grained
#'   column such as a publication date can be used instead of the year alone).
#' @param verbose Logical; print analysis details to the console as the
#'   function runs (data loading, the random-effects model, heterogeneity,
#'   required information size, the cumulative analysis and the boundaries).
#'   Default \code{TRUE}. Warnings are raised regardless of \code{verbose}.
#' @param boundary_route Character string, one of \code{"design"} (default)
#'   or \code{"analysis"}, selecting which of RTSA's two retrospective
#'   boundary-computation routes to use. \code{"design"} computes the alpha and
#'   beta boundaries on the observed information-fraction timeline with no
#'   further inflation, and the formal endpoint is DARIS itself.
#'   \code{"analysis"} first solves an inflation factor \code{design_R} and
#'   recomputes the boundaries on the timeline scaled by it, so the formal
#'   endpoint becomes \code{design_R * DARIS}, and the efficacy boundaries change
#'   as well as the futility ones. This changes more than the futility numbers:
#'   the formal endpoint, the route-endpoint verdict (\code{final_reached}) and
#'   every decision field in \code{results} refer to the route endpoint. See
#'   \dQuote{Route endpoint versus DARIS} and \dQuote{Retrospective boundary
#'   timeline} under Details.
#' @param legacy_fallback Logical, default \code{TRUE}. Governs what happens
#'   if the compiled RTSA-derived boundary engine fails to produce a result
#'   for the requested design (e.g. no root bracket exists for an unusual
#'   information-fraction schedule). When \code{TRUE} (the default),
#'   \code{tsa_cor()} falls back to a legacy, R-only approximate engine, with
#'   an immediate warning and a visible banner in \code{print()}/\code{summary()}
#'   output, and marks the result (\code{beta_engine$engine ==
#'   "legacy_r_fallback"}) so the fallback is never silent -- but it IS a
#'   fallback: a caller who wraps the call in \code{suppressWarnings()} will not
#'   see the warning, and the returned boundaries are not RTSA-comparable when
#'   this happens. Set \code{legacy_fallback = FALSE} for strict fail-closed
#'   behaviour: an engine failure then stops \code{tsa_cor()} with an error
#'   instead of silently substituting the approximate engine, appropriate when
#'   the result will be reported as RTSA-equivalent and an unnoticed fallback
#'   would be worse than a hard stop. Must be a single \code{TRUE} or
#'   \code{FALSE}.
#' @param projection_stat Character string, one of \code{"median"} (default)
#'   or \code{"mean"}, selecting the summary statistic used to turn the
#'   OBSERVED per-study contributions into \dQuote{typical} values for the
#'   retrospective projection described under \dQuote{Estimated additional
#'   studies/participants}: the information per study (\code{1/se_z^2}, giving
#'   the additional-studies estimate) and the information per participant
#'   (\code{1/se_z^2} divided by the study's \code{n_subjects}, giving the
#'   additional-participants estimate). \code{"median"} is the default because
#'   it is more robust to a single unusually large or small study;
#'   \code{"mean"} is offered as an alternative when that robustness is not
#'   wanted (e.g. a deliberately evenly-sized set of studies). This has no
#'   effect on any other quantity returned by \code{tsa_cor()} -- it only feeds
#'   the \code{projection} element of the return value and the corresponding
#'   printed/plotted text.
#' @param info_per_participant_basis Character string, one of
#'   \code{"per_study"} (default) or \code{"pooled"}, selecting how the
#'   \dQuote{historical rate} of information per participant used for the
#'   additional-PARTICIPANTS projection is obtained. \code{"per_study"}: the
#'   \code{projection_stat} (median by default, or mean) of each study's own
#'   information per participant (\code{1/se_z^2} divided by its
#'   \code{n_subjects}). \code{"pooled"}: the ratio of sums (total information /
#'   total participants), which weights larger studies more and equals the slope
#'   of observed cumulative information against cumulative participants. It does
#'   not affect the additional-STUDIES estimate, which always uses
#'   \code{projection_stat}. Both variants are always returned in
#'   \code{projection} as a sensitivity check. (Unlike a time-to-event analysis
#'   with zero-event studies, every study has \code{n_subjects > 3} and a finite
#'   positive standard error, so no study is ever excluded from these rates.)
#' @param re_inference Character string selecting how inference on the pooled
#'   \emph{random-effects} effect, and on every cumulative look, is carried
#'   out; matching is case-insensitive. The heterogeneity-variance estimator
#'   (\code{method}) and the pooled point estimate are the same under every
#'   option -- only the standard error, test statistic, p-value and 95\% CI
#'   change. One of:
#'   \describe{
#'     \item{\code{"standard"}}{(default) the usual Wald-type inference:
#'       standard normal reference distribution, variance \code{1 / sum(w)}
#'       with random-effects weights \code{w = 1 / (se_z^2 + tau^2)}.}
#'     \item{\code{"hksj"} (alias \code{"knha"})}{Hartung-Knapp-Sidik-Jonkman
#'       (Knapp-Hartung) adjustment: the variance is multiplied by
#'       \code{q = sum(w * (z - pooled)^2) / (k - 1)} and the test and CI use
#'       a t distribution with \code{k - 1} degrees of freedom
#'       (\code{metafor::rma(test = "knha")}). \code{q} may be below 1, in
#'       which case the adjusted CI can be \emph{narrower} than the standard
#'       one.}
#'     \item{\code{"hksj_adhoc"} (alias \code{"knha_adhoc"})}{HKSJ with the
#'       ad hoc correction that \code{q} is never allowed to be below 1: the
#'       variance is multiplied by \code{max(1, q)} (t distribution,
#'       \code{k - 1} df), so the standard error is never smaller than under
#'       \code{"standard"}. Also known as the modified/truncated HKSJ method.}
#'   }
#'   The cumulative Z-curve then uses the normal-equivalent of the HKSJ t
#'   statistic (the value with the same two-sided p-value on the normal scale),
#'   so that it stays on the scale of the monitoring boundaries; it is
#'   \code{NA} at the first look, where HKSJ is undefined. The Diversity
#'   D^2, the adjustment factor, DARIS, the information fractions and hence
#'   the boundaries do not depend on \code{re_inference} (they do depend on
#'   \code{method}). A warning about the instability of early HKSJ looks is
#'   raised whenever a non-standard option is actually used, and the plot
#'   caption, \code{print()} and \code{summary()} name the option. See
#'   \dQuote{Random-effects inference} under Details.
#'
#' @details
#' \strong{Fisher z scale and information.} Each study contributes
#' \eqn{z_i = \mathrm{atanh}(r_i)} with standard error \eqn{se_i} (see
#' \code{se_source}); the pooled effect, its standard error, the cumulative
#' Z-curve, tau^2, I^2, Q and the Diversity D^2 all live on this scale, and
#' the pooled Fisher z and its confidence limits are back-transformed with
#' \code{tanh()} for reporting. Statistical information is the inverse
#' variance, \code{1 / se_z^2}, which for Pearson correlations is close to
#' \eqn{n - 3}. With the anticipated effect \eqn{z_0 = \mathrm{atanh}(r_0)}
#' the required (allocation-free) information is
#' \deqn{I = (z_{1-\alpha/2} + z_{power})^2 / z_0^2,}
#' the analogue of the Schoenfeld formula used for hazard ratios, and the
#' diversity-adjusted required information size is
#' \eqn{DARIS = I / (1 - D^2)}.
#'
#' \strong{Participants instead of events.} tsahr, the sister package for
#' hazard ratios, draws the Z-curve against cumulative \emph{events}; the
#' natural counterpart for correlations is the cumulative number of
#' \emph{participants} (\code{n_subjects}). The theoretical sample-size
#' equivalent of an information target \eqn{I} is
#' \eqn{c I + 3}, i.e. the familiar single-study sample-size formula
#' \eqn{n = c (z_{1-\alpha/2} + z_{power})^2 / z_0^2 + 3}, with the
#' variance factor \eqn{c} of the correlation type (1 for Pearson; 1.06 or
#' \eqn{1 + \rho^2/2} for Spearman, see \code{spearman_variance}). It is
#' reported as \code{RIS_participants} (from \eqn{I}) and
#' \code{DARIS_participants} (from DARIS). Because a meta-analysis of \eqn{k}
#' studies loses 3 units of information per study, whereas this equivalent
#' loses them only once, it is a single-study-equivalent reference, not an
#' exact count of the participants needed across studies; the
#' observed-information criterion (\code{info_accrued} against DARIS) is the
#' one used for all "reached" verdicts.
#'
#' \strong{Circularity caution:} using the observed pooled effect
#' (\code{target_r = NA}) to determine the required information size is
#' circular -- it tends to make the required information size small whenever
#' the pooled correlation is large and precise, which can make the TSA
#' boundary collapse to the conventional boundary almost immediately. For a
#' publication-quality TSA, set \code{target_r} to a value fixed independently
#' of (and ideally before looking at) the meta-analysis result. Because the
#' required information grows as \eqn{1 / \mathrm{atanh}(r_0)^2}, target
#' correlations close to 0 require very large amounts of information.
#'
#' \strong{Random-effects caveat:} the cumulative Z-curve is estimated from a
#' random-effects model whose between-study variance is re-estimated at every
#' step. The canonical Lan-DeMets/O'Brien-Fleming theory assumes a fixed,
#' canonical information process with independent Brownian-motion increments;
#' the displayed monitoring boundaries should therefore be regarded as an
#' approximation (as in the official Copenhagen Trial Unit TSA software), not
#' an exact result.
#'
#' \strong{Random-effects inference (\code{re_inference}).} \code{re_inference}
#' changes how the random-effects pooled effect and each cumulative look are
#' tested, not how heterogeneity is estimated: \code{method} still chooses the
#' tau^2 estimator, and the pooled Fisher z and tau^2 at every look are
#' identical under all three options. For \code{"hksj"} and \code{"hksj_adhoc"}
#' the standard error, p-value and 95\% CI of the pooled correlation
#' (\code{res_re}, printed output, summary table and plot subtitle) and of every
#' row of \code{cumulative} are the HKSJ ones, computed look by look from the
#' studies accrued so far.
#'
#' \emph{Z-curve scale.} The monitoring boundaries and the conventional
#' boundary are defined on the standard normal scale, whereas an HKSJ statistic
#' follows a t distribution with \code{k - 1} df at look \code{k}. The
#' cumulative Z is therefore the normal-equivalent of the HKSJ t statistic: the
#' value with the same two-sided p-value on the normal scale
#' (\code{sign(estimate) * qnorm(1 - p / 2)}). Consequently
#' \code{|Z| >= qnorm(1 - alpha/2)} exactly when the HKSJ p-value is below
#' alpha, and the conventional (naive) boundary keeps its usual meaning. The raw
#' t statistic is kept in \code{cumulative$zval}, the degrees of freedom in
#' \code{cumulative$re_df} and the variance multiplier in
#' \code{cumulative$re_scale} (these columns exist only for non-standard
#' options). The t distribution is not what the alpha-spending boundaries were
#' derived for, so this is a further approximation on top of the one described
#' under \dQuote{Random-effects caveat}.
#'
#' \emph{Early looks.} HKSJ needs at least two studies: the HKSJ statistic is
#' undefined at \code{k = 1}, so \code{cumulative$Z} is \code{NA} at the first
#' look for \code{"hksj"}/\code{"hksj_adhoc"} -- shown as \code{NA} in the
#' printed cumulative tables and simply not drawn on the plotted Z-curve
#' (\code{cumulative$se}, \code{$pval} etc. still hold the standard, z-based
#' values at that look; only Z is withheld). A look at which the scale factor is
#' zero or not finite (e.g. all estimates identical) likewise keeps the standard
#' z-based Z; both cases have \code{re_df = Inf} and \code{re_scale = 1}. At the
#' second look there is \code{k - 1 = 1} degree of freedom, so early HKSJ looks
#' are very conservative and erratic; a \code{warning()} says so (with
#' \code{call. = FALSE}) whenever a non-standard \code{re_inference} is actually
#' used, regardless of \code{verbose}. Early cumulative HKSJ values should not
#' be interpreted as directly comparable in magnitude with conventional normal
#' Z-statistics.
#'
#' \emph{What is not changed.} The Diversity D^2, the adjustment factor,
#' DARIS, the information fractions and hence the alpha/beta boundaries are
#' always computed from the standard random-effects and equal-effects variances
#' and do not depend on \code{re_inference}, as do I^2, tau^2 and Q (they do
#' depend on \code{method}). Only the Z-curve (and the decisions that compare
#' it with the boundaries) and the reported pooled inference change. With
#' \code{target_r = NA} the anticipated correlation is the pooled point
#' estimate, which is also unaffected.
#'
#' \emph{Degenerate data.} If the HKSJ scale factor of the full data set is zero
#' or not finite, \code{tsa_cor()} warns and falls back to \code{"standard"};
#' \code{parameters$re_inference} then reads \code{"standard"} while
#' \code{parameters$re_inference_requested} keeps what was asked for.
#'
#' \strong{Retrospective boundary timeline:} the observed cumulative Z-curve
#' continues through every included study, but the formal alpha and beta
#' boundaries follow the RTSA retrospective convention: observed looks are
#' retained only while \code{info_fraction < 1}, followed by one synthetic
#' final-analysis point at \code{t = 1} (HARIS/DARIS). Boundaries are not
#' continued through studies occurring after DARIS. The synthetic point is
#' stored in \code{boundary_timeline}; the observed-study rows in
#' \code{cumulative} have boundary values set to \code{NA} at and after DARIS.
#' Formal crossing/futility decisions are evaluated only through the first
#' observed look reaching the route endpoint (DARIS for
#' \code{boundary_route = "design"}; \code{design_R * DARIS} for
#' \code{"analysis"}), using the definitive boundary at that endpoint.
#'
#' \strong{Decision fields: "at any formal look" versus "at the definitive
#' look".} \code{crossed_tsa} and \code{entered_futility_region} are
#' \code{TRUE} if the cumulative Z-curve crossed the efficacy boundary
#' (respectively lay inside the futility region) at any look up to and
#' including the definitive one. A study series that crossed efficacy at an
#' interim look keeps \code{crossed_tsa = TRUE} even if the definitive look then
#' fell back below the boundary. \code{final_crossed_efficacy},
#' \code{final_non_efficacy} (= \code{!final_crossed_efficacy}) and
#' \code{final_entered_futility_region} refer ONLY to the first look reaching
#' the route endpoint (\code{results$final_tsa_look}) and are \code{NA} when
#' that endpoint has not been reached. At the definitive look the futility
#' boundary equals the final efficacy boundary (as in RTSA's design pass, where
#' the two are calibrated to meet there; e.g. about 2.1-2.2 rather than 1.959964
#' for a two-sided alpha of 0.05 with many looks), so
#' \code{final_entered_futility_region} is the complement of
#' \code{final_crossed_efficacy} (both are TRUE only at \code{|Z|} exactly equal
#' to the boundary); it merely says that the Z-curve did not reach the final
#' efficacy boundary and is not a formal interim futility stop. Neither family
#' of fields is a recommendation to stop. \code{crossed_conventional} is
#' deliberately evaluated over the FULL cumulative Z-curve, including studies
#' added after DARIS, unlike \code{crossed_tsa}, which is restricted to the
#' formal decision horizon: it answers \dQuote{did the naive, uncorrected
#' cumulative P-value ever drop below alpha?}, as a deliberate contrast that
#' illustrates the repeated-testing inflation TSA guards against.
#'
#' \strong{Route endpoint versus DARIS.} \code{boundary_route = "design"}
#' (default) is a port of \code{RTSA::boundaries(type = "design")}: alpha and
#' beta boundaries are computed on the observed information-fraction timeline
#' (looks below \code{t = 1} plus one synthetic \code{t = 1} point) with no
#' further inflation, and the formal endpoint is DARIS itself.
#' \code{boundary_route = "analysis"} is a port of \code{RTSA::RTSA(type =
#' "analysis", design = NULL)}: a design pass first solves an inflation factor
#' (\code{design_R}); the boundaries returned are then recomputed on the
#' timeline scaled by \code{design_R}, and the formal endpoint becomes
#' \code{design_R * DARIS}. DARIS and the route endpoint are reported
#' separately: \code{results$daris_reached} and
#' \code{information_size$DARIS_info_threshold_n} refer to DARIS itself,
#' whereas \code{results$final_reached}, \code{information_size$route_endpoint_info}
#' and \code{information_size$route_endpoint_n} refer to the route endpoint
#' (identical to DARIS for \code{"design"}).
#'
#' \strong{Fallbacks.} \code{settings$route_used} (\code{"design"},
#' \code{"analysis"} or \code{"legacy"}), \code{settings$fallback_used},
#' \code{settings$fallback_route} and \code{settings$fallback_reason} record
#' what actually produced the boundaries: \code{"design"} means
#' \code{boundary_route = "analysis"} failed and the design-route result is
#' returned; \code{"legacy"} means the compiled engine failed and the legacy,
#' approximate R-only engine was used. Both are announced by warnings. For
#' confirmatory or RTSA-parity work use \code{legacy_fallback = FALSE}.
#'
#' \strong{Numerical diagnostics and RTSA quirks kept.} The compiled engine
#' warns when a boundary search converged only within a loose tolerance, when
#' an integration grid collapsed to a degenerate interval, and -- with a
#' separate, more alarming warning -- when an integration interval was
#' REVERSED (lower wall above the upper wall), a state RTSA's own code would
#' have stopped on and which invalidates the boundaries from that look onward.
#' Diagnostics are reported for the converged passes, not for the transient
#' candidate information scales tried inside the root searches.
#'
#' \emph{Placeholder boundaries: RTSA's tolerance, kept.} At looks where the
#' cumulative alpha spend is below RTSA's absolute search tolerance of 1e-9
#' (very early looks, small information fractions) RTSA does not solve for a
#' boundary but reports the placeholder value 20, and tsacor reproduces this
#' deliberately (\code{TSA_boundary_upper == 20}). It is not a computed
#' boundary: the true boundary there is finite, and a Z-curve above the
#' placeholder is not counted as crossing. The same tolerance also makes the
#' first boundary that IS solved after such looks slightly inaccurate, with an
#' error shrinking at later looks. Both effects are RTSA's own; tsacor keeps
#' RTSA's tolerance so that its bounds match RTSA's, and the placeholder
#' stretches the vertical axis of \code{plot()}. A Z-curve that starts with a
#' very small information fraction (a target correlation that is small
#' compared with the precision of the first studies) is the typical trigger.
#'
#' \emph{Final-look-only beta spend.} When every interim look is suppressed
#' (small information fractions) the final look carries the whole beta spend.
#' RTSA has an exact-float shortcut for that case that fires or not depending on
#' the last bit of beta; tsacor computes \code{beta = 1 - power}, which at power
#' 0.80 lands on the wrong side of that bit and would lose the compiled engine
#' to the legacy fallback. The shortcut is deliberately not ported, so the
#' result no longer depends on the last bit of beta.
#'
#' \strong{Estimated additional studies/participants (retrospective
#' projection).} Whenever the route's target has NOT been reached in the
#' observed data (\code{results$daris_reached == FALSE} for
#' \code{boundary_route = "design"}; \code{results$final_reached == FALSE} for
#' \code{"analysis"}), \code{tsa_cor()} additionally projects how many more
#' studies, and roughly how many more participants, would be needed to reach it,
#' based directly on the OBSERVED study-level information increments already in
#' the data, not on any new assumption about the size of future studies. The
#' target information \eqn{I_{required}} differs by route and is never
#' conflated: DARIS (\code{information_size$DARIS_info}) for \code{"design"};
#' \code{design_R * DARIS} (\code{information_size$route_endpoint_info}) for
#' \code{"analysis"}. The shortfall \eqn{I_{required}} minus the accrued
#' information is translated in two separate ways, which answer different
#' questions:
#' \describe{
#'   \item{Additional participants (primary)}{The shortfall divided by a
#'     \dQuote{historical rate} of information per participant, chosen with
#'     \code{info_per_participant_basis}: \code{"per_study"} uses the
#'     \code{projection_stat} of each study's own \code{1/se_z^2} divided by its
#'     \code{n_subjects}; \code{"pooled"} uses total information divided by total
#'     participants. The result is continuous and is rounded up only for
#'     display. It is returned in
#'     \code{projection$additional_participants_estimated}; both variants are
#'     always returned (\code{additional_participants_study_level} and
#'     \code{additional_participants_pooled}) so the choice can be checked as a
#'     sensitivity analysis.}
#'   \item{Additional studies (secondary)}{The shortfall divided by a single
#'     \dQuote{typical future study} information increment (the
#'     \code{projection_stat} of each study's own \code{1/se_z^2}), rounded up to
#'     the nearest whole study, so it is always a natural number of at least 1
#'     whenever there is a genuine shortfall. It is not affected by
#'     \code{info_per_participant_basis}.}
#' }
#' The two are deliberately not chained: the participants figure is NOT the
#' number of studies times a participants-per-study increment, because rounding
#' the studies up would inflate it. That chained quantity (whole typical studies
#' times the typical study size) is still returned, as
#' \code{projection$participants_from_whole_studies}, but it answers \dQuote{how
#' many participants would the minimum whole number of typical studies bring?},
#' not \dQuote{how many participants are needed?}.
#'
#' When the route's own target has not been reached, three separate quantities
#' are reported (printed output, \code{summary()} and the plot caption; for
#' \code{boundary_route = "design"} the target is DARIS, for \code{"analysis"}
#' the analysis-route endpoint):
#' \describe{
#'   \item{Theoretical additional participants}{The direct, deterministic
#'     difference between the single-study-equivalent participant target of the
#'     route endpoint (\eqn{c \cdot} route endpoint \eqn{\cdot} DARIS + 3) and the
#'     participants accrued, floored at 0
#'     (\code{projection$additional_participants_theoretical}).}
#'   \item{Estimated additional participants (historical rate)}{The participants
#'     projection above, with \eqn{I_{required}} set to the route's target
#'     information (\code{projection$additional_participants_estimated}).
#'     \code{plot()} draws it as a second reference line, at accrued participants
#'     plus these additional participants
#'     (\code{projection$target_participants_historical_rate}), next to the
#'     theoretical line; see \code{show_historical_daris} in
#'     \code{\link{plot.tsa_cor}}.}
#'   \item{Estimated additional studies required}{The studies projection above.}
#' }
#' The two participants figures rest on different assumptions and can disagree
#' -- for example, when the studies supply less information per participant than
#' the single-study formula assumes, the historical-rate figure is larger, and it
#' can be positive when the theoretical figure is already 0. Both are
#' approximations.
#'
#' \emph{Caveat: fixed-effect information scale.} All projections are linear
#' extrapolations on the fixed-effect, study-level inverse-variance information
#' scale (\code{1/se_z^2}) that \code{tsa_cor()} compares with DARIS. They
#' deliberately do not model how random-effects weights or the between-study
#' variance (tau^2) would change as further studies are added, nor any change in
#' the size, precision or design of future studies; the mathematics are not
#' adjusted for random effects. Read every projected number of participants or
#' studies as indicative, not as a required quantity. This projection is
#' deliberately NOT called \dQuote{number of studies required} anywhere in the
#' output: that phrasing reads as deterministic, and it is not. Figures are
#' labelled \dQuote{Estimated ...}, and the printed output always carries a note
#' that the projection assumes future studies contribute information at
#' approximately the observed historical rate and is not a formal guarantee.
#' When the endpoint HAS already been reached, the projected fields
#' (\code{n_additional_studies}, \code{additional_participants_estimated} and
#' their variants) are \code{NA}: there is nothing left to project, and nothing
#' is printed or plotted for it.
#'
#' \strong{Order dependence and few studies.} TSA is a cumulative analysis, so
#' its result depends on the order of the studies: use \code{order_by} (or
#' supply the rows chronologically) and check the ordering warnings. With fewer
#' than 10 studies the heterogeneity estimate, \eqn{D^2} and hence DARIS are
#' unstable, and \code{tsa_cor()} warns; with a heavy-tailed or very
#' heterogeneous set of studies \eqn{D^2} may reach its numerical cap (see
#' \code{heterogeneity$D2_was_capped} under Value).
#'
#' @return An object of class \code{"tsa_cor"}: a list containing
#'   \describe{
#'     \item{\code{data}}{the (possibly re-ordered) input data with the added
#'       columns \code{z_fisher} (Fisher z of \code{r}), \code{se_z_ci}
#'       and \code{se_z_n} (the two standard-error versions; \code{se_z_ci}
#'       is \code{NA} without usable CI columns) and \code{se_z} (the one
#'       used, per \code{se_source}).}
#'     \item{\code{parameters}}{the design settings, including
#'       \code{target_r}, \code{r_anticipated}, \code{z_anticipated},
#'       \code{cor_type}, \code{se_source}, \code{spearman_variance},
#'       \code{var_factor}, \code{ci_level}, \code{method},
#'       \code{method_requested}, \code{re_inference} and
#'       \code{re_inference_requested} (the normalised options actually used
#'       and what the caller passed; they differ for aliases such as
#'       \code{"CO"} or \code{"knha"}, and when HKSJ falls back to
#'       \code{"standard"}). With a non-standard \code{re_inference},
#'       \code{res_re} is the corresponding \code{metafor::rma} fit,
#'       \code{cumulative} gains the columns \code{re_scale}, \code{re_df} and
#'       \code{zval}, and the summary table gains a \dQuote{Random-effects
#'       inference} row.}
#'     \item{\code{res_re}, \code{res_fe}}{the fitted random-effects and
#'       equal-effects \code{metafor::rma} objects (on the Fisher z scale).}
#'     \item{\code{pooled}}{the pooled random-effects correlation:
#'       \code{z}, \code{r}, \code{r_lb}, \code{r_ub} (back-transformed with
#'       \code{tanh()}) and \code{pval}.}
#'     \item{\code{heterogeneity}}{\code{Q}, \code{df}, \code{I2},
#'       \code{tau2}, \code{D2} (capped at 99.9\%), the uncapped
#'       \code{D2_raw}, \code{D2_was_capped} and the adjustment factor
#'       \code{AF}.}
#'     \item{\code{information_size}}{\code{z_alpha}, \code{z_beta},
#'       \code{info_required}, \code{RIS_participants}, \code{DARIS_info},
#'       \code{DARIS_participants}, \code{DARIS_info_threshold_n} (estimated
#'       cumulative participants at which the observed information reached
#'       DARIS, by interpolation), \code{route_endpoint_info},
#'       \code{route_endpoint_n}, \code{route_endpoint_participants_theoretical},
#'       \code{var_factor}, and the circularity flags
#'       \code{circularity_warning} (\code{target_r} unspecified) and
#'       \code{circularity_severe} (additionally, accrued participants exceed
#'       three times the resulting DARIS participant-equivalent).}
#'     \item{\code{cumulative}}{the cumulative analysis data frame: cumulative
#'       pooled Fisher z (\code{estimate}) with \code{se}, \code{Z},
#'       \code{ci.lb}/\code{ci.ub}, tau^2 etc. from \code{metafor::cumul()},
#'       the back-transformed \code{r_estimate}, \code{r_ci_lb} and
#'       \code{r_ci_ub}, \code{cum_n}, \code{info_accrued},
#'       \code{info_fraction} and the boundary columns.}
#'     \item{\code{boundary_timeline}}{the formal sequential boundary schedule
#'       including the synthetic \code{t = 1} final-analysis point.}
#'     \item{\code{results}}{\code{crossed_conventional} (evaluated over the
#'       FULL cumulative Z-curve, unlike \code{crossed_tsa}),
#'       \code{crossed_tsa}, \code{entered_futility_region}, the
#'       definitive-look fields, \code{daris_reached}, \code{final_reached},
#'       \code{final_tsa_look}, \code{participants_accrued} and
#'       \code{info_accrued_final}.}
#'     \item{\code{settings}}{\code{boundary_route}, \code{route_used},
#'       \code{legacy_fallback}, \code{fallback_used}, \code{fallback_route},
#'       \code{fallback_reason}, \code{route_endpoint},
#'       \code{route_endpoint_info} and \code{used_legacy_engine}.}
#'     \item{\code{projection}}{the retrospective projection (see Details):
#'       \code{method} (the \code{projection_stat} used),
#'       \code{info_per_participant_basis}, \code{I_required},
#'       \code{info_accrued}, \code{participants_accrued},
#'       \code{additional_info_required}, \code{central_info_increment},
#'       \code{central_participant_increment},
#'       \code{central_info_per_participant} (the historical rate used),
#'       \code{study_level_info_per_participant},
#'       \code{pooled_info_per_participant}, \code{n_additional_studies},
#'       \code{additional_participants_estimated} (continuous, information
#'       divided by information per participant, on the chosen basis),
#'       \code{additional_participants_study_level},
#'       \code{additional_participants_pooled},
#'       \code{participants_from_whole_studies},
#'       \code{target_participants_historical_rate} (accrued plus estimated
#'       additional participants), \code{additional_participants_theoretical},
#'       \code{n_studies} and \code{note}. The projected quantities are
#'       \code{NA} once the route's endpoint has been reached.}
#'     \item{\code{summary_table}}{a data frame with character
#'       \code{Parameter} and \code{Value} columns (\code{Value} is character so
#'       logical rows print as TRUE/FALSE/NA); the expansions of the short
#'       abbreviations used in \code{Parameter} are in
#'       \code{attr(summary_table, "abbreviations")}.}
#'   }
#'   Use \code{plot()}, \code{summary()}, or \code{print()} on the result.
#'
#' @references
#' Wetterslev J, Thorlund K, Brok J, Gluud C. "Estimating required
#' information size by quantifying diversity in random-effects model
#' meta-analyses." BMC Med Res Methodol. 2009;9:86.
#'
#' DerSimonian R, Laird N. "Meta-analysis in clinical trials." Control Clin
#' Trials. 1986;7:177-188.
#'
#' Hartung J, Knapp G. "On tests of the overall treatment effect in
#' meta-analysis with normally distributed responses." Stat Med.
#' 2001;20:1771-1782.
#'
#' Sidik K, Jonkman JN. "A simple confidence interval for meta-analysis."
#' Stat Med. 2002;21:3153-3164.
#'
#' Fisher RA. "On the 'probable error' of a coefficient of correlation
#' deduced from a small sample." Metron. 1921;1:3-32.
#'
#' Fieller EC, Hartley HO, Pearson ES. "Tests for rank correlation
#' coefficients. I." Biometrika. 1957;44:470-481.
#'
#' Bonett DG, Wright TA. "Sample size requirements for estimating Pearson,
#' Kendall and Spearman correlations." Psychometrika. 2000;65:23-28.
#'
#' @examples
#' \donttest{
#' path <- tsacor_example_data()
#' res <- tsa_cor(path, target_r = 0.10, order_by = "Year")
#' summary(res)
#' plot(res)
#' }
#'
#' \dontrun{
#' path <- tsacor_example_data()
#'
#' ## Sensitivity to the tau^2 estimator: D^2, DARIS and the boundaries change
#' res_reml <- tsa_cor(path, target_r = 0.10, order_by = "Year",
#'                     method = "REML", verbose = FALSE)
#' res_reml$heterogeneity$D2
#' res_reml$information_size$DARIS_participants
#'
#' ## Hartung-Knapp-Sidik-Jonkman inference and the analysis boundary route
#' res_hk <- tsa_cor(path, target_r = 0.10, order_by = "Year",
#'                   re_inference = "hksj", boundary_route = "analysis",
#'                   verbose = FALSE)
#' }
#'
#' @importFrom stats qnorm
#' @export
tsa_cor <- function(data,
                  alpha_two_sided = 0.05,
                  power = 0.80,
                  target_r = NA_real_,
                  cor_type = c("pearson", "spearman"),
                  se_source = c("ci", "n"),
                  spearman_variance = c("fieller", "bonett_wright"),
                  ci_level = 0.95,
                  method = "DL",
                  order_by = NULL,
                  verbose = TRUE,
                  boundary_route = c("design", "analysis"),
                  legacy_fallback = TRUE,
                  projection_stat = c("median", "mean"),
                  info_per_participant_basis = c("per_study", "pooled"),
                  re_inference = "standard") {

  ## Accept a couple of common shorthands/aliases for cor_type, normalising
  ## case and surrounding whitespace before match.arg(): "r" -> "pearson",
  ## "rho" -> "spearman" (so "R", " Rho ", "PEARSON", etc. all work too).
  if (is.character(cor_type) && length(cor_type) == 1L && !is.na(cor_type)) {
    cor_type <- tolower(trimws(cor_type))
    if (identical(cor_type, "rho")) cor_type <- "spearman"
    if (identical(cor_type, "r"))   cor_type <- "pearson"
  }
  cor_type <- match.arg(cor_type, c("pearson", "spearman"))
  se_source <- match.arg(se_source)
  spearman_variance <- match.arg(spearman_variance)
  boundary_route <- match.arg(boundary_route)
  projection_stat <- match.arg(projection_stat)
  info_per_participant_basis <- match.arg(info_per_participant_basis)
  if (!is.logical(legacy_fallback) || length(legacy_fallback) != 1L ||
      is.na(legacy_fallback))
    stop("legacy_fallback must be a single TRUE or FALSE")
  vcat <- function(...) if (verbose) cat(...)
  cor_label <- .tsacor_cor_label(cor_type)

  ## --- Random-effects (tau^2) estimator method --------------------------
  ## Validated up front, before any data handling. "FE" is intentionally NOT
  ## offered: the equal-effects comparator used for the Diversity (D^2)
  ## adjustment is always fitted internally with method = "FE", and allowing
  ## it for res_re would make D2 degenerate (0 by construction). "GENQ" /
  ## "GENQM" need a user-supplied `weights` argument that tsa_cor() does not
  ## collect. "CO" / "VC" are metafor's aliases of the Hedges estimator; they
  ## are normalised to "HE" here so behaviour does not depend on the metafor
  ## version (older releases reject a bare "CO").
  method_requested <- method
  if (is.character(method) && length(method) == 1L && !is.na(method) &&
      method %in% c("CO", "VC")) {
    method <- "HE"
  }
  valid_methods <- c("DL", "HE", "HS", "HSk", "SJ", "ML", "REML",
                      "EB", "PM", "PMM")
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("method must be a single character string; one of: ",
         paste(valid_methods, collapse = ", "), ".")
  }
  if (method %in% c("GENQ", "GENQM")) {
    stop("method = \"", method, "\" is not currently supported by tsa_cor(): ",
         "metafor's generalized-Q-statistic estimators require a ",
         "user-supplied `weights` argument to metafor::rma(), which ",
         "tsa_cor() does not currently collect or pass through. Supported ",
         "methods are: ", paste(valid_methods, collapse = ", "), ".")
  }
  if (!(method %in% valid_methods)) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "),
         " (the random-effects heterogeneity-variance estimators supported ",
         "by metafor::rma() that work without additional arguments this ",
         "package does not currently collect; see ?tsa_cor).")
  }

  ## --- Random-effects inference -----------------------------------------
  ## "standard" (default), "hksj" (alias "knha") or "hksj_adhoc" (alias
  ## "knha_adhoc"); case-insensitive. The user's own string is kept in
  ## `re_inference_requested`; `re_inference` holds the canonical value.
  re_inference_requested <- re_inference
  re_inference <- .tsacor_normalise_re_inference(re_inference)

  ## --- Scalar design-parameter validation -------------------------------
  if (!is.numeric(alpha_two_sided) || length(alpha_two_sided) != 1L ||
      !is.finite(alpha_two_sided) || alpha_two_sided <= 0 || alpha_two_sided >= 1) {
    stop("alpha_two_sided must be a single finite value strictly between 0 and 1.")
  }
  if (!is.numeric(power) || length(power) != 1L ||
      !is.finite(power) || power <= 0 || power >= 1) {
    stop("power must be a single finite value strictly between 0 and 1.")
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a single finite value strictly between 0 and 1.")
  }

  ## -----------------------------------------------------------------
  ## 1. Load / validate data
  ## -----------------------------------------------------------------
  if (is.character(data)) {
    data <- as.data.frame(readxl::read_excel(data))
  } else {
    data <- as.data.frame(data)
  }
  ## Excel headers commonly contain spaces; normalise them to underscores.
  ## This can MERGE two distinct headers into one name, after which a column
  ## lookup silently resolves to whichever came first -- refuse that case.
  names_before <- names(data)
  names(data) <- gsub(" ", "_", names(data))
  if (anyDuplicated(names(data))) {
    clashing <- unique(names(data)[duplicated(names(data))])
    stop("Column names are not unique after spaces are replaced with ",
         "underscores: ", paste(sQuote(clashing), collapse = ", "),
         ". Original column name(s) involved: ",
         paste(sQuote(names_before[names(data) %in% clashing]),
               collapse = ", "),
         ". Rename the columns in the source data so they remain ",
         "distinct once spaces become underscores.")
  }
  if (!is.null(order_by) && is.character(order_by) &&
      length(order_by) == 1L && !is.na(order_by)) {
    order_by <- gsub(" ", "_", order_by)
  }

  required_cols <- c("Study", "r", "n_subjects")
  if (identical(se_source, "ci")) required_cols <- c(required_cols, "lbound", "ubound")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "),
         if (identical(se_source, "ci") && any(c("lbound", "ubound") %in% missing_cols))
           paste0(" (the confidence-interval columns are needed because ",
                  "se_source = \"ci\"; use se_source = \"n\" to derive the ",
                  "standard errors from n_subjects instead)") else "")
  }
  if (any(is.na(data$Study)) ||
      any(!nzchar(trimws(as.character(data$Study))))) {
    stop("Study must contain non-missing, non-empty identifiers.")
  }
  if (anyDuplicated(data$Study)) {
    stop("Study names must be unique (duplicate found in 'Study' column).")
  }
  ## --- Minimum number of studies ---------------------------------------
  if (nrow(data) < 2L) {
    stop("tsa_cor() requires at least two studies.")
  }
  if (nrow(data) < 10L) {
    warning("Only ", nrow(data), " studies were supplied. Heterogeneity/D2 ",
            "estimates (and therefore DARIS and the monitoring boundaries) can ",
            "be highly unstable with few studies; the Copenhagen TSA manual ",
            "cautions about this below roughly 10 studies.")
  }

  ## --- Type / range sanity checks --------------------------------------
  ## Type first, then finiteness: is.finite() on a character or factor column
  ## returns all-FALSE, which would otherwise be misreported as "NA/NaN/Inf".
  numeric_cols <- c("r", "n_subjects")
  if (identical(se_source, "ci")) numeric_cols <- c(numeric_cols, "lbound", "ubound")
  not_numeric <- vapply(data[numeric_cols], function(col) !is.numeric(col),
                         logical(1))
  if (any(not_numeric)) {
    stop("Column(s) must be numeric, but are not: ",
         paste(sprintf("%s (%s)", numeric_cols[not_numeric],
                       vapply(data[numeric_cols[not_numeric]],
                              function(col) class(col)[1], character(1))),
               collapse = ", "),
         ". Check for text, footnote markers, or blank-but-not-empty cells ",
         "in the source data.")
  }
  if (any(!is.finite(data$r))) {
    stop("r must be finite for every study (found NA/NaN/Inf).")
  }
  if (any(abs(data$r) >= 1)) {
    stop("r must lie strictly between -1 and 1 for every study (Fisher's z ",
         "is infinite at |r| = 1).")
  }
  if (any(!is.finite(data$n_subjects))) {
    stop("n_subjects must be finite for every study (found NA/NaN/Inf).")
  }
  if (any(abs(data$n_subjects - round(data$n_subjects)) > sqrt(.Machine$double.eps))) {
    stop("n_subjects must be whole numbers.")
  }
  if (any(data$n_subjects <= 3)) {
    stop("n_subjects must be greater than 3 for every study (the variance ",
         "of Fisher's z is 1/(n - 3)).")
  }

  ## Standard errors of Fisher's z, both versions, study by study.
  data$z_fisher <- atanh(data$r)
  data$se_z_n   <- .tsacor_se_from_n(data$r, data$n_subjects, cor_type,
                                    spearman_variance)
  data$se_z_ci  <- rep(NA_real_, nrow(data))
  ci_cols_ok <- all(c("lbound", "ubound") %in% names(data)) &&
    is.numeric(data$lbound) && is.numeric(data$ubound) &&
    all(is.finite(data$lbound)) && all(is.finite(data$ubound)) &&
    all(abs(data$lbound) < 1) && all(abs(data$ubound) < 1) &&
    all(data$lbound < data$ubound)
  if (identical(se_source, "ci")) {
    ## Explicit, specific diagnostics for the CI columns actually used.
    if (any(!is.finite(data$lbound)) || any(!is.finite(data$ubound))) {
      stop("lbound and ubound must be finite for every study (found NA/NaN/Inf).")
    }
    if (any(abs(data$lbound) >= 1) || any(abs(data$ubound) >= 1)) {
      stop("lbound and ubound must lie strictly between -1 and 1 for every study.")
    }
    if (any(data$lbound >= data$ubound)) {
      stop("lbound must be smaller than ubound for every study.")
    }
  }
  if (ci_cols_ok) {
    data$se_z_ci <- .tsacor_se_from_ci(data$lbound, data$ubound, ci_level)
  }
  if (identical(se_source, "ci")) {
    data$se_z <- data$se_z_ci
    outside <- data$r < data$lbound | data$r > data$ubound
    if (any(outside)) {
      warning("r lies outside its own confidence interval [lbound, ubound] for ",
              "study/studies: ",
              paste(data$Study[outside], collapse = ", "),
              ". Check the columns; the standard errors are still derived ",
              "from the interval width.", call. = FALSE)
    }
    ## On the Fisher z scale a symmetric (Wald-type) interval is centred on z.
    ## A large departure means the interval was not built on the z scale
    ## (e.g. a bootstrap interval); the width-based standard error is then
    ## only an approximation.
    z_mid_dev <- abs((atanh(data$lbound) + atanh(data$ubound)) / 2 - data$z_fisher) /
      data$se_z_ci
    if (any(z_mid_dev > 0.25)) {
      warning("The confidence interval of study/studies ",
              paste(data$Study[z_mid_dev > 0.25], collapse = ", "),
              " is not centred on Fisher's z of r (the midpoint differs by ",
              "more than 0.25 standard errors), so it does not look like a ",
              "Fisher-z interval; the standard error taken from its width is ",
              "only approximate. Consider se_source = \"n\".", call. = FALSE)
    }
  } else {
    data$se_z <- data$se_z_n
  }
  if (any(!is.finite(data$se_z)) || any(data$se_z <= 0)) {
    stop("The standard error of Fisher's z must be finite and strictly ",
         "positive for every study.")
  }

  ## --- target_r ----------------------------------------------------------
  if (length(target_r) != 1L) {
    stop("target_r must be a single finite numeric value, or NA.")
  }
  if (!is.na(target_r) && (!is.numeric(target_r) || !is.finite(target_r))) {
    stop("target_r must be a single finite numeric value, or NA.")
  }
  if (!is.na(target_r) && abs(target_r) >= 1) {
    stop("target_r must lie strictly between -1 and 1.")
  }
  if (!is.na(target_r) && target_r == 0) {
    stop("target_r cannot equal 0: atanh(0)=0 makes the required information infinite.")
  }
  if (!is.na(target_r) && abs(target_r) < 0.10) {
    ## Not an error -- small correlations are not invalid, just increasingly
    ## information-hungry: required information is proportional to
    ## 1/[atanh(target_r)]^2, which grows rapidly as target_r -> 0.
    warning("target_r (", target_r, ") is very close to the null value of 0; ",
            "the required information size increases rapidly as atanh(target_r) ",
            "approaches zero. Confirm that this represents a meaningful ",
            "target correlation.", call. = FALSE)
  }

  ## Note: TSA is inherently order-dependent (it is a CUMULATIVE analysis).
  ## If order_by names a column the data are explicitly sorted by it
  ## (ascending); otherwise the SUPPLIED row order is used as-is and assumed
  ## to be chronological -- the package has no way to verify this.
  if (!is.null(order_by)) {
    if (!order_by %in% names(data)) {
      stop("order_by = '", order_by, "' is not a column in data.")
    }
    ob_col <- data[[order_by]]
    if (any(is.na(ob_col))) {
      warning("order_by column '", order_by, "' contains NA values; ",
              "those rows will sort to one end under R's default NA handling.")
    }
    if (!is.numeric(ob_col) && !inherits(ob_col, c("Date", "POSIXct", "POSIXt"))) {
      warning("order_by column '", order_by, "' is not numeric or a Date/POSIXct; ",
              "sorting will use R's default ordering for its type (e.g. lexical ",
              "for character), which may not reflect chronological order unless ",
              "e.g. formatted as 'YYYY-MM-DD'.")
    }
    if (anyDuplicated(ob_col[!is.na(ob_col)])) {
      warning("order_by column '", order_by, "' contains tied values; TSA is ",
              "order-dependent, so the relative order of tied studies (broken ",
              "by R's stable sort, i.e. their original row order among ties) ",
              "may affect the cumulative TSA. Consider a finer-grained ",
              "order_by column (e.g. publication date instead of year alone) ",
              "if the exact ordering of tied studies matters.")
    }
    data <- data[order(ob_col), , drop = FALSE]
    rownames(data) <- NULL
    vcat("Studies sorted by '", order_by, "' (ascending) for the cumulative analysis.\n", sep = "")
  } else {
    vcat("Study order used for sequential analysis (assumed chronological",
         " -- pass order_by = <column name> to sort explicitly):\n")
  }
  if (verbose) print(data$Study)
  vcat("\n")

  vcat(sprintf("Loaded %d studies (%s; Fisher z scale). Total participants: %.0f\n",
               nrow(data), cor_label, sum(data$n_subjects)))
  vcat(sprintf("Standard errors of Fisher's z taken from: %s\n",
               if (identical(se_source, "ci"))
                 sprintf("the reported %.0f%% confidence intervals", ci_level * 100)
               else "the sample sizes (n_subjects)"))
  if (verbose && ci_cols_ok) {
    ratio_ci_n <- data$se_z_ci / data$se_z_n
    vcat(sprintf(paste0("Median ratio of CI-based to n-based SE(z): %.3f ",
                        "(range %.3f-%.3f)\n"),
                 stats::median(ratio_ci_n), min(ratio_ci_n), max(ratio_ci_n)))
  }
  vcat("\n")

  ## -----------------------------------------------------------------
  ## 2. Conventional (overall) meta-analysis on Fisher's z
  ## -----------------------------------------------------------------
  ## `res_std` is the STANDARD random-effects fit (z-based inference). It
  ## always exists because the Diversity D^2 / DARIS / boundary calculations,
  ## I^2, tau^2 and Q are defined on it and must not depend on `re_inference`
  ## (the HKSJ adjustment rescales the variance of the pooled estimate, which
  ## would otherwise leak into D^2 = (Var_random - Var_fixed) / Var_random).
  res_std <- metafor::rma(yi = z_fisher, sei = se_z, data = data, method = method)
  res_fe  <- metafor::rma(yi = z_fisher, sei = se_z, data = data, method = "FE")

  ## `res_re`, the model reported back (pooled correlation, CI, p-value),
  ## carries the requested inference. Same tau^2 estimator and point estimate;
  ## only se/vb/test/CI differ. "hksj" -> metafor's test = "knha". "hksj_adhoc"
  ## -> variance multiplier max(1, q): for q >= 1 that is exactly test =
  ## "knha"; for q < 1 it is the unscaled variance with a t reference
  ## distribution, i.e. test = "t".
  if (identical(re_inference, "standard")) {
    res_re <- res_std
  } else {
    q_pooled <- .tsacor_hksj_q(data$z_fisher, data$se_z, res_std$tau2)
    if (!is.finite(q_pooled) || q_pooled <= 0) {
      warning("re_inference = \"", re_inference_requested, "\" needs a positive, ",
              "finite Hartung-Knapp scale factor, but it is ",
              format(q_pooled), " for these data (e.g. identical study estimates); ",
              "falling back to re_inference = \"standard\".", call. = FALSE)
      re_inference <- "standard"
      res_re <- res_std
    } else {
      test_arg <- if (identical(re_inference, "hksj_adhoc") && q_pooled < 1) "t" else "knha"
      res_re <- metafor::rma(yi = z_fisher, sei = se_z, data = data,
                             method = method, test = test_arg)
    }
  }
  ## Early-look caveat for HKSJ, raised as a plain warning() regardless of
  ## `verbose`, and only when HKSJ inference is actually used.
  if (re_inference %in% c("hksj", "hksj_adhoc")) {
    warning("HKSJ inference (re_inference = \"", re_inference_requested,
            "\") can be unstable at early cumulative looks because the ",
            "degrees of freedom are k-1: the HKSJ statistic is undefined ",
            "at k=1 (\"NA\" at the first look) and is based on only one ",
            "degree of freedom at k=2. Early cumulative HKSJ values should ",
            "not be interpreted as directly comparable in magnitude with ",
            "conventional normal Z-statistics.", call. = FALSE)
  }

  pooled_z  <- as.numeric(res_re$b)
  pooled_r  <- tanh(pooled_z)
  pooled_lb <- tanh(as.numeric(res_re$ci.lb))
  pooled_ub <- tanh(as.numeric(res_re$ci.ub))

  if (verbose) {
    if (identical(re_inference, "standard")) {
      cat(sprintf("=== Random-effects (%s) meta-analysis of Fisher's z ===\n",
                  .tsacor_method_label(method)))
    } else {
      cat(sprintf("=== Random-effects (%s) meta-analysis of Fisher's z; inference: %s ===\n",
                  .tsacor_method_label(method),
                  .tsacor_re_inference_label(re_inference)))
    }
    print(res_re)
    cat(sprintf("\nPooled %s (random effects, back-transformed): %.3f  95%% CI: %.3f-%.3f\n\n",
                cor_label, pooled_r, pooled_lb, pooled_ub))
  }

  ## -----------------------------------------------------------------
  ## 3. Heterogeneity / Diversity (D^2)
  ##    Wetterslev J, Thorlund K, Brok J, Gluud C. BMC Med Res Methodol.
  ##    2009;9:86. D^2 = (Var_random - Var_fixed) / Var_random.
  ## -----------------------------------------------------------------
  ## Taken from the STANDARD fit (res_std) so they do not depend on
  ## `re_inference`; for re_inference = "standard" res_std and res_re are the
  ## same object.
  Q    <- res_std$QE
  df   <- res_std$k - 1
  I2   <- res_std$I2
  tau2 <- res_std$tau2
  var_random <- res_std$vb[1, 1]
  var_fixed  <- res_fe$vb[1, 1]

  D2_raw <- max(0, (var_random - var_fixed) / var_random)
  ## D2 is bounded in [0,1) by definition. With very few studies and extreme
  ## heterogeneity it can approach 1 closely enough that 1/(1-D2) becomes
  ## numerically unstable/explosive, so it is capped defensively and a warning
  ## is raised. The 99.9% cap is a purely NUMERICAL safeguard against division
  ## by (near) zero, not a statistically justified correction. D2_raw is kept
  ## (uncapped) so a capped value is auditable rather than silent.
  D2_was_capped <- D2_raw >= 0.999
  D2 <- if (D2_was_capped) 0.999 else D2_raw
  if (D2_was_capped) {
    warning("Diversity D2 is at or very near its theoretical upper bound (100%), ",
            "indicating extreme heterogeneity relative to the number of studies. ",
            "D2 has been capped at 99.9% to avoid a numerically unstable/explosive ",
            "heterogeneity adjustment factor; interpret the required information ",
            "size and DARIS with caution in this scenario.", call. = FALSE)
  }
  AF <- 1 / (1 - D2)

  vcat("=== Heterogeneity ===\n")
  vcat(sprintf("Q = %.2f (df = %d), p = %.4f\n", Q, df, res_std$QEp))
  vcat(sprintf("I^2 = %.1f%%   tau^2 = %.4f (Fisher z scale)\n", I2, tau2))
  vcat(sprintf("Diversity D^2 = %.1f%%   Adjustment factor (1/(1-D2)) = %.3f\n\n", D2 * 100, AF))

  ## -----------------------------------------------------------------
  ## 4. Required information size (Fisher z scale)
  ## -----------------------------------------------------------------
  z_alpha <- stats::qnorm(1 - alpha_two_sided / 2)
  z_beta  <- stats::qnorm(power)

  if (is.na(target_r)) {
    r_anticipated <- pooled_r
    z_anticipated <- pooled_z
    effect_source <- "OBSERVED pooled correlation (random-effects model) -- see circularity caution in ?tsa_cor"
    warning("target_r was not specified: the OBSERVED pooled correlation from ",
            "this meta-analysis is being used to calculate the required ",
            "information size. This is circular (see ?tsa_cor, 'Details') and is ",
            "appropriate for exploratory use only -- for a publication-quality ",
            "TSA, set target_r to a pre-specified, meaningful correlation, ",
            "e.g. target_r = 0.20.", call. = FALSE)
  } else {
    r_anticipated <- target_r
    z_anticipated <- atanh(target_r)
    effect_source <- "user-specified target_r (pre-specified)"
  }
  if (!is.finite(z_anticipated) || z_anticipated == 0) {
    stop("The anticipated correlation is exactly 0 (the observed pooled ",
         "correlation is 0, or target_r = 0), which makes the required ",
         "information infinite. Specify a non-zero target_r.")
  }

  var_factor <- .tsacor_variance_factor(cor_type, spearman_variance, r_anticipated)

  info_required      <- (z_alpha + z_beta)^2 / z_anticipated^2
  RIS_participants   <- var_factor * info_required + 3
  DARIS_info         <- info_required * AF
  DARIS_participants <- var_factor * DARIS_info + 3

  if (verbose) {
    cat("=== Required Information Size (Fisher z scale) ===\n")
    cat("Effect size source:", effect_source, "\n")
    cat(sprintf("Anticipated %s (used for RIS calculation): %.3f (Fisher z = %.4f)\n",
                cor_label, r_anticipated, z_anticipated))
    cat(sprintf("alpha (2-sided) = %.3f, power = %.0f%%, variance factor c = %.3f\n",
                alpha_two_sided, power * 100, var_factor))
    cat(sprintf("Required statistical information (inverse-variance units): %.4f\n", info_required))
    cat(sprintf("Required sample size (RIS, no heterogeneity adj.; single-study equivalent, c*I + 3): %.0f\n",
                ceiling(RIS_participants)))
    cat(sprintf("Diversity-Adjusted Required Information (DARIS, information units): %.4f\n", DARIS_info))
    cat(sprintf("DARIS translated to an equivalent number of participants (c*DARIS + 3): %.0f\n\n",
                ceiling(DARIS_participants)))
    cat(sprintf("Total participants accrued across included studies: %.0f (%.1f%% of DARIS participant-equivalent)\n\n",
                sum(data$n_subjects), 100 * sum(data$n_subjects) / DARIS_participants))
  }

  ## Circularity flags -- two DIFFERENT facts:
  ##   - circularity_warning: the RIS was derived from the OBSERVED pooled
  ##     effect (target_r = NA). Circular regardless of how much accrued.
  ##   - circularity_severe: circular AND the accrued participants dwarf the
  ##     resulting DARIS participant-equivalent, the runaway case in which the
  ##     boundary collapses to the conventional one almost immediately.
  circularity_warning <- is.na(target_r)
  circularity_severe  <- circularity_warning &&
    sum(data$n_subjects) / DARIS_participants > 3
  if (verbose && circularity_severe) {
    cat("*** NOTE: accrued participants greatly exceed the DARIS because the RIS was\n")
    cat("    calculated from the observed (very large, very precise) pooled effect.\n")
    cat("    This is circular and will make the TSA boundary collapse almost\n")
    cat("    immediately to the conventional boundary. Consider re-running with\n")
    cat("    a pre-specified 'target_r' for a more standard, protocol-driven TSA. ***\n\n")
  } else if (verbose && circularity_warning) {
    cat("*** NOTE: 'target_r' was not specified, so the required information size\n")
    cat("    was calculated from the OBSERVED pooled effect. This is circular: the\n")
    cat("    required information size depends on the result it is being used to\n")
    cat("    evaluate. Set a pre-specified 'target_r' for a standard,\n")
    cat("    protocol-driven TSA. ***\n\n")
  }

  ## -----------------------------------------------------------------
  ## 5. Cumulative (sequential) meta-analysis
  ## -----------------------------------------------------------------
  ## The cumulative point estimates and tau^2 come from the STANDARD fit; for
  ## re_inference != "standard" the look-by-look inference is then replaced by
  ## its HKSJ counterpart (see .tsacor_cumulative_re_inference()).
  cumul_re <- metafor::cumul(res_std, order = seq_len(nrow(data)))
  cumul_df <- as.data.frame(cumul_re)

  cumul_df$Study <- data$Study
  cumul_df$cum_n <- cumsum(data$n_subjects)
  if (identical(re_inference, "standard")) {
    cumul_df$Z <- cumul_df$estimate / cumul_df$se
  } else {
    ## sets se/zval/pval/ci.lb/ci.ub, re_scale, re_df and Z (the
    ## normal-equivalent of the HKSJ t statistic; see ?tsa_cor)
    cumul_df <- .tsacor_cumulative_re_inference(cumul_df, yi = data$z_fisher,
                                              sei = data$se_z,
                                              re_inference = re_inference,
                                              method = method)
  }
  ## Back-transformed cumulative pooled correlation and its CI.
  cumul_df$r_estimate <- tanh(cumul_df$estimate)
  cumul_df$r_ci_lb    <- tanh(cumul_df$ci.lb)
  cumul_df$r_ci_ub    <- tanh(cumul_df$ci.ub)

  ## info_accrued: cumulative STUDY-LEVEL inverse-variance information,
  ## sum(1/se_z^2). This is a REPORTED-DATA surrogate for statistical
  ## information, not the exact Fisher information of the cumulative
  ## random-effects estimator (whose tau^2 is re-estimated at every look); the
  ## package nonetheless uses it as its DARIS/monitoring information scale,
  ## following common TSA practice -- see ?tsa_cor.
  cumul_df$info_accrued  <- cumsum(1 / data$se_z^2)
  cumul_df$info_fraction <- cumul_df$info_accrued / DARIS_info
  cumul_df$info_fraction_participants <- cumul_df$cum_n / DARIS_participants

  if (verbose) {
    print(cumul_df[, c("Study", "estimate", "se", "Z", "r_estimate", "cum_n",
                        "info_fraction", "info_fraction_participants")])
    cat("\n")
    cat("*** IMPORTANT CAVEAT: the cumulative Z-curve above is from a RANDOM-\n")
    cat("    EFFECTS model, whose between-study variance (tau^2) is RE-ESTIMATED\n")
    cat("    at every step. The Lan-DeMets/O'Brien-Fleming monitoring boundaries\n")
    cat("    strictly assume a fixed, CANONICAL information process with\n")
    cat("    independent Brownian-motion increments. A random-effects cumulative\n")
    cat("    Z-curve with tau^2 re-estimated at each look does not exactly satisfy\n")
    cat("    those assumptions, so applying the boundaries here is a widely-used\n")
    cat("    APPROXIMATION (as in the official Copenhagen Trial Unit TSA\n")
    cat("    software), not an exact result. ***\n\n")
  }

  ## -----------------------------------------------------------------
  ## 6. Trial sequential monitoring boundaries (alpha- and beta-spending)
  ##    Computed via the RTSA-derived compiled recursive-integration engine
  ##    (src/rtsa_core.h, driven by R/rtsa_engine.R), with a legacy R-only
  ##    fallback (R/obf_boundaries.R). The engine works on information
  ##    fractions only, so it is exactly the one used by tsahr.
  ## -----------------------------------------------------------------
  info_fracs <- cumul_df$info_fraction

  ## design-pass boundary timeline on the observed information fractions:
  ##   * retain observed interim looks strictly before DARIS (t < 1);
  ##   * append one definitive final-analysis point at t = 1 (HARIS/DARIS);
  ##   * do not continue the alpha or beta boundary through studies that
  ##     occur after the required information size has been reached.
  ##
  ## This is deliberately separate from `cumul_df`, which continues to
  ## contain every observed study and its cumulative Z-score.  Thus the
  ## evidence curve can extend beyond DARIS while the formal monitoring
  ## boundaries terminate at the RTSA-style final information point.
  ##
  ## ** 0.2.7.11: RTSA-derived compiled engine. **
  ## Alpha (efficacy) and beta (non-binding futility) boundaries are now
  ## computed by src/rtsa_core.h -- a C++ port of RTSA's own
  ## alpha_boundary()/beta_boundary() recursion -- driven through RTSA's
  ## boundaries(side = 2, futility = "non-binding", type = "design")
  ## orchestration (.rtsa_design_bounds() in R/rtsa_engine.R): alpha bounds
  ## on the (t < 1, 1) timeline, then the two-pass information-scale root
  ## search for the futility bounds against THAT timeline's own final
  ## efficacy wall (the alpha recursion's value at t = 1, e.g. 2.13 -- NOT
  ## qnorm(1 - alpha/2) = 1.96, which is what 0.2.6.x-0.2.7.10 substituted
  ## and which shifted design_R and with it every futility bound).
  ## See NEWS.md (0.2.7.11) and inst/REVERSE_ENGINEERING_RTSA.md for the
  ## numerical evidence.
  ##
  ## ** 0.2.7.13: ** this design pass is now ALWAYS run first, regardless of
  ## `boundary_route`, because its root (design_R) is what the "analysis"
  ## route is calibrated against. `legacy_fallback` (default TRUE) controls
  ## what happens if the compiled engine cannot produce a result (e.g. no
  ## root bracket exists for an unusual schedule): TRUE falls back to the
  ## pre-0.2.7.11 R-only approximate engine, with an impossible-to-miss
  ## warning and console banner, and marks the result
  ## (beta_engine$engine == "legacy_r_fallback"); FALSE fails closed --
  ## tsa_cor() stops with an error instead of silently substituting a
  ## non-RTSA-comparable engine. See ?tsa_cor, "legacy_fallback".
  boundary_timing_design <- sort(unique(c(info_fracs[info_fracs < 1], 1)))

  ## ** 0.2.7.14: ** explicit, programmatic record of any fallback, returned in
  ## `settings` (fallback_used / fallback_route / fallback_reason / route_used)
  ## so calling code never has to parse warning text to tell the cases apart:
  ##   fallback_route "none"   -- the requested route ran as requested;
  ##   fallback_route "design" -- boundary_route = "analysis" failed and the
  ##                              RTSA-derived DESIGN-route result is returned;
  ##   fallback_route "legacy" -- the RTSA-derived engine failed and the
  ##                              legacy, approximate R-only engine is used.
  fallback_used   <- FALSE
  fallback_route  <- "none"
  fallback_reason <- NA_character_
  route_used      <- boundary_route

  rtsa_fit <- tryCatch(
    .rtsa_design_bounds(boundary_timing_design, alpha = alpha_two_sided,
                        beta = 1 - power),
    error = function(e) e
  )
  used_legacy <- inherits(rtsa_fit, "error")
  if (used_legacy) {
    fb_msg <- paste0(
      "*** WARNING: THE RTSA-DERIVED BOUNDARY ENGINE FAILED (",
      conditionMessage(rtsa_fit), "). The alpha and futility boundaries in ",
      "this result were computed with the LEGACY, APPROXIMATE R-only engine ",
      "(pre-0.2.7.11) and are NOT comparable with RTSA. Do not report them ",
      "as RTSA-equivalent; check the design (information fractions, alpha, ",
      "power) or report the problem. ***"
    )
    if (!isTRUE(legacy_fallback)) {
      stop(paste0(
        "The RTSA-derived boundary engine failed (",
        conditionMessage(rtsa_fit), "), and legacy_fallback = FALSE means ",
        "tsa_cor() will not silently substitute the legacy, approximate ",
        "R-only engine. Set legacy_fallback = TRUE to allow that fallback ",
        "(with a warning), or address the underlying issue (check the ",
        "design: information fractions, alpha, power)."
      ), call. = FALSE)
    }
    warning(fb_msg, call. = FALSE, immediate. = TRUE)
    if (verbose) cat("\n", fb_msg, "\n\n", sep = "")
    legacy <- .tsacor_legacy_boundaries(info_fracs, boundary_timing_design,
                                       alpha_two_sided, 1 - power)
    alpha_bounds_design <- legacy$alpha_bounds_design
    beta_pre_daris      <- legacy$beta_pre_daris
    beta_engine <- c(legacy$beta_engine,
                     list(engine = "legacy_r_fallback",
                          engine_error = conditionMessage(rtsa_fit)))
    ## The legacy engine has no analysis-route counterpart; a request for
    ## boundary_route = "analysis" is honoured as closely as possible by
    ## falling back to the (also legacy) design-route result, noted below.
    if (boundary_route == "analysis" && verbose) {
      cat("Note: boundary_route = \"analysis\" was requested, but the legacy\n",
          "  fallback engine has no analysis-route equivalent; using its\n",
          "  design-route result instead.\n", sep = "")
    }
    route_endpoint      <- 1
    beta_final          <- utils::tail(alpha_bounds_design, 1)
    beta_bounds_design  <- c(beta_pre_daris, beta_final)
    boundary_timing     <- boundary_timing_design
    fallback_used       <- TRUE
    fallback_route      <- "legacy"
    fallback_reason     <- conditionMessage(rtsa_fit)
    route_used          <- "legacy"
  } else if (boundary_route == "design") {
    alpha_bounds_design <- rtsa_fit$alpha_ubound
    beta_pre_daris      <- rtsa_fit$beta_ubound[boundary_timing_design < 1]
    beta_engine <- c(rtsa_fit, list(engine = "rtsa_design_cpp",
                                    boundary = rtsa_fit$beta_ubound,
                                    warp_root = rtsa_fit$root))
    route_endpoint <- 1
    ## Final-look futility value (changed in 0.2.7.12): equal to the final
    ## efficacy bound, exactly as in RTSA's design pass, where the root
    ## search makes the futility bound meet the efficacy bound at t = 1.
    beta_final          <- utils::tail(alpha_bounds_design, 1)
    beta_bounds_design  <- c(beta_pre_daris, beta_final)
    boundary_timing     <- boundary_timing_design
  } else {
    ## boundary_route == "analysis": RTSA::RTSA(type = "analysis",
    ## design = NULL). The design pass above already gives design_R
    ## (rtsa_fit$root); build the analysis-pass timeline (observed looks
    ## capped/extended to design_R) exactly as RTSA's own RTSA() does, and
    ## run the analysis-route recursion (.rtsa_analysis_bounds()) against
    ## it. If THIS pass fails, the same legacy_fallback contract applies,
    ## but there is no legacy analysis-route engine to fall back to, so
    ## the fallback is the already-computed design-route result (with a
    ## warning explaining the substitution) rather than the pre-0.2.7.11
    ## approximate engine.
    design_R <- rtsa_fit$root
    t_ext <- if (max(info_fracs) < design_R) {
      c(info_fracs, design_R)
    } else if (max(info_fracs) > design_R) {
      c(info_fracs[info_fracs < design_R], design_R)
    } else {
      info_fracs
    }
    ana_fit <- tryCatch(
      .rtsa_analysis_bounds(t_ext, design_R, alpha = alpha_two_sided,
                            beta = 1 - power),
      error = function(e) e
    )
    if (inherits(ana_fit, "error")) {
      fb_msg2 <- paste0(
        "*** WARNING: THE RTSA ANALYSIS-ROUTE ENGINE FAILED (",
        conditionMessage(ana_fit), "). boundary_route = \"analysis\" was ",
        "requested, but tsa_cor() is falling back to its \"design\"-route ",
        "result instead (still the RTSA-derived compiled engine, just the ",
        "other route -- NOT the legacy R-only engine). ***"
      )
      if (!isTRUE(legacy_fallback)) {
        stop(paste0(
          "The RTSA analysis-route boundary engine failed (",
          conditionMessage(ana_fit), "), and legacy_fallback = FALSE means ",
          "tsa_cor() will not silently fall back to the design-route result. ",
          "Set legacy_fallback = TRUE to allow that fallback (with a ",
          "warning), pass boundary_route = \"design\" directly, or address ",
          "the underlying issue."
        ), call. = FALSE)
      }
      warning(fb_msg2, call. = FALSE, immediate. = TRUE)
      if (verbose) cat("\n", fb_msg2, "\n\n", sep = "")
      alpha_bounds_design <- rtsa_fit$alpha_ubound
      beta_pre_daris      <- rtsa_fit$beta_ubound[boundary_timing_design < 1]
      beta_engine <- c(rtsa_fit, list(engine = "rtsa_design_cpp",
                                      boundary = rtsa_fit$beta_ubound,
                                      warp_root = rtsa_fit$root,
                                      analysis_route_error = conditionMessage(ana_fit)))
      route_endpoint      <- 1
      beta_final          <- utils::tail(alpha_bounds_design, 1)
      beta_bounds_design  <- c(beta_pre_daris, beta_final)
      boundary_timing     <- boundary_timing_design
      fallback_used       <- TRUE
      fallback_route      <- "design"
      fallback_reason     <- conditionMessage(ana_fit)
      route_used          <- "design"
    } else {
      route_endpoint      <- design_R
      boundary_timing     <- ana_fit$timing
      alpha_bounds_design <- ana_fit$alpha_ubound
      beta_bounds_design  <- ana_fit$beta_ubound
      beta_pre_daris      <- beta_bounds_design[boundary_timing < route_endpoint]
      beta_engine <- c(ana_fit, list(engine = "rtsa_analysis_cpp",
                                     boundary = beta_bounds_design,
                                     design_R = design_R))
    }
  }

  ## -----------------------------------------------------------------
  ## 6b. Reconcile the participants-scale DARIS reference with the
  ##     observed-information "reached" verdict.
  ##
  ##     DARIS_participants (Section 4) is a THEORETICAL translation of the
  ##     required information into a single-study-equivalent number of
  ##     participants (c * DARIS + 3). The observed-information criterion
  ##     (info_accrued vs DARIS_info) is what `final_reached`/`daris_reached`
  ##     are based on. The cumulative-participants point at which the
  ##     observed information first reaches the target is therefore located by
  ##     linear interpolation between the two bracketing looks; because no
  ##     study actually occurred exactly there it is an ESTIMATE, deliberately
  ##     NOT called "DARIS participants" (that label is reserved for the
  ##     theoretical quantity). plot.tsa_cor() shows BOTH, separately labelled.
  ##
  ##     Two DISTINCT thresholds are kept apart:
  ##       * DARIS itself (info_fraction = 1) -- `daris_reached` /
  ##         `DARIS_info_threshold_n`;
  ##       * the ROUTE ENDPOINT (route_endpoint * DARIS_info; = DARIS for the
  ##         design route, design_R * DARIS for the analysis route) --
  ##         `final_reached` / `route_endpoint_n_est`. THIS is the formal
  ##         analysis endpoint the decision layer is evaluated at.
  ## -----------------------------------------------------------------
  final_reached <- max(info_fracs) >= route_endpoint
  daris_reached <- max(info_fracs) >= 1
  analysis_endpoint <- identical(route_used, "analysis")
  route_endpoint_info <- route_endpoint * DARIS_info
  ## theoretical single-study-equivalent participants of the route endpoint
  route_endpoint_participants_theoretical <- var_factor * route_endpoint_info + 3

  DARIS_info_threshold_n <- .tsacor_participants_at_fraction(cumul_df, 1)
  route_endpoint_n_est <- if (route_endpoint == 1) {
    DARIS_info_threshold_n
  } else {
    .tsacor_participants_at_fraction(cumul_df, route_endpoint)
  }
  endpoint_name <- if (analysis_endpoint) {
    sprintf("analysis-route endpoint (%.3f x DARIS)", route_endpoint)
  } else {
    "DARIS"
  }

  ## Map only genuine pre-endpoint observed looks back to the cumulative study
  ## table. The synthetic t = 1 HARIS point is stored separately in
  ## `boundary_timeline`; it is not an observed study and therefore must not
  ## be inserted into `cumul_df`.
  boundary_z <- rep(NA_real_, length(info_fracs))
  futility_z <- rep(NA_real_, length(info_fracs))
  pre_daris <- info_fracs < route_endpoint
  if (any(pre_daris)) {
    boundary_z[pre_daris] <- alpha_bounds_design[match(info_fracs[pre_daris],
                                                       boundary_timing)]
    futility_z[pre_daris] <- beta_bounds_design[match(info_fracs[pre_daris],
                                                       boundary_timing)]
  }

  cumul_df$TSA_boundary_upper <- boundary_z
  cumul_df$TSA_boundary_lower <- -boundary_z
  cumul_df$TSA_futility_upper <- futility_z
  cumul_df$TSA_futility_lower <- -futility_z

  ## X-coordinate of the synthetic t = 1 boundary. When the observed
  ## cumulative information has actually reached the route endpoint, the
  ## formal final boundary terminates at that OBSERVED-information threshold
  ## (interpolated participants), not at the theoretical participant
  ## equivalent, which remains a separate reference on the plot. If the
  ## endpoint has not been reached, RTSA's retrospective design endpoint
  ## remains the theoretical one.
  boundary_endpoint_n <- if (final_reached && is.finite(route_endpoint_n_est)) {
    route_endpoint_n_est
  } else {
    route_endpoint_participants_theoretical
  }

  boundary_timeline <- data.frame(
    info_fraction = boundary_timing,
    cum_n = c(
      vapply(boundary_timing[boundary_timing < route_endpoint], function(tt) {
        idx <- which(info_fracs == tt)[1]
        cumul_df$cum_n[idx]
      }, numeric(1)),
      boundary_endpoint_n
    ),
    TSA_boundary_upper = alpha_bounds_design,
    TSA_boundary_lower = -alpha_bounds_design,
    TSA_futility_upper = beta_bounds_design,
    TSA_futility_lower = -beta_bounds_design,
    synthetic = boundary_timing == route_endpoint,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("=== Trial sequential monitoring boundaries (alpha/beta spending) ===\n")
    if (analysis_endpoint) {
      cat(sprintf(paste0("Formal final boundary endpoint (%s = %.4f information ",
                         "units): %.1f cumulative participants\n"),
                  endpoint_name, route_endpoint_info, boundary_endpoint_n))
    } else {
      cat(sprintf("Formal final boundary endpoint (DARIS information reached): %.1f cumulative participants\n",
                  boundary_endpoint_n))
    }
    print(cumul_df[, c("Study", "info_fraction", "Z", "TSA_boundary_upper", "TSA_futility_upper")])
    cat("\n")
  }

  ## -----------------------------------------------------------------
  ## 7c. Formal decision logic: restrict to looks up to and including
  ##     DARIS being first reached.
  ##
  ##     Once DARIS is reached, the planned monitoring process is complete.
  ##     The first DARIS-reaching observed look is therefore compared with
  ##     the definitive t = 1 boundary, while later observed studies are
  ##     not treated as additional final looks.
  ##
  ##     The formal "crossed_tsa" / "entered_futility_region" verdicts are
  ##     therefore evaluated only through the FIRST look at which DARIS
  ##     was reached (or through all looks, if DARIS was never reached).
  ##     The full cumulative Z-curve, including any studies added after
  ##     DARIS, is still returned/plotted in full -- only the formal
  ##     sequential decision is restricted, not what is shown.
  final_tsa_look <- if (final_reached) which(info_fracs >= route_endpoint)[1] else length(info_fracs)
  decision_idx <- seq_len(final_tsa_look)

  ## For the formal decision, the first look reaching the route endpoint
  ## (1 for the design route; design_R for the analysis route) is compared
  ## with the definitive boundary at that endpoint, even though that
  ## boundary is displayed at the synthetic HARIS endpoint rather than on
  ## the observed post-endpoint study rows.
  decision_boundary_upper <- boundary_timeline$TSA_boundary_upper[
    match(pmin(info_fracs[decision_idx], route_endpoint), boundary_timeline$info_fraction)
  ]
  decision_futility_upper <- boundary_timeline$TSA_futility_upper[
    match(pmin(info_fracs[decision_idx], route_endpoint), boundary_timeline$info_fraction)
  ]

  crossed_tsa <- any(abs(cumul_df$Z[decision_idx]) >= decision_boundary_upper,
                     na.rm = TRUE)
  ## `crossed_conventional` is DELIBERATELY evaluated over the FULL
  ## cumulative Z-curve (unlike crossed_tsa/entered_futility_region,
  ## which are restricted to decision_idx). This is intentional: it
  ## answers "did the naive, uncorrected cumulative P-value ever drop
  ## below 0.05?", including at studies added after DARIS was reached,
  ## as a deliberate contrast against the properly-scoped crossed_tsa --
  ## illustrating exactly the repeated-testing inflation risk that TSA
  ## exists to guard against. If you instead want "did the formal
  ## sequential analysis reach conventional significance at/before the
  ## DARIS decision look" (the same horizon as crossed_tsa), use
  ## `any(abs(cumul_df$Z[decision_idx]) >= z_alpha, na.rm = TRUE)`
  ## instead. Decide explicitly before changing this -- both readings
  ## are defensible, but they answer different questions.
  ## 0.2.8.12: na.rm = TRUE guards against cumul_df$Z[1] == NA at k = 1
  ## under HKSJ inference (see .tsacor_cumulative_re_inference()); with a
  ## single look otherwise NA, any(NA) would propagate to NA here.
  crossed_conventional <- any(abs(cumul_df$Z) >= z_alpha, na.rm = TRUE)
  ## `entered_futility_region` uses the same decision_idx / definitive
  ## t=1-boundary comparison as crossed_tsa (see comment above
  ## decision_idx). At the first DARIS-reaching look specifically, this
  ## means comparing against the FINAL futility boundary (which, as in
  ## RTSA's design pass, equals the final efficacy boundary), not an interim futility
  ## boundary -- so entered_futility_region == TRUE at that look does
  ## NOT mean "TSA recommends stopping for futility now"; the printed
  ## summary label says "(not a formal stopping decision)" for exactly
  ## this reason. That caveat lives only in the printed/summary text,
  ## not in the boolean itself -- a caller reading
  ## `results$entered_futility_region` programmatically (bypassing the
  ## printed output) will not see it. See ?tsa_cor, "Value", for the
  ## corresponding caveat in the documented return value.
  entered_futility_region <- any(abs(cumul_df$Z[decision_idx]) <= decision_futility_upper,
                                  na.rm = TRUE)

  ## ** 0.2.7.14: DEFINITIVE-LOOK fields (fix of the 0.2.7.13 semantics). **
  ##
  ## `crossed_tsa` and `entered_futility_region` above are "at ANY formal
  ## look up to the route endpoint" quantities. 0.2.7.13 defined
  ## `final_non_efficacy <- !crossed_tsa`, which is NOT "the definitive look
  ## did not cross efficacy": a trial that crossed at an interim look and
  ## then fell back below the boundary at the definitive look would report
  ## FALSE. The fields below refer to the definitive look ONLY --
  ## `final_tsa_look`, the first look reaching the route endpoint -- and are
  ## NA when that endpoint has not been reached (there is then no
  ## definitive look to report on).
  ##
  ## Because the futility and efficacy boundaries meet at the definitive
  ## look, `final_entered_futility_region` there is the complement of
  ## `final_crossed_efficacy` (both would be TRUE only at |Z| exactly equal
  ## to the boundary). Neither is a recommendation to stop early.
  definitive <- .tsacor_definitive_look(cumul_df$Z, final_reached, final_tsa_look,
                                       decision_boundary_upper,
                                       decision_futility_upper)
  final_crossed_efficacy        <- definitive$final_crossed_efficacy
  final_entered_futility_region <- definitive$final_entered_futility_region
  final_non_efficacy            <- definitive$final_non_efficacy

  ## -----------------------------------------------------------------
  ## 6d. Retrospective projection: estimated additional studies/participants
  ##     needed to reach the route's own target information, based directly
  ##     on the OBSERVED study-level information increments (see "Estimated
  ##     additional studies/participants" under ?tsa_cor, Details).
  ##
  ##     I_required is computed per route and never conflated:
  ##       - design route:   I_required = DARIS            (DARIS_info)
  ##       - analysis route: I_required = design_R x DARIS (route_endpoint_info)
  ##     `central_info_increment` / `central_participant_increment` are the
  ##     `projection_stat` (median by default, mean if requested) of each
  ##     included study's own information (1/se_z^2) and sample size -- i.e. a
  ##     single "typical future study" estimated from the studies already in
  ##     the data, not from any new assumption. `n_additional_studies` is
  ##     always rounded UP to a natural number (>= 1 whenever there is a
  ##     genuine positive shortfall).
  ## -----------------------------------------------------------------
  participants_accrued <- sum(data$n_subjects)
  info_accrued_final   <- cumul_df$info_accrued[nrow(cumul_df)]

  per_study_info_increment <- 1 / data$se_z^2
  per_study_participants   <- data$n_subjects
  n_studies_total <- length(per_study_participants)
  central_fun <- if (identical(projection_stat, "mean")) mean else stats::median
  central_info_increment        <- central_fun(per_study_info_increment)
  central_participant_increment <- central_fun(per_study_participants)

  ## Information per participant, study by study. Additional PARTICIPANTS are
  ## the primary projected quantity and come straight from the information
  ## shortfall (continuous, not through a whole number of studies); the number
  ## of studies is a separate, secondary translation of the same shortfall.
  ## Because n_subjects > 3 and se_z is finite and positive for every study
  ## (validated above), every ratio is finite and positive.
  per_study_info_per_participant <- per_study_info_increment / per_study_participants
  study_level_info_per_participant <- central_fun(per_study_info_per_participant)
  ## ratio of sums (dominated by the larger studies), i.e. the slope of
  ## observed cumulative information against cumulative participants
  pooled_info_per_participant <- sum(per_study_info_increment) / sum(per_study_participants)
  ## The historical rate actually used for the participants projection; both
  ## variants are always returned so the choice can be checked.
  central_info_per_participant <- if (identical(info_per_participant_basis, "pooled")) {
    pooled_info_per_participant
  } else {
    study_level_info_per_participant
  }
  info_per_participant_basis_label <- if (identical(info_per_participant_basis, "pooled")) {
    "pooled ratio: total information / total participants"
  } else {
    sprintf("%s of the study-level information per participant", projection_stat)
  }
  ## short form for the summary-table row labels ("IPP" is spelled out in the
  ## abbreviations legend attached to summary_df below)
  info_per_participant_basis_short <- if (identical(info_per_participant_basis, "pooled")) {
    "pooled"
  } else {
    projection_stat
  }

  I_required <- route_endpoint_info
  additional_info_required_raw <- I_required - info_accrued_final

  n_additional_studies_est            <- NA_integer_
  additional_participants_estimated   <- NA_real_
  additional_participants_study_level <- NA_real_
  additional_participants_pooled      <- NA_real_
  participants_from_whole_studies     <- NA_real_
  target_participants_historical_rate <- NA_real_
  if (!final_reached && is.finite(additional_info_required_raw) &&
      additional_info_required_raw > 0) {
    ## Primary: participants = remaining information / information per
    ## participant (continuous; rounded up only for display).
    if (is.finite(study_level_info_per_participant) && study_level_info_per_participant > 0) {
      additional_participants_study_level <- additional_info_required_raw /
        study_level_info_per_participant
    }
    if (is.finite(pooled_info_per_participant) && pooled_info_per_participant > 0) {
      additional_participants_pooled <- additional_info_required_raw /
        pooled_info_per_participant
    }
    additional_participants_estimated <- if (identical(info_per_participant_basis, "pooled")) {
      additional_participants_pooled
    } else {
      additional_participants_study_level
    }
    ## Cumulative-participants position at which the route's target
    ## information would be reached at the historical rate (drawn by plot()
    ## as the "historical rate" reference line).
    if (is.finite(additional_participants_estimated)) {
      target_participants_historical_rate <- participants_accrued +
        additional_participants_estimated
    }
    ## Secondary: whole number of typical studies (rounded UP, >= 1 here
    ## because a genuine shortfall exists) and the participants those whole
    ## studies would bring -- a different question from the one above.
    if (is.finite(central_info_increment) && central_info_increment > 0) {
      n_additional_studies_est <- max(1L, ceiling(additional_info_required_raw /
                                                    central_info_increment))
      if (is.finite(central_participant_increment)) {
        participants_from_whole_studies <- n_additional_studies_est *
          central_participant_increment
      }
    }
  }

  ## THEORETICAL additional participants: the direct, deterministic
  ## difference between the single-study-equivalent participant target of the
  ## route endpoint (c * route_endpoint * DARIS + 3) and the participants
  ## accrued. Reported next to the historical-rate estimate above; the two
  ## rest on different assumptions and can disagree.
  additional_participants_theoretical <- if (!final_reached) {
    max(ceiling(route_endpoint_participants_theoretical - participants_accrued), 0)
  } else {
    NA_real_
  }

  projection_note <- paste0(
    "Projection assumes future studies contribute information at ",
    "approximately the observed historical rate (", info_per_participant_basis_label,
    "); it is not a formal guarantee of the number of future studies or ",
    "participants required. It is a linear extrapolation on the fixed-effect, ",
    "study-level inverse-variance information scale that is compared with ",
    "DARIS, and it does not model how random-effects weights or the ",
    "between-study variance (\u03c4\u00b2) would change as further studies ",
    "are added; treat it as indicative only."
  )

  projection <- list(
    method = projection_stat,
    info_per_participant_basis = info_per_participant_basis,
    I_required = I_required,
    info_accrued = info_accrued_final,
    participants_accrued = participants_accrued,
    additional_info_required = if (!final_reached) max(additional_info_required_raw, 0) else NA_real_,
    central_info_increment = central_info_increment,
    central_participant_increment = central_participant_increment,
    central_info_per_participant = central_info_per_participant,
    study_level_info_per_participant = study_level_info_per_participant,
    pooled_info_per_participant = pooled_info_per_participant,
    n_additional_studies = n_additional_studies_est,
    additional_participants_estimated = additional_participants_estimated,
    additional_participants_study_level = additional_participants_study_level,
    additional_participants_pooled = additional_participants_pooled,
    participants_from_whole_studies = participants_from_whole_studies,
    target_participants_historical_rate = target_participants_historical_rate,
    additional_participants_theoretical = additional_participants_theoretical,
    n_studies = n_studies_total,
    note = projection_note
  )

  if (verbose) {
    if (final_reached && final_tsa_look < nrow(cumul_df)) {
      cat(sprintf(paste0("Note: %s was reached at study #%d of %d ('%s'). Formal TSA\n",
                          "  boundary-crossing/futility decisions below are evaluated only\n",
                          "  through that look (studies added afterward are still shown in\n",
                          "  the returned data and plot, but are not treated as additional\n",
                          "  formal '%s' analyses -- see ?tsa_cor).\n"),
                  endpoint_name, final_tsa_look, nrow(cumul_df),
                  cumul_df$Study[final_tsa_look],
                  if (analysis_endpoint) sprintf("t=%.3f", route_endpoint) else "t=1"))
    }
    cat(sprintf("Cumulative Z-curve crossed the conventional (P<0.05) boundary: %s\n",
                ifelse(crossed_conventional, "YES", "NO")))
    cat(sprintf("Cumulative Z-curve crossed the TSA monitoring boundary       : %s\n",
                ifelse(crossed_tsa, "YES", "NO")))
    cat(sprintf("Cumulative Z-curve entered the non-binding futility region  : %s\n",
                ifelse(entered_futility_region, "YES", "NO")))
    if (final_reached) {
      cat(sprintf("Definitive look (%s) crossed the efficacy boundary%s: %s\n",
                  if (analysis_endpoint) "route endpoint" else "DARIS",
                  if (analysis_endpoint) "" else "         ",
                  ifelse(final_crossed_efficacy, "YES", "NO")))
    }
    cat(sprintf("Required information size (DARIS) reached                    : %s\n",
                ifelse(daris_reached, "YES", "NO")))
    if (analysis_endpoint) {
      cat(sprintf("Analysis-route endpoint (%.3f x DARIS) reached               : %s\n",
                  route_endpoint, ifelse(final_reached, "YES", "NO")))
    }
    cat(sprintf("  Theoretical DARIS participant-equivalent (c*DARIS + 3)       : %.0f\n",
                ceiling(DARIS_participants)))
    if (analysis_endpoint) {
      cat(sprintf(paste0("  Theoretical participant-equivalent of the analysis-route endpoint\n",
                          "  (%.3f x DARIS)                                            : %.0f\n"),
                  route_endpoint, ceiling(route_endpoint_participants_theoretical)))
    }
    if (daris_reached) {
      cat(sprintf(paste0("  Estimated cumulative participants at which DARIS information\n",
                          "  was reached (interpolated, not an observed look)          : %.0f\n"),
                  ceiling(DARIS_info_threshold_n)))
      if (abs(DARIS_info_threshold_n - DARIS_participants) > 0.01 * DARIS_participants) {
        cat(paste0("  (These may differ because the observed study-level information\n",
                   "   per participant differs from the 1/c assumed by the theoretical\n",
                   "   participant-equivalent, and because that equivalent counts the\n",
                   "   3 lost units of information only once -- in either direction, not\n",
                   "   necessarily because information accrued faster.)\n"))
      }
    }
    cat("\n")
    ## Estimated additional studies/participants (6d) -- printed only when the
    ## route's own target has NOT been reached; design and analysis routes use
    ## different labels/figures.
    if (!final_reached) {
      target_lab <- if (analysis_endpoint) "analysis-route endpoint" else "DARIS"
      cat(sprintf("Participants accrued = %.0f\n\n", participants_accrued))
      if (analysis_endpoint) {
        cat(sprintf(
          paste0("Analysis-route endpoint (%.3f x DARIS) = %s information units): ",
                 "%s cumulative participants (theoretical)\n"),
          route_endpoint,
          formatC(route_endpoint_info, format = "f", digits = 4, big.mark = ","),
          formatC(boundary_endpoint_n, format = "f", digits = 0, big.mark = ",")))
        cat(sprintf("Additional information required: %s\n",
                    formatC(projection$additional_info_required, format = "f",
                            digits = 3, big.mark = ",")))
      }
      cat(sprintf("Theoretical additional participants to %s (c*I + 3 equivalent): %.0f\n",
                  target_lab, additional_participants_theoretical))
      ## Where (in cumulative participants) the target would be reached at the
      ## historical information/participant rate -- the same position plot()
      ## draws as the historical-rate reference line.
      if (is.finite(target_participants_historical_rate)) {
        cat(sprintf(paste0("%s: ~%s cumulative participants\n",
                           "  (%s information units per participant)\n"),
                    if (analysis_endpoint) "Historical information/participant-rate projection"
                    else "DARIS (historical rate)",
                    formatC(ceiling(target_participants_historical_rate), format = "d", big.mark = ","),
                    formatC(central_info_per_participant, format = "f", digits = 4)))
      }
      if (is.na(additional_participants_estimated)) {
        cat(sprintf("Estimated additional participants to %s (historical rate): cannot be estimated\n",
                    target_lab))
      } else {
        cat(sprintf("Estimated additional participants to %s (historical rate): ~%s\n",
                    target_lab,
                    formatC(ceiling(additional_participants_estimated), format = "d", big.mark = ",")))
      }
      if (is.na(n_additional_studies_est)) {
        cat("Estimated additional studies required: cannot be estimated\n")
        cat("  (no usable historical per-study information increment to project from)\n")
      } else {
        cat(sprintf("Estimated additional studies required: %d\n",
                    n_additional_studies_est))
      }
      cat(strwrap(paste0("Note: ", projection_note), width = 78,
                  prefix = "", initial = ""), sep = "\n")
      cat("\n")
    }
    cat(strwrap(paste0(
      "Note: The \u03c4\u00b2 estimator may have limited influence on the pooled ",
      "average effect-size when the evidence base is substantial, but it can ",
      "materially influence heterogeneity-dependent quantities, prediction ",
      "intervals, DARIS, and the timing of TSA conclusions\u2014particularly ",
      "when cumulative information is near the DARIS threshold."
    ), width = 78, prefix = "", initial = ""), sep = "\n")
    cat("\n")
  }

  ## -----------------------------------------------------------------
  ## 7. Summary table
  ## -----------------------------------------------------------------
  ## Decision rows distinguish "at ANY formal look" from "at the DEFINITIVE
  ## look"; the analysis-route endpoint rows are added only when that route
  ## actually ran. Parameter labels are kept short (abbreviations expanded
  ## once in the "Abbreviations" legend returned as attr(summary_table,
  ## "abbreviations"), printed by summary.tsa_cor(), and in ?tsa_cor, Value).
  sum_par <- c(sprintf("Pooled %s (RE, observed)", cor_label), "95% CI lower", "95% CI upper",
               sprintf("Anticipated %s (for RIS)", cor_label),
               "Pooled Fisher z (RE)",
               "I2 (%)", "tau2 (z scale)", "Diversity D2 (%)", "Adjustment factor",
               "Variance factor c",
               "Required info (inv-var units)",
               "RIS, participants (c*I + 3)",
               "DARIS (info units)",
               "DARIS participant-equiv. (c*DARIS + 3)",
               "Participants accrued",
               "Info accrued (observed inv-var)",
               "% of DARIS reached",
               "Participants at DARIS reached (est.)")
  sum_val <- c(round(pooled_r, 3), round(pooled_lb, 3), round(pooled_ub, 3),
               round(r_anticipated, 3),
               round(pooled_z, 4),
               round(I2, 1), round(tau2, 4), round(D2 * 100, 1), round(AF, 3),
               round(var_factor, 3),
               round(info_required, 4),
               ceiling(RIS_participants),
               round(DARIS_info, 4),
               ceiling(DARIS_participants),
               participants_accrued,
               round(info_accrued_final, 4),
               round(100 * info_accrued_final / DARIS_info, 1),
               ifelse(is.na(DARIS_info_threshold_n), NA_real_,
                      ceiling(DARIS_info_threshold_n)))
  if (analysis_endpoint) {
    sum_par <- c(sum_par, sprintf(
      "Participants at AR endpoint (%.3fxDARIS) reached (est.)",
      route_endpoint))
    sum_val <- c(sum_val, ifelse(is.na(route_endpoint_n_est), NA_real_,
                                 ceiling(route_endpoint_n_est)))
  }
  ## Estimated additional studies/participants rows -- appended only when the
  ## route's own target has NOT been reached, matching the printed output.
  if (!analysis_endpoint && !daris_reached) {
    sum_par <- c(sum_par,
                 "Add'l participants to DARIS (theoretical)",
                 "DARIS reached, hist. rate (est. participants)",
                 sprintf("Add'l participants to DARIS (hist. rate; %s IPP)",
                         info_per_participant_basis_short),
                 sprintf("Add'l studies to DARIS (%s-based proj.)",
                        projection_stat))
    sum_val <- c(sum_val,
                 additional_participants_theoretical,
                 ifelse(is.finite(target_participants_historical_rate),
                        ceiling(target_participants_historical_rate), NA_real_),
                 ifelse(is.na(additional_participants_estimated), NA_real_,
                        ceiling(additional_participants_estimated)),
                 ifelse(is.na(n_additional_studies_est), NA_real_, n_additional_studies_est))
  } else if (analysis_endpoint && !final_reached) {
    sum_par <- c(sum_par,
                 sprintf("Add'l info to AR endpoint (%.3fxDARIS)",
                        route_endpoint),
                 "Add'l participants to AR endpoint (theoretical)",
                 "AR endpoint reached, hist. rate (est. participants)",
                 sprintf("Add'l participants to AR endpoint (hist. rate; %s IPP)",
                         info_per_participant_basis_short),
                 sprintf("Add'l studies to AR endpoint (%s-based proj.)",
                        projection_stat))
    sum_val <- c(sum_val,
                 round(projection$additional_info_required, 4),
                 additional_participants_theoretical,
                 ifelse(is.finite(target_participants_historical_rate),
                        ceiling(target_participants_historical_rate), NA_real_),
                 ifelse(is.na(additional_participants_estimated), NA_real_,
                        ceiling(additional_participants_estimated)),
                 ifelse(is.na(n_additional_studies_est), NA_real_, n_additional_studies_est))
  }
  sum_par <- c(sum_par,
               "Crossed conventional boundary",
               "Crossed TSA boundary (any formal look)",
               "Entered futility region (any look; not a stop decision)",
               "Definitive look crossed efficacy (NA if not reached)",
               "Definitive look: non-efficacy (NA if not reached)")
  ## `Value` is a CHARACTER column: building it with c() of numbers and
  ## logicals would silently coerce everything to numeric, so the TRUE/FALSE
  ## decision rows would print as 1/0. Numbers keep their rounding; logicals
  ## print as TRUE/FALSE/NA; NA stays NA.
  sum_val <- c(.tsacor_format_numeric(sum_val),
               .tsacor_format_logical(c(crossed_conventional, crossed_tsa,
                                       entered_futility_region,
                                       final_crossed_efficacy, final_non_efficacy)))
  summary_df <- data.frame(Parameter = sum_par, Value = sum_val,
                           stringsAsFactors = FALSE)
  ## Name the inference option right below the pooled correlation and its CI
  ## when it is not the standard one (the default table is unchanged).
  if (!identical(re_inference, "standard")) {
    summary_df <- rbind(
      summary_df[1:3, , drop = FALSE],
      data.frame(Parameter = "Random-effects inference",
                 Value = .tsacor_re_inference_label(re_inference),
                 stringsAsFactors = FALSE),
      summary_df[-(1:3), , drop = FALSE])
    rownames(summary_df) <- NULL
  }
  ## The single place the abbreviations are spelled out, returned as an
  ## attribute so it travels with the table and is printed by summary.tsa_cor().
  attr(summary_df, "abbreviations") <- c(
    "RE"          = "random effects",
    "RIS"         = "Required Information Size",
    "DARIS"       = "Diversity-Adjusted RIS",
    "AR endpoint" = "analysis-route endpoint",
    "c"           = "variance factor of Fisher's z (Var(z) = c/(n-3); 1 for Pearson r)",
    "hist. rate"  = "historical participant/information rate",
    "IPP"         = "information per participant (pooled: total information / total participants; otherwise the per-study statistic named, e.g. median)",
    "proj."       = "projection",
    "inv-var"     = "inverse-variance",
    "est."        = "estimated",
    "Add'l"       = "Additional")

  out <- list(
    data = data,
    call = match.call(),
    parameters = list(alpha_two_sided = alpha_two_sided, power = power,
                       target_r = target_r, r_anticipated = r_anticipated,
                       z_anticipated = z_anticipated,
                       cor_type = cor_type, se_source = se_source,
                       spearman_variance = spearman_variance,
                       var_factor = var_factor, ci_level = ci_level,
                       method = method,
                       method_requested = method_requested,
                       re_inference = re_inference,
                       re_inference_requested = re_inference_requested),
    res_re = res_re,
    res_fe = res_fe,
    pooled = list(z = pooled_z, r = pooled_r, r_lb = pooled_lb,
                  r_ub = pooled_ub, pval = as.numeric(res_re$pval)),
    heterogeneity = list(Q = Q, df = df, I2 = I2, tau2 = tau2, D2 = D2,
                          D2_raw = D2_raw, D2_was_capped = D2_was_capped, AF = AF),
    beta_engine = beta_engine,
    information_size = list(z_alpha = z_alpha, z_beta = z_beta,
                             info_required = info_required,
                             RIS_participants = RIS_participants,
                             DARIS_info = DARIS_info,
                             DARIS_participants = DARIS_participants,
                             DARIS_info_threshold_n = DARIS_info_threshold_n,
                             route_endpoint_info = route_endpoint_info,
                             route_endpoint_n = route_endpoint_n_est,
                             route_endpoint_participants_theoretical =
                               route_endpoint_participants_theoretical,
                             var_factor = var_factor,
                             circularity_warning = circularity_warning,
                             circularity_severe = circularity_severe),
    cumulative = cumul_df,
    boundary_timeline = boundary_timeline,
    results = list(crossed_conventional = crossed_conventional,
                   crossed_tsa = crossed_tsa,
                   entered_futility_region = entered_futility_region,
                   final_crossed_efficacy = final_crossed_efficacy,
                   final_non_efficacy = final_non_efficacy,
                   final_entered_futility_region = final_entered_futility_region,
                   final_reached = final_reached,
                   daris_reached = daris_reached,
                   final_tsa_look = final_tsa_look,
                   participants_accrued = participants_accrued,
                   info_accrued_final = info_accrued_final),
    settings = list(boundary_route = boundary_route,
                    legacy_fallback = legacy_fallback,
                    route_endpoint = route_endpoint,
                    route_endpoint_info = route_endpoint_info,
                    route_used = route_used,
                    fallback_used = fallback_used,
                    fallback_route = fallback_route,
                    fallback_reason = fallback_reason,
                    used_legacy_engine = used_legacy),
    projection = projection,
    summary_table = summary_df
  )
  class(out) <- "tsa_cor"
  out
}
