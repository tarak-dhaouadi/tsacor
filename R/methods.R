#' Print a tsa_cor object
#'
#' @param x An object of class \code{"tsa_cor"}.
#' @param ... Currently unused.
#' @return An object of class \code{"tsa_cor"}, namely the same object
#'   \code{x} supplied to the method, returned invisibly. The method
#'   prints a concise summary of the Trial Sequential Analysis results
#'   and returns the original object unchanged.
#' @export
print.tsa_cor <- function(x, ...) {
  cor_label <- .tsacor_cor_label(x$parameters$cor_type)
  cat(sprintf("Trial Sequential Analysis (%s correlations, Fisher z)\n", cor_label))
  cat("----------------------------------------------------------\n")
  cat(sprintf("Studies: %d | Participants accrued: %.0f\n",
              nrow(x$data), x$results$participants_accrued))
  cat(sprintf("Pooled %s (random effects): %.3f [95%% CI: %.3f, %.3f]\n",
              cor_label, x$pooled$r, x$pooled$r_lb, x$pooled$r_ub))
  re_inf <- x$parameters$re_inference
  if (!is.null(re_inf) && !identical(re_inf, "standard")) {
    cat(sprintf("Random-effects inference: %s\n", .tsacor_re_inference_label(re_inf)))
  }
  cat(sprintf("Anticipated %s (RIS calc): %.3f\n", cor_label, x$parameters$r_anticipated))
  cat(sprintf("Theoretical DARIS participant-equivalent: %.0f\n",
              ceiling(x$information_size$DARIS_participants)))
  analysis_route <- identical(x$settings$route_used, "analysis")
  daris_reached  <- if (is.null(x$results$daris_reached)) x$results$final_reached else x$results$daris_reached
  cat(sprintf("Crossed TSA boundary: %s | Entered futility region: %s | DARIS information reached: %s\n",
              ifelse(x$results$crossed_tsa, "YES", "NO"),
              ifelse(x$results$entered_futility_region, "YES", "NO"),
              ifelse(daris_reached, "YES", "NO")))
  if (analysis_route) {
    cat(sprintf("Boundary route: RTSA analysis | Analysis-route endpoint (%.3f x DARIS) reached: %s\n",
                x$settings$route_endpoint,
                ifelse(x$results$final_reached, "YES", "NO")))
  }
  if (isTRUE(x$results$final_reached)) {
    cat(sprintf("Definitive look crossed the efficacy boundary: %s\n",
                ifelse(isTRUE(x$results$final_crossed_efficacy), "YES", "NO")))
  }
  if (identical(x$settings$fallback_route, "design")) {
    cat("\n*** NOTE: boundary_route = \"analysis\" FAILED; the results shown are the\n")
    cat("    DESIGN-route (RTSA-derived) result -- see settings$fallback_reason. ***\n")
  }
  if (identical(x$beta_engine$engine, "legacy_r_fallback")) {
    cat("\n*** WARNING: boundaries computed with the LEGACY, APPROXIMATE fallback\n")
    cat("    engine (the RTSA-derived engine failed) -- NOT comparable with RTSA. ***\n")
  }
  cat("\nUse summary() for the full results table, or plot() for the TSA chart.\n")
  invisible(x)
}

#' Summarise a tsa_cor object
#'
#' Prints \code{object$summary_table} (\code{Parameter}/\code{Value}). The
#' \code{Parameter} column uses short abbreviations (e.g. \code{"DARIS"},
#' \code{"AR endpoint"}, \code{"RE"}) to keep the table readable, and an
#' \dQuote{Abbreviations:} line spelling them out is printed directly below
#' the table (from \code{attr(object$summary_table, "abbreviations")}).
#'
#' @param object An object of class \code{"tsa_cor"}.
#' @param ... Currently unused.
#' @return The underlying summary data.frame (invisibly printed).
#' @export
summary.tsa_cor <- function(object, ...) {
  print(object$summary_table, row.names = FALSE)
  ## The Parameter column uses short abbreviations to keep the table
  ## readable; spell them out once here (see ?tsa_cor, "Value").
  abbr <- attr(object$summary_table, "abbreviations")
  if (!is.null(abbr) && length(abbr)) {
    cat("\nAbbreviations: ",
        paste(sprintf("%s = %s", names(abbr), abbr), collapse = "; "),
        ".\n", sep = "")
  }
  if (identical(object$settings$fallback_route, "design")) {
    cat("\n*** NOTE: boundary_route = \"analysis\" FAILED; the results shown are the\n")
    cat("    DESIGN-route (RTSA-derived) result -- see settings$fallback_reason. ***\n")
  }
  if (identical(object$beta_engine$engine, "legacy_r_fallback")) {
    cat("\n*** WARNING: boundaries computed with the LEGACY, APPROXIMATE fallback\n")
    cat("    engine (the RTSA-derived engine failed) -- NOT comparable with RTSA. ***\n")
  }
  if (isTRUE(object$information_size$circularity_warning)) {
    cat("\nNOTE: target_r was not specified, so the observed pooled correlation was used\n")
    cat("for the required information size. This is circular -- see ?tsa_cor.\n")
    if (isTRUE(object$information_size$circularity_severe)) {
      cat("Accrued participants also greatly exceed the resulting DARIS, so the TSA\n")
      cat("boundary will collapse to the conventional boundary almost immediately.\n")
    }
  }
  invisible(object$summary_table)
}

## Internal helper: the second subtitle line of the TSA plot -- pooled
## random-effects correlation (back-transformed from Fisher's z) with its 95%
## CI, the p-value of the pooled effect, tau^2 (Fisher z scale) and I^2.
## Everything comes from the fitted random-effects model (x$res_re, the
## metafor::rma() object) and x$heterogeneity, i.e. the same numbers
## print()/summary() report; nothing is recomputed here. "2" is written as the
## Unicode superscript two (\u00b2). Missing pieces print as "NA" instead of
## erroring. Not exported.
.tsacor_pooled_subtitle <- function(x) {
  num <- function(v) if (is.null(v) || length(v) != 1L) NA_real_ else as.numeric(v)
  re   <- x$res_re
  b    <- num(re$b)
  lo   <- num(re$ci.lb)
  hi   <- num(re$ci.ub)
  pv   <- num(re$pval)
  tau2 <- num(x$heterogeneity$tau2)
  I2   <- num(x$heterogeneity$I2)
  sym  <- if (identical(x$parameters$cor_type, "spearman")) "\u03c1" else "r"
  p_txt <- if (is.na(pv)) "p = NA" else if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)
  sprintf("Pooled %s = %.2f [95%% CI: %.2f, %.2f] | %s | Tau\u00b2 = %.4f | I\u00b2 = %.1f%%",
          sym, tanh(b), tanh(lo), tanh(hi), p_txt, tau2, I2)
}

#' Plot a tsa_cor object
#'
#' Produces the standard Trial Sequential Analysis chart for a meta-analysis
#' of correlations: cumulative Z-curve (against the cumulative number of
#' participants), O'Brien-Fleming-type alpha (efficacy) and beta (futility)
#' spending boundaries, the conventional (naive) significance boundary, the
#' theoretical Diversity-Adjusted Required Information Size (DARIS)
#' participant-equivalent reference line, and -- whenever DARIS has actually
#' been reached (see \code{?tsa_cor}, "DARIS reached" criterion) -- a second
#' reference line marking the estimated cumulative-participants point at
#' which the observed accrued statistical information reached DARIS. The
#' two lines are shown and labelled separately, since they are different
#' quantities that need not coincide (see \code{?tsa_cor}).
#'
#' The subtitle has two lines: the model and design summary (random-effects
#' model, Diversity D^2, anticipated correlation, alpha and power), and the
#' pooled random-effects correlation with its 95\% CI (back-transformed from
#' Fisher's z), the p-value of the pooled effect, tau^2 (Fisher z scale) and
#' I^2 (the same values \code{print()} and \code{summary()} report).
#'
#' @param x An object of class \code{"tsa_cor"}.
#' @param legend Logical; show the boundary-type legend at the bottom of
#'   the plot. Default \code{TRUE}.
#' @param caption Logical; show the methods caption below the plot.
#'   Default \code{TRUE}. When \code{tsa_cor()} was run with a non-standard
#'   \code{re_inference} (\code{"hksj"}/\code{"knha"} or \code{"hksj_adhoc"}),
#'   the caption gains a line right under the first "Methods" line naming
#'   that inference option; it is absent for the default \code{"standard"}.
#'   When the route's own endpoint (DARIS for
#'   \code{boundary_route = "design"}; the analysis-route endpoint for
#'   \code{"analysis"}) has not yet been reached, the caption gains final
#'   lines with the projection: "Theoretical additional participants to
#'   DARIS", "Estimated additional participants to DARIS (historical rate)"
#'   and "Estimated additional studies required" (for the analysis route,
#'   "DARIS" reads "analysis-route endpoint") -- see \code{?tsa_cor},
#'   "Estimated additional studies/participants".
#' @param caption_size Font size for the methods caption text. Default
#'   \code{8}.
#' @param caption_face Font face for the methods caption text: one of
#'   \code{"italic"} (default, matching the previous fixed styling) or
#'   \code{"plain"}. Also accepts any other value \code{ggplot2::element_text()}
#'   understands for \code{face} (e.g. \code{"bold"}, \code{"bold.italic"}).
#' @param show_theoretical_daris Logical; show the theoretical DARIS
#'   participant-equivalent reference line and its label (see Details). Default
#'   \code{TRUE}. Set to \code{FALSE} to hide it -- e.g. when it would
#'   clutter the plot, or when only the observed-information "DARIS
#'   information reached" marker is of interest. Has no effect on the
#'   underlying DARIS calculation or on the "DARIS reached" verdict,
#'   only on what is drawn.
#' @param daris_label_size Font size for the theoretical "DARIS
#'   participant-equivalent" label. Default \code{3.2}.
#' @param daris_label_x,daris_label_y Position (in data coordinates: x =
#'   cumulative participants, y = Z-score) for the theoretical DARIS
#'   participant-equivalent label. Default \code{NULL} uses the built-in
#'   position (just right of its vertical line); it is near the top of the
#'   plot when the final Z-score is negative (or zero) and near the bottom
#'   when it is positive. The same default rule applies to the other three
#'   DARIS-related labels (historical-rate, \dQuote{DARIS information
#'   reached} and analysis-route endpoint), which keep their relative
#'   stacking order.
#' @param info_threshold_label_size Font size for the "DARIS information
#'   reached" label (only shown when DARIS has actually been reached).
#'   Default \code{3.2}.
#' @param info_threshold_label_x,info_threshold_label_y Position (in data
#'   coordinates) for the "DARIS information reached" label. Default
#'   \code{NULL} uses the built-in position.
#' @param participants_label_size Font size for the "Participants accrued"
#'   label. Default \code{3.2}.
#' @param participants_label_x,participants_label_y Position (in data
#'   coordinates) for the "Participants accrued" label. Default \code{NULL}
#'   uses the built-in position (right-aligned at the last data point): near
#'   the bottom of the plot when the final Z-score is negative (or zero), and
#'   near the top, over the curve, when it is positive.
#' @param endpoint_label_size Font size for the "Analysis-route endpoint
#'   (Design_R x DARIS) reached" label (only shown when
#'   \code{boundary_route = "analysis"}; it is also drawn, worded
#'   "not yet reached; theoretical ~ N participants", at the theoretical
#'   position when the endpoint has not been reached). Default \code{NULL}
#'   uses the same size as \code{info_threshold_label_size} (\code{3.2}
#'   unless changed).
#' @param endpoint_label_x,endpoint_label_y Position (in data coordinates: x
#'   = cumulative participants, y = Z-score) for the "Analysis-route endpoint
#'   (Design_R x DARIS) reached" label. Default \code{NULL} uses the
#'   built-in position (just right of its vertical line, below the "DARIS
#'   information reached" label).
#' @param show_historical_daris Logical; when the route's own target has not
#'   been reached, draw an extra reference line at the target position
#'   projected from the historical information-per-participant rate (accrued
#'   participants plus the estimated additional participants,
#'   \code{projection$target_participants_historical_rate}). For
#'   \code{boundary_route = "design"} the target is DARIS and the line
#'   ("DARIS (historical rate)") sits next to the theoretical
#'   DARIS line; for \code{"analysis"} the target is the analysis-route
#'   endpoint (design_R x DARIS) and the line ("Historical
#'   information/participant-rate projection") sits next to the theoretical
#'   endpoint line; it is a projection of where that same endpoint would be
#'   reached, not a second definition of the endpoint.
#'   Default \code{TRUE}. Only affects what is drawn.
#' @param historical_label_size Font size for the historical-rate label.
#'   Default \code{NULL} uses \code{daris_label_size}.
#' @param historical_label_x,historical_label_y Position (in data
#'   coordinates) for the historical-rate label. Default \code{NULL} uses
#'   the built-in position (just right of its vertical line; design route
#'   below the theoretical DARIS label, analysis route below the
#'   analysis-route endpoint label).
#' @param xmax_mult Positive number; multiplier applied to the largest x
#'   value that must fit in the plot (accrued participants, theoretical DARIS,
#'   DARIS information marker, historical-rate projection, analysis-route
#'   endpoint, and the last x of the formal boundaries) to obtain the upper
#'   limit of the x-axis. Default
#'   \code{1.15}, i.e. 15\% of free space to the right. Use a larger value
#'   (e.g. \code{1.4}) to leave more room for labels, or \code{1} to end the
#'   axis exactly at the largest element. Values below \code{1} crop the
#'   right-hand part of the plot.
#' @param alpha_col Color for the alpha (efficacy) boundary line. Default
#'   \code{"firebrick"}.
#' @param beta_col Color for the beta (futility) boundary line. Default
#'   \code{"blue"}.
#' @param naive_col Color for the naive/conventional significance boundary
#'   line. Default \code{"darkgreen"}.
#' @param z_col Color for the cumulative Z-score line/points. Default
#'   \code{"black"}.
#' @param ... Currently unused.
#' @return A \code{ggplot} object (invisibly), also drawn on the current
#'   graphics device / returned for further customisation, e.g.
#'   \code{ggplot2::ggsave()}.
#' @export
plot.tsa_cor <- function(x, legend = TRUE, caption = TRUE,
                         caption_size = 8, caption_face = "italic",
                         show_theoretical_daris = TRUE,
                         daris_label_size = 3.2,
                         daris_label_x = NULL, daris_label_y = NULL,
                         info_threshold_label_size = 3.2,
                         info_threshold_label_x = NULL, info_threshold_label_y = NULL,
                         participants_label_size = 3.2,
                         participants_label_x = NULL, participants_label_y = NULL,
                         endpoint_label_size = NULL,
                         endpoint_label_x = NULL, endpoint_label_y = NULL,
                         show_historical_daris = TRUE,
                         historical_label_size = NULL,
                         historical_label_x = NULL, historical_label_y = NULL,
                         xmax_mult = 1.15,
                         alpha_col = "firebrick", beta_col = "blue",
                         naive_col = "darkgreen", z_col = "black",
                         ...) {

  cumul_df <- x$cumulative

  if (!is.numeric(xmax_mult) || length(xmax_mult) != 1L ||
      !is.finite(xmax_mult) || xmax_mult <= 0) {
    stop("`xmax_mult` must be a single positive number (default 1.15).",
         call. = FALSE)
  }
  DARIS_participants <- x$information_size$DARIS_participants
  DARIS_info_threshold_n <- x$information_size$DARIS_info_threshold_n
  final_reached <- x$results$final_reached
  ## DARIS (t = 1) and the route endpoint are separate quantities.
  ## For boundary_route = "design" they coincide (route endpoint == DARIS)
  ## and the plot is unchanged; for "analysis" the formal endpoint is
  ## design_R * DARIS and gets its own, separately labelled marker.
  daris_reached <- if (is.null(x$results$daris_reached)) final_reached else x$results$daris_reached
  analysis_route <- identical(x$settings$route_used, "analysis")
  route_endpoint <- if (is.null(x$settings$route_endpoint)) 1 else x$settings$route_endpoint
  route_endpoint_n <- x$information_size$route_endpoint_n
  show_endpoint_marker <- analysis_route && isTRUE(final_reached) &&
    !is.null(route_endpoint_n) && !is.na(route_endpoint_n)
  ## The formal boundaries always terminate at the route endpoint
  ## (observed-information estimate if reached, otherwise the theoretical
  ## participant-equivalent c * route_endpoint * DARIS + 3). When design_R > 1
  ## the endpoint typically lies beyond the observed information, so it is NOT
  ## reached; the theoretical position is then drawn and labelled "not yet
  ## reached".
  endpoint_theoretical_participants <- x$information_size$route_endpoint_participants_theoretical
  if (is.null(endpoint_theoretical_participants)) endpoint_theoretical_participants <- NA_real_
  show_endpoint_theoretical <- analysis_route && !show_endpoint_marker &&
    is.finite(endpoint_theoretical_participants) &&
    !(isTRUE(all.equal(route_endpoint, 1)) && isTRUE(show_theoretical_daris))

  ## Two distinct quantities are shown on the plot, deliberately NOT
  ## conflated into a single "DARIS participants" figure:
  ##
  ##  1. DARIS_participants: the THEORETICAL participant-equivalent of the
  ##     required information (c * DARIS + 3, the single-study sample-size
  ##     equivalent on the Fisher z scale). Always drawn.
  ##  2. DARIS_info_threshold_n: an ESTIMATE (interpolated between
  ##     looks) of the cumulative-participants point at which the OBSERVED
  ##     accrued inverse-variance information actually reached DARIS_info
  ##     -- only meaningful, and only drawn, when DARIS has been reached
  ##     (see tsa_cor(), Section 6b). This is why the printed "DARIS
  ##     reached" verdict can never contradict what is plotted: whenever
  ##     it says YES, this second marker is shown at (or before) the
  ##     accrued-participants point; whenever it says NO, only the theoretical
  ##     line (1) is shown, clearly labelled as not yet reached.
  show_info_threshold_marker <- daris_reached && !is.na(DARIS_info_threshold_n)

  ## (route's own target not reached): an extra reference line at the target
  ## position projected from the historical information-per-participant rate,
  ## i.e. accrued participants + estimated additional participants. Design
  ## route: DARIS; analysis route: the analysis-route endpoint (design_R x
  ## DARIS). Distinct from the theoretical line (1) and the
  ## observed-information markers.
  historical_target_participants <- if (is.null(x$projection)) NA_real_ else
    x$projection$target_participants_historical_rate
  if (is.null(historical_target_participants) || length(historical_target_participants) != 1L)
    historical_target_participants <- NA_real_
  historical_target_reached <- if (analysis_route) isTRUE(final_reached) else isTRUE(daris_reached)
  show_historical_daris <- isTRUE(show_historical_daris) && !historical_target_reached &&
    is.finite(historical_target_participants)

  z_alpha <- x$information_size$z_alpha
  D2 <- x$heterogeneity$D2
  AF <- x$heterogeneity$AF
  alpha_two_sided <- x$parameters$alpha_two_sided
  power <- x$parameters$power
  r_anticipated <- x$parameters$r_anticipated
  cor_type <- if (is.null(x$parameters$cor_type)) "pearson" else x$parameters$cor_type
  cor_label <- .tsacor_cor_label(cor_type)
  cor_sym <- if (identical(cor_type, "spearman")) "\u03c1" else "r"

  ## Boundary plotting uses the dedicated RTSA-style timeline.  The
  ## observed cumulative table continues through all studies, whereas the
  ## formal boundaries terminate at the t = 1 information target.  If DARIS
  ## was reached in the observed data, that endpoint is plotted at the
  ## interpolated observed-information participant coordinate
  ## (DARIS_info_threshold_n); the theoretical DARIS_participants line remains
  ## separate.
  boundary_line <- x$boundary_timeline[, c("cum_n", "TSA_boundary_upper",
                                           "TSA_boundary_lower",
                                           "TSA_futility_upper",
                                           "TSA_futility_lower")]

  finite_bounds <- c(boundary_line$TSA_boundary_upper[is.finite(boundary_line$TSA_boundary_upper)],
                      abs(boundary_line$TSA_boundary_lower[is.finite(boundary_line$TSA_boundary_lower)]))
  y_abs_max <- max(abs(cumul_df$Z), finite_bounds, na.rm = TRUE)
  y_limit <- y_abs_max * 1.15

  ## Default vertical placement of the labels depends on where the Z-curve
  ## ends. The four DARIS-related labels (theoretical DARIS, historical-rate
  ## projection, "DARIS information reached", analysis-route endpoint) sit in
  ## the UPPER part of the plot when the Z-curve is negative (or zero) and in
  ## the LOWER part when it is positive, i.e. always on the side away from the
  ## curve; the "Participants accrued" label takes the opposite side (over
  ## the curve when it is positive). `daris_sign` = +1 (upper) / -1 (lower);
  ## the user-supplied *_label_y arguments always override these defaults.
  z_last <- {
    zz <- cumul_df$Z[!is.na(cumul_df$Z)]
    if (length(zz)) zz[length(zz)] else NA_real_
  }
  z_positive <- is.finite(z_last) && z_last > 0
  daris_sign <- if (z_positive) -1 else 1

  boundary_line$TSA_boundary_upper <- pmin(boundary_line$TSA_boundary_upper, y_limit)
  boundary_line$TSA_boundary_lower <- pmax(boundary_line$TSA_boundary_lower, -y_limit)

  participants_accrued <- x$results$participants_accrued
  x_max <- max(c(cumul_df$cum_n,
                  if (show_theoretical_daris) DARIS_participants else NA,
                  if (show_info_threshold_marker) DARIS_info_threshold_n else NA,
                  if (show_endpoint_marker) route_endpoint_n else NA,
                  if (show_endpoint_theoretical) endpoint_theoretical_participants else NA,
                  if (show_historical_daris) historical_target_participants else NA,
                  ## the boundaries' own last x (the route endpoint)
                  ## must always fit inside the axis range
                  boundary_line$cum_n[is.finite(boundary_line$cum_n)]),
               na.rm = TRUE) * xmax_mult

  alpha_lines <- data.frame(
    cum_n = rep(boundary_line$cum_n, 2),
    y = c(boundary_line$TSA_boundary_upper, boundary_line$TSA_boundary_lower),
    side = rep(c("upper", "lower"), each = nrow(boundary_line)),
    type = "Alpha boundaries"
  )
  beta_lines <- data.frame(
    cum_n = rep(boundary_line$cum_n, 2),
    y = c(boundary_line$TSA_futility_upper, boundary_line$TSA_futility_lower),
    side = rep(c("upper", "lower"), each = nrow(boundary_line)),
    type = "Non-binding futility boundaries"
  )
  naive_lines <- data.frame(
    cum_n = rep(c(0, x_max), 2),
    y = rep(c(z_alpha, -z_alpha), each = 2),
    side = rep(c("upper", "lower"), each = 2),
    type = "Naive boundaries"
  )
  ## the first-look Z is NA under HKSJ re_inference (undefined at
  ## k = 1; see .tsacor_cumulative_re_inference()). That row is dropped here
  ## explicitly, rather than relying only on the layers' na.rm = TRUE, so
  ## the point/line are never drawn regardless of ggplot2 version behaviour.
  z_curve_df <- cumul_df[!is.na(cumul_df$Z), , drop = FALSE]
  z_line <- data.frame(
    cum_n = z_curve_df$cum_n,
    y = z_curve_df$Z,
    side = "z",
    type = "Z scores"
  )

  plot_lines <- rbind(alpha_lines, beta_lines, naive_lines, z_line)
  line_types  <- c("Alpha boundaries" = "solid", "Non-binding futility boundaries" = "dashed",
                    "Naive boundaries" = "dashed", "Z scores" = "solid")
  line_colors <- c("Alpha boundaries" = alpha_col, "Non-binding futility boundaries" = beta_col,
                    "Naive boundaries" = naive_col, "Z scores" = z_col)
  line_widths <- c("Alpha boundaries" = 0.8, "Non-binding futility boundaries" = 0.8,
                    "Naive boundaries" = 0.5, "Z scores" = 0.8)

  cum_n <- y <- type <- side <- Z <- NULL  # avoid R CMD check NOTE for NSE

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = plot_lines,
      ggplot2::aes(x = cum_n, y = y, color = type, linetype = type,
                   linewidth = type, group = interaction(type, side)),
      na.rm = TRUE
    ) +
    ## z_curve_df already excludes the NA first look under HKSJ re_inference
    ## (see above); na.rm = TRUE is kept as a belt-and-braces guard against
    ## any other NA that might reach this layer.
    ggplot2::geom_point(data = z_curve_df, ggplot2::aes(x = cum_n, y = Z, color = "Z scores"),
                         size = 2, na.rm = TRUE) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = x_max, y = 0, yend = 0),
                           color = "grey60", linewidth = 0.3) +
    ggplot2::annotate("text",
                       x = if (is.null(participants_label_x)) max(cumul_df$cum_n) else participants_label_x,
                       y = if (is.null(participants_label_y)) -daris_sign * y_limit * 0.92 else participants_label_y,
                       label = paste0("Participants accrued = ", participants_accrued),
                       hjust = 1, vjust = 0, size = participants_label_size, color = "steelblue4") +
    ggplot2::scale_color_manual(name = NULL, values = line_colors) +
    ggplot2::scale_linetype_manual(name = NULL, values = line_types) +
    ggplot2::scale_linewidth_manual(values = line_widths, guide = "none") +
    ggplot2::labs(
      title = sprintf("Trial Sequential Analysis of Correlations (%s, Fisher z)", cor_label),
      subtitle = paste0(
        sprintf(
          "Random-effects model | Diversity D\u00b2 = %.0f%% | Anticipated %s = %.2f | alpha=%.0f%%, power=%.0f%%",
          D2 * 100, cor_sym, r_anticipated, alpha_two_sided * 100, power * 100),
        "\n", .tsacor_pooled_subtitle(x)),
      x = "Cumulative number of participants",
      y = "Cumulative Z-score"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::coord_cartesian(xlim = c(0, x_max), ylim = c(-y_limit, y_limit), clip = "off")

  ## Theoretical DARIS participant-equivalent reference line (Section 4 of
  ## tsa_cor()). Optional -- toggled off via show_theoretical_daris = FALSE
  ## -- since some users only want the observed-information marker below,
  ## or find two reference lines cluttered. Hiding it only affects the
  ## plot; the underlying DARIS calculation and "DARIS reached" verdict
  ## are unchanged either way.
  if (show_theoretical_daris) {
    p <- p +
      ggplot2::geom_vline(xintercept = DARIS_participants, color = "black",
                           linetype = "dotted", linewidth = 0.6) +
      ggplot2::annotate("text",
                         x = if (is.null(daris_label_x)) DARIS_participants else daris_label_x,
                         y = if (is.null(daris_label_y)) daris_sign * y_limit * 0.92 else daris_label_y,
                         label = paste0("Theoretical DARIS participant-equivalent ~ ", ceiling(DARIS_participants)),
                         hjust = -0.05, vjust = 0, size = daris_label_size)
  }

  ## Route's own target not reached: target position at the historical
  ## information-per-participant rate (accrued participants + estimated
  ## additional participants). Drawn in addition to -- never instead of -- the theoretical
  ## line(s); the label says it is projected. Design route: DARIS (label y =
  ## 0.75, below the theoretical DARIS label; the observed-information
  ## marker that shares that height needs DARIS to have been reached, so
  ## they never coexist). Analysis route: the analysis-route endpoint
  ## (label y = 0.41, below the endpoint label at 0.58 and the "DARIS
  ## information reached" label at 0.75, which can coexist with it).
  if (show_historical_daris) {
    hist_label <- if (analysis_route) {
      ## deliberately NOT worded "analysis-route endpoint ...", which
      ## would read as a second, competing definition of the endpoint; the
      ## theoretical endpoint line keeps that name, this one is a projection.
      paste0("Historical information/participant-rate projection ~ ",
             ceiling(historical_target_participants), " participants")
    } else {
      paste0("DARIS (historical rate) ~ ",
             ceiling(historical_target_participants), " participants (projected)")
    }
    p <- p +
      ggplot2::geom_vline(xintercept = historical_target_participants, color = "darkorange3",
                           linetype = "dotdash", linewidth = 0.6) +
      ggplot2::annotate("text",
                         x = if (is.null(historical_label_x)) historical_target_participants
                             else historical_label_x,
                         y = if (is.null(historical_label_y)) {
                               daris_sign * y_limit * (if (analysis_route) 0.41 else 0.75)
                             } else historical_label_y,
                         label = hist_label,
                         hjust = -0.05, vjust = 0,
                         size = if (is.null(historical_label_size)) daris_label_size
                                else historical_label_size,
                         color = "darkorange3")
  }

  ## Second reference marker: the ESTIMATED cumulative-participants point at
  ## which the observed accrued statistical information reached DARIS_info
  ## (see tsa_cor(), Section 7b). Added separately, and only when DARIS has
  ## actually been reached, so it is never confused with the theoretical
  ## participant-equivalent line above -- both are shown, distinctly labelled,
  ## per the package's documented "observed inverse-variance information"
  ## criterion.
  if (show_info_threshold_marker) {
    p <- p +
      ggplot2::geom_vline(xintercept = DARIS_info_threshold_n, color = "grey35",
                           linetype = "dashed", linewidth = 0.6) +
      ggplot2::annotate("text",
                         x = if (is.null(info_threshold_label_x)) DARIS_info_threshold_n
                             else info_threshold_label_x,
                         y = if (is.null(info_threshold_label_y)) daris_sign * y_limit * 0.75
                             else info_threshold_label_y,
                         label = paste0("DARIS information reached ~ ",
                                        ceiling(DARIS_info_threshold_n), " participants (est.)"),
                         hjust = -0.05, vjust = 0, size = info_threshold_label_size,
                         color = "grey35")
  }

  ## Analysis route only: the formal analysis endpoint (design_R x DARIS
  ## information), estimated by interpolation between looks. Drawn
  ## separately from the DARIS marker above so the two are never conflated.
  if (show_endpoint_marker) {
    p <- p +
      ggplot2::geom_vline(xintercept = route_endpoint_n, color = "purple4",
                           linetype = "longdash", linewidth = 0.6) +
      ggplot2::annotate("text",
                         x = if (is.null(endpoint_label_x)) route_endpoint_n
                             else endpoint_label_x,
                         y = if (is.null(endpoint_label_y)) daris_sign * y_limit * 0.58
                             else endpoint_label_y,
                         label = paste0("Analysis-route endpoint (",
                                        sprintf("%.3f", route_endpoint),
                                        " x DARIS) reached ~ ",
                                        ceiling(route_endpoint_n), " participants (est.)"),
                         hjust = -0.05, vjust = 0,
                         size = if (is.null(endpoint_label_size)) info_threshold_label_size
                                else endpoint_label_size,
                         color = "purple4")
  }

  ## Analysis route, endpoint NOT yet reached in the observed data (typically
  ## design_R > 1): mark where the formal boundaries terminate, i.e. the
  ## theoretical participant-equivalent of design_R x DARIS. Analogous to the
  ## theoretical DARIS line; the label says explicitly that it is not reached.
  if (show_endpoint_theoretical) {
    p <- p +
      ggplot2::geom_vline(xintercept = endpoint_theoretical_participants, color = "purple4",
                           linetype = "longdash", linewidth = 0.6) +
      ggplot2::annotate("text",
                         x = if (is.null(endpoint_label_x)) endpoint_theoretical_participants
                             else endpoint_label_x,
                         y = if (is.null(endpoint_label_y)) daris_sign * y_limit * 0.58
                             else endpoint_label_y,
                         label = paste0("Analysis-route endpoint (",
                                        sprintf("%.3f", route_endpoint),
                                        " x DARIS) not yet reached; theoretical ~ ",
                                        ceiling(endpoint_theoretical_participants), " participants"),
                         hjust = -0.05, vjust = 0,
                         size = if (is.null(endpoint_label_size)) info_threshold_label_size
                                else endpoint_label_size,
                         color = "purple4")
  }

  if (caption) {
    se_source <- x$parameters$se_source
    se_txt <- if (identical(se_source, "n")) {
      if (identical(cor_type, "spearman")) {
        sprintf("sample sizes (%s variance)",
                if (identical(x$parameters$spearman_variance, "bonett_wright"))
                  "Bonett-Wright" else "Fieller")
      } else {
        "sample sizes (1/(n-3))"
      }
    } else {
      "reported confidence intervals"
    }
    methods_caption <- sprintf(
      paste0("Methods: Random-effects (%s) model of Fisher's z (%s); SE(z) from %s\n",
             "Alpha spending: O'Brien-Fleming-type (asOF); ",
             "Non-binding futility: RTSA-reconstructed recursive-integration ",
             "engine, O'Brien-Fleming-type beta-spending (bsOF)\n",
             "alpha = %.0f%% (two-sided), power = %.0f%% | Diversity D\u00b2 = %.0f%%, Adjustment factor = %.2f"),
      .tsacor_method_label(if (is.null(x$parameters$method)) "DL" else x$parameters$method),
      cor_label, se_txt, alpha_two_sided * 100, power * 100, D2 * 100, AF)
    ## Name the random-effects inference option on its own line, directly
    ## under the first "Methods:" line, when it is not "standard".
    re_inf <- x$parameters$re_inference
    if (!is.null(re_inf) && !identical(re_inf, "standard")) {
      re_line <- if (identical(re_inf, "hksj_adhoc")) {
        paste0("Random-effects inference: HKSJ with ad hoc correction (variance scale max(1, q); ",
               "t distribution, k-1 df); Z = normal-equivalent of the t-statistic")
      } else {
        paste0("Random-effects inference: HKSJ (Hartung-Knapp-Sidik-Jonkman; ",
               "t distribution, k-1 df); Z = normal-equivalent of the t-statistic")
      }
      cap_split <- strsplit(methods_caption, "\n", fixed = TRUE)[[1]]
      methods_caption <- paste(c(cap_split[1], re_line, cap_split[-1]), collapse = "\n")
    }
    if (analysis_route) {
      methods_caption <- paste0(methods_caption, sprintf(
        "\nBoundary route: RTSA analysis (formal endpoint = %.3f x DARIS information)",
        route_endpoint))
    }
    ## Estimated additional participants/studies (6d in tsa_cor()) -- appended
    ## as their own caption lines, below everything else, only when the
    ## route's own target has not been reached (design: DARIS; analysis: the
    ## analysis-route endpoint): the theoretical participant figure AND the
    ## historical-rate estimate (remaining information / information per
    ## participant, NOT whole studies x participants per study), plus the
    ## studies estimate.
    projection <- x$projection
    if (!is.null(projection)) {
      show_projection_caption <- if (analysis_route) !isTRUE(final_reached) else !isTRUE(daris_reached)
      fmt_n <- function(v) formatC(ceiling(v), format = "d", big.mark = ",")
      n_studies_cap <- projection$n_additional_studies
      est_n_cap <- projection$additional_participants_estimated
      cap_lines <- character(0)
      if (show_projection_caption) {
        target_txt <- if (analysis_route) "analysis-route endpoint" else "DARIS"
        theo_n_cap <- projection$additional_participants_theoretical
        n_parts <- c(
          if (!is.null(theo_n_cap) && !is.na(theo_n_cap))
            sprintf("Theoretical additional participants to %s: %s",
                    target_txt, fmt_n(theo_n_cap)),
          if (!is.null(est_n_cap) && !is.na(est_n_cap))
            sprintf("Estimated additional participants to %s (historical rate): ~%s",
                    target_txt, fmt_n(est_n_cap)))
        cap_lines <- c(
          if (length(n_parts)) paste(n_parts, collapse = ", "),
          if (!is.null(n_studies_cap) && !is.na(n_studies_cap))
            sprintf("Estimated additional studies required: %d", n_studies_cap))
      }
      if (length(cap_lines)) {
        methods_caption <- paste0(methods_caption, "\n",
                                  paste(cap_lines, collapse = "\n"))
      }
    }
    p <- p + ggplot2::labs(caption = methods_caption) +
      ggplot2::theme(plot.caption = ggplot2::element_text(
        hjust = 0, size = caption_size, face = caption_face))
  }

  if (legend) {
    p <- p + ggplot2::theme(legend.position = "bottom")
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  print(p)
  invisible(p)
}
