## Functional tests for tsa_cor(), against the bundled r_meta.xlsx example
## (20 studies, Pearson r, 2005-2024). Reference values below were computed
## independently in Python (DerSimonian-Laird on Fisher's z; see PR notes) and
## are checked at a loose tolerance since they are cross-implementation
## checks, not bit-for-bit ports of tsa_cor()'s own arithmetic. The compiled
## RTSA-derived boundary engine itself is checked bit-for-bit against RTSA's
## own output in test-rtsa-live-reference.R and friends; this file is about
## the correlation-specific statistics layer built on top of it.

.path <- function() tsacor_example_data()

.fit <- function(...) {
  suppressWarnings(suppressMessages(tsa_cor(.path(), verbose = FALSE, ...)))
}

## Same as .fit() but WITHOUT suppressing warnings, for tests that assert on
## the warnings themselves (expect_warning() cannot see warnings that .fit()
## has already swallowed).
.fit_w <- function(...) {
  suppressMessages(tsa_cor(.path(), verbose = FALSE, ...))
}

test_that("tsacor_example_data() points at a real, readable file with the expected columns", {
  p <- tsacor_example_data()
  expect_true(file.exists(p))
  d <- as.data.frame(readxl::read_excel(p))
  expect_identical(nrow(d), 20L)
  expect_true(all(c("Study", "Year", "Ethnicity", "Age", "r", "lbound",
                     "ubound", "n_subjects") %in% names(d)))
  expect_true(all(d$r > d$lbound & d$r < d$ubound))
  expect_equal(sum(d$n_subjects), 5355)
})

test_that("data validation rejects malformed inputs before any fitting", {
  d <- as.data.frame(readxl::read_excel(.path()))
  expect_error(tsa_cor(d[, setdiff(names(d), "n_subjects")], target_r = 0.10, verbose = FALSE),
               "Missing required column")
  bad <- d; bad$r[1] <- 1.4
  expect_error(tsa_cor(bad, target_r = 0.10, verbose = FALSE),
               "strictly between -1 and 1")
  bad2 <- d; bad2$n_subjects[1] <- 3
  expect_error(tsa_cor(bad2, target_r = 0.10, verbose = FALSE),
               "greater than 3")
  bad3 <- d; bad3$n_subjects[1] <- 10.5
  expect_error(tsa_cor(bad3, target_r = 0.10, verbose = FALSE),
               "whole numbers")
  bad4 <- d; bad4$lbound[1] <- bad4$ubound[1] + 0.01
  expect_error(tsa_cor(bad4, target_r = 0.10, verbose = FALSE),
               "lbound must be smaller than ubound")
  bad5 <- d; bad5$Study[2] <- bad5$Study[1]
  expect_error(tsa_cor(bad5, target_r = 0.10, verbose = FALSE),
               "unique")
  bad6 <- d; bad6$r <- as.character(bad6$r)
  expect_error(tsa_cor(bad6, target_r = 0.10, verbose = FALSE),
               "must be numeric")
  expect_error(tsa_cor(d, target_r = 0, verbose = FALSE), "cannot equal 0")
  expect_error(tsa_cor(d, target_r = 1, verbose = FALSE),
               "strictly between -1 and 1")
  expect_error(tsa_cor(d[1, , drop = FALSE], target_r = 0.10, verbose = FALSE),
               "at least two studies")
  expect_error(tsa_cor(d, target_r = 0.10, method = "GENQ", verbose = FALSE),
               "not currently supported")
  expect_error(tsa_cor(d, target_r = 0.10, method = "bogus", verbose = FALSE),
               "method must be one of")
  expect_no_warning(tsa_cor(d, target_r = 0.10, verbose = FALSE, order_by = "Year"))
})

test_that("colliding column names after space-to-underscore normalisation are refused", {
  d <- as.data.frame(readxl::read_excel(.path()))
  d$`n subjects` <- d$n_subjects # a second column that collides once
  expect_error(tsa_cor(d, target_r = 0.10, verbose = FALSE),        # normalised
               "not unique after spaces")
})

test_that("Fisher z, se_z_n and se_z_ci are computed correctly and se_source selects between them", {
  d <- as.data.frame(readxl::read_excel(.path()))
  res_ci <- .fit(target_r = 0.10, se_source = "ci")
  res_n  <- .fit(target_r = 0.10, se_source = "n")

  expect_equal(res_ci$data$z_fisher, atanh(d$r), tolerance = 1e-10)
  expect_equal(res_n$data$se_z_n, 1 / sqrt(d$n_subjects - 3), tolerance = 1e-10)
  zcrit <- stats::qnorm(0.975)
  expect_equal(res_ci$data$se_z_ci,
               (atanh(d$ubound) - atanh(d$lbound)) / (2 * zcrit), tolerance = 1e-10)
  ## se_source picks the column actually used
  expect_identical(res_ci$data$se_z, res_ci$data$se_z_ci)
  expect_identical(res_n$data$se_z, res_n$data$se_z_n)
  ## both variants are always present regardless of which is used
  expect_true(all(is.finite(res_ci$data$se_z_n)))
  expect_true(all(is.finite(res_n$data$se_z_ci)))
  expect_identical(res_ci$parameters$var_factor, 1) # Pearson: c = 1
})

test_that("Spearman variance models change se_z_n (and only se_z_n) as expected", {
  d <- as.data.frame(readxl::read_excel(.path()))
  fieller <- .fit(target_r = 0.10, cor_type = "spearman",
                  se_source = "n", spearman_variance = "fieller")
  bw <- .fit(target_r = 0.10, cor_type = "spearman",
             se_source = "n", spearman_variance = "bonett_wright")
  expect_equal(fieller$data$se_z_n, sqrt(1.06 / (d$n_subjects - 3)), tolerance = 1e-10)
  expect_equal(bw$data$se_z_n, sqrt((1 + d$r^2 / 2) / (d$n_subjects - 3)), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fieller$data$se_z_n, bw$data$se_z_n)))
  ## se_z_ci is unaffected by spearman_variance (it only touches se_source="n")
  expect_equal(fieller$data$se_z_ci, bw$data$se_z_ci)
})

test_that("target_r = NA uses the observed pooled effect and warns about circularity", {
  expect_warning(res <- .fit_w(target_r = NA_real_), "circular")
  expect_true(res$information_size$circularity_warning)
  expect_equal(res$parameters$r_anticipated, res$pooled$r, tolerance = 1e-8)
})

test_that("target_r close to zero warns but is not an error", {
  expect_warning(.fit_w(target_r = 0.05), "close to the null value")
})

test_that("random-effects meta-analysis and heterogeneity match an independent DerSimonian-Laird computation", {
  res <- .fit(target_r = 0.10, method = "DL")
  expect_equal(res$pooled$r, 0.4997, tolerance = 1e-3)
  expect_equal(res$pooled$r_lb, 0.4572, tolerance = 1e-3)
  expect_equal(res$pooled$r_ub, 0.5398, tolerance = 1e-3)
  expect_equal(res$heterogeneity$Q, 72.177, tolerance = 1e-2)
  expect_equal(res$heterogeneity$I2, 73.68, tolerance = 0.2)
  expect_equal(res$heterogeneity$D2, 0.7462, tolerance = 1e-3)
  expect_equal(res$heterogeneity$AF, 3.9406, tolerance = 1e-3)
  expect_false(res$heterogeneity$D2_was_capped)
})

test_that("required information size and DARIS follow the stated formulas for a Pearson target", {
  res <- .fit(target_r = 0.10)
  z_alpha <- stats::qnorm(0.975); z_beta <- stats::qnorm(0.80)
  info_req <- (z_alpha + z_beta)^2 / atanh(0.10)^2
  expect_equal(res$information_size$info_required, info_req, tolerance = 1e-8)
  expect_equal(res$information_size$RIS_participants, info_req + 3, tolerance = 1e-8)
  DARIS <- info_req * res$heterogeneity$AF
  expect_equal(res$information_size$DARIS_info, DARIS, tolerance = 1e-6)
  expect_equal(res$information_size$DARIS_participants, DARIS + 3, tolerance = 1e-6)
  expect_equal(res$information_size$var_factor, 1)
})

test_that("the variance factor c multiplies participant-equivalents for Spearman targets", {
  res_p <- .fit(target_r = 0.10, cor_type = "pearson")
  res_s <- .fit(target_r = 0.10, cor_type = "spearman", se_source = "n",
                spearman_variance = "fieller")
  expect_equal(res_s$information_size$var_factor, 1.06)
  ## same required information (c doesn't enter info_required), different
  ## participant-equivalents (c does enter those)
  expect_equal(res_s$information_size$info_required, res_p$information_size$info_required,
               tolerance = 1e-8)
  expect_gt(res_s$information_size$RIS_participants, res_p$information_size$RIS_participants)
})

test_that("cumulative analysis is monotonic in information and reaches its DARIS target midway with target_r = 0.10", {
  res <- .fit(target_r = 0.10, order_by = "Year")
  cd <- res$cumulative
  expect_identical(nrow(cd), 20L)
  expect_true(all(diff(cd$info_accrued) > 0))
  expect_true(all(diff(cd$cum_n) > 0))
  expect_equal(cd$r_estimate, tanh(cd$estimate), tolerance = 1e-10)
  expect_equal(tail(cd$info_fraction, 1), tail(cd$info_accrued, 1) / res$information_size$DARIS_info,
               tolerance = 1e-10)
  ## DARIS is crossed partway through this 20-study series, not at the very
  ## first or very last look
  reach <- which(cd$info_fraction >= 1)[1]
  expect_true(is.finite(reach) && reach > 1 && reach < 20)
  expect_true(res$results$daris_reached)
})

test_that("boundary_route = 'design' terminates the formal boundaries at DARIS itself", {
  res <- .fit(target_r = 0.10, boundary_route = "design")
  expect_identical(res$settings$route_used, "design")
  expect_equal(res$settings$route_endpoint, 1)
  expect_equal(tail(res$boundary_timeline$info_fraction, 1), 1, tolerance = 1e-10)
  expect_true(tail(res$boundary_timeline$synthetic, 1))
})

test_that("boundary_route = 'analysis' produces a route endpoint at design_R x DARIS", {
  res <- .fit(target_r = 0.10, boundary_route = "analysis")
  expect_identical(res$settings$route_used, "analysis")
  expect_gte(res$settings$route_endpoint, 1)
  expect_equal(res$information_size$route_endpoint_info,
               res$settings$route_endpoint * res$information_size$DARIS_info,
               tolerance = 1e-6)
})

test_that("a much smaller target correlation than observed is not reached and triggers a projection", {
  ## The pooled r is about 0.50, so DARIS is reached for any target of about
  ## 0.08 or more; a small target such as 0.05 needs far more information
  ## than the 20 studies supply.
  res <- .fit(target_r = 0.05)
  expect_false(res$results$daris_reached)
  expect_false(res$results$final_reached)
  expect_true(is.na(res$results$final_crossed_efficacy))
  expect_gt(res$information_size$DARIS_participants, sum(res$data$n_subjects))
  expect_true(is.finite(res$projection$additional_participants_estimated))
  expect_true(is.finite(res$projection$additional_participants_theoretical))
  expect_true(is.finite(res$projection$n_additional_studies))
  expect_gte(res$projection$n_additional_studies, 1L)
})

test_that("re_inference is validated, case-insensitive and aliased, and changes only the inference layer", {
  expect_error(.fit(target_r = 0.10, re_inference = "bogus"), "re_inference must be one of")
  expect_identical(.fit(target_r = 0.10, re_inference = "KNHA")$parameters$re_inference, "hksj")
  expect_identical(.fit(target_r = 0.10, re_inference = "Hksj_adhoc")$parameters$re_inference,
                   "hksj_adhoc")

  rs <- .fit(target_r = 0.10, re_inference = "standard")
  rh <- .fit(target_r = 0.10, re_inference = "hksj")
  expect_identical(rh$res_re$test, "knha")
  expect_equal(as.numeric(rh$res_re$b), as.numeric(rs$res_re$b), tolerance = 1e-10)
  expect_equal(rh$heterogeneity, rs$heterogeneity)
  expect_equal(rh$information_size$DARIS_info, rs$information_size$DARIS_info)
  expect_equal(rh$boundary_timeline, rs$boundary_timeline)
  ## HKSJ is undefined at the very first cumulative look
  expect_true(is.na(rh$cumulative$Z[1]))
  expect_false(is.na(rh$cumulative$Z[2]))
})

test_that("print/summary/plot run without error and plot is a ggplot", {
  res <- .fit(target_r = 0.10)
  expect_output(print(res), "Trial Sequential Analysis")
  expect_output(summary(res), "Abbreviations")
  p <- plot(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot label placement follows the sign of the final Z-score", {
  res <- .fit(target_r = 0.10)
  ## annotate("text", ...) stores `label` in aes_params (not in the layer
  ## data), while x/y live in the layer data.
  label_y <- function(p, pattern) {
    for (ly in p$layers) {
      lab <- ly$aes_params$label
      if (is.null(lab) && is.data.frame(ly$data) && "label" %in% names(ly$data))
        lab <- ly$data$label
      if (length(lab) >= 1L && is.character(lab) && grepl(pattern, lab[1])) {
        yy <- ly$data$y
        if (is.null(yy)) yy <- ly$aes_params$y
        return(as.numeric(yy[1]))
      }
    }
    NA_real_
  }
  place <- function(r) {
    p <- suppressWarnings(plot(r))
    c(daris = label_y(p, "^Theoretical DARIS"),
      accrued = label_y(p, "^Participants accrued"))
  }
  pos <- res
  zz <- pos$cumulative$Z
  pos$cumulative$Z <- abs(zz)
  neg <- res
  neg$cumulative$Z <- -abs(zz)
  yp <- place(pos)
  yn <- place(neg)
  ## Z-curve positive: DARIS label low, participants label high
  expect_lt(yp[["daris"]], 0)
  expect_gt(yp[["accrued"]], 0)
  ## Z-curve negative: DARIS label high, participants label low (unchanged)
  expect_gt(yn[["daris"]], 0)
  expect_lt(yn[["accrued"]], 0)
  ## user-supplied positions still override the defaults
  p <- suppressWarnings(plot(pos, daris_label_y = 1.5, participants_label_y = -1.5))
  expect_equal(label_y(p, "^Theoretical DARIS"), 1.5)
  expect_equal(label_y(p, "^Participants accrued"), -1.5)
})

test_that("order_by sorts ascending and affects the cumulative trajectory", {
  d <- as.data.frame(readxl::read_excel(.path()))
  reversed <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  res_default <- suppressWarnings(suppressMessages(
    tsa_cor(reversed, target_r = 0.10, verbose = FALSE)))
  res_ordered <- suppressWarnings(suppressMessages(
    tsa_cor(reversed, target_r = 0.10, verbose = FALSE, order_by = "Year")))
  expect_identical(res_ordered$data$Year, sort(d$Year))
  ## same final pooled estimate (order-invariant), different intermediate path
  expect_equal(tail(res_default$cumulative$estimate, 1),
               tail(res_ordered$cumulative$estimate, 1), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(res_default$cumulative$Z, res_ordered$cumulative$Z)))
})

test_that("cor_type accepts the 'r'/'rho' shorthands, case- and whitespace-insensitively", {
  ref <- .fit(target_r = 0.10, cor_type = "pearson")
  r1  <- .fit(target_r = 0.10, cor_type = "r")
  r2  <- .fit(target_r = 0.10, cor_type = "R")
  r3  <- .fit(target_r = 0.10, cor_type = " Pearson ")
  for (alt in list(r1, r2, r3)) {
    expect_identical(alt$parameters$cor_type, "pearson")
    expect_equal(alt$pooled$r, ref$pooled$r, tolerance = 1e-10)
  }

  ref_s <- .fit(target_r = 0.10, cor_type = "spearman", se_source = "n")
  s1 <- .fit(target_r = 0.10, cor_type = "rho", se_source = "n")
  s2 <- .fit(target_r = 0.10, cor_type = " RHO ", se_source = "n")
  for (alt in list(s1, s2)) {
    expect_identical(alt$parameters$cor_type, "spearman")
    expect_equal(alt$pooled$r, ref_s$pooled$r, tolerance = 1e-10)
  }

  expect_error(.fit(target_r = 0.10, cor_type = "kendall"), "should be one of")
})
