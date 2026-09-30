## 0.2.7.15: regression tests for the strict-convergence regression of
## 0.2.7.13/14 and for schedules beyond the handful of looks the earlier tests
## used.
##
## In 0.2.7.13/14, searchfunc() threw whenever a beta search could not
## converge. During the design route's information-scale root search that
## happens routinely at the UPPER edge of a bracketing window (the futility
## spend asked for exceeds the probability mass still alive, i.e. the futility
## bound would lie beyond the efficacy wall), so the bracketing failed
## ("no root bracket") for many schedules -- including a real 37-look
## (40-study) schedule -- and tsa_cor() fell back to the legacy engine.
## Reference values below come from the compiled core driven through the same
## orchestration (tools/cpp_vs_py.py, section 7); the 37-look values were
## additionally checked digit-for-digit against a real tsa_cor() 40-study
## printout (efficacy 5.429581 / 2.166182; futility 0.2044558 / 1.9381223 /
## 2.0492984).

user37 <- c(0.02538682, 0.05287274, 0.07337192, 0.10764900, 0.12546251,
            0.16341297, 0.18678300, 0.21645033, 0.23176335, 0.26724431,
            0.29263113, 0.32011705, 0.34061623, 0.37489331, 0.39270682,
            0.43065727, 0.45402731, 0.48369464, 0.49900766, 0.53448862,
            0.55987544, 0.58736136, 0.60786054, 0.64213762, 0.65995113,
            0.69790158, 0.72127162, 0.75093894, 0.76625197, 0.80173292,
            0.82711975, 0.85460567, 0.87510485, 0.90938192, 0.92719543,
            0.96514589, 0.98851592)

test_that("real 37-look (40-study, target correlation 0.94) schedule: design route calibrates", {
  des <- tsacor:::.rtsa_design_bounds(user37, alpha = 0.05, beta = 0.20)
  expect_equal(length(des$timing), 38L)
  expect_equal(des$root, 1.227781, tolerance = 2e-5)
  expect_equal(des$rm_bs, 14L)
  expect_true(all(is.na(des$beta_ubound[1:14])))
  expect_equal(des$alpha_ubound[c(6, 37, 38)],
               c(5.429581, 2.166182, 2.170275), tolerance = 1e-5)
  expect_equal(des$beta_ubound[c(15, 36, 37)],
               c(0.2044558, 1.938122, 2.049299), tolerance = 1e-5)
  expect_equal(des$beta_ubound[38], des$alpha_ubound[38], tolerance = 1e-6)
  ## RTSA-parity of the alpha tolerance (0.2.7.22): RTSA's alpha search uses an
  ## ABSOLUTE tolerance of 1e-9, so a look whose cumulative spend is below that
  ## is not solved but reported as the placeholder 20 (looks 1-5 here), and the
  ## first solved look after them (look 6, 5.4296) is a little off the value a
  ## tighter tolerance would give (5.4227, error 7e-3). tsacor deliberately keeps
  ## RTSA's tolerance and placeholder; this pins that.
  expect_equal(des$alpha_ubound[1:5], rep(20, 5))
  expect_lt(des$alpha_ubound[6], 20)
})

test_that("many-look and awkward schedules calibrate (design route)", {
  cases <- list(
    list(name = "50 even looks", t = seq(0.02, 1, length.out = 50),
         root = 1.229871, rm = 19L, wall = 2.163619),
    list(name = "100 even looks", t = seq(0.01, 1, length.out = 100),
         root = 1.233670, rm = 39L, wall = 2.185481),
    list(name = "dense at start",
         t = c(seq(0.001, 0.1, length.out = 30), seq(0.15, 0.95, length.out = 10)),
         root = 1.207230, rm = 33L, wall = 2.113442),
    list(name = "dense at end",
         t = c(seq(0.05, 0.8, length.out = 10), seq(0.9, 0.999, length.out = 30)),
         root = 1.217509, rm = 4L, wall = 2.206697),
    list(name = "tiny first fraction", t = c(1e-4, 0.05, 0.2, 0.4, 0.7, 0.95),
         root = 1.161691, rm = 3L, wall = 2.085270),
    list(name = "many early looks",
         t = c(seq(0.01, 0.3, length.out = 20), 0.5, 0.7, 0.9),
         root = 1.175011, rm = 20L, wall = 2.070314)
  )
  for (cs in cases) {
    des <- tsacor:::.rtsa_design_bounds(cs$t, alpha = 0.05, beta = 0.20)
    info <- paste("schedule:", cs$name)
    expect_equal(des$root, cs$root, tolerance = 2e-5, info = info)
    expect_identical(des$rm_bs, cs$rm, info = info)
    n <- length(des$alpha_ubound)
    expect_equal(des$alpha_ubound[n], cs$wall, tolerance = 1e-5, info = info)
    ## the calibration's defining property: futility meets efficacy at t = 1
    expect_equal(des$beta_ubound[n], des$alpha_ubound[n], tolerance = 1e-6, info = info)
  }
})

test_that("an unreachable beta target is 'beyond the wall' (negative gap), not an error", {
  t100 <- seq(0.01, 1, length.out = 100)
  ub <- tsacor:::.rtsa_alpha_cpp(t100, side = 2L, alpha = 0.05)$alpha_ubound
  delta <- abs(stats::qnorm(0.025) + stats::qnorm(0.20))

  ## x = 1.25 is above the root (1.2337): 0.2.7.13/14 threw here
  lb <- tsacor:::.rtsa_beta_cpp(t100, ub, 0.20, delta, warp_root = 1.25, warn = FALSE)
  expect_gt(lb$unreachable_look, 0L)
  expect_lt(ub[100] - lb$za[100], 0)
  ## from the unreachable look onward `za` is the documented NON-PHYSICAL
  ## sentinel zb + 1 (never a boundary); before it, ordinary bounds
  idx <- lb$unreachable_look:100
  expect_equal(lb$za[idx], ub[idx] + 1)
  expect_true(all(lb$za[seq_len(lb$unreachable_look - 1L)] < ub[seq_len(lb$unreachable_look - 1L)]))

  ## at the root the pass is reachable (and the converged-pass check accepts it)
  lb0 <- tsacor:::.rtsa_beta_cpp(t100, ub, 0.20, delta, rm_bs = 0L,
                                warp_root = 1.10, warn = FALSE)
  expect_identical(lb0$unreachable_look, 0L)
  expect_gt(ub[100] - lb0$za[100], 0)

  ## an INTERIOR unreachable look also gives a negative gap (dense-at-end schedule)
  t_de <- c(seq(0.05, 0.8, length.out = 10), seq(0.9, 0.999, length.out = 30), 1)
  ub_de <- tsacor:::.rtsa_alpha_cpp(t_de, side = 2L, alpha = 0.05)$alpha_ubound
  lb_de <- tsacor:::.rtsa_beta_cpp(t_de, ub_de, 0.20, delta, warp_root = 1.25,
                                  warn = FALSE)
  expect_gt(lb_de$unreachable_look, 0L)
  expect_lt(lb_de$unreachable_look, length(t_de))
  expect_lt(ub_de[length(t_de)] - lb_de$za[length(t_de)], 0)
})

test_that(".rtsa_check_converged_pass() rejects unreachable and non-root passes", {
  ok <- list(za = c(0.5, 2.1), unreachable_look = 0L)
  expect_no_error(tsacor:::.rtsa_check_converged_pass(ok, 2.1, "test pass"))
  expect_error(
    tsacor:::.rtsa_check_converged_pass(list(za = c(0.5, 3.1), unreachable_look = 2L),
                                       2.1, "test pass"),
    "unreachable futility target"
  )
  expect_error(
    tsacor:::.rtsa_check_converged_pass(list(za = c(0.5, 1.5), unreachable_look = 0L),
                                       2.1, "test pass"),
    "does not meet the final efficacy bound"
  )
})

test_that(".rtsa_slide_root() reports the last real engine error, not just 'no root'", {
  expect_error(
    tsacor:::.rtsa_slide_root(function(x) stop("engine boom"), start = 0.95,
                             step = 0.02, max_iter = 3L),
    "last engine error: engine boom"
  )
  ## a plain absence of sign change adds no engine-error suffix
  expect_error(
    tsacor:::.rtsa_slide_root(function(x) 1, start = 0.95, step = 0.02,
                             max_iter = 3L),
    "did not converge"
  )
})

test_that("analysis route: unreachable final look is clamped (RTSA rule); interior is an error", {
  skip_if_not_installed("testthat", "3.2.0")
  t_ext <- c(0.4, 0.8, 1.2)
  fake <- function(unreach) {
    function(...) {
      list(za = c(-20, 1, 3.5), as_cum = c(0, 0.1, 0.2), as_incr = c(0, 0.1, 0.1),
           grid_collapses = 0L, grid_reversed = 0L, slow_searches = 0L,
           unreachable_look = unreach)
    }
  }
  testthat::local_mocked_bindings(rtsa_beta_boundary_cpp = fake(3L), .package = "tsacor")
  res <- tsacor:::.rtsa_analysis_bounds(t_ext, design_R = 1.2, alpha = 0.05, beta = 0.20)
  expect_true(res$final_beyond_wall)
  expect_equal(res$beta_ubound[3], res$alpha_ubound[3])      # RTSA's final clamp
  expect_true(is.na(res$beta_ubound[1]))                     # sentinel suppressed

  testthat::local_mocked_bindings(rtsa_beta_boundary_cpp = fake(2L), .package = "tsacor")
  expect_error(
    tsacor:::.rtsa_analysis_bounds(t_ext, design_R = 1.2, alpha = 0.05, beta = 0.20),
    "reaches the efficacy wall at look 2 of 3"
  )
})

test_that("0.2.7.16: engine helpers reject Inf / -Inf / NaN, not just NA", {
  expect_error(tsacor:::.rtsa_alpha_cpp(c(0.5, Inf), side = 2L, alpha = 0.05),
               "finite and strictly positive")
  expect_error(tsacor:::.rtsa_alpha_cpp(c(0.5, -Inf), side = 2L, alpha = 0.05),
               "finite and strictly positive")
  expect_error(tsacor:::.rtsa_alpha_cpp(c(0.5, NaN), side = 2L, alpha = 0.05),
               "finite and strictly positive")
  expect_error(tsacor:::.rtsa_design_bounds(c(0.3, Inf), 0.05, 0.20),
               "finite and strictly positive")
  expect_error(tsacor:::.rtsa_design_bounds(c(0.3, NaN, 0.8), 0.05, 0.20),
               "finite and strictly positive")
  expect_error(tsacor:::.rtsa_analysis_bounds(c(0.4, Inf), 1.2, 0.05, 0.20),
               "finite and strictly positive")
  expect_error(tsacor:::.rtsa_beta_cpp(c(0.5, Inf), c(2, 2), 0.2, 2.8), "finite")
})

test_that("look-spacing diagnostic: warns below 0.25% of the required information, else silent", {
  expect_no_warning(tsacor:::.rtsa_check_look_spacing(user37, "test schedule"))
  ## the FIRST look's own size is irrelevant (only increments between looks)
  expect_no_warning(tsacor:::.rtsa_check_look_spacing(c(1e-4, 0.05, 0.4, 1), "test schedule"))
  w <- tryCatch(
    tsacor:::.rtsa_check_look_spacing(c(0.3, 0.3005, 0.6, 1), "test schedule"),
    warning = function(w) conditionMessage(w))
  expect_match(w, "less than 0.25%", fixed = TRUE)
  expect_match(w, "at look 2 of 4", fixed = TRUE)
  expect_match(w, "RTSA itself refuses", fixed = TRUE)
  ## 0.25% is described as an empirical threshold for this implementation's
  ## default grid, not as an RTSA rule
  expect_match(w, "empirical warning threshold", fixed = TRUE)
  expect_match(w, "not an RTSA rule", fixed = TRUE)
  ## evidence claim kept conservative: "numerically unreliable", not "silently wrong"
  expect_match(w, "numerically unreliable", fixed = TRUE)
  expect_false(grepl("silently", w, fixed = TRUE))
  ## the engine entry point applies it (whatever the calibration then does)
  msgs <- character()
  withCallingHandlers(
    try(tsacor:::.rtsa_design_bounds(c(0.3, 0.3005, 0.6, 0.9), 0.05, 0.20), silent = TRUE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("less than 0.25%", msgs, fixed = TRUE)))
})

test_that("0.2.7.22: knife-edge -- calibration does not depend on the last bit of beta (power 0.80 etc.)", {
  ## Low-information schedule: every interim look is suppressed (rm_bs = 30) and
  ## the final look carries the whole beta spend. RTSA's `spend == beta -> za = 0`
  ## shortcut used to fire for beta = 1 - 0.8 = 0.19999999999999996 (RTSA itself is
  ## called with a literal 0.2) and made the root search fail. Roots are ~1.00003
  ## (final wall ~1.960039: with only the last look effective, ~ the fixed-sample
  ## value).
  sched <- seq(0.009, 0.26, length.out = 30)
  roots <- c(`0.80` = 1.0000404, `0.85` = 1.0000382, `0.90` = 1.0000358,
             `0.95` = 1.0000328, `0.99` = 1.0000284)
  for (pw in c(0.80, 0.85, 0.90, 0.95, 0.99)) {
    exp_root <- roots[[sprintf("%.2f", pw)]]
    for (b in c(1 - pw, round(1 - pw, 12), (1 - pw) + 1e-12, (1 - pw) - 1e-12)) {
      info <- sprintf("power %.2f, beta = %.17g", pw, b)
      des <- tsacor:::.rtsa_design_bounds(sched, alpha = 0.05, beta = b)
      expect_identical(des$rm_bs, 30L, info = info)
      expect_true(all(is.na(des$beta_ubound[1:30])), info = info)
      expect_equal(des$root, exp_root, tolerance = 2e-6, info = info)
      expect_equal(des$alpha_ubound[31], 1.960039, tolerance = 1e-6, info = info)
      expect_equal(des$beta_ubound[31], des$alpha_ubound[31], tolerance = 1e-6, info = info)
    }
  }
})


