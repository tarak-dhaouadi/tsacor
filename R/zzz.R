## usethis namespace: start
#' @useDynLib tsacor, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
## usethis namespace: end
NULL

# z_fisher and se_z are column names of the internal `data` frame built by
# tsa_cor(), referenced via non-standard evaluation inside
# metafor::rma(yi = z_fisher, sei = se_z, data = data, ...) -- exactly as in
# stats::lm() formulas. This declaration tells R CMD check these are not
# undefined global variables.
utils::globalVariables(c("z_fisher", "se_z"))
