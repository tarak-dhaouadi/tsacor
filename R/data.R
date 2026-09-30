#' Path to a bundled example correlation meta-analysis dataset
#'
#' Returns the file path to the example dataset bundled with the package,
#' suitable for trying out \code{\link{tsa_cor}}. It is an \code{.xlsx} sheet
#' with one row per study (20 studies, 2005-2024) and the columns
#' \code{Study}, \code{Year}, \code{Ethnicity}, \code{Age}, \code{r} (a
#' Pearson correlation or Spearman rho), \code{lbound} and \code{ubound}
#' (the lower and upper limits of its 95\% confidence interval) and
#' \code{n_subjects} (the number of subjects of the study). \code{tsa_cor()}
#' uses \code{Study}, \code{r}, \code{lbound}, \code{ubound} and
#' \code{n_subjects}; \code{Year} can be passed to \code{order_by}.
#'
#' @param dataset Which example to return. Currently only \code{"r_meta"}
#'   (default; 20 studies) is available.
#'
#' @return A character string giving the path to the example .xlsx file.
#' @examples
#' path <- tsacor_example_data()
#' d <- readxl::read_excel(path)
#' head(d)
#' @export
tsacor_example_data <- function(dataset = c("r_meta")) {
  dataset <- match.arg(dataset)
  system.file("extdata", paste0(dataset, ".xlsx"), package = "tsacor")
}
