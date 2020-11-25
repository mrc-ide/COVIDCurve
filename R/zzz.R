#------------------------------------------------
#' @title COVIDCurve: Simple Inference of Age-Specific IFRs
#' @description In many settings, individual-based (i.e. "line-list") data are often unavailable or
#' unreliable. Here we provide a model framework that relies on aggregated death counts and seroprevalence
#' data to infer the infection fatality ratio, or the probability of death givin infection. Models are
#' fit within a Bayesian framework and account for onset-outcome delays and serologic test characteristics.
#'
#' @docType package
#' @name COVIDCurve
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib COVIDCurve, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("COVIDCurve", libpath)  # nocov
}
