#------------------------------------------------
#' @title COVIDCurve: Simple Inference of IFRs and I0 for Aggregate COVID Death Counts
#'
#' @description Assumes an incidence curve can be mapped on to a onset-to-death curve in order to infer
#'              the IFR or I0 of an epidemic.
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
