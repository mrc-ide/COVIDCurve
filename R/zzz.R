#------------------------------------------------
#' @title COVIDCurve: Simple Inference of Age-Specific IFRs
#'
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
