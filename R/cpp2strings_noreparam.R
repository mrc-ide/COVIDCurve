#' @title Liftover User Param Df to Prior Function for Dr. Jacoby
# TODO -- can't inherit params from R6 class so copy and paste

make_user_logprior_noreparam <- function(modinf) {
  paramdf <- modinf$paramdf
  Infxnparams <- modinf$Infxnparams
  IFRparams <- modinf$IFRparams

  Infxnparams <- paramdf[paramdf$name %in% Infxnparams, ]
  IFRparams <- paramdf[paramdf$name %in% IFRparams, ]

  #..................
  # gaussian for infxnpts
  #..................
  Infxnextractparams <- sapply(Infxnparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makenormpriors <- mapply(function(param, d1, d2){
    paste0("R::dnorm(", param, ",", d1, ",", d2, ",", "true) +")
  }, param = Infxnparams$name, d1 = Infxnparams$dsc1, d2 = Infxnparams$dsc2)

  #..................
  # beta for IFRparams
  #..................
  IFRextractparams <- sapply(IFRparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makebetapriors <- mapply(function(param, d1, d2){
    paste0("R::dbeta(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = IFRparams$name, d1 = IFRparams$dsc1, d2 = IFRparams$dsc2)

  #..................
  # bring together
  #..................
  extractparams <- c(Infxnextractparams, IFRextractparams)
  makebetapriors[length(makebetapriors)] <- sub("\\+", ";", makebetapriors[length(makebetapriors)])
  priors <- c("double ret =", makenormpriors, makebetapriors)

  ret <- c("SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {",
           extractparams,
           priors,
           "if (!std::isfinite(ret)) { const double OVERFLO_DOUBLE = DBL_MAX/100.0; ret = -OVERFLO_DOUBLE; }",
           "return Rcpp::wrap(ret);",
           "}"
           )

  out <- capture.output(cat(ret))
  return(out)
}

#' @title Liftover User Param Df to Likelihood Function for Dr. Jacoby
# TODO -- can't inherit params from R6 class so copy and paste

make_user_loglike_noreparam <- function(modinf) {

  paramdf <- modinf$paramdf
  Infxnparams <- modinf$Infxnparams
  IFRparams <- modinf$IFRparams

  #..................
  # extract misc
  #..................
  extmisc <- "std::vector<double> pa = Rcpp::as< std::vector<double> >(misc[\"pa\"]); std::vector<double> pgmms = Rcpp::as< std::vector<double> >(misc[\"pgmms\"]); bool level = misc[\"level\"]; std::vector<double> node_x = Rcpp::as< std::vector<double> >(misc[\"knots\"]);"

  #..................
  # extract inputs
  #..................
  params <- sapply(paramdf$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  #..................
  # storage items
  #..................
  storageitems <- "int agelen = pa.size(); std::vector<double>ma(agelen); std::vector<double> node_y(node_x.size());"
  #..................
  # Fill storage items
  #..................
  Infxnparamschar <- rep(NA, length(Infxnparams))
  for (i in 1:length(Infxnparams)){
    Infxnparamschar[i] <- paste(paste0("node_y[", i-1, "]"), "=", Infxnparams[i], ";")
  }

  IFRparamschar <- rep(NA, length(IFRparams))
  for (i in 1:length(IFRparams)){
    IFRparamschar[i] <- paste(paste0("ma[", i-1, "]"), "=", IFRparams[i], ";")
  }

  #..................
  # end liftover
  #..................
  loglike <- readLines("src/NatCubic_Spline_Growth_ExpDeaths.cpp")
  loglike <- loglike[(grep("// end liftover", loglike)+1):(grep(" // return as Rcpp list", loglike)-1)]
  # remove comment lines
  commlines <- grep("//", loglike)
  loglike <- loglike[!(1:length(loglike) %in% commlines)]

  ret <- c("SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {",
           extmisc,
           params,
           storageitems,
           Infxnparamschar,
           IFRparamschar,
           loglike,
           "return Rcpp::wrap(loglik);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}


