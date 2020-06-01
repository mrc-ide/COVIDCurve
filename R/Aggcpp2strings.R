#' @title Create the logprior from the IFRmodel R6 Class for DrJacoby Inference
#' @param IFRmodel R6 class; Internal model object for COVIDCurve
#' @param reparamIFR logical; Whether IFRs should be reparameterized or inferred seperately
#' @noRd

make_user_Agg_logprior <- function(IFRmodel, reparamIFR) {
  #..................
  # assertsions
  #..................
  assert_custom_class(IFRmodel, "IFRmodel")
  #..................
  # setup
  #..................
  paramdf <- IFRmodel$paramdf
  Infxnparams <- IFRmodel$Infxnparams
  IFRparams <- IFRmodel$IFRparams
  Seroparams <- IFRmodel$Seroparams

  Infxnparams <- paramdf[paramdf$name %in% Infxnparams, ]
  IFRparams <- paramdf[paramdf$name %in% IFRparams, ]
  Seroparams <- paramdf[paramdf$name %in% Seroparams, ]

  #..................
  # gaussian for infxnpts
  #..................
  Infxnextractparams <- sapply(Infxnparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeinfxnpriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(", param, ",", d1, ",", d2, ",", "true) +")
  }, param = Infxnparams$name, d1 = Infxnparams$dsc1, d2 = Infxnparams$dsc2)

  if (reparamIFR) {
    #..................
    # account for reparam
    #..................
    assert_non_null(IFRmodel$maxMa, message = "Reparameterization requires expected max mortaltiy group to be indicated (i.e. maxMa)")
    maxMa <- IFRmodel$maxMa
    scalars <- IFRparams$name[IFRparams$name != maxMa]

  }

  #..................
  # beta for IFRparams
  #..................
  IFRextractparams <- sapply(IFRparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeifrpriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = IFRparams$name, d1 = IFRparams$dsc1, d2 = IFRparams$dsc2)

  #......................
  # Serology priors
  #......................
  Seroextractparams <- sapply(Seroparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })
  makeSerobetapriors <- mapply(function(param, d1, d2){
    paste0("R::dbeta(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = Seroparams$name[!Seroparams$name %in% c("sero_rate", "sero_day")],
  d1 = Seroparams$dsc1[!Seroparams$name %in% c("sero_rate", "sero_day")],
  d2 = Seroparams$dsc2[!Seroparams$name %in% c("sero_rate", "sero_day")])
  serorateprior <- paste0("R::dunif(sero_rate,", Seroparams$dsc1[Seroparams$name == "sero_rate"], ",", Seroparams$dsc2[Seroparams$name == "sero_rate"], ",true) +")
  serodateprior <- paste0("R::dunif(sero_day,", Seroparams$dsc1[Seroparams$name == "sero_day"], ",", Seroparams$dsc2[Seroparams$name == "sero_day"], ",true) +")

  #..................
  # bring together
  #..................
  extractparams <- c(Infxnextractparams, IFRextractparams, Seroextractparams)
  if (reparamIFR) {
    priors <- c("double ret =", makeinfxnpriors, makeifrpriors, makeSerobetapriors, serorateprior, serodateprior,
                paste0(length(scalars), "*log(", maxMa, ");"))
  } else {
    priors <- c("double ret =", makeinfxnpriors, makeifrpriors, makeSerobetapriors, serorateprior, serodateprior)
    priors[length(priors)] <- sub("\\) \\+$", ");", priors[length(priors)]) # trailing + sign to a semicolon
  }

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



#' @title Create the loglikelihood from the IFRmodel R6 Class for DrJacoby Inference
#' @param IFRmodel R6 class; Internal model object for COVIDCurve
#' @param reparamIFR logical; Whether IFRs should be reparameterized or inferred seperately
#' @noRd

make_user_Agg_loglike <- function(IFRmodel, reparamIFR) {
  #..................
  # assertsions
  #..................
  assert_custom_class(IFRmodel, "IFRmodel")
  #..................
  # setup
  #..................
  paramdf <- IFRmodel$paramdf
  Infxnparams <- IFRmodel$Infxnparams
  IFRparams <- IFRmodel$IFRparams

  #..................
  # extract misc
  #..................
  extmisc <- "std::vector<double> pa = Rcpp::as< std::vector<double> >(misc[\"pa\"]); std::vector<double> pgmms = Rcpp::as< std::vector<double> >(misc[\"pgmms\"]); bool level = misc[\"level\"]; std::vector<double> node_x = Rcpp::as< std::vector<double> >(misc[\"knots\"]); int popN = misc[\"popN\"];"

  #..................
  # extract inputs
  #..................
  params <- sapply(paramdf$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })
  # convert sero_date
  params <- gsub("double sero_day", "double sero_day_raw", params)
  params <- c(params, "int sero_day = std::floor(sero_day_raw);")


  #..................
  # storage items
  #..................
  storageitems <- "int agelen = pa.size(); std::vector<double>ma(agelen); std::vector<double> node_y(node_x.size());"


  #..................
  # liftover reparam vars to Mas
  #.................
  if (reparamIFR) {
    assert_non_null(IFRmodel$maxMa, message = "Reparameterization requires expected max mortaltiy group to be indicated (i.e. maxMa)")
    infxnparamdf <- paramdf[paramdf$name %in% IFRparams, ]
    maxMa <- IFRmodel$maxMa
    scalars <- infxnparamdf$name[infxnparamdf$name != maxMa]
    mavec <- rep(NA, (length(scalars)+1))
    for (i in 1:length(scalars)) {
      mavec[i] <- paste0("ma[", i-1, "] = ", scalars[i], "*", maxMa, ";")
    }
    mavec[length(mavec)] <- paste0("ma[", length(mavec)-1, "] = ", maxMa, ";")
  } else {
    mavec <- rep(NA, length(IFRparams))
    for (i in 1:length(IFRparams)){
      mavec[i] <- paste0("ma[", i-1, "]", " = ", IFRparams[i], ";")
    }
  }

  #..................
  # Fill in Infxn Pts
  #..................
  Infxnparamschar <- rep(NA, length(Infxnparams))
  for (i in 1:length(Infxnparams)){
    Infxnparamschar[i] <- paste0("node_y[", i-1, "]", " = ", Infxnparams[i], ";")
  }

  #..................
  # get loglike
  #..................
  loglike <- readLines(system.file("src/NatCubic_AggExpDeaths_loglike_cubicspline.cpp", package = "COVIDCurve", mustWork = TRUE))
  loglike_start <- grep("// Deaths Section", loglike)
  loglike_end <- grep("// return as Rcpp list", loglike)
  loglike <- loglike[loglike_start:loglike_end]

  # remove comments which cause issues when coercing to string in this format
  commlines <- grep("//", loglike)
  loglike <- loglike[! 1:length(loglike) %in% commlines]

  # end loglike and now out
  ret <- c("SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {",
           extmisc,
           params,
           storageitems,
           Infxnparamschar,
           mavec,
           loglike,
           "return Rcpp::wrap(loglik);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}


