#' @title Liftover User Param Df to Prior Function for Dr. Jacoby
# TODO -- can't inherit params from R6 class so copy and paste
make_user_LineList_logprior_reparam <- function(modinf) {
  #..................
  # assertsions
  #..................
  assert_custom_class(modinf, "Inference-LineList-Model")
  #..................
  # setup
  #..................
  paramdf <- modinf$paramdf
  IFRparams <- modinf$IFRparams
  death_mod <- modinf$DeathModparam
  death_sod <- modinf$DeathSodparam
  recov_mor <- modinf$RecovMorparam
  recov_sor <- modinf$RecovSorparam

  IFRparams <- paramdf[paramdf$name %in% IFRparams, ]
  death_mod <- paramdf[paramdf$name %in% death_mod, ]
  death_sod <- paramdf[paramdf$name %in% death_sod, ]
  recov_mor <- paramdf[paramdf$name %in% recov_mor, ]
  recov_sor <- paramdf[paramdf$name %in% recov_sor, ]

  #..................
  # norm for gamma means
  #..................
  gmmean <- rbind.data.frame(death_mod, recov_mor)
  gmmeanextractparams <- sapply(gmmean$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makenormpriors <- mapply(function(param, d1, d2){
    paste0("R::dnorm(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = gmmean$name, d1 = gmmean$dsc1, d2 = gmmean$dsc2)


  #..................
  # norm for gamma coeff var
  #..................
  gmsod <- rbind.data.frame(death_sod, recov_sor)
  gmsodextractparams <- sapply(gmsod$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makelognormpriors <- mapply(function(param, d1, d2){
    paste0("R::dlnorm(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = gmsod$name, d1 = gmsod$dsc1, d2 = gmsod$dsc2)



  #..................
  # account for reparam
  #..................
  maxMa <- IFRparams$name[which(IFRparams$init == max(IFRparams$init))]
  scalars <- IFRparams$name[IFRparams$name != maxMa]

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
  extractparams <- c(gmmeanextractparams, gmsodextractparams, IFRextractparams)
  priors <- c("double ret =", makenormpriors, makelognormpriors, makebetapriors,
              paste0(length(scalars), "* log(", maxMa, ");"))

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

#-----------------------------------------------------------------------------------------------------
#' @title Liftover User Param Df to Likelihood Function for Dr. Jacoby
# TODO -- can't inherit params from R6 class so copy and paste

make_user_LineList_loglike_reparam <- function(modinf) {

  # Note, the IFR params are the grouping parameters -- i.e. in the data, the age-specific mortaltiy
  # must be refelcted in the IFR params

  #..................
  # assertsions
  #..................
  assert_custom_class(modinf, "Inference-LineList-Model")
  #..................
  # setup/storage
  #..................
  paramdf <- modinf$paramdf
  mastore <- "std::vector<int> ifrparams = Rcpp::as< std::vector<int> >(misc[\"IFRparams\"]); std::vector<double> ma(ifrparams.size());"


  #..................
  # extract inputs
  #..................
  deathdata <- "std::vector<int> death = Rcpp::as< std::vector<int> >(data[\"death_interval\"]);"
  deathgroup <- "std::vector<int> recov = Rcpp::as< std::vector<int> >(data[\"recovery_interval\"]);"
  recovdata <- "std::vector<int> deathgroup = Rcpp::as< std::vector<int> >(data[\"death_group\"]);"
  recovgroup <- "std::vector<int> recovgroup = Rcpp::as< std::vector<int> >(data[\"recovery_group\"]);"


  extractIFRparams <- sapply(paramdf$name[paramdf$name %in% modinf$IFRparams], function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  extractmod <- paste("double mod = ", paste0("params[\"",  modinf$DeathModparam, "\"];"))
  extractsod <- paste("double sod = ", paste0("params[\"",  modinf$DeathSodparam, "\"];"))
  extractmor <- paste("double mor = ", paste0("params[\"",  modinf$RecovMorparam, "\"];"))
  extractsor <- paste("double sor = ", paste0("params[\"",  modinf$RecovSorparam, "\"];"))


  #..................
  # liftover reparam vars to Mas
  #.................
  infxnparamdf <- paramdf[paramdf$name %in% modinf$IFRparams, ]
  maxMa <- infxnparamdf$name[which(infxnparamdf$init == max(infxnparamdf$init))]
  scalars <- infxnparamdf$name[infxnparamdf$name != maxMa]
  mavec <- rep(NA, (length(scalars)+1))
  for (i in 1:length(scalars)) {
    mavec[i] <- paste0("ma[", i-1, "] = ", scalars[i], "*", maxMa, ";")
  }
  mavec[length(mavec)] <- paste0("ma[", length(mavec)-1, "] = ", maxMa, ";")

  loglike <- "double loglik = 0.0;

  for (int i = 0; i < ifrparams.size(); i++) {
    for (int j = 0; j < death.size(); j++) {
      if (ifrparams[i] == deathgroup[j]) {

        double l1 = R::pgamma(death[j] + 1, 1/(sod*sod), mod*sod*sod, true, false);
        double l2 = R::pgamma(death[j], 1/(sod*sod), mod*sod*sod, true, false);
        loglik += log(ma[i]) + log(l1 - l2);

      }
    }
  }

  for (int i = 0; i < ifrparams.size(); i++) {
    for (int j = 0; j < recov.size(); j++) {
      if (ifrparams[i] == recovgroup[j]) {

        double l1 = R::pgamma(recov[j] + 1, 1/(sor*sor), mor*sor*sor, true, false);
        double l2 = R::pgamma(recov[j], 1/(sor*sor), mor*sor*sor, true, false);
        loglik += log(1 - ma[i]) + log(l1 - l2);

      }
    }
  }

  if (!isfinite(loglik)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    loglik = -OVERFLO_DOUBLE;
  }
  "
  loglike <- unlist(stringr::str_split(loglike, "\n", simplify = F))


  #..................
  # ret
  #..................
  datainput <- c(deathdata, deathgroup, recovdata, recovgroup)
  extractparams <- c(extractIFRparams, extractmod, extractsod, extractmor, extractsor)
  ret <- c("SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {",
           datainput,
           extractparams,
           mastore,
           mavec,
           loglike,
           "return Rcpp::wrap(loglik);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}


