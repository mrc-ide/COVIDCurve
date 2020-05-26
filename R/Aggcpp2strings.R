#' @title Liftover User Param Df to Prior Function for Dr. Jacoby
# TODO -- can't inherit params from R6 class so copy and paste

make_user_Agg_logprior <- function(modinf, reparamIFR) {
  #..................
  # assertsions
  #..................
  assert_custom_class(modinf, "Inference-Aggregate-Model")
  #..................
  # setup
  #..................
  paramdf <- modinf$paramdf
  Infxnparams <- modinf$Infxnparams
  IFRparams <- modinf$IFRparams
  Seroparams <- modinf$Seroparams

  Infxnparams <- paramdf[paramdf$name %in% Infxnparams, ]
  IFRparams <- paramdf[paramdf$name %in% IFRparams, ]
  Seroparams <- paramdf[paramdf$name %in% Seroparams, ]

  #..................
  # gaussian for infxnpts
  #..................
  Infxnextractparams <- sapply(Infxnparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makenormpriors <- mapply(function(param, d1, d2){
    paste0("R::dnorm(", param, ",", d1, ",", d2, ",", "true) +")
  }, param = Infxnparams$name, d1 = Infxnparams$dsc1, d2 = Infxnparams$dsc2)

  if (reparamIFR) {
    #..................
    # account for reparam
    #..................
    maxMa <- IFRparams$name[which(IFRparams$init == max(IFRparams$init))]
    scalars <- IFRparams$name[IFRparams$name != maxMa]

  }

  #..................
  # beta for IFRparams
  #..................
  IFRextractparams <- sapply(IFRparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makebetapriors <- mapply(function(param, d1, d2){
    paste0("R::dbeta(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = IFRparams$name, d1 = IFRparams$dsc1, d2 = IFRparams$dsc2)

  #......................
  # Serology priors
  #......................
  Seroextractparams <- sapply(Seroparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })
  makeSerobetapriors <- mapply(function(param, d1, d2){
    paste0("R::dbeta(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = Seroparams$name[Seroparams$name != "sero_date"],
  d1 = Seroparams$dsc1[Seroparams$name != "sero_date"],
  d2 = Seroparams$dsc2[Seroparams$name != "sero_date"])
  serodateprior <- paste0("R::dnorm(sero_date,", Seroparams$dsc1[Seroparams$name == "sero_date"], ",", Seroparams$dsc2[Seroparams$name == "sero_date"], ",true) +")

  #..................
  # bring together
  #..................
  extractparams <- c(Infxnextractparams, IFRextractparams, Seroextractparams)
  if (reparamIFR) {
    priors <- c("double ret =", makenormpriors, makebetapriors, makeSerobetapriors, serodateprior,
                paste0(length(scalars), "* log(", maxMa, ");"))
  } else {
    priors <- c("double ret =", makenormpriors, makebetapriors, makeSerobetapriors, serodateprior)
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



#---------------------------------------------------------------------------------------------------------------------

#' @title Liftover User Param Df to Likelihood Function for Dr. Jacoby
# TODO -- can't inherit params from R6 class so copy and paste

make_user_Agg_loglike <- function(modinf, reparamIFR) {
  #..................
  # assertsions
  #..................
  assert_custom_class(modinf, "Inference-Aggregate-Model")
  #..................
  # setup
  #..................
  paramdf <- modinf$paramdf
  Infxnparams <- modinf$Infxnparams
  IFRparams <- modinf$IFRparams

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
  params <- gsub("double sero_date", "double sero_day_count", params)
  params <- c(params, "int sero_day = std::floor(sero_day_count);")


  #..................
  # storage items
  #..................
  storageitems <- "int agelen = pa.size(); std::vector<double>ma(agelen); std::vector<double> node_y(node_x.size());"


  #..................
  # liftover reparam vars to Mas
  #.................
  if (reparamIFR) {
    infxnparamdf <- paramdf[paramdf$name %in% IFRparams, ]
    maxMa <- infxnparamdf$name[which(infxnparamdf$init == max(infxnparamdf$init))]
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
  # check src file for commented cpp file
  loglike <-   loglike <- "int n_knots = node_x.size();
  int n_dat = node_x[n_knots-1] - node_x[0];
  std::vector<double> denom(n_knots-1);
  for (int i = 0; i < denom.size(); i++) {
    denom[i] = node_x[i+1] - node_x[i];
  }
  std::vector<double> z(n_knots-1);
  z[0] = 0;
  std::vector<double> m(n_knots-2);
  for (int i = 1; i <= m.size(); i++) {
    m[i-1] = (3/denom[i])*(node_y[i+1] - node_y[i]) - (3/denom[i-1])*(node_y[i] - node_y[i-1]);
  }
  std::vector<double> g(n_knots-1);
  std::vector<double> k(n_knots-1);
  g[0] = 1;
  k[0] = 0;
  for (int i = 1; i < (n_knots-1); i++) {
    g[i] = 2*(node_x[i+1] - node_x[i-1]) - (denom[i-1])*(k[i-1]);
    k[i] = denom[i]/g[i];
    z[i] = (m[i-1] - denom[i-1]*z[i-1])/g[i];
  }
  std::vector<double> sp1(n_knots-1);
  std::vector<double> sp3(n_knots-1);
  std::vector<double> sp2(n_knots);
  sp2[n_knots-1] = 0;
  for (int i = (n_knots-2); i >= 0; i--) {
    sp2[i] = z[i] - k[i]*sp2[i+1];
    sp1[i] = (node_y[i+1] - node_y[i])/(denom[i]) - (denom[i]*(sp2[i+1] + 2*sp2[i]))/3;
    sp3[i] = (sp2[i+1] - sp2[i])/(3*denom[i]);
  }
  std::vector<double> infxn_spline(n_dat);
  int node_j = 0;
  for (int i = 0; i < n_dat; i++) {
    infxn_spline[i] = node_y[node_j] +
      sp1[node_j] * (i - node_x[node_j]) +
      sp2[node_j] * pow((i - node_x[node_j]), 2) +
      sp3[node_j] * pow((i - node_x[node_j]), 3);
    if ((node_x[0] + i) >= node_x[node_j+1]) {
      node_j++;
    }
  }
  std::vector<double> cumm_infxn_spline(infxn_spline.size());
  cumm_infxn_spline[0] = exp(infxn_spline[0]);
  for (int i = 1; i < infxn_spline.size(); i++) {
    cumm_infxn_spline[i] = exp(infxn_spline[i]) + exp(infxn_spline[i-1]);
  }
  std::vector<double> auc(infxn_spline.size());
  for (int i = 0; i < infxn_spline.size(); i++) {
    for (int j = i+1; j < (infxn_spline.size() + 1); j++) {
      int delta = j - i - 1;
      auc[j-1] += exp(infxn_spline[i]) * (pgmms[delta + 1] - pgmms[delta]);
    }
  }
  double death_loglik = 0.0;
  if (level) {
    std::vector<int> obsd = Rcpp::as< std::vector<int> >(data[\"obs_deaths\"]);
    double aucsum = 0;
    for (int i = 0; i < auc.size(); i++) {
      aucsum += auc[i];
    }
    std::vector<double>expd(agelen);
    for (int a = 0; a < agelen; a++) {
      expd[a] = aucsum * pa[a] * ma[a];
    }
    for (int a = 0; a < agelen; a++) {
      death_loglik += R::dpois(obsd[a], expd[a], true);
    }
  } else {
    std::vector<int> raw = Rcpp::as< std::vector<int> >(data[\"obs_deaths\"]);
    std::vector<std::vector<int>> obsd(infxn_spline.size(), std::vector<int>(agelen));
    int iter = 0;
    for (int i = 0; i < infxn_spline.size(); i++) {
      for (int j = 0; j < agelen; j++) {
        obsd[i][j] = raw[iter];
        iter++;
      }
    }
    std::vector<std::vector<double>> expd(infxn_spline.size(), std::vector<double>(agelen));
    for (int  i = 0; i < infxn_spline.size(); i++) {
      for (int a = 0; a < agelen; a++) {
        expd[i][a] = auc[i] * pa[a] * ma[a];
      }
    }
    for (int  i = 0; i < infxn_spline.size(); i++) {
      for (int a = 0; a < agelen; a++) {
        death_loglik += R::dpois(obsd[i][a], expd[i][a], true);
      }
    }
  }
  std::vector<double> fpr(cumm_infxn_spline.size());
  for (int i = 0; i < cumm_infxn_spline.size(); i++) {
    fpr[i] = (1-spec) * (1 - cumm_infxn_spline[i]/popN);
  }
  std::vector<double> sero_con(infxn_spline.size());
  for (int i = 0; i < infxn_spline.size(); i++) {
    for (int j = i+1; j < (infxn_spline.size() + 1); j++) {
      double delta = j - i - 1;
      sero_con[j-1] += sens*(1-exp(((-delta+1)/sero_rate)));
    }
  }
  double datpos = data[\"obs_serologyrate\"];
  double pos = sero_con[sero_day - 1] + fpr[sero_day - 1]*popN;
  int posint = std::round(pos);
  double sero_loglik = R::dbinom(posint, popN, datpos, true);
  double loglik = death_loglik + sero_loglik;
  if (!std::isfinite(loglik)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    loglik = -OVERFLO_DOUBLE;
  }"

  loglike <- unlist(stringr::str_split(loglike, "\n", simplify = F))

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


