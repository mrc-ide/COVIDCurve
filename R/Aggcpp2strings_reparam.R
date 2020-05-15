#' @title Liftover User Param Df to Prior Function for Dr. Jacoby
# TODO -- can't inherit params from R6 class so copy and paste

make_user_Agg_logprior_reparam <- function(modinf) {
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
  extractparams <- c(Infxnextractparams, IFRextractparams)
  priors <- c("double ret =", makenormpriors, makebetapriors,
              paste0(length(scalars), "* log(", maxMa, ");"))

  ret <- c("SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {",
           extractparams,
           priors,
           "if (!std::isfinite(ret)) { const double OVERFLO_DOUBLE = DBL_MAX/100.0; ret = -OVERFLO_DOUBLE; }",
           #"Rcout << \"The LP is \" << ret << std::endl;",
           "return Rcpp::wrap(ret);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}



#---------------------------------------------------------------------------------------------------------------------

#' @title Liftover User Param Df to Likelihood Function for Dr. Jacoby
# TODO -- can't inherit params from R6 class so copy and paste

make_user_Agg_loglike_reparam <- function(modinf) {
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
  # liftover reparam vars to Mas
  #.................
  infxnparamdf <- paramdf[paramdf$name %in% IFRparams, ]
  maxMa <- infxnparamdf$name[which(infxnparamdf$init == max(infxnparamdf$init))]
  scalars <- infxnparamdf$name[infxnparamdf$name != maxMa]
  mavec <- rep(NA, (length(scalars)+1))
  for (i in 1:length(scalars)) {
    mavec[i] <- paste0("ma[", i-1, "] = ", scalars[i], "*", maxMa, ";")
  }
  mavec[length(mavec)] <- paste0("ma[", length(mavec)-1, "] = ", maxMa, ";")


  #..................
  # Fill in Infxn Pts
  #..................
  Infxnparamschar <- rep(NA, length(Infxnparams))
  for (i in 1:length(Infxnparams)){
    Infxnparamschar[i] <- paste(paste0("node_y[", i-1, "]"), "=", Infxnparams[i], ";")
  }


  #..................
  # end liftover
  #..................
  loglike <- "int n_knots = node_x.size();
  int n3 = node_x[n_knots - 1] - node_x[0];

  std::vector<double> h(n_knots-1);
  for (int i = 0; i < h.size(); i++) {
    h[i] = node_x[i+1] - node_x[i];
  }
  std::vector<double> alpha(n_knots-1);
  for (int i = 1; i < h.size(); i++) {
    alpha[i] = (3/h[i])*(node_y[i+1] - node_y[i]) - (3/h[i-1])*(node_y[i] - node_y[i-1]);
  }

  std::vector<double> l(n_knots);
  l[0] = 1;
  l[(l.size()-1)] = 1;
  std::vector<double> u(n_knots-1);
  u[0] = 0;
  std::vector<double> z(n_knots);
  z[0] = 0;
  z[(z.size()-1)] = 0;
  std::vector<double> c(n_knots);
  c[(c.size()-1)] = 0;

  for (int i = 1; i < h.size(); i++) {
    l[i] = 2 * (node_x[i+1] - node_x[i-1]) - h[i-1]*u[i-1];
    u[i] = h[i]/l[i];
    z[i] = (alpha[i] - h[i-1]*z[i-1])/l[i];
  }

  std::vector<double> b(n_knots-1);
  std::vector<double> d(n_knots-1);
  for (int i = ((n_knots-1)-1); i >= 0; i--) {
    c[i] = z[i] - u[i]*c[i+1];
    b[i] = (node_y[i+1] - node_y[i])/(h[i]) - (h[i] * (c[i+1] + 2*c[i]))/3;
    d[i] = (c[i+1] - c[i])/(3*h[i]);
  }

  std::vector<double> infxn_spline(n3);
  infxn_spline[0] = node_y[0];
  int node_j = 0;
  for (int i = 1; i < n3; ++i) {

    infxn_spline[i] = node_y[node_j] +
      b[node_j] * (i - node_x[node_j]) +
      c[node_j] * pow((i - node_x[node_j]), 2) +
      d[node_j] * pow((i - node_x[node_j]), 3);

    if ((node_x[0] + i) >= node_x[node_j+1]) {
      node_j++;
    }
  }

  std::vector<double>auc(infxn_spline.size());
  for (int i = 0; i < infxn_spline.size(); i++) {
    for (int j = i+1; j < (infxn_spline.size() + 1); j++) {
      int delta = j - i - 1;
      auc[j-1] += exp(infxn_spline[i]) * (pgmms[delta + 1] - pgmms[delta]);
    }
  }

  double loglik = 0.0;
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
      loglik += R::dpois(obsd[a], expd[a], true);
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
        loglik += R::dpois(obsd[i][a], expd[i][a], true);
      }
    }
  }

  if (!isfinite(loglik)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    loglik = -OVERFLO_DOUBLE;
  }
  "
  loglike <- unlist(stringr::str_split(loglike, "\n", simplify = F))

  ret <- c("SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {",
           extmisc,
           params,
           storageitems,
           Infxnparamschar,
           mavec,
           loglike,
           #"Rcout << \"The LL is \" << loglik << std::endl;",
           "return Rcpp::wrap(loglik);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}


