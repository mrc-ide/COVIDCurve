cpp_loglike <- "SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {

  // fixed parameters
  // from Verity et al. 2020 Lancet ID
  double r = 0.14;
  double m_od = 18.8;
  double s_od = 0.45;

  // unpack data
  std::vector<double> dat = Rcpp::as< std::vector<double> >(data[\"obs_death\"]);

  // extract parameters
  double ma = params[\"ma\"];
  double I0 = params[\"I0\"];

  // extract misc objects
  std::vector<double> pa = Rcpp::as< std::vector<double> >(misc[\"pa\"]);
  int min_day = Rcpp::as<int>(misc[\"min_day\"]);
  int curr_day = Rcpp::as<int>(misc[\"curr_day\"]);

  // counter iterms
  int agelen = pa.size();

  // how do I integrate
  std::vector<double> expc;

  // get log-likelihood over all days
  double ret = 0.0;
  for(int d = min_day; d < curr_day; d++){
    for(int a = 0; a < agelen; a++){
    ret += R::dpois(expc[d][a], dat[d][a], true);
    }
  }

  // account for reparameterization
  // ret = ret + (agelen - 1) * log(ma)

  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"



# define cpp logprior function
cpp_logprior <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {

  // extract parameters
  //double ma = params[\"ma\"];
  double I0 = params[\"I0\"];

  // fixed/derived parameters
  double sigma_mean = 0.7;
  double sigma_var = 1;
  double lnorm_sigma2 = log(sigma_var/(sigma_mean*sigma_mean) + 1.0);
  double lnorm_mu = log(sigma_mean) - lnorm_sigma2/2.0;

  // calculate logprior
  // ma prior is dunif 1 to 10 or -log(9)
  double ret = -log(9) + R::dlnorm(I0, lnorm_mu, sqrt(lnorm_sigma2), true);

  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"
