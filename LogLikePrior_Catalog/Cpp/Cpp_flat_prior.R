# define cpp logprior function
cpp_flat_prior <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {

  // extract parameters
  double ma2 = params[\"ma2\"];
  double r1 = params[\"r1\"];

  // fixed/derived parameters
  double r1_mean = 0.0;
  double r1_var = 5.0;
  double lnorm_r12 = log(r1_var/(r1_mean*r1_mean) + 1.0);
  double lnorm_mu = log(r1_mean) - lnorm_r12/2.0;

  // calculate logprior
  double ret = -log(9.0) + R::dlnorm(r1, lnorm_mu, sqrt(lnorm_r12), true);

  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"


