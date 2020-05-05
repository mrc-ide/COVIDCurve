# define cpp logprior function
flat_prior <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {

  double ret = 0.0;

  // return as SEXP
  return Rcpp::wrap(ret);
}"

