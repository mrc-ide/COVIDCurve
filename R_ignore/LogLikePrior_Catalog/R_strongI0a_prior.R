# define cpp logprior function
r_strongI0_prior <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {

  // extract parameters
  double ma2 = params[\"ma2\"];
  double r1 = params[\"r1\"];
  double y1 = params[\"y1\"];
  double y2 = params[\"y2\"];
  double y3 = params[\"y3\"];
  double y4 = params[\"y4\"];
  double y5 = params[\"y5\"];
  //double y6 = params[\"y6\"];

  double ret = R::dunif(r1, 0.0, 1.0, true ) +
               R::dunif(ma2, 0.0, 1.0, true ) +
               R::dnorm(y1, 4.15, 0.25, true ) +
               R::dnorm(y2, 5.90, 0.25, true ) +
               R::dnorm(y3, 7.70, 0.25, true ) +
               R::dnorm(y4, 8.42, 0.25, true ) +
               R::dnorm(y5, 8.51, 0.25, true ) +
               //R::dnorm(y6, 8.52, 0.25, true ) +
               log(ma2); // account for reparam

  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"


