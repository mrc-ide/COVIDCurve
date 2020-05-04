# define cpp logprior function
flat_prior <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {

  // extract parameters
  double ma2 = params[\"ma2\"];
  double r1 = params[\"r1\"];
  double y1 = params[\"y1\"];
  double y2 = params[\"y2\"];
  double y3 = params[\"y3\"];

  double ret = R::dunif(ma2, 0.0, 1.0, true ) +
               R::dunif(r1, 0.0, 1.0, true ) +
               R::dunif(y1, 0.0, 25.0, true ) +
               R::dunif(y2, 0.0, 25.0, true ) +
               R::dunif(y3, 0.0, 5.0, true );

  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"

