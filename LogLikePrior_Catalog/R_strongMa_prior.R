# define cpp logprior function
r_strongMa_prior <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {

  // extract parameters
  double ma2 = params[\"ma2\"];
  double r1 = params[\"r1\"];
  double y1 = params[\"y1\"];
  double y2 = params[\"y2\"];
  double y3 = params[\"y3\"];
  double y4 = params[\"y4\"];
  double y5 = params[\"y5\"];

  double ret = R::dbeta(ma2, 500, 500, true ) +
               R::dbeta(r1, 200, 800, true ) +
               R::dunif(y1, 0.0, 10.0, true ) +
               R::dunif(y2, 0.0, 10.0, true ) +
               R::dunif(y3, 0.0, 10.0, true ) +
               R::dunif(y4, 0.0, 10.0, true ) +
               R::dunif(y5, 0.0, 10.0, true );

  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"


