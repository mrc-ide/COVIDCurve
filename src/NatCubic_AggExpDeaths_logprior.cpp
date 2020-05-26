#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------
// Log-prior for Aggregate Expected Deaths with a Natural Cubic Spline for the Incidence Curve and a Gamma distribution for the onset-to-death course
// [[Rcpp::export]]

Rcpp::List NatCubic_SplineGrowth_logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {
  double y1 = params["y1"];
  double y2 = params["y2"];
  double y3 = params["y3"];
  double r1 = params["r1"];
  double ma2 = params["ma2"];
  double sens = params["sens"];
  double spec = params["spec"];
  double sero_rate = params["sero_rate"];
  double sero_date = params["sero_date"];
  double ret = R::dnorm(y1,3,5,true) +
    R::dnorm(y2,3,5,true) +
    R::dnorm(y3,3,5,true) +
    R::dbeta(r1,200,800,true) +
    R::dbeta(ma2,500,500,true) +
    R::dbeta(sens,800,200,true) +
    R::dbeta(spec,950,50,true) +
    R::dbeta(sero_rate,40,960,true) +
    R::dnorm(sero_date,35,8,true) +
    2* log(ma2);

  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
    }
  // return as Rcpp list
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("logprior") = ret);
  return out;
}






