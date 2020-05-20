#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------
// Log-Likelihood for line list data
// [[Rcpp::export]]
Rcpp::List LineList_loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {

  // extract free parameters
  double ma1 = params["ma1"];
  double ma2 = params["ma2"];
  double mod = params["mod"];
  double sod = params["sod"];
  double mor = params["mor"];
  double sor = params["sor"];

  // read in data
  std::vector<int> death = Rcpp::as< std::vector<int> >(data["death_interval"]);
  std::vector<int> recov = Rcpp::as< std::vector<int> >(data["recovery_interval"]);
  std::vector<int> deathgroup = Rcpp::as< std::vector<int> >(data["death_group"]);
  std::vector<int> recovgroup = Rcpp::as< std::vector<int> >(data["recovery_group"]);

  // extract fixed and storage items
  std::vector<int> ifrparams = Rcpp::as< std::vector<int> >(misc["IFRparams"]);
  std::vector<double> ma(ifrparams.size());
  ma[0] = ma1;
  ma[1] = ma2;

  // get loglikelihood
  double loglik = 0.0;

  // loop through deaths
  for (int i = 0; i < ifrparams.size(); i++) {
    for (int j = 0; j < death.size(); j++) {
      if (ifrparams[i] == deathgroup[j]) {

        double l1 = R::pgamma(death[j] + 1, 1/(sod*sod), mod*sod*sod, true, false);
        double l2 = R::pgamma(death[j], 1/(sod*sod), mod*sod*sod, true, false);
        loglik += log(ma[i]) + log(l1 - l2);

      }
    }
  }

  // loop through recovery
  for (int i = 0; i < ifrparams.size(); i++) {
    for (int j = 0; j < recov.size(); j++) {
      if (ifrparams[i] == recovgroup[j]) {

        double l1 = R::pgamma(recov[j] + 1, 1/(sor*sor), mor*sor*sor, true, false);
        double l2 = R::pgamma(recov[j], 1/(sor*sor), mor*sor*sor, true, false);
        loglik += log(1 - ma[i]) + log(l1 - l2);

      }
    }
  }

  // catch underflow
  if (!std::isfinite(loglik)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    loglik = -OVERFLO_DOUBLE;
  }

  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("LogLik") = loglik);
  return ret;
}
