#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------
// Log-Likelihood for Aggregate Expected Deaths with a Natural Cubic Spline for the Incidence Curve and a Gamma distribution for the onset-to-death course
// [[Rcpp::export]]
Rcpp::List NatCubic_SplineGrowth_loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {

    // extract misc items
  std::vector<double> pa = Rcpp::as< std::vector<double> >(misc["pa"]);
  std::vector<double> pgmms = Rcpp::as< std::vector<double> >(misc["pgmms"]);
  bool level = misc["level"];
  std::vector<double> node_x = Rcpp::as< std::vector<double> >(misc["knots"]);

  // extract serology items
  int popN = misc["popN"];
  double sens = params["sens"];
  double spec = params["spec"];
  double sero_rate = params["sero_rate"];
  double sero_day_count = params["sero_day"];
  int sero_day = std::floor(sero_day_count);

  // extract free parameters
  double y1 = params["y1"];
  double y2 = params["y2"];
  double y3 = params["y3"];
  double y4 = params["y4"];
  double y5 = params["y5"];
  double y6 = params["y6"];
  double ma2 = params["ma2"];
  double r1 = params["r1"];
  double ma1 = ma2 * r1;

  // storage items
  int agelen = pa.size();
  std::vector<double>ma(agelen);
  std::vector<double> node_y(node_x.size());
  // fill storage
  ma[0] = ma1;
  ma[1] = ma2;
  node_y[0] = y1;
  node_y[1] = y2;
  node_y[2] = y3;
  node_y[3] = y4;
  node_y[4] = y5;
  node_y[5] = y6;

  // end liftover
  //........................................................
  // Deaths Section
  //........................................................
  //.............................
  // natural cubic spline
  //.............................
  int n_knots = node_x.size();
  int n_dat = node_x[n_knots-1] - node_x[0];
  // NB, we want N spline functions (interpolants) and so we have n+1 knots
  // as part of simplifying the linear spline, function we write this denom: x_{i+1} - x_i
  std::vector<double> denom(n_knots-1);
  for (int i = 0; i < denom.size(); i++) {
    denom[i] = node_x[i+1] - node_x[i];
  }

  // NB, our intercepts are the y_coordinates of the n+1 knots -- ends of knots will serve as boundaries
  // initialize our knot 2nd derive/linear slope
  std::vector<double> z(n_knots-1);
  z[0] = 0;
  // now store piecewise slopes of second deriv so we can calculate zi's later and make sure we are smoothing through knot
  // -2 b/c first knot and last knot set to 0
  std::vector<double> m(n_knots-2);
  for (int i = 1; i <= m.size(); i++) {
    m[i-1] = (3/denom[i])*(node_y[i+1] - node_y[i]) - (3/denom[i-1])*(node_y[i] - node_y[i-1]);
  }

  // need g and k for first deriv calc of our knots, z
  std::vector<double> g(n_knots-1);
  std::vector<double> k(n_knots-1);
  g[0] = 1;
  k[0] = 0;
  for (int i = 1; i < (n_knots-1); i++) {
    // now we can get our zs and sp2s
    g[i] = 2*(node_x[i+1] - node_x[i-1]) - (denom[i-1])*(k[i-1]);
    k[i] = denom[i]/g[i];
    z[i] = (m[i-1] - denom[i-1]*z[i-1])/g[i];
  }
  // initialize our three "slopes"
  std::vector<double> sp1(n_knots-1);
  std::vector<double> sp3(n_knots-1);
  std::vector<double> sp2(n_knots);
  sp2[n_knots-1] = 0;
  // finally loop through and get our "slopes" for our interpolants
  for (int i = (n_knots-2); i >= 0; i--) {
    sp2[i] = z[i] - k[i]*sp2[i+1];
    sp1[i] = (node_y[i+1] - node_y[i])/(denom[i]) - (denom[i]*(sp2[i+1] + 2*sp2[i]))/3 ;
    sp3[i] = (sp2[i+1] - sp2[i])/(3*denom[i]);
  }

  // create infection spline
  std::vector<double> infxn_spline(n_dat);
  int node_j = 0;
  for (int i = 0; i < n_dat; i++) {
    // update curve
    infxn_spline[i] = node_y[node_j] +
      sp1[node_j] * (i - node_x[node_j]) +
      sp2[node_j] * pow((i - node_x[node_j]), 2) +
      sp3[node_j] * pow((i - node_x[node_j]), 3);

    // update node_j
    if ((node_x[0] + i) >= node_x[node_j+1]) {
      node_j++;
    }
  }

  // convert cumulative infection spline into daily infection spline
  std::vector<double> cumm_infxn_spline(infxn_spline.size());
  cumm_infxn_spline[0] = exp(infxn_spline[0]);
  for (int i = 1; i < cumm_infxn_spline.size(); i++) {
    cumm_infxn_spline[i] = exp(infxn_spline[i]) + exp(infxn_spline[i-1]);
  }

  // loop through days and TOD integral
  std::vector<double> auc(infxn_spline.size());
  for (int i = 0; i < infxn_spline.size(); i++) {
    for (int j = i+1; j < (infxn_spline.size() + 1); j++) {
      int delta = j - i - 1;
      auc[j-1] += exp(infxn_spline[i]) * (pgmms[delta + 1] - pgmms[delta]);
    }
  }

  // Expectation
  double death_loglik = 0.0;
  // True is for Cumulative Calculation
  if (level) {

    // extract data
    std::vector<int> obsd = Rcpp::as< std::vector<int> >(data["obs_deaths"]);

    // sum up to current day
    double aucsum = 0;
    for (int i = 0; i < auc.size(); i++) {
      aucsum += auc[i];
    }
    // get exp deaths per age group
    std::vector<double>expd(agelen);
    for (int a = 0; a < agelen; a++) {
      expd[a] = aucsum * pa[a] * ma[a];
    }
    // get log-likelihood over all days
    for (int a = 0; a < agelen; a++) {
      death_loglik += R::dpois(obsd[a], expd[a], true);
    }

  } else {
  // False is for Time-Series Calculation
    // get data in right format
    std::vector<int> raw = Rcpp::as< std::vector<int> >(data["obs_deaths"]);
    std::vector<std::vector<int>> obsd(infxn_spline.size(), std::vector<int>(agelen));
    int iter = 0;
    for (int i = 0; i < infxn_spline.size(); i++) {
      for (int j = 0; j < agelen; j++) {
        obsd[i][j] = raw[iter];
        iter++;
      }
    }

    // get exp deaths per age group
    std::vector<std::vector<double>> expd(infxn_spline.size(), std::vector<double>(agelen));
    for (int  i = 0; i < infxn_spline.size(); i++) {
      for (int a = 0; a < agelen; a++) {
        expd[i][a] = auc[i] * pa[a] * ma[a];
      }
    }

    // get log-likelihood over all days
    for (int  i = 0; i < infxn_spline.size(); i++) {
      for (int a = 0; a < agelen; a++) {
        death_loglik += R::dpois(obsd[i][a], expd[i][a], true);
      }
    }
  }

  //........................................................
  // Serology Section
  //........................................................
  // account for false positives
  std::vector<double> fpr(cumm_infxn_spline.size());
  for (int i = 0; i < cumm_infxn_spline.size(); i++) {
    fpr[i] = (1-spec) * (1 - cumm_infxn_spline[i]/popN);
  }

  // account for serology delay
  std::vector<double> sero_con(infxn_spline.size());
  for (int i = 0; i < infxn_spline.size(); i++) {
    for (int j = i+1; j < (infxn_spline.size() + 1); j++) {
      double delta = j - i - 1;
      sero_con[j-1] += sens*(1-exp(((-delta+1)/sero_rate)));
    }
  }
  double datpos = data["obs_serologyrate"];
  // -1 for day to being 1-based to a 0-based call
  double pos = sero_con[sero_day - 1] + fpr[sero_day - 1]*popN;
  int posint = std::round(pos);
  double sero_loglik = R::dbinom(posint, popN, datpos, true);
  // bring together
  double loglik = death_loglik + sero_loglik;

  // catch underflow
  if (!std::isfinite(loglik)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    loglik = -OVERFLO_DOUBLE;
  }

  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("auc") = auc,
                                      Rcpp::Named("pos") = pos,
                                      Rcpp::Named("sero_con") = sero_con,
                                      Rcpp::Named("death_loglik") = death_loglik,
                                      Rcpp::Named("sero_loglik") = sero_loglik,
                                      Rcpp::Named("LogLik") = loglik,
                                      Rcpp::Named("agefat") = ma,
                                      Rcpp::Named("infxn_spline") = infxn_spline,
                                      Rcpp::Named("cumm_infxn_spline") = cumm_infxn_spline,
                                      Rcpp::Named("sp1") = sp1,
                                      Rcpp::Named("sp2") = sp2,
                                      Rcpp::Named("sp3") = sp3);
  return ret;
}

