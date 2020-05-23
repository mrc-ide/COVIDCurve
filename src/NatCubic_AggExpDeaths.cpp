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
  double sens =  misc["sens"];
  double sero_rate =  misc["sero_rate"];
  int sero_day =  misc["sero_day"];
  double spec =  misc["spec"];

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
  int n3 = node_x[n_knots - 1] - node_x[0];

  // need alpha and h for second deriv
  std::vector<double> h(n_knots-1);
  for (int i = 0; i < h.size(); i++) {
    h[i] = node_x[i+1] - node_x[i];
  }
  std::vector<double> alpha(n_knots-1);
  for (int i = 1; i < h.size(); i++) {
    alpha[i] = (3/h[i])*(node_y[i+1] - node_y[i]) - (3/h[i-1])*(node_y[i] - node_y[i-1]);
  }

  // need items for first deriv and constraints
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
  // other coeffs now
  std::vector<double> b(n_knots-1);
  std::vector<double> d(n_knots-1);
  for (int i = ((n_knots-1)-1); i >= 0; i--) {
    c[i] = z[i] - u[i]*c[i+1];
    b[i] = (node_y[i+1] - node_y[i])/(h[i]) - (h[i] * (c[i+1] + 2*c[i]))/3;
    d[i] = (c[i+1] - c[i])/(3*h[i]);
  }

  // create CUMULATIVE infection spline
  std::vector<double> cumm_infxn_spline(n3);
  cumm_infxn_spline[0] = node_y[0];
  int node_j = 0;
  for (int i = 1; i < n3; ++i) {

    // update curve
    cumm_infxn_spline[i] = node_y[node_j] +
      b[node_j] * (i - node_x[node_j]) +
      c[node_j] * pow((i - node_x[node_j]), 2) +
      d[node_j] * pow((i - node_x[node_j]), 3);

    // update node_j
    if ((node_x[0] + i) >= node_x[node_j+1]) {
      node_j++;
    }
  }

  // convert cumulative infection spline into daily infection spline
  std::vector<double> infxn_spline(cumm_infxn_spline.size());
  infxn_spline[0] = exp(cumm_infxn_spline[0]);
  for (int i = 1; i < cumm_infxn_spline.size(); i++) {
    infxn_spline[i] = exp(cumm_infxn_spline[i]) - exp(cumm_infxn_spline[i-1]);
  }

  // loop through days and TOD integral
  std::vector<double> auc(infxn_spline.size());
  for (int i = 0; i < infxn_spline.size(); i++) {
    for (int j = i+1; j < (infxn_spline.size() + 1); j++) {
      int delta = j - i - 1;
      auc[j-1] += infxn_spline[i] * (pgmms[delta + 1] - pgmms[delta]);
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
  std::vector<double> sero_folks(infxn_spline.size()); //TEMP
  for (int i = 0; i < infxn_spline.size(); i++) {
    for (int j = i+1; j < (infxn_spline.size() + 1); j++) {
      double delta = j - i - 1;
      sero_con[j-1] += sens*(exp(1/sero_rate)-1)*exp((-delta+1)/sero_rate);
      sero_folks[j-1] += infxn_spline[i] * sens*(exp(1/sero_rate)-1)*exp((-delta+1)/sero_rate);
    }
  }
  double datpos = data["obs_serologyrate"];
  double pos = sero_con[sero_day] + fpr[sero_day]*popN;
  pos = round(pos);
  double sero_loglik = R::dbinom(pos, popN, datpos, true);
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
                                      Rcpp::Named("sero_folks") = sero_folks,
                                      Rcpp::Named("death_loglik") = death_loglik,
                                      Rcpp::Named("sero_loglik") = sero_loglik,
                                      Rcpp::Named("LogLik") = loglik,
                                      Rcpp::Named("agefat") = ma,
                                      Rcpp::Named("infxn_spline") = infxn_spline,
                                      Rcpp::Named("cumm_infxn_spline") = cumm_infxn_spline);
  return ret;
}

