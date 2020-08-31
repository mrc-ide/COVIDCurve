#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------
// Log-Likelihood for Aggregate Expected Deaths with a Natural Cubic Spline for the Incidence Curve and a Gamma distribution for the onset-to-death course
// [[Rcpp::export]]
Rcpp::List natcubspline_loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {

  // extract misc items
  std::vector<double> rho = Rcpp::as< std::vector<double> >(misc["rho"]);
  std::vector<int> demog = Rcpp::as< std::vector<int> >(misc["demog"]);
  int n_knots = misc["n_knots"];
  int rcensor_day = misc["rcensor_day"];
  int days_obsd = misc["days_obsd"];
  int n_sero_obs = misc["n_sero_obs"];
  std::vector<int> sero_survey_start = Rcpp::as< std::vector<int> >(misc["sero_survey_start"]);
  std::vector<int> sero_survey_end = Rcpp::as< std::vector<int> >(misc["sero_survey_end"]);
  int max_seroday_obsd = misc["max_seroday_obsd"];

  // extract serology items
  double sens = params["sens"];
  double spec = params["spec"];
  double sero_rate = params["sero_rate"];
  // death delay params
  double mod = params["mod"];
  double sod = params["sod"];

  // storage items
  int stratlen = rho.size();
  std::vector<double>ma(stratlen);
  std::vector<double>ne(stratlen);
  std::vector<double> node_x(n_knots);
  std::vector<double> node_y(n_knots);
  // fill storage with parameters
  node_x[0] = 1;
  node_x[1] = params["x1"];
  node_x[2] = params["x2"];
  node_x[3] = params["x3"];
  node_x[4] = params["x4"];
  node_y[0] = params["y1"];
  node_y[1] = params["y2"];
  node_y[2] = params["y3"];
  node_y[3] = params["y4"];
  node_y[4] = params["y5"];
  ma[0] = params["ma1"];
  ma[1] = params["ma2"];
  ma[2] = params["ma3"];
  ne[0] = params["ne1"];
  ne[1] = params["ne2"];
  ne[2] = params["ne3"];

  //........................................................
  // Lookup Items
  //........................................................
  // rescale Ne by attack rate
  double nedom = 0.0;
  for (int i = 0; i < stratlen; i++) {
    ne[i] = ne[i] * rho[i];
    nedom += ne[i];
  }
  // Ne as a prop vector
  for (int i = 0; i < stratlen; i++) {
    ne[i] = ne[i]/nedom;
  }

  // gamma look up table
  std::vector<double> pgmms(days_obsd + 1);
  for (int i = 0; i < (days_obsd+1); i++) {
    pgmms[i] = R::pgamma(i, 1/pow(sod,2), mod*pow(sod,2), true, false);
  }

  //........................................................
  // Deaths Section
  //........................................................
  //.............................
  // knots catch
  //.............................
  const double OVERFLO_DOUBLE = DBL_MAX/100.0;
  double loglik = -OVERFLO_DOUBLE;
  bool nodex_pass = true;
  for (int i = 1; i < node_x.size(); i++) {
    if (node_x[i] <= node_x[i-1]) {
      nodex_pass = false;
    }
  }
  if (nodex_pass) {
    //.............................
    // natural cubic spline
    //.............................
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
    std::vector<double> infxn_spline(days_obsd);
    int node_j = 0;
    infxn_spline[0] = node_y[0];
    for (int i = 1; i < days_obsd; i++) {
      // update curve and have i+1 to account for fact days are 1-based
      infxn_spline[i] = node_y[node_j] +
        sp1[node_j] * ((i+1) - node_x[node_j]) +
        sp2[node_j] * pow(((i+1) - node_x[node_j]), 2) +
        sp3[node_j] * pow(((i+1) - node_x[node_j]), 3);


      // for all interpolants except (potentially) the last knot
      if (node_j < (node_x.size()-2)) {
        // update node_j
        if ((node_x[0] + i) >= node_x[node_j+1]) {
          node_j++;
        }
      }
    }

    // exponentiate infxn spline out of log space
    for (int i = 0; i < days_obsd; i++) {
      infxn_spline[i] = exp(infxn_spline[i]);
    }

    //.............................
    // popN/size catch
    //.............................
    // check if startified infections exceed stratified population denominator
    bool popN_pass = true;
    std::vector<double> cum_infxn_check(stratlen);
    for (int i = 0; i < days_obsd; i++) {
      for (int a = 0; a < stratlen; a++) {
        cum_infxn_check[a] += ne[a] * infxn_spline[i];
      }
    }

    for (int a = 0; a < stratlen; a++) {
      if (cum_infxn_check[a] > demog[a]) {
        popN_pass = false;
      }
    }

    if (popN_pass) {

      // loop through days and TOD integral
      std::vector<double> auc(days_obsd);
      for (int i = 0; i < days_obsd; i++) {
        for (int j = i+1; j < (days_obsd + 1); j++) {
          int delta = j - i - 1;
          auc[j-1] += infxn_spline[i] * (pgmms[delta + 1] - pgmms[delta]);
        }
      }

      //.............................
      // L1 Deaths Shape Expectation
      //.............................
      // items for expectation
      double L1deathshape_loglik = 0.0;
      double texpd = 0.0;
      std::vector<std::vector<double>> strata_expdailyd(days_obsd, std::vector<double>(stratlen));
      std::vector<double> expd(days_obsd);
      // read in observed deaths
      std::vector<int> obsd = Rcpp::as< std::vector<int> >(data["obs_deaths"]);
      // get log-likelihood over all days
      for (int  i = 0; i < days_obsd; i++) {
        for (int a = 0; a < stratlen; a++) {
          // store strata deaths for L2
          strata_expdailyd[i][a] = auc[i] * ne[a] * ma[a];
          // get exp deaths per day by "summing out" strata
          expd[i] += strata_expdailyd[i][a];
        }
        // a+1 to account for 1-based dates
        if ((i+1) < rcensor_day) {
          if (obsd[i] != -1) {
            L1deathshape_loglik += R::dpois(obsd[i], expd[i], true);
          }
        }
        // store total expected deaths for L2 likelihood
        texpd += expd[i];
      }

      //.............................
      // L2 Deaths proportions Expectation
      //.............................
      // items for expectation
      double L2deathprop_loglik = 0.0;
      std::vector<double> strata_expd(stratlen);
      for (int a = 0; a < stratlen; a++) {
        for (int  i = 0; i < days_obsd; i++) {
          if ((i+1) < rcensor_day) {
            strata_expd[a] += strata_expdailyd[i][a];
          }
        }
      }

      // extract observed data
      std::vector<double> paobsd = Rcpp::as< std::vector<double> >(data["prop_strata_obs_deaths"]);
      for (int a = 0; a < stratlen; a++) {
        L2deathprop_loglik += R::dbinom(round(strata_expd[a]), round(texpd), paobsd[a], true);
      }

      //........................................................
      // Serology Section
      //........................................................
      // get cumulative hazard for each day up to the latest serology observation date
      // i.e. cumulative hazard of seroconversion on given day look up table
      std::vector<double> cum_hazard(max_seroday_obsd);
      for (int d = 0; d < max_seroday_obsd; d++) {
        cum_hazard[d] = 1-exp((-(d+1)/sero_rate));
      }

      // seroconversion by strata look up table
      std::vector<std::vector<double>> sero_con_num_full(max_seroday_obsd, std::vector<double>(stratlen));
      // loop through and split infection curve by strata and by number of seroconversion study dates
      for (int a = 0; a < stratlen; a++) {
        for (int i = 0; i < max_seroday_obsd; i++) {
          // go to the "end" of the day
          for (int j = i+1; j < (max_seroday_obsd + 1); j++) {
            int time_elapsed = j - i - 1;
            sero_con_num_full[j-1][a] += infxn_spline[i] * ne[a] * cum_hazard[time_elapsed];
          }
        }
      }

      // get average over serostudy data
      std::vector<std::vector<double>> sero_con_num(n_sero_obs, std::vector<double>(stratlen));
      for (int i = 0; i < n_sero_obs; i++) {
        for (int j = 0; j < stratlen; j++) {
          for (int k = sero_survey_start[i]; k <= sero_survey_end[i]; k++) {
            // days are 1 based
            sero_con_num[i][j] += sero_con_num_full[k-1][j];
          }
          // now get average
          sero_con_num[i][j] =  sero_con_num[i][j]/(sero_survey_end[i] - sero_survey_start[i] + 1);
        }
      }

      //.............................
      // L3 Serology Proportions Expectation
      //.............................
      double L3sero_loglik = 0.0;
      // unpack serology observed data
      std::vector<double> datpos_raw = Rcpp::as< std::vector<double> >(data["obs_serologypos"]);
      std::vector<double> datn_raw = Rcpp::as< std::vector<double> >(data["obs_serologyn"]);
      // recast datpos
      std::vector<std::vector<double>> datpos(n_sero_obs, std::vector<double>(stratlen));
      std::vector<std::vector<double>> datn(n_sero_obs, std::vector<double>(stratlen));
      int seroiter = 0;
      for (int i = 0; i < n_sero_obs; i++) {
        for (int j = 0; j < stratlen; j++) {
          datpos[i][j] = datpos_raw[seroiter];
          datn[i][j] = datn_raw[seroiter];
          seroiter++;
        }
      }
      // loop through sero likelihood
      for (int i = 0; i < n_sero_obs; i++) {
        for (int j = 0; j < stratlen; j++) {
          if (datpos[i][j] != -1 | datn[i][j] != -1 ) {
            // Gelman Estimator for numerical stability
            double obs_prev = sens*(sero_con_num[i][j]/demog[j]) + (1-spec)*(1 - (sero_con_num[i][j]/demog[j]));
            L3sero_loglik += R::dbinom(datpos[i][j], datn[i][j], obs_prev, true);
          }
        }
      }
      // bring together
      loglik = L1deathshape_loglik + L2deathprop_loglik + L3sero_loglik;

      // catch underflow
      if (!std::isfinite(loglik)) {
        loglik = -OVERFLO_DOUBLE;
      }

      // end cumulative vs. popN check
    }
    // end node_x check
  }

  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("LogLik") = loglik);
  return ret;
}
