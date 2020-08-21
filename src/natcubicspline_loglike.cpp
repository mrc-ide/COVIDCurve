#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------
// Log-Likelihood for Aggregate Expected Deaths with a Natural Cubic Spline for the Incidence Curve and a Gamma distribution for the onset-to-death course
// [[Rcpp::export]]
Rcpp::List natcubspline_loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {

  // extract misc items
  // items needed for death spline
  int n_knots = misc["n_knots"];
  int rcensor_day = misc["rcensor_day"];
  int days_obsd = misc["days_obsd"];

  // age and region strata lengths
  int agestratlen = misc["agestratlen"];
  int rgnstratlen = misc["rgnstratlen"];

  // items needed serology
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

  // extract population demographics and put in right format
  std::vector<double> count_marg_Rdemog = Rcpp::as< std::vector<double> >(misc["countmarginal_Rgn_demog"]);
  std::vector<double> count_marg_Ademog = Rcpp::as< std::vector<double> >(misc["countmarginal_Age_demog"]);

  std::vector<int> rawdemog = Rcpp::as< std::vector<int> >(misc["demog"]);
  std::vector<std::vector<int>> demog(agestratlen, std::vector<int>(rgnstratlen));
  // Ages 1,2,3, 1,2,3; region A,A,A,  B,B,B
  int iter = 0;
  for (int r = 0; r < rgnstratlen; r++) {
    for (int a = 0; a < agestratlen; a++) {
      demog[a][r] = rawdemog[iter];
      iter++;
    }
  }

  // storage items
  std::vector<double> ma(agestratlen);
  std::vector<double> Ane(agestratlen);
  std::vector<double> Rne(rgnstratlen);
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
  Ane[0] = params["Ane1"];
  Ane[1] = params["Ane2"];
  Ane[2] = params["Ane3"];
  Rne[0] = params["Rne1"];
  //Rne[0] = params["Rne2"];

  // gamma look up table
  std::vector<double> pgmms(days_obsd + 1);
  for (int i = 0; i < (days_obsd+1); i++) {
    pgmms[i] = R::pgamma(i, 1/pow(sod,2), mod*pow(sod,2), true, false);
  }


  //........................................................
  // L1 -- Deaths Shape Section
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
    // check if stratified infections exceed stratified population denominator
    bool popN_pass = true;
    std::vector<std::vector<double>> cum_infxn_check(agestratlen, std::vector<double>(rgnstratlen));
    for (int i = 0; i < days_obsd; i++) {
      for (int r = 0; r < rgnstratlen; r++) {
        for (int a = 0; a < agestratlen; a++) {
          cum_infxn_check[a][r] += Ane[a] * Rne[r] * infxn_spline[i];
        }
      }
    }

    for (int r = 0; r < rgnstratlen; r++) {
      for (int a = 0; a < agestratlen; a++) {
        if (cum_infxn_check[a][r] > demog[a][r]) {
          popN_pass = false;
        }
      }
    }

    if (popN_pass) {
      // loop through days and TOD integral to get time-delayed area under the death curve
      std::vector<double> auc(days_obsd);
      for (int i = 0; i < days_obsd; i++) {
        for (int j = i+1; j < (days_obsd + 1); j++) {
          int delta = j - i - 1;
          auc[j-1] += infxn_spline[i] * (pgmms[delta + 1] - pgmms[delta]);
        }
      }

      // get the global number of daily deaths and then calculate the global expected deaths per days
      std::vector<double> gdexpd(days_obsd);
      std::vector<std::vector<std::vector<double>>> strata_death(days_obsd, std::vector<std::vector<double>>(agestratlen, std::vector<double>(rgnstratlen)));
      double gdeaths = 0.0;
      for (int i = 0; i < days_obsd; i++) {
        for (int r = 0; r < rgnstratlen; r++) {
          for (int a = 0; a < agestratlen; a++) {
            gdexpd[i] += auc[i] * Ane[a] * Rne[r] * ma[a];
            // also store stratified deaths for marginal later
            strata_death[i][a][r] = auc[i] * Ane[a] * Rne[r] * ma[a];
          }
        }
        // also store global cumulative deaths for marginal later
        gdeaths += gdexpd[i];
      }
      // Expectation
      double L1death_loglik = 0.0;
      // extract observed data
      std::vector<int> gobsd = Rcpp::as< std::vector<int> >(data["obs_deaths"]);
      for (int  i = 0; i < days_obsd; i++) {
        // +1 for 1-based days
        if ((i+1) < rcensor_day) {
          if (gobsd[i] != -1) {
            L1death_loglik += R::dpois(gobsd[i], gdexpd[i], true);
          }
        }
      }



      //........................................................
      // L2 -- Marginal Age Deaths (Proportions)
      //........................................................
      std::vector<double> cum_age_marg_deaths(agestratlen);
      for (int a = 0; a < agestratlen; a++) {
        for (int i = 0; i < days_obsd; i++) {
          // +1 for 1-based days
          if ((i+1) < rcensor_day) {
            for (int r = 0; r < rgnstratlen; r++) {
              cum_age_marg_deaths[a] += strata_death[i][a][r];
            }
          }
        }
      }

      // Expectation
      double L2cum_age_marg_deaths_loglik = 0.0;
      // extract observed data
      std::vector<double> paobsd = Rcpp::as< std::vector<double> >(data["prop_age_obs_deaths"]);
      for (int a = 0; a < agestratlen; a++) {
        L2cum_age_marg_deaths_loglik += R::dbinom(round(cum_age_marg_deaths[a]), round(gdeaths), paobsd[a], true);
      }



      //........................................................
      // L3 -- Marginal Region Deaths (Proportions)
      //........................................................
      std::vector<double> cum_rgn_marg_deaths(rgnstratlen);
      for (int r = 0; r < rgnstratlen; r++) {
        // +1 for 1-based days
        for (int i = 0; i < days_obsd; i++) {
          if ((i+1) < rcensor_day) {
            for (int a = 0; a < agestratlen; a++) {
              cum_rgn_marg_deaths[r] += strata_death[i][a][r];
            }
          }
        }
      }

      // Expectation
      double L3cum_rgn_marg_deaths_loglik = 0.0;
      // extract observed data
      std::vector<double> probsd = Rcpp::as< std::vector<double> >(data["prop_rgn_obs_deaths"]);
      for (int r = 0; r < rgnstratlen; r++) {
        L3cum_rgn_marg_deaths_loglik += R::dbinom(round(cum_rgn_marg_deaths[r]), round(gdeaths), probsd[r], true);
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

      // // seroconversion by strata look up table
      std::vector<std::vector<std::vector<double>>> sero_con_num_full(max_seroday_obsd, std::vector<std::vector<double>>(agestratlen, std::vector<double>(rgnstratlen)));
      // // loop through and split infection curve by strata and by number of seroconversion study dates
      for (int a = 0; a < agestratlen; a++) {
        for (int r = 0; r < rgnstratlen; r++){
          for (int i = 0; i < max_seroday_obsd; i++) {
            // go to the "end" of the day
            for (int j = i+1; j < (max_seroday_obsd + 1); j++) {
              int time_elapsed = j - i - 1;
              sero_con_num_full[j-1][a][r] += infxn_spline[i] * Ane[a] * Rne[r] * cum_hazard[time_elapsed];
            }
          }
        }
      }

      // get average over serostudy data
      std::vector<std::vector<std::vector<double>>> sero_con_num(n_sero_obs, std::vector<std::vector<double>>(agestratlen, std::vector<double>(rgnstratlen)));
      for (int i = 0; i < n_sero_obs; i++) {
        for (int a = 0; a < agestratlen; a++) {
          for (int r = 0; r < rgnstratlen; r++) {
            for (int k = sero_survey_start[i]; k <= sero_survey_end[i]; k++) {
              // days are 1 based
              sero_con_num[i][a][r] += sero_con_num_full[k-1][a][r];
            }
            // now get average (+1 because days are 1-based)
            sero_con_num[i][a][r] =  sero_con_num[i][a][r]/(sero_survey_end[i] - sero_survey_start[i] + 1);
          }
        }
      }


      //........................................................
      // L4 -- Marginal Age Seroprevalences
      //........................................................
      // unpack serology observed data
      double L4_agesero_loglik = 0.0;
      std::vector<int> agedatpos_raw = Rcpp::as< std::vector<int> >(data["age_obs_serologypos"]);
      std::vector<int> agedatn_raw = Rcpp::as< std::vector<int> >(data["age_obs_serologyn"]);
      // recast datpos
      std::vector<std::vector<int>> age_datpos(n_sero_obs, std::vector<int>(agestratlen));
      std::vector<std::vector<int>> age_datn(n_sero_obs, std::vector<int>(agestratlen));
      int ageseroiter = 0;
      for (int i = 0; i < n_sero_obs; i++) {
        for (int a = 0; a < agestratlen; a++) {
          age_datpos[i][a] = agedatpos_raw[ageseroiter];
          age_datn[i][a] = agedatn_raw[ageseroiter];
          ageseroiter++;
        }
      }
      // get marginal age_sero_con
      std::vector<std::vector<double>> age_sero_con_num(n_sero_obs, std::vector<double>(agestratlen));
      for (int i = 0; i < n_sero_obs; i++) {
        for (int a = 0; a < agestratlen; a++) {
          for (int r = 0; r < rgnstratlen; r++) {
            age_sero_con_num[i][a] += sero_con_num[i][a][r];
          }
        }
      }

      // loop through sero likelihood
      for (int i = 0; i < n_sero_obs; i++) {
        for (int a = 0; a < agestratlen; a++) {
          if (age_datpos[i][a] != -1 | age_datn[i][a] != -1 ) {
            // Gelman Estimator for numerical stability
            double obs_prev = sens*(age_sero_con_num[i][a]/count_marg_Ademog[a]) + (1-spec)*(1 - (age_sero_con_num[i][a]/count_marg_Ademog[a]));
            L4_agesero_loglik += R::dbinom(age_datpos[i][a], age_datn[i][a], obs_prev, true);
          }
        }
      }




      //........................................................
      // L5 -- Marginal Regional Seroprevalences
      //........................................................
      // unpack serology observed data
      double L5_rgnsero_loglik = 0.0;
      std::vector<int> rgndatpos_raw = Rcpp::as< std::vector<int> >(data["rgn_obs_serologypos"]);
      std::vector<int> rgndatn_raw = Rcpp::as< std::vector<int> >(data["rgn_obs_serologyn"]);
      // recast datpos
      std::vector<std::vector<int>> rgn_datpos(n_sero_obs, std::vector<int>(rgnstratlen));
      std::vector<std::vector<int>> rgn_datn(n_sero_obs, std::vector<int>(rgnstratlen));
      int rgnseroiter = 0;
      for (int i = 0; i < n_sero_obs; i++) {
        for (int r = 0; r < rgnstratlen; r++) {
          rgn_datpos[i][r] = rgndatpos_raw[rgnseroiter];
          rgn_datn[i][r] = rgndatn_raw[rgnseroiter];
          rgnseroiter++;
        }
      }
      // get marginal rgn_sero_con
      std::vector<std::vector<double>> rgn_sero_con_num(n_sero_obs, std::vector<double>(rgnstratlen));
      for (int i = 0; i < n_sero_obs; i++) {
        for (int r = 0; r < rgnstratlen; r++) {
          for (int a = 0; a < agestratlen; a++) {
            rgn_sero_con_num[i][r] += sero_con_num[i][a][r];
          }
        }
      }

      // loop through sero likelihood
      for (int i = 0; i < n_sero_obs; i++) {
        for (int r = 0; r < rgnstratlen; r++) {
          if (rgn_datpos[i][r] != -1 | rgn_datn[i][r] != -1 ) {
            // Gelman Estimator for numerical stability
            double obs_prev = sens*(rgn_sero_con_num[i][r]/count_marg_Rdemog[r]) + (1-spec)*(1 - (rgn_sero_con_num[i][r]/count_marg_Rdemog[r]));
            L5_rgnsero_loglik += R::dbinom(rgn_datpos[i][r], rgn_datn[i][r], obs_prev, true);
          }
        }
      }




      //........................................................
      // // bring together
      //........................................................
      loglik = L1death_loglik + L2cum_age_marg_deaths_loglik + L3cum_rgn_marg_deaths_loglik + L4_agesero_loglik + L5_rgnsero_loglik;
      // catch underflow
      if (!std::isfinite(loglik)) {
        loglik = -OVERFLO_DOUBLE;
      }
      // end cumulative vs. popN check
    }
    // end node_x check
  }

  // // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("LogLik") = loglik);
  return ret;
}
