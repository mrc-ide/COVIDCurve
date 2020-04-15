# define cpp loglike function for cumulative data
cpp_tod_log_like_timeseries <- "SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {

  // unpack data
  std::vector<int> raw = Rcpp::as< std::vector<int> >(data[\"obs_deaths\"]);

  // extract misc
  std::vector<double> pa =  Rcpp::as< std::vector<double> >(misc[\"pa\"]);
  int min_day = Rcpp::as<int>(misc[\"min_day\"]);
  int curr_day = Rcpp::as<int>(misc[\"curr_day\"]);
  int nsteps =  Rcpp::as<int>(misc[\"nsteps\"]);

  // extract parameters
  double ma2 = params[\"ma2\"];
  double r1 = params[\"r1\"];
  double ma1 = ma2 * r1;

  // fixed parameters
  // from Verity et al. 2020 Lancet ID
  double r = 0.14;
  double m_od = 18.8;
  double s_od = 0.45;
  double alpha = 1/(s_od*s_od);
  double beta = (m_od*s_od*s_od); // R does 1/beta
  double I0 = 2.0;

  // storage items
  int agelen = pa.size();
  std::vector<double>ma(agelen);
  ma[0] = ma1;
  ma[1] = ma2;

  // get data in right format
  std::vector<std::vector<int>> obsd(curr_day, std::vector<int>(agelen));
  int iter = 0;
  for (int i = 0; i < curr_day; i++) {
    for (int j = 0; j < agelen; j++) {
      obsd[i][j] = raw[iter];
      iter++;
    }
  }
  // get expected deaths
  std::vector<std::vector<double>> expd(curr_day, std::vector<double>(agelen));
  for (int day = 0; day < curr_day; day++) {
    // integration items
    double auc = 0.0; // reset for a new day
    double lower_tail = -nsteps/100;
    double upper_tail = day + 1; // 1-based

    // get our steps from lower tail to upper tail
    std::vector<double> steps(nsteps+1);
    double step_width = (upper_tail - lower_tail)/nsteps; // as long as nsteps is large relative to up - low, this approximates the right width
    steps[0] = lower_tail;
    for (int i = 1; i <= nsteps; i++) {
      steps[i] = steps[i-1] + step_width;
      if(steps[i] > upper_tail){
        steps[i] = upper_tail; // corner case
      }
    }

    // integrate w/ trapezoidal method
    for (int i = 0; i <= nsteps; i++) {
      auc += ((steps[i+1] - steps[i]) / 2) *
                  (I0 * exp(r*steps[i]) *
                    (R::pgamma((upper_tail - steps[i] + 1), alpha, beta, true, false) -
                     R::pgamma((upper_tail - steps[i]), alpha, beta, true, false)
                   ) +
                   I0 * exp(r*steps[i+1]) *
                    (R::pgamma((upper_tail - steps[i+1] + 1), alpha, beta, true, false) -
                     R::pgamma((upper_tail - steps[i+1]), alpha, beta, true, false)
                    )
                  );
    } // end integral

    // get exp deaths per age group
    for (int a = 0; a < agelen; a++) {
      expd[day][a] = auc * pa[a] * ma[a];
    }
  } // end loop through days

  // get log-likelihood over all days
  double loglik = 0.0;
  for (int day = 0; day < curr_day; day++) {
    for (int a = 0; a < agelen; a++) {
      loglik += R::dpois(obsd[day][a], expd[day][a], true);
    }
  }

  // account for reparameterization
  loglik += log(ma2);

  // catch underflow
    if (!isfinite(loglik)) {
      const double OVERFLO_DOUBLE = DBL_MAX/100.0;
      loglik = -OVERFLO_DOUBLE;
    }


  // return as SEXP
  return Rcpp::wrap(loglik);
}"

# define cpp logprior function
cpp_tod_log_prior_timeseries <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {

  // extract parameters
  double ma2 = params[\"ma2\"];
  double r1 = params[\"r1\"];

  // fixed/derived parameters
  double r1_mean = 0.0;
  double r1_var = 5.0;
  double lnorm_r12 = log(r1_var/(r1_mean*r1_mean) + 1.0);
  double lnorm_mu = log(r1_mean) - lnorm_r12/2.0;

  // calculate logprior
  double ret = -log(9.0) + R::dlnorm(r1, lnorm_mu, sqrt(lnorm_r12), true);

  // catch underflow
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }

  // return as SEXP
  return Rcpp::wrap(ret);
}"


