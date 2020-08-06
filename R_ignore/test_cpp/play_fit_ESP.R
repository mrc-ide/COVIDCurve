#########################################################################
# Purpose: Deploy
#
# Author: Nicholas F. Brazeau
#########################################################################
cpp_loglike <- "SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {

  // extract misc items
  std::vector<double> rho = Rcpp::as< std::vector<double> >(misc[\"rho\"]);
  std::vector<int> demog = Rcpp::as< std::vector<int> >(misc[\"demog\"]);
  int n_knots = misc[\"n_knots\"];
  int n_sero_obs = misc[\"n_sero_obs\"];
  int rcensor_day = misc[\"rcensor_day\"];
  int days_obsd = misc[\"days_obsd\"];

  // extract serology items
  double sens = params[\"sens\"];
  double spec = params[\"spec\"];
  double sero_rate = params[\"sero_rate\"];
  double sero_day1_raw = params[\"sero_day1\"];
  double sero_day2_raw = params[\"sero_day2\"];
  std::vector<double> sero_day_raw(n_sero_obs);
  sero_day_raw[0] = sero_day1_raw;
  sero_day_raw[1] = sero_day2_raw;
  std::vector<int> sero_day(n_sero_obs);
  for (int i = 0; i < n_sero_obs; i++) {
    sero_day[i] =  std::floor(sero_day_raw[i]);
  }

  // extract free spline parameters
  double x1 = params[\"x1\"];
  double x2 = params[\"x2\"];
  double x3 = params[\"x3\"];
  double x4 = params[\"x4\"];
  double y1 = params[\"y1\"];
  double y2 = params[\"y2\"];
  double y3 = params[\"y3\"];
  double y4 = params[\"y4\"];
  double y5 = params[\"y5\"];
  // extract IFR parameters
  double ma1 = params[\"ma1\"];
  double ma2 = params[\"ma2\"];
  double ma3 = params[\"ma3\"];
  double ma4 = params[\"ma4\"];
  double ma5 = params[\"ma5\"];
  double ma6 = params[\"ma6\"];
  double ma7 = params[\"ma7\"];
  double ma8 = params[\"ma8\"];
  double ma9 = params[\"ma9\"];
  double ma10 = params[\"ma10\"];
  // extract noise parameters
  double ne1 = params[\"ne1\"];
  double ne2 = params[\"ne2\"];
  double ne3 = params[\"ne3\"];
  double ne4 = params[\"ne4\"];
  double ne5 = params[\"ne5\"];
  double ne6 = params[\"ne6\"];
  double ne7 = params[\"ne7\"];
  double ne8 = params[\"ne8\"];
  double ne9 = params[\"ne9\"];
  double ne10 = params[\"ne10\"];
  // death delay params
  double mod = params[\"mod\"];
  double sod = params[\"sod\"];

  // storage items
  int stratlen = rho.size();
  std::vector<double>ma(stratlen);
  std::vector<double>ne(stratlen);
  std::vector<double> node_x(n_knots);
  std::vector<double> node_y(n_knots);
  // fill storage
  node_x[0] = 1;
  node_x[1] = x1 * x4;
  node_x[2] = x2 * x4;
  node_x[3] = x3 * x4;
  node_x[4] = x4;
  node_y[0] = y1 * y3;
  node_y[1] = y2 * y3;
  node_y[2] = y3;
  node_y[3] = y4 * y3;
  node_y[4] = y5 * y3;
  ma[0] = ma1;
  ma[1] = ma2 * ma1;
  ma[2] = ma3* ma1;
  ma[3] = ma4* ma1;
  ma[4] = ma5* ma1;
  ma[5] = ma6* ma1;
  ma[6] = ma7* ma1;
  ma[7] = ma8* ma1;
  ma[8] = ma9* ma1;
  ma[9] = ma10* ma1;
  // break more correlation
  mod = mod * 1/spec;
  ne1 = ne1 * 1/spec;

  sero_rate = mod*sero_rate;

  ne[0] = ne1;
  ne[1] = ne2 * ne1;
  ne[2] = ne3 * ne1;
  ne[3] = ne4 * ne1;
  ne[4] = ne5 * ne1;
  ne[5] = ne6 * ne1;
  ne[6] = ne7 * ne1;
  ne[7] = ne8 * ne1;
  ne[8] = ne9 * ne1;
  ne[9] = ne10 * ne1;

  //........................................................
  // Lookup Items
  //........................................................
  // rescale ne by attack rate
  for (int i = 0; i < stratlen; i++) {
    ne[i] = ne[i] * rho[i];
  }

  // get popN for catch
  int popN = 0;
  for (int i = 0; i < stratlen; i++) {
    popN += demog[i];
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
    // initialize our three \"slopes\"
    std::vector<double> sp1(n_knots-1);
    std::vector<double> sp3(n_knots-1);
    std::vector<double> sp2(n_knots);
    sp2[n_knots-1] = 0;
    // finally loop through and get our \"slopes\" for our interpolants
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

    // get cumulative infection spline
    double cum_infxn_check = 0.0;
    for (int i = 0; i < days_obsd; i++) {
      for (int a = 0; a < stratlen; a++) {
        cum_infxn_check += ne[a] * infxn_spline[i];
      }
    }

    // check if total infections exceed population denominator
    if (cum_infxn_check <= popN) {

      // loop through days and TOD integral
      std::vector<double> auc(days_obsd);
      for (int i = 0; i < days_obsd; i++) {
        for (int j = i+1; j < (days_obsd + 1); j++) {
          int delta = j - i - 1;
          auc[j-1] += infxn_spline[i] * (pgmms[delta + 1] - pgmms[delta]);
        }
      }

      // Expectation
      double death_loglik = 0.0;
      // get data in right format
      std::vector<int> raw = Rcpp::as< std::vector<int> >(data[\"obs_deaths\"]);
      std::vector<std::vector<int>> obsd(days_obsd, std::vector<int>(stratlen));
      int iter = 0;
      for (int i = 0; i < days_obsd; i++) {
        for (int j = 0; j < stratlen; j++) {
          obsd[i][j] = raw[iter];
          iter++;
        }
      }


      std::vector<std::vector<double>> expd(days_obsd, std::vector<double>(stratlen));
      // get log-likelihood over all days
      for (int  i = 0; i < days_obsd; i++) {
        for (int a = 0; a < stratlen; a++) {
          // get exp deaths per age group
          expd[i][a] = auc[i] * ne[a] * ma[a];
          // a+1 to account for 1-based dates
          if ((a+1) < rcensor_day) {
            if (obsd[i][a] != -1) {
              death_loglik += R::dpois(obsd[i][a], expd[i][a], true);
            }
          }
        }
      }

      //........................................................
      // Serology Section
      //........................................................
      // account for serology delay -- cumulative hazard of seroconversion on given day look up table
      // days are 1-based
      std::vector<std::vector<double>> sero_con_num(n_sero_obs, std::vector<double>(stratlen));
      for (int i = 0; i < n_sero_obs; i++) {
        // get cumulative hazard for each study date
        std::vector<double> cum_hazard(sero_day[i]);
        for (int d = 0; d < sero_day[i]; d++) {
          cum_hazard[d] = 1-exp((-(d+1)/sero_rate));
        }
        // loop through and split infection curve by strata and by number of seroconversion study dates
        for (int j = 0; j < stratlen; j++) {
          // go to the \"end\" of the day
          for (int d = 0; d < sero_day[i]; d++) {
            int time_elapsed = sero_day[i] - d - 1;
            sero_con_num[i][j] += infxn_spline[d] * ne[j] * cum_hazard[time_elapsed];
          }
        }
      }

      // update now for sensitivity and false positives; -1 for day to being 1-based to a 0-based call
      double sero_loglik = 0.0;
      std::vector<double> datpos_raw = Rcpp::as< std::vector<double> >(data[\"obs_serology\"]);
      // recast datpos
      std::vector<std::vector<double>> datpos(n_sero_obs, std::vector<double>(stratlen));
      int seroiter = 0;
      for (int i = 0; i < n_sero_obs; i++) {
        for (int j = 0; j < stratlen; j++) {
          datpos[i][j] = datpos_raw[seroiter];
          seroiter++;
        }
      }

      std::vector<std::vector<double>> obs_prev(n_sero_obs, std::vector<double>(stratlen));
      for (int i = 0; i < n_sero_obs; i++) {
        for (int j = 0; j < stratlen; j++) {
          // Rogan-Gladen Estimator
          obs_prev[i][j] = sens*(sero_con_num[i][j]/demog[j]) + (1-spec)*(1 - (sero_con_num[i][j]/demog[j]));
          int posint = round(obs_prev[i][j] * demog[j]);
          sero_loglik += R::dbinom(posint, demog[j], datpos[i][j], true);
        }
      }

      // bring together
      loglik = death_loglik + sero_loglik;

      // catch underflow
      if (!std::isfinite(loglik)) {
        loglik = -OVERFLO_DOUBLE;
      }

      // end cumulative vs. popN check
    }
    // end node_x check
  }

  // return as SEXP
  return Rcpp::wrap(loglik);
}"


set.seed(1234)
library(drjacoby)
library(tidyverse)


# onset to deaths
tod_paramsdf <- tibble::tibble(name = c("mod", "sod"),
                               min  = c(10,     0.01),
                               init = c(14,     0.7),
                               max =  c(30,     3.00),
                               dsc1 = c(2.657,  -0.236),
                               dsc2 = c(0.05,   0.05))

sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1", "sero_day2"),
                                min =   c(0.83,    0.50,     0,          125,        146),
                                init =  c(0.85,    0.99,     1,        125,        146),
                                max =   c(1.00,    1.00,     2,          125,        146),
                                dsc1 =  c(123.5,   156.5,    1,        118,        139),
                                dsc2 =  c(30.5,    0.5,      2,        132,        153))

ifr_paramsdf <- tibble::tibble(name = c("ma1", "ma2",  "ma3", "ma4", "ma5",  "ma6", "ma7",  "ma8", "ma9",  "ma10"),
                               min  = rep(0, 10),
                               init = rep(0.5, 10),
                               max = rep(1, 10),
                               dsc1 = rep(0, 10),
                               dsc2 = rep(1, 10))

noise_paramsdf <- tibble::tibble(name = c("ne1", "ne2",  "ne3", "ne4", "ne5",  "ne6", "ne7",  "ne8", "ne9",  "ne10"),
                                 min  = rep(0, 10),
                                 init = rep(1, 10),
                                 max = rep(10, 10),
                                 dsc1 = rep(0, 10),
                                 dsc2 = rep(10, 10))

infxn_paramsdf <- tibble::tibble(name = paste0("y", 1:5),
                                 min  = c(0,   0,   0,  0,   0),
                                 init = c(0.5, 0.5, 10, 0.5, 0.5),
                                 max  = c(1,   1,   17, 1,   1),
                                 dsc1 = c(0,   0,   0,  0,   0),
                                 dsc2 = c(1,   1,   17, 1,   1))

knot_paramsdf <- tibble::tibble(name = paste0("x", 1:4),
                                min  = c(0,    0.33, 0.66, 175),
                                init = c(0.05, 0.40, 0.75, 185),
                                max =  c(0.33, 0.66, 0.99, 205),
                                dsc1 = c(0,    0.33, 0.66, 175),
                                dsc2 = c(0.33, 0.66, 0.99, 205))


df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf, tod_paramsdf)

#......................
# format data
#......................
dat <- readRDS("data/derived/ESP/ESP_agebands.RDS")
num_mas <- 10
dictkey <- tibble::tibble(groupvar = as.character(unlist(unique(dat$seroprevMCMC[, "ageband"]))), "Strata" = paste0("ma", 1:num_mas))
colnames(dictkey) <- c("ageband", "Strata")
# deaths
dat$deaths <- dplyr::left_join(dat$deathsMCMC, dictkey) %>%
  dplyr::select(c("ObsDay", "Strata", "Deaths"))  %>%
  dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
  dplyr::arrange(ObsDay, Strata) %>%
  dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat

# seroprev
seroprev_day_lftvr <- tibble::tibble(ObsDaymin = unique(dat$seroprevMCMC$ObsDaymin),
                                     ObsDaymax = unique(dat$seroprevMCMC$ObsDaymax),
                                     SeroDay = c("sero_day1", "sero_day2"))

dat$obs_serology <- dplyr::left_join(dat$seroprevMCMC, dictkey) %>%
  dplyr::left_join(., seroprev_day_lftvr) %>%
  dplyr::select(c("SeroDay", "Strata", "SeroPrev")) %>%
  dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
  dplyr::arrange(SeroDay, Strata) %>%
  dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat



demog <- dat$prop_pop %>%
  dplyr::left_join(., dictkey) %>%
  dplyr::select(c("Strata", "popN")) %>%
  dplyr::group_by(Strata) %>%
  dplyr::summarise(popN = round(sum(popN))) %>%
  dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
  dplyr::arrange(Strata) %>%
  dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat

#......................
# pull out bits from model
#......................
misc_list = list(rho = rep(1, 10),
                 demog = demog$popN,
                 rcensor_day = .Machine$integer.max,
                 days_obsd = 205,
                 n_knots = 5,
                 n_sero_obs = 2)

datinput <- list("obs_deaths" = dat$deaths$Deaths,
                 "obs_serology" = dat$obs_serology$SeroPrev)


#..................
# get items for model
#..................
bvec <- seq(5, 1.1, length.out = 50)

flatprior <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc){ return Rcpp::wrap(0.0); }"


modout <- drjacoby::run_mcmc(data = datinput,
                             df_params = df_params[,1:4],
                             loglike = cpp_loglike,
                             logprior = flatprior,
                             misc = misc_list,
                             chains = 3,
                             burnin = 1e3,
                             samples = 1e3,
                             rungs = 50,
                             GTI_pow = bvec,
                             silent = F)


modout <- readRDS("~/Documents/GitHub/reestimate_covidIFR_analysis/results/temp/corr.RDS")

plot_mc_acceptance(modout)
summary(modout$output$loglikelihood)
summary(modout$output$logprior)
plot_par(modout, "spec")
plot_par(modout, "sero_rate")
plot_par(modout, "ne1")


plot_par(modout, "ma1", rung = 1)
plot_par(modout, "ma2", rung = 1)
plot_par(modout, "ma3", rung = 1)
plot_par(modout, "ma4", rung = 1)
plot_par(modout, "ma5", rung = 1)
plot_par(modout, "ma6", rung = 1)


plot_par(modout, "sens")
plot_par(modout, "spec")

plot_par(modout, "mod")
plot_par(modout, "sero_rate")

plot_par(modout, "sero_day1")
plot_par(modout, "sero_day2")

plot_par(modout, "y1", rung = 1)
plot_par(modout, "y2", rung = 1)
plot_par(modout, "y3", rung = 1)
plot_par(modout, "y4", rung = 1)
plot_par(modout, "y5", rung = 1)
plot_par(modout, "x1", rung = 1)
plot_par(modout, "x2", rung = 1)
plot_par(modout, "x3", rung = 1)
plot_par(modout, "x4", rung = 1)
plot_par(modout, "ne1", rung = 1)
plot_par(modout, "ne2", rung = 1)
plot_par(modout, "ne3", rung = 1)

summary(modout$output$loglikelihood)
summary(modout$output$logprior)
modout$output[modout$output$loglikelihood == max(modout$output$loglikelihood), ]



plot_cor(modout, "ne1", "ne2", rung = 1)
plot_cor(modout, "ne2", "ne3", rung = 1)

plot_cor(modout, "y3", "spec", rung = 1)
plot_cor(modout, "y3", "ma3", rung = 1)
plot_cor(modout, "ma1", "ma3", rung = 1)



plot_mc_acceptance(modout)
drjacoby::plot_rung_loglike(modout)
drjacoby::plot_rung_loglike(modout, x_axis_type = 1, y_axis_type = 2, phase = "burnin")
drjacoby::plot_rung_loglike(modout, x_axis_type = 1, y_axis_type = 3, phase = "sampling")
drjacoby::plot_rung_loglike(modout, x_axis_type = 2, y_axis_type = 2)
drjacoby::plot_rung_loglike(modout, x_axis_type = 2, y_axis_type = 3)


# plot out
jpeg("~/Desktop/posterior_curve_draws.jpg", width = 11, height = 8, units = "in", res = 500)
library(ggplot2)
liftover <- data.frame(param = c("ma1", "ma2", "ma3"),
                       Strata = c("ma1", "ma2", "ma3"))

fatalitydataplot <- fatalitydata %>%
  dplyr::left_join(liftover, ., by = "Strata")

plot1 <- ggplot() +
  geom_pointrange(data = ifr, aes(x = param, ymin = LCI, ymax = UCI, y = median, color = param)) +
  geom_hline(data = fatalitydataplot, aes(yintercept  = IFR, group = param), color = "#3182bd", linetype = "dashed", size = 1.1) +
  facet_wrap(.~param) +
  scale_color_viridis_d() +
  theme_bw()

plot2 <- curve$plotObj +
  geom_line(data = infxns, aes(x = time, y = infxns), color = "#3182bd")

cowplot::plot_grid(plot1, plot2, ncol = 1, nrow = 2)

graphics.off()



#......................
# get deaths posterior pred check
#......................
postdeaths <- COVIDCurve::posterior_check_infxns_to_death(IFRmodel_inf = modout,
                                                          dwnsmpl = 1e2,
                                                          by_chain = FALSE)
postdeaths.plotObj <- postdeaths %>%
  dplyr::select(c("time", dplyr::starts_with("deaths"))) %>%
  tidyr::gather(., key = "strata", value = "deaths", 2:ncol(.)) %>%
  ggplot() +
  geom_line(aes(x= time, y = deaths, group = strata, color = strata), size = 1.2) +
  scale_color_viridis_d()

jpeg("~/Desktop/posterior_check.jpg", width = 11, height = 8, units = "in", res = 500)
postdeaths.plotObj +
  geom_line(data = dat$AggDat, aes(x=ObsDay, y = Deaths, group = Strata), color = "#bdbdbd", size = 0.75) +
  theme_bw() +
  ggtitle("Posterior Predictive Check", subtitle = "Grey Lines are Simulated Data, Viridis Lines are Draws from Posterior")
graphics.off()






#............................................................
# look at rungs
#...........................................................
COVIDCurve::get_cred_intervals(IFRmodel = mod1,
                               mcmcout = modout,
                               whichrung = "rung50",
                               what = "IFRparams",
                               by_chain = F)

modout$output %>%
  dplyr::mutate(
    logpost = loglikelihood + logprior
  ) %>%
  dplyr::filter(stage == "sampling" &
                  rung == "rung50") %>%
  ggplot() +
  geom_histogram(aes(x  = logpost))

rungdf <- modout$output %>%
  dplyr::mutate(
    logpost = loglikelihood + logprior
  ) %>%
  dplyr::filter(stage == "sampling" &
                  rung == "rung50")

minloglike <- min(rungdf$logpost)
# is this neg infinity
if (minloglike < -1.796e306 & minloglike > -1.8e306) {
  sum(rungdf$logpost == minloglike)/nrow(rungdf)
}

#............................................................
# are the problem ones the rungs that are inf?
#...........................................................
rungdf %>%
  dplyr::filter(logpost > -1.8e306) %>%
  dplyr::select_at(c("chain", "r1", "r2", "ma3")) %>%
  tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
  dplyr::group_by_at("param") %>%
  dplyr::summarise(
    min = min(est),
    LCI = quantile(est, 0.025),
    median = median(est),
    mean = mean(est),
    UCI = quantile(est, 0.975),
    max = max(est),
    ESS = coda::effectiveSize(coda::as.mcmc(est)),
    GewekeZ = coda::geweke.diag(coda::as.mcmc(est))[[1]],
    GewekeP = dnorm(GewekeZ)
  )











#......................
# cumulative
#......................
dat <- COVIDCurve::Aggsim_infxn_2_death(
  casefat = casefat,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 150,
  level = "Cumulative",
  infections = infxns$infxns
)

mod1 <- make_modinf_agg$new()
mod1$set_level("Cumulative")
mod1$set_data(dat)
mod1$set_IFRparams(c("r1", "r2", "ma3"))
mod1$set_Infxnparams(c("y1", "y2", "y3", "y4", "y5", "y6"))
mod1$set_paramdf(df_params)
mod1$set_pa(c(1/3, 1/3, 1/3))
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_knots(c(1, 30, 60, 90, 120, 150))
r_mcmc_out.cumm <- COVIDCurve::run_modinf_agg(modinf = mod1, reparamIFR = T, rungs = 10)
plot_par(r_mcmc_out.cumm, "y1")



# sanity
