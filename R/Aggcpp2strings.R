#' @title Create the logprior from the IFRmodel R6 Class for DrJacoby Inference
#' @param IFRmodel R6 class; Internal model object for COVIDCurve
#' @param reparamIFR logical; Whether IFRs should be reparameterized or inferred seperately
#' @param reparamKnots logical; Whether infection knots (i.e. the x-coordinates of the infection spline) should be reparameterized or inferred seperately
#' @param reparamInfxn logical; Whether infection curve (i.e. the  y-coordinates infection spline) should be reparameterized or inferred seperately
#' @noRd

make_user_Agg_logprior <- function(IFRmodel, reparamIFR, reparamInfxn, reparamKnots) {
  #..................
  # assertsions
  #..................
  assert_custom_class(IFRmodel, "IFRmodel")
  #..................
  # setup
  #..................
  paramdf <- IFRmodel$paramdf
  IFRparams <- paramdf[paramdf$name %in% IFRmodel$IFRparams, ]
  Knotparams <- paramdf[paramdf$name %in% IFRmodel$Knotparams, ]
  Infxnparams <- paramdf[paramdf$name %in% IFRmodel$Infxnparams, ]
  Seroparams <- paramdf[paramdf$name %in% IFRmodel$Seroparams, ]

  if (reparamKnots) {
    #..................
    # account for knot reparam -- xpos
    #..................
    assert_non_null(IFRmodel$relKnot, message = "Reparameterization requires relative knot to be indicated (i.e. relKnot)")
    relKnot <- IFRmodel$relKnot
    knotscalars <- Knotparams$name[Knotparams$name != relKnot]
  }

  if (reparamInfxn) {
    #..................
    # account for Infection Point reparam -- Ypos
    #..................
    assert_non_null(IFRmodel$relInfxn, message = "Reparameterization requires relative infection point to be indicated (i.e. relInfxn)")
    relInfxn <- IFRmodel$relInfxn
    infxnscalars <- Infxnparams$name[Infxnparams$name != relInfxn]
  }

  if (reparamIFR) {
    #..................
    # account for IFR reparam
    #..................
    assert_non_null(IFRmodel$maxMa, message = "Reparameterization requires expected max mortaltiy group to be indicated (i.e. maxMa)")
    maxMa <- IFRmodel$maxMa
    ifrscalars <- IFRparams$name[IFRparams$name != maxMa]
  }

  #..................
  # priors for knots -- Xpos
  #..................
  Knotextractparams <- sapply(Knotparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeknotpriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(", param, ",", d1, ",", d2, ",", "true) +")
  }, param = Knotparams$name, d1 = Knotparams$dsc1, d2 = Knotparams$dsc2)

  #..................
  # priors for infxnpts -- Ypos
  #..................
  Infxnextractparams <- sapply(Infxnparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeinfxnpriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(", param, ",", d1, ",", d2, ",", "true) +")
  }, param = Infxnparams$name, d1 = Infxnparams$dsc1, d2 = Infxnparams$dsc2)

  #..................
  # priors for IFRparams
  #..................
  IFRextractparams <- sapply(IFRparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeifrpriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = IFRparams$name, d1 = IFRparams$dsc1, d2 = IFRparams$dsc2)

  #......................
  # Serology priors
  #......................
  Seroextractparams <- sapply(Seroparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })
  makeSerobetapriors <- mapply(function(param, d1, d2){
    paste0("R::dbeta(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = Seroparams$name[!Seroparams$name %in% c("sero_rate", "sero_day")],
  d1 = Seroparams$dsc1[!Seroparams$name %in% c("sero_rate", "sero_day")],
  d2 = Seroparams$dsc2[!Seroparams$name %in% c("sero_rate", "sero_day")])
  serorateprior <- paste0("R::dunif(sero_rate,", Seroparams$dsc1[Seroparams$name == "sero_rate"], ",", Seroparams$dsc2[Seroparams$name == "sero_rate"], ",true) +")
  serodateprior <- paste0("R::dunif(sero_day,", Seroparams$dsc1[Seroparams$name == "sero_day"], ",", Seroparams$dsc2[Seroparams$name == "sero_day"], ",true) +")

  #..................
  # bring together
  #..................
  extractparams <- c(IFRextractparams, Knotextractparams, Infxnextractparams, Seroextractparams)

  switch(paste0(reparamIFR, "-", reparamKnots, "-", reparamInfxn),
         "TRUE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerobetapriors, serorateprior, serodateprior,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ");"))
         },

         "TRUE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerobetapriors, serorateprior, serodateprior,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ");"))
         },

         "TRUE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerobetapriors, serorateprior, serodateprior,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ");"))
         },

         "FALSE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerobetapriors, serorateprior, serodateprior,
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ");"))
         },

         "TRUE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerobetapriors, serorateprior, serodateprior,
                       paste0(length(ifrscalars), "*log(", maxMa, ");"))
         },

         "FALSE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerobetapriors, serorateprior, serodateprior,
                       paste0(length(knotscalars), "*log(", relKnot, ");"))
         },

         "FALSE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerobetapriors, serorateprior, serodateprior,
                       paste0(length(infxnscalars), "*log(", relInfxn, ");"))
         },

         "FALSE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerobetapriors, serorateprior, serodateprior)
           priors[length(priors)] <- sub("\\) \\+$", ");", priors[length(priors)]) # trailing + sign to a semicolon
         },

         {
           stop("Prior option not correctly specified during make priors")
         }
  )

  ret <- c("SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {",
           extractparams,
           priors,
           "if (!std::isfinite(ret)) { const double OVERFLO_DOUBLE = DBL_MAX/100.0; ret = -OVERFLO_DOUBLE; }",
           "return Rcpp::wrap(ret);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}



#' @title Create the loglikelihood from the IFRmodel R6 Class for DrJacoby Inference
#' @inheritParams make_user_Agg_logprior
#' @noRd

make_user_Agg_loglike <- function(IFRmodel, reparamIFR, reparamInfxn, reparamKnots) {
  #..................
  # assertsions
  #..................
  assert_custom_class(IFRmodel, "IFRmodel")
  #..................
  # setup
  #..................
  paramdf <- IFRmodel$paramdf
  Knotparams <- IFRmodel$Knotparams
  Infxnparams <- IFRmodel$Infxnparams
  IFRparams <- IFRmodel$IFRparams

  #..................
  # extract misc
  #..................
  extmisc <- "std::vector<double> pa = Rcpp::as< std::vector<double> >(misc[\"pa\"]); std::vector<double> pgmms = Rcpp::as< std::vector<double> >(misc[\"pgmms\"]); bool level = misc[\"level\"]; int popN = misc[\"popN\"]; int rcensor_day = misc[\"rcensor_day\"]; int days_obsd = misc[\"days_obsd\"]; int n_knots = misc[\"n_knots\"];"

  #..................
  # extract inputs
  #..................
  params <- sapply(paramdf$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })
  # convert sero_date
  params <- gsub("double sero_day", "double sero_day_raw", params)
  params <- c(params, "int sero_day = std::floor(sero_day_raw);")


  #..................
  # storage items
  #..................
  storageitems <- "int agelen = pa.size(); std::vector<double>ma(agelen); std::vector<double> node_x_raw(n_knots); std::vector<double> node_x(n_knots); std::vector<double> node_y(n_knots);"

  #..................
  # liftover knotreparam vars for Knots -- Infxn Xpositions
  # NB, "raw" here because we take in a double and need to convert it to an integer day later
  #.................
  if (reparamKnots) {
    assert_non_null(IFRmodel$relKnot, message = "Reparameterization requires relative knot to be indicated (i.e. relKnot)")
    knotparamdf <- paramdf[paramdf$name %in% Knotparams, ]
    relKnot <- IFRmodel$relKnot
    knotscalars <- knotparamdf$name[knotparamdf$name != relKnot]
    node_xvec <- rep(NA, length(Knotparams))
    # determine which relative position in our knot vector
    relnodex_pos <- which(relKnot == Knotparams)
    nodex_counter <- 1
    for (i in 1:length(node_xvec)) {
      if (i == relnodex_pos) {
        node_xvec[i] <- paste0("node_x_raw[", i, "] = ", relKnot, ";")
      } else {
        node_xvec[i] <- paste0("node_x_raw[", i, "] = ", knotscalars[nodex_counter], "*", relKnot, ";")
        nodex_counter <- nodex_counter + 1
      }
    }
  } else {
    node_xvec <- rep(NA, length(Knotparams))
    for (i in 1:length(Knotparams)){
      node_xvec[i] <- paste0("node_x_raw[", i, "]", " = ", Knotparams[i], ";")
    }
  }
  # account for internal knot at position 1
  node_xvec <- c("node_x_raw[0] = 1.0;", node_xvec)

  #..................
  # liftover infxnreparam vars to Infxn Ypositions
  #.................
  if (reparamInfxn) {
    assert_non_null(IFRmodel$relInfxn, message = "Reparameterization requires relative infection points to be indicated (i.e. relInfxn)")
    infxnparamdf <- paramdf[paramdf$name %in% Infxnparams, ]
    relInfxn <- IFRmodel$relInfxn
    infxnscalars <- infxnparamdf$name[infxnparamdf$name != relInfxn]
    node_yvec <- rep(NA, length(Infxnparams))
    # determine which relative position in our node-y vector
    relnodey_pos <- which(relInfxn == Infxnparams)
    nodey_counter <- 1
    for (i in 1:length(node_yvec)) {
      if (i == relnodey_pos) {
        node_yvec[i] <- paste0("node_y[", i-1, "] = ", relInfxn, ";")
      } else {
        node_yvec[i] <- paste0("node_y[", i-1, "] = ", infxnscalars[nodey_counter], "*", relInfxn, ";")
        nodey_counter <- nodey_counter + 1
      }
    }
  } else {
    node_yvec <- rep(NA, length(Infxnparams))
    for (i in 1:length(Infxnparams)){
      node_yvec[i] <- paste0("node_y[", i-1, "]", " = ", Infxnparams[i], ";")
    }
  }

  #..................
  # liftover ifrreparam vars to Mas
  #.................
  if (reparamIFR) {
    assert_non_null(IFRmodel$maxMa, message = "Reparameterization requires expected max mortaltiy group to be indicated (i.e. maxMa)")
    ifrparamdf <- paramdf[paramdf$name %in% IFRparams, ]
    maxMa <- IFRmodel$maxMa
    ifrscalars <- ifrparamdf$name[ifrparamdf$name != maxMa]
    mavec <- rep(NA, length(IFRparams))
    # determine which relative position in our ma vector
    mamax_pos <- which(maxMa == IFRparams)
    ma_counter <- 1
    for (i in 1:length(mavec)) {
      if (i == mamax_pos) {
        mavec[i] <- paste0("ma[", i-1, "] = ", maxMa, ";")
      } else {
        mavec[i] <- paste0("ma[", i-1, "] = ", ifrscalars[ma_counter], "*", maxMa, ";")
        ma_counter <- ma_counter + 1
      }
    }
  } else {
    mavec <- rep(NA, length(IFRparams))
    for (i in 1:length(IFRparams)){
      mavec[i] <- paste0("ma[", i-1, "]", " = ", IFRparams[i], ";")
    }
  }

  #..................
  # discretize knots (i.e. day is a discrete time)
  #..................
  node_xvec.discretize <- "for (int i = 0; i < node_x.size(); i++) { node_x[i] = std::ceil(node_x_raw[i]); }"

  #..................
  # get loglike
  #..................
  # TODO temp for debugging otherwise won't work for devtools::load_all() quickly
  #loglike <- readLines("~/Documents/GitHub/COVIDCurve/src/NatCubic_AggExpDeaths_loglike_cubicspline.cpp")
  #loglike <- readLines(system.file("src/NatCubic_AggExpDeaths_loglike_cubicspline.cpp", package = "COVIDCurve", mustWork = TRUE))
  #loglike_start <- grep("// Deaths Section", loglike)
  #loglike_end <- grep("// return as Rcpp list", loglike)
  #loglike <- loglike[loglike_start:loglike_end]
  loglike <- "const double OVERFLO_DOUBLE = DBL_MAX/100.0;
  double loglik = -OVERFLO_DOUBLE;
  bool nodex_pass = true;
  for (int i = 1; i < node_x.size(); i++) {
    if (node_x[i] <= node_x[i-1]) {
      nodex_pass = false;
    }
  }
  if (nodex_pass) {
    std::vector<double> denom(n_knots-1);
    for (int i = 0; i < denom.size(); i++) {
      denom[i] = node_x[i+1] - node_x[i];
    }
    std::vector<double> z(n_knots-1);
    z[0] = 0;
    std::vector<double> m(n_knots-2);
    for (int i = 1; i <= m.size(); i++) {
      m[i-1] = (3/denom[i])*(node_y[i+1] - node_y[i]) - (3/denom[i-1])*(node_y[i] - node_y[i-1]);
    }
    std::vector<double> g(n_knots-1);
    std::vector<double> k(n_knots-1);
    g[0] = 1;
    k[0] = 0;
    for (int i = 1; i < (n_knots-1); i++) {
      g[i] = 2*(node_x[i+1] - node_x[i-1]) - (denom[i-1])*(k[i-1]);
      k[i] = denom[i]/g[i];
      z[i] = (m[i-1] - denom[i-1]*z[i-1])/g[i];
    }
    std::vector<double> sp1(n_knots-1);
    std::vector<double> sp3(n_knots-1);
    std::vector<double> sp2(n_knots);
    sp2[n_knots-1] = 0;
    for (int i = (n_knots-2); i >= 0; i--) {
      sp2[i] = z[i] - k[i]*sp2[i+1];
      sp1[i] = (node_y[i+1] - node_y[i])/(denom[i]) - (denom[i]*(sp2[i+1] + 2*sp2[i]))/3 ;
      sp3[i] = (sp2[i+1] - sp2[i])/(3*denom[i]);
    }
    std::vector<double> infxn_spline(days_obsd);
    int node_j = 0;
    infxn_spline[0] = node_y[0];
    for (int i = 1; i < days_obsd; i++) {
      infxn_spline[i] = node_y[node_j] +
        sp1[node_j] * ((i+1) - node_x[node_j]) +
        sp2[node_j] * pow(((i+1) - node_x[node_j]), 2) +
        sp3[node_j] * pow(((i+1) - node_x[node_j]), 3);
      if (node_j < (node_x.size()-2)) {
        if ((node_x[0] + i) >= node_x[node_j+1]) {
          node_j++;
        }
      }
    }
    for (int i = 0; i < infxn_spline.size(); i++) {
      infxn_spline[i] = exp(infxn_spline[i]);
    }
    std::vector<double> cumm_infxn_spline(infxn_spline.size());
    cumm_infxn_spline[0] = infxn_spline[0];
    for (int i = 1; i < cumm_infxn_spline.size(); i++) {
      cumm_infxn_spline[i] = infxn_spline[i] + cumm_infxn_spline[i-1];
    }
    if (cumm_infxn_spline[cumm_infxn_spline.size() - 1] <= popN) {
      std::vector<double> auc(infxn_spline.size());
      for (int i = 0; i < infxn_spline.size(); i++) {
        for (int j = i+1; j < (infxn_spline.size() + 1); j++) {
          int delta = j - i - 1;
          auc[j-1] += infxn_spline[i] * (pgmms[delta + 1] - pgmms[delta]);
        }
      }
      double death_loglik = 0.0;
      if (level) {
        std::vector<int> obsd = Rcpp::as< std::vector<int> >(data[\"obs_deaths\"]);
        double aucsum = 0;
        for (int i = 0; i < auc.size(); i++) {
          aucsum += auc[i];
        }
        std::vector<double>expd(agelen);
        for (int a = 0; a < agelen; a++) {
          expd[a] = aucsum * pa[a] * ma[a];
        }
        for (int a = 0; a < agelen; a++) {
          if ((a+1) < rcensor_day) {
            if (obsd[a] != -1) {
              death_loglik += R::dpois(obsd[a], expd[a], true);
            }
          }
        }
      } else {
        std::vector<int> raw = Rcpp::as< std::vector<int> >(data[\"obs_deaths\"]);
        std::vector<std::vector<int>> obsd(infxn_spline.size(), std::vector<int>(agelen));
        int iter = 0;
        for (int i = 0; i < infxn_spline.size(); i++) {
          for (int j = 0; j < agelen; j++) {
            obsd[i][j] = raw[iter];
            iter++;
          }
        }
        std::vector<std::vector<double>> expd(infxn_spline.size(), std::vector<double>(agelen));
        for (int  i = 0; i < infxn_spline.size(); i++) {
          for (int a = 0; a < agelen; a++) {
            expd[i][a] = auc[i] * pa[a] * ma[a];
          }
        }
        for (int  i = 0; i < infxn_spline.size(); i++) {
          for (int a = 0; a < agelen; a++) {
            if ((a+1) < rcensor_day) {
              if (obsd[i][a] != -1) {
                death_loglik += R::dpois(obsd[i][a], expd[i][a], true);
              }
            }
          }
        }
      }
      std::vector<double> cum_hazard(sero_day);
      cum_hazard[0] = 0.0;
      for (int i = 1; i < sero_day; i++) {
        cum_hazard[i] = (1-exp((-i/sero_rate)));
      }
      double sero_con_num = 0.0;
      for (int i = 0; i < sero_day; i++) {
        int time_elapsed = sero_day - i - 1;
        sero_con_num += infxn_spline[i] * cum_hazard[time_elapsed];
      }
      double obs_prev = (sero_con_num/popN) * (spec + (sens-1)) - (spec-1);
      int posint = round(obs_prev * popN);
      double datpos = data[\"obs_serologyrate\"];
      double sero_loglik = R::dbinom(posint, popN, datpos, true);
      loglik = death_loglik + sero_loglik;
      if (!std::isfinite(loglik)) {
        loglik = -OVERFLO_DOUBLE;
      }
    }
  }
"
  # remove comments which cause issues when coercing to string in this format
  #commlines <- grep("//", loglike)
  #loglike <- loglike[! 1:length(loglike) %in% commlines]
  loglike <- gsub("\\\n", "", loglike)
  # end loglike and now out
  ret <- c("SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {",
           extmisc,
           params,
           storageitems,
           node_xvec,
           node_xvec.discretize,
           node_yvec,
           mavec,
           loglike,
           "return Rcpp::wrap(loglik);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}


