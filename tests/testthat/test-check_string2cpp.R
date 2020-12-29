context("mod fit with cpp background")
test_that("R writes Cpp function correctly", {
  set.seed(1234)

  #............................................................
  # simulate data
  #...........................................................
  # sigmoidal function
  infxns <- data.frame(time = 1:200)
  sig <- function(x){1 / (1 +  exp(-x))}
  timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
  infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                              min = -25,
                                              max = 50)

  # make up fatality data
  fatalitydata <- data.frame(Strata = c("ma1", "ma2", "ma3"),
                             IFR = c(0.05, 0.2, 0.5),
                             Rho = c(1, 1, 1))
  # make up demography data
  demog <- data.frame(Strata = c("ma1", "ma2", "ma3"),
                      popN = c(1500000, 2250000, 1250000))
  # pick serology date
  sero_days <- c(110, 135)

  #..................
  # run sim
  #..................
  dat <- COVIDCurve::Agesim_infxn_2_death(
    fatalitydata = fatalitydata,
    demog = demog,
    m_od = 19.8,
    s_od = 0.85,
    curr_day = 200,
    infections = infxns$infxns,
    simulate_seroreversion = FALSE,
    sens = 0.85,
    spec = 0.95,
    sero_delay_rate = 18.3
  )


  # liftover proprtion deaths
  totdeaths <- sum(dat$StrataAgg_TimeSeries_Death$Deaths)
  prop_strata_obs_deaths <- dat$StrataAgg_TimeSeries_Death %>%
    dplyr::group_by(Strata) %>%
    dplyr::summarise(Deaths = sum(Deaths),
                     PropDeaths = Deaths/totdeaths) %>%
    dplyr::select(c("Strata", "PropDeaths"))

  # liftover obs serology
  obs_serology <- dat$StrataAgg_Seroprev %>%
    dplyr::group_by(Strata) %>%
    dplyr::filter(ObsDay %in% sero_days) %>%
    dplyr::mutate(
      SeroPos = round(ObsPrev * testedN),
      SeroN = testedN ) %>%
    dplyr::rename(
      SeroPrev = ObsPrev) %>%
    dplyr::mutate(SeroStartSurvey = sero_days - 5,
                  SeroEndSurvey = sero_days + 5) %>%
    dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev")) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(SeroStartSurvey, Strata) %>%
    dplyr::mutate(SeroLCI = SeroPrev - 0.01,
                  SeroUCI = SeroPrev + 0.01) # make up some tight CIs

  datinput <- list(obs_deaths = dat$Agg_TimeSeries_Death,
                   prop_deaths = prop_strata_obs_deaths,
                   obs_serology = obs_serology)


  #..................
  # make model
  #..................
  ifr_paramsdf <- tibble::tibble(name = c("ma1", "ma2",  "ma3"),
                                 min  = rep(0, 3),
                                 init = rep(0.1, 3),
                                 max = rep(0.4, 3),
                                 dsc1 = rep(0, 3),
                                 dsc2 = rep(0.4, 3))

  infxn_paramsdf <- tibble::tibble(name = paste0("y", 1:5),
                                   min  = rep(0, 5),
                                   init = c(rep(0.5, 4), 5),
                                   max =  c(rep(1, 4), 10),
                                   dsc1 = rep(0, 5),
                                   dsc2 = c(rep(1, 4), 10))

  knot_paramsdf <- tibble::tibble(name = paste0("x", 1:4),
                                  min  = c(0,    0.33, 0.66, 175),
                                  init = c(0.05, 0.40, 0.75, 185),
                                  max =  c(0.33, 0.66, 0.99, 200),
                                  dsc1 = c(0,    0.33, 0.66, 175),
                                  dsc2 = c(0.33, 0.66, 0.99, 200))
  sero_paramsdf <- tibble::tibble(name =  c("sens", "spec"),
                                  min =   c(0.83,     0.8),
                                  init =  c(0.85,     0.95),
                                  max =   c(0.87,     1.00),
                                  dsc1 =  c(8500,     950),
                                  dsc2 =  c(1500,     50))

  noise_paramsdf <- tibble::tibble(name = c("ne1", "ne2", "ne3"),
                                   min  = rep(0.5, 3),
                                   init = rep(1, 3),
                                   max = rep(1.5, 3),
                                   dsc1 = rep(1, 3),
                                   dsc2 = rep(0.05, 3))

  # onset to deaths
  tod_paramsdf <- tibble::tibble(name = c("mod", "sod", "sero_con_rate"),
                                 min  = c(18,     0,     16),
                                 init = c(19,     0.8,   18),
                                 max =  c(20,     1,     21),
                                 dsc1 = c(19.8,  2550,  18.3),
                                 dsc2 = c(0.1,    450,   0.1))



  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sero_paramsdf, noise_paramsdf, tod_paramsdf)

  #......................
  # make mode
  #......................
  mod1 <- make_IFRmodel_age$new()
  mod1$set_MeanTODparam("mod")
  mod1$set_CoefVarOnsetTODparam("sod")
  mod1$set_IFRparams(c("ma1", "ma2", "ma3"))
  mod1$set_maxMa("ma3")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y5")
  mod1$set_Serotestparams(c("sens", "spec", "sero_con_rate"))
  mod1$set_Noiseparams(c("ne1", "ne2", "ne3"))
  mod1$set_data(datinput)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rcensor_day(.Machine$integer.max)

  #............................................................
  # run MCMC
  #...........................................................
  modout <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod1,
                                         reparamIFR = TRUE,
                                         reparamInfxn = TRUE,
                                         reparamKnot = TRUE,
                                         burnin = 1e3,
                                         samples = 1e3,
                                         chains = 1,
                                         silent = FALSE)


  # internal cpp function
  agestring_cpp_function_wrapper <- function(params, infmodel) {
    #......................
    # inputs needed for cpp function from sim data
    #......................
    misc_list = list(rcensor_day = .Machine$integer.max,
                     days_obsd = infmodel$inputs$IFRmodel$maxObsDay,
                     n_knots = length(infmodel$inputs$IFRmodel$Knotparams) + 1, # +1 because we set an internal knot for pos 1
                     n_sero_obs = length(unique(infmodel$inputs$IFRmodel$data$obs_serology$SeroStartSurvey)),
                     sero_survey_start = unique(infmodel$inputs$IFRmodel$data$obs_serology$SeroStartSurvey),
                     sero_survey_end = unique(infmodel$inputs$IFRmodel$data$obs_serology$SeroEndSurvey),
                     max_seroday_obsd = max(infmodel$inputs$IFRmodel$data$obs_serology$SeroEndSurvey),
                     demog = infmodel$inputs$IFRmodel$demog$popN,
                     account_serorev = infmodel$inputs$account_seroreversion)

    datin <- list(obs_deaths = dat$Agg_TimeSeries_Death$Deaths,
                  prop_strata_obs_deaths = prop_strata_obs_deaths$PropDeaths,
                  obs_serologypos = obs_serology$SeroPos,
                  obs_serologyn = obs_serology$SeroN)

    #......................
    # internall Cpp function controlling likelihood
    #......................
    fitcurve_string <- COVIDCurve:::make_user_Age_loglike(IFRmodel = infmodel$inputs$IFRmodel,
                                                          account_serorev = infmodel$inputs$account_seroreversion,
                                                          binomial_likelihood = infmodel$inputs$binomial_likelihood,
                                                          reparamIFR = FALSE,
                                                          reparamKnots = FALSE,
                                                          reparamInfxn = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
    # pull out pieces I need
    fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "double loglik = -OVERFLO_DOUBLE;", n = 2)[,1]
    fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
    fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
    fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "loglik = -OVERFLO_DOUBLE;", n = 2)[,1]
    fitcurve_curve <- stringr::str_replace(fitcurve_curve, "  if \\(popN_pass\\) \\{", "")
    fitcurve_string <- paste(fitcurve_start, "double loglik = -OVERFLO_DOUBLE;",
                             fitcurve_curve,
                             "loglik = -OVERFLO_DOUBLE; }",
                             "Rcpp::List ret = Rcpp::List::create(loglik); return ret;}",
                             collapse = "")
    Rcpp::cppFunction(fitcurve_string)

    paramsin <- unlist(params[c(infmodel$inputs$IFRmodel$modparam,
                                infmodel$inputs$IFRmodel$sodparam,
                                infmodel$inputs$IFRmodel$IFRparams,
                                infmodel$inputs$IFRmodel$Infxnparams,
                                infmodel$inputs$IFRmodel$Knotparams,
                                infmodel$inputs$IFRmodel$Serotestparams,
                                infmodel$inputs$IFRmodel$Noiseparams)])

    full_mod_out <- loglike(params = paramsin,
                            param_i = 1,
                            data = datin,
                            misc = misc_list)
    names(full_mod_out) <- c("loglik")

    return(full_mod_out)

  }
  #......................
  # get params for max likelihood
  #......................
  maxparams <- modout$mcmcout$output[modout$mcmcout$output$loglikelihood == max(modout$mcmcout$output$loglikelihood), ]
  out <- agestring_cpp_function_wrapper(params = maxparams, infmodel = modout)
  #......................
  # test equal
  #......................
  testthat::expect_equal(out$loglik, max(modout$mcmcout$output$loglikelihood))

})
