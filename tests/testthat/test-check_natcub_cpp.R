context("mod fit with cpp background")

test_that("likelihood accurate", {
  set.seed(1234)
  # sigmoidal function
  infxns <- data.frame(time = 1:200)
  sig <- function(x){1 / (1 +  exp(-x))}
  timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
  infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                              min = -25,
                                              max = 50)
  sum(infxns$infxns < 0)

  # make up fatality data
  fatalitydata <- tibble::tibble(Strata = c("ma1", "ma2", "ma3"),
                                 IFR = c(0.01, 0.05, 0.1),
                                 Rho = 1)
  demog <- tibble::tibble(Strata = c("ma1", "ma2", "ma3"),
                          popN = c(1500000, 2250000, 1250000))

  # pick serology date
  sero_days <- c(135, 160)

  #..............................................................
  # AGGREGATE
  #..............................................................
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
    dplyr::mutate(SeroStartSurvey = c(130, 155),
                  SeroEndSurvey = c(140, 165)) %>%
    dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev")) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(SeroStartSurvey, Strata)

  #......................
  # binomial
  #......................
  datinput_binomial <- list(obs_deaths = dat$Agg_TimeSeries_Death$Deaths,
                            prop_strata_obs_deaths = prop_strata_obs_deaths$PropDeaths,
                            obs_serologypos = obs_serology$SeroPos,
                            obs_serologyn = obs_serology$SeroN)
  #......................
  # logit
  #......................
  obs_serology <- obs_serology %>%
    dplyr::mutate(SeroMu = log(SeroPrev/(1-SeroPrev)),
                  SeroSE = 0.05) # pick some SE
  datinput_logit <- list(obs_deaths = dat$Agg_TimeSeries_Death$Deaths,
                         prop_strata_obs_deaths = prop_strata_obs_deaths$PropDeaths,
                         obs_serologymu = obs_serology$SeroMu,
                         obs_serologyse = obs_serology$SeroSE)

  #..................
  # inputs
  #..................
  # misc list
  knots <- c(30, 60, 90, 120)
  misc_list = list(rcensor_day = .Machine$integer.max,
                   days_obsd = max(dat$Agg_TimeSeries_Death$ObsDay),
                   n_sero_obs = 2,
                   max_seroday_obsd = 140,
                   sero_survey_start = c(105, 130),
                   sero_survey_end = c(115, 140),
                   n_knots = length(knots)+1,
                   demog = demog$popN,
                   account_serorev = FALSE)


  # liftover to Rcpp list
  morelikely.paramsin <- c("mod" = 19.8, "sod" = 0.85,
                           "ma1" = 0.01, "ma2" = 0.05, "ma3" = 0.1,
                           "x1" = 30, "x2" = 60, "x3" = 90, "x4" = 120,
                           "y1" = 2.8, "y2" = 5.7, "y3" = 7.7, "y4" = 8.4, "y5" = 8.5,
                           "ne1" = 0.33, "ne2" = 0.33, "ne3" = 0.33,
                           "sens" = 0.85, "spec" = 0.95, "sero_con_rate" = 18,
                           "sero_rev_shape" = NA, "sero_rev_scale" = NA) # just to make visible to cpp

  #......................
  # binomial
  #......................
  # near truth
  morelikely <- COVIDCurve:::natcubspline_loglike_binomial(params = morelikely.paramsin,
                                                           param_i = 1,
                                                           data = datinput_binomial,
                                                           misc = misc_list)


  # random
  lesslikely.paramsin <-  c("mod" = 19.8, "sod" = 0.85,
                            "ma1" = 0.1, "ma2" = 0.39, "ma3" = 0.98,
                            "x1" = 22.25, "x2" = 54.47, "x3" = 109.9, "x4" = 145.58,
                            "y1" = 2.98, "y2" = 4.52, "y3" = 6.74, "y4" = 7.82, "y5" = 7.88,
                            "ne1" = 0.1, "ne2" = 0.4, "ne3" = 0.5,
                            "sens" = 0.85, "spec" = 0.99, "sero_con_rate" = 10,
                            "sero_rev_shape" = NA, "sero_rev_scale" = NA) # just to make visible to cpp


  lesslikely <- COVIDCurve:::natcubspline_loglike_binomial(params = lesslikely.paramsin,
                                                           param_i = 1,
                                                           data = datinput_binomial,
                                                           misc = misc_list)
  # check
  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

  #......................
  # logit
  #......................
  # near truth
  morelikely <- COVIDCurve:::natcubspline_loglike_logit(params = morelikely.paramsin,
                                                        param_i = 1,
                                                        data = datinput_logit,
                                                        misc = misc_list)


  lesslikely <- COVIDCurve:::natcubspline_loglike_logit(params = lesslikely.paramsin,
                                                        param_i = 1,
                                                        data = datinput_logit,
                                                        misc = misc_list)
  # check
  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)



})
