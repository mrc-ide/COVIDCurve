context("mod fit with cpp background")

test_that("serology likelihood accurate", {
  set.seed(1234)
  library(magrittr)
  library(tidyverse)
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
  dat <- COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    demog = demog,
    m_od = 19.2,
    s_od = 0.79,
    curr_day = 200,
    infections = infxns$infxns,
    simulate_seroreversion = TRUE,
    sero_rev_shape = 4.75,
    sero_rev_scale = 272,
    sens = 0.85,
    spec = 0.95,
    sero_delay_rate = 15
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

  datinput <- list(obs_deaths = dat$Agg_TimeSeries_Death$Deaths,
                   prop_strata_obs_deaths = prop_strata_obs_deaths$PropDeaths,
                   obs_serologypos = obs_serology$SeroPos,
                   obs_serologyn = obs_serology$SeroN)

  #..................
  # inputs
  #..................
  # misc list
  days_obsd <- 150
  knots <- c(30, 60, 90, 120)
  misc_list = list(rho = fatalitydata$Rho,
                   rcensor_day = .Machine$integer.max,
                   days_obsd = days_obsd,
                   n_sero_obs = 2,
                   max_seroday_obsd = 140,
                   sero_survey_start = c(105, 130),
                   sero_survey_end = c(115, 140),
                   n_knots = length(knots)+1,
                   demog = demog$popN,
                   account_serorev = FALSE)


  # liftover to Rcpp list
  morelikely.paramsin <- c("mod" = 17.8, "sod" = 0.45,
                           "ma1" = 0.01, "ma2" = 0.05, "ma3" = 0.1,
                           "x1" = 30, "x2" = 60, "x3" = 90, "x4" = 120,
                           "y1" = 2.8, "y2" = 5.7, "y3" = 7.7, "y4" = 8.4, "y5" = 8.5,
                           "ne1" = 0.33, "ne2" = 0.33, "ne3" = 0.33,
                           "sens" = 0.85, "spec" = 0.95, "sero_rate" = 10,
                           "sero_rev_scale" = 272, "sero_rev_shape" = 4.75)

  # truth
  morelikely <- COVIDCurve:::natcubspline_loglike(params = morelikely.paramsin,
                                                  param_i = 1,
                                                  data = datinput,
                                                  misc = misc_list)

  morelikely


  # random
  lesslikely.paramsin <-  c("mod" = 17.8, "sod" = 0.45,
                            "ma1" = 0.1, "ma2" = 0.39, "ma3" = 0.98,
                            "x1" = 22.25, "x2" = 54.47, "x3" = 109.9, "x4" = 145.58,
                            "y1" = 2.98, "y2" = 4.52, "y3" = 6.74, "y4" = 7.82, "y5" = 7.88,
                            "ne1" = 0.1, "ne2" = 0.4, "ne3" = 0.5,
                            "sens" = 0.85, "spec" = 0.99, "sero_rate" = 10,
                            "sero_day1" = 110, "sero_day2" = 135)

  lesslikely <- COVIDCurve:::natcubspline_loglike(params = lesslikely.paramsin,
                                                  param_i = 1,
                                                  data = datinput,
                                                  misc = misc_list)


  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

})
