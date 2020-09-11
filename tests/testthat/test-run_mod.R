context("IFR model runs")
test_that("R wrapper works", {
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
  #sero_days <- c(150)
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
    m_od = 19.2,
    s_od = 0.79,
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
                                 init = rep(0.5, 3),
                                 max = rep(1, 3),
                                 dsc1 = rep(0, 3),
                                 dsc2 = rep(1, 3))

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
                                  init =  c(0.85,     0.90),
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
                                 init = c(19,     0.79,  18),
                                 max =  c(20,     1,     21),
                                 dsc1 = c(19.26,  2370,  18.3),
                                 dsc2 = c(0.1,    630,   0.1))

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

  #..................
  # run model
  #..................
  modout <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod1,
                                         reparamIFR = TRUE,
                                         reparamInfxn = TRUE,
                                         reparamKnot = TRUE,
                                         burnin = 1e2,
                                         samples = 1e2,
                                         chains = 1,
                                         rungs = 1,
                                         thinning = 0,
                                         silent = FALSE)

  testthat::expect_is(modout, "IFRmodel_inf")
  testthat::expect_is(modout$mcmcout, "drjacoby_output")

})
