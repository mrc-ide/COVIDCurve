context("R function that writes cpp is working")
test_that("R wrapper works", {

  set.seed(1234)
  library(drjacoby)
  # sigmoidal function
  infxns <- data.frame(time = 1:200)
  sig <- function(x){1 / (1 +  exp(-x))}
  timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
  infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                              min = -25,
                                              max = 50)
  sum(infxns$infxns < 0)

  # make up fatality data
  fatalitydata <- data.frame(strata = c("ma1", "ma2", "ma3"),
                             ifr = c(0.05, 0.2, 0.5),
                             rho = 1/3)
  # pick serology date
  sero_days <- c(135, 160)

  #..................
  # run sim
  #..................
  dat <- COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 200,
    level = "Time-Series",
    infections = infxns$infxns,
    simulate_seroprevalence = TRUE,
    sens = 0.85,
    spec = 0.99,
    sero_delay_rate = 10,
    popN = 5e6
  )


  datinput <- list(obs_deaths = dat$AggDat,
                   obs_serologyrate = dat$seroprev$ObsPrev[sero_days])

  #..................
  # make model
  #..................
  ifr_paramsdf <- tibble::tibble(name = c("r1", "r2",  "ma3"),
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
  sero_paramsdf <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1", "sero_day2"),
                                  min =   c(0.83,     0.97,   10,         130,         150),
                                  init =  c(0.85,     0.99,   10,         135,         160),
                                  max =   c(0.87,     1.00,   10,         140,         170),
                                  dsc1 =  c(8500,     990,    5,          130,         150),
                                  dsc2 =  c(1500,     10,     15,         140,         170))

  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sero_paramsdf)

  #......................
  # make mode
  #......................
  mod1 <- make_IFRmodel_agg$new()
  mod1$set_MeanOnset(18.8)
  mod1$set_CoefVarOnset(0.45)
  mod1$set_level("Time-Series")
  mod1$set_data(datinput)
  mod1$set_IFRparams(c("r1", "r2", "ma3"))
  mod1$set_maxMa("ma3")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y5")
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams(c("sero_day1", "sero_day2"))
  mod1$set_popN(5e6)
  mod1$set_paramdf(df_params)
  mod1$set_rho(c(1/3, 1/3, 1/3))
  mod1$set_rcensor_day(.Machine$integer.max)

  #..................
  # run model
  #..................
  modout <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod1,
                                         reparamIFR = TRUE,
                                         reparamInfxn = TRUE,
                                         reparamKnot = TRUE,
                                         burnin = 1e1,
                                         samples = 1e1,
                                         chains = 1)

  testthat::expect_is(modout, "IFRmodel_inf")
  testthat::expect_is(modout$mcmcout, "drjacoby_output")

})
