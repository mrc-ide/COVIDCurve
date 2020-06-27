context("mod fit with cpp background")
test_that("serology likelihood accurate", {
  # shape of the curve is essentially the same, but pin is wrong for less likely scenario

  set.seed(1234)
  library(drjacoby)
  # sigmoidal function
  infxns <- data.frame(time = 1:150)
  sig <- function(x){1 / (1 +  exp(-x))}
  timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
  infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                              min = -25,
                                              max = 50)
  sum(infxns$infxns < 0)

  # make up fatality data
  fatalitydata <- data.frame(Strata = c("ma1", "ma2", "ma3"),
                             IFR = c(0.05, 0.2, 0.5),
                             Rho = c(0.4, 0.4, 0.2),
                             Ne = c(0.1, 0.4, 0.5))
  # make up demography data
  demog <- data.frame(Strata = c("ma1", "ma2", "ma3"),
                      popN = c(1500000, 2250000, 1250000))
  # pick serology date
  sero_days <- c(110, 135)

  #..................
  # run sim
  #..................
  dat <- COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    demog = demog,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 150,
    level = "Time-Series",
    infections = infxns$infxns,
    simulate_seroprevalence = TRUE,
    sens = 0.85,
    spec = 0.99,
    sero_delay_rate = 10
  )

  #......................
  # input data
  #......................
  obs_serology <- dat$seroprev %>%
    dplyr::filter(event_obs_day %in% sero_days) %>%
    dplyr::arrange(event_obs_day, Strata)

  datinput <- list(obs_deaths = dat$AggDat$Deaths,
                   obs_serology = obs_serology$ObsPrev)

  #..................
  # inputs
  #..................
  # misc list
  days_obsd <- 150
  knots <- c(30, 60, 90, 120)
  day <- 1:(days_obsd+1)
  gamma_lookup <- stats::pgamma((day-1),
                                shape = 1/0.45^2, scale = 18.8*0.45^2)

  misc_list = list(rho = fatalitydata$Rho,
                   pgmms = gamma_lookup,
                   level = FALSE,
                   popN = 5e6,
                   rcensor_day = .Machine$integer.max,
                   days_obsd = days_obsd,
                   n_sero_obs = 2,
                   n_knots = length(knots)+1,
                   demog = demog$popN)

  # liftover to Rcpp list

  morelikely.paramsin <- c("r1" = 0.05, "r2" = 0.2, "ma3" = 0.5,
                           "x1" = 30, "x2" = 60, "x3" = 90, "x4" = 120,
                           "y1" = 2.8, "y2" = 5.7, "y3" = 7.7, "y4" = 8.4, "y5" = 8.5,
                           "ne1" = 0.1, "ne2" = 0.45,
                           "sens" = 0.85, "spec" = 0.99, "sero_rate" = 10,
                           "sero_day1" = 110, "sero_day2" = 135)


  # truth
  morelikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike_cubicspline(params = morelikely.paramsin,
                                                                       param_i = 1,
                                                                       data = datinput,
                                                                       misc = misc_list)

  # random
  lesslikely.paramsin <- c("r1" = 0.1, "r2" = 0.39, "ma3" = 0.98,
                           "x1" = 22.25, "x2" = 54.47, "x3" = 109.9, "x4" = 145.58,
                           "y1" = 2.98, "y2" = 4.52, "y3" = 6.74, "y4" = 7.82, "y5" = 7.88,
                           "ne1" = 0.1, "ne2" = 0.45,
                           "sens" = 0.85, "spec" = 0.99, "sero_rate" = 10,
                           "sero_day1" = 110, "sero_day2" = 135)

  lesslikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike_cubicspline(params = lesslikely.paramsin,
                                                                       param_i = 1,
                                                                       data = datinput,
                                                                       misc = misc_list)


  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

})
