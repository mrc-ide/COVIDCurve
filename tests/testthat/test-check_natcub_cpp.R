context("aggregate fits for cpp")
test_that("natcub cpp likelihood works", {

  #............................................................
  # N.B. By default, we are not reparameterizing
  #...........................................................

  set.seed(1234)
  #..................
  # sim data
  #..................
  # sigmoidal function
  infxns <- data.frame(time = 1:150)
  sig <- function(x){1 / (1 +  exp(-x))}
  timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
  infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                              min = -25,
                                              max = 50)
  # make up fatality data
  fatalitydata <- data.frame(strata = c("r1", "r2", "ma3"),
                             ifr = c(0.05, 0.2, 0.5),
                             pa = 1/3)
  # pick serology date
  sero_day <- 135

  simdat <- COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 150,
    level = "Time-Series",
    infections = infxns$infxns,
    simulate_seroprevalence = TRUE,
    spec = 0.95,
    sens = 0.80,
    sero_delay_rate = 10,
    popN = 5e6
  )

  #..................
  # inputs
  #..................
  # params in
  # misc list
  days_obsd <- 150
  rcensor_day <- 130
  knots <- c(30, 60, 90, 120)
  day <- 1:(days_obsd+1)
  gamma_lookup <- stats::pgamma((day-1),
                                shape = 1/0.45^2, scale = 18.8*0.45^2)

  misc_list = list(pa = fatalitydata$pa,
                   pgmms = gamma_lookup,
                   level = FALSE,
                   popN = 5e6,
                   rcensor_day = rcensor_day,
                   days_obsd = days_obsd,
                   n_knots = length(knots)+1)

  # liftover to Rcpp list

  morelikely.paramsin <- c("r1" = 0.05, "r2" = 0.2, "ma3" = 0.5,
                           "x1" = 30, "x2" = 60, "x3" = 90, "x4" = 120,
                           "y1" = 2.8, "y2" = 5.7, "y3" = 7.7, "y4" = 8.4, "y5" = 8.5,
                           "sens" = 0.8, "spec" = 0.95, "sero_rate" = 10, "sero_day" = 135)
  datin <- list("obs_deaths" = simdat$AggDat$Deaths,
                "obs_serologyrate" = simdat$seroprev$SeroRateFP[sero_day])

  # truth
  morelikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike_cubicspline(params = morelikely.paramsin,
                                                                                 param_i = 1,
                                                                                 data = datin,
                                                                                 misc = misc_list)

  # random
  lesslikely.paramsindup <- c("r1" = 0.1, "r2" = 0.39, "ma3" = 0.98,
                           "x1" = 22.25, "x2" = 54.47, "x3" = 109.9, "x4" = 145.58,
                           "y1" = 2.98, "y2" = 4.52, "y3" = 6.74, "y4" = 7.82, "y5" = 7.88,
                           "sens" = 0.8, "spec" = 0.95, "sero_rate" = 10, "sero_day" = 125)
  lesslikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike_cubicspline(params = lesslikely.paramsin,
                                                                       param_i = 1,
                                                                       data = datin,
                                                                       misc = misc_list)


  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

})
