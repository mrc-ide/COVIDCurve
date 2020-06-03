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
  pa <- rep(1/3, 3)
  fatalitydata <- data.frame(strata = c("r1", "r2", "ma3"),
                             ifr = c(0.05, 0.2, 0.5),
                             pa = pa)
  dat <- COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 150,
    level = "Time-Series",
    infections = infxns$infxns,
    simulate_seroprevalence = TRUE,
    specificity = 1.0,
    sensitivity = 1.0,
    sero_delay_rate = 10,
    popN = 5e5
  )

  #..................
  # inputs
  #..................
  # params in
  # misc list
  days_obsd <- 150
  rcensor_day <- .Machine$integer.max
  knots <- c(1, 30, 60, 90, 120, 150)
  day <- 1:(days_obsd+1)
  gamma_lookup <- stats::pgamma((day-1),
                                shape = 1/0.45^2, scale = 18.8*0.45^2)

  misc_list = list(pa = pa,
                   pgmms = gamma_lookup,
                   level = FALSE,
                   popN = 5e5,
                   rcensor_day = rcensor_day,
                   days_obsd = days_obsd,
                   n_knots = length(knots)+1)

  # liftover to Rcpp list

  morelikely.paramsin <- c("r1" = 0.05, "r2" = 0.2, "ma3" = 0.5,
                           "x1" = 30, "x2" = 60, "x3" = 90, "x4" = 120,  "x5" = 150,
                           "y1" = 17, "y2" = 304, "y3" = 2230, "y4" = 4529, "y5" = 4984,  "y6" = 5014,
                           "sens" = 1.0, "spec" = 1.0, "sero_rate" = 10, "sero_day" = 135.1)
  sero_day <- 135
  datin <- list("obs_deaths" = dat$AggDat$Deaths,
                "obs_serologyrate" = dat$seroprev$SeroRateFP[sero_day])

  # truth
  morelikely <- COVIDCurve:::playgroundNatCubic_SplineGrowth_loglike_cubicspline(params = morelikely.paramsin,
                                                                                 param_i = 1,
                                                                                 data = datin,
                                                                                 misc = misc_list)
  morelikely$LogLik
  morelikely$auc
  morelikely$sp1
  morelikely$sp2
  morelikely$sp3

  plot(morelikely$infxn_spline)

  # random
  lesslikely.paramsin <- c("r1" = 0.5, "r2" = 0.5, "ma3" = 0.5,
                           "x1" = 1, "x2" = 30, "x3" = 60, "x4" = 90, "x5" = 120,  "x6" = 150,
                           "y1" = 50, "y2" = 150, "y3" = 1000, "y4" = 5000, "y5" = 5000, "y6" = 5000,
                           "sens" = 0.95, "spec" = 0.95, "sero_rate" = 5, "sero_day" = 120.6)
  lesslikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike_cubicspline(params = lesslikely.paramsin,
                                                                       param_i = 1,
                                                                       data = datin,
                                                                       misc = misc_list)

  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

})
