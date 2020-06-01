context("aggregate fits for cpp")
test_that("natcub cpp likelihood works", {
  set.seed(48)
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
                             ifr = c(0.05, 0.1, 0.5),
                             pa = pa)
  dat <- COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 150,
    level = "Time-Series",
    infections = infxns$infxns,
    simulate_seroprevalence = TRUE,
    specificity = 0.95,
    sensitivity = 0.8,
    sero_delay_rate = 10,
    popN = 5e5
  )

  #..................
  # inputs
  #..................
  # params in
  # misc list
  knots <- c(1, 30, 60, 90, 120, 150)
  day <- knots[1]:(knots[length(knots)]+1)
  gamma_lookup <- stats::pgamma((day-1),
                                shape = 1/0.45^2, scale = 18.8*0.45^2)

  misc_list = list(pa = pa,
                   pgmms = gamma_lookup,
                   knots = knots,
                   level = FALSE,
                   popN = 5e5)

  # liftover to Rcpp list

  morelikely.paramsin <- c("r1" = 0.1, "r2" = 0.2, "ma3" = 0.5,
                           "y1" = 3.95, "y2" = 5.83, "y3" = 7.70, "y4" = 8.41, "y5" = 8.50, "y6" = 8.52,
                           "sens" = 0.8, "spec" = 0.95, "sero_rate" = 10, "sero_day" = 135.1)
  sero_day <- 135
  datin <- list("obs_deaths" = dat$AggDat$Deaths,
                "obs_serologyrate" = dat$seroprev$SeroRateFP[sero_day])

  # truth
  morelikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike_cubicspline(params = morelikely.paramsin,
                                                                       param_i = 1,
                                                                       data = datin,
                                                                       misc = misc_list)

  # random
  lesslikely.paramsin <- c("r1" = 0.5, "r2" = 0.5, "ma3" = 0.5,
                           "y1" = 4, "y2" = 4, "y3" = 4, "y4" = 4, "y5" = 4, "y6" = 4,
                           "sens" = 0.95, "spec" = 0.95, "sero_rate" = 0.95, "sero_day" = 45.6)
  lesslikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike_cubicspline(params = lesslikely.paramsin,
                                                                       param_i = 1,
                                                                       data = datin,
                                                                       misc = misc_list)

  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

})
