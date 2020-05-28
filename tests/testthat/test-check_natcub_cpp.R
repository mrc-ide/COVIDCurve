context("aggregate fits for cpp")
test_that("natcub cpp likelihood works", {
  set.seed(48)
  #..................
  # sim data
  #..................
  # sigmoidal function
  infxns <- data.frame(time = 1:199)
  sig <- function(x){1 / (1 +  exp(-x))}
  timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
  infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                              min = -25,
                                              max = 50)
  pa <- c(0.5, 0.5)
  casefat <- data.frame(age = c("ma1", "ma2"),
                        ifr = c(0.1, 0.5),
                        pa = pa)
  dat <- COVIDCurve::Aggsim_infxn_2_death(
    casefat = casefat,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 200,
    level = "Time-Series",
    infections = infxns$infxns,
    simulate_seroprevalence = TRUE,
    sero_spec = 0.95,
    sero_sens = 0.8,
    sero_delay_rate = 10,
    popN = 5e5
  )

  #..................
  # inputs
  #..................
  # params in
  # misc list
  knots <- c(1, 40, 80, 120, 160, 200)
  day <- knots[1]:knots[length(knots)]
  gamma_lookup <- stats::pgamma((day-1),
                                shape = 1/0.45^2, scale = 18.8*0.45^2)

  misc_list = list(pa = pa,
                   pgmms = gamma_lookup,
                   knots = knots,
                   level = FALSE,
                   popN = 5e5)

  # liftover to Rcpp list

  morelikely.paramsin <- c("r1" = 0.2, "ma2" = 0.5, "y1" = 3.95,
                           "y2" = 5.73, "y3" = 7.70, "y4" = 8.42,
                           "y5" = 8.51, "y6" = 8.52,
                           "sens" = 0.8, "spec" = 0.95, "sero_rate" = 10, "sero_day" = 35.1)
  sero_day <- 35
  datin <- list("obs_deaths" = dat$AggDat$Deaths,
                "obs_serologyrate" = dat$seroprev$SeroRateFP[sero_day])
  morelikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike(params = morelikely.paramsin,
                                                           param_i = 1,
                                                           data = datin,
                                                           misc = misc_list)

  lesslikely.paramsin <- c("r1" = 0.1, "ma2" = 0.5, "y1" = 4, "y2" = 4, "y3" = 4, "y4" = 4, "y5" = 4, "y6" = 4,
                           "sens" = 0.95, "spec" = 0.95, "sero_rate" = 0.95, "sero_day" = 45.6)
  lesslikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike(params = lesslikely.paramsin,
                                                           param_i = 1,
                                                           data = datin,
                                                           misc = misc_list)

  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

})
