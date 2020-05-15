test_that("cpp likelihood works", {
  set.seed(48)
  #..................
  # sim data
  #..................
  # sigmoidal function
  infxns <- data.frame(time = 1:49)
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
    curr_day = 50,
    level = "Time-Series",
    expgrowth = F,
    infections = infxns$infxns
  )
  #..................
  # inputs
  #..................
  # params in

  # misc list
  knots <- c(1, 10, 20, 30, 40, 50)
  day <- knots[1]:knots[length(knots)]
  gamma_lookup <- stats::pgamma((day-1),
                                shape = 1/0.45^2, scale = 18.8*0.45^2)

  misc_list = list(pa = pa,
                   pgmms = gamma_lookup,
                   knots = knots,
                   level = FALSE)

  # liftover to Rcpp list

  morelikely.paramsin <- c("r1" = 0.2, "ma2" = 0.5, "y1" = 3.95, "y2" = 5.60, "y3" = 7.54, "y4" = 8.41, "y5" = 8.51, "y6" = 8.52)
  datin <- list("obs_deaths" = dat$Deaths)
  morelikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike(params = morelikely.paramsin,
                                                           param_i = 1,
                                                           data = datin,
                                                           misc = misc_list)

  lesslikely.paramsin <- c("r1" = 0.1, "ma2" = 0.5, "y1" = 4, "y2" = 4, "y3" = 4, "y4" = 4, "y5" = 4, "y6" = 4)
  lesslikely <- COVIDCurve:::NatCubic_SplineGrowth_loglike(params = lesslikely.paramsin,
                                                           param_i = 1,
                                                           data = datin,
                                                           misc = misc_list)

  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

})
