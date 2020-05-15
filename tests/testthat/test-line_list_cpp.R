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
  dat <- COVIDCurve::LineListsim_infxn_2_death(
    casefat = casefat,
    m_od = 18.8,
    s_od = 0.45,
    m_or = 24.7,
    s_or = 0.35,
    hospprob = 0.005,
    curr_day = 50,
    infections = infxns$infxns,
    expgrowth = F) %>%
    dplyr::filter(hosp == 1) %>%
    dplyr::select(-c("hosp"))

  #..................
  # inputs
  #..................
  misc_list = list(IFRparams = c(1, 2))

  #..................
  # make data list
  #..................
  death_onset_day <- dat %>%
    dplyr::filter(Outcome == "Death") %>%
    dplyr::select(c("OnsetDay")) %>%
    unlist(.) %>%
    unname(.)

  death_event_day <- dat %>%
    dplyr::filter(Outcome == "Death") %>%
    dplyr::select(c("EventDay")) %>%
    unlist(.) %>%
    unname(.)

  death_group <- dat %>%
    dplyr::filter(Outcome == "Death") %>%
    dplyr::select(c("AgeGroup")) %>%
    unlist(.) %>%
    unname(.)
  death_group <- factor(death_group)

  recovery_onset_day <- dat %>%
    dplyr::filter(Outcome == "Recovery") %>%
    dplyr::select(c("OnsetDay")) %>%
    unlist(.) %>%
    unname(.)

  recovery_event_day <- dat %>%
    dplyr::filter(Outcome == "Recovery") %>%
    dplyr::select(c("EventDay")) %>%
    unlist(.) %>%
    unname(.)

  recovery_group <- dat %>%
    dplyr::filter(Outcome == "Recovery") %>%
    dplyr::select(c("AgeGroup")) %>%
    unlist(.) %>%
    unname(.)
  recovery_group <- factor(recovery_group)

  datin <- list(death_interval = death_event_day - death_onset_day,
                death_group = as.numeric(death_group),
                recovery_interval = recovery_event_day -recovery_onset_day,
                recovery_group = as.numeric(recovery_group))



  # liftover to Rcpp list

  morelikely.paramsin <- c("ma1" = 0.1, "ma2" = 0.5, "mod" = 18.8, "sod" = 0.45, "mor" = 24.7, "sor" = 0.35)

  morelikely <- COVIDCurve:::LineList_loglike(params = morelikely.paramsin,
                                              param_i = 1,
                                              data = datin,
                                              misc = misc_list)

  lesslikely.paramsin <- c("ma1" = 0.3, "ma2" = 0.8, "mod" = 13.8, "sod" = 1.45, "mor" = 14.7, "sor" = 5.35)
  lesslikely <- COVIDCurve:::LineList_loglike(params = lesslikely.paramsin,
                                              param_i = 1,
                                              data = datin,
                                              misc = misc_list)

  testthat::expect_gt(object = morelikely$LogLik, expected = lesslikely$LogLik)

})
