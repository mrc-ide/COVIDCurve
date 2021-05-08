context("posterior summaries working as expected")

test_that("check get cred intervals function as dataframes", {
  library(COVIDCurve)
  data("covidcurve_modfit")
  covidcurve_modfit
  ifr <- COVIDCurve::get_cred_intervals(IFRmodel_inf = covidcurve_modfit, whichrung = paste0("rung", 1),
                                        what = "IFRparams", by_chain = F)
  sero <- COVIDCurve::get_cred_intervals(IFRmodel_inf = covidcurve_modfit, whichrung = paste0("rung", 1),
                                         what = "Serotestparams", by_chain = F)
  knotspost <- COVIDCurve::get_cred_intervals(IFRmodel_inf = covidcurve_modfit,  whichrung = paste0("rung", 1),
                                              what = "Knotparams", by_chain = F)
  infxn <- COVIDCurve::get_cred_intervals(IFRmodel_inf = covidcurve_modfit,  whichrung = paste0("rung", 1),
                                          what = "Infxnparams", by_chain = F)
  ddelays <- COVIDCurve::get_cred_intervals(IFRmodel_inf = covidcurve_modfit,  whichrung = paste0("rung", 1),
                                            what = "DeathDelayparams", by_chain = F)
  neparams <- COVIDCurve::get_cred_intervals(IFRmodel_inf = covidcurve_modfit,  whichrung = paste0("rung", 1),
                                             what = "Noiseparams", by_chain = F)

  # tests
  testthat::expect_is(ifr, "data.frame")
  testthat::expect_is(sero, "data.frame")
  testthat::expect_is(knotspost, "data.frame")
  testthat::expect_is(infxn, "data.frame")
  testthat::expect_is(ddelays, "data.frame")
  testthat::expect_is(neparams, "data.frame")

})
