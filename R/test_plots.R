context("plotting works as expected")

test_that("infection plot works", {
  data(covidcurve_modfit)
  curve <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf = covidcurve_modfit,
                                                          whichrung = paste0("rung", 1),
                                                          by_chain = F,
                                                          by_strata = F,
                                                          dwnsmpl = 1e2)

  testthat::expect_is(curve$plotObj, "ggplot")
  testthat::expect_is(curve$curvedata, "data.frame")
})


test_that("seroprevalence dataframe for plot works", {
  data(covidcurve_modfit)
  serocurve <- COVIDCurve::draw_posterior_sero_curves(IFRmodel_inf = covidcurve_modfit,
                                                      whichrung = paste0("rung", 1),
                                                      by_chain = F,
                                                      dwnsmpl = 1e2)
  testthat::expect_is(serocurve, "data.frame")
})


test_that("posterior deaths dataframe for plot works", {
  data(covidcurve_modfit)
  postdeaths <- COVIDCurve::posterior_check_infxns_to_death(IFRmodel_inf = covidcurve_modfit,
                                                            dwnsmpl = 1e2,
                                                            by_chain = FALSE)
  testthat::expect_is(postdeaths, "data.frame")

})
