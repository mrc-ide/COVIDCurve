context("mod fit with cpp background")
test_that("serology likelihood accurate", {
  # shape of the curve is essentially the same, but pin is wrong for less likely scenario

  set.seed(1234)
  library(drjacoby)


  #............................................................
  # simulate data
  #...........................................................
  # sigmoidal function
  infxns <- data.frame(time = 1:150)
  sig <- function(x){1 / (1 +  exp(-x))}
  timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
  infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                              min = -25,
                                              max = 50)

  # make up fatality data
  fatalitydata <- data.frame(Strata = c("ma1", "ma2", "ma3"),
                             IFR = c(0.05, 0.2, 0.5),
                             Rho = c(1, 1, 1),
                             Ne = c(1, 1, 1))
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
    infections = infxns$infxns,
    simulate_seroprevalence = TRUE,
    sens = 0.85,
    spec = 0.95,
    sero_delay_rate = 10
  )

  obs_serology <- dat$SeroPrev %>%
    dplyr::group_by(Strata) %>%
    dplyr::filter(event_obs_day %in% sero_days) %>%
    dplyr::rename(
      SeroDay = event_obs_day,
      SeroPrev = ObsPrev) %>%
    dplyr::select(c("SeroDay", "Strata", "SeroPrev")) %>%
    dplyr::mutate(SeroDay = ifelse(SeroDay == 135, "sero_day1", "sero_day2")) %>%
    dplyr::arrange(SeroDay) %>%
    dplyr::ungroup(.)

  datinput <- list(obs_deaths = dat$AggDeath,
                   obs_serology = obs_serology)

  #..................
  # make model
  #..................
  ifr_paramsdf <- tibble::tibble(name = c("ma1", "ma2",  "ma3"),
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
                                  min =   c(0.83,     0.8,    0,           135,        160),
                                  init =  c(0.85,     0.95,   0.9,         135,        160),
                                  max =   c(0.87,     1.00,   1,           135,        160),
                                  dsc1 =  c(8500,     950,    90,          130,        150),
                                  dsc2 =  c(1500,     50,     10,          140,        170))

  noise_paramsdf <- tibble::tibble(name = c("ne1", "ne2", "ne3"),
                                   min  = rep(1, 3),
                                   init = rep(1, 3),
                                   max = rep(1, 3),
                                   dsc1 = rep(0, 3),
                                   dsc2 = rep(10, 3))

  # onset to deaths
  tod_paramsdf <- tibble::tibble(name = c("mod", "sod"),
                                 min  = c(10,     0.01),
                                 init = c(14,     0.7),
                                 max =  c(20,     1.00),
                                 dsc1 = c(2.657,  -0.236),
                                 dsc2 = c(0.01,   0.01))


  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sero_paramsdf, noise_paramsdf, tod_paramsdf)

  #......................
  # make mode
  #......................
  mod1 <- make_IFRmodel_agg$new()
  mod1$set_MeanTODparam("mod")
  mod1$set_CoefVarOnsetTODparam("sod")
  mod1$set_IFRparams(c("ma1", "ma2", "ma3"))
  mod1$set_maxMa("ma3")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y5")
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams(c("sero_day1", "sero_day2"))
  mod1$set_Noiseparams(c("ne1", "ne2", "ne3"))
  mod1$set_data(datinput)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rho(c(1/3, 1/3, 1/3))
  mod1$set_rcensor_day(.Machine$integer.max)

  #............................................................
  # run MCMC
  #...........................................................
  modout <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod1,
                                         reparamIFR = TRUE,
                                         reparamInfxn = TRUE,
                                         reparamKnot = TRUE,
                                         reparamSeroRate = TRUE,
                                         burnin = 1e2,
                                         samples = 1e2,
                                         chains = 1,
                                         silent = FALSE)



  aggstring_cpp_function_wrapper <- function(params, infmodel) {
    #......................
    # inputs needed for cpp function from sim data
    #......................
    misc_list = list(rho = infmodel$inputs$IFRmodel$rho,
                     demog = infmodel$inputs$IFRmodel$demog$popN,
                     rcensor_day = infmodel$inputs$IFRmodel$rcensor_day,
                     days_obsd = infmodel$inputs$IFRmodel$maxObsDay,
                     n_knots = length(infmodel$inputs$IFRmodel$Knotparams)+1,
                     n_sero_obs = length(infmodel$inputs$IFRmodel$Serodayparams))

    datin <- list("obs_deaths" = infmodel$inputs$IFRmodel$data$obs_deaths$Deaths,
                  "obs_serology" = infmodel$inputs$IFRmodel$data$obs_serology$SeroPrev)

    #......................
    # internall Cpp function controlling likelihood
    #......................
    fitcurve_string <- COVIDCurve:::make_user_Agg_loglike(IFRmodel = infmodel$inputs$IFRmodel,
                                                          reparamIFR = FALSE,
                                                          reparamKnots = FALSE,
                                                          reparamInfxn = FALSE,
                                                          reparamSeroRate = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
    # pull out pieces I need
    fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "double loglik = -OVERFLO_DOUBLE;", n = 2)[,1]
    fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)

    fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
    fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "loglik = -OVERFLO_DOUBLE;", n = 2)[,1]
    fitcurve_curve <- sub("if \\(cum_infxn_check <= popN\\) \\{", "", fitcurve_curve)
    fitcurve_string <- paste(fitcurve_start, "double loglik = -OVERFLO_DOUBLE;",
                             fitcurve_curve,
                             "loglik = -OVERFLO_DOUBLE; }",
                             "Rcpp::List ret = Rcpp::List::create(loglik, sero_con_num, death_loglik, sero_loglik, infxn_spline); return ret;}",
                             collapse = "")
    Rcpp::cppFunction(fitcurve_string)

    paramsin <- unlist(params[c(infmodel$inputs$IFRmodel$modparam,
                                infmodel$inputs$IFRmodel$sodparam,
                                infmodel$inputs$IFRmodel$IFRparams,
                                infmodel$inputs$IFRmodel$Infxnparams,
                                infmodel$inputs$IFRmodel$Knotparams,
                                infmodel$inputs$IFRmodel$Serotestparams,
                                infmodel$inputs$IFRmodel$Serodayparams,
                                infmodel$inputs$IFRmodel$Noiseparams)])

    full_mod_out <- loglike(params = paramsin,
                            param_i = 1,
                            data = datin,
                            misc = misc_list)
    names(full_mod_out) <- c("loglik", "sero_con_num", "death_loglik", "sero_loglik", "infxn_spline")

    return(full_mod_out)

  }
  #......................
  # get params for max likelihood
  #......................
  maxparams <- modout$mcmcout$output[modout$mcmcout$output$loglikelihood == max(modout$mcmcout$output$loglikelihood), ]
  out <- aggstring_cpp_function_wrapper(params = maxparams, infmodel = modout)

  #......................
  # test equal
  #......................
  testthat::expect_equal(out$loglik, max(modout$mcmcout$output$loglikelihood))



})
