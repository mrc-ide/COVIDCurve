#' @title Run Aggregate Model
#' @details Wraps the Metropolic-Coupled MCMC Framework from Dr. Jacoby
#' @inheritParams drjacoby::run_mcmc
#' @export

run_IFRmodel_agg <- function(IFRmodel, reparamIFR = TRUE, reparamInfxn = TRUE, reparamKnots = TRUE,
                             burnin = 1e3, samples = 1e3, chains = 3,
                             rungs = 1, GTI_pow = 3, coupling_on = TRUE,
                             pb_markdown = FALSE, silent = TRUE) {
  #..................
  # assertions
  #..................
  assert_custom_class(IFRmodel, "IFRmodel")
  assert_logical(reparamIFR)
  assert_logical(reparamInfxn)
  assert_logical(reparamKnots)
  assert_numeric(burnin)
  assert_numeric(samples)
  assert_numeric(chains)
  assert_numeric(rungs)
  assert_numeric(GTI_pow)
  assert_logical(coupling_on)
  assert_logical(pb_markdown)
  assert_logical(silent)
  assert_non_null(IFRmodel$level)
  assert_non_null(IFRmodel$data)
  assert_non_null(IFRmodel$IFRparams)
  assert_non_null(IFRmodel$Infxnparams)
  assert_non_null(IFRmodel$Knotparams)
  assert_non_null(IFRmodel$paramdf)
  assert_non_null(IFRmodel$pa)
  assert_non_null(IFRmodel$Seroparams)
  assert_non_null(IFRmodel$popN)
  assert_non_null(IFRmodel$mod)
  assert_non_null(IFRmodel$sod)
  assert_non_null(IFRmodel$gamma_lookup)
  assert_non_null(IFRmodel$maxObsDay)

  #............................................................
  # "Warm-Up" MCMC
  #...........................................................
  warmdf_params <- rbind.data.frame(list("x", 1, 1, 1))
  names(warmdf_params) <- c("name", "min", "max", "init")
  warmloglike <- "SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) { double ret = -1.0; return Rcpp::wrap(ret);}"
  warmlogprior <- "SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) { double ret = -1.0; return Rcpp::wrap(ret);}"
  warmup <- drjacoby::run_mcmc(data = list("dat" = c(1)),
                               df_params = warmdf_params,
                               misc = list(),
                               loglike = warmloglike,
                               logprior = warmlogprior,
                               burnin = 1,
                               samples = 1,
                               chains = 1,
                               rungs = 1,
                               silent = T)


  #..............................................................
  # unpack object
  #..............................................................

  #..................
  # Get loglike and logprior
  #..................
  if (reparamIFR) {
    assert_non_null(IFRmodel$maxMa, message = "If performing reparameterization, must set a maximum Ma in the R6 class object")
  }

  if (reparamInfxn) {
    assert_non_null(IFRmodel$relInfxn, message = "If performing reparameterization, must set a relative infection point in the R6 class object")
  }

  if (reparamKnots) {
    assert_non_null(IFRmodel$relKnot, message = "If performing reparameterization, must set a relative knot point in the R6 class object")
  }

  logpriorfunc <- COVIDCurve:::make_user_Agg_logprior(IFRmodel, reparamIFR = reparamIFR, reparamInfxn = reparamInfxn, reparamKnots = reparamKnots)
  loglikfunc <- COVIDCurve:::make_user_Agg_loglike(IFRmodel, reparamIFR = reparamIFR, reparamInfxn = reparamInfxn, reparamKnots = reparamKnots)

  #..................
  # make misc
  #..................
  misc_list = list(pa = IFRmodel$pa,
                   pgmms = IFRmodel$gamma_lookup,
                   level = ifelse(IFRmodel$level == "Cumulative", TRUE, FALSE),
                   popN = IFRmodel$popN,
                   rcensor_day = IFRmodel$rcensor_day,
                   days_obsd = IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel$Knotparams) + 1) # +1 because we set an internal knot for pos 1

  #..................
  # make data list
  #..................
  if (IFRmodel$level == "Time-Series"){
    data_list <- split(IFRmodel$data$obs_deaths$Deaths, factor(IFRmodel$data$obs_deaths$ObsDay))
    data_list <- unname(unlist(data_list))
    data_list <- list(obs_deaths = data_list,
                      obs_serologyrate = IFRmodel$data$obs_serologyrate)

  } else if (IFRmodel$level == "Cumulative") {
    data_list <- list(obs_deaths = unname(IFRmodel$data$obs_deaths$Deaths),
                      obs_serologyrate = IFRmodel$data$obs_serologyrate)
  }

  #..................
  # make df param
  #..................
  df_params <-  IFRmodel$paramdf[, 1:4]

  #..............................................................
  # Dr Jacoby
  #..............................................................

  mcmcout <- drjacoby::run_mcmc(data = data_list,
                                df_params = df_params,
                                misc = misc_list,
                                loglike = loglikfunc,
                                logprior = logpriorfunc,
                                burnin = burnin,
                                samples = samples,
                                chains = chains,
                                rungs = rungs,
                                coupling_on = coupling_on,
                                GTI_pow = GTI_pow,
                                pb_markdown = pb_markdown,
                                silent = silent
  )

  if (reparamIFR) {
    #..................
    # account for ifr reparam
    #..................
    IFRparams <- IFRmodel$paramdf[IFRmodel$paramdf$name %in% IFRmodel$IFRparams, ]
    maxMa <- IFRmodel$maxMa
    scalars <- IFRparams$name[IFRparams$name != maxMa]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, maxMa]})
  }

  if (reparamKnots) {
    #..................
    # account for knots (infxn X position) reparam
    #..................
    Knotparams <- IFRmodel$paramdf[IFRmodel$paramdf$name %in% IFRmodel$Knotparams, ]
    relKnot <- IFRmodel$relKnot
    scalars <- Knotparams$name[Knotparams$name != relKnot]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, relKnot]})
  }

  if (reparamInfxn) {
    #..................
    # account for Infxn Y position reparam
    #..................
    Infxnparams <- IFRmodel$paramdf[IFRmodel$paramdf$name %in% IFRmodel$Infxnparams, ]
    relInfxn <- IFRmodel$relInfxn
    scalars <- Infxnparams$name[Infxnparams$name != relInfxn]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, relInfxn]})
  }


  # store input along with Dr.Jacoby output for later use
  inputs <- list(
    IFRmodel = IFRmodel,
    reparamIFR = reparamIFR,
    reparamInfxn = reparamInfxn,
    reparamKnots = reparamKnots,
    burnin = burnin,
    samples = samples,
    chains = chains)
  if (rungs > 1) {
    inputs <- append(inputs, list(rungs = rungs,
                                  GTI_pow = GTI_pow,
                                  coupling_on = coupling_on))
  }

  ret <- list(
    inputs = inputs,
    mcmcout = mcmcout
  )
  class(ret) <- c("IFRmodel_inf")
  return(ret)
}



