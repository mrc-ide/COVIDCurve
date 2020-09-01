#' @title Run Aggregate Model
#' @details Wraps the Metropolic-Coupled MCMC Framework from Dr. Jacoby
#' @inheritParams drjacoby::run_mcmc
#' @param IFRmodel R6 class; Internal model object for COVIDCurve
#' @param reparamIFR logical; Whether IFRs should be reparameterized or inferred seperately
#' @param reparamKnots logical; Whether infection knots (i.e. the x-coordinates of the infection spline) should be reparameterized or inferred seperately
#' @param reparamInfxn logical; Whether infection curve (i.e. the  y-coordinates infection spline) should be reparameterized or inferred seperately
#' @param reparamDelays logical; Whether the mean offset-to-death and the seroconversion rate delay should be reparameterized or inferred seperately
#' @param reparamNe logical; Whether "noise scalar effects" should be reparameterized or inferred seperately (if TRUE, considered relateve to Ne1)
#' @param thinning integer; The regular sequence to count by to thin MCMC posterior chain (iterations are kept as: \code{seq(from = thinning, to = (burnin+samples), by = thinning)}).
#' @export

run_IFRmodel_agg <- function(IFRmodel, reparamIFR = TRUE, reparamInfxn = TRUE, reparamKnots = TRUE, reparamDelays = TRUE, reparamNe = TRUE,
                             burnin = 1e3, samples = 1e3, chains = 3, thinning = 0,
                             rungs = 1, GTI_pow = 3, coupling_on = TRUE,
                             cluster = NULL, pb_markdown = FALSE, silent = TRUE) {
  #..................
  # assertions
  #..................
  assert_custom_class(IFRmodel, "IFRmodel")
  assert_logical(reparamIFR)
  assert_logical(reparamInfxn)
  assert_logical(reparamKnots)
  assert_logical(reparamDelays)
  assert_logical(reparamNe)
  assert_numeric(burnin)
  assert_numeric(samples)
  assert_numeric(chains)
  assert_numeric(rungs)
  assert_numeric(GTI_pow)
  assert_logical(coupling_on)
  assert_logical(pb_markdown)
  assert_logical(silent)
  assert_non_null(IFRmodel$data)
  assert_non_null(IFRmodel$IFRparams)
  assert_non_null(IFRmodel$Infxnparams)
  assert_non_null(IFRmodel$Knotparams)
  assert_non_null(IFRmodel$paramdf)
  assert_non_null(IFRmodel$rho)
  assert_non_null(IFRmodel$Serotestparams)
  assert_non_null(IFRmodel$Noiseparams)
  assert_non_null(IFRmodel$modparam)
  assert_non_null(IFRmodel$sodparam)
  assert_non_null(IFRmodel$maxObsDay)
  assert_non_null(IFRmodel$demog)
  assert_eq(as.character(IFRmodel$demog$Strata),
            as.character(IFRmodel$data$prop_deaths$Strata[1:length(IFRmodel$IFRparams)]),
            message = "Strata within the demography data-frame must be in the same order as the strata in the observed deaths data frame")
  assert_eq(as.character(IFRmodel$data$obs_serology$Strata[1:length(IFRmodel$IFRparams)]),
            as.character(IFRmodel$data$prop_deaths$Strata[1:length(IFRmodel$IFRparams)]),
            message = "Strata within the observed serology data-frame must be in the same order as the strata in the observed deaths data frame")
  assert_eq(as.character(IFRmodel$data$obs_serology$Strata[1:length(IFRmodel$IFRparams)]),
            as.character(IFRmodel$demog$Strata),
            message = "Strata within the observed serology data-frame must be in the same order as the strata in thedemography data frame")

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
  # catches
  #..............................................................

  if (reparamIFR) {
    assert_non_null(IFRmodel$maxMa, message = "If performing reparameterization, must set a maximum Ma in the R6 class object")
  }

  if (reparamInfxn) {
    assert_non_null(IFRmodel$relInfxn, message = "If performing reparameterization, must set a relative infection point in the R6 class object")
  }

  if (reparamKnots) {
    assert_non_null(IFRmodel$relKnot, message = "If performing reparameterization, must set a relative knot point in the R6 class object")
  }

  # catch seroreversion
  check_serorev <- is.na(unlist(IFRmodel$paramdf[IFRmodel$paramdf$name %in% c("sero_rev_scale", "sero_rev_shape"),
                                                 c("min", "init", "max", "dsc1", "dsc2")]))
  if (all(check_serorev)) {
    account_serorev <- FALSE

  } else if (any(check_serorev)) {
    stop("You have inputted a mixture of NAs and values for Seroreversion. This is not allowed.")
  } else {
    account_serorev <- TRUE
  }

  #..................
  # use R to write Cpp prior and likelihood
  #..................
  logpriorfunc <- COVIDCurve:::make_user_Agg_logprior(IFRmodel, account_serorev = account_serorev,
                                                      reparamIFR = reparamIFR, reparamInfxn = reparamInfxn, reparamKnots = reparamKnots, reparamDelays = reparamDelays, reparamNe = reparamNe)
  loglikfunc <- COVIDCurve:::make_user_Agg_loglike(IFRmodel, account_serorev = account_serorev,
                                                   reparamIFR = reparamIFR, reparamInfxn = reparamInfxn, reparamKnots = reparamKnots, reparamDelays = reparamDelays, reparamNe = reparamNe)

  #..................
  # make misc
  #..................
  misc_list = list(rho = IFRmodel$rho,
                   rcensor_day = IFRmodel$rcensor_day,
                   days_obsd = IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel$Knotparams) + 1, # +1 because we set an internal knot for pos 1
                   n_sero_obs = length(unique(IFRmodel$data$obs_serology$SeroStartSurvey)),
                   sero_survey_start = unique(IFRmodel$data$obs_serology$SeroStartSurvey),
                   sero_survey_end = unique(IFRmodel$data$obs_serology$SeroEndSurvey),
                   max_seroday_obsd = max(IFRmodel$data$obs_serology$SeroEndSurvey),
                   demog = IFRmodel$demog$popN,
                   account_serorev = account_serorev)
  #..................
  # make data list
  #..................
  data_list <- list(obs_deaths = IFRmodel$data$obs_deaths$Deaths,
                    prop_strata_obs_deaths = IFRmodel$data$prop_deaths$PropDeaths,
                    obs_serologypos = IFRmodel$data$obs_serology$SeroPos,
                    obs_serologyn = IFRmodel$data$obs_serology$SeroN)

  #..................
  # make df param
  #..................
  df_params <-  IFRmodel$paramdf[, 1:4]
  # catch if no serorevesion
  if (!account_serorev) {
    df_params <- df_params %>%
      dplyr::filter(!name %in% c("sero_rev_scale", "sero_rev_shape"))
  }


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
                                silent = silent,
                                cluster = cluster
  )

  # apply thinning
  if (thinning > 0) {
    keepiters <- seq(from = thinning, to = (burnin+samples), by = thinning)
    mcmcout$output <- mcmcout$output %>%
      dplyr::group_by(chain, rung) %>%
      dplyr::filter(iteration %in% keepiters) %>%
      dplyr::ungroup(.)
  }

  #......................
  # account for reparameterizations
  #......................

  if (reparamIFR) {
    # account for ifr reparam
    IFRparams <- IFRmodel$IFRparams
    maxMa <- IFRmodel$maxMa
    scalars <- IFRparams[IFRparams != maxMa]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, maxMa]})
  }

  if (reparamKnots) {
    # account for knots (infxn X position) reparam
    Knotparams <- IFRmodel$Knotparams
    relKnot <- IFRmodel$relKnot
    scalars <- Knotparams[Knotparams != relKnot]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, relKnot]})
  }

  if (reparamInfxn) {
    # account for Infxn Y position reparam
    Infxnparams <- IFRmodel$Infxnparams
    relInfxn <- IFRmodel$relInfxn
    scalars <- Infxnparams[Infxnparams != relInfxn]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, relInfxn]})
  }

  if (reparamDelays) {
    # reparameterize serology rate relative to specificity
    mcmcout$output[, "sero_rate"] <- unname(unlist(1/mcmcout$output[, "spec"] * mcmcout$output[, "sero_rate"]))
  }

  if (reparamNe) {
    # reparameterize noise parameters
    Noiseparams <- IFRmodel$Noiseparams
    liftovercols <- Noiseparams[2:length(Noiseparams)]
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, Noiseparams[1]]})
  }

  # store input along with Dr.Jacoby output for later use
  inputs <- list(
    IFRmodel = IFRmodel,
    reparamIFR = reparamIFR,
    reparamInfxn = reparamInfxn,
    reparamKnots = reparamKnots,
    account_seroreversion = account_serorev,
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



