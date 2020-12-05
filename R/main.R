#' @title Run Aggregate Model
#' @details Wraps the Metropolic-Coupled MCMC Framework from Dr. Jacoby
#' @inheritParams drjacoby::run_mcmc
#' @param binomial_likelihood logical; Whether the binomial or the logit likelihood should be used
#' @param IFRmodel R6 class; Internal model object for COVIDCurve
#' @param reparamIFR logical; Whether IFRs should be reparameterized or inferred separately
#' @param reparamKnots logical; Whether infection knots (i.e. the x-coordinates of the infection spline) should be reparameterized or inferred separately
#' @param reparamInfxn logical; Whether infection curve (i.e. the  y-coordinates infection spline) should be reparameterized or inferred separately
#' @param thinning integer; The regular sequence to count by to thin MCMC posterior chain (iterations are kept as: \code{seq(from = thinning, to = (burnin+samples), by = thinning)}).
#' @export

run_IFRmodel_age <- function(IFRmodel,
                             binomial_likelihood = TRUE,
                             reparamIFR = TRUE, reparamInfxn = TRUE, reparamKnots = TRUE,
                             thinning = 0,
                             burnin = 1e3,
                             samples = 1e4,
                             rungs = 1,
                             chains = 5,
                             coupling_on = TRUE,
                             GTI_pow = 1.0,
                             beta_manual = NULL,
                             cluster = NULL,
                             pb_markdown = FALSE,
                             silent = FALSE) {
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
  assert_numeric_or_NULL(beta_manual)
  assert_logical(coupling_on)
  assert_logical(pb_markdown)
  assert_logical(silent)
  assert_non_null(IFRmodel$data)
  assert_non_null(IFRmodel$IFRparams)
  assert_non_null(IFRmodel$Infxnparams)
  assert_non_null(IFRmodel$Knotparams)
  assert_non_null(IFRmodel$paramdf)
  assert_non_null(IFRmodel$Serotestparams)
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
  # catch missing data
  if (binomial_likelihood) {
    if( any(c(is.na(IFRmodel$data$obs_serology$SeroPos), is.na(IFRmodel$data$obs_serology$SeroN))) ) {
      stop("Cannot have missing values of SeroPos or SeroN when considering the binomial likelihood")
    }
  } else {
    if( any(c(is.na(IFRmodel$data$obs_serology$SeroPrev), is.na(IFRmodel$data$obs_serology$SeroUCI), is.na(IFRmodel$data$obs_serology$SeroLCI))) ) {
      stop("Cannot have missing values of SeroPrev, SeroUCI, or SeroLCI when considering the logit likelihood")
    }
  }
  if (any(is.na(IFRmodel$data$prop_deaths$PropDeaths))) {
    stop("Cannot have missing proportions of cumulative deaths")
  }

  if (any(is.na(IFRmodel$data$obs_deaths$Deaths))) {
    warning("Missing daily deaths -- will skip over in likelihood")
  }
  # catch only one age group
  if (length(IFRmodel$inputs$IFRmodel$IFRparams) == 1) {
    stop("One strata currently not supported. Must specify more than strate to be estimated.")
  }

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
  if (length(IFRmodel$IFRparams) != 1) {
    if (length(IFRmodel$IFRparams) != length(IFRmodel$Noiseparams)) {
      stop("You must specificy the same number of IFR params and Noise params to estimate. Having
           a mix of no noise params and noise params among IFR params is not currently supported")
    }
  }

  # catch seroreversion
  if( all( c("sero_rev_shape", "sero_rev_scale") %in% IFRmodel$Serotestparams) )  {
    account_serorev <- TRUE
  } else if (sum( c("sero_rev_shape", "sero_rev_scale") %in% IFRmodel$Serotestparams)  == 1) {
    stop("Must specify both the Weibull Shape and Scale paramter for seroreversion to be considered. Specifying only one parameter is not an option")
  } else {
    account_serorev <- FALSE
  }

  #..................
  # use R to write Cpp prior and likelihood
  # and make data list
  # dependent on binomial or logit likelihood
  #..................
  logpriorfunc <- COVIDCurve:::make_user_Age_logprior(IFRmodel,
                                                      account_serorev = account_serorev,
                                                      reparamIFR = reparamIFR, reparamInfxn = reparamInfxn, reparamKnots = reparamKnots)
  loglikfunc <- COVIDCurve:::make_user_Age_loglike(IFRmodel,
                                                   binomial_likelihood = binomial_likelihood,
                                                   account_serorev = account_serorev,
                                                   reparamIFR = reparamIFR, reparamInfxn = reparamInfxn, reparamKnots = reparamKnots)



  #......................
  # catch user missing data -- liftover for Cpp
  #......................
  IFRmodel$data$obs_deaths$Deaths[is.na(IFRmodel$data$obs_deaths$Deaths)] <- -1
  #..................
  # make data list
  #..................
  if (binomial_likelihood) {
    data_list <- list(obs_deaths = IFRmodel$data$obs_deaths$Deaths,
                      prop_strata_obs_deaths = IFRmodel$data$prop_deaths$PropDeaths,
                      obs_serologypos = IFRmodel$data$obs_serology$SeroPos,
                      obs_serologyn = IFRmodel$data$obs_serology$SeroN)
  } else {
    # logit-normal transformation for likelihood
    # machine min tolerance is to protect against zeroes
    IFRmodel$data$obs_serology <- IFRmodel$data$obs_serology %>%
      dplyr::mutate(SeroSE = (COVIDCurve:::logit(SeroUCI + .Machine$double.xmin) - COVIDCurve:::logit(SeroLCI + .Machine$double.xmin)) / (2*1.96),
                    SeroMu = COVIDCurve:::logit(SeroPrev + .Machine$double.xmin))

    data_list <- list(obs_deaths = IFRmodel$data$obs_deaths$Deaths,
                      prop_strata_obs_deaths = IFRmodel$data$prop_deaths$PropDeaths,
                      obs_serologymu = IFRmodel$data$obs_serology$SeroMu,
                      obs_serologyse = IFRmodel$data$obs_serology$SeroSE)
  }


  #..................
  # make misc
  #..................
  misc_list = list(rcensor_day = IFRmodel$rcensor_day,
                   days_obsd = IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel$Knotparams) + 1, # +1 because we set an internal knot for pos 1
                   n_sero_obs = length(unique(IFRmodel$data$obs_serology$SeroStartSurvey)),
                   sero_survey_start = unique(IFRmodel$data$obs_serology$SeroStartSurvey),
                   sero_survey_end = unique(IFRmodel$data$obs_serology$SeroEndSurvey),
                   max_seroday_obsd = max(IFRmodel$data$obs_serology$SeroEndSurvey),
                   demog = IFRmodel$demog$popN,
                   account_serorev = account_serorev)

  #..................
  # make df param
  #..................
  # columns 5 & 6 used for prior distributions
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
                                beta_manual = beta_manual,
                                pb_markdown = pb_markdown,
                                silent = silent,
                                cluster = cluster)

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


  # store MCMC input along with Dr.Jacoby output for later use
  inputs <- list(
    IFRmodel = IFRmodel,
    reparamIFR = reparamIFR,
    reparamInfxn = reparamInfxn,
    reparamKnots = reparamKnots,
    account_seroreversion = account_serorev,
    binomial_likelihood = binomial_likelihood,
    burnin = burnin,
    samples = samples,
    chains = chains)

  if (rungs > 1) {
    inputs <- append(inputs, list(rungs = rungs,
                                  GTI_pow = GTI_pow,
                                  beta_manual = beta_manual,
                                  coupling_on = coupling_on))
  }

  # out
  ret <- list(
    inputs = inputs,
    mcmcout = mcmcout
  )
  class(ret) <- c("IFRmodel_inf")
  return(ret)
}



