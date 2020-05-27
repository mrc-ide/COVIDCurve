#' @title Run Aggregate Model
#' @details Wraps the Metropolic-Coupled MCMC Framework from Dr. Jacoby
#' @inheritParams drjacoby::run_mcmc
#'
#' @export

run_modinf_agg <- function(modinf, reparamIFR = T,
                           burnin = 1e3, samples = 1e3, chains = 3,
                           rungs = 1, GTI_pow = 3, coupling_on = T,
                           pb_markdown = F, silent = T) {
  #..................
  # assertions
  #..................
  assert_custom_class(modinf, "Inference-Aggregate-Model")
  assert_logical(reparamIFR)
  assert_numeric(burnin)
  assert_numeric(samples)
  assert_numeric(chains)
  assert_numeric(rungs)
  assert_numeric(GTI_pow)
  assert_logical(coupling_on)
  assert_logical(pb_markdown)
  assert_logical(silent)
  assert_non_null(modinf$level)
  assert_non_null(modinf$data)
  assert_non_null(modinf$IFRparams)
  assert_non_null(modinf$Infxnparams)
  assert_non_null(modinf$paramdf)
  assert_non_null(modinf$knots)
  assert_non_null(modinf$pa)
  assert_non_null(modinf$Seroparams)
  assert_non_null(modinf$popN)

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
    assert_non_null(modinf$maxMa, message = "If performing reparameterization, must set a maximum Ma in the R6 class object")
  }


  logpriorfunc <- COVIDCurve:::make_user_Agg_logprior(modinf, reparamIFR = reparamIFR)
  loglikfunc <- COVIDCurve:::make_user_Agg_loglike(modinf, reparamIFR = reparamIFR)



  #..................
  # make misc
  #..................
  misc_list = list(pa = modinf$pa,
                   pgmms = modinf$gamma_lookup,
                   knots = modinf$knots,
                   level = ifelse(modinf$level == "Cumulative", TRUE, FALSE),
                   popN = modinf$popN)

  #..................
  # make data list
  #..................
  if (modinf$level == "Time-Series"){
    data_list <- split(modinf$data$obs_deaths$Deaths, factor(modinf$data$obs_deaths$ObsDay))
    data_list <- unname(unlist(data_list))
    data_list <- list(obs_deaths = data_list,
                      obs_serologyrate = modinf$data$obs_serologyrate)

  } else if (modinf$level == "Cumulative") {
    data_list <- list(obs_deaths = unname(modinf$data$obs_deaths$Deaths),
                      obs_serologyrate = modinf$data$obs_serologyrate)
  }

  #..................
  # make df param
  #..................
  df_params <-  modinf$paramdf[, 1:4]

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
    # account for reparam
    #..................
    IFRparams <- modinf$paramdf[modinf$paramdf$name %in% modinf$IFRparams, ]
    maxMa <- modinf$maxMa
    scalars <- IFRparams$name[IFRparams$name != maxMa]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, maxMa]})
  }

  return(mcmcout)
}

















#............................................................
# line list
#...........................................................
#' @title Run Line-List Model
#' @details Wraps the Metropolic-Coupled MCMC Framework from Dr. Jacoby
#' @inheritParams drjacoby::run_mcmc
#'
#' @export


run_modinf_linelist <- function(modinf, reparamIFR = T,
                                burnin = 1e3, samples = 1e3, chains = 3,
                                rungs = 1, GTI_pow = 3, coupling_on = T,
                                pb_markdown = F, silent = T) {
  #..................
  # assertions
  #..................
  assert_custom_class(modinf, "Inference-LineList-Model")
  assert_logical(reparamIFR)
  assert_numeric(burnin)
  assert_numeric(samples)
  assert_numeric(chains)
  assert_numeric(rungs)
  assert_numeric(GTI_pow)
  assert_logical(coupling_on)
  assert_logical(pb_markdown)
  assert_logical(silent)
  assert_non_null(modinf$data)
  assert_non_null(modinf$IFRparams)
  assert_non_null(modinf$DeathModparam)
  assert_non_null(modinf$DeathSodparam)
  assert_non_null(modinf$RecovMorparam)
  assert_non_null(modinf$RecovSorparam)
  assert_non_null(modinf$paramdf)


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
    assert_same_length(max(modinf$paramdf$init[modinf$paramdf$name %in% modinf$IFRparams]), 1,
                       message = "One IFR-Param must be considered the max
                                  (i.e. in your paramdf, one IFR must have highest init value
                                  for other IFRs to be scaled towards)")

    logpriorfunc <- COVIDCurve:::make_user_LineList_logprior_reparam(modinf)
    loglikfunc <- COVIDCurve:::make_user_LineList_loglike_reparam(modinf)

  } else {
    logpriorfunc <- COVIDCurve:::make_user_LineList_logprior_noreparam(modinf)
    loglikfunc <- COVIDCurve:::make_user_LineList_loglike_noreparam(modinf)

  }


  #..................
  # make misc
  #..................
  misc_list = list(IFRparams = as.numeric(factor(modinf$IFRparams)))

  #..................
  # make data list
  #..................
  death_onset_day <- modinf$data %>%
    dplyr::filter(Outcome == "Death") %>%
    dplyr::select(c("OnsetDay")) %>%
    unlist(.) %>%
    unname(.)

  death_event_day <- modinf$data %>%
    dplyr::filter(Outcome == "Death") %>%
    dplyr::select(c("EventDay")) %>%
    unlist(.) %>%
    unname(.)

  death_group <- modinf$data %>%
    dplyr::filter(Outcome == "Death") %>%
    dplyr::select(c("AgeGroup")) %>%
    unlist(.) %>%
    unname(.)
  death_group <- factor(death_group)

  recovery_onset_day <- modinf$data %>%
    dplyr::filter(Outcome == "Recovery") %>%
    dplyr::select(c("OnsetDay")) %>%
    unlist(.) %>%
    unname(.)

  recovery_event_day <- modinf$data %>%
    dplyr::filter(Outcome == "Recovery") %>%
    dplyr::select(c("EventDay")) %>%
    unlist(.) %>%
    unname(.)

  recovery_group <- modinf$data %>%
    dplyr::filter(Outcome == "Recovery") %>%
    dplyr::select(c("AgeGroup")) %>%
    unlist(.) %>%
    unname(.)
  recovery_group <- factor(recovery_group)

  data_list <- list(death_interval = death_event_day - death_onset_day,
                    death_group = as.numeric(death_group),
                    recovery_interval = recovery_event_day -recovery_onset_day,
                    recovery_group = as.numeric(recovery_group)
  )


  #..................
  # make df param
  #..................
  df_params <-  modinf$paramdf[, 1:4]

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
    # account for reparam
    #..................
    IFRparams <- modinf$paramdf[modinf$paramdf$name %in% modinf$IFRparams, ]
    maxMa <- IFRparams$name[which(IFRparams$init == max(IFRparams$init))]
    scalars <- IFRparams$name[IFRparams$name != maxMa]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, maxMa]})
  }

  return(mcmcout)
}


