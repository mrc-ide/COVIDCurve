#' @title Run Model
#' @param level character; Must be either Time-Series or Cumulative
#' @details Wraps the Metropolic-Coupled MCMC Framework from Dr. Jacoby
#' @inheritParams drjacoby::run_mcmc
#'
#' @export


run_modinf <- function(modinf, reparamIFR = T,
                       burnin = 1e3, samples = 1e3, chains = 3,
                       rungs = 1, GTI_pow = 3, coupling_on = T,
                       pb_markdown = F, silent = T) {
  #..................
  # assertions
  #..................
  assert_custom_class(modinf, "Inference-Model")
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
  assert_non_null(modinf$level)
  assert_non_null(modinf$IFRparams)
  assert_non_null(modinf$Infxnparams)
  assert_non_null(modinf$paramdf)
  assert_non_null(modinf$knots)
  assert_non_null(modinf$pa)


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

    logpriorfunc <- COVIDCurve:::make_user_logprior_reparam(modinf)
    loglikfunc <- COVIDCurve:::make_user_loglike_reparam(modinf)

  } else {
    logpriorfunc <- COVIDCurve:::make_user_logprior_noreparam(modinf)
    loglikfunc <- COVIDCurve:::make_user_loglike_noreparam(modinf)

  }


  #..................
  # make misc
  #..................
  misc_list = list(pa = modinf$pa,
                   pgmms = modinf$gamma_lookup,
                   knots = modinf$knots,
                   level = ifelse(modinf$level == "Cumulative", TRUE, FALSE))

  #..................
  # make data list
  #..................
  if (modinf$level == "Time-Series"){
    data_list <- split(modinf$data$Deaths, factor(modinf$data$ObsDay))
    data_list <- unname(unlist(data_list))
    data_list <- list(obs_deaths = data_list)

  } else if (modinf$level == "Cumulative") {
    data_list <- list(obs_deaths = unname(modinf$data$Deaths))
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
    maxMa <- IFRparams$name[which(IFRparams$init == max(IFRparams$init))]
    scalars <- IFRparams$name[IFRparams$name != maxMa]

    liftovercols <- colnames(mcmcout$output) %in% scalars
    liftovercols.list <- mcmcout$output[, liftovercols]
    liftovercols.list <- lapply(colnames(liftovercols.list), function(x){liftovercols.list[,x]})
    mcmcout$output[, liftovercols] <- sapply(liftovercols.list, function(x) {x * mcmcout$output[, maxMa]})
  }

  return(mcmcout)
}
