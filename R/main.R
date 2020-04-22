#' @title Wrap the Metropolic-Coupled MCMC Framework from Dr. Jacoby
#' @param level character; Must be either Time-Series or Cumulative
#' @param LogLike global function;
#' @param LogPrior global function;
#' @param misclist list; Required list object that must contain the following items:
#'   elements:
#'   \itemize{
#'     \item \code{min_day} - the presumed minimum day of the epidemic (required only for time-series likelihoods)
#'     \item \code{curr_day} - the presumed minimum day of the epidemic.
#'     \item \code{pa} - A vector of the attack rates for the given age bands
#'   }
#' @inheritParams drjacoby::run_mcmc
#'
#' @export


wrap_drjacoby_mcmc <- function(data, df_params, misclist, LogLike, LogPrior, level,
                               burnin = 1e3, samples = 1e3, chains = 3,
                               rungs = 1, GTI_pow = 3, coupling_on = T,
                               pb_markdown = F, silent = T) {
  #..............................................................
  # Assertions that are specific to this project
  #..............................................................
  assert_in(x = level, y = c("Time-Series", "Cumulative"))

  if (level == "Time-Series") {
    #TODO specific assertions
    # split in to days
    dat <- split(data, factor(data$obs_day))
    dat <- purrr::map(dat, "day_deaths")
    # list for Dr. Jacoby
    data_list <- list(obs_deaths = dat)

  } else if (level == "Cumulative") {
    #TODO specific assertions
    assert_vector(data)
    sapply(data, function(x){assert_numeric(x)})
    # put in list format for Dr. Jacoby
    data_list <- list(obs_deaths = unname(data))
  }

  #..............................................................
  # Dr Jacoby
  #..............................................................

  mcmcout <- drjacoby::run_mcmc(data = data_list,
                                df_params = df_params,
                                misc = misclist,
                                loglike = LogLike,
                                logprior = LogPrior,
                                burnin = burnin,
                                samples = samples,
                                chains = chains,
                                rungs = rungs,
                                coupling_on = coupling_on,
                                GTI_pow = GTI_pow,
                                pb_markdown = pb_markdown,
                                silent = silent
                                )

  return(mcmcout)
}
