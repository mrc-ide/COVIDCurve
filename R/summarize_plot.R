#' IFR-Inference-Aggregate-Model, overload summary
#' @noRd
#' @export

summary.IFRmodel <- function(object, ...){
  cat(crayon::cyan("*~*~*~*~ IFR Inference Model ~*~*~*~*\n"))
  cat(crayon::green("IFR strata params: "), paste(object$IFRparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Point Params: "), paste(object$Infxnparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Knot Params: "), paste(object$Knotparams, collapse = ", "), "\n")
  cat(crayon::magenta("Serology Test Parameters: "), paste(object$Serotestparams, collapse = ", "), "\n")
  cat(crayon::magenta("Serology Day Parameters: "), paste(object$Serodayparams, collapse = ", "), "\n")
  cat(crayon::yellow("Model Type: "), object$level, "\n")
  cat(crayon::yellow("Total Population Size: "), object$popN, "\n")
  cat(crayon::yellow("Prob. of Infection Given Strata: "),  paste(round(object$rho, 2), collapse = ", "), "\n")
  cat(crayon::yellow("Mean Delay of Onset-to-Death: "), object$mod, "\n")
  cat(crayon::yellow("Coef. Var. Delay of Onset-to-Death: "), object$sod, "\n")
}


#' IFR-Inference-Aggregate-Model, overload print
#' @noRd
#' @export

print.IFRmodel <- function(x, ...){
  cat(crayon::cyan("*~*~*~*~ IFR Inference Model ~*~*~*~*\n"))
  cat(crayon::green("IFR strata params: "), paste(x$IFRparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Point Params: "), paste(x$Infxnparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Knot Params: "), paste(x$Knotparams, collapse = ", "), "\n")
  cat(crayon::magenta("Serology Parameters: "), paste(x$Serotestparams, collapse = ", "), "\n")
  cat(crayon::magenta("Serology Parameters: "), paste(x$Serodayparams, collapse = ", "), "\n")
  cat(crayon::yellow("Model Type: "), x$level, "\n")
  cat(crayon::yellow("Total Population Size: "), x$popN, "\n")
  cat(crayon::yellow("Prob. of Infection Given Strata: "),  paste(round(x$rho, 2), collapse = ", "), "\n")
  cat(crayon::yellow("Mean Delay of Onset-to-Death: "), x$mod, "\n")
  cat(crayon::yellow("Coef. Var. Delay of Onset-to-Death: "), x$sod, "\n")
}

#' @title Get Credible Intervals for Parameters from Sampling Iterations
#' @param IFRmodel_inf R6 class; The result of the IFR Model MCMC run along with the model input.
#' @param what character; Specify which parameters you want: "IFRparams", "Infxnparams", "Serotestparams", or "Serodayparams"
#' @param whichrung character; Specify which rung to sample from (default is rung1)
#' @param by_chain logical; Whether or not credible intervals should be reported with respect to individual chains (TRUE) or not.
#' @importFrom magrittr %>%
#' @export

get_cred_intervals <- function(IFRmodel_inf, what, whichrung = "rung1", by_chain = TRUE) {
  assert_custom_class(IFRmodel_inf$inputs$IFRmodel, "IFRmodel")
  assert_custom_class(IFRmodel_inf$mcmcout, "drjacoby_output")
  assert_custom_class(IFRmodel_inf, "IFRmodel_inf")
  assert_in(what, c("IFRparams", "Knotparams", "Infxnparams", "Serotestparams", "Serodayparams"))
  assert_string(whichrung)
  assert_logical(by_chain)

  # grouping vars
  switch(paste0(what, "-", by_chain),

         "IFRparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel_inf$inputs$IFRmodel$IFRparams)
         },
         "IFRparams-FALSE"={
           groupingvar <- "param"
           params <-  c("iteration", IFRmodel_inf$inputs$IFRmodel$IFRparams)
         },

         "Knotparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel_inf$inputs$IFRmodel$Knotparams)
         },
         "Knotparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel_inf$inputs$IFRmodel$Knotparams)
         },

         "Infxnparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel_inf$inputs$IFRmodel$Infxnparams)
         },
         "Infxnparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel_inf$inputs$IFRmodel$Infxnparams)
         },

         "Serotestparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel_inf$inputs$IFRmodel$Serotestparams)
         },
         "Serotestparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel_inf$inputs$IFRmodel$Serotestparams)
         }

         "Serodayparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel_inf$inputs$IFRmodel$Serodayparams)
         },
         "Serodayparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel_inf$inputs$IFRmodel$Serodayparams)
         }
  )

  IFRmodel_inf$mcmcout$output %>%
    dplyr::filter(stage == "sampling" & rung == whichrung) %>%
    dplyr::select_at(params) %>%
    tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
    dplyr::group_by_at(groupingvar) %>%
    dplyr::summarise(
      min = min(est),
      LCI = quantile(est, 0.025),
      median = median(est),
      mean = mean(est),
      UCI = quantile(est, 0.975),
      max = max(est),
      ESS = coda::effectiveSize(coda::as.mcmc(est)),
      GewekeZ = coda::geweke.diag(coda::as.mcmc(est))[[1]],
      GewekeP = dnorm(GewekeZ)
    )
}




#' @title Draw posterior results from the Cubic Spline
#' @details Given sampling iterations with posterior-log-likes greater than or equal to a specific threshold, posterior results for the linear spline are generated. Assumed that the spline was fit in "un-transformed" space
#' @inheritParams get_cred_intervals
#' @param dwnsmpl integer; Number of posterior results to draw -- weighted by posterior prob
#' @importFrom magrittr %>%
#' @export

draw_posterior_infxn_points_cubic_splines <- function(IFRmodel_inf, whichrung = "rung1", dwnsmpl, by_chain = TRUE) {
  assert_custom_class(IFRmodel_inf$inputs$IFRmodel, "IFRmodel")
  assert_custom_class(IFRmodel_inf, "IFRmodel_inf")
  assert_custom_class(IFRmodel_inf$mcmcout, "drjacoby_output")
  assert_pos_int(dwnsmpl)
  assert_string(whichrung)
  assert_logical(by_chain)
  #......................
  # fitler to sampling and by rung
  #......................
  mcmcout.nodes <-  IFRmodel_inf$mcmcout$output %>%
    dplyr::filter(stage == "sampling" & rung == whichrung)

  #......................
  # sample by CI limit and make infxn curves
  #......................
  mcmcout.nodes <- mcmcout.nodes %>%
    dplyr::mutate(logposterior = loglikelihood + logprior)
  # Log-Sum-Exp trick
  convert_post_probs <- function(logpost) {
    exp(logpost - (log(sum(exp(logpost - max(logpost)))) + max(logpost)))
  }
  probs <- convert_post_probs(mcmcout.nodes$logposterior)
  # downsample
  dwnsmpl_rows <- sample(1:nrow(mcmcout.nodes), size = dwnsmpl,
                         prob = probs)
  dwnsmpl_rows <- sort(dwnsmpl_rows)
  mcmcout.nodes <- mcmcout.nodes[dwnsmpl_rows, ]

  #......................
  # get natural cubic gradients
  #......................
  # internal function, liftover cpp likelihood to get infxn curve
  # NOTE, this is extremely sensitive to the placements of the Cpp source file and therefore, is not generalizable
  fitcurve_string <- COVIDCurve:::make_user_Agg_loglike(IFRmodel = IFRmodel_inf$inputs$IFRmodel,
                                                       reparamIFR = FALSE,
                                                       reparamKnots = FALSE,
                                                       reparamInfxn = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
  # pull out pieces I need
  fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
  fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "std::vector\\<double\\> cumm_infxn_spline\\(infxn_spline.size\\(\\)\\);", n = 2)[,1]
  fitcurve_string <- paste(fitcurve_start, fitcurve_curve, "Rcpp::List ret = Rcpp::List::create(infxn_spline); return ret;}", collapse = "")
  Rcpp::cppFunction(fitcurve_string)

  #......................
  # inputs needed for cpp function
  #......................
  misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
                   pgmms = IFRmodel_inf$inputs$IFRmodel$gamma_lookup,
                   level = ifelse(IFRmodel_inf$inputs$IFRmodel$level == "Cumulative", TRUE, FALSE),
                   popN = IFRmodel_inf$inputs$IFRmodel$popN,
                   rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
                   days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams)+1,
                   n_sero_obs = length(IFRmodel_inf$inputs$IFRmodel$Serodayparams))

  datin <- list("obs_deaths" = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
                "obs_serologyrate" = IFRmodel_inf$inputs$IFRmodel$data$obs_serologyrate)


  #......................
  # split, run, recombine
  #......................
  cpp_function_wrapper <- function(params, data, misc) {
    paramsin <- unlist(params[c(IFRmodel_inf$inputs$IFRmodel$IFRparams,
                                IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                                IFRmodel_inf$inputs$IFRmodel$Knotparams,
                                IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                                IFRmodel_inf$inputs$IFRmodel$Serodayparams)])
    infxns <- unlist(loglike(params = paramsin,
                             param_i = 1,
                             data = datin,
                             misc = misc_list))
    ret <- data.frame(time = 1:length(infxns),
                      infxns = infxns)
    return(ret)

  }

  mcmcout.node.rows <- split(mcmcout.nodes, 1:nrow(mcmcout.nodes))
  mcmcout.nodes$infxncurves <- purrr::map(mcmcout.node.rows, cpp_function_wrapper,
                                          data = datin, misc = misc_list)

  #......................
  # tidy
  #......................
  # keep params around for convenience
  if (by_chain) {
    plotdat <- mcmcout.nodes %>%
      dplyr::select(c("chain",
                      IFRmodel_inf$inputs$IFRmodel$IFRparams,
                      IFRmodel_inf$inputs$IFRmodel$Knotparams,
                      IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                      IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                      IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                      "infxncurves")) %>%
      dplyr::group_by(chain) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      dplyr::ungroup(chain) %>%
      tidyr::unnest(cols = "infxncurves")

    plotObj <- ggplot2::ggplot() +
      ggplot2::geom_line(data = plotdat, mapping =  ggplot2::aes(time, infxns, group = sim), alpha = 0.25,
                         lwd = 0.5, color = "#d9d9d9") +
      ggplot2::xlab("Time") +  ggplot2::ylab("Num. Infxns")  +
      ggplot2::labs(title = "Posterior Draws of the Infection Curve") +
      ggplot2::facet_wrap(. ~ chain) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title =  ggplot2::element_text(family = "Helvetica", face = "bold", vjust = 0.5,  hjust = 0.5, size = 18),
        plot.subtitle =  ggplot2::element_text(family = "Helvetica", face = "bold", vjust = 0.5,  hjust = 0.5, size = 18),
        axis.title =  ggplot2::element_text(family = "Helvetica", face = "bold", hjust = 0.5, vjust = 0.5, size = 16),
        axis.text.x =  ggplot2::element_text(family = "Helvetica", angle = 45, hjust = 0.5, vjust = 0.5, size = 15),
        axis.text.y =  ggplot2::element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 15),
        panel.background =  ggplot2::element_blank(),
        plot.background =  ggplot2::element_blank(),
        axis.line =  ggplot2::element_line(color = "#000000", size = 1.2),
        legend.position = "none")

    if (IFRmodel_inf$inputs$IFRmodel$rcensor_day < IFRmodel_inf$inputs$IFRmodel$maxObsDay) {
      plotObj <- plotObj +
        ggplot2::geom_vline(xintercept = IFRmodel_inf$inputs$IFRmodel$rcensor_day, linetype = "dashed", size = 1.2, color = "#de2d26")
    }


  } else {
    # keep params around for convenience
    plotdat <- mcmcout.nodes %>%
      dplyr::select(c(IFRmodel_inf$inputs$IFRmodel$IFRparams,
                      IFRmodel_inf$inputs$IFRmodel$Knotparams,
                      IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                      IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                      IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                      "infxncurves")) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      tidyr::unnest(cols = "infxncurves")

    plotObj <-  ggplot2::ggplot() +
      ggplot2::geom_line(data = plotdat, mapping =  ggplot2::aes(time, infxns, group = sim), alpha = 0.25,
                         lwd = 0.5, color = "#d9d9d9") +
      ggplot2::xlab("Time") +  ggplot2::ylab("Num. Infxns")  +
      ggplot2::labs(title = "Posterior Draws of the Infection Curve") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(family = "Helvetica", face = "bold", vjust = 0.5,  hjust = 0.5, size = 18),
        plot.subtitle = ggplot2::element_text(family = "Helvetica", face = "bold", vjust = 0.5,  hjust = 0.5, size = 18),
        axis.title = ggplot2::element_text(family = "Helvetica", face = "bold", hjust = 0.5, vjust = 0.5, size = 16),
        axis.text.x = ggplot2::element_text(family = "Helvetica", angle = 45, hjust = 0.5, vjust = 0.5, size = 15),
        axis.text.y = ggplot2::element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 15),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "#000000", size = 1.2),
        legend.position = "none")

    if (IFRmodel_inf$inputs$IFRmodel$rcensor_day < IFRmodel_inf$inputs$IFRmodel$maxObsDay) {
      plotObj <- plotObj +
        ggplot2::geom_vline(xintercept = IFRmodel_inf$inputs$IFRmodel$rcensor_day, linetype = "dashed", size = 1.2, color = "#de2d26")
    }
  }

  #......................
  # out
  #......................
  ret <- list(
    plotdata = plotdat,
    plotObj = plotObj
  )
  return(ret)
}

#' @title Posterior Check for Infections to Death
#' @details Given sampling iterations with posterior-log-likes greater than or equal to a specific threshold, posterior results for the linear spline are generated. Assumed that the spline was fit in "un-transformed" space
#' @inheritParams draw_posterior_infxn_points_cubic_splines
#' @importFrom magrittr %>%
#' @export

posterior_check_infxns_to_death <- function(IFRmodel_inf, whichrung = "rung1",
                                            dwnsmpl, by_chain = FALSE) {

  postdat <- COVIDCurve::draw_posterior_infxn_points_cubic_splines(IFRmodel_inf, whichrung = whichrung,
                                                                   dwnsmpl = dwnsmpl, by_chain = by_chain)$plotdat
  # set up function to draw posterior deaths
  postdat.sims <- split(postdat, factor(postdat$sim))
  draw_post_deaths <- function(postdatsim){
    gmmlkup <- stats::dgamma(postdatsim$time, shape = 1/(IFRmodel_inf$inputs$IFRmodel$sod^2),
                             scale = IFRmodel_inf$inputs$IFRmodel$mod*IFRmodel_inf$inputs$IFRmodel$sod^2, log = FALSE)

    # exp deaths day and strata
    exp_death.day <- rep(0, length = IFRmodel_inf$inputs$IFRmodel$maxObsDay)
    # spread infxns out to day when they may or may not die
    for (i in 1:nrow(postdatsim)) {
      for (j in (i+1):(nrow(postdatsim) + 1)) {
        delta <- j - i
        exp_death.day[j-1] <- exp_death.day[j-1] + postdatsim$infxns[i] * gmmlkup[delta]
      }
    }

    exp_death.day.strata <- exp_death.day %*% t(IFRmodel_inf$inputs$IFRmodel$pa) # account for pa
    exp_death.day.strata <- exp_death.day.strata * postdatsim[, IFRmodel_inf$inputs$IFRmodel$IFRparams]

    #......................
    # tidy up and out
    #......................
    out <- cbind.data.frame(time = 1:nrow(exp_death.day.strata), exp_death.day.strata)
    colnames(out)[2:ncol(out)] <- paste0("deaths_", IFRmodel_inf$inputs$IFRmodel$IFRparams)
    return(out)
  }
  # get post deaths
  postdat.curves <- postdat %>%
    dplyr::select("sim", IFRmodel_inf$inputs$IFRmodel$IFRparams) %>%
    dplyr::filter(!duplicated(.)) %>%
    dplyr::mutate(
      post_deaths = purrr::map(postdat.sims, draw_post_deaths)
    ) %>%
    tidyr::unnest(cols = post_deaths)

  return(postdat.curves)
}





