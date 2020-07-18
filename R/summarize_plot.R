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
  cat(crayon::yellow("Total Population Size: "), sum(object$demog$popN), "\n")
}


#' IFR-Inference-Aggregate-Model, overload print
#' @noRd
#' @export

print.IFRmodel <- function(x, ...){
  cat(crayon::cyan("*~*~*~*~ IFR Inference Model ~*~*~*~*\n"))
  cat(crayon::green("IFR strata params: "), paste(x$IFRparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Point Params: "), paste(x$Infxnparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Knot Params: "), paste(x$Knotparams, collapse = ", "), "\n")
  cat(crayon::magenta("Serology Test Parameters: "), paste(x$Serotestparams, collapse = ", "), "\n")
  cat(crayon::magenta("Serology Day Parameters: "), paste(x$Serodayparams, collapse = ", "), "\n")
  cat(crayon::yellow("Total Population Size: "), sum(x$demog$popN), "\n")
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
  assert_in(what, c("IFRparams", "Knotparams", "Infxnparams", "Serotestparams", "Serodayparams", "Noiseparams", "TODparams"))
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
         },

         "Serodayparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel_inf$inputs$IFRmodel$Serodayparams)
         },
         "Serodayparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel_inf$inputs$IFRmodel$Serodayparams)
         },

         "Noiseparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel_inf$inputs$IFRmodel$Noiseparams)
         },
         "Noiseparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel_inf$inputs$IFRmodel$Noiseparams)
         },

         "TODparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", c(IFRmodel_inf$inputs$IFRmodel$modparam, IFRmodel_inf$inputs$IFRmodel$sodparam))
         },
         "TODparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", c(IFRmodel_inf$inputs$IFRmodel$modparam, IFRmodel_inf$inputs$IFRmodel$sodparam))
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

#' @title Log-Sum-Exp trick
#' @param logpost numeric; log-posterior
#' @noRD
#' @export

convert_post_probs <- function(logpost) {
  exp(logpost - (log(sum(exp(logpost - max(logpost)))) + max(logpost)))
}

#' @title Draw posterior results from the Infection Curve (Cubic Spline)
#' @details Given sampling iterations with posterior-log-likes greater than or equal to a specific threshold, posterior results for the linear spline are generated. Assumed that the spline was fit in "un-transformed" space
#' @inheritParams get_cred_intervals
#' @param dwnsmpl integer; Number of posterior results to draw -- weighted by posterior prob
#' @importFrom magrittr %>%
#' @export

draw_posterior_infxn_cubic_splines <- function(IFRmodel_inf, whichrung = "rung1", dwnsmpl, by_chain = TRUE) {
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
                                                        reparamInfxn = FALSE,
                                                        reparamSpec = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
  # pull out pieces I need
  fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
  fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "double cum_infxn_check = 0.0;", n = 2)[,1]
  fitcurve_string <- paste(fitcurve_start, fitcurve_curve,
                           "std::vector<std::vector<double>> infxn_spline_strata(days_obsd, std::vector<double>(stratlen));
                            for (int i = 0; i < days_obsd; i++) {
                              for (int a = 0; a < stratlen; a++) {
                                infxn_spline_strata[i][a] =  ne[a] * infxn_spline[i];
                              }
                           }",
                           "Rcpp::List ret = Rcpp::List::create(infxn_spline_strata); return ret;}",
                           collapse = "")
  Rcpp::cppFunction(fitcurve_string)

  #......................
  # inputs needed for cpp function
  #......................
  misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
                   demog = IFRmodel_inf$inputs$IFRmodel$demog$popN,
                   rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
                   days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams)+1,
                   n_sero_obs = length(IFRmodel_inf$inputs$IFRmodel$Serodayparams))

  datin <- list("obs_deaths" = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
                "obs_serology" = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroPrev)


  #......................
  # split, run, recombine
  #......................
  cpp_function_wrapper <- function(params, data, misc) {
    paramsin <- unlist(params[c(IFRmodel_inf$inputs$IFRmodel$modparam,
                                IFRmodel_inf$inputs$IFRmodel$sodparam,
                                IFRmodel_inf$inputs$IFRmodel$IFRparams,
                                IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                                IFRmodel_inf$inputs$IFRmodel$Knotparams,
                                IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                                IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                                IFRmodel_inf$inputs$IFRmodel$Noiseparams)])

    # run efficient cpp code
    infxns <- loglike(params = paramsin,
                      param_i = 1,
                      data = datin,
                      misc = misc_list)[[1]]
    infxns <- infxns %>%
      do.call("rbind.data.frame", .)

    # out
    colnames(infxns) <- paste0("infxns_", IFRmodel_inf$inputs$IFRmodel$IFRparams)
    ret <- cbind.data.frame(time = 1:nrow(infxns),
                            infxns)
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
                      IFRmodel_inf$inputs$IFRmodel$modparam,
                      IFRmodel_inf$inputs$IFRmodel$sodparam,
                      IFRmodel_inf$inputs$IFRmodel$IFRparams,
                      IFRmodel_inf$inputs$IFRmodel$Knotparams,
                      IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                      IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                      IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                      IFRmodel_inf$inputs$IFRmodel$Noiseparams,
                      "infxncurves")) %>%
      dplyr::group_by(chain) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      dplyr::ungroup(chain) %>%
      tidyr::unnest(cols = "infxncurves") %>%
      dplyr::mutate(totinfxns = rowSums(dplyr::select(., dplyr::starts_with("infxns_"))))

    plotObj <- ggplot2::ggplot() +
      ggplot2::geom_line(data = plotdat, mapping =  ggplot2::aes(x = time, y = totinfxns, group = sim), alpha = 0.25,
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
      dplyr::select(c(IFRmodel_inf$inputs$IFRmodel$modparam,
                      IFRmodel_inf$inputs$IFRmodel$sodparam,
                      IFRmodel_inf$inputs$IFRmodel$IFRparams,
                      IFRmodel_inf$inputs$IFRmodel$Knotparams,
                      IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                      IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                      IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                      IFRmodel_inf$inputs$IFRmodel$Noiseparams,
                      "infxncurves")) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      tidyr::unnest(cols = "infxncurves") %>%
      dplyr::mutate(totinfxns = rowSums(dplyr::select(., dplyr::starts_with("infxns_"))))

    plotObj <-  ggplot2::ggplot() +
      ggplot2::geom_line(data = plotdat, mapping =  ggplot2::aes(x = time, y = totinfxns, group = sim), alpha = 0.25,
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
#' @inheritParams draw_posterior_infxn_cubic_splines
#' @importFrom magrittr %>%
#' @export

posterior_check_infxns_to_death <- function(IFRmodel_inf, whichrung = "rung1",
                                            dwnsmpl, by_chain = FALSE) {
  #......................
  # draw infection curves
  #......................
  postdat <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf, whichrung = whichrung,
                                                            dwnsmpl = dwnsmpl, by_chain = by_chain)$plotdat
  # split infection data for looping through later
  postdat.sims <- split(postdat, factor(postdat$sim))

  #......................
  # draw deaths from infections
  #......................
  # make eff cpp function for drawing deaths
  src <- "Rcpp::List get_post_deaths(Rcpp::NumericVector infxns, Rcpp::NumericVector ifr, int daylen, int stratlen, double mod, double sod) {

      // get gamma look up table
      std::vector<double> gmmlkup(daylen + 1);
      for (int i = 0; i < (daylen+1); i++) {
        gmmlkup[i] = R::dgamma(i, 1/pow(sod,2), mod*pow(sod,2), false);
      }
      std::vector<double> infxns_raw = Rcpp::as< std::vector<double> >(infxns);
      std::vector<std::vector<double>> infxns_strata(daylen, std::vector<double>(stratlen));
      int iter = 0;
      // recast infxns to matrix -- r unlists by row (so column1, column2, ...)
      for (int j = 0; j < stratlen; j++) {
        for (int i = 0; i < daylen; i++) {
          infxns_strata[i][j] = infxns_raw[iter];
          iter++;
        }
      }
      std::vector<std::vector<double>> exp_death_day_strata(daylen, std::vector<double>(stratlen));
      for (int i = 0; i < daylen; i++) {
        for (int j = 0; j < stratlen; j++) {
        exp_death_day_strata[i][j] = 0;
        }
      }
      for (int i = 0; i < daylen; i++) {
        for (int j = i+1; j < (daylen + 1); j++) {
          int delta = j - i - 1;
          for (int a = 0; a < stratlen; a++) {
            exp_death_day_strata[j-1][a] += infxns_strata[i][a] * ifr[a] * gmmlkup[delta];
          }
        }
      }
    // return as Rcpp list
    Rcpp::List ret = Rcpp::List::create(Rcpp::Named(\"exp_death_day_strata\") = exp_death_day_strata);
    return ret;
    }"
  # source cpp function
  Rcpp::cppFunction(src)


  # wrapper function for call cpp function
  draw_post_deaths_wrapper <- function(postdatsim,
                                       IFRmodel_inf){
    # run cpp function
    exp_death_day_strata <- get_post_deaths(infxns = unlist(postdatsim[, paste0("infxns_", IFRmodel_inf$inputs$IFRmodel$IFRparams)]),
                                            ifr = unlist(unique(postdatsim[, IFRmodel_inf$inputs$IFRmodel$IFRparams])),
                                            mod = unlist(unique(postdatsim[, IFRmodel_inf$inputs$IFRmodel$modparam])),
                                            sod = unlist(unique(postdatsim[, IFRmodel_inf$inputs$IFRmodel$sodparam])),
                                            stratlen = length(unique(postdatsim[, IFRmodel_inf$inputs$IFRmodel$IFRparams])),
                                            daylen = IFRmodel_inf$inputs$IFRmodel$maxObsDay)[[1]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("deaths_", IFRmodel_inf$inputs$IFRmodel$IFRparams)) %>%
      dplyr::mutate(time = 1:nrow(.)) %>%
      dplyr::select(c("time", dplyr::everything()))
    # out
    return(exp_death_day_strata)
  }

  #......................
  # Run draw post deaths for each simulation iteration
  #......................
  postdat.curves <- postdat %>%
    dplyr::select("sim", IFRmodel_inf$inputs$IFRmodel$IFRparams) %>%
    dplyr::filter(!duplicated(.)) %>%
    dplyr::mutate(
      post_deaths = purrr::map(.x = postdat.sims,
                               .f = draw_post_deaths_wrapper,
                               IFRmodel_inf = IFRmodel_inf)
    ) %>%
    tidyr::unnest(cols = post_deaths)

  return(postdat.curves)
}


#' @title Draw the Putative Observed Seroprevalnces using the Rogan-Gladen Estimator
#' @details Given sampling iterations with posterior-log-likes greater than or equal to a specific threshold, posterior results for the linear spline are generated. Assumed that the spline was fit in "un-transformed" space
#' @inheritParams get_cred_intervals
#' @param dwnsmpl integer; Number of posterior results to draw -- weighted by posterior prob
#' @importFrom magrittr %>%
#' @export

draw_posterior_RGobserved_seroprevalences <- function(IFRmodel_inf, whichrung = "rung1", dwnsmpl, by_chain = TRUE) {
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
                                                        reparamInfxn = FALSE,
                                                        reparamSpec = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
  # pull out pieces I need
  fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
  fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
  fitcurve_curve <- stringr::str_replace(fitcurve_curve, "if \\(cum_infxn_check <= popN\\) \\{", "")
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "double sero_loglik = 0.0;", n = 2)[,1]
  fitcurve_string <- paste(fitcurve_start, fitcurve_curve,
                           "std::vector<std::vector<double>> obs_prev_mat(n_sero_obs, std::vector<double>(stratlen));
                            for (int i = 0; i < n_sero_obs; i++) {
                              for (int j = 0; j < stratlen; j++) {
                                obs_prev_mat[i][j] = (sero_con_num[i][j]/demog[j]) * (spec + (sens-1)) - (spec-1);
                              }
                            }",
                           "Rcpp::List ret = Rcpp::List::create(obs_prev_mat); return ret;}",
                           collapse = "")
  Rcpp::cppFunction(fitcurve_string)

  #......................
  # inputs needed for cpp function
  #......................
  misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
                   demog = IFRmodel_inf$inputs$IFRmodel$demog$popN,
                   rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
                   days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams)+1,
                   n_sero_obs = length(IFRmodel_inf$inputs$IFRmodel$Serodayparams))

  datin <- list("obs_deaths" = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
                "obs_serology" = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroPrev)


  #......................
  # split, run, recombine
  #......................
  cpp_function_wrapper <- function(params, data, misc) {
    paramsin <- unlist(params[c(IFRmodel_inf$inputs$IFRmodel$modparam,
                                IFRmodel_inf$inputs$IFRmodel$sodparam,
                                IFRmodel_inf$inputs$IFRmodel$IFRparams,
                                IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                                IFRmodel_inf$inputs$IFRmodel$Knotparams,
                                IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                                IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                                IFRmodel_inf$inputs$IFRmodel$Noiseparams)])


    seroprev_strata <- loglike(params = paramsin,
                               param_i = 1,
                               data = datin,
                               misc = misc_list)[[1]]
    seroprev <- seroprev_strata %>%
      do.call("rbind.data.frame", .)

    colnames(seroprev) <- paste0("seroprev_", IFRmodel_inf$inputs$IFRmodel$IFRparams)
    ret <- cbind.data.frame(SeroDay = IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                            seroprev)
    return(ret)

  }

  mcmcout.node.rows <- split(mcmcout.nodes, 1:nrow(mcmcout.nodes))
  mcmcout.nodes$seroprev <- purrr::map(mcmcout.node.rows, cpp_function_wrapper,
                                       data = datin, misc = misc_list)

  #......................
  # tidy
  #......................
  # keep params around for convenience
  if (by_chain) {
    dat <- mcmcout.nodes %>%
      dplyr::select(c("chain",
                      IFRmodel_inf$inputs$IFRmodel$modparam,
                      IFRmodel_inf$inputs$IFRmodel$sodparam,
                      IFRmodel_inf$inputs$IFRmodel$IFRparams,
                      IFRmodel_inf$inputs$IFRmodel$Knotparams,
                      IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                      IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                      IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                      IFRmodel_inf$inputs$IFRmodel$Noiseparams,
                      "seroprev")) %>%
      dplyr::group_by(chain) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      dplyr::ungroup(chain) %>%
      tidyr::unnest(cols = "seroprev")


  } else {
    # keep params around for convenience
    dat <- mcmcout.nodes %>%
      dplyr::select(c(IFRmodel_inf$inputs$IFRmodel$modparam,
                      IFRmodel_inf$inputs$IFRmodel$sodparam,
                      IFRmodel_inf$inputs$IFRmodel$IFRparams,
                      IFRmodel_inf$inputs$IFRmodel$Knotparams,
                      IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                      IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                      IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                      IFRmodel_inf$inputs$IFRmodel$Noiseparams,
                      "seroprev")) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      tidyr::unnest(cols = "seroprev")
  }

  #......................
  # out
  #......................
  return(dat)
}


#' @title Draw the Inferred Seropevalence Curves both Adjusted and Unadjusted for Specificity and Sensitivity with the Rogan-Gladen Estimator
#' @details Given sampling iterations with posterior-log-likes greater than or equal to a specific threshold, posterior results for the linear spline are generated. Assumed that the spline was fit in "un-transformed" space
#' @inheritParams get_cred_intervals
#' @param dwnsmpl integer; Number of posterior results to draw -- weighted by posterior prob
#' @importFrom magrittr %>%
#' @export

draw_posterior_sero_curves <- function(IFRmodel_inf, whichrung = "rung1", dwnsmpl, by_chain = TRUE) {
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
                                                        reparamInfxn = FALSE,
                                                        reparamSpec = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
  # pull out pieces I need
  fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
  fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
  fitcurve_curve <- stringr::str_replace(fitcurve_curve, "if \\(cum_infxn_check <= popN\\) \\{", "")
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "std::vector<std::vector<double>> sero_con_num\\(n_sero_obs, std::vector<double>\\(stratlen\\)\\);", n = 2)[,1]
  fitcurve_string <- paste(fitcurve_start, fitcurve_curve,
                           "std::vector<std::vector<double>> full_sero_con_num(days_obsd, std::vector<double>(stratlen));
                            std::vector<std::vector<double>> RGfull_sero_con_num(days_obsd, std::vector<double>(stratlen));
                            // get cumulative hazard for study period
                            std::vector<double> cum_hazard(days_obsd);
                            for (int d = 0; d < days_obsd; d++) {
                              cum_hazard[d] = 1-exp((-(d+1)/sero_rate));
                            }

                            for (int i = 0; i < days_obsd; i++) {
                              for (int j = 0; j < stratlen; j++) {
                                // loop through and split infection curve by strata and by number of seroconversion study period
                                // note this cumulative, so loop through previous days
                                for (int d = 0; d < days_obsd; d++) {
                                  int time_elapsed = days_obsd - d - 1;
                                  full_sero_con_num[i][j] += infxn_spline[d] * ne[j] * cum_hazard[time_elapsed];
                                }
                              }
                            }

                            // correct for spec/sens
                            for (int i = 0; i < days_obsd; i++) {
                              for (int j = 0; j < stratlen; j++) {
                                // Rogan-Gladen Estimator
                                double obs_prev = (full_sero_con_num[i][j]/demog[j]) * (spec + (sens-1)) - (spec-1);
                                RGfull_sero_con_num[i][j] = round(obs_prev * demog[j]);
                              }
                            }",
                           "Rcpp::List ret = Rcpp::List::create(full_sero_con_num,  RGfull_sero_con_num); return ret;}",
                           collapse = "")
  Rcpp::cppFunction(fitcurve_string)

  #......................
  # inputs needed for cpp function
  #......................
  misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
                   demog = IFRmodel_inf$inputs$IFRmodel$demog$popN,
                   rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
                   days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams)+1,
                   n_sero_obs = length(IFRmodel_inf$inputs$IFRmodel$Serodayparams))

  datin <- list("obs_deaths" = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
                "obs_serology" = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroPrev)


  #......................
  # split, run, recombine
  #......................
  cpp_function_wrapper <- function(params, data, misc) {
    paramsin <- unlist(params[c(IFRmodel_inf$inputs$IFRmodel$modparam,
                                IFRmodel_inf$inputs$IFRmodel$sodparam,
                                IFRmodel_inf$inputs$IFRmodel$IFRparams,
                                IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                                IFRmodel_inf$inputs$IFRmodel$Knotparams,
                                IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                                IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                                IFRmodel_inf$inputs$IFRmodel$Noiseparams)])

    seroprev_lists <- loglike(params = paramsin,
                              param_i = 1,
                              data = datin,
                              misc = misc_list)

    inf_sero_con_num <- seroprev_lists[[1]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("inf_seroprev_", IFRmodel_inf$inputs$IFRmodel$IFRparams)) %>%
      dplyr::mutate(ObsDay = sort(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$ObsDay))) %>%
      dplyr::select(c("ObsDay", dplyr::everything()))

    RG_sero_con_num <- seroprev_lists[[2]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("RG_inf_seroprev_", IFRmodel_inf$inputs$IFRmodel$IFRparams)) %>%
      dplyr::mutate(ObsDay = sort(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$ObsDay))) %>%
      dplyr::select(c("ObsDay", dplyr::everything()))


    ret <- list(inf_sero_con_num = inf_sero_con_num,
                RG_sero_con_num = RG_sero_con_num)
    return(ret)

  }

  mcmcout.node.rows <- split(mcmcout.nodes, 1:nrow(mcmcout.nodes))
  mcmcout.nodes$seroprev <- purrr::map(mcmcout.node.rows, cpp_function_wrapper,
                                       data = datin, misc = misc_list)

  #......................
  # tidy
  #......................
  # keep params around for convenience
  if (by_chain) {
    dat <- mcmcout.nodes %>%
      dplyr::select(c("chain",
                      IFRmodel_inf$inputs$IFRmodel$modparam,
                      IFRmodel_inf$inputs$IFRmodel$sodparam,
                      IFRmodel_inf$inputs$IFRmodel$IFRparams,
                      IFRmodel_inf$inputs$IFRmodel$Knotparams,
                      IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                      IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                      IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                      IFRmodel_inf$inputs$IFRmodel$Noiseparams,
                      "seroprev")) %>%
      dplyr::group_by(chain) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      dplyr::ungroup(chain) %>%
      tidyr::unnest(cols = "seroprev")


  } else {
    # keep params around for convenience
    dat <- mcmcout.nodes %>%
      dplyr::select(c(IFRmodel_inf$inputs$IFRmodel$modparam,
                      IFRmodel_inf$inputs$IFRmodel$sodparam,
                      IFRmodel_inf$inputs$IFRmodel$IFRparams,
                      IFRmodel_inf$inputs$IFRmodel$Knotparams,
                      IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                      IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                      IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                      IFRmodel_inf$inputs$IFRmodel$Noiseparams,
                      "seroprev")) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      tidyr::unnest(cols = "seroprev")
  }

  #......................
  # out
  #......................
  return(dat)
}
