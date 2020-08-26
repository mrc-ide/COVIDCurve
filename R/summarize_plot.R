#' IFR-Inference-Aggregate-Model, overload summary
#' @noRd
#' @export

summary.IFRmodel <- function(object, ...){
  cat(crayon::cyan("*~*~*~*~ IFR Inference Model ~*~*~*~*\n"))
  cat(crayon::green("IFR strata params: "), paste(object$IFRparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Point Params: "), paste(object$Infxnparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Knot Params: "), paste(object$Knotparams, collapse = ", "), "\n")
  cat(crayon::magenta("Serology Test Parameters: "), paste(object$Serotestparams, collapse = ", "), "\n")
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
  cat(crayon::yellow("Total Population Size: "), sum(x$demog$popN), "\n")
}

#' @title Get Credible Intervals for Parameters from Sampling Iterations
#' @param IFRmodel_inf R6 class; The result of the IFR Model MCMC run along with the model input.
#' @param what character; Specify which parameters you want: "IFRparams", "Infxnparams", "Serotestparams", or "Noiseparams"
#' @param whichrung character; Specify which rung to sample from (default is rung1)
#' @param by_chain logical; Whether or not credible intervals should be reported with respect to individual chains (TRUE) or not.
#' @importFrom magrittr %>%
#' @export

get_cred_intervals <- function(IFRmodel_inf, what, whichrung = "rung1", by_chain = TRUE) {
  assert_custom_class(IFRmodel_inf$inputs$IFRmodel, "IFRmodel")
  assert_custom_class(IFRmodel_inf$mcmcout, "drjacoby_output")
  assert_custom_class(IFRmodel_inf, "IFRmodel_inf")
  assert_in(what, c("IFRparams", "Knotparams", "Infxnparams", "Serotestparams", "Noiseparams"))
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

         "Noiseparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel_inf$inputs$IFRmodel$Noiseparams)
         },
         "Noiseparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel_inf$inputs$IFRmodel$Noiseparams)
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




#' @title Draw posterior results from the Infection Curve (Cubic Spline)
#' @details Given sampling iterations with posterior-log-likes greater than or equal to a specific threshold, posterior results for the linear spline are generated. Assumed that the spline was fit in "un-transformed" space
#' @inheritParams get_cred_intervals
#' @param dwnsmpl integer; Number of posterior results to draw -- weighted by posterior prob
#' @importFrom magrittr %>%
#' @export

draw_posterior_infxn_cubic_splines <- function(IFRmodel_inf, whichrung = "rung1", dwnsmpl,
                                               by_chain = TRUE, by_strata = FALSE) {
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
                                                        reparamInfxn = FALSE,
                                                        reparamDelays = FALSE,
                                                        reparamNe = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
  # pull out pieces I need
  fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
  fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "bool popN_pass = true;", n = 2)[,1]
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
  # misc list
  misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
                   rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
                   days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams) + 1, # +1 because we set an internal knot for pos 1
                   n_sero_obs = length(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroStartSurvey)),
                   sero_survey_start = unique(IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroStartSurvey),
                   sero_survey_end = unique(IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroEndSurvey),
                   max_seroday_obsd = max(IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroEndSurvey),
                   demog = IFRmodel_inf$inputs$IFRmodel$demog$popN)
  # data in
  datin <- list(obs_deaths = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
                prop_strata_obs_deaths = IFRmodel_inf$inputs$IFRmodel$data$prop_deaths$PropDeaths,
                obs_serologypos = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroPos,
                obs_serologyn = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroN)

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
  switch(paste0(by_chain, "-", by_strata),
         # by chain and by strata
         "TRUE-TRUE"={
           plotdat <- mcmcout.nodes %>%
             dplyr::select(c("chain",
                             IFRmodel_inf$inputs$IFRmodel$modparam,
                             IFRmodel_inf$inputs$IFRmodel$sodparam,
                             IFRmodel_inf$inputs$IFRmodel$IFRparams,
                             IFRmodel_inf$inputs$IFRmodel$Knotparams,
                             IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                             IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                             IFRmodel_inf$inputs$IFRmodel$Noiseparams,
                             "infxncurves")) %>%
             dplyr::group_by(chain) %>%
             dplyr::mutate(sim = 1:dplyr::n()) %>%
             dplyr::ungroup(chain) %>%
             tidyr::unnest(cols = "infxncurves")
           # downsize for what is needed in ggplot vs. returrn
           plotdat_sm <- plotdat %>%
             dplyr::select(c("chain", "time", "sim", dplyr::starts_with("infxns_"))) %>%
             magrittr::set_colnames(gsub("infxns_", "", colnames(.))) %>%
             tidyr::gather(., key = "Strata", value = "infxns", 4:ncol(.))
           if( !is.null(IFRmodel_inf$inputs$IFRmodel$IFRdictkey) ) {
             # set up dictionary key to be ageband, regions, etc.
             IFRdictkey <- IFRmodel_inf$inputs$IFRmodel$IFRdictkey
             colnames(IFRdictkey)[which(colnames(IFRdictkey) != "Strata")] <- "stratavar"
             plotdat_sm <- plotdat_sm %>%
               dplyr::left_join(., IFRdictkey, by = "Strata")
           } else {
             plotdat_sm <- plotdat_sm %>%
               dplyr::rename(stratavar = Strata)
           }


           plotObj <- ggplot2::ggplot() +
             ggplot2::geom_line(data = plotdat_sm,
                                mapping =  ggplot2::aes(x = time, y = infxns, group = sim), alpha = 0.25,
                                lwd = 0.5, color = "#d9d9d9") +
             ggplot2::xlab("Time") +  ggplot2::ylab("Num. Infxns")  +
             ggplot2::labs(title = "Posterior Draws of the Infection Curve") +
             ggplot2::facet_wrap(stratavar ~ chain) +
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
         },
         # by chain but not by strata
         "TRUE-FALSE"={
           plotdat <- mcmcout.nodes %>%
             plotdat <- mcmcout.nodes %>%
               dplyr::select(c("chain",
                               IFRmodel_inf$inputs$IFRmodel$modparam,
                               IFRmodel_inf$inputs$IFRmodel$sodparam,
                               IFRmodel_inf$inputs$IFRmodel$IFRparams,
                               IFRmodel_inf$inputs$IFRmodel$Knotparams,
                               IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                               IFRmodel_inf$inputs$IFRmodel$Serotestparams,
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
         },
         # not by chain but by strata
         "FALSE-TRUE"={
           plotdat <- mcmcout.nodes %>%
             dplyr::select(c("chain",
                             IFRmodel_inf$inputs$IFRmodel$modparam,
                             IFRmodel_inf$inputs$IFRmodel$sodparam,
                             IFRmodel_inf$inputs$IFRmodel$IFRparams,
                             IFRmodel_inf$inputs$IFRmodel$Knotparams,
                             IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                             IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                             IFRmodel_inf$inputs$IFRmodel$Noiseparams,
                             "infxncurves")) %>%
             dplyr::group_by(chain) %>%
             dplyr::mutate(sim = 1:dplyr::n()) %>%
             dplyr::ungroup(chain) %>%
             tidyr::unnest(cols = "infxncurves")
           # downsize for what is needed in ggplot vs. returrn
           plotdat_sm <- plotdat %>%
             dplyr::select(c("time", "sim", dplyr::starts_with("infxns_"))) %>%
             magrittr::set_colnames(gsub("infxns_", "", colnames(.))) %>%
             tidyr::gather(., key = "Strata", value = "infxns", 3:ncol(.))
           if( !is.null(IFRmodel_inf$inputs$IFRmodel$IFRdictkey) ) {
             # set up dictionary key to be ageband, regions, etc.
             IFRdictkey <- IFRmodel_inf$inputs$IFRmodel$IFRdictkey
             colnames(IFRdictkey)[which(colnames(IFRdictkey) != "Strata")] <- "stratavar"
             plotdat_sm <- plotdat_sm %>%
               dplyr::left_join(., IFRdictkey, by = "Strata")
           } else {
             plotdat_sm <- plotdat_sm %>%
               dplyr::rename(stratavar = Strata)
           }

           plotObj <- ggplot2::ggplot() +
             ggplot2::geom_line(data = plotdat_sm,
                                mapping =  ggplot2::aes(x = time, y = infxns, group = sim), alpha = 0.25,
                                lwd = 0.5, color = "#d9d9d9") +
             ggplot2::xlab("Time") +  ggplot2::ylab("Num. Infxns")  +
             ggplot2::labs(title = "Posterior Draws of the Infection Curve") +
             ggplot2::facet_wrap(stratavar ~ .) +
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

         },

         # not by chain, not by strata
         "FALSE-FALSE"={
           plotdat <- mcmcout.nodes %>%
             dplyr::select(c(IFRmodel_inf$inputs$IFRmodel$modparam,
                             IFRmodel_inf$inputs$IFRmodel$sodparam,
                             IFRmodel_inf$inputs$IFRmodel$IFRparams,
                             IFRmodel_inf$inputs$IFRmodel$Knotparams,
                             IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                             IFRmodel_inf$inputs$IFRmodel$Serotestparams,
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

         })

  # check to see if censoring was applied
  if (IFRmodel_inf$inputs$IFRmodel$rcensor_day < IFRmodel_inf$inputs$IFRmodel$maxObsDay) {
    plotObj <- plotObj +
      ggplot2::geom_vline(xintercept = IFRmodel_inf$inputs$IFRmodel$rcensor_day, linetype = "dashed", size = 1.2, color = "#de2d26")
  }

  #......................
  # out
  #......................
  ret <- list(
    curvedata = plotdat,
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
                                                            dwnsmpl = dwnsmpl, by_chain = by_chain)$curvedata
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



#' ' @title Draw the Inferred Seropevalence Curves both Adjusted and Unadjusted for Specificity and Sensitivity with the Rogan-Gladen Estimator
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
  # get serocurves
  #......................
  # internal function, liftover cpp likelihood to get infxn curve
  # NOTE, this is extremely sensitive to the placements of the Cpp source file and therefore, is not generalizable
  fitcurve_string <- COVIDCurve:::make_user_Agg_loglike(IFRmodel = IFRmodel_inf$inputs$IFRmodel,
                                                        reparamIFR = FALSE,
                                                        reparamKnots = FALSE,
                                                        reparamInfxn = FALSE,
                                                        reparamDelays = FALSE,
                                                        reparamNe = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
  # pull out pieces I need
  fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
  fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
  fitcurve_curve <- stringr::str_replace(fitcurve_curve, "  if \\(popN_pass\\) \\{", "")
  fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "std::vector<double> cum_hazard\\(max_seroday_obsd\\);", n = 2)[,1]

  # rewriting the sero_con_num_full vector here to be all days observed, not just serology days
  fitcurve_string <- paste(fitcurve_start, fitcurve_curve,
                           "std::vector<double> cum_hazard(days_obsd);
                           for (int d = 0; d < days_obsd; d++) {
                             cum_hazard[d] = 1-exp((-(d+1)/sero_rate));
                           }
                           std::vector<std::vector<double>> sero_con_num_full(days_obsd, std::vector<double>(stratlen));
                            for (int a = 0; a < stratlen; a++) {
                              for (int i = 0; i < days_obsd; i++) {
                                for (int j = i+1; j < (days_obsd + 1); j++) {
                                  int time_elapsed = j - i - 1;
                                  sero_con_num_full[j-1][a] += infxn_spline[i] * ne[a] * cum_hazard[time_elapsed];
                                }
                              }
                            }
                           std::vector<std::vector<double>> crude_seroprev(days_obsd, std::vector<double>(stratlen));
                           std::vector<std::vector<double>> RG_seroprev(days_obsd, std::vector<double>(stratlen));
                            for (int i = 0; i < days_obsd; i++) {
                              for (int j = 0; j < stratlen; j++) {
                                crude_seroprev[i][j] = (sero_con_num_full[i][j]/demog[j]);
                                RG_seroprev[i][j] = sens*crude_seroprev[i][j] + (1-spec)*(1-crude_seroprev[i][j]);
                              }
                            }",
                           "Rcpp::List ret = Rcpp::List::create(sero_con_num_full, crude_seroprev,  RG_seroprev); return ret;}",
                           collapse = "")
  Rcpp::cppFunction(fitcurve_string)

  #......................
  # inputs needed for cpp function
  #......................
  # misc list
  misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
                   rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
                   days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
                   n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams) + 1, # +1 because we set an internal knot for pos 1
                   n_sero_obs = length(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroStartSurvey)),
                   sero_survey_start = unique(IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroStartSurvey),
                   sero_survey_end = unique(IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroEndSurvey),
                   max_seroday_obsd = max(IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroEndSurvey),
                   demog = IFRmodel_inf$inputs$IFRmodel$demog$popN)
  # data in
  datin <- list(obs_deaths = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
                prop_strata_obs_deaths = IFRmodel_inf$inputs$IFRmodel$data$prop_deaths$PropDeaths,
                obs_serologypos = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroPos,
                obs_serologyn = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroN)

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
                                IFRmodel_inf$inputs$IFRmodel$Noiseparams)])

    seroprev_lists <- loglike(params = paramsin,
                              param_i = 1,
                              data = datin,
                              misc = misc_list)

    crude_seroprev <- seroprev_lists[[1]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("serocounts_", IFRmodel_inf$inputs$IFRmodel$IFRparams)) %>%
      dplyr::mutate(ObsDay = sort(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$ObsDay))) %>%
      dplyr::select(c("ObsDay", dplyr::everything()))

    crude_seroprev <- seroprev_lists[[2]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("crude_pd_seroprev_", IFRmodel_inf$inputs$IFRmodel$IFRparams)) %>%
      dplyr::mutate(ObsDay = sort(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$ObsDay))) %>%
      dplyr::select(c("ObsDay", dplyr::everything()))

    RG_seroprev <- seroprev_lists[[3]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("RG_pd_seroprev_", IFRmodel_inf$inputs$IFRmodel$IFRparams)) %>%
      dplyr::mutate(ObsDay = sort(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$ObsDay))) %>%
      dplyr::select(c("ObsDay", dplyr::everything()))


    ret <- dplyr::left_join(crude_seroprev, RG_seroprev, by = "ObsDay")
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
