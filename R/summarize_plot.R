#' IFR-Inference-Aggregate-Model, overload summary
#' @noRd
#' @export

summary.IFRmodel <- function(object, ...){
  cat(crayon::cyan("*~*~*~*~ IFR Inference Model ~*~*~*~*\n"))
  cat(crayon::green("IFR strata params: "), paste(round(object$IFRparams, 2), collapse = ", "), "\n")
  cat(crayon::blue("Infection Point Params: "), paste(round(object$Infxnparams, 2), collapse = ", "), "\n")
  cat(crayon::blue("Infection Knot Positions: "), paste(round(object$knots, 2), collapse = ", "), "\n")
  cat(crayon::magenta("Serology Parameters: "), paste(round(object$Seroparams, 2), collapse = ", "), "\n")
  cat(crayon::yellow("Model Type: "), object$level, "\n")
  cat(crayon::yellow("Total Population Size: "), object$popN, "\n")
  cat(crayon::yellow("Prob. of Infection Given Strata: "),  paste(round(object$pa, 2), collapse = ", "), "\n")
  cat(crayon::yellow("Mean Delay of Onset-to-Death: "), object$mod, "\n")
  cat(crayon::yellow("Coef. Var. Delay of Onset-to-Death: "), object$sod, "\n")
}


#' IFR-Inference-Aggregate-Model, overload print
#' @noRd
#' @export

print.IFRmodel <- function(object, ...){
  cat(crayon::cyan("*~*~*~*~ IFR Inference Model ~*~*~*~*\n"))
  cat(crayon::green("IFR strata params: "), paste(object$IFRparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Point Params: "), paste(object$Infxnparams, collapse = ", "), "\n")
  cat(crayon::blue("Infection Knot Positions: "), paste(object$knots, collapse = ", "), "\n")
  cat(crayon::magenta("Serology Parameters: "), paste(object$Seroparams, collapse = ", "), "\n")
  cat(crayon::yellow("Model Type: "), object$level, "\n")
  cat(crayon::yellow("Total Population Size: "), object$popN, "\n")
  cat(crayon::yellow("Prob. of Infection Given Strata: "),  paste(round(object$pa, 2), collapse = ", "), "\n")
  cat(crayon::yellow("Mean Delay of Onset-to-Death: "), object$mod, "\n")
  cat(crayon::yellow("Coef. Var. Delay of Onset-to-Death: "), object$sod, "\n")
}

#' @title Get Credible Intervals for Parameters from Sampling Iterations
#' @param IFRmodel R6 class; Internal model object for COVIDCurve
#' @param mcmcout IFR Model MCMC Output; The result of the IFR Model MCMC run. Object will inherits classes from COVIDCurve ("IFRmodel_inf") and DrJacoby ("drjacoby_output")
#' @param what character; Specify which parameters you want: "IFRparams", "Infxnparams", or "Seroparams"
#' @param whichrung character; Specify which rung to sample from (default is rung1)
#' @param by_chain logical; Whether or not credible intervals should be reported with respect to individual chains (TRUE) or not.
#' @importFrom magrittr %>%
#' @export

get_cred_intervals <- function(IFRmodel, mcmcout, what, whichrung = "rung1", by_chain = TRUE) {
  assert_custom_class(IFRmodel, "IFRmodel")
  assert_custom_class(mcmcout, "IFRmodel_inf")
  assert_in(what, c("IFRparams", "Infxnparams", "Seroparams"))
  assert_logical(by_chain)

  # grouping vars
  switch(paste0(what, "-", by_chain),
         "IFRparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel$IFRparams)
         },
         "IFRparams-FALSE"={
           groupingvar <- "param"
           params <-  c("iteration", IFRmodel$IFRparams)
         },
         "Infxnparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel$Infxnparams)
         },
         "Infxnparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel$Infxnparams)
         },

         "Seroparams-TRUE"={
           groupingvar <- c("chain", "param")
           params <- c("chain", IFRmodel$Seroparams)
         },
         "Seroparams-FALSE"={
           groupingvar <- "param"
           params <- c("iteration", IFRmodel$Seroparams)
         }
  )

  mcmcout$output %>%
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
#' @param IFRmodel R6 class; Internal model object for COVIDCurve
#' @param mcmcout IFR Model MCMC Output; The result of the IFR Model MCMC run. Object will inherits classes from COVIDCurve ("IFRmodel_inf") and DrJacoby ("drjacoby_output")
#' @param by_chain logical; Whether or not credible intervals should be reported with respect to individual chains (TRUE) or not.
#' @param whichrung character; Specify which rung to sample from (default is rung1)
#' @param CIquant numeric; Quantile to draw from the posterior curve
#' @importFrom magrittr %>%
#' @export

draw_posterior_infxn_points_cubic_splines <- function(IFRmodel, mcmcout, whichrung = "rung1", CIquant, by_chain = TRUE) {
  assert_custom_class(IFRmodel, "IFRmodel")
  assert_custom_class(mcmcout, "IFRmodel_inf")
  assert_numeric(CIquant)
  assert_bounded(CIquant, left = 0, right = 1)
  assert_logical(by_chain)
  #......................
  # fitler to sampling and by rung
  #......................
  mcmcout.nodes <-  mcmcout$output %>%
    dplyr::filter(stage == "sampling" & rung == whichrung)

  #......................
  # get natural cubic gradients
  #......................
  # internal functions, not generalizable
  # inputs
  knots <- IFRmodel$knots
  denom <- (knots - dplyr::lag(knots))[2:(length(knots))]

  # function
  make_cubic_infxn_curve <- function(mcmcout_row, IFRmodel, knots, denom){
    # get y param for this row
    node_y <- unlist( mcmcout_row[, IFRmodel$Infxnparams] )
    # get m
    m <- rep(NA, length(knots)-2)
    for (i in 2:(length(m)+1)) {
      m[i-1] <- (3/denom[i])*(node_y[i+1] - node_y[i]) - (3/denom[i-1])*(node_y[i] - node_y[i-1]);
    }

    # get z, g, k
    z <- rep(0, length(knots)-1)
    g <- rep(1, length(knots)-1)
    k <- rep(0, length(knots)-1)
    for (i in 2:(length(knots)-2)) {
      g[i] = 2*(knots[i+1] - knots[i-1]) - (denom[i-1])*(k[i-1])
      k[i] = denom[i]/g[i]
      z[i] = (m[i-1] - denom[i-1]*z[i-1])/g[i]
    }

    # finally loop through and get our "slopes" for our interpolants
    sp3 <- sp1 <- rep(NA, length(knots)-1)
    sp2 <- rep(0, length(knots))
    for(i in (length(knots)-1):1) {
      sp2[i] = z[i] - k[i]*sp2[i+1]
      sp1[i] = (node_y[i+1] - node_y[i])/(denom[i]) - (denom[i]*(sp2[i+1] + 2*sp2[i]))/3
      sp3[i] = (sp2[i+1] - sp2[i])/(3*denom[i])
    }

    # now recreate infection spline
    infxn_spline <- rep(NA, knots[length(knots)]- knots[1] + 1)
    infxn_spline[1] <- node_y[1]
    node_j <- 1
    for (i in 2:length(infxn_spline)) {
      infxn_spline[i] = node_y[node_j] +
        sp1[node_j] * (i - knots[node_j]) +
        sp2[node_j] * (i - knots[node_j])^2 +
        sp3[node_j] * (i - knots[node_j])^3

      # update node j if needed
      if ((knots[1] + i - 1) >= knots[node_j + 1]) {
        node_j <- node_j +1
      }
    }
    out <- data.frame(time = 1:length(infxn_spline),
                      infxns = infxn_spline)
    return(out)
  }

  #......................
  # sample by CI limit and make infxn curves
  #......................
  mcmcout.nodes <- mcmcout.nodes %>%
    dplyr::mutate(logposterior = loglikelihood + logprior)

  # filter
  upperci <- quantile(mcmcout.nodes$logposterior,
                      probs = CIquant)
  mcmcout.nodes <- mcmcout.nodes %>%
    dplyr::filter(logposterior >= upperci)

  # split, run, recombine
  mcmcout.node.rows <- split(mcmcout.nodes, 1:nrow(mcmcout.nodes))
  mcmcout.nodes$infxncurves <- furrr::future_map(mcmcout.node.rows, make_cubic_infxn_curve,
                                                 knots = knots, denom = denom, IFRmodel = IFRmodel)

  #......................
  # tidy
  #......................
  # keep IFR params around for convenience
  if (by_chain) {
    plotdat <- mcmcout.nodes %>%
      dplyr::select(c("chain", IFRmodel$IFRparams, "infxncurves")) %>%
      dplyr::group_by(chain) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      dplyr::ungroup(chain) %>%
      tidyr::unnest(cols = "infxncurves")

    plotObj <- ggplot2::ggplot() +
      ggplot2::geom_line(data = plotdat, mapping =  ggplot2::aes(time, infxns, group = sim), alpha = 0.25,
                         lwd = 0.5, color = "#d9d9d9") +
      ggplot2::geom_vline(xintercept = IFRmodel$knots, color = "#cb181d", lwd = 0.25, linetype = "dashed", alpha = 0.5) +
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


  } else {
    # keep IFR params around for convenience
    plotdat <- mcmcout.nodes %>%
      dplyr::select(c(IFRmodel$IFRparams, "infxncurves")) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      tidyr::unnest(cols = "infxncurves")

    plotObj <-  ggplot2::ggplot() +
      ggplot2::geom_line(data = plotdat, mapping =  ggplot2::aes(time, infxns, group = sim), alpha = 0.25,
                         lwd = 0.5, color = "#d9d9d9") +
      ggplot2::geom_vline(xintercept = IFRmodel$knots, color = "#cb181d", lwd = 0.25, linetype = "dashed", alpha = 0.5) +
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

posterior_check_infxns_to_death <- function(IFRmodel, mcmcout, whichrung = "rung1",
                                            CIquant, by_chain = FALSE) {
  postdat <- COVIDCurve::draw_posterior_infxn_points_cubic_splines(IFRmodel, mcmcout, whichrung = whichrung,
                                                                   CIquant, by_chain = by_chain)$plotdat
  # set up function to draw posterior deaths
  postdat.sims <- split(postdat, factor(postdat$sim))
  draw_post_deaths <- function(postdatsim){
    gmmlkup <- stats::dgamma(postdatsim$time, shape = 1/(IFRmodel$sod^2), scale = IFRmodel$mod*IFRmodel$sod^2, log = F)

    # exp deaths day and strata
    exp_death.day <- rep(0, length = max(IFRmodel$knots) - IFRmodel$knots[1] + 1)
    # spread infxns out to day when they may or may not die
    for (i in 1:nrow(postdatsim)) {
      for (j in (i+1):(nrow(postdatsim) + 1)) {
        delta <- j - i
        exp_death.day[j-1] <- exp_death.day[j-1] + postdatsim$infxns[i] * gmmlkup[delta]
      }
    }

    exp_death.day.strata <- exp_death.day %*% t(IFRmodel$pa) # account for pa
    exp_death.day.strata <- exp_death.day.strata * postdatsim[, IFRmodel$IFRparams]

    #......................
    # tidy up and out
    #......................
    out <- cbind.data.frame(time = 1:nrow(exp_death.day.strata), exp_death.day.strata)
    colnames(out)[2:ncol(out)] <- paste0("deaths_", IFRmodel$IFRparams)
    return(out)
  }
  # get post deaths
  postdat.curves <- postdat %>%
    dplyr::select("sim", IFRmodel$IFRparams) %>%
    dplyr::filter(!duplicated(.)) %>%
    dplyr::mutate(
      post_deaths = furrr::future_map(postdat.sims, draw_post_deaths)
    ) %>%
    tidyr::unnest(cols = post_deaths)

  return(postdat.curves)
}





