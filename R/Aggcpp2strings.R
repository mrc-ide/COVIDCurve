#' @title Create the logprior from the IFRmodel R6 Class for DrJacoby Inference
#' @param IFRmodel R6 class; Internal model object for COVIDCurve
#' @param reparamIFR logical; Whether IFRs should be reparameterized or inferred seperately
#' @param reparamKnots logical; Whether infection knots (i.e. the x-coordinates of the infection spline) should be reparameterized or inferred seperately
#' @param reparamInfxn logical; Whether infection curve (i.e. the  y-coordinates infection spline) should be reparameterized or inferred seperately
#' @param reparamDelays logical; Whether the numerous correlation serology paratmers should be reparameterized (mean offset-to-death is scaled by 1/specificity, attack rate noise vector is scaled by 1/specificity, and the seroconversion rate delay is recast as function of the mean offset-to-death) or inferred seperately
#' @param reparamNe logical; Whether "noise scalar effects" should be reparameterized or inferred seperately (if TRUE, considered relateve to Ne1)
#' @noRd

make_user_Agg_logprior <- function(IFRmodel, reparamIFR, reparamInfxn, reparamKnots, reparamDelays, reparamNe) {
  #..................
  # assertsions
  #..................
  assert_custom_class(IFRmodel, "IFRmodel")
  assert_logical(reparamIFR)
  assert_logical(reparamInfxn)
  assert_logical(reparamKnots)
  assert_logical(reparamDelays)
  assert_logical(reparamNe)

  #..................
  # setup
  #..................
  paramdf <- IFRmodel$paramdf
  IFRparams <- paramdf[paramdf$name %in% IFRmodel$IFRparams, ]
  Knotparams <- paramdf[paramdf$name %in% IFRmodel$Knotparams, ]
  Infxnparams <- paramdf[paramdf$name %in% IFRmodel$Infxnparams, ]
  Serotestparams <- paramdf[paramdf$name %in% IFRmodel$Serotestparams, ]
  Noiseparams <- paramdf[paramdf$name %in% IFRmodel$Noiseparams, ]
  TODparams <- paramdf[paramdf$name %in% c(IFRmodel$modparam, IFRmodel$sodparam), ]

  if (reparamKnots) {
    #..................
    # account for knot reparam -- xpos
    #..................
    assert_non_null(IFRmodel$relKnot, message = "Reparameterization requires relative knot to be indicated (i.e. relKnot)")
    relKnot <- IFRmodel$relKnot
    knotscalars <- Knotparams$name[Knotparams$name != relKnot]
  }

  if (reparamInfxn) {
    #..................
    # account for Infection Point reparam -- Ypos
    #..................
    assert_non_null(IFRmodel$relInfxn, message = "Reparameterization requires relative infection point to be indicated (i.e. relInfxn)")
    relInfxn <- IFRmodel$relInfxn
    infxnscalars <- Infxnparams$name[Infxnparams$name != relInfxn]
  }

  if (reparamIFR) {
    #..................
    # account for IFR reparam
    #..................
    assert_non_null(IFRmodel$maxMa, message = "Reparameterization requires expected max mortaltiy group to be indicated (i.e. maxMa)")
    maxMa <- IFRmodel$maxMa
    ifrscalars <- IFRparams$name[IFRparams$name != maxMa]
  }

  #..................
  # priors for knots -- Xpos
  #..................
  Knotextractparams <- sapply(Knotparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeknotpriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(", param, ",", d1, ",", d2, ",", "true) +")
  }, param = Knotparams$name, d1 = Knotparams$dsc1, d2 = Knotparams$dsc2)

  #..................
  # priors for infxnpts -- Ypos
  #..................
  Infxnextractparams <- sapply(Infxnparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeinfxnpriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(", param, ",", d1, ",", d2, ",", "true) +")
  }, param = Infxnparams$name, d1 = Infxnparams$dsc1, d2 = Infxnparams$dsc2)

  #..................
  # priors for IFRparams
  #..................
  IFRextractparams <- sapply(IFRparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeifrpriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = IFRparams$name, d1 = IFRparams$dsc1, d2 = IFRparams$dsc2)


  #..................
  # priors for Time Onset to Death
  #..................
  TODextractparams <- sapply(TODparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  modparams <- TODparams[grepl("mod", TODparams$name), ]
  maketodpriors_mod <- paste0("R::dnorm(", modparams$name, ",", modparams$dsc1, ",", modparams$dsc2, ",", "true) +")
  sodparams <- TODparams[grepl("sod", TODparams$name), ]
  maketodpriors_sod <- paste0("R::dbeta(", sodparams$name, ",", sodparams$dsc1, ",", sodparams$dsc2, ",", "true) +")
  maketodpriors <- c(maketodpriors_mod, maketodpriors_sod)

  #......................
  # Serology test priors
  #......................
  Serotestextractparams <- sapply(Serotestparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makeSerotestpriors <- c(
    paste0("R::dnorm(sero_rate,", Serotestparams$dsc1[Serotestparams$name == "sero_rate"], ",", Serotestparams$dsc2[Serotestparams$name == "sero_rate"], ", true) +"),
    paste0("R::dbeta(sens,", Serotestparams$dsc1[Serotestparams$name == "sens"], ",", Serotestparams$dsc2[Serotestparams$name == "sens"], ", true) +"),
    paste0("R::dbeta(spec,", Serotestparams$dsc1[Serotestparams$name == "spec"], ",", Serotestparams$dsc2[Serotestparams$name == "spec"], ", true) +")
  )


  #..................
  # priors for Noiseparams
  #..................
  Noiseextractparams <- sapply(Noiseparams$name, function(param){
    paste0("double ", param, " = params[\"",  param, "\"];")
  })

  makenoisepriors <- mapply(function(param, d1, d2){
    paste0("R::dunif(",param, ",", d1, ",", d2, ",", "true) +")
  }, param = Noiseparams$name, d1 = Noiseparams$dsc1, d2 = Noiseparams$dsc2)


  #..................
  # bring together
  #..................
  extractparams <- c(IFRextractparams, Knotextractparams, Infxnextractparams, Serotestextractparams, Noiseextractparams, TODextractparams)
  # account for reparam
  switch(paste0(reparamIFR, "-", reparamKnots, "-", reparamInfxn, "-", reparamDelays, "-", reparamNe),
         "TRUE-TRUE-TRUE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec) +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "TRUE-TRUE-TRUE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec);")
           )
         },

         "TRUE-TRUE-TRUE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "TRUE-TRUE-TRUE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ");")
           )
         },

         "TRUE-TRUE-FALSE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec) +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "TRUE-TRUE-FALSE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec);")
           )
         },

         "TRUE-TRUE-FALSE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "TRUE-TRUE-FALSE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(knotscalars), "*log(", relKnot, ");")
           )
         },

         "TRUE-FALSE-TRUE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec) +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "TRUE-FALSE-TRUE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec);")
           )
         },

         "TRUE-FALSE-TRUE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "TRUE-FALSE-TRUE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ");")
           )
         },

         "TRUE-FALSE-FALSE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec) +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "TRUE-FALSE-FALSE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec);")
           )
         },

         "TRUE-FALSE-FALSE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ") +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "TRUE-FALSE-FALSE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(ifrscalars), "*log(", maxMa, ");"))
         },

         "FALSE-TRUE-TRUE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec) +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "FALSE-TRUE-TRUE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec);")
           )
         },

         "FALSE-TRUE-TRUE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "FALSE-TRUE-TRUE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0(length(infxnscalars), "*log(", relInfxn, ");"),
           )
         },

         "FALSE-TRUE-FALSE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec) +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "FALSE-TRUE-FALSE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec);")
           )
         },

         "FALSE-TRUE-FALSE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(knotscalars), "*log(", relKnot, ") +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "FALSE-TRUE-FALSE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(knotscalars), "*log(", relKnot, ");")
           )
         },

         "FALSE-FALSE-TRUE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec) +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "FALSE-FALSE-TRUE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec);")
           )
         },

         "FALSE-FALSE-TRUE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(infxnscalars), "*log(", relInfxn, ") +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "FALSE-FALSE-TRUE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0(length(infxnscalars), "*log(", relInfxn, ");"),
           )
         },

         "FALSE-FALSE-FALSE-TRUE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec) +"),
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "FALSE-FALSE-FALSE-TRUE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0("log(", IFRmodel$modparam, ") + 2*log(1/spec);"),
           )
         },

         "FALSE-FALSE-FALSE-FALSE-TRUE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       paste0((length(Noiseextractparams)-1), "*log(", Noiseparams$name[1], ");")
           )
         },

         "FALSE-FALSE-FALSE-FALSE-FALSE" = {
           priors <- c("double ret =", makeifrpriors, makeknotpriors, makeinfxnpriors, makeSerotestpriors, makenoisepriors, maketodpriors,
                       ";")
         },


         {
           stop("Prior option not correctly specified during make priors")
         }
  )

  # out
  ret <- c("SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {",
           extractparams,
           priors,
           "if (!std::isfinite(ret)) { const double OVERFLO_DOUBLE = DBL_MAX/100.0; ret = -OVERFLO_DOUBLE; }",
           "return Rcpp::wrap(ret);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}



#' @title Create the loglikelihood from the IFRmodel R6 Class for DrJacoby Inference
#' @inheritParams make_user_Agg_logprior
#' @noRd

make_user_Agg_loglike <- function(IFRmodel, reparamIFR, reparamInfxn, reparamKnots, reparamDelays, reparamNe) {
  #..................
  # assertions
  #..................
  assert_custom_class(IFRmodel, "IFRmodel")
  assert_logical(reparamIFR)
  assert_logical(reparamInfxn)
  assert_logical(reparamKnots)
  assert_logical(reparamDelays)
  assert_logical(reparamNe)

  #..................
  # setup
  #..................
  paramdf <- IFRmodel$paramdf
  Knotparams <- IFRmodel$Knotparams
  Infxnparams <- IFRmodel$Infxnparams
  IFRparams <- IFRmodel$IFRparams
  Noiseparams <- IFRmodel$Noiseparams

  #..................
  # extract misc and storage items
  #..................
  extmisc <- "std::vector<double> rho = Rcpp::as< std::vector<double> >(misc[\"rho\"]); int stratlen = rho.size(); std::vector<int> demog = Rcpp::as< std::vector<int> >(misc[\"demog\"]); int rcensor_day = misc[\"rcensor_day\"]; int days_obsd = misc[\"days_obsd\"]; int n_knots = misc[\"n_knots\"]; int n_sero_obs = misc[\"n_sero_obs\"]; std::vector<int> sero_survey_start = Rcpp::as< std::vector<int> >(misc[\"sero_survey_start\"]); std::vector<int> sero_survey_end = Rcpp::as< std::vector<int> >(misc[\"sero_survey_end\"]); int max_seroday_obsd = misc[\"max_seroday_obsd\"];"
  storageitems <- "std::vector<double> node_x(n_knots); std::vector<double> node_y(n_knots); std::vector<double>ma(stratlen); std::vector<double>ne(stratlen);"

  #......................
  # extract fixed parameters
  #......................
  fixedparams <- "double sens = params[\"sens\"]; double spec = params[\"spec\"]; double sero_rate = params[\"sero_rate\"]; double mod = params[\"mod\"]; double sod = params[\"sod\"];"

  #..................
  # Extract and potentially liftover knotreparam vars for Knots -- Infxn Xpositions
  # NB, "raw" here because we take in a double and need to convert it to an integer day later
  #.................
  if (reparamKnots) {
    assert_non_null(IFRmodel$relKnot, message = "Reparameterization requires relative knot to be indicated (i.e. relKnot)")
    knotparamdf <- paramdf[paramdf$name %in% Knotparams, ]
    relKnot <- IFRmodel$relKnot
    knotscalars <- knotparamdf$name[knotparamdf$name != relKnot]
    node_xvec <- rep(NA, length(Knotparams))
    # determine which relative position in our knot vector
    relnodex_pos <- which(relKnot == Knotparams)
    nodex_counter <- 1
    for (i in 1:length(node_xvec)) {
      if (i == relnodex_pos) {
        node_xvec[i] <- paste0("node_x[", i, "] = params[\"", relKnot, "\"];")
      } else {
        node_xvec[i] <- paste0("node_x[", i, "] = params[\"", knotscalars[nodex_counter], "\"] * params[\"", relKnot, "\"];")
        nodex_counter <- nodex_counter + 1
      }
    }
  } else {
    node_xvec <- rep(NA, length(Knotparams))
    for (i in 1:length(Knotparams)){
      node_xvec[i] <- paste0("node_x[", i, "] = params[\"", Knotparams[i], "\"];")
    }
  }
  # account for internal knot at position 1
  node_xvec <- c("node_x[0] = 1.0;", node_xvec)

  #..................
  # Extract and potentially liftover infxnreparam vars to Infxn Ypositions
  #.................
  if (reparamInfxn) {
    assert_non_null(IFRmodel$relInfxn, message = "Reparameterization requires relative infection points to be indicated (i.e. relInfxn)")
    infxnparamdf <- paramdf[paramdf$name %in% Infxnparams, ]
    relInfxn <- IFRmodel$relInfxn
    infxnscalars <- infxnparamdf$name[infxnparamdf$name != relInfxn]
    node_yvec <- rep(NA, length(Infxnparams))
    # determine which relative position in our node-y vector
    relnodey_pos <- which(relInfxn == Infxnparams)
    nodey_counter <- 1
    for (i in 1:length(node_yvec)) {
      if (i == relnodey_pos) {
        node_yvec[i] <- paste0("node_y[", i-1, "] = params[\"", relInfxn, "\"];")
      } else {
        node_yvec[i] <- paste0("node_y[", i-1, "] = params[\"", infxnscalars[nodey_counter], "\"] * params[\"", relInfxn, "\"];")
        nodey_counter <- nodey_counter + 1
      }
    }
  } else {
    node_yvec <- rep(NA, length(Infxnparams))
    for (i in 1:length(Infxnparams)){
      node_yvec[i] <- paste0("node_y[", i-1, "] = params[\"", Infxnparams[i],  "\"];")
    }
  }

  #..................
  # Extract and liftover ifrreparam vars to Mas
  #.................
  if (reparamIFR) {
    assert_non_null(IFRmodel$maxMa, message = "Reparameterization requires expected max mortaltiy group to be indicated (i.e. maxMa)")
    ifrparamdf <- paramdf[paramdf$name %in% IFRparams, ]
    maxMa <- IFRmodel$maxMa
    ifrscalars <- ifrparamdf$name[ifrparamdf$name != maxMa]
    mavec <- rep(NA, length(IFRparams))
    # determine which relative position in our ma vector
    mamax_pos <- which(maxMa == IFRparams)
    ma_counter <- 1
    for (i in 1:length(mavec)) {
      if (i == mamax_pos) {
        mavec[i] <- paste0("ma[", i-1, "] = params[\"", maxMa, "\"];")
      } else {
        mavec[i] <- paste0("ma[", i-1, "] = params[\"", ifrscalars[ma_counter], "\"] * params[\"", maxMa, "\"];")
        ma_counter <- ma_counter + 1
      }
    }
  } else {
    mavec <- rep(NA, length(IFRparams))
    for (i in 1:length(IFRparams)){
      mavec[i] <- paste0("ma[", i-1, "] = params[\"", IFRparams[i], "\"];")
    }
  }

  #......................
  # extract out noise params and potentially liftover
  #......................
  if (reparamNe) {
    noisevec_p1 <- paste0("ne[0] = params[\"",  Noiseparams[1], "\"];")
    noisevec_p2 <- mapply(function(x, y){
      paste0("ne", "[", y-1, "] ", " = params[\"",  x, "\"] * ne[0];")
    }, x = paramdf$name[paramdf$name %in% Noiseparams[2:length(Noiseparams)]], y = 2:length(Noiseparams))
    # bring together
    noisevec <- c(noisevec_p1, noisevec_p2)

  } else {
    noisevec <- mapply(function(x, y){
      paste0("ne", "[", y-1, "] ", " = params[\"",  x, "\"];")
    }, x = paramdf$name[paramdf$name %in% Noiseparams[1:length(Noiseparams)]], y = 1:length(Noiseparams))
  }



  #......................
  # liftover for highly correlated serology parameters, including:
  #  mean offset (gamma) for sero rate (exp) and Spec
  #......................
  if (reparamDelays) {
    seroscorrlftovr <- "sero_rate = sero_rate * 1/spec;"
  } else {
    seroscorrlftovr <- ""
  }


  #..................
  # get loglike
  #..................
  loglike <- readLines(system.file("covidcurve", "natcubicspline_loglike.cpp", package = "COVIDCurve", mustWork = TRUE))
  loglike_start <- grep("// Lookup Items", loglike)
  loglike_end <- grep("// return as Rcpp list", loglike)
  loglike <- loglike[loglike_start:loglike_end]
  # remove comments which cause issues when coercing to string in this format
  commlines <- grep("//", loglike)
  loglike <- loglike[! 1:length(loglike) %in% commlines]
  loglike <- gsub("\\\n", "", loglike)

  # end loglike and now out
  ret <- c("SEXP loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {",
           extmisc,
           storageitems,
           fixedparams,
           node_xvec,
           node_yvec,
           mavec,
           noisevec,
           seroscorrlftovr,
           loglike,
           "return Rcpp::wrap(loglik);",
           "}"
  )

  out <- capture.output(cat(ret))
  return(out)
}


