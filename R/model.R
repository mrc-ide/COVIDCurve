#' @title Make Model Inference Object-Oriented Class
#' @description A simple R6 class with only public features. Purpose is to provide a model framework
#' for inference of IFRs from Age-Specific Death Data.
#' @name make_IFRmodel_age
#' @export
NULL

#..............................................................
# declare
#..............................................................
make_IFRmodel_age <- R6::R6Class(classname = "IFRmodel",
                                 public = list(
                                   data = NULL,
                                   maxObsDay = NULL,
                                   modparam = NULL,
                                   sodparam = NULL,
                                   IFRparams = NULL,
                                   maxMa = NULL,
                                   Knotparams = NULL,
                                   relKnot = NULL,
                                   Infxnparams = NULL,
                                   relInfxn = NULL,
                                   Noiseparams = NULL,
                                   paramdf = NULL,
                                   Serotestparams = NULL,
                                   demog = NULL,
                                   rcensor_day = NULL,
                                   IFRdictkey = NULL,

                                   initialize = function(data = NULL, maxObsDay = NULL,
                                                         modparam = NULL, sodparam = NULL,
                                                         IFRparams = NULL, maxMa = NULL,
                                                         Infxnparams = NULL, relInfxn = NULL,
                                                         Knotparams = NULL, relKnot = NULL,
                                                         Serotestparams = NULL,
                                                         Noiseparams = NULL,
                                                         paramdf = NULL,
                                                         demog = NULL,
                                                         rcensor_day = NULL, IFRdictkey = NULL
                                   ) {
                                     #......................
                                     # assertions and checks
                                     #......................
                                     items <- c(data, IFRparams, Infxnparams, Knotparams, Serotestparams, Noiseparams, paramdf, modparam, sodparam, demog)
                                     if ( !all(sapply(items, is.null)) ) { # all user to input all options at init, otherwise initialize empty
                                       # assert data
                                       assert_list(data)
                                       assert_in(names(data), c("obs_deaths", "prop_deaths", "obs_serology"))
                                       # L1
                                       assert_dataframe(data$obs_deaths)
                                       assert_in(x = colnames(data$obs_deaths), y = c("ObsDay", "Deaths"))
                                       assert_numeric(data$obs_deaths$ObsDay)
                                       assert_increasing(data$obs_deaths$ObsDay)
                                       assert_numeric(data$obs_deaths$Deaths)
                                       # L2
                                       assert_dataframe(data$prop_deaths)
                                       assert_in(x = colnames(data$prop_deaths), y = c("Strata", "PropDeaths"))
                                       assert_bounded(data$prop_deaths$PropDeaths, left = 0, right = 1)
                                       assert_eq(sum(data$prop_deaths$PropDeaths), 1,
                                                 message = c("Proportion of deaths by strata must sum to 1"))
                                       assert_in(data$prop_deaths$Strata, IFRparams)

                                       # L3
                                       assert_dataframe(data$obs_serology)
                                       assert_in(colnames(data$obs_serology), c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev", "SeroLCI", "SeroUCI"))
                                       assert_in(data$obs_serology$Strata, IFRparams)
                                       sapply(data$obs_serology$SeroPos, assert_pos_int_or_NA)
                                       sapply(data$obs_serology$SeroN, assert_pos_int_or_NA)
                                       sapply(data$obs_serology$SeroPrev, assert_bounded_or_NA, left = 0, right = 1)
                                       sapply(data$obs_serology$SeroLCI, assert_numeric_or_NA)
                                       sapply(data$obs_serology$SeroUCI, assert_numeric_or_NA)
                                       assert_pos_int(data$obs_serology$SeroStartSurvey)
                                       assert_pos_int(data$obs_serology$SeroEndSurvey)

                                       #assert params
                                       assert_string(modparam)
                                       assert_string(sodparam)
                                       assert_string(IFRparams)
                                       assert_unique(IFRparams)
                                       assert_string(Infxnparams)
                                       assert_unique(Infxnparams)
                                       assert_string(Knotparams)
                                       assert_unique(Knotparams)
                                       assert_string(Serotestparams)
                                       assert_unique(Serotestparams)
                                       assert_in(Serotestparams, c("sens", "spec", "sero_con_rate", "sero_rev_shape", "sero_rev_scale"))

                                       # consider Ne if multiple strata
                                       if (length(IFRparams) > 1) {
                                         assert_string(Noiseparams)
                                         assert_unique(Noiseparams)
                                         assert_same_length(Noiseparams, IFRparams)
                                       }


                                       # assert paramdf
                                       assert_dataframe(paramdf)
                                       assert_in(x = colnames(paramdf), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                       assert_string(paramdf$name)
                                       assert_unique(paramdf$name)

                                       # consider Ne if multiple strata
                                       if (length(IFRparams) > 1) {
                                       assert_in(paramdf$name, c(modparam, sodparam, IFRparams, Infxnparams, Knotparams, Serotestparams, Noiseparams))
                                       assert_in(c(modparam, sodparam, IFRparams, Infxnparams, Knotparams, Serotestparams, Noiseparams), paramdf$name)
                                       } else {
                                         assert_in(paramdf$name, c(modparam, sodparam, IFRparams, Infxnparams, Knotparams, Serotestparams))
                                         assert_in(c(modparam, sodparam, IFRparams, Infxnparams, Knotparams, Serotestparams), paramdf$name)
                                       }

                                       sapply(paramdf$min, assert_numeric_or_NA)
                                       sapply(paramdf$init, assert_numeric_or_NA)
                                       sapply(paramdf$max, assert_numeric_or_NA)
                                       sapply(paramdf$dsc1, assert_numeric_or_NA)
                                       sapply(paramdf$dsc2, assert_numeric_or_NA)

                                       # demography
                                       assert_dataframe(demog)
                                       assert_in(colnames(demog), c("Strata", "popN"))
                                       assert_in(demog$Strata, IFRparams)
                                       assert_string(demog$strata)
                                       assert_pos_int(demog$popN)

                                       # censoring
                                       if (!is.null(rcensor_day)) {
                                         assert_pos_int(rcensor_day)
                                       }
                                     }

                                     # user dictionary key
                                     if (!is.null(IFRdictkey)) {
                                       assert_dataframe(IFRdictkey)
                                     }


                                     #......................
                                     # fill in
                                     #......................
                                     if (!is.null(data)) {
                                       if (unique(length(data$obs_deaths$ObsDay)) == 1) {
                                         stop("Time-Series data required but only one observation")
                                       }
                                       if (min(data$obs_deaths$ObsDay) != 1) {
                                         stop("Time-Series Data must start on day 1")
                                       }

                                       self$data <- data
                                       self$maxObsDay <- max(data$obs_deaths$ObsDay)
                                     }

                                     self$IFRparams <- IFRparams
                                     if (!is.null(maxMa)) {
                                       assert_in(maxMa, IFRparams)
                                       self$maxMA <- maxMa
                                     }

                                     # consider Ne if multiple strata
                                     if (length(IFRparams) > 1) {
                                       self$Noiseparams <- Noiseparams
                                     }

                                     self$Knotparams <- Knotparams
                                     if (!is.null(relKnot)) {
                                       assert_in(relKnot, Knotparams)
                                       self$relKnot <- relKnot
                                     }

                                     self$Infxnparams <- Infxnparams
                                     if (!is.null(relInfxn)) {
                                       assert_in(relInfxn, Infxnparams)
                                       self$relInfxn <- relInfxn
                                     }

                                     self$modparam <- modparam
                                     self$sodparam <- sodparam
                                     self$Serotestparams <- Serotestparams
                                     self$demog <- demog
                                     self$paramdf <- paramdf

                                     if (!is.null(rcensor_day)) {
                                       self$rcensor_day <- rcensor_day
                                     } else {
                                       self$rcensor_day <- .Machine$integer.max # user has elected to not censor (i.e. default no censoring)
                                     }

                                     self$IFRdictkey <- IFRdictkey
                                   },

                                   #......................
                                   # set functions
                                   #......................
                                   set_IFRparams = function(val) {
                                     assert_string(val)
                                     assert_unique(val)
                                     self$IFRparams <- val
                                   },

                                   set_maxMa = function(val) {
                                     if (length(self$IFRparams) == 0) {
                                       stop("Must specify IFR parameters before specifying the maximum mortality strata (maxMa)")
                                     }
                                     assert_string(val)
                                     assert_in(val, self$IFRparams)
                                     self$maxMa <- val
                                   },

                                   set_Noiseparams = function(val) {
                                     assert_string(val)
                                     assert_unique(val)
                                     if (is.null(self$IFRparams)) {
                                       stop("Must specificy IFR parmaeters before Noise Effect parameters")
                                     }
                                     if (length(self$IFRparams) <= 1) {
                                       stop("Cannot specify Noise Effect parameters when only one strata is being considered")
                                     }
                                     assert_same_length(val, self$IFRparams)
                                     self$Noiseparams <- val
                                   },


                                   set_Knotparams = function(val) {
                                     assert_string(val)
                                     assert_unique(val)
                                     self$Knotparams <- val
                                   },
                                   set_relKnot = function(val) {
                                     if (length(self$Knotparams) == 0) {
                                       stop("Must specify Knotparams before specifying the relative knot point (relKnot)")
                                     }
                                     assert_string(val)
                                     assert_in(val, self$Knotparams)
                                     self$relKnot <- val
                                   },

                                   set_Infxnparams = function(val) {
                                     assert_string(val)
                                     assert_unique(val)
                                     self$Infxnparams <- val
                                   },
                                   set_relInfxn = function(val) {
                                     if (length(self$Infxnparams) == 0) {
                                       stop("Must specify Infxnparams before specifying the relative infection point (relInfxn)")
                                     }
                                     assert_string(val)
                                     assert_in(val, self$Infxnparams)
                                     self$relInfxn <- val
                                   },

                                   set_MeanTODparam = function(val) {
                                     assert_string(val)
                                     self$modparam <- val
                                   },

                                   set_CoefVarOnsetTODparam = function(val) {
                                     assert_string(val)
                                     self$sodparam <- val
                                   },

                                   set_Serotestparams = function(val) {
                                     assert_string(val)
                                     assert_unique(val)
                                     assert_in(val, c("sens", "spec", "sero_con_rate", "sero_rev_shape", "sero_rev_scale"),
                                               message = "Serology test parameters currently limited to specifitiy (spec),
                                                          sensitivity (sens), lambda of the exponentially distributed onset of infection
                                                          to seroconversion (sero_con_rate), and the scale and shape of the Weibull
                                                          distribution for seroreversion after seroconversion (sero_rev_scale; sero_rev_shape).
                                                          Seroreversion parameters can be excluded if seroreversion should be ignored")
                                     self$Serotestparams <- val
                                   },

                                  # set data and paramdf
                                   set_data = function(val) {
                                     if (is.null(self$IFRparams)) {
                                       stop("Must specify IFR parameters before specifying data")
                                     }
                                     assert_list(val)
                                     assert_in(names(val), c("obs_deaths", "prop_deaths", "obs_serology"))
                                     # L1
                                     assert_dataframe(val$obs_deaths)
                                     assert_in(x = colnames(val$obs_deaths), y = c("ObsDay", "Deaths"))
                                     assert_numeric(val$obs_deaths$ObsDay)
                                     assert_increasing(val$obs_deaths$ObsDay)
                                     assert_numeric(val$obs_deaths$Deaths)
                                     # L2
                                     assert_dataframe(val$prop_deaths)
                                     assert_in(x = colnames(val$prop_deaths), y = c("Strata", "PropDeaths"))
                                     assert_bounded(val$prop_deaths$PropDeaths, left = 0, right = 1)
                                     assert_in(val$prop_deaths$Strata, self$IFRparams)
                                     assert_eq(sum(val$prop_deaths$PropDeaths), 1,
                                               message = c("Proportion of deaths by strata must sum to 1"))
                                     # L3
                                     assert_dataframe(val$obs_serology)
                                     assert_in(colnames(val$obs_serology), c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev", "SeroLCI", "SeroUCI"))
                                     assert_in(val$obs_serology$Strata, self$IFRparams)
                                     sapply(val$obs_serology$SeroPos, assert_pos_int_or_NA)
                                     sapply(val$obs_serology$SeroN, assert_pos_int_or_NA)
                                     sapply(val$obs_serology$SeroPrev, assert_bounded_or_NA, left = 0, right = 1)
                                     sapply(val$obs_serology$SeroLCI, assert_numeric_or_NA)
                                     sapply(val$obs_serology$SeroUCI, assert_numeric_or_NA)
                                     assert_pos_int(val$obs_serology$SeroStartSurvey)
                                     assert_pos_int(val$obs_serology$SeroEndSurvey)

                                     # extra checks
                                     if (unique(length(val$obs_deaths$ObsDay)) == 1) {
                                       stop("Time-Series data has only one observed day specified")
                                     }
                                     if (min(val$obs_deaths$ObsDay) != 1) {
                                       stop("Time-Series Data must start on day 1")
                                     }
                                     self$data <- val
                                     self$maxObsDay <- max(val$obs_deaths$ObsDay)
                                   },

                                   set_paramdf = function(val) {
                                     if (length(self$IFRparams) == 0 | length(self$Knotparams) == 0 | length(self$Infxnparams) == 0 | length(self$Serotestparams) == 0 | length(self$modparam) == 0 | length(self$sodparam) == 0) {
                                       stop("Must specify modparam, sodparam, IFRparams, Knotparams, Infxnparams, Serotestparams, and Noiseparams before specifying the param dataframe")
                                     }
                                     assert_dataframe(val)
                                     assert_in(x = colnames(val), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                     assert_string(val$name)
                                     assert_unique(val$name)
                                     assert_in(val$name, c(self$modparam, self$sodparam, self$IFRparams, self$Knotparams, self$Infxnparams, self$Serotestparams, self$Noiseparams))
                                     sapply(val$min, assert_numeric_or_NA)
                                     sapply(val$init, assert_numeric_or_NA)
                                     sapply(val$max, assert_numeric_or_NA)
                                     sapply(val$dsc1, assert_numeric_or_NA)
                                     sapply(val$dsc2, assert_numeric_or_NA)
                                     self$paramdf <- val
                                   },


                                   set_rcensor_day = function(val) {
                                     assert_pos_int(val)
                                     self$rcensor_day <- val
                                   },

                                   set_demog = function(val) {
                                     if (is.null(self$IFRparams)) {
                                       stop("Must specificy IFR parameters prior to specifying demography data")
                                     }
                                     assert_dataframe(val)
                                     assert_in(colnames(val), c("Strata", "popN"))
                                     assert_string(val$Strata)
                                     assert_in(val$Strata, self$IFRparams)
                                     assert_pos_int(val$popN)
                                     self$demog <- val
                                   },

                                   set_IFRdictkey = function(val) {
                                     assert_dataframe(val)
                                     self$IFRdictkey <- val
                                   }

                                 )
)







