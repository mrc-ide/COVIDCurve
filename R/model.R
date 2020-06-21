#' @title Make Model Inference Object-Oriented Class
#' @description A simple R6 class with only public features. Purpose is to provide a model framework
#' for inference of IFRs and the incidence curve from Aggregate Death Data.
#' @name make_IFRmodel_agg
#' @section Public Variables:
#'  \tabular{ll}{
#'   \strong{Variable} \tab \strong{Return} \cr
#'   \code{data} \tab  Input data for inference \cr
#'   }
#' @section Constructor Arguments:
#' \tabular{lll}{
#' \strong{Argument} \tab \strong{Type} \tab \strong{Details} \cr
#' \code{data} \tab dataframe \tab Dataframe with .... \cr
#' }
#'
#' @export
NULL

#..............................................................
# declare
#..............................................................
make_IFRmodel_agg <- R6::R6Class(classname = "IFRmodel",
                                 public = list(
                                   data = NULL,
                                   maxObsDay = NULL,
                                   level = NULL,
                                   IFRparams = NULL,
                                   maxMa = NULL,
                                   Knotparams = NULL,
                                   relKnot = NULL,
                                   Infxnparams = NULL,
                                   relInfxn = NULL,
                                   paramdf = NULL,
                                   rho = NULL,
                                   mod = NULL,
                                   sod = NULL,
                                   gamma_lookup = NULL,
                                   # sero items
                                   Serodayparams = NULL,
                                   Serotestparams = NULL,
                                   popN = NULL,
                                   # censoring items
                                   rcensor_day = NULL,

                                   initialize = function(data = NULL, maxObsDay = NULL, level = NULL, rho = NULL,
                                                         mod = NULL, sod = NULL, gamma_lookup = NULL,
                                                         IFRparams = NULL, maxMa = NULL,
                                                         Infxnparams = NULL, relInfxn = NULL,
                                                         Knotparams = NULL, relKnot = NULL,
                                                         Serotestparams = NULL, Serodayparams = NULL, popN = NULL,
                                                         rcensor_day = NULL,
                                                         paramdf = NULL) {
                                     #......................
                                     # assertions and checks
                                     #......................
                                     items <- c(data, level, IFRparams, Infxnparams, Knotparams, Serotestparams, Serodayparams, paramdf, rho, mod, sod, popN)
                                     if ( !all(sapply(items, is.null)) ) { # if user tries to input things, assert otherwise initialize empty -- N.B., we initialize gamma_lookup later based on knots
                                       #assert level
                                       assert_in(x = level, y = c("Time-Series", "Cumulative"))
                                       # assert data
                                       assert_list(data)
                                       assert_in(names(data), c("obs_deaths", "obs_serologyrate"))
                                       assert_dataframe(data$obs_deaths)
                                       assert_in(x = colnames(data$obs_deaths), y = c("ObsDay", "Deaths"))
                                       assert_numeric(data$obs_deaths$ObsDay)
                                       assert_increasing(data$obs_deaths$ObsDay)
                                       assert_numeric(data$obs_deaths$Deaths)
                                       assert_numeric(data$obs_serologyrate)
                                       assert_bounded(data$obs_serologyrate, left = 0, right = 1)
                                       #assert params
                                       assert_string(IFRparams)
                                       assert_unique(IFRparams)
                                       assert_string(Infxnparams)
                                       assert_unique(Infxnparams)
                                       assert_string(Knotparams)
                                       assert_unique(Knotparams)
                                       assert_string(Serotestparams)
                                       assert_unique(Serotestparams)
                                       assert_in(Serotestparams, c("sens", "spec", "sero_rate"))
                                       assert_string(Serodayparams)
                                       assert_unique(Serodayparams)
                                       # assert paramdf
                                       assert_dataframe(paramdf)
                                       assert_in(x = colnames(paramdf), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                       assert_string(paramdf$name)
                                       assert_unique(paramdf$name)
                                       assert_in(paramdf$name, c(IFRparams, Infxnparams, Knotparams, Serotestparams, Serodayparams))
                                       assert_in(c(IFRparams, Infxnparams, Knotparams, Serotestparams, Serodayparams), paramdf$name)
                                       assert_numeric(paramdf$init)
                                       assert_numeric(paramdf$min)
                                       assert_numeric(paramdf$max)
                                       assert_numeric(paramdf$dsc1)
                                       assert_numeric(paramdf$dsc2)

                                       # OTD
                                       assert_numeric(mod)
                                       assert_numeric(sod)

                                       # censoring
                                       if (!is.null(rcensor_day)) {
                                         assert_pos_int(rcensor_day)
                                         assert_gr(rcensor_day, max(paramdf$max[paramdf$name %in% Serodayparams]),
                                                   message = "Day of censoring must be after maximum serology date")
                                       }

                                       # rho
                                       assert_numeric(rho)
                                       assert_same_length(rho, IFRparams)
                                       assert_eq(sum(rho), 1)
                                     }

                                     # fill in
                                     self$data <- data
                                     if (!is.null(data)) {
                                       if (level == "Time-Series") {
                                         if (unique(length(data$obs_deaths$ObsDay)) == 1) {
                                           stop("Time-Series specified but only one observation")
                                         }
                                         if (min(data$obs_deaths$ObsDay) != 1) {
                                           stop("Time-Series Data must start on day 1")
                                         }
                                       }

                                       day <- 1:(max(data$obs_deaths$ObsDay) + 1)
                                       self$gamma_lookup <- stats::pgamma((day-1), shape = 1/sod^2, scale = mod*sod^2)
                                       self$maxObsDay <- max(data$obs_deaths$ObsDay)
                                     }

                                     self$level <- level
                                     self$IFRparams <- IFRparams
                                     if (!is.null(maxMa)) {
                                       assert_in(maxMa, IFRparams)
                                       self$maxMA <- maxMa
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
                                     self$Serotestparams <- Serotestparams
                                     self$Serodayparams <- Serodayparams
                                     self$popN <- popN
                                     self$paramdf <- paramdf
                                     self$rho <- rho
                                     self$mod <- mod
                                     self$sod <- sod

                                     if (!is.null(rcensor_day)) {
                                       self$rcensor_day <- rcensor_day
                                     } else {
                                       self$rcensor_day <- .Machine$integer.max # user has elected to not censor (i.e. default no censoring)
                                     }
                                   },

                                   #......................
                                   # set functions
                                   #......................
                                   set_level = function(val) {
                                     assert_in(x = val, y = c("Time-Series", "Cumulative"))
                                     self$level <- val
                                   },

                                   set_data = function(val) {
                                     if (is.null(self$mod) | is.null(self$sod)) {
                                       stop("Must specify a Mean and Coeff. of Variation Delay of Onset to Death (mod & sod) before specifying data")
                                     }
                                     if (is.null(self$level)) {
                                       stop("Must specify a level before specifying data")
                                     }
                                     assert_list(val)
                                     assert_in(names(val), c("obs_deaths", "obs_serologyrate"))
                                     assert_dataframe(val$obs_deaths)
                                     assert_in(x = c("ObsDay", "Deaths"), y = colnames(val$obs_deaths))
                                     assert_numeric(val$obs_deaths$ObsDay)
                                     assert_increasing(val$obs_deaths$ObsDay)
                                     assert_numeric(val$obs_deaths$Deaths)
                                     assert_numeric(val$obs_serologyrate)
                                     assert_bounded(val$obs_serologyrate, left = 0, right = 1)

                                     if (self$level == "Time-Series") {
                                       if (unique(length(val$obs_deaths$ObsDay)) == 1) {
                                         stop("Time-Series specified but only one observation")
                                       }
                                       if (min(val$obs_deaths$ObsDay) != 1) {
                                         stop("Time-Series Data must start on day 1")
                                       }
                                     }
                                     self$data <- val
                                     day <- 1:(max(val$obs_deaths$ObsDay) + 1)
                                     self$gamma_lookup <- stats::pgamma((day-1), shape = 1/self$sod^2, scale = self$mod*self$sod^2)
                                     self$maxObsDay <- max(val$obs_deaths$ObsDay)
                                   },

                                   set_IFRparams = function(val) {
                                     assert_string(val)
                                     assert_unique(val)
                                     self$IFRparams <- val
                                   },

                                   set_maxMa = function(val) {
                                     if (length(self$IFRparams) == 0) {
                                       stop("Must specify IFRparams before specifying the maximum mortality strata (maxMa)")
                                     }
                                     assert_string(val)
                                     assert_in(val, self$IFRparams)
                                     self$maxMa <- val
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

                                   set_Serotestparams = function(val) {
                                     assert_string(val)
                                     assert_unique(val)
                                     assert_in(val, c("sens", "spec", "sero_rate"),
                                               message = "Serology parameters currently limited to specifitiy (spec),
                                             sensitivity (sens), and serology rate (sero_rate)")
                                     self$Serotestparams <- val
                                   },

                                   set_Serodayparams = function(val) {
                                     assert_string(val)
                                     assert_unique(val)
                                     self$Serodayparams <- val
                                   },

                                   set_paramdf = function(val) {
                                     if (length(self$IFRparams) == 0 | length(self$Knotparams) == 0 | length(self$Infxnparams) == 0 | length(self$Serotestparams) == 0 | length(self$Serodayparams) == 0) {
                                       stop("Must specify IFRparams, Knotparams, Infxnparams, Serotestparams, and Serodayparams before specifying the param dataframe")
                                     }
                                     assert_dataframe(val)
                                     assert_in(x = colnames(val), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                     assert_string(val$name)
                                     assert_unique(val$name)
                                     assert_in(val$name, c(self$IFRparams, self$Knotparams, self$Infxnparams, self$Serotestparams, self$Serodayparams))
                                     assert_numeric(val$init)
                                     assert_numeric(val$min)
                                     assert_numeric(val$max)
                                     assert_numeric(val$dsc1)
                                     assert_numeric(val$dsc2)
                                     self$paramdf <- val
                                   },

                                   set_MeanOnset = function(val) {
                                     assert_numeric(val)
                                     self$mod <- val
                                   },

                                   set_CoefVarOnset = function(val) {
                                     assert_numeric(val)
                                     self$sod <- val
                                   },

                                   set_rcensor_day = function(val) {
                                     if (is.null(self$paramdf)) {
                                       stop("Must specificy parameter dataframe prior to specifying day to right censor from")
                                     }
                                     if (is.null(self$Serodayparams)) {
                                       stop("Must specificy serology day parameters prior to specifying day to right censor from")
                                     }
                                     assert_pos_int(val)
                                     assert_gr(val, max(self$paramdf$max[self$paramdf$name %in% self$Serodayparams]),
                                               message = "Day of censoring must be after maximum serology date")
                                     self$rcensor_day <- val
                                   },

                                   set_rho = function(val) {
                                     assert_numeric(val)
                                     assert_same_length(val, self$IFRparams)
                                     assert_eq(sum(val), 1)
                                     self$rho <- val
                                   },

                                   set_popN = function(val) {
                                     assert_pos_int(val)
                                     self$popN <- val
                                   }
                                 )
)







