#' @title
#' @description

#' @title Object-Oriented Class for Model Inference
#'
#' @description A simple R6 class with only public features.
#' @name Distribution
#' @section Constructor Arguments:
#' \tabular{lll}{
#' \strong{Argument} \tab \strong{Type} \tab \strong{Details} \cr
#' \code{data} \tab dataframe \tab Dataframe with .... \cr
#' }
#' @section Public Variables:
#'  \tabular{ll}{
#'   \strong{Variable} \tab \strong{Return} \cr
#'   \code{data} \tab  Input data for ifnerence \cr
#'   }
#'
#' @export
NULL

#..............................................................
# declare
#..............................................................
make_modinf <- R6::R6Class(classname = "Inference-Model",
                           public = list(
                             data = NULL,
                             level = NULL,
                             IFRparams = NULL,
                             Infxnparams = NULL,
                             paramdf = NULL,
                             knots = NULL,
                             pa = NULL,
                             mod = NULL,
                             sod = NULL,
                             gamma_lookup = NULL,

                             initialize = function(data = NULL, level = NULL, IFRparams = NULL, Infxnparams = NULL,
                                                   paramdf = NULL, knots = NULL, pa = NULL, mod = NULL, sod = NULL, gamma_lookup = NULL) {
                               items <- c(data, level, IFRparams, Infxnparams, paramdf, knots, pa, mod, sod)
                               if ( !all(sapply(items, is.null)) ) { # if use tries to input things, assert otherwise initialize empty
                                 #assert level
                                 assert_in(x = level, y = c("Time-Series", "Cumulative"))
                                 # assert data
                                 assert_dataframe(data)
                                 assert_in(x = colnames(data), y = c("ObsDay", "AgeGroup", "Deaths"))
                                 sapply(data$ObsDay, assert_numeric)
                                 assert_increasing(data$ObsDay)
                                 sapply(data$AgeGroup, assert_factor)
                                 sapply(data$Deaths, assert_numeric)
                                 #assert params
                                 sapply(IFRparams, assert_string)
                                 sapply(Infxnparams, assert_string)
                                 # assert paramdf
                                 assert_dataframe(paramdf)
                                 assert_in(x = colnames(paramdf), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                 sapply(paramdf$name, assert_string)
                                 assert_in(paramdf$name, c(self$IFRparams, self$Infxnparams))
                                 assert_in(paramdf$name, c(self$IFRparams, self$Infxnparams))

                                 assert_in(self$IFRparams, paramdf$name)
                                 assert_in(self$Infxnparams, paramdf$name)

                                 sapply(paramdf$init, assert_numeric)
                                 sapply(paramdf$min, assert_numeric)
                                 sapply(paramdf$max, assert_numeric)
                                 sapply(paramdf$dsc1, assert_numeric)
                                 sapply(paramdf$dsc2, assert_numeric)
                                 # OTD
                                 assert_numeric(mod)
                                 assert_numeric(sod)
                                 # knots
                                 sapply(knots, assert_numeric)
                                 assert_same_length(knots, Infxnparams)
                                 # pa
                                 sapply(pa, assert_numeric)
                                 assert_same_length(pa, IFRparams)
                                 assert_eq(sum(pa), 1)
                               }


                               # fill in
                               self$data <- data
                               self$level <- level
                               self$IFRparams <- IFRparams
                               self$Infxnparams <- Infxnparams
                               self$paramdf <- paramdf
                               self$knots <- knots
                               self$pa <- pa
                               self$mod <- mod
                               self$sod <- sod
                               if (!is.null(self$knots)) {
                                 day <- self$knots[1]:self$knots[length(self$knots)]
                                 self$gamma_lookup <- stats::pgamma((day-1), shape = 1/self$sod^2, scale = self$mod*self$sod^2)
                               } else {
                                 self$gamma_lookup <- gamma_lookup
                               }

                             },

                             set_level = function(val) {
                               assert_in(x = val, y = c("Time-Series", "Cumulative"))
                               self$level <- val
                             },

                             set_data = function(val) {
                               if (is.null(self$level)) {
                                 stop("Must specify a level before specifying data")
                               }
                               assert_dataframe(val)
                               assert_in(x = colnames(val), y = c("ObsDay", "AgeGroup", "Deaths"))
                               sapply(val$ObsDay, assert_numeric)
                               assert_increasing(val$ObsDay)
                               sapply(val$AgeGroup, assert_factor)
                               sapply(val$Deaths, assert_numeric)

                               if (self$level == "Time-Series") {
                                 if (unique(length(val$ObsDay)) == 1) {
                                   stop("Time-Series specified but only one observation")
                                 }
                               }
                               self$data <- val
                             },

                             set_IFRparams = function(val) {
                               sapply(val, assert_string)
                               self$IFRparams <- val
                             },

                             set_Infxnparams = function(val) {
                               sapply(val, assert_string)
                               self$Infxnparams <- val
                             },

                             set_paramdf = function(val) {
                               if (length(self$IFRparams) == 0 | length(self$Infxnparams) ==0) {
                                 stop("Must specify IFRparams and Infxnparams before specifying the param dataframe")
                               }
                               assert_dataframe(val)
                               assert_in(x = colnames(val), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                               sapply(val$name, assert_string)
                               assert_in(val$name, c(self$IFRparams, self$Infxnparams))
                               sapply(val$init, assert_numeric)
                               sapply(val$min, assert_numeric)
                               sapply(val$max, assert_numeric)
                               sapply(val$dsc1, assert_numeric)
                               sapply(val$dsc2, assert_numeric)
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

                             set_knots = function(val) {
                               if (is.null(self$mod) | is.null(self$sod)) {
                                 stop("Must specificy the mean and coefficient of variation for the onset-to-death distribution prior to specifying knots")
                               }
                               sapply(val, assert_numeric)
                               assert_same_length(val, self$Infxnparams)
                               self$knots <- val

                               # get gamma look up table
                               day <- self$knots[1]:self$knots[length(self$knots)]
                               self$gamma_lookup <- stats::pgamma((day-1), shape = 1/self$sod^2, scale = self$mod*self$sod^2)
                             },

                             set_pa = function(val) {
                               sapply(val, assert_numeric)
                               assert_same_length(val, self$IFRparams)
                               assert_eq(sum(val), 1)
                               self$pa <- val
                             }
                           )
)
