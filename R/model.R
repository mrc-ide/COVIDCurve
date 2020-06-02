#' @title Make Model Inference Object-Oriented Class
#' @description A simple R6 class with only public features. Purpose is to provide a model framework
#' for inference of IFRs and the incidence curve from Aggregate Death Data.
#' @name make_modinf_agg
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
make_modinf_agg <- R6::R6Class(classname = "IFRmodel",
                               public = list(
                                 data = NULL,
                                 level = NULL,
                                 IFRparams = NULL,
                                 maxMa = NULL,
                                 Infxnparams = NULL,
                                 relInfxn = NULL,
                                 paramdf = NULL,
                                 knots = NULL,
                                 pa = NULL,
                                 mod = NULL,
                                 sod = NULL,
                                 gamma_lookup = NULL,
                                 # sero items
                                 Seroparams = NULL,
                                 popN = NULL,

                                 initialize = function(data = NULL, level = NULL, knots = NULL, pa = NULL,
                                                       mod = NULL, sod = NULL, gamma_lookup = NULL,
                                                       IFRparams = NULL, maxMa = NULL,
                                                       Infxnparams = NULL, relInfxn = NULL,
                                                       Seroparams = NULL, popN = NULL,
                                                       paramdf = NULL) {
                                   #......................
                                   # assertions and checks
                                   #......................
                                   items <- c(data, level, IFRparams, Infxnparams, paramdf, knots, pa, mod, sod,
                                              Seroparams, popN)
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
                                     assert_string(Seroparams)
                                     assert_unique(Seroparams)
                                     assert_in(Seroparams, c("sens", "spec", "sero_day", "sero_rate"))
                                     # assert paramdf
                                     assert_dataframe(paramdf)
                                     assert_in(x = colnames(paramdf), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                     assert_string(paramdf$name)
                                     assert_in(paramdf$name, c(IFRparams, Infxnparams,  Seroparams))
                                     assert_in(c(IFRparams, Infxnparams, Seroparams), paramdf$name)

                                     assert_numeric(paramdf$init)
                                     assert_numeric(paramdf$min)
                                     assert_numeric(paramdf$max)
                                     assert_numeric(paramdf$dsc1)
                                     assert_numeric(paramdf$dsc2)
                                     # OTD
                                     assert_numeric(mod)
                                     assert_numeric(sod)
                                     # knots
                                     assert_numeric(knots)
                                     assert_same_length(knots, Infxnparams)
                                     # pa
                                     assert_numeric(pa)
                                     assert_same_length(pa, IFRparams)
                                     assert_eq(sum(pa), 1)
                                   }

                                   # fill in
                                   self$data <- data
                                   self$level <- level
                                   self$IFRparams <- IFRparams
                                   if (!is.null(maxMa)) {
                                     self$maxMA <- maxMa
                                   }
                                   self$Infxnparams <- Infxnparams
                                   if (!is.null(relInfxn)) {
                                     self$relInfxn <- relInfxn
                                   }
                                   self$Seroparams <- Seroparams
                                   self$popN <- popN
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

                                 #......................
                                 # set functions
                                 #......................
                                 set_level = function(val) {
                                   assert_in(x = val, y = c("Time-Series", "Cumulative"))
                                   self$level <- val
                                 },

                                 set_data = function(val) {
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
                                     if (unique(length(val$ObsDay)) == 1) {
                                       stop("Time-Series specified but only one observation")
                                     }
                                   }
                                   self$data <- val
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

                                 set_Seroparams = function(val) {
                                   assert_string(val)
                                   assert_unique(val)
                                   assert_in(val, c("sens", "spec", "sero_rate", "sero_day"),
                                             message = "Serology parameters currently limited to specifitiy (spec),
                                             sensitivity (sens), serology rate (sero_rate) and date of serology rate (sero_day)")
                                   self$Seroparams <- val
                                 },

                                 set_paramdf = function(val) {
                                   if (length(self$IFRparams) == 0 | length(self$Infxnparams) == 0 | length(self$Seroparams) ==0) {
                                     stop("Must specify IFRparams, Infxnparams, and Seroparams before specifying the param dataframe")
                                   }
                                   assert_dataframe(val)
                                   assert_in(x = colnames(val), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                   assert_string(val$name)
                                   assert_in(val$name, c(self$IFRparams, self$Infxnparams, self$Seroparams))
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

                                 set_knots = function(val) {
                                   if (is.null(self$mod) | is.null(self$sod)) {
                                     stop("Must specificy the mean and coefficient of variation for the onset-to-death distribution prior to specifying knots")
                                   }
                                   assert_numeric(val)
                                   assert_same_length(val, self$Infxnparams)
                                   self$knots <- val

                                   # get gamma look up table
                                   day <- self$knots[1]:self$knots[length(self$knots)]
                                   self$gamma_lookup <- stats::pgamma((day-1), shape = 1/self$sod^2, scale = self$mod*self$sod^2)
                                 },

                                 set_pa = function(val) {
                                   assert_numeric(val)
                                   assert_same_length(val, self$IFRparams)
                                   assert_eq(sum(val), 1)
                                   self$pa <- val
                                 },

                                 set_popN = function(val) {
                                   assert_pos_int(val)
                                   self$popN <- val
                                 }
                               )
)







