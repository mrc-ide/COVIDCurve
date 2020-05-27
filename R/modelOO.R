#' @title Object-Oriented Class for Model Inference with Aggregate
#'
#' @description A simple R6 class with only public features.
#' @name make_modinf_agg
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
make_modinf_agg <- R6::R6Class(classname = "Inference-Aggregate-Model",
                               public = list(
                                 data = NULL,
                                 level = NULL,
                                 IFRparams = NULL,
                                 maxMa = NULL,
                                 Infxnparams = NULL,
                                 paramdf = NULL,
                                 knots = NULL,
                                 pa = NULL,
                                 mod = NULL,
                                 sod = NULL,
                                 gamma_lookup = NULL,
                                 # sero items
                                 Seroparams = NULL,
                                 popN = NULL,

                                 initialize = function(data = NULL, level = NULL, IFRparams = NULL, maxMa = NULL, Infxnparams = NULL,
                                                       paramdf = NULL, knots = NULL, pa = NULL, mod = NULL, sod = NULL, gamma_lookup = NULL,
                                                       Seroparams = NULL, popN = NULL) {
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
                                     assert_string(Infxnparams)
                                     assert_string(Seroparams)
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

                                 # set functions
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
                                   self$Infxnparams <- val
                                 },

                                 set_Seroparams = function(val) {
                                   assert_string(val)
                                   assert_in(val, c("sens", "spec", "sero_rate", "sero_date"),
                                             message = "Serology parameters currently limited to specifitiy (spec),
                                             sensitivity (sens), serology rate (sero_rate) and date of serology rate (sero_date)")
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









#............................................................
# line list
#...........................................................
#' @title Object-Oriented Class for Model Inference with LineList
#'
#' @description A simple R6 class with only public features.
#' @name make_modinf_linelist
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
make_modinf_linelist <- R6::R6Class(classname = "Inference-LineList-Model",
                                    public = list(
                                      data = NULL,
                                      DeathModparam = NULL,
                                      DeathSodparam = NULL,
                                      RecovMorparam = NULL,
                                      RecovSorparam = NULL,
                                      IFRparams = NULL,
                                      paramdf = NULL,
                                      # init
                                      initialize = function(data = NULL, DeathModparam = NULL, DeathSodparam = NULL,
                                                            RecovMorparam = NULL, RecovSorparam = NULL, IFRparams = NULL, paramdf = NULL){

                                        items <- c(data, DeathModparam, DeathSodparam, RecovMorparam, RecovSorparam, IFRparams, paramdf)
                                        if ( !all(sapply(items, is.null)) ) { # if use tries to input things, assert otherwise initialize empty

                                          # assert data
                                          assert_dataframe(data)
                                          assert_in(x = colnames(data), y = c("AgeGroup", "OnsetDay", "Outcome", "EventDay"))
                                          sapply(data$AgeGroup, assert_factor)
                                          sapply(data$OnsetDay, assert_numeric)
                                          assert_in(data$outcome, c("Death", "Recovery"))
                                          sapply(data$EventDay, assert_numeric)
                                          # asssert gamma shape params and ifr params
                                          assert_string(DeathModparam)
                                          assert_string(DeathSodparam)
                                          assert_string(RecovMorparam)
                                          assert_string(RecovSorparam)
                                          sapply(IFRparams, assert_string)
                                          # assert paramdf
                                          assert_dataframe(paramdf)
                                          assert_in(x = colnames(paramdf), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                          sapply(paramdf$name, assert_string)
                                          assert_in(paramdf$name, c(self$IFRparams, self$DeathModparam, self$DeathSodparam, self$RecovMorparam, self$RecovSorparam))
                                          sapply(paramdf$init, assert_numeric)
                                          sapply(paramdf$min, assert_numeric)
                                          sapply(paramdf$max, assert_numeric)
                                          sapply(paramdf$dsc1, assert_numeric)
                                          sapply(paramdf$dsc2, assert_numeric)
                                          # assert paramdf and data match/IFR params
                                          assert_in(unique(data$AgeGroup), self$IFRparams,
                                                    message = "Grouping variable in the data (i.e. data$AgeGroup)
                                                            must be the same names as the IFR parameter names")

                                        }
                                        # fill in
                                        self$data <- data
                                        self$DeathModparam <- DeathModparam
                                        self$DeathSodparam <- DeathSodparam
                                        self$RecovMorparam <- RecovMorparam
                                        self$RecovSorparam <- RecovSorparam
                                        self$IFRparams <- IFRparams
                                        self$paramdf <- paramdf


                                      },

                                      # set functions
                                      set_data = function(val) {
                                        assert_dataframe(val)
                                        assert_in(x = colnames(val), y = c("AgeGroup", "OnsetDay", "Outcome", "EventDay"))
                                        sapply(val$AgeGroup, assert_factor)
                                        sapply(val$OnsetDay, assert_numeric)
                                        assert_in(val$Outcome, c("Death", "Recovery"))
                                        sapply(val$EventDay, assert_numeric)

                                        self$data <- val
                                      },

                                      set_DeathModparam = function(val) {
                                        assert_string(val)
                                        self$DeathModparam <- val
                                      },

                                      set_DeathSodparam = function(val) {
                                        assert_string(val)
                                        self$DeathSodparam <- val
                                      },

                                      set_RecovMorparam = function(val) {
                                        assert_string(val)
                                        self$RecovMorparam <- val
                                      },

                                      set_RecovSorparam = function(val) {
                                        assert_string(val)
                                        self$RecovSorparam <- val
                                      },

                                      set_IFRparams = function(val) {
                                        sapply(val, assert_string)
                                        # assert paramdf and data match/IFR params
                                        assert_non_null(self$data, message = "Must set data before setting IFRparams")
                                        assert_in(unique(self$data$AgeGroup), val,
                                                  message = "Grouping variable in the data (i.e. data$AgeGroup)
                                                            must be the same names as the IFR parameter names")
                                        self$IFRparams <- val
                                      },

                                      set_paramdf = function(val) {
                                        if (length(self$IFRparams) == 0 |
                                            length(self$set_DeathModparam) == 0 |
                                            length(self$set_DeathSodparam) == 0 |
                                            length(self$set_RecovMorparam) == 0 |
                                            length(self$set_RecovSorparam) == 0) {
                                          stop("Must specify IFRparams and Gamma Shape Paramerts for Death and Recovery before specifying the param dataframe")
                                        }
                                        assert_dataframe(val)
                                        assert_in(x = colnames(val), y = c("name", "init", "min", "max", "dsc1", "dsc2"))
                                        sapply(val$name, assert_string)
                                        assert_in(val$name, c(self$IFRparams, self$DeathModparam, self$DeathSodparam, self$RecovMorparam, self$RecovSorparam))
                                        sapply(val$init, assert_numeric)
                                        sapply(val$min, assert_numeric)
                                        sapply(val$max, assert_numeric)
                                        sapply(val$dsc1, assert_numeric)
                                        sapply(val$dsc2, assert_numeric)
                                        self$paramdf <- val
                                      }
                                    )
)

