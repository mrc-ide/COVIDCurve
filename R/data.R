#' Simulated Data from Under the Model
#'
#' This data was simulated using the model simulator. It follows the
#' `runningmodel.Rmd` vignette.
#'
#' @format A list with three dataframes: (1) StrataAgg_TimeSeries_Death: stratified daily deaths by age group;
#'        (2) Agg_TimeSeries_Death: summed daily deaths (summed across age groups);
#'        (3) StrataAgg_Seroprev: daily seroprevalences by age group.
#' \describe{
#'   This file was generated with the following code:
#'   \dontrun{
#'            set.seed(1234)
#'            # users are required to input their own "infection" curve,
#'            # here we will use a simple sigmoidal function with some noise
#'            infxns <- data.frame(time = 1:200)
#'            sig <- function(x){1 / (1 +  exp(-x))}
#'            timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
#'            infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
#'                                             min = -25,
#'                                             max = 50)
#'
#'          # make up IFR values and attack rates within the three "ma" age groups
#'          fatalitydata <- tibble::tibble(Strata = c("ma1", "ma2", "ma3"),
#'                                IFR = c(0.01, 0.05, 0.1),
#'                                Rho = 1) # assuming uniform attack rate wrt demographic size
#'          # population sizes
#'          demog <- tibble::tibble(Strata = c("ma1", "ma2", "ma3"),
#'                         popN = c(1500000, 2250000, 1250000))
#'
#'         #..............................
#'         # Running the Simulation Function
#'         #...............................
#'         covidcurve_simdata <- COVIDCurve::Agesim_infxn_2_death(
#'                       fatalitydata = fatalitydata,
#'                       demog = demog,
#'                       m_od = 19.8,
#'                       s_od = 0.85,
#'                       curr_day = 200,
#'                       infections = infxns$infxns,
#'                       simulate_seroreversion = FALSE,
#'                       sens = 0.85,
#'                       spec = 0.95,
#'                       sero_delay_rate = 18.3
#'                       )
#'   }
#' }
#'
#'
"covidcurve_simdata"

