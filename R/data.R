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



#' Model Run on Simulated
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
#'            # tidy up data for model input
#'               # here we are calculating the total proportion of deaths
#'               # within each age group
#'               # proportion deaths
#'               prop_deaths <- simdat$StrataAgg_TimeSeries_Death %>%
#'                 dplyr::group_by(Strata) %>%
#'                 dplyr::summarise(deaths = sum(Deaths)) %>%
#'                 dplyr::ungroup(.) %>%
#'                 dplyr::mutate(PropDeaths = deaths/sum(simdat$Agg_TimeSeries_Death$Deaths)) %>%
#'                 dplyr::select(-c("deaths"))
#'
#'               # Serologic studies are typically conducted over a few days
#'               # here we have assumed that two time points with 5 days on
#'               # either side of our "midpoint".
#'               # For each study and we will take the average over that time
#'
#'               sero_days <- c(135, 160)
#'               sero_days <- lapply(sero_days, function(x){seq(from = (x-5), to = (x+5), by = 1)})
#'               obs_serology <- simdat$StrataAgg_Seroprev %>%
#'                 dplyr::group_by(Strata) %>%
#'                 dplyr::filter(ObsDay %in% unlist(sero_days)) %>%
#'                 dplyr::mutate(serodaynum = sort(rep(1:length(sero_days), 11))) %>%
#'                 dplyr::mutate(
#'                   SeroPos = ObsPrev * testedN,
#'                   SeroN = testedN ) %>%
#'                 dplyr::group_by(Strata, serodaynum) %>%
#'                 dplyr::summarise(SeroPos = mean(SeroPos),
#'                                  SeroN = mean(SeroN)) %>% # seroN doesn't change
#'                 dplyr::mutate(SeroStartSurvey = sapply(sero_days, median) - 5,
#'                               SeroEndSurvey = sapply(sero_days, median) + 5,
#'                               SeroPos = round(SeroPos),
#'                               SeroPrev = SeroPos/SeroN,
#'                               SeroLCI = NA,
#'                               SeroUCI = NA) %>%
#'                 dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev", "SeroLCI", "SeroUCI")) %>%
#'                 dplyr::ungroup(.) %>%
#'                 dplyr::arrange(SeroStartSurvey, Strata)
#'
#'               ```
#'
#'               ## Building Model Object
#'               We will now make our model fitting object.
#'
#'               ```{r}
#'               #......................
#'               # put binomial data in correct format
#'               #......................
#'               inputdata <- list(obs_deaths = simdat$Agg_TimeSeries_Death,
#'                                 prop_deaths = prop_deaths,
#'                                 obs_serology = obs_serology)
#'
#'               #......................
#'               # Prior distrubitions and sampling distributions
#'               #......................
#'               # sens/spec
#'               sens_spec_tbl <- tibble::tibble(name =  c("sens",  "spec"),
#'                                               min =   c(0.5,      0.5),
#'                                               init =  c(0.85,     0.95),
#'                                               max =   c(1,        1),
#'                                               dsc1 =  c(850.5,    950.5),
#'                                               dsc2 =  c(150.5,    50.5))
#'
#'               # delay priors
#'               tod_paramsdf <- tibble::tibble(name = c("mod", "sod", "sero_con_rate"),
#'                                              min  = c(18,     0,     16),
#'                                              init = c(19,     0.85,  18),
#'                                              max =  c(20,     1,     21),
#'                                              dsc1 = c(19.8,   2550,  18.3),
#'                                              dsc2 = c(0.1,    450,   0.1))
#'
#'               # make param dfs
#'               ifr_paramsdf <- tibble::tibble(name = c("ma1", "ma2",  "ma3"),
#'                                              min  = rep(0, 3),
#'                                              init = rep(0.2, 3),
#'                                              max = rep(0.4, 3),
#'                                              dsc1 = rep(0, 3),
#'                                              dsc2 = rep(0.4, 3))
#'
#'               infxn_paramsdf <- tibble::tibble(name = paste0("y", 1:5),
#'                                                min  = rep(0, 5),
#'                                                init = rep(2, 5),
#'                                                max =  rep(10, 5),
#'                                                dsc1 = rep(0, 5),
#'                                                dsc2 = rep(10, 5))
#'
#'               knot_paramsdf <- tibble::tibble(name = paste0("x", 1:4),
#'                                               min  = c(2,  50, 100, 150),
#'                                               init = c(5,  75, 125, 175),
#'                                               max =  c(49, 99, 149, 200),
#'                                               dsc1 = c(2,  50, 100, 150),
#'                                               dsc2 = c(49, 99, 149, 200))
#'
#'
#'               noise_paramsdf <- tibble::tibble(name = c("ne1", "ne2", "ne3"),
#'                                                min  = rep(0.5, 3),
#'                                                init = rep(1, 3),
#'                                                max = rep(1.5, 3),
#'                                                dsc1 = rep(1, 3),
#'                                                dsc2 = rep(0.05, 3))
#'
#'               # bring together
#'               df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, noise_paramsdf, knot_paramsdf, sens_spec_tbl, tod_paramsdf)
#'
#'               #......................
#'               # make model for serorev and regular
#'               #......................
#'               # reg
#'               mod1 <- COVIDCurve::make_IFRmodel_age$new()
#'               mod1$set_MeanTODparam("mod")
#'               mod1$set_CoefVarOnsetTODparam("sod")
#'               mod1$set_IFRparams(paste0("ma", 1:3))
#'               mod1$set_maxMa("ma3")
#'               mod1$set_Knotparams(paste0("x", 1:4))
#'               mod1$set_Infxnparams(paste0("y", 1:5))
#'               mod1$set_Noiseparams(c(paste0("ne", 1:3)))
#'               mod1$set_Serotestparams(c("sens", "spec", "sero_con_rate"))
#'               mod1$set_data(inputdata)
#'               mod1$set_demog(demog)
#'               mod1$set_paramdf(df_params)
#'               mod1$set_rcensor_day(.Machine$integer.max) # no censoring
#' #'            modout <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod1,
#'                                                      reparamIFR = FALSE,
#'                                                      reparamInfxn = FALSE,
#'                                                      reparamKnot = FALSE,
#'                                                      burnin = 1e4,
#'                                                      samples = 1e4,
#'                                                      chains = 3,
#'                                                      rungs = 1,
#'                                                      thinning = 0,
#'                                                      silent = TRUE)
#'                       )
#'   }
#' }

"covidcurve_modfit"

