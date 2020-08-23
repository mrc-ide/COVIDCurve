
set.seed(1234)
devtools::load_all()
library(drjacoby)
# sigmoidal function
infxns <- data.frame(time = 1:150)
sig <- function(x){1 / (1 +  exp(-x))}
timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                            min = -25,
                                            max = 50)

# make up fatality data
fatalitydata <- data.frame(Strata = c("ma1", "ma2", "ma3"),
                           IFR = c(0.05, 0.2, 0.5),
                           Rho = c(1, 1, 1),
                           Ne = c(1, 1, 1))
# make up demography data
demog <- data.frame(Strata = c("ma1", "ma2", "ma3"),
                    popN = c(1500000, 2250000, 1250000))
# pick serology date
sero_days <- c(110, 135)

#..................
# run sim
#..................
dat <- COVIDCurve::Aggsim_infxn_2_death(
  fatalitydata = fatalitydata,
  demog = demog,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 150,
  infections = infxns$infxns,
  simulate_seroprevalence = TRUE,
  sens = 0.85,
  spec = 0.95,
  sero_delay_rate = 10
)

#......................
# input data
#......................
# major liftover
obs_deaths <- dat$AggDeath %>%
  dplyr::group_by(ObsDay) %>%
  dplyr::summarise(deaths = sum(Deaths)) %>%
  dplyr::pull(deaths)



prop_age_obs_deaths <- dat$AggDeath %>%
  dplyr::group_by(Strata) %>%
  dplyr::summarise(deaths = sum(Deaths)) %>%
  dplyr::pull(deaths)/sum(dat$AggDeath$Deaths)

prop_rgn_obs_deaths <- dat$AggDeath %>%
  #dplyr::group_by(Strata) %>%
  dplyr::summarise(deaths = sum(Deaths)) %>%
  dplyr::pull(deaths)/sum(dat$AggDeath$Deaths)


obs_serology_age <- dat$AggSeroPrev %>%
  dplyr::filter(event_obs_day %in% sero_days) %>%
  dplyr::arrange(event_obs_day, Strata) %>%
  dplyr::mutate(
    n_pos = ObsPrev * popN,
    n_test = popN
  )

obs_serology_rgn <- dat$AggSeroPrev %>%
  dplyr::filter(event_obs_day %in% sero_days) %>%
  dplyr::mutate(
    n_pos = ObsPrev * popN,
    n_test = popN) %>%
  dplyr::group_by(event_obs_day) %>%
  dplyr::summarise(
    n_pos = sum(n_pos),
    n_test = sum(n_test)
  )

datinput <- list(obs_deaths = obs_deaths,
                 prop_age_obs_deaths = prop_age_obs_deaths,
                 prop_rgn_obs_deaths = prop_rgn_obs_deaths,
                 age_obs_serologypos = round(obs_serology_age$n_pos),
                 age_obs_serologyn = round(obs_serology_age$n_test),
                 rgn_obs_serologypos = round(obs_serology_rgn$n_pos),
                 rgn_obs_serologyn = round(obs_serology_rgn$n_test))

#..................
# inputs
#..................
# misc list
days_obsd <- 150
knots <- c(30, 60, 90, 120)
misc_list = list(rcensor_day = .Machine$integer.max,
                 days_obsd = days_obsd,
                 n_sero_obs = 2,
                 max_seroday_obsd = 140,
                 sero_survey_start = c(105, 130),
                 sero_survey_end = c(115, 140),
                 n_knots = length(knots)+1,
                 demog = demog$popN,
                 agestratlen = 3,
                 rgnstratlen = 1,
                 countmarginal_Rgn_demog = 5e6,
                 countmarginal_Age_demog = demog$popN,
                 total_observed_deaths = sum(dat$AggDeath$Deaths))

# liftover to Rcpp list
morelikely.paramsin <- c("mod" = 17.8, "sod" = 0.45,
                         "ma1" = 0.05, "ma2" = 0.2, "ma3" = 0.5,
                         "x1" = 30, "x2" = 60, "x3" = 90, "x4" = 120,
                         "y1" = 2.8, "y2" = 5.7, "y3" = 7.7, "y4" = 8.4, "y5" = 8.5,
                         "Ane1" = 0.33, "Ane2" = 0.33, "Ane3" = 0.33, "Rne1" = 1,
                         "sens" = 0.85, "spec" = 0.95, "sero_rate" = 10)


# truth
morelikely <- COVIDCurve:::natcubspline_loglike(params = morelikely.paramsin,
                                                param_i = 1,
                                                data = datinput,
                                                misc = misc_list)
morelikely



