#########################################################################
# Purpose: Deploy
#
# Author: Nicholas F. Brazeau
#########################################################################
devtools::load_all()
library(drjacoby)
library(tidyverse)
# sigmoidal function
infxns <- data.frame(time = 1:200)
sig <- function(x){1 / (1 +  exp(-x))}
timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                            min = -25,
                                            max = 50)
sum(infxns$infxns < 0)

# make up fatality data
fatalitydata <- tibble::tibble(Strata = c("ma1", "ma2", "ma3"),
                               IFR = c(0.01, 0.05, 0.1),
                               Rho = 1)
demog <- tibble::tibble(Strata = c("ma1", "ma2", "ma3"),
                        popN = c(1500000, 2250000, 1250000))

#..............................................................
# AGGREGATE
#..............................................................
#..................
# run sim
#..................
dat <- COVIDCurve::Agesim_infxn_2_death(
  fatalitydata = fatalitydata,
  demog = demog,
  m_od = 19.66,
  s_od = 0.90,
  curr_day = 200,
  infections = infxns$infxns,
  simulate_seroreversion = FALSE,
  sens = 0.85,
  spec = 0.95,
  sero_delay_rate = 18.3,
  smplfrac = 1e-3
)

# liftover proprtion deaths
totdeaths <- sum(dat$StrataAgg_TimeSeries_Death$Deaths)
prop_strata_obs_deaths <- dat$StrataAgg_TimeSeries_Death %>%
  dplyr::group_by(Strata) %>%
  dplyr::summarise(Deaths = sum(Deaths),
                   PropDeaths = Deaths/totdeaths) %>%
  dplyr::select(c("Strata", "PropDeaths"))

# pick serology date
#sero_days <- c(150)
sero_days <- c(135, 180)

# sero_days full
sero_days <- lapply(sero_days, function(x){seq(from = (x-5), to = (x+5), by = 1)})
obs_serology <- dat$StrataAgg_Seroprev %>%
  dplyr::group_by(Strata) %>%
  dplyr::filter(ObsDay %in% unlist(sero_days)) %>%
  dplyr::mutate(serodaynum = sort(rep(1:length(sero_days), 11))) %>%
  dplyr::mutate(
    SeroPos = ObsPrev * testedN,
    SeroN = testedN ) %>%
  dplyr::group_by(Strata, serodaynum) %>%
  dplyr::summarise(SeroPos = mean(SeroPos),
                   SeroN = mean(SeroN)) %>% # seroN doesn't change
  dplyr::mutate(SeroStartSurvey = sapply(sero_days, median) - 5,
                SeroEndSurvey = sapply(sero_days, median) + 5,
                SeroPos = round(SeroPos),
                SeroPrev = SeroPos/SeroN) %>%
  dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev")) %>%
  dplyr::ungroup(.) %>%
  dplyr::arrange(SeroStartSurvey, Strata) %>%
  dplyr::mutate(SeroLCI = SeroPrev - 0.01,
                SeroUCI = SeroPrev + 0.01) # make up some tight CIs

datinput <- list(obs_deaths = dat$Agg_TimeSeries_Death,
                 prop_deaths = prop_strata_obs_deaths,
                 obs_serology = obs_serology)

#..................
# make model
#..................
ifr_paramsdf <- tibble::tibble(name = c("ma1", "ma2",  "ma3"),
                               min  = rep(0, 3),
                               init = rep(0.2, 3),
                               max = rep(0.4, 3),
                               dsc1 = rep(0, 3),
                               dsc2 = rep(0.4, 3))

infxn_paramsdf <- tibble::tibble(name = paste0("y", 1:5),
                                 min  = rep(0, 5),
                                 init = c(rep(0.5, 4), 5),
                                 max =  c(rep(1, 4), 10),
                                 dsc1 = rep(0, 5),
                                 dsc2 = c(rep(1, 4), 10))

# knot_paramsdf <- tibble::tibble(name = paste0("x", 1:4),
#                                 min  = c(0,    0.33, 0.66, 175),
#                                 init = c(0.05, 0.40, 0.75, 185),
#                                 max =  c(0.33, 0.66, 0.99, 200),
#                                 dsc1 = c(0,    0.33, 0.66, 175),
#                                 dsc2 = c(0.33, 0.66, 0.99, 200))

knot_paramsdf <- tibble::tibble(name = paste0("x", 1:4),
                                min  = c(rep(0,3), 175),
                                init = c(rep(0.2,3), 185),
                                max =  c(rep(1,3), 200),
                                dsc1 = c(rep(0,3), 175),
                                dsc2 = c(rep(1,3), 200))

sero_paramsdf <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.83,     0.8),
                                init =  c(0.85,     0.90),
                                max =   c(0.87,     1.00),
                                dsc1 =  c(8500,     950),
                                dsc2 =  c(1500,     50))


noise_paramsdf <- tibble::tibble(name = c("ne1", "ne2", "ne3"),
                                 min  = rep(0.5, 3),
                                 init = rep(1, 3),
                                 max = rep(1.5, 3),
                                 dsc1 = rep(1, 3),
                                 dsc2 = rep(0.05, 3))

# onset to deaths
tod_paramsdf <- tibble::tibble(name = c("mod", "sod", "sero_con_rate"),
                               min  = c(18,     0,     16),
                               init = c(19,     0.90,  18),
                               max =  c(20.5,   1,     21),
                               dsc1 = c(19.66,  2700,  18.3),
                               dsc2 = c(0.1,    300,   0.1))

df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sero_paramsdf, noise_paramsdf, tod_paramsdf)

#......................
# make mode
#......................
mod1 <- COVIDCurve::make_IFRmodel_age$new()
mod1$set_MeanTODparam("mod")
mod1$set_CoefVarOnsetTODparam("sod")
mod1$set_IFRparams(c("ma1", "ma2", "ma3"))
mod1$set_maxMa("ma3")
mod1$set_Knotparams(paste0("x", 1:4))
mod1$set_relKnot("x4")
mod1$set_Infxnparams(paste0("y", 1:5))
mod1$set_relInfxn("y5")
mod1$set_Serotestparams(c("sens", "spec", "sero_con_rate"))
mod1$set_Noiseparams(c("ne1", "ne2", "ne3"))
mod1$set_data(datinput)
mod1$set_demog(demog)
mod1$set_paramdf(df_params)
mod1$set_rcensor_day(.Machine$integer.max)

#..................
# run model
#..................
start <- Sys.time()
modout <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod1,
                                       reparamIFR = TRUE,
                                       reparamInfxn = TRUE,
                                       reparamKnot = TRUE,
                                       burnin = 1e4,
                                       samples = 1e4,
                                       chains = 3,
                                       rungs = 10,
                                       thinning = 10,
                                       silent = FALSE)
Sys.time() - start
modout


#............................................................
# examine outputs
#...........................................................
(ifr <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, whichrung = paste0("rung", 1),
                                       what = "IFRparams", by_chain = F))
(sero <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, whichrung = paste0("rung", 1),
                                        what = "Serotestparams", by_chain = F))
(knotspost <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                             what = "Knotparams", by_chain = F))
(infxn <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                         what = "Infxnparams", by_chain = F))

(ddelays <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                         what = "DeathDelayparams", by_chain = F))


summary(modout$mcmcout$output$loglikelihood)
summary(modout$mcmcout$output$logprior)
modout$mcmcout$output[modout$mcmcout$output$loglikelihood == max(modout$mcmcout$output$loglikelihood), ]


#............................................................
# curves and posteriors
#...........................................................
#......................
# ifrs and infxns
#......................
curve <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf = modout,
                                                        whichrung = paste0("rung", 1),
                                                        by_chain = F,
                                                        by_strata = F,
                                                        dwnsmpl = 1e3)
# tidy up and make plots
liftover <- data.frame(param = c("ma1", "ma2", "ma3"),
                       Strata = c("ma1", "ma2", "ma3"))

fatalitydataplot <- fatalitydata %>%
  dplyr::left_join(liftover, ., by = "Strata")

plot1 <- ggplot() +
  geom_pointrange(data = ifr, aes(x = param, ymin = LCI, ymax = UCI, y = median, color = param)) +
  geom_hline(data = fatalitydataplot, aes(yintercept  = IFR, group = param), color = "#3182bd", linetype = "dashed", size = 1.1) +
  facet_wrap(.~param) +
  scale_color_viridis_d() +
  theme_bw()

plot2 <- curve$plotObj +
  geom_line(data = infxns, aes(x = time, y = infxns), color = "#3182bd")

cowplot::plot_grid(plot1, plot2, ncol = 1, nrow = 2)


#......................
# serology
#......................
serocurve <- COVIDCurve::draw_posterior_sero_curves(IFRmodel_inf = modout,
                                                    whichrung = paste0("rung", 1),
                                                    by_chain = F,
                                                    dwnsmpl = 1e2)

serocurvedat <- serocurve %>%
  dplyr::select(c("sim", "ObsDay", dplyr::starts_with("RG_pd_"),
                  dplyr::starts_with("crude_pd_"))) %>%
  tidyr::pivot_longer(., cols = -c("sim", "ObsDay"),
                      names_to = "seroprev_strata_lvl", values_to = "seroprev") %>%
  dplyr::mutate(seroprevlvl = ifelse(stringr::str_detect(seroprev_strata_lvl, "RG_"), "RG Corr.", "Crude"),
                param = stringr::str_extract(seroprev_strata_lvl, "ma[0-9]+"))

SeroPrevObs <- tibble::tibble(obsdaymin = sapply(sero_days, median) - 5,
                              obsdaymax = sapply(sero_days, median) + 5)
serocurvedat %>%
  ggplot() +
  geom_line(aes(x = ObsDay, y = seroprev, color = seroprevlvl), alpha = 0.5) +
  geom_rect(data = SeroPrevObs, aes(xmin = obsdaymin, xmax = obsdaymax, ymin = -Inf, ymax = Inf),
            fill = "#d9d9d9", alpha = 0.4) +
  facet_wrap(.~param) +
  scale_color_manual("Seroprev. \n Adjustment", values = c("#FFD301", "#246BCF"),
                     labels = c("Inferred 'Truth'", "Inferred 'Observed' - \n Rogan-Gladen Corrected")) +
  labs(caption = "Grey box is observed seroprevalence across study period") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.90, hjust= 1, face = "bold"))

#......................
# get deaths posterior pred check
#......................
postdeaths <- COVIDCurve::posterior_check_infxns_to_death(IFRmodel_inf = modout,
                                                          dwnsmpl = 1e2,
                                                          by_chain = FALSE)
postdeaths.plotObj <- postdeaths %>%
  dplyr::select(c("time", dplyr::starts_with("deaths"))) %>%
  tidyr::gather(., key = "strata", value = "deaths", 2:ncol(.)) %>%
  ggplot() +
  geom_line(aes(x= time, y = deaths, group = strata, color = strata), size = 1.2) +
  scale_color_viridis_d()

postdeaths.plotObj +
  geom_line(data = dat$Agg_TimeSeries_Death, aes(x=ObsDay, y = Deaths, group = Strata), color = "#bdbdbd", size = 0.75) +
  theme_bw() +
  ggtitle("Posterior Predictive Check", subtitle = "Grey Lines are Simulated Data, Viridis Lines are Draws from Posterior")


