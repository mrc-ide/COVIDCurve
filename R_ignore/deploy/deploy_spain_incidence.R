#########################################################################
# Purpose: Deploy
#
# Author: Nicholas F. Brazeau
#########################################################################
set.seed(48)
library(drjacoby)
devtools::load_all()
library(squire)
# Simulate Data
country <- "Spain"
start_date <- as.Date("2020-02-01")
base_R0 <- 3
daysobs <- 136
R0_change <- c(seq(base_R0, 1, -0.5), 0.9, 1.4, 1.5)/base_R0
R0_change <- tail(R0_change, -1)
tt_R0 <- c(0, 40, 50, 65, 70, 80, 130, 132)
date_R0_change <- start_date + tail(tt_R0, -1)
get <- squire::run_explicit_SEEIR_model(country,
                                        R0 = c(base_R0, base_R0 * R0_change),
                                        tt_R0 = tt_R0,
                                        replicates = 1,
                                        day_return = TRUE,
                                        time_period = daysobs,
                                        dt = 0.25)

infxns <- format_output(get, var_select = "E")

infxns <- data.frame(time = 1:(daysobs+1),
                     infxns = infxns$y)
plot(infxns$infxns)
cumsum(infxns$infxns)[131]/sum(get_population(country)$n)
# divide by 10 but same shape
infxns$infxns <- infxns$infxns/10

# make up fatality data
fatalitydata <- data.frame(strata = paste0("ma", 1:10),
                           ifr = c(1e-26, 1e-26, 1e-26, 1e-26, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2),
                           pa = 1/10)


#..............................................................
# AGGREGATE
#..............................................................
#..................
# run sim
#..................
dat <- COVIDCurve::Aggsim_infxn_2_death(
  fatalitydata = fatalitydata,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = (daysobs+1),
  level = "Time-Series",
  infections = infxns$infxns,
  simulate_seroprevalence = TRUE,
  sens = 0.85,
  spec = 0.99,
  sero_delay_rate = 10,
  popN = sum(get_population(country)$n)
)

save(dat, get, infxns, file = "sandbox/squire_spain.rda")
load("sandbox/squire_spain.rda")
# pick serology date
sero_day <- 125
datin <- list(obs_deaths = dat$AggDat,
              obs_serologyrate = dat$seroprev$ObsPrev[sero_day])

#..................
# make model
#..................
ifr_paramsdf <- tibble::tibble(name = paste0("ma", 1:10),
                               min  = rep(0, 10),
                               init = rep(0.5, 10),
                               max = rep(1, 10),
                               dsc1 = rep(0, 10),
                               dsc2 = rep(1, 10))

infxndf_params <- tibble::tibble(name = paste0("y", 1:5),
                                 min  = c(0,   0,   0,  0,   0),
                                 init = c(0.5, 0.5, 9,  0.5, 0.5),
                                 max =  c(1,   1,   14, 1,   1),
                                 dsc1 = c(0,   0,    0,  0,   0),
                                 dsc2 = c(1,   1,    14, 1,   1))
# knotdf_params <- tibble::tibble(name = paste0("x", 1:4),
#                                 min  = c(0,    0.30, 0.50, 125),
#                                 init = c(0.05, 0.50, 0.60, 128),
#                                 max =  c(0.5,  0.75, 0.99, 137),
#                                 dsc1 = c(0,    0.30, 0.50, 125),
#                                 dsc2 = c(0.5,  0.75, 0.99, 137))
knotdf_params <- tibble::tibble(name = paste0("x", 1:4),
                                min  = c(40, 70, 90, 135),
                                init = c(40, 70, 90, 135),
                                max =  c(42, 72, 92, 137),
                                dsc1 = c(40, 70, 90, 135),
                                dsc2 = c(42, 72, 92, 137))

sero_paramsdf <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day"),
                                min =   c(0.83,     0.97,   10,         120),
                                init =  c(0.85,     0.99,   10,         125),
                                max =   c(0.87,     1.00,   10,         130),
                                dsc1 =  c(8500,     9900,    5,         125),
                                dsc2 =  c(1500,     100,     15,        0.1))

df_params <- rbind.data.frame(ifr_paramsdf, infxndf_params, knotdf_params, sero_paramsdf)

#......................
# make mode
#......................
mod1 <- make_IFRmodel_agg$new()
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_level("Time-Series")
mod1$set_data(datin)
mod1$set_IFRparams(paste0("ma", 1:10))
mod1$set_maxMa("ma10")
mod1$set_Knotparams(paste0("x", 1:4))
mod1$set_relKnot("x4")
mod1$set_Infxnparams(paste0("y", 1:5))
mod1$set_relInfxn("y3")
mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_day"))
mod1$set_popN(sum(squire::get_population("Spain")$n))
mod1$set_paramdf(df_params)
mod1$set_pa(rep(1/10, 10))
mod1$set_rcensor_day(.Machine$integer.max)

#..................
# run model
#..................
modout <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod1,
                                       reparamIFR = TRUE,
                                       reparamInfxn = TRUE,
                                       reparamKnot = FALSE,
                                       burnin = 1e3,
                                       samples = 1e3,
                                       rungs = 10,
                                       GTI_pow = 4.0,
                                       chains = 3)
summary(modout$mcmcout$output$loglikelihood[modout$mcmcout$output$rung == "rung1" & modout$mcmcout$output$stage == "sampling"])
plot(modout$mcmcout$output$loglikelihood[modout$mcmcout$output$rung == "rung1" & modout$mcmcout$output$stage == "sampling"])
summary(modout$mcmcout$output$logprior[modout$mcmcout$output$rung == "rung1" & modout$mcmcout$output$stage == "sampling"])
plot(modout$mcmcout$output$logprior[modout$mcmcout$output$rung == "rung1" & modout$mcmcout$output$stage == "sampling"])

plot_par(modout$mcmcout, "ma9", chain = 2)
plot_par(modout$mcmcout, "ma10")
plot_par(modout$mcmcout, "sens")
plot_par(modout$mcmcout, "spec")
plot_par(modout$mcmcout, "sero_day")
plot_par(modout$mcmcout, "x1")
plot_par(modout$mcmcout, "x2")
plot_par(modout$mcmcout, "x3")
plot_par(modout$mcmcout, "x4")
plot_par(modout$mcmcout, "y1")
plot_par(modout$mcmcout, "y2")
plot_par(modout$mcmcout, "y3")
plot_par(modout$mcmcout, "y4")
plot_par(modout$mcmcout, "y5")


plot_mc_acceptance(modout$mcmcout)
drjacoby::plot_rung_loglike(modout$mcmcout)
drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 1, y_axis_type = 2)
drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 1, y_axis_type = 3)
drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 2, y_axis_type = 2)
drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 2, y_axis_type = 3)
(ifr <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, whichrung = paste0("rung", 1),
                                       what = "IFRparams", by_chain = T))
(knotspost <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                             what = "Knotparams", by_chain = F))
(infxn <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                         what = "Infxnparams", by_chain = F))
(sero <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                        what = "Seroparams", by_chain = F))


curve <- COVIDCurve::draw_posterior_infxn_points_cubic_splines(IFRmodel_inf = modout,
                                                               whichrung = paste0("rung", 1),
                                                               by_chain = TRUE,
                                                               CIquant = 0.95)
# plot out
jpeg("~/Desktop/posterior_curve_draws.jpg", width = 11, height = 8, units = "in", res = 500)
library(ggplot2)

fatalitydataplot <- fatalitydata %>%
  dplyr::rename(param = strata) %>%
  dplyr::mutate(param = factor(param, levels = paste0("ma", 1:10)))
ifr <- ifr %>%
  dplyr::mutate(param = factor(param, levels = paste0("ma", 1:10)))

plot1 <- ggplot() +
  geom_hline(data = fatalitydataplot, aes(yintercept  = ifr, group = param), color = "#3182bd", linetype = "dashed", size = 1.1) +
  geom_pointrange(data = ifr, aes(x = param, ymin = LCI, ymax = UCI, y = median, group = param)) +
  facet_wrap(.~param, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank())


plot2 <- curve$plotObj +
  geom_line(data = infxns, aes(x = time, y = infxns), color = "#3182bd")

cowplot::plot_grid(plot1, plot2, ncol = 1, nrow = 2)

graphics.off()

# save out
params <- rbind.data.frame(ifr, infxn, sero)
saveRDS(params,
        file = "sandbox/infparams.RDS")


#......................
# get deaths posterior pred check
#......................
postdeaths <- COVIDCurve::posterior_check_infxns_to_death(IFRmodel_inf = modout,
                                                          CIquant = 0.95,
                                                          by_chain = FALSE)
postdeaths.plotObj <- postdeaths %>%
  dplyr::select(c("time", dplyr::starts_with("deaths"))) %>%
  tidyr::gather(., key = "strata", value = "deaths", 2:ncol(.)) %>%
  ggplot() +
  geom_line(aes(x= time, y = deaths, group = strata, color = strata), size = 1.2) +
  scale_color_viridis_d()

jpeg("~/Desktop/posterior_check.jpg", width = 11, height = 8, units = "in", res = 500)
postdeaths.plotObj +
  geom_line(data = dat$obs_deaths, aes(x=ObsDay, y = Deaths, group = Strata), color = "#bdbdbd", size = 0.75) +
  theme_bw() +
  ggtitle("Posterior Predictive Check", subtitle = "Grey Lines are Simulated Data, Viridis Lines are Draws from Posterior")
graphics.off()




plot_par(modout$mcmcout, "r1", rung = 1)
plot_par(modout$mcmcout, "r2", rung = 1)
plot_par(modout$mcmcout, "ma3", rung = 1)
plot_par(modout$mcmcout, "y1", rung = 1)
plot_par(modout$mcmcout, "y2", rung = 1)
plot_par(modout$mcmcout, "y3", rung = 1)
plot_par(modout$mcmcout, "y4", rung = 1)
plot_par(modout$mcmcout, "y5", rung = 1)
plot_par(modout$mcmcout, "y6", rung = 1)
plot_par(modout$mcmcout, "sens", rung = 1)
plot_par(modout$mcmcout, "spec", rung = 1)
plot_par(modout$mcmcout, "sero_day", rung = 1)

plot_mc_acceptance(modout$mcmcout)
plot_rung_loglike(modout$mcmcout)
plot_rung_loglike(modout$mcmcout, y_axis_type = 2)
plot_rung_loglike(modout$mcmcout, y_axis_type = 3)


plot_cor(modout$mcmcout, "r1", "y4")
plot_cor(modout$mcmcout, "ma3", "y4")



#............................................................
# look at rungs
#...........................................................
COVIDCurve::get_cred_intervals(IFRmodel = mod1,
                               mcmcout = modout,
                               whichrung = "rung50",
                               what = "IFRparams",
                               by_chain = F)

modout$output %>%
  dplyr::mutate(
    logpost = loglikelihood + logprior
  ) %>%
  dplyr::filter(stage == "sampling" &
                  rung == "rung50") %>%
  ggplot() +
  geom_histogram(aes(x  = logpost))

rungdf <- modout$output %>%
  dplyr::mutate(
    logpost = loglikelihood + logprior
  ) %>%
  dplyr::filter(stage == "sampling" &
                  rung == "rung50")

minloglike <- min(rungdf$logpost)
# is this neg infinity
if (minloglike < -1.796e306 & minloglike > -1.8e306) {
  sum(rungdf$logpost == minloglike)/nrow(rungdf)
}

#............................................................
# are the problem ones the rungs that are inf?
#...........................................................
rungdf %>%
  dplyr::filter(logpost > -1.8e306) %>%
  dplyr::select_at(c("chain", "r1", "r2", "ma3")) %>%
  tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
  dplyr::group_by_at("param") %>%
  dplyr::summarise(
    min = min(est),
    LCI = quantile(est, 0.025),
    median = median(est),
    mean = mean(est),
    UCI = quantile(est, 0.975),
    max = max(est),
    ESS = coda::effectiveSize(coda::as.mcmc(est)),
    GewekeZ = coda::geweke.diag(coda::as.mcmc(est))[[1]],
    GewekeP = dnorm(GewekeZ)
  )











#......................
# cumulative
#......................
dat <- COVIDCurve::Aggsim_infxn_2_death(
  casefat = casefat,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 150,
  level = "Cumulative",
  infections = infxns$infxns
)

mod1 <- make_modinf_agg$new()
mod1$set_level("Cumulative")
mod1$set_data(dat)
mod1$set_IFRparams(c("r1", "r2", "ma3"))
mod1$set_Infxnparams(c("y1", "y2", "y3", "y4", "y5", "y6"))
mod1$set_paramdf(df_params)
mod1$set_pa(c(1/3, 1/3, 1/3))
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_knots(c(1, 30, 60, 90, 120, 150))
r_mcmc_out.cumm <- COVIDCurve::run_modinf_agg(modinf = mod1, reparamIFR = T, rungs = 10)
plot_par(r_mcmc_out.cumm, "y1")



# sanity
