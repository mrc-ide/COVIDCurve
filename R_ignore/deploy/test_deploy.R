#########################################################################
# Purpose: Deploy
#
# Author: Nicholas F. Brazeau
#########################################################################
set.seed(1234)
library(drjacoby)
devtools::load_all()
# sigmoidal function
infxns <- data.frame(time = 1:150)
sig <- function(x){1 / (1 +  exp(-x))}
timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                            min = -25,
                                            max = 50)
sum(infxns$infxns < 0)

# make up fatality data
fatalitydata <- data.frame(strata = c("ma1", "ma2", "ma3"),
                           ifr = c(0.05, 0.2, 0.5),
                           pa = 1/3)
# pick serology date
sero_day <- 75

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
  curr_day = 150,
  level = "Time-Series",
  infections = infxns$infxns,
  simulate_seroprevalence = TRUE,
  sensitivity = 0.8,
  specificity = 0.95,
  sero_delay_rate = 10,
  popN = 5e6
)


dat <- list(obs_deaths = dat$AggDat,
            obs_serologyrate = dat$seroprev$SeroRateFP[sero_day])

#..................
# make model
#..................
ifr_paramsdf <- tibble::tibble(name = c("r1", "r2",  "ma3"),
                            min  = rep(0, 3),
                            init = rep(0.5, 3),
                            max = rep(1, 3),
                            dsc1 = rep(0, 3),
                            dsc2 = rep(0, 3))
infxn_paramsdf <- tibble::tibble(name = paste0("y", 1:5),
                                 min  = rep(0, 5),
                                 init = c(0.01, rep(0.5, 3), 1e3),
                                 max =  c(0.05, rep(1, 3),   1e4),
                                 dsc1 = c(0.01, rep(0, 4)),
                                 dsc2 = c(0.05, rep(1, 3),   1e4))
knot_paramsdf <- tibble::tibble(name = paste0("x", 1:4),
                                 min  = c(0,    0.33, 0.66, 120),
                                 init = c(0.05, 0.40, 0.75, 135),
                                 max =  c(0.33, 0.66, 0.99, 150),
                                 dsc1 = c(0,    0.33, 0.66, 120),
                                 dsc2 = c(0.33, 0.66, 0.99, 150))
sero_paramsdf <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day"),
                                min =   c(0.75,    0.92,   10,          75),
                                init =  c(0.8,     0.95,   10,          75),
                                max =   c(0.9,     0.98,   10,          75),
                                dsc1 =  c(8000,     9500,    5,           70),
                                dsc2 =  c(2000,     500,     15,          80))

df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sero_paramsdf)

#......................
# make mode
#......................
mod1 <- make_IFRmodel_agg$new()
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(c("r1", "r2", "ma3"))
mod1$set_maxMa("ma3")
mod1$set_Knotparams(paste0("x", 1:4))
mod1$set_relKnot("x4")
mod1$set_Infxnparams(paste0("y", 1:5))
mod1$set_relInfxn("y5")
mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_day"))
mod1$set_popN(5e6)
mod1$set_paramdf(df_params)
mod1$set_pa(c(1/3, 1/3, 1/3))
mod1$set_rcensor_day(130)
#..................
# run model
#..................
start <- Sys.time()
modout <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod1,
                                       reparamIFR = TRUE,
                                       reparamInfxn = TRUE,
                                       reparamKnot = TRUE,
                                       burnin = 1e4,
                                       samples = 1e3,
                                       chains = 10)
Sys.time() - start
modout
plot_par(modout$mcmcout, "r1")
plot_par(modout$mcmcout, "sens")
plot_par(modout$mcmcout, "spec")
plot_par(modout$mcmcout, "sero_day")




plot_mc_acceptance(modout$mcmcout)
drjacoby::plot_rung_loglike(modout$mcmcout)
drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 1, y_axis_type = 2)
drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 1, y_axis_type = 3)
drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 2, y_axis_type = 2)
drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 2, y_axis_type = 3)
(ifr <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, whichrung = paste0("rung", 1),
                                       what = "IFRparams", by_chain = F))
(knotspost <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                             what = "Knotparams", by_chain = F))
(infxn <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                         what = "Infxnparams", by_chain = F))
(sero <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,  whichrung = paste0("rung", 1),
                                        what = "Seroparams", by_chain = F))


curve <- COVIDCurve::draw_posterior_infxn_points_cubic_splines(IFRmodel_inf = modout,
                                                               whichrung = paste0("rung", 1),
                                                               by_chain = F,
                                                               CIquant = 0.9)
# plot out
jpeg("~/Desktop/posterior_curve_draws.jpg", width = 11, height = 8, units = "in", res = 500)
library(ggplot2)
liftover <- data.frame(param = c("r1", "r2", "ma3"),
                       strata = c("ma1", "ma2", "ma3"))

fatalitydataplot <- fatalitydata %>%
  dplyr::left_join(liftover, ., by = "strata")

plot1 <- ggplot() +
  geom_pointrange(data = ifr, aes(x = param, ymin = LCI, ymax = UCI, y = median, color = param)) +
  geom_hline(data = fatalitydataplot, aes(yintercept  = ifr, group = param), color = "#3182bd", linetype = "dashed", size = 1.1) +
  facet_wrap(.~param) +
  scale_color_viridis_d() +
  theme_bw()

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
                                                          CIquant = 0.9,
                                                          by_chain = FALSE)
postdeaths.plotObj <- postdeaths %>%
  dplyr::select(c("time", dplyr::starts_with("deaths"))) %>%
  tidyr::gather(., key = "strata", value = "deaths", 2:ncol(.)) %>%
  ggplot() +
  geom_line(aes(x= time, y = deaths, group = strata, color = strata), size = 1.2) +
  scale_color_viridis_d()

jpeg("~/Desktop/posterior_check.jpg", width = 11, height = 8, units = "in", res = 500)
postdeaths.plotObj +
  geom_line(data = dat$obs_deaths, aes(x=ObsDay, y = Deaths, group = Strata), color = "#bdbdbd") +
  theme_bw() +
  ggtitle("Posterior Predictive Check", subtitle = "Grey Lines are Simulated Data, Viridis Lines are Draws from Posterior")
graphics.off()




plot_par(modout$mcmcout, "r1", rung = 1)
plot_par(modout, "r2", rung = 1)
plot_par(modout, "ma3", rung = 1)
plot_par(modout, "y1", rung = 1)
plot_par(modout, "y2", rung = 1)
plot_par(modout, "y3", rung = 1)
plot_par(modout, "y4", rung = 1)
plot_par(modout, "y5", rung = 1)
plot_par(modout, "y6", rung = 1)
plot_par(modout, "sens", rung = 1)
plot_par(modout, "spec", rung = 1)
plot_par(modout, "sero_day", rung = 1)

plot_mc_acceptance(modout)
plot_rung_loglike(modout)
plot_rung_loglike(modout, y_axis_type = 2)
plot_rung_loglike(modout, y_axis_type = 3)


plot_cor(modout, "r1", "y4")
plot_cor(modout, "ma3", "y4")



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
