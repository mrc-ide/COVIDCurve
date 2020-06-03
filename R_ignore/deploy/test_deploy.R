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
  sensitivity = 1,
  specificity = 1,
  sero_delay_rate = 10,
  popN = 5e6
)


dat <- list(obs_deaths = dat$AggDat,
            obs_serologyrate = dat$seroprev$SeroRateFP[sero_day])

#..................
# make model
#..................
# deathdf_params <- tibble::tibble(name = c("r1", "r2", "ma3", paste0("y", 1:nrow(infxns))),
#                                  min =  c(0,    0,    0,      infxns$infxns-2),
#                                  init = c(0.5,  0.5,  0.5,    infxns$infxns),
#                                  max =  c(1,    1,    1,      infxns$infxns+2),
#                                  dsc1 = c(0,    0,    0,      infxns$infxns-2),
#                                  dsc2 = c(1,    1,    1,      infxns$infxns+2))


# deathdf_params <- tibble::tibble(name = c("r1", "r2", "ma3", "y1", "y2", "y3",  "y4",  "y5",   "y6"),
#                                  min =  c(0,    0,    0,      15,   300,  2228,  4527,  4982,   5012),
#                                  init = c(0.5,  0.5,  0.5,    17,   303,  2230,  4529,  4984,   5014),
#                                  max =  c(1,    1,    1,      18,   306,  2232,  4531,  4986,   5016),
#                                  dsc1 = c(0,    0,    0,      15,   300,  2228,  4527,  4982,   5012),
#                                  dsc2 = c(1,    1,    1,      18,   306,  2232,  4531,  4986,   5016))
# deathdf_params <- tibble::tibble(name = c("r1", "r2", "ma3", "y1", "y2", "y3",  "y4",  "y5",   "y6"),
#                                  min =  c(0,    0,    0,      2.6, 5.5,  7.69,  8.2,  8.3,   8.3),
#                                  init = c(0.5,  0.5,  0.5,    2.8, 5.7,  7.71,  8.4,  8.5,   8.5),
#                                  max =  c(1,    1,    1,      3,   5.9,  7.73,  8.6,  8.7,   8.7),
#                                  dsc1 = c(0,    0,    0,      2.6, 5.5,  7.69,  8.2,  8.3,   8.3),
#                                  dsc2 = c(1,    1,    1,      3,   5.9,  7.73,  8.6,  8.7,   8.7))

# deathdf_params <- tibble::tibble(name = c("r1", "r2", "ma3", paste0("y", 1:6)),
#                                  min =  c(0,    0,    0,     rep(0, 6)),
#                                  init = c(0.5,  0.5,  0.5,   rep(5, 6)),
#                                  max =  c(1,    1,    1,     rep(12, 6)),
#                                  dsc1 = c(0,    0,    0,     rep(0, 6)),
#                                  dsc2 = c(1,    1,    1,     rep(12, 6)))

# deathdf_params <- tibble::tibble(name = c("r1", "r2",  "ma3",  paste0("y", 1:15)),
#                                  min =  c(0.58,  0.78,  0.88,  rep(0, 15)),
#                                  init = c(0.6,   0.8,   0.9,   rep(5, 15)),
#                                  max =  c(0.63,  0.82,  0.92,  rep(1e4, 15)),
#                                  dsc1 = c(0,     0,     0,     rep(0,15)),
#                                  dsc2 = c(1,     1,     1,     rep(1e4,15)))

deathdf_params <- tibble::tibble(name = c("r1", "r2",  "ma3",  paste0("y", 1:5), "y6",  paste0("x", 1:4), "x5"),
                                 min =  c(0,     0,     0,     rep(0, 5),         2500,  rep(0, 4),         120),
                                 init = c(0.5,   0.5,   0.5,   rep(2, 5),         5000,  rep(0.5, 4),       130),
                                 max =  c(1,     1,     1,     rep(5, 5),         10000, rep(1, 4),         150),
                                 dsc1 = c(0,     0,     0,     rep(0, 5),         2500,  rep(0, 4),         120),
                                 dsc2 = c(1,     1,     1,     rep(5, 5),         10000, rep(1, 4),         150)
)

# deathdf_params <- tibble::tibble(name = c("r1", "r2",  "ma3",  paste0("y", 1:5), "y6"),
#                                  min =  c(0,     0,     0,     rep(0, 5),         2500),
#                                  init = c(0.5,   0.5,   0.5,   rep(2, 5),         5000),
#                                  max =  c(1,     1,     1,     rep(5, 5),         10000),
#                                  dsc1 = c(0,     0,     0,     rep(0, 5),         2500),
#                                  dsc2 = c(1,     1,     1,     rep(5, 5),         10000))
#
# knotdf_params <- tibble::tibble(name = c("x1", "x2",  "x4",  "x4", "x5", "x6"),
#                                 min =  c(0,     0,     0,      1,     1,  1),
#                                 init = c(0.5,   0.5,   0.5,    1,     1,  1),
#                                 max =  c(1,     1,     1,     1,     1,   1),
#                                 dsc1 = c(0,     0,     0,    1,     1,   1),
#                                 dsc2 = c(1,     1,     1,     1,     1,  1))


# serodf_params <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day"),
#                                 min =   c(0.75,    0.91,   10,          75),
#                                 init =  c(0.8,     0.95,   10,          75),
#                                 max =   c(0.85,    0.99,   10,          75),
#                                 dsc1 =  c(800,     950,    5,           65),
#                                 dsc2 =  c(200,     50,     15,          85))


serodf_params <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day"),
                                min =   c(1,    1,   10,          75),
                                init =  c(1,     1,   10,          75),
                                max =   c(1,    1,   10,          75),
                                dsc1 =  c(999,     999,    5,           65),
                                dsc2 =  c(1,     1,     15,          85))

#df_params <- rbind.data.frame(deathdf_params, knotdf_params, serodf_params)
df_params <- rbind.data.frame(deathdf_params, serodf_params)
#knots <- round(seq(from = 1, to = 150, length.out = sum(grepl("y[0-9]+", df_params$name))))
mod1 <- make_IFRmodel_agg$new()
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(c("r1", "r2", "ma3"))
mod1$set_maxMa("ma3")
mod1$set_Knotparams(paste0("x", 1:5))
mod1$set_relKnot("x5")
mod1$set_Infxnparams(paste0("y", 1:6))
mod1$set_relInfxn("y6")
mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_day"))
mod1$set_popN(5e6)
mod1$set_paramdf(df_params)
mod1$set_pa(c(1/3, 1/3, 1/3))
mod1$set_rcensor_day(130)
#..................
# run model
#..................
pout <- capture.output(COVIDCurve::run_IFRmodel_agg(IFRmodel = mod1,
                                                    reparamIFR = TRUE,
                                                    reparamInfxn = TRUE,
                                                    reparamKnot = TRUE,
                                                    burnin = 1e2,
                                                    samples = 1e2,
                                                    chains = 1,
                                                    rungs = 1))


r_mcmc_out.ts <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod1,
                                              reparamIFR = TRUE,
                                              reparamInfxn = TRUE,
                                              reparamKnot = TRUE,
                                              burnin = 1e3,
                                              samples = 1e3,
                                              chains = 10,
                                              rungs = 1)

r_mcmc_out.ts
(ifr <- COVIDCurve::get_cred_intervals(IFRmodel = mod1, mcmcout = r_mcmc_out.ts, whichrung = paste0("rung", 1),
                                       what = "IFRparams", by_chain = F))
(knotspost <- COVIDCurve::get_cred_intervals(IFRmodel = mod1, mcmcout = r_mcmc_out.ts,  whichrung = paste0("rung", 1),
                                             what = "Knotparams", by_chain = F))
(infxn <- COVIDCurve::get_cred_intervals(IFRmodel = mod1, mcmcout = r_mcmc_out.ts,  whichrung = paste0("rung", 1),
                                         what = "Infxnparams", by_chain = F))
(sero <- COVIDCurve::get_cred_intervals(IFRmodel = mod1, mcmcout = r_mcmc_out.ts,  whichrung = paste0("rung", 1),
                                        what = "Seroparams", by_chain = F))


curve <- COVIDCurve::draw_posterior_infxn_points_cubic_splines(IFRmodel = mod1,
                                                               whichrung = paste0("rung", 1),
                                                               mcmcout = r_mcmc_out.ts,
                                                               by_chain = F,
                                                               CIquant = 0.95)
# plot out
jpeg("~/Desktop/posterior_curve_draws.jpg", width = 11, height = 8, units = "in", res = 500)

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
postdeaths <- COVIDCurve::posterior_check_infxns_to_death(IFRmodel = mod1,
                                                          mcmcout = r_mcmc_out.ts,
                                                          CIquant = 0.9,
                                                          by_chain = FALSE)
postdeaths.plotObj <- postdeaths %>%
  dplyr::select(c("time", dplyr::starts_with("deaths"))) %>%
  tidyr::gather(., key = "strata", value = "deaths", 2:ncol(.)) %>%
  ggplot() +
  geom_line(aes(x= time, y = deaths, group = strata, color = strata), size = 1.2) +
  scale_color_viridis_d()

postdeaths.plotObj +
  geom_line(data = dat$obs_deaths, aes(x=ObsDay, y = Deaths, group = AgeGroup), color = "#bdbdbd") +
  theme_bw() +
  ggtitle("Posterior Predictive Check", subtitle = "Grey Lines are Simulated Data, Viridis Lines are Draws from Posterior")





plot_par(r_mcmc_out.ts, "r1", rung = 1)
plot_par(r_mcmc_out.ts, "r2", rung = 1)
plot_par(r_mcmc_out.ts, "ma3", rung = 1)
plot_par(r_mcmc_out.ts, "y1", rung = 1)
plot_par(r_mcmc_out.ts, "y2", rung = 1)
plot_par(r_mcmc_out.ts, "y3", rung = 1)
plot_par(r_mcmc_out.ts, "y4", rung = 1)
plot_par(r_mcmc_out.ts, "y5", rung = 1)
plot_par(r_mcmc_out.ts, "y6", rung = 1)
plot_par(r_mcmc_out.ts, "sens", rung = 1)
plot_par(r_mcmc_out.ts, "spec", rung = 1)
plot_par(r_mcmc_out.ts, "sero_day", rung = 1)

plot_mc_acceptance(r_mcmc_out.ts)
plot_rung_loglike(r_mcmc_out.ts)
plot_rung_loglike(r_mcmc_out.ts, y_axis_type = 2)
plot_rung_loglike(r_mcmc_out.ts, y_axis_type = 3)


plot_cor(r_mcmc_out.ts, "r1", "y4")
plot_cor(r_mcmc_out.ts, "ma3", "y4")



#............................................................
# look at rungs
#...........................................................
COVIDCurve::get_cred_intervals(IFRmodel = mod1,
                               mcmcout = r_mcmc_out.ts,
                               whichrung = "rung50",
                               what = "IFRparams",
                               by_chain = F)

r_mcmc_out.ts$output %>%
  dplyr::mutate(
    logpost = loglikelihood + logprior
  ) %>%
  dplyr::filter(stage == "sampling" &
                  rung == "rung50") %>%
  ggplot() +
  geom_histogram(aes(x  = logpost))

rungdf <- r_mcmc_out.ts$output %>%
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
