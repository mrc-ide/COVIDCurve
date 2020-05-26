#########################################################################
# Purpose: Deploy
#
# Author: Nicholas F. Brazeau
#########################################################################
set.seed(1234)
library(drjacoby)
devtools::load_all()
# sigmoidal function
infxns <- data.frame(time = 1:149)
sig <- function(x){1 / (1 +  exp(-x))}
timevec <- seq(from = -5, to = 7, length.out = nrow(infxns))
infxns$infxns <- sig(timevec) * 5e3 + runif(n = nrow(infxns),
                                            min = -25,
                                            max = 50)
sum(infxns$infxns < 0)

# make up casefat
casefat <- data.frame(age = c("ma1", "ma2", "ma3"),
                      ifr = c(0.1, 0.2, 0.5),
                      pa = 0.5)
# pick serology date
sero_day <- 75

#..............................................................
# AGGREGATE
#..............................................................
#..................
# run sim
#..................
dat <- COVIDCurve::Aggsim_infxn_2_death(
  casefat = casefat,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 150,
  level = "Time-Series",
  infections = infxns$infxns,
  simulate_seroprevalence = TRUE,
  sero_sens = 0.8,
  sero_spec = 0.95,
  sero_delay_rate = 10,
  popN = 5e6
)


dat <- list(obs_deaths = dat$AggDat,
            obs_serologyrate = dat$seroprev$SeroRateFP[sero_day])

#..................
# make model
#..................
deathdf_params <- tibble::tibble(name = c("r1", "r2", "ma3", "y1", "y2", "y3", "y4", "y5", "y6"),
                                 min =  c(0,    0,    0,      1e-5, 6,    7,    9,    10,   10),
                                 max =  c(1,    1,    1,      5,    11,   13,   15,   15,   15),
                                 init = c(0.15, 0.2,  0.4,    2,    8.5,  10,   12,   12.5, 12.5),
                                 dsc1 = c(200,  400,  500,    2.5,  8.5,  10,   12,   12.5, 13),
                                 dsc2 = c(800,  600,  500,    5,    5,    5,    5,    5,    5))

serodf_params <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_date"),
                                 min =  c(0.75,    0.90,   0.1,        65),
                                 max =  c(0.85,    0.99,   0.1,        85),
                                 init = c(0.8,     0.95,   0.1,        75),
                                 dsc1 = c(800,     950,    100,        75),
                                 dsc2 = c(200,     50,     900,        3))

df_params <- rbind.data.frame(deathdf_params, serodf_params)

mod1 <- make_modinf_agg$new()
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(c("r1", "r2", "ma3"))
mod1$set_Infxnparams(c("y1", "y2", "y3", "y4", "y5", "y6"))
mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_date"))
mod1$set_popN(5e6)
mod1$set_paramdf(df_params)
mod1$set_pa(c(1/3, 1/3, 1/3))
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_knots(c(1, 30, 60, 90, 120, 150))

#..................
# run model
#..................
r_mcmc_out.ts <- COVIDCurve::run_modinf_agg(modinf = mod1, reparamIFR = T, rungs = 10)
r_mcmc_out.ts
plot_par(r_mcmc_out.ts, "r1")
plot_par(r_mcmc_out.ts, "r2")
plot_par(r_mcmc_out.ts, "ma3")
plot_par(r_mcmc_out.ts, "y1")
plot_par(r_mcmc_out.ts, "y2")
plot_par(r_mcmc_out.ts, "y3")
plot_par(r_mcmc_out.ts, "y4")
plot_par(r_mcmc_out.ts, "y5")
plot_par(r_mcmc_out.ts, "y6")
plot_par(r_mcmc_out.ts, "sero_date")




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



#..............................................................
# LINE LIST
#..............................................................
dat <- COVIDCurve::LineListsim_infxn_2_death(
  casefat = casefat,
  m_od = 18.8,
  s_od = 0.45,
  m_or = 24.7,
  s_or = 0.35,
  hospprob = 0.005,
  curr_day = 150,
  infections = infxns$infxns) %>%
  dplyr::filter(hosp == 1) %>%
  dplyr::select(-c("hosp"))

df_params <- tibble::tibble(name = c("ma1", "ma2", "ma3", "mod", "sod", "mor", "sor"),
                            min =  c(0,     0,     0,    10,    0,     20,     0),
                            max =  c(1,     1,     1,    30,    1,     40,     1),
                            init = c(0.15, 0.2,   0.4,   15,   0.5,    30,    0.5),
                            dsc1 = c(200, 400, 500,     19,    -3,    25,     -3),
                            dsc2 = c(800, 600, 500,     5,     1,      5,     1)
)

mod1 <- make_modinf_linelist$new()
mod1$set_data(dat)
mod1$set_DeathModparam("mod")
mod1$set_DeathSodparam("sod")
mod1$set_RecovMorparam("mor")
mod1$set_RecovSorparam("sor")
mod1$set_IFRparams(c("ma1", "ma2", "ma3"))
mod1$set_paramdf(df_params)

#..................
# run model
#..................
r_mcmc_out <- COVIDCurve::run_modinf_linelist(modinf = mod1, reparamIFR = T, rungs = 1)
r_mcmc_out
plot_par(r_mcmc_out, "ma1")
plot_par(r_mcmc_out, "ma2")
plot_par(r_mcmc_out, "ma3")
plot_par(r_mcmc_out, "mod")
plot_par(r_mcmc_out, "sod")
plot_par(r_mcmc_out, "mor")
plot_par(r_mcmc_out, "sor")















# sanity
