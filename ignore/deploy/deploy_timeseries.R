#########################################################################
# Purpose: To simulate under the time series model and perform inference
#
# Author: Nicholas F. Brazeau
#
# Date: April 15 2020
#########################################################################
library(CurveAware)
library(drjacoby)
#..................
# Run simulator
#..................
# make up a case-fatality
casefat <- data.frame(age = c("1:60", "60:100"),
                      cfr = c(0.1, 0.5),
                      pa = 1/2)

simdat <- CurveAware::sim_infxn_2_death_timeseries(
  casefat = casefat,
  I0 = 2,
  r = 0.14,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 100
)



#..................
# Input for Curve Aware
#..................
source("LogLikePrior_Catalog/R_timeseries_likelihood.R")
source("LogLikePrior_Catalog/R_strongI0_prior.R")

# params
df_params <- data.frame(name = c("r1", "ma2", "I0"),
                        min = c(0, 0, 0),
                        init = c(2, 0.1, 2),
                        max = c(10, 1, 10))
misclist <- list(min_day = min(simdat$obs_day), curr_day = max(simdat$obs_day), pa = c(0.5, 0.5))

r_mcmc_out <- CurveAware::wrap_drjacoby_mcmc(
  data = simdat,
  level = "Time-Series",
  reparameterization = T,
  df_params = df_params,
  misc = misclist,
  LogLike = r_timeseries_likelihood,
  LogPrior = r_strongI0_prior,
  burnin = 1e3,
  samples = 1e3,
  chains = 3,
  pb_markdown = TRUE
)

