#########################################################################
# Purpose: To simulate under the time series model and perform inference
#
# Author: Nicholas F. Brazeau
#
# Date: April 15 2020
#########################################################################
#devtools::install_github("mrc-ide/drjacoby", ref = "develop")
library(drjacoby)
library(tidyverse)
source("R/sim_expd_timeseries.R")

#..................
# Run simulator
#..................
# make up a case-fatality
casefat <- data.frame(age = c("1:60", "60:100"),
                      cfr = c(0.1, 0.5),
                      pa = 1/2)

# define key parameters
I0 <- 2
min_day <- 1
curr_day <- 50

# simulate data
dat <- sim_infxn_2_death_timeseries(
  casefat = casefat,
  I0 = I0,
  r = 0.14,
  m_od = 18.8,
  s_od = 0.45,
  min_day = min_day,
  curr_day = curr_day
)

# split in to days
dat <- split(dat, factor(dat$obs_day))
dat <- purrr::map(dat, "day_deaths")

# list for Dr. Jacoby
data_list <- list(obs_deaths = dat)

#..................
# Run Dr. Jacoby w/ R-LL
#..................
source("R/R_likelihood_timeseries.R")

# params
df_params <- rbind.data.frame(list("I0", I0, I0, I0),  # fixed parameter
                              list("r1", 0, 100, 2),
                              list("ma2", 0, 1, 0.1))

names(df_params) <- c("name", "min", "max", "init")

# create list of misc elements to pass to MCMC
misc_list <- list(min_day = min_day,
                  curr_day = curr_day,
                  pa = casefat$pa)

# define MCMC parameters
burnin <- 1e2
samples <- 1e2
chains <- 1

# MCMC
r_mcmc_out <- run_mcmc(data = data_list,
                       df_params = df_params,
                       misc = misc_list,
                       loglike = r_tod_log_like_timeseries,
                       logprior = r_tod_log_prior_timeseries,
                       burnin = burnin,
                       samples = samples,
                       chains = chains,
                       pb_markdown = FALSE)


# append Dr. Jacoby output with reparameterized posteriors
r_mcmc_out$output$ma1 <- r_mcmc_out$output$ma2 * r_mcmc_out$output$r1

# parameter plots
plot_par(r_mcmc_out, "ma1")
plot_par(r_mcmc_out, "ma2")




# plots
library(patchwork)
sapply(c(1,2), function(x){
  plotObj <- plot_par(r_mcmc_out, show = paste0("ma", x),
                      phase = "sampling", display = F)
  plotObj[[1]][[2]] <- plotObj[[1]][[2]] +
    geom_vline(xintercept = casefat$cfr[x], color = "red", size = 2)
  plotObj <- (plotObj[[1]][[1]]) / (plotObj[[1]][[2]] | plotObj[[1]][[3]])

  jpeg(paste0("~/Desktop/scalar_ma", x, "_.jpg"), height = 11, width = 8, units = "in", res = 500)
  plot(plotObj)
  graphics.off()
})


jpeg("~/Desktop/corrma.jpg", height = 11, width = 8, units = "in", res = 500)
plot_cor(r_mcmc_out, "ma2", "ma1")
graphics.off()


#..................
# Run Dr. Jacoby w/ Cpp-LL
#..................
source("R/cpp_likelihood_timeseries.R")
raw <- list(obs_deaths = unname(unlist(data_list)))

# MCMC
# NB -- I am very slow ~ 30 minutes per chain
r_mcmc_out <- run_mcmc(data = raw,
                       df_params = df_params,
                       misc = list(min_day = 0, curr_day = 50, pa = c(0.5, 0.5), nsteps = 1e4),
                       loglike = cpp_tod_log_like_timeseries,
                       logprior = cpp_tod_log_prior_timeseries,
                       burnin = 1e3,
                       samples = 1e3,
                       chains = 3,
                       pb_markdown = TRUE)
# append Dr. Jacoby output with reparameterized posteriors
r_mcmc_out$output$ma1 <- r_mcmc_out$output$ma2 * r_mcmc_out$output$r1
plot_par(r_mcmc_out, "ma1")




