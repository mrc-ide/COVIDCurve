#########################################################################
# Purpose: To simulate under the time series model and perform inference
#
# Author: Nicholas F. Brazeau
#
# Date: April 15 2020
#########################################################################
devtools::install_github("mrc-ide/drjacoby", ref = "develop")
library(drjacoby)
library(tidyverse)
source("R/sim_expd_cumulative.R")

#..................
# Run simulator
#..................
# make up a case-fatality
casefat <- data.frame(age = c("1:60", "60:100"),
                      cfr = c(0.1, 0.5),
                      pa = 1/2)

dat <- sim_infxn_2_death_cumulative(
  casefat = casefat,
  I0 = 2,
  r = 0.14,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 50
)


# list for Dr. Jacoby
data_list <- list(obs_deaths = unname(dat))

#..................
# Run Dr. Jacoby w/ R-LL
#..................
source("R/R_likelihood_cumulative.R")

# params
df_params <- data.frame(name = c("r1", "ma2"),
                        min = c(0, 0),
                        max = c(100, 1),
                        init = c(2, 0.1))
# MCMC
# Note the misc list, where current day and pa must be consistent with your simulation
r_mcmc_out <- run_mcmc(data = data_list,
                       df_params = df_params,
                       misc = list(curr_day = 50, pa = c(0.5, 0.5)),
                       loglike = r_tod_log_like_cumulative,
                       logprior = r_tod_log_prior_cumulative,
                       burnin = 1e3,
                       samples = 1e3,
                       chains = 3,
                       pb_markdown = TRUE)

# append Dr. Jacoby output with reparameterized posteriors
r_mcmc_out$output$ma1 <- r_mcmc_out$output$ma2 * r_mcmc_out$output$r1

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
source("R/cpp_likelihood_cumulative.R")

# MCMC
r_mcmc_out <- run_mcmc(data = data_list,
                       df_params = df_params,
                       misc = list(curr_day = 50, pa = c(0.5, 0.5), nsteps = 1e3),
                       loglike = cpp_tod_log_like_cumulative,
                       logprior = cpp_tod_log_prior_cumulative,
                       burnin = 1e3,
                       samples = 1e3,
                       chains = 3,
                       pb_markdown = TRUE)

# append Dr. Jacoby output with reparameterized posteriors
r_mcmc_out$output$ma1 <- r_mcmc_out$output$ma2 * r_mcmc_out$output$r1

# plots
plot_par(r_mcmc_out, "ma1")

