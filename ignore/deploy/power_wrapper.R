#########################################################################
# Purpose: Goal here is to  wrap the cumulative
#          simulator and the time-series simulator and run Dr. Jacoby
#          to estimate the "power" of each method
#
# Author: Nicholas F. Brazeau
#
# Date: April 14 2020
#########################################################################
library(tidyverse)
source("R/sim_expd_timeseries.R")
source("R/sim_expd_cumulative.R")
source("R/R_likelihood_timeseries.R")
source("R/R_likelihood_cumulative.R")
devtools::install_github("mrc-ide/drjacoby", ref = "develop")
library(drjacoby)

run_sims <- function(seed){

  casefat <- data.frame(age = c("1:60", "60:100"),
                        cfr = c(0.1, 0.5),
                        pa = 1/2)

set.seed(seed)
#..............................................................
# Time Series
#..............................................................

  time_simdat <- sim_infxn_2_death_timeseries(
    casefat = casefat,
    I0 = 2,
    r = 0.14,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 50)

  # split in to days
  time_dat <- split(time_simdat, factor(time_simdat$obs_day))
  time_dat <- purrr::map(time_dat, "day_deaths")
  # list for Dr. Jacoby
  data_list <- list(obs_deaths = time_dat)
  #..................
  # Run Dr. Jacoby
  #..................
  df_params <- data.frame(name = c("r1", "ma2", "I0"),
                          min = c(0, 0, 0),
                          max = c(10, 1, 4),
                          init = c(1, 0.1, 2))
  # create list of misc elements to pass to MCMC
  misc_list <- list(min_day = min(time_simdat$obs_day),
                    curr_day = max(time_simdat$obs_day),
                    pa = casefat$pa)

  r_mcmc_out_time <- run_mcmc(data = data_list,
                              df_params = df_params,
                              misc = misc_list,
                              loglike = r_tod_log_like_timeseries,
                              logprior = r_tod_log_prior_timeseries,
                              burnin = 1e3,
                              samples = 1e3,
                              chains = 3,
                              pb_markdown = TRUE)
  #..................
  # LiftOver Dr. Jacoby
  #..................
  r_mcmc_out_time$output$ma1 <- r_mcmc_out_time$output$ma2 * r_mcmc_out_time$output$r1
  
  #..................
  # out
  #..................
  ageband <- data.frame(param = c("ma1", "ma2", "I0"), trueest = c(casefat$cfr, 2))

  timeout <- r_mcmc_out_time$output %>%
    dplyr::filter(stage == "sampling") %>%
    dplyr::select(c("chain", "ma1", "ma2", "I0")) %>%
    tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
    dplyr::group_by(chain, param) %>%
    dplyr::summarise(
      min = min(est),
      LCI = quantile(est, 0.025),
      median = median(est),
      mean = mean(est),
      UCI = quantile(est, 0.975),
      max = max(est)
    ) %>%
    dplyr::left_join(., ageband, by = "param") %>%
    dplyr::ungroup(chain, param)


set.seed(seed)
#..............................................................
# Cumulative Day Data
#..............................................................
  dat <- sim_infxn_2_death_cumulative(
    casefat = casefat,
    I0 = 2,
    r = 0.14,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 50)

  data_list <- list(obs_deaths = unname(dat))
  #..................
  # Run Dr. Jacoby
  #..................
  df_params <- data.frame(name = c("r1", "ma2", "I0"),
                          min = c(0, 0, 0),
                          max = c(10, 1, 4),
                          init = c(1, 0.1, 2))

  r_mcmc_out_cumul <- run_mcmc(data = data_list,
                               df_params = df_params,
                               misc = list(curr_day = 50, pa = c(0.5, 0.5)),
                               loglike = r_tod_log_like_cumulative,
                               logprior = r_tod_log_prior_cumulative,
                               burnin = 1e3,
                               samples = 1e3,
                               chains = 3,
                               pb_markdown = TRUE)
  #..................
  # LiftOver Dr. Jacoby
  #..................
  r_mcmc_out_cumul$output$ma1 <- r_mcmc_out_cumul$output$ma2 * r_mcmc_out_cumul$output$r1
  #..................
  # out
  #..................
  cumulout <- r_mcmc_out_cumul$output %>%
    dplyr::filter(stage == "sampling") %>%
    dplyr::select(c("chain", "ma1", "ma2", "I0")) %>%
    tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
    dplyr::group_by(chain, param) %>%
    dplyr::summarise(
      min = min(est),
      LCI = quantile(est, 0.025),
      median = median(est),
      mean = mean(est),
      UCI = quantile(est, 0.975),
      max = max(est)
    ) %>%
    dplyr::left_join(., ageband, by = "param") %>%
    dplyr::ungroup(chain, param)



#..............................................................
# return
#.............................................................
compar <- dplyr::left_join(timeout, cumulout, by = c("chain", "param"))
colnames(compar) <- gsub("\\.x", "_time", colnames(compar))
colnames(compar) <- gsub("\\.y", "_cumul", colnames(compar))

ret <- list(compar = compar,
            cumul = r_mcmc_out_cumul,
            time = r_mcmc_out_time)

return(ret)
}

seeds <- sort(round(runif(5, 1, 100)))

ret <- parallel::mclapply(seeds, run_sims)



compar_df <- purrr::map(ret, "compar")

compar_df <- compar_df %>%
  dplyr::bind_rows(., .id = "run")



plotObj <- compar_df  %>%
  dplyr::mutate(run = factor(run, levels = 1:5),
                true_est = trueest_time) %>%
  ggplot() +
  geom_errorbar(aes(x = median_cumul, ymin = min_time, ymax = max_time), color = "#969696", size = 1.2) +
  geom_errorbarh(aes(xmin = min_cumul, xmax = max_cumul, y = median_time, height = 0), color = "#525252", size = 1.2) +
  geom_point(aes(x = median_cumul, y = median_time, color = run)) +
  scale_color_viridis_d("Run") +
  geom_hline(aes(yintercept = true_est), color = "red", linetype = "dashed", size = 0.5, alpha = 0.5) +
  geom_vline(aes(xintercept = true_est), color = "red", linetype = "dashed", size = 0.5, alpha = 0.5) +
  facet_wrap(~param, scales = "free") +
  xlab("Cumulative Median Est (Min-Max CI)") + ylab("Time-Series Median Est (Min-Max CI)") +
  labs(caption = "Runs are using \"same\" data.") +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0),
        legend.position = "bottom")

plotObj


jpeg("~/Desktop/pw_nfb.jpg", width = 11, height = 8, units = "in", res = 500)
plotObj
graphics.off()

for(i in 1:10){
  Sys.sleep(3)
  beepr::beep()
}