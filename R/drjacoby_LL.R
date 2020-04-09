#' @title Likelihood for Curve Aware
r_tod_log_like <- function(params, param_i, data, misc) {
  # assume r is fixed
  r <- 0.14
  curr_day <- misc$curr_day

  # alpha and beta for time from infection to time to death
  # from lancet id paper
  m_od <- 18.8
  s_od <- 0.45
  alpha <- 1/(s_od*s_od)
  beta <- 1/(m_od*s_od*s_od)

  # free params
  I0 <- params[1]
  ma1 <- params[2]
  ma2 <- params[3]
  ma3 <- params[4]
  ma4 <- params[5]
  ma5 <- params[6]
  ma6 <- params[7]
  ma7 <- params[8]
  ma8 <- params[9]
  ma9 <- params[10]
  ma <- c(ma1, ma2, ma3, ma4, ma5, ma6, ma7, ma8, ma9)

  #..................
  # expected deaths
  #..................
  integrand <- function(t, talpha = alpha, tbeta = beta, gr = r){ return(
    I0 * exp(gr*t) * pgamma(curr_day - t, talpha, tbeta)) }
  integral <- integrate(integrand, lower = -Inf, curr_day)

  # total exp deaths
  exp.deaths <- misc$pa * ma * integral$value

  #..................
  # poisson
  #..................
  ret <- sum(dpois(x = data$obs_deaths, lambda = exp.deaths, log = T))
  return(ret)

}

#' @title Prior for Curve Aware
r_tod_log_prior <- function(params, param_i, misc) {
  I0 <- params[1]
  ma1 <- params[2]
  ma2 <- params[3]
  ma3 <- params[4]
  ma4 <- params[5]
  ma5 <- params[6]
  ma6 <- params[7]
  ma7 <- params[8]
  ma8 <- params[9]
  ma9 <- params[10]
  ma <- c(ma1, ma2, ma3, ma4, ma5, ma6, ma7, ma8, ma9)
  # pa * ma
  ma <- ma * misc$pa
  # flat prior
  ret <- dunif(I0, min = 0, max = 10, log = TRUE) +
    sum( sapply(ma, function(x){dunif(x, min = 0, max = 1, log = TRUE)}) )
  return(ret)
}


#..............................................................
# Run Dr. Jacoby
#..............................................................
#devtools::install_github("mrc-ide/drjacoby")
devtools::install_github("mrc-ide/drjacoby", ref = "develop")
library(drjacoby)

# define parameters dataframe
df_params <- data.frame(name = c("I0",
                                 "ma1", "ma2", "ma3", "ma4",
                                 "ma5", "ma6", "ma7", "ma8", "ma9"),
                        min = c(1, rep(0, 9)),
                        max = c(3, rep(1, 9)),
                        init = c(2, rep(0.5, 9)))

#obs_deaths <- list(obs_deaths = unname(as.numeric(obs_deaths)))
dat <- readRDS("~/Desktop/temp.rds")
casefat <- dat$case_fatality_ratio_by_age
data_list <- list(obs_deaths =  as.numeric(dat$observed_deaths))

r_mcmc_out <- run_mcmc(data = data_list,
                       df_params = df_params,
                       misc = list(curr_day = 50, pa = 1/9),
                       loglike = r_tod_log_like,
                       logprior = r_tod_log_prior,
                       burnin = 1e3,
                       samples = 1e3,
                       chains = 3,
                       pb_markdown = TRUE)



library(tidyverse)
library(patchwork)
postdata <- r_mcmc_out$output
cred_intervals <- postdata %>%
  dplyr::select(c("chain", "iteration", "I0", dplyr::starts_with("ma"))) %>%
  tidyr::gather(., key = "param", value = "est", 3:ncol(.)) %>%
  dplyr::group_by(chain, param) %>%
  dplyr::summarise(
    min = min(est),
    LCI = quantile(est, 0.025),
    median = median(est),
    mean = mean(est),
    UCI = quantile(est, 0.975),
    max = max(est)
  ) %>%
  dplyr::mutate_if(is.numeric, round, 2)

truth <- data.frame(param = c("I0", paste0("ma", 1:9)),
                    est = c(2, casefat$cfr))
ggplot() +
  geom_point(data = truth, aes(x = param, y = est)) +
  geom_pointrange(data = cred_intervals, aes(x=param, y = mean, ymin = LCI, ymax = UCI,
                                             color = factor(chain)),
                  alpha = 0.3) +
  theme_bw()


truthnoI <- truth %>%
  dplyr::filter(param != "I0")

cred_intervalsnoI <- cred_intervals %>%
  dplyr::filter(param != "I0")
ggplot() +
  geom_point(data = truthnoI, aes(x = param, y = est)) +
  geom_pointrange(data = cred_intervalsnoI, aes(x=param, y = mean, ymin = LCI, ymax = UCI,
                                                color = factor(chain)),
                  alpha = 0.3) +
  theme_bw()









plot_par(r_mcmc_out, show = "I0", phase = "sampling")
ma1 <- plot_par(r_mcmc_out, show = "ma1", phase = "sampling")
ma1$Plot_ma1[[2]] <- ma1$Plot_ma1[[2]] + geom_vline(xintercept = casefat$cfr[1], color = "red")
(ma1$Plot_ma1[[1]]) / (ma1$Plot_ma1[[2]] | ma1$Plot_ma1[[3]])
plot_credible(r_mcmc_out, show = paste0("ma", 1:9)) +
  geom_point(data = casefat, aes(x = age, y = cfr), color = "red")

plot_par(r_mcmc_out, show = "ma2", phase = "sampling")
plot_par(r_mcmc_out, show = "ma3", phase = "sampling")
plot_par(r_mcmc_out, show = "ma4", phase = "sampling")
plot_par(r_mcmc_out, show = "ma5", phase = "sampling")
ma5 <- plot_par(r_mcmc_out, show = "ma5", phase = "sampling")
ma5$Plot_ma5[[2]] <- ma5$Plot_ma5[[2]] + geom_vline(xintercept = casefat$cfr[5], color = "red")
(ma5$Plot_ma5[[1]]) / (ma5$Plot_ma5[[2]] | ma5$Plot_ma5[[3]])

plot_par(r_mcmc_out, show = "ma6", phase = "sampling")
plot_par(r_mcmc_out, show = "ma7", phase = "sampling")
plot_par(r_mcmc_out, show = "ma8", phase = "sampling")
ma8 <- plot_par(r_mcmc_out, show = "ma8", phase = "sampling")
ma8$Plot_ma8[[2]] <- ma8$Plot_ma8[[2]] + geom_vline(xintercept = casefat$cfr[8], color = "red")
(ma8$Plot_ma8[[1]]) / (ma8$Plot_ma8[[2]] | ma8$Plot_ma8[[3]])

plot_par(r_mcmc_out, show = "ma9", phase = "sampling")


casefat$cfr
drjacoby::plot_cor(x = r_mcmc_out, parameter1 = "ma1", parameter2 = "ma2")
drjacoby::plot_cor(x = r_mcmc_out, parameter1 = "ma3", parameter2 = "ma2")
drjacoby::plot_cor(x = r_mcmc_out, parameter1 = "ma4", parameter2 = "ma2")
drjacoby::plot_cor(x = r_mcmc_out, parameter1 = "ma5", parameter2 = "ma2")
drjacoby::plot_cor(x = r_mcmc_out, parameter1 = "ma6", parameter2 = "ma2")
drjacoby::plot_cor(x = r_mcmc_out, parameter1 = "ma7", parameter2 = "ma2")
drjacoby::plot_cor(x = r_mcmc_out, parameter1 = "ma8", parameter2 = "ma2")
drjacoby::plot_cor(x = r_mcmc_out, parameter1 = "ma9", parameter2 = "ma2")







# end
