#..............................................................
# Run Sim and set up data for Dr. Jacoby
#..............................................................
source("R/sim_infxns_2_death.R")
# make up casefat
casefat <- data.frame(age = paste0(seq(0, 80, by = 10), ":", seq(10, 90, 10)),
                      cfr = c(0.09, 0.025, 0.054, 0.047,
                              0.051, 0.13, 0.097, 0.021, 0.15),
                      pa = 1/9)

dat <- sim_infxn_2_death(
  casefat = casefat,
  I0 = 2,
  r = 0.14,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 50
)


data_list <- list(obs_deaths = unname(dat))

#..............................................................
# Run Dr. Jacoby
#..............................................................
#devtools::install_github("mrc-ide/drjacoby")
devtools::install_github("mrc-ide/drjacoby", ref = "develop")
library(drjacoby)
source("R/drjacoby_LL.R")

df_params <- data.frame(name = c("ma9", paste0("r", 1:8), "I0"),
                        min = c(0, rep(0, 8), 1),
                        max = c(1, rep(100, 8), 100),
                        init = c(0.5, rep(2, 8), 5))

r_mcmc_out <- run_mcmc(data = data_list,
                       df_params = df_params,
                       misc = list(curr_day = 50, pa = 1/9),
                       loglike = r_tod_log_like,
                       logprior = r_tod_log_prior,
                       burnin = 1e3,
                       samples = 1e3,
                       chains = 3,
                       rungs = 1,
                       pb_markdown = TRUE)

#..............................................................
# append Dr. Jacoby output with reparameterized posteriors
#..............................................................
reparcols <- grepl("r[0-9]", colnames(r_mcmc_out$output))
params_repar <- r_mcmc_out$output[, c(reparcols | colnames(r_mcmc_out$output) == "ma9")]
reparcols <- grepl("r[0-9]", colnames(params_repar))
params_repar <- apply(params_repar[, reparcols], 2,
                      function(x){params_repar$ma9 * x})
colnames(params_repar) <- paste0("ma", 1:sum(reparcols))
r_mcmc_out$output <- cbind.data.frame(r_mcmc_out$output, params_repar)

#..............................................................
# now plot
#..............................................................
library(tidyverse)
library(patchwork)
postdata <- r_mcmc_out$output
cred_intervals <- postdata %>%
  dplyr::select(c("chain", "iteration", dplyr::starts_with("ma"))) %>%
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

truth <- data.frame(param = c(paste0("ma", 1:9)),
                    est = c(casefat$cfr))
credintervals <- ggplot() +
  geom_point(data = truth, aes(x = param, y = est), size = 3) +
  geom_pointrange(data = cred_intervals, aes(x=param, y = mean, ymin = LCI, ymax = UCI,
                                             color = chain), size = 1.1,
                  alpha = 0.3) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(family = "Helvetica", face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica", face = "bold", size = 15),
        legend.position = "none")
credintervals

sapply(1:9, function(x){
  plotObj <- plot_par(r_mcmc_out, show = paste0("ma", x),
                      phase = "sampling", display = F)
  plotObj[[1]][[2]] <- plotObj[[1]][[2]] +
    geom_vline(xintercept = casefat$cfr[x], color = "red", size = 2)
  plotObj <- (plotObj[[1]][[1]]) / (plotObj[[1]][[2]] | plotObj[[1]][[3]])

  jpeg(paste0("~/Desktop/scalar_ma", x, "_.jpg"), height = 11, width = 8, units = "in", res = 500)
  plot(plotObj)
  graphics.off()
})






plot_par(r_mcmc_out, show = "r1", phase = "sampling")
