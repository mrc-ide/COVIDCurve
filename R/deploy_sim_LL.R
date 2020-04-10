#..............................................................
# Run Sim and set up data for Dr. Jacoby
#..............................................................
source("R/sim_infxns_2_death.R")
# casefat <- data.frame(age = c("0:60", "60:100"),
#                       cfr = c(0.1, 0.5),
#                       pa = c(0.5, 0.5))

casefat <- data.frame(age = "0:100", cfr = 0.1, pa = 1)
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
# define parameters dataframe
# df_params <- data.frame(name = c("I0",
#                                  "ma1", "ma2", "ma3", "ma4",
#                                  "ma5", "ma6", "ma7", "ma8", "ma9"),
#                         min = c(1, rep(0, 9)),
#                         max = c(3, rep(1, 9)),
#                         init = c(2, rep(0.5, 9)))

df_params <- data.frame(name = c("I0", "ma1"),
                        min = c(1, 0),
                        max = c(10, 1),
                        init = c(5, 0.5))

r_mcmc_out <- run_mcmc(data = data_list,
                       df_params = df_params,
                       misc = list(curr_day = 50, pa = 1),
                       loglike = r_tod_log_like,
                       logprior = r_tod_log_prior,
                       burnin = 1e3,
                       samples = 1e3,
                       chains = 3,
                       rungs = 1,
                       pb_markdown = TRUE)


