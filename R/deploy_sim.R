source("R/sim_infxns_2_death.R")
#..............................................................
# Model parameters
#..............................................................
# infxn rate
I0 <- t0 <- 2 # number infected at time 0 (day 0)
curr_day <- 50 # current day for the model
r <- 0.14
# CFR from lancet id
# casefat <- data.frame(age = paste0(seq(0, 80, by = 10), ":", seq(10, 90, 10)),
#                       cfr = c(0.000016, 0.0000695, 0.000309, 0.000844,
#                               0.00161, 0.00595, 0.0193, 0.0428, 0.078))

casefat <- data.frame(age = paste0(seq(0, 80, by = 10), ":", seq(10, 90, 10)),
                      cfr = rexp(9, 0.05))
casefat$cfr <- casefat$cfr/sum(casefat$cfr)

# alpha and beta for time from infection to time to death
# from lancet id paper
m_od <- 18.8
s_od <- 0.45
alpha <- 1/(s_od*s_od)
beta <- 1/(m_od*s_od*s_od)

#..............................................................
# Simulate
#..............................................................

simdeaths <- sim_infxn_2_death(t0 = I0,
                               r = r,
                               alpha = alpha,
                               beta = beta,
                               casefat = casefat,
                               curr_day = curr_day)
saveRDS(simdeaths, "~/Desktop/temp.rds")
