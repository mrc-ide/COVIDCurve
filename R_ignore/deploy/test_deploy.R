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
plot(infxns$infxns)
sum(infxns$infxns < 0)

#..............................................................
# Run Sim
#..............................................................
# make up casefat
casefat <- data.frame(age = c("0", "1", "2"),
                      cfr = c(0.1, 0.2, 0.5),
                      pa = 0.5)

dat <- COVIDCurve::sim_infxn_2_death(
  casefat = casefat,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 150,
  level = "Time-Series",
  infections = infxns$infxns,
  expgrowth = F
)

#..............................................................
# Make Model
#..............................................................
df_params <- tibble::tibble(name = c("r1", "r2", "ma3", paste0("y", 1:6)),
                            min = c(0, 0, 0, rep(0, 6)),
                            max = c(1, 1, 1, rep(10, 6)),
                            init = c(0.15, 0.2, 0.4, rep(4, 6)),
                            dsc1 = c(200, 400, 500, rep(3, 6)),
                            dsc2 = c(800, 600, 500, rep(5, 6)))

mod1 <- make_modinf$new()
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(c("r1", "r2", "ma3"))
mod1$set_Infxnparams(c("y1", "y2", "y3", "y4", "y5", "y6"))
mod1$set_paramdf(df_params)
mod1$set_pa(c(1/3, 1/3, 1/3))
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_knots(c(1, 30, 60, 90, 120, 150))

#..............................................................
# Run Model
#..............................................................

r_mcmc_out <- COVIDCurve::run_modinf(modinf = mod1, reparamIFR = T)
r_mcmc_out
plot_par(r_mcmc_out, "r1")














# sanity
