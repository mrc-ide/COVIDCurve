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
set.seed(44046)
r_mcmc_out <- COVIDCurve::run_modinf_linelist(modinf = mod1, reparamIFR = T, rungs = 1)
beepr::beep()

testcovid <- function(x){
  set.seed(x)
  cat(c(x, "\n"))
  ret <- COVIDCurve::run_modinf_linelist(modinf = mod1, reparamIFR = T, rungs = 1)
  return(ret)
}
testcovid(44)

COVIDCurve::run_modinf_linelist(modinf = mod1, reparamIFR = T, rungs = 1)

r_mcmc_out
plot_par(r_mcmc_out, "ma1")
plot_par(r_mcmc_out, "ma2")
plot_par(r_mcmc_out, "ma3")
plot_par(r_mcmc_out, "mod")
plot_par(r_mcmc_out, "sod")
plot_par(r_mcmc_out, "mor")
plot_par(r_mcmc_out, "sor")


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
  expgrowth = F
)

#..................
# make model
#..................
df_params <- tibble::tibble(name = c("r1", "r2", "ma3", paste0("y", 1:6)),
                            min = c(0, 0, 0, rep(0, 6)),
                            max = c(1, 1, 1, rep(10, 6)),
                            init = c(0.15, 0.2, 0.4, rep(4, 6)),
                            dsc1 = c(200, 400, 500, rep(3, 6)),
                            dsc2 = c(800, 600, 500, rep(5, 6)))

mod1 <- make_modinf_agg$new()
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(c("r1", "r2", "ma3"))
mod1$set_Infxnparams(c("y1", "y2", "y3", "y4", "y5", "y6"))
mod1$set_paramdf(df_params)
mod1$set_pa(c(1/3, 1/3, 1/3))
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_knots(c(1, 30, 60, 90, 120, 150))

#..................
# run model
#..................
set.seed(11347)
r_mcmc_out <- COVIDCurve::run_modinf_agg(modinf = mod1, reparamIFR = T, rungs = 1)
r_mcmc_out
plot_par(r_mcmc_out, "r1")















# sanity
