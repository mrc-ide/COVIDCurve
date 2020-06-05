####################################################################################
## Purpose:
##
## Author: Nick Brazeau
##
## Date: 05 June, 2020
####################################################################################
library(tidyverse)

#......................
# spanish data
#......................
esp <- readr::read_csv("https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1") %>%
  dplyr::filter(countryterritoryCode == "ESP") %>%
  dplyr::mutate(dateRep = lubridate::dmy(dateRep),
    ObsDay = as.numeric(dateRep - min(dateRep)) + 1)
esp <- esp %>%
  dplyr::select(c("ObsDay", "deaths")) %>%
  dplyr::rename(Deaths = deaths) %>%
  dplyr::arrange(ObsDay)
plot(esp$Deaths)


#............................................................
# pull out Rcpp function
#...........................................................
#......................
# pull out dummy model
#......................
dat <- list(obs_deaths = esp,
            obs_serologyrate = 0.05)

ifr_paramsdf <- tibble::tibble(name = c("r1"),
                               min  = 0,
                               init = 0.5,
                               max = 1,
                               dsc1 = 0,
                               dsc2 = 1)
infxn_paramsdf <- tibble::tibble(name = paste0("y", 1:5),
                                 min  = rep(0, 5),
                                 init = c(0.01, rep(0.5, 3), 1e3),
                                 max =  c(0.05, rep(1, 3),   1e4),
                                 dsc1 = c(0.01, rep(0, 4)),
                                 dsc2 = c(0.05, rep(1, 3),   1e4))
knot_paramsdf <- tibble::tibble(name = paste0("x", 1:4),
                                min  = c(0,    0.33, 0.66, 120),
                                init = c(0.05, 0.40, 0.75, 135),
                                max =  c(0.33, 0.66, 0.99, 150),
                                dsc1 = c(0,    0.33, 0.66, 120),
                                dsc2 = c(0.33, 0.66, 0.99, 150))
sero_paramsdf <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day"),
                                min =   c(0.78,    0.93,   10,          70),
                                init =  c(0.8,     0.95,   10,          75),
                                max =   c(0.82,     0.97,   10,          80),
                                dsc1 =  c(8000,     9500,    5,           70),
                                dsc2 =  c(2000,     500,     15,          80))

df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sero_paramsdf)

#......................
# make mode
#......................
mod1 <- make_IFRmodel_agg$new()
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(c("r1"))
mod1$set_maxMa("r1")
mod1$set_Knotparams(paste0("x", 1:4))
mod1$set_relKnot("x4")
mod1$set_Infxnparams(paste0("y", 1:5))
mod1$set_relInfxn("y5")
mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_day"))
mod1$set_popN(sum(squire::get_population("Spain")$n))
mod1$set_paramdf(df_params)
mod1$set_pa(1)
mod1$set_rcensor_day(130)
#......................
# fit curve string
#......................
fitcurve_string <- COVIDCurve:::make_user_Agg_loglike(IFRmodel = mod1,
                                                      reparamIFR = FALSE,
                                                      reparamKnots = FALSE,
                                                      reparamInfxn = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, don't need it)
# pull out pieces I need
fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "bool spline_test = true;", n = 2)[,1]
#fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "std::vector\\<double\\> cumm_infxn_spline\\(infxn_spline.size\\(\\)\\);", n = 2)[,1]
fitcurve_string <- paste(fitcurve_start, fitcurve_curve, "Rcpp::List ret = Rcpp::List::create(infxn_spline); return ret;}", collapse = "")
Rcpp::cppFunction(fitcurve_string)



#......................
# inputs needed for cpp function
#......................
misc_list = list(pa = 1,
                 pgmms = c(1,1,1),
                 level = TRUE,
                 popN = 1,
                 rcensor_day = 1,
                 days_obsd = max(dat$obs_deaths$ObsDay),
                 n_knots = nrow(knot_paramsdf)+1)

extparams <- c("sens" = 1, "spec" = 1, "sero_day" = 1, "r1" = 1)
#............................................................
# run function w/ parameters in that matter
#...........................................................
xparams <- c(40, 60, 80, 100)
names(xparams) <- paste0("x", 1:4)
yparams <- c(0, 0, 0, 107, 743, 184)
names(yparams) <- paste0("y", 1:5)

loglike(params = c(extparams, xparams, yparams),
        param_i = 1,
        data = datin,
        misc = misc_list)


























# sanity

