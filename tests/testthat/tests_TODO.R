devtools::load_all()
#............................................................
#
#...........................................................
IFRmodel_inf <- readRDS("../reestimate_covidIFR_analysis/results/ModFits/ESP_age_GTI3_rung50_burn1000_smpl1000.RDS")
mcmcout.nodes <-  IFRmodel_inf$mcmcout$output %>%
  dplyr::filter(stage == "sampling" & rung == "rung1" & chain == "chain4")


maxparams <- mcmcout.nodes[which(mcmcout.nodes$loglikelihood == max(mcmcout.nodes$loglikelihood)),]
maxparams <- maxparams[!names(maxparams) %in% c("chain", "rung", "iteration", "stage", "logprior", "loglikelihood")]


#............................................................
# posterior infxn curve
#...........................................................
fitcurve_string <- COVIDCurve:::make_user_Agg_loglike(IFRmodel = IFRmodel_inf$inputs$IFRmodel,
                                                      reparamIFR = FALSE,
                                                      reparamKnots = FALSE,
                                                      reparamInfxn = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
# pull out pieces I need
fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "double cum_infxn_check = 0.0;", n = 2)[,1]
fitcurve_string <- paste(fitcurve_start, fitcurve_curve,
                         "std::vector<std::vector<double>> infxn_spline_strata(days_obsd, std::vector<double>(stratlen));
                            for (int i = 0; i < days_obsd; i++) {
                              for (int a = 0; a < stratlen; a++) {
                                infxn_spline_strata[i][a] =  ne[a] * infxn_spline[i];
                              }
                           }",
                         "Rcpp::List ret = Rcpp::List::create(infxn_spline_strata); return ret;}",
                         collapse = "")
Rcpp::cppFunction(fitcurve_string)

#......................
# inputs needed for cpp function
#......................
misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
                 demog = IFRmodel_inf$inputs$IFRmodel$demog$popN,
                 rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
                 days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
                 n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams)+1,
                 n_sero_obs = length(IFRmodel_inf$inputs$IFRmodel$Serodayparams))

datin <- list("obs_deaths" = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
              "obs_serology" = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroPrev)


#......................
# split, run, recombine
#......................
cpp_function_wrapper <- function(params, data, misc) {
  paramsin <- unlist(params[c(IFRmodel_inf$inputs$IFRmodel$modparam,
                              IFRmodel_inf$inputs$IFRmodel$sodparam,
                              IFRmodel_inf$inputs$IFRmodel$IFRparams,
                              IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                              IFRmodel_inf$inputs$IFRmodel$Knotparams,
                              IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                              IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                              IFRmodel_inf$inputs$IFRmodel$Noiseparams)])

  # run efficient cpp code
  infxns <- loglike(params = paramsin,
                    param_i = 1,
                    data = datin,
                    misc = misc_list)[[1]]
  infxns <- infxns %>%
    do.call("rbind.data.frame", .)

  # out
  colnames(infxns) <- paste0("infxns_", IFRmodel_inf$inputs$IFRmodel$IFRparams)
  ret <- cbind.data.frame(time = 1:nrow(infxns),
                          infxns)
  return(ret)

}


demog <- IFRmodel_inf$inputs$IFRmodel$demog

maxCurveout <- cpp_function_wrapper(params = maxparams, data = datin, misc = misc_list)
maxCurveout <- maxCurveout %>%
  tidyr::gather(., key = "Strata", value = "infxns", 2:ncol(.)) %>%
  dplyr::mutate(Strata = gsub("infxns_", "", Strata))

maxCurveout %>%
  dplyr::group_by(Strata) %>%
  dplyr::summarise(
    totinfxns = sum(infxns)
  ) %>%
  dplyr::left_join(., demog, by = "Strata") %>%
  dplyr::mutate(
    IR = totinfxns/popN
  )




#............................................................
# posterior sero curve
#...........................................................
# internal function, liftover cpp likelihood to get infxn curve
# NOTE, this is extremely sensitive to the placements of the Cpp source file and therefore, is not generalizable
fitcurve_string <- COVIDCurve:::make_user_Agg_loglike(IFRmodel = IFRmodel_inf$inputs$IFRmodel,
                                                      reparamIFR = FALSE,
                                                      reparamKnots = FALSE,
                                                      reparamInfxn = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
# pull out pieces I need
fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "const double OVERFLO_DOUBLE = DBL_MAX/100.0;", n = 2)[,1]
fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
fitcurve_curve <- stringr::str_replace(fitcurve_curve, "if \\(cum_infxn_check <= popN\\) \\{", "")
fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "std::vector<std::vector<double>> sero_con_num\\(n_sero_obs, std::vector<double>\\(stratlen\\)\\);", n = 2)[,1]
fitcurve_string <- paste(fitcurve_start, fitcurve_curve,
                         "std::vector<std::vector<double>> full_sero_con_num(days_obsd, std::vector<double>(stratlen));
                            std::vector<std::vector<double>> RGfull_sero_con_num(days_obsd, std::vector<double>(stratlen));
                            // get cumulative hazard for study period
                            std::vector<double> cum_hazard(days_obsd);
                            for (int d = 0; d < days_obsd; d++) {
                              cum_hazard[d] = 1-exp((-(d+1)/sero_rate));
                            }

                            for (int i = 0; i < days_obsd; i++) {
                              for (int j = 0; j < stratlen; j++) {
                                // loop through and split infection curve by strata and by number of seroconversion study period
                                // note this cumulative, so loop through previous days
                                for (int d = 0; d <= i; d++) {
                                  int time_elapsed = days_obsd - d - 1;
                                  full_sero_con_num[i][j] += infxn_spline[d] * ne[j] * cum_hazard[time_elapsed];
                                }
                              }
                            }

                            // correct for spec/sens
                            for (int i = 0; i < days_obsd; i++) {
                              for (int j = 0; j < stratlen; j++) {
                                // Rogan-Gladen Estimator
                                double obs_prev = (full_sero_con_num[i][j]/demog[j]) * (spec + (sens-1)) - (spec-1);
                                RGfull_sero_con_num[i][j] = round(obs_prev * demog[j]);
                              }
                            }",
                         "Rcpp::List ret = Rcpp::List::create(full_sero_con_num,  RGfull_sero_con_num); return ret;}",
                         collapse = "")
Rcpp::cppFunction(fitcurve_string)

#......................
# inputs needed for cpp function
#......................
misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
                 demog = IFRmodel_inf$inputs$IFRmodel$demog$popN,
                 rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
                 days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
                 n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams)+1,
                 n_sero_obs = length(IFRmodel_inf$inputs$IFRmodel$Serodayparams))

datin <- list("obs_deaths" = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
              "obs_serology" = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroPrev)


#......................
# split, run, recombine
#......................
cpp_function_wrapper <- function(params, data, misc) {
  paramsin <- unlist(params[c(IFRmodel_inf$inputs$IFRmodel$modparam,
                              IFRmodel_inf$inputs$IFRmodel$sodparam,
                              IFRmodel_inf$inputs$IFRmodel$IFRparams,
                              IFRmodel_inf$inputs$IFRmodel$Infxnparams,
                              IFRmodel_inf$inputs$IFRmodel$Knotparams,
                              IFRmodel_inf$inputs$IFRmodel$Serotestparams,
                              IFRmodel_inf$inputs$IFRmodel$Serodayparams,
                              IFRmodel_inf$inputs$IFRmodel$Noiseparams)])

  seroprev_lists <- loglike(params = paramsin,
                            param_i = 1,
                            data = datin,
                            misc = misc_list)

  inf_sero_con_num <- seroprev_lists[[1]] %>%
    do.call("rbind.data.frame", .) %>%
    magrittr::set_colnames(paste0("inf_seroprev_", IFRmodel_inf$inputs$IFRmodel$IFRparams)) %>%
    dplyr::mutate(ObsDay = sort(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$ObsDay))) %>%
    dplyr::select(c("ObsDay", dplyr::everything()))

  RG_sero_con_num <- seroprev_lists[[2]] %>%
    do.call("rbind.data.frame", .) %>%
    magrittr::set_colnames(paste0("RG_inf_seroprev_", IFRmodel_inf$inputs$IFRmodel$IFRparams)) %>%
    dplyr::mutate(ObsDay = sort(unique(IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$ObsDay))) %>%
    dplyr::select(c("ObsDay", dplyr::everything()))


  ret <- dplyr::left_join(inf_sero_con_num, RG_sero_con_num, by = "ObsDay")
  return(ret)

}

mcmcout.node.rows <- split(mcmcout.nodes, 1:nrow(mcmcout.nodes))
mcmcout.nodes$seroprev <- purrr::map(mcmcout.node.rows, cpp_function_wrapper,
                                     data = datin, misc = misc_list)



# DONE below

# #............................................................
# # return full loglike
# #...........................................................
# # TODO
# IFRmodel_inf <- readRDS("~/Desktop/temp.rds")
# mcmcout.nodes <-  IFRmodel_inf$mcmcout$output
#
# maxparams <- mcmcout.nodes[which(mcmcout.nodes$loglikelihood == max(mcmcout.nodes$loglikelihood)),]
# maxparams <- maxparams[!names(maxparams) %in% c("chain", "rung", "iteration", "stage", "logprior", "loglikelihood")]
#
#
# # internal function, liftover cpp likelihood to get infxn curve
# # NOTE, this is extremely sensitive to the placements of the Cpp source file and therefore, is not generalizable
# fitcurve_string <- COVIDCurve:::make_user_Agg_loglike(IFRmodel = IFRmodel_inf$inputs$IFRmodel,
#                                                       reparamIFR = FALSE,
#                                                       reparamKnots = FALSE,
#                                                       reparamInfxn = FALSE) #NOTE, must be false because we re-parameterized the posterior already if reparameterization was requested (and if not, not needed)
# # pull out pieces I need
# fitcurve_start <- stringr::str_split_fixed(fitcurve_string, "double loglik = -OVERFLO_DOUBLE;", n = 2)[,1]
# fitcurve_start <- sub("SEXP", "Rcpp::List", fitcurve_start)
#
# fitcurve_curve <- stringr::str_split_fixed(fitcurve_string, "if \\(nodex_pass\\) \\{", n = 2)[,2]
# fitcurve_curve <- stringr::str_split_fixed(fitcurve_curve, "loglik = -OVERFLO_DOUBLE;", n = 2)[,1]
# fitcurve_curve <- sub("if \\(cum_infxn_check <= popN\\) \\{", "", fitcurve_curve)
# fitcurve_string <- paste(fitcurve_start, "double loglik = -OVERFLO_DOUBLE;",
#                          fitcurve_curve,
#                          "loglik = -OVERFLO_DOUBLE; }",
#                          "Rcpp::List ret = Rcpp::List::create(loglik, sero_con_num, death_loglik, sero_loglik, infxn_spline); return ret;}",
#                          collapse = "")
# Rcpp::cppFunction(fitcurve_string)
#
# #......................
# # inputs needed for cpp function
# #......................
# misc_list = list(rho = IFRmodel_inf$inputs$IFRmodel$rho,
#                  demog = IFRmodel_inf$inputs$IFRmodel$demog$popN,
#                  rcensor_day = IFRmodel_inf$inputs$IFRmodel$rcensor_day,
#                  days_obsd = IFRmodel_inf$inputs$IFRmodel$maxObsDay,
#                  n_knots = length(IFRmodel_inf$inputs$IFRmodel$Knotparams)+1,
#                  n_sero_obs = length(IFRmodel_inf$inputs$IFRmodel$Serodayparams))
#
# datin <- list("obs_deaths" = IFRmodel_inf$inputs$IFRmodel$data$obs_deaths$Deaths,
#               "obs_serology" = IFRmodel_inf$inputs$IFRmodel$data$obs_serology$SeroPrev)
# # wrapper function
# cpp_function_wrapper <- function(params, data, misc) {
#   paramsin <- unlist(params[c(IFRmodel_inf$inputs$IFRmodel$modparam,
#                               IFRmodel_inf$inputs$IFRmodel$sodparam,
#                               IFRmodel_inf$inputs$IFRmodel$IFRparams,
#                               IFRmodel_inf$inputs$IFRmodel$Infxnparams,
#                               IFRmodel_inf$inputs$IFRmodel$Knotparams,
#                               IFRmodel_inf$inputs$IFRmodel$Serotestparams,
#                               IFRmodel_inf$inputs$IFRmodel$Serodayparams,
#                               IFRmodel_inf$inputs$IFRmodel$Noiseparams)])
#
#   full_mod_out <- loglike(params = paramsin,
#                             param_i = 1,
#                             data = datin,
#                             misc = misc_list)
#   names(full_mod_out) <- c("loglik", "sero_con_num", "death_loglik", "sero_loglik", "infxn_spline")
#
#   return(full_mod_out)
#
# }
# out <- cpp_function_wrapper(params = maxparams, data = datin, misc = misc_list)
# out$loglik
# max(IFRmodel_inf$mcmcout$output$loglikelihood)
