#########################################################################
# Purpose: pull sims together
#
# Author: Nicholas F. Brazeau
#
# Date: April 20 2020
#########################################################################
library(tidyverse)

#..............................................................
# read in
#..............................................................
ret <- readr::read_tsv("ignore/deploy/cluster_deploy/paramset/mastermap_map.txt") %>%
  dplyr::mutate(modelit = paste0(model, "_iter_", iter)) %>%
  dplyr::filter(lvl == "Time-Series")

filepaths <- tibble::tibble(
  modelit = gsub(".rds", "", basename(list.files("~/Desktop/curvawr_mcmc_results/", full.names = T))),
  path = list.files("~/Desktop/curvawr_mcmc_results/", full.names = T)
)
ret <- dplyr::left_join(ret, filepaths, by = "modelit") %>%
  dplyr::filter(!is.na(path))

#..............................................................
# temp debug
#..............................................................

ret$drjacoby <- purrr::map(ret$path, function(x){readRDS(x)})
ret$ma1_summ <- "ma1"
ret$ma2_summ <- "ma2"
ret$I0_summ <- "I0"

get_summ_results <- function(djobj, param){
  switch(param,
    "ma1"={
      djobj$output %>%
        dplyr::filter(stage == "sampling") %>%
        dplyr::select(c("chain", "ma1")) %>%
        tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
        dplyr::group_by(chain, param) %>%
        dplyr::mutate(ESS = coda::effectiveSize(coda::as.mcmc(est)),
                      gelmansz = coda::geweke.diag(coda::as.mcmc(est))$z) %>%
        dplyr::summarise(
          min = min(est),
          LCI = quantile(est, 0.025),
          median = median(est),
          mean = mean(est),
          UCI = quantile(est, 0.975),
          max = max(est),
          ESS = unique(ESS),
          gelmansz = unique(gelmansz)
        ) %>%
        dplyr::mutate(trueest = 0.1) %>%
        dplyr::ungroup(chain, param) %>%
        dplyr::mutate(CI_contains_truth = ifelse(LCI <= trueest & UCI >= trueest, T, F))

    },
    "ma2"={
      djobj$output %>%
        dplyr::filter(stage == "sampling") %>%
        dplyr::select(c("chain", "ma2")) %>%
        tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
        dplyr::group_by(chain, param) %>%
        dplyr::mutate(ESS = coda::effectiveSize(coda::as.mcmc(est)),
                      gelmansz = coda::geweke.diag(coda::as.mcmc(est))$z) %>%
        dplyr::summarise(
          min = min(est),
          LCI = quantile(est, 0.025),
          median = median(est),
          mean = mean(est),
          UCI = quantile(est, 0.975),
          max = max(est),
          ESS = unique(ESS),
          gelmansz = unique(gelmansz)
        ) %>%
        dplyr::mutate(trueest = 0.5) %>%
        dplyr::ungroup(chain, param) %>%
        dplyr::mutate(CI_contains_truth = ifelse(LCI <= trueest & UCI >= trueest, T, F))

    },
    "I0"={
      djobj$output %>%
        dplyr::filter(stage == "sampling") %>%
        dplyr::select(c("chain", "I0")) %>%
        tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
        dplyr::group_by(chain, param) %>%
        dplyr::mutate(ESS = coda::effectiveSize(coda::as.mcmc(est)),
                      gelmansz = coda::geweke.diag(coda::as.mcmc(est))$z) %>%
        dplyr::summarise(
          min = min(est),
          LCI = quantile(est, 0.025),
          median = median(est),
          mean = mean(est),
          UCI = quantile(est, 0.975),
          max = max(est),
          ESS = unique(ESS),
          gelmansz = unique(gelmansz)
        ) %>%
        dplyr::mutate(trueest = 2) %>%
        dplyr::ungroup(chain, param) %>%
        dplyr::mutate(CI_contains_truth = ifelse(LCI <= trueest & UCI >= trueest, T, F))

    }
  )

}

ret$ma1_summ <- purrr::map2(ret$drjacoby, ret$ma1_summ, get_summ_results)
ret$ma2_summ <- purrr::map2(ret$drjacoby, ret$ma2_summ, get_summ_results)
ret$I0_summ <- purrr::map2(ret$drjacoby, ret$I0_summ, get_summ_results)

#..............................................................
# tidy up results
#..............................................................
ret.tidy <- ret %>%
  dplyr::select(-c("modelit", "path", "drjacoby")) %>%
  tidyr::gather(., key = "paramsumms", value = "df", 9:11) %>%
  tidyr::unnest(cols = "df")


ret.summary <- ret.tidy %>%
  dplyr::mutate(sampling = 1e3) %>%
  dplyr::group_by(model, lvl, curr_day, prior, burnin, sampling, rungs, gti) %>%
  dplyr::summarise(
    nchains = sum(param == "ma1"),
    meanESS = mean(ESS),
    mean_CI_truth = mean(CI_contains_truth)
  )

ret.summary %>%
  dplyr::mutate_if(is.numeric, round, 2) %>%
  write_csv(., "~/Desktop/time_series_summary.csv")
