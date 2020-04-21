#!/usr/bin/env Rscript

#..............................................................
# Purpose of this script is to run the curve aware mcmc on the LL cluster
#
# Note, this funciton is not generalizable -- purpose is
# project specfici
#..............................................................
#remotes::install_github("nickbrazeau/CurveAware")
library(CurveAware)
library(optparse)

option_list=list(
  make_option(c("-m", "--mastermap"),
              type = "character",
              default = NULL,
              help = paste("start parameters path"),
              metavar = "character"),
  make_option(c("-O", "--output"),
              type = "character",
              default = NULL,
              help = "Output directory path.",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$mastermap)){
  print_help(opt_parser)
  stop("Missing input argument", call. = FALSE)
}

if(is.null(opt$output)){
  print_help(opt_parser)
  stop("Missing output argument", call. = FALSE)
}

#..............................................................
# setup
#..............................................................
library(CurveAware)
library(drjacoby)
# make up a case-fatality
casefat <- data.frame(age = c("1:60", "60:100"),
                      cfr = c(0.1, 0.5),
                      pa = 1/2)

#..............................................................
# extract priors and likelihood
#..............................................................
params <- readRDS(file = opt$mastermap)
params <- unlist(params)

switch(params["lvl"],
       "Cumulative"={
         source("LogLikePrior_Catalog/R_cumulative_likelihood.R")
         ll <- r_cumulative_likelihood
       },
       "Time-Series"={
         source("LogLikePrior_Catalog/R_timeseries_likelihood.R")
         ll <- r_timeseries_likelihood
       },
       {
         stop("Likelihood not found")
       }
)


switch(params["prior"],
       "stronger"={
         source("LogLikePrior_Catalog/R_strongerI0_prior.R")
         lp <- r_strongerI0_prior
       },
       "intermed"={
         source("LogLikePrior_Catalog/R_intermedI0_prior.R")
         lp <- r_intermedI0_prior
       },
       "strong"={
         source("LogLikePrior_Catalog/R_strongI0_prior.R")
         lp <- r_strongI0_prior
       },
       {
         stop("Prior not found")
       }
)
#..............................................................
# unpack and run sim and mcmc
#..............................................................
switch(params["lvl"],
       "Cumulative"={
         simdat <- CurveAware::sim_infxn_2_death_cumulative(
           casefat = casefat,
           I0 = 2,
           r = 0.14,
           m_od = 18.8,
           s_od = 0.45,
           curr_day = unname(as.numeric(params["curr_day"]))
         )
         # dr jacoby
         df_params <- data.frame(name = c("r1", "ma2", "I0"),
                                 min = c(0, 0, 0),
                                 init = c(2, 0.1, 2),
                                 max = c(10, 1, 10))
         misclist <- list(curr_day = unname(as.numeric(params["curr_day"])), pa = c(0.5, 0.5))

         r_mcmc_out <- CurveAware::wrap_drjacoby_mcmc(
           data = simdat,
           level = params["lvl"],
           reparameterization = T,
           df_params = df_params,
           misc = misclist,
           LogLike = ll,
           LogPrior = lp,
           burnin = unname(as.numeric(params["burnin"])),
           samples = 1e3,
           chains = 3,
           rungs = unname(as.numeric(params["rungs"])),
           GTI_power = unname(as.numeric(params["gti"])),
           pb_markdown = T
         )

       },
       "Time-Series"={
         simdat <- CurveAware::sim_infxn_2_death_timeseries(
           casefat = casefat,
           I0 = 2,
           r = 0.14,
           m_od = 18.8,
           s_od = 0.45,
           curr_day = unname(as.numeric(params["curr_day"]))
         )

         # params
         df_params <- data.frame(name = c("r1", "ma2", "I0"),
                                 min = c(0, 0, 0),
                                 init = c(2, 0.1, 2),
                                 max = c(10, 1, 10))
         misclist <- list(min_day = min(simdat$obs_day),
                          curr_day = max(simdat$obs_day), pa = c(0.5, 0.5))

         r_mcmc_out <- CurveAware::wrap_drjacoby_mcmc(
           data = simdat,
           level = params["lvl"],
           reparameterization = T,
           df_params = df_params,
           misc = misclist,
           LogLike = ll,
           LogPrior = lp,
           burnin = unname(as.numeric(params["burnin"])),
           samples = 1e3,
           chains = 3,
           rungs = unname(as.numeric(params["rungs"])),
           GTI_power = unname(as.numeric(params["gti"])),
           pb_markdown = F
         )
       },
       {
         stop("Level not found")
       }
)







#..............................................................
# out
#..............................................................
saveRDS(object = r_mcmc_out,
        file = opt$output)
