#' @title Simualte Expected Deaths
#' @param t0 integer; number of infected individuals at the beginning of the epidemic
#' @param r double; the expontential growth rate of the epidemic
#' @param alpha double; the mean of the onset of infection to death (gamma distribution)
#' @param beta double; the coefficient of variation of the onset of infection to death (gamma distribution)
#' @param casefat
sim_infxn_2_death <- function(t0, r, # for exp growth
                              alpha, beta, # for gamma distribution
                              casefat, # case fatality and age bands
                              curr_day){

  # draw number of infections
  exp.integrand <- function(t, gr = r) { return(t0*exp(gr*t)) }
  num_infxn <- integrate(exp.integrand, lower = 0, upper = curr_day)
  num_infxn <- round(num_infxn$value) # round to nearest person

  # assume that the time infection is exp distrib
  t_probs <- sapply(1:curr_day, exp.integrand)
  t_probs <- t_probs/sum(t_probs) # standardize
  t_infxn <- sample(1:curr_day, size = num_infxn, prob = t_probs, replace = T)

  # infections occurring equally over all age groups, "multinomial draw"
  casefat$n_infxn_t <- as.vector( rmultinom(n = 1, size = num_infxn, prob = rep(1, nrow(casefat))/nrow(casefat) ) )

  # now that we have number infected across groups, determine whether each case dies
  case_fatal <- rep(casefat$cfr, times = casefat$n_infxn_t)
  case_fatal <- sapply(case_fatal, function(x){return(as.logical(rbinom(1, 1, x)))})
  # do we observe the death in this interval
  death_times <- rgamma(num_infxn, alpha, beta)
  # total deaths we observe at this time point
  die_in_interval <- (death_times + t_infxn) < curr_day
  # record deaths at time t
  obs_deaths_t <- rep(F, times = num_infxn)
  obs_deaths_t[die_in_interval & case_fatal] <- T
  names(obs_deaths_t) <- rep(casefat$age, times = casefat$n_infxn_t)
  obs_deaths_t <- table(factor(names(obs_deaths_t)[obs_deaths_t], levels = casefat$age))

  ret <- list(
    observed_deaths = obs_deaths_t,
    case_fatality_ratio_by_age = casefat
  )

  return(ret)
}
