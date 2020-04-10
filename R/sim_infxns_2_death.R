#' @title Simualte Expected Deaths
#' @param I0 integer; number of infected individuals at the beginning of the epidemic
#' @param r double; the expontential growth rate of the epidemic
#' @param m_od double; the mean of the onset of infection to death (gamma distribution)
#' @param s_od double; the coefficient of variation of the onset of infection to death (gamma distribution)
#' @param casefat dataframe; column names age, cfr, and pa that correspond to age-band, case-fatality ratio, and proportion of population in given age band, respectively
#' @param curr_day numeric; current day of epidemic

sim_infxn_2_death <- function(casefat, I0, r, curr_day, m_od = 18.8, s_od = 0.45){

  # get expected number of infections up to today
  total_expected_inf <- casefat$pa * I0/r * exp(r*curr_day)

  # get expected deaths without censoring
  total_expected_death <- total_expected_inf * casefat$cfr

  # get actual number of deaths without censoring
  total_death_uncensored <- rpois(2, lambda = total_expected_death)

  # distribute infections through time
  t_infection <- lapply(total_death_uncensored, function(x) { log(runif(x))/r + curr_day })

  # draw times of death
  t_death <- mapply(function(d_u, t_i){
    t_d <- t_i + rgamma(d_u, shape = 1/s_od^2, scale = m_od*s_od^2)
    return(t_d)
  }, d_u = total_death_uncensored, t_i = t_infection)

  # get final number of deaths with censoring
  total_death <- sapply(t_death, function(x){sum(x <= curr_day)})

  # out
  names(total_death) <- casefat$age
  return(total_death)

}










# sanity
