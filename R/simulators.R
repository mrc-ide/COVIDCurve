#' @title Simulate Cumulative Expected Deaths for a Given Day
#' @param I0 integer; number of infected individuals at the beginning of the epidemic
#' @param r double; the expontential growth rate of the epidemic
#' @param m_od double; the mean of the onset of infection to death (gamma distribution)
#' @param s_od double; the coefficient of variation of the onset of infection to death (gamma distribution)
#' @param casefat dataframe; column names age, cfr, and pa that correspond to age-band, case-fatality ratio, and the attack rate in given age band, respectively
#' @param curr_day numeric; current day of epidemic
#' @importFrom magrittr %>%
#' @export

sim_infxn_2_death_cumulative <- function(casefat, I0, r, curr_day, m_od = 18.8, s_od = 0.45){

  # get expected number of infections up to today
  total_expected_inf <- casefat$pa * I0/r * exp(r*curr_day)

  # get expected deaths without censoring
  total_expected_death <- total_expected_inf * casefat$cfr

  # get actual number of deaths without censoring
  total_death_uncensored <- rpois(length(total_expected_death),
                                  lambda = total_expected_death)

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


#' @title Simulate Expected Deaths in a Time-Series
#' @inheritParams sim_infxn_2_death_cumulative
#' @importFrom magrittr %>%
#' @export

sim_infxn_2_death_timeseries <- function(casefat, I0, r, curr_day, m_od = 18.8, s_od = 0.45){

  # get expected number of infections up to today
  total_expected_inf <- casefat$pa * I0/r * exp(r*curr_day)

  # get expected deaths without censoring
  total_expected_death <- total_expected_inf * casefat$cfr

  # get actual number of deaths without censoring
  total_death_uncensored <- rpois(length(total_expected_death),
                                  lambda = total_expected_death)

  # distribute infections through time
  t_infection <- lapply(total_death_uncensored, function(x) { log(runif(x))/r + curr_day })

  # draw times of death
  t_death <- mapply(function(d_u, t_i){
    t_d <- t_i + rgamma(d_u, shape = 1/s_od^2, scale = m_od*s_od^2)
    return(t_d)
  }, d_u = total_death_uncensored, t_i = t_infection)

  # expand out the death grid
  casefat_adddeaths <- function(datrow, dths){
    n <- length(dths)
    datrow <- datrow[rep(1, times = n), ]
    datrow$t_death <- dths
    return(datrow)
  }
  case_line_list <- split(casefat[,c("age", "cfr")], 1:nrow(casefat))
  case_line_list <- mapply(casefat_adddeaths, datrow = case_line_list,
                           dths = t_death, SIMPLIFY = F) %>%
    do.call("rbind.data.frame", .)


  # Tidy up so that we observe deaths on a daily time step
  agelvls <- unique(as.character(casefat$age)) # protect against data.frame, string as factor = F
  lowertail <- floor(min(unlist(t_death)))
  death_line_list <- case_line_list %>%
    dplyr::mutate(age = factor(age, levels = agelvls)) %>%  # need this for later summarize
    dplyr::select(c("age", "t_death")) %>%
    dplyr::mutate(tod_fct = cut(t_death, breaks = c((lowertail-1), lowertail:curr_day),
                                labels = lowertail:curr_day,
                                left = F, right = T)) %>%
    dplyr::filter(!is.na(tod_fct)) %>% # drop "future" deaths
    dplyr::group_by(tod_fct, age, .drop = F) %>%
    dplyr::summarise(
      day_deaths = n()) %>%
    dplyr::ungroup(tod_fct, age) %>%
    dplyr::group_by(age) %>%
    dplyr::mutate(obs_day = as.numeric( levels(tod_fct)) # char to numeric -- fine
    ) %>%
    dplyr::ungroup(age)  # for export

  # out
  return(death_line_list)
}

