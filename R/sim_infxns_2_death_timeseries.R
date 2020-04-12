#' @title Simualte Expected Deaths
#' @param I0 integer; number of infected individuals at the beginning of the epidemic
#' @param r double; the expontential growth rate of the epidemic
#' @param m_od double; the mean of the onset of infection to death (gamma distribution)
#' @param s_od double; the coefficient of variation of the onset of infection to death (gamma distribution)
#' @param casefat dataframe; column names age, cfr, and pa that correspond to age-band, case-fatality ratio, and proportion of population in given age band, respectively
#' @param min_day numeric; day recording of epidemic started (after true start of "time 0")
#' @param curr_day numeric; current day of epidemic
#' @importFrom magrittr %>%

sim_infxn_2_death_timeseries <- function(casefat, I0, r, min_day, curr_day, m_od = 18.8, s_od = 0.45){

  # get expected number of infections up to today
  total_expected_inf <- casefat$pa * I0/r * exp(r*curr_day)

  # get expected deaths without censoring
  total_expected_death <- total_expected_inf * casefat$cfr

  # get actual number of deaths without censoring
  total_death_uncensored <- rpois(length(total_expected_death),
                                  lambda = total_expected_death)

  # expand out the death grid
  casefat_expand <- function(datrow, n){
    datrow <- datrow[rep(1, times = n), ]
    return(datrow)
  }
  case_line_list <- split(casefat[,c("age", "cfr")], 1:nrow(casefat))
  case_line_list <- mapply(casefat_expand, datrow = case_line_list,
                           n = total_death_uncensored, SIMPLIFY = F) %>%
    do.call("rbind.data.frame", .)

  case_line_list$t_infection <- unlist(lapply(total_death_uncensored, function(x) { log(runif(x))/r + curr_day }))
  case_line_list$t_death <- case_line_list$t_infection +
                                rgamma(n = nrow(case_line_list), shape = 1/s_od^2, scale = m_od*s_od^2)


  #..................
  # Tidy up so that we observe deaths on a daily time step
  #..................
  agelvls <- unique(as.character(casefat$age)) # protect against data.frame, string as factor = F
  death_line_list <- case_line_list %>%
    dplyr::mutate(age = factor(age, levels = agelvls)) %>%  # need this for later summarize
    dplyr::select(c("age", "t_death")) %>%
    dplyr::mutate(tod_fct = cut(t_death, breaks = c(1:curr_day), right = F)) %>%
    dplyr::filter(!is.na(tod_fct)) %>% # exceeds current day, so "future" death
    dplyr::group_by(tod_fct) %>%
    dplyr::mutate(obs_day = floor(min(t_death)),
                  obs_day = factor(obs_day, levels = 1:curr_day) # need this for later summarize
                  ) %>%
    dplyr::ungroup(tod_fct) %>%
    dplyr::group_by(obs_day, age, .drop = F) %>%
    dplyr::summarise(
      day_deaths = n()) %>%
    dplyr::ungroup(obs_day, age) %>%
    dplyr::group_by(age) %>%
    dplyr::mutate(obs_deaths = cumsum(day_deaths),
                  obs_day = as.numeric(obs_day) # factor levels match day, only had to factorize so .drop=F worked
                  ) %>%
    dplyr::ungroup(age) %>%  # for export
   dplyr::filter(obs_day >= min_day)

  # out
  return(death_line_list)
}
