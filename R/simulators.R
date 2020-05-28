#' @title Simulate Seroprevalence Study
#' @inheritParams LineListsim_infxn_2_death
#' @param sero_spec double; Sensitivity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param sero_sens double; Specificity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param popN double; Population Size (only considered if simulate_seroprevalence is set to TRUE)
#' @param sero_delay_rate double; Rate of time from infection to seroconversion, assumed to be exponentially distributed (only considered if simulate_seroprevalence is set to TRUE)
#' internal function, not exported

sim_seroprev <- function(infxns,
                         sero_spec,
                         sero_sens,
                         sero_delay_rate,
                         popN,
                         casefat,
                         min_day,
                         curr_day) {
  # store prevalence which is function of infections before delay to seroconversion
  dz_pos <- cumsum(colSums(infxns))
  dz_neg <- popN - dz_pos
  fps <- dz_neg - (dz_neg * sero_spec)

  # recast infections for seroprev delay
  sero.df <- as.data.frame(infxns)
  colnames(sero.df)[1:ncol(sero.df)] <- min_day:(curr_day - 1)
  sero.df <- sero.df %>%
    tidyr::gather(., key = "day", value = "infxncount", 1:ncol(.)) %>%
    dplyr::filter(infxncount != 0) # can't observe if not infxn

  df_expand <- function(datrow){
    datrow <- datrow[rep(1, times = datrow$infxncount), ]
    return(datrow)
  }

  sero_line_list <- split(sero.df, 1:nrow(sero.df))
  sero_line_list <- lapply(sero_line_list, df_expand) %>%
    dplyr::bind_rows(.) %>%
    tibble::as_tibble(.) %>%
    dplyr::select(-c("infxncount"))

  # draw time to seroconversion
  draw_tosc <- function(day, sero_delay_rate){
    as.numeric(day) + rexp(n = 1, rate = sero_delay_rate)
  }
  sero_line_list$tosc <- sapply(sero_line_list$day, draw_tosc, sero_delay_rate = sero_delay_rate)

  # Tidy up so that we observe deaths on a daily time step

  sero_line_list.agg <- sero_line_list %>%
    dplyr::mutate(event_obs_day = cut(tosc, breaks = c((min_day-1), min_day:(curr_day-1)),
                                      labels = min_day:(curr_day-1))) %>%
    dplyr::filter(!is.na(event_obs_day)) %>%   # drop "future" deaths
    dplyr::group_by(event_obs_day) %>%
    dplyr::mutate(
      pInfxns = sero_sens,
      TestposDzpos = purrr::map_dbl(pInfxns, function(x) rbinom(1, 1, x)) ) %>%
    dplyr::group_by(event_obs_day, .drop = F) %>%
    dplyr::summarise(nInfxns = dplyr::n(),
                     day_seros = sum(TestposDzpos)) %>%
    dplyr::mutate(day_seros = cumsum(day_seros),
                  day_seros_fp = day_seros + fps)

  #..................
  # out
  #..................
  sero_line_list.agg %>%
    dplyr::mutate(event_obs_day = as.numeric(as.character(event_obs_day))) %>% # protect against factor and min_day > 1
    dplyr::rename(
      ObsDay = event_obs_day,
      SeroCount = day_seros,
      SeroCountFP = day_seros_fp) %>%
    dplyr::mutate(
      SeroRate = SeroCount/popN,
      SeroRateFP = SeroCountFP/popN)

}


#' @title Simulate Aggregate Expected Deaths
#' @inheritParams LineListsim_infxn_2_death
#' @inheritParams sim_seroprev
#' @importFrom magrittr %>%
#' @export

Aggsim_infxn_2_death <- function(casefat, infections, m_od = 18.8, s_od = 0.45,
                                 min_day = 1, curr_day, level,
                                 simulate_seroprevalence = TRUE,
                                 sero_spec, sero_sens, popN, sero_delay_rate){

  #..................
  # Assertions that are specific to this project
  #..................
  assert_single_numeric(m_od)
  assert_single_numeric(s_od)
  assert_dataframe(casefat)
  assert_single_int(min_day)
  assert_single_int(curr_day)
  assert_single_string(level)
  assert_in(x = level, y = c("Time-Series", "Cumulative"))
  assert_vector(infections)
  assert_same_length(infections, min_day:(curr_day - 1))
  assert_logical(simulate_seroprevalence)
  if (simulate_seroprevalence) {
    assert_numeric(sero_spec)
    assert_bounded(sero_spec, left = 0, right = 1)
    assert_numeric(sero_sens)
    assert_bounded(sero_sens, left = 0, right = 1)
    assert_numeric(sero_delay_rate)
    assert_pos_int(popN)
  }
  #..................
  # Run Infxns and Deaths
  #..................
  # get  number of infections for each day
  expected_inf.day <- rpois(length(infections), lambda = infections)

  # Split infxns by the Age Prop.
  expected_inf.age.day <- matrix(NA, nrow = length(casefat$pa), ncol = length(expected_inf.day))
  for (i in 1:length(expected_inf.day)) {
    expected_inf.age.day[,i] <- rmultinom(n = 1, size = expected_inf.day[i], prob = casefat$pa)
  }

  # draw deaths among Infected
  t_death <- matrix(NA, nrow = nrow(expected_inf.age.day), ncol = ncol(expected_inf.age.day))
  for (i in 1:nrow(expected_inf.age.day)) {
    for(j in 1:ncol(expected_inf.age.day)) {
      t_death[i,j] <- rbinom(n = 1, size = expected_inf.age.day[i,j], prob = casefat$ifr[i])
    }
  }

  # expand out the death grid
  t_death.df <- cbind.data.frame(age = casefat$age, t_death)
  colnames(t_death.df)[2:ncol(t_death.df)] <- min_day:(curr_day - 1)
  t_death.df <- t_death.df %>%
    tidyr::gather(., key = "day", value = "deathcount", 2:ncol(.))

  df_expand <- function(datrow){
    datrow <- datrow[rep(1, times = datrow$deathcount), ]
    return(datrow)
  }
  death_line_list <- split(t_death.df, 1:nrow(t_death.df))
  death_line_list <- lapply(death_line_list, df_expand) %>%
    do.call("rbind.data.frame", .) %>%
    dplyr::select(-c("deathcount"))
  # note, we have a discrete day + a continuous time
  death_line_list$tod <- as.numeric(death_line_list$day) + rgamma(nrow(death_line_list), shape = 1/s_od^2, scale = m_od*s_od^2)

  # Tidy up so that we observe deaths on a daily time step
  agelvls <- unique(as.character(casefat$age)) # protect against data.frame, string as factor = F
  death_line_list <- death_line_list %>%
    dplyr::mutate(age = factor(age, levels = agelvls)) %>%  # need this for later summarize
    dplyr::select(c("age", "tod")) %>%
    dplyr::mutate(obs_day = cut(tod, breaks = c((min_day-1), min_day:(curr_day-1)),
                                labels = min_day:(curr_day-1))) %>%
    dplyr::filter(!is.na(obs_day)) %>% # drop "future" deaths
    dplyr::group_by(obs_day, age, .drop = F) %>%
    dplyr::summarise(
      day_deaths = dplyr::n()) %>%
    dplyr::ungroup(obs_day, age)

  #..................
  # run seroprev
  #..................
  if (simulate_seroprevalence) {
    seroprev <- sim_seroprev(infxns = expected_inf.age.day, sero_spec = sero_spec, sero_sens = sero_sens, sero_delay_rate = sero_delay_rate,
                             popN = popN, casefat = casefat, min_day = min_day, curr_day = curr_day)
  }

  #..................
  # out
  #..................
  if (level == "Cumulative") {
    death_line_list <- death_line_list %>%
      dplyr::mutate(obs_day = as.numeric(as.character(obs_day)),
                    ObsDay = max(obs_day)) %>% # protect against factor and min_day > 1
      dplyr::group_by(ObsDay, age) %>%
      dplyr::summarise(deaths = sum(day_deaths)) %>%
      dplyr::rename(
        AgeGroup = age,
        Deaths = deaths)
  } else {
    death_line_list <- death_line_list %>%
      dplyr::ungroup(age) %>%
      dplyr::mutate(obs_day = as.numeric(as.character(obs_day))) %>% # protect against factor and min_day > 1
      dplyr::rename(
        ObsDay = obs_day,
        AgeGroup = age,
        Deaths = day_deaths)
  }

  if (simulate_seroprevalence) {
    ret <- list(
      AggDat = death_line_list,
      seroprev = seroprev)
    return(ret)
  } else {
    return(death_line_list)
  }
}






#............................................................
# Line List
#...........................................................

#' @title Simulate Line-List Outcomes
#' @param infections integer vector; The infections for each day up to the current day (must be supplied if \param{expgrowth} is false).
#' @param m_od double; The mean of the onset of infection to death (gamma distribution).
#' @param s_od double; The coefficient of variation of the onset of infection to death (gamma distribution).
#' @param m_or double; The mean of the onset of infection to recovery (gamma distribution).
#' @param s_or double; The coefficient of variation of the onset of infection to recovery (gamma distribution).
#' @param casefat dataframe; The column names: age, ifr, and pa correspond to age-band, infection-fatality ratio, and the attack rate, respectively.
#' @param hospprob numeric; Probability of hospitalization given infection.
#' @param min_day numeric; First day epidemic was observed.
#' @param curr_day numeric; Current day of epidemic (considered up to but not including this day).
#' @param simulate_seroprevalence logical; Boolean indicating whether or not the seroprevalence should also be simulated
#' @inheritParams sim_seroprev
#' @importFrom magrittr %>%
#' @export

LineListsim_infxn_2_death <- function(casefat, infections,
                                      m_od = 18.8, s_od = 0.45,
                                      m_or = 24.7, s_or = 0.45,
                                      hospprob,
                                      min_day = 1, curr_day,
                                      simulate_seroprevalence = FALSE) {

  #..................
  # Assertions that are specific to this project
  #..................
  assert_single_numeric(m_od)
  assert_single_numeric(s_od)
  assert_single_numeric(m_or)
  assert_single_numeric(s_or)
  assert_single_numeric(hospprob)
  assert_dataframe(casefat)
  assert_single_int(min_day)
  assert_single_int(curr_day)
  assert_vector(infections)
  assert_same_length(infections, min_day:(curr_day - 1))
  assert_logical(simulate_seroprevalence)

  #..................
  # Run Line List Infection to Deaths
  #..................
  # get  number of infections for each day
  expected_inf.day <- rpois(length(infections), lambda = infections)

  # Split infxns by the Age Prop.
  expected_inf.age.day <- matrix(NA, nrow = length(casefat$pa), ncol = length(expected_inf.day))
  for (i in 1:length(expected_inf.day)) {
    expected_inf.age.day[,i] <- rmultinom(n = 1, size = expected_inf.day[i], prob = casefat$pa)
  }

  # recast infections
  infxn.df <- cbind.data.frame(age = casefat$age, expected_inf.age.day)
  colnames(infxn.df)[2:ncol(infxn.df)] <- min_day:(curr_day - 1)
  infxn.df <- infxn.df %>%
    tidyr::gather(., key = "day", value = "infxncount", 2:ncol(.)) %>%
    dplyr::filter(infxncount != 0) # can't observe if not infxn

  df_expand <- function(datrow){
    datrow <- datrow[rep(1, times = datrow$infxncount), ]
    return(datrow)
  }

  infxn_line_list <- split(infxn.df, 1:nrow(infxn.df))
  infxn_line_list <- lapply(infxn_line_list, df_expand) %>%
    dplyr::bind_rows(.) %>%
    tibble::as_tibble(.) %>%
    dplyr::select(-c("infxncount"))

  # draw time of outcome among Infected
  draw_out <- function(infxn_line, casefat, m_od = m_od, s_od = s_od, m_or = m_or, s_or = s_or){
    if (rbinom(n = 1, size = 1, prob = casefat$ifr[infxn_line$age == casefat$age])) {
      ret <- list(
        outcome = "Death",
        todr = as.numeric(infxn_line$day) + rgamma(1, shape = 1/s_od^2, scale = m_od*s_od^2)
      )
    } else {
      ret <- list(
        outcome = "Recovery",
        todr = as.numeric(infxn_line$day) + rgamma(1, shape = 1/s_or^2, scale = m_or*s_or^2)
      )
    }
    return(ret)
  }
  infxn_line_list$out <- lapply(split(infxn_line_list, 1:nrow(infxn_line_list)),
                                draw_out, casefat = casefat,
                                m_od = m_od, s_od = s_od,
                                m_or = m_or, s_or = s_or)
  infxn_line_list$outcome <- purrr::map_chr(infxn_line_list$out, "outcome")
  infxn_line_list$todr <- purrr::map_dbl(infxn_line_list$out, "todr")

  # draw hospitilization
  infxn_line_list$hosp <- rbinom(n = nrow(infxn_line_list), size = 1, prob = hospprob)

  # Tidy up so that we observe deaths on a daily time step
  agelvls <- unique(as.character(casefat$age)) # protect against data.frame, string as factor = F
  infxn_line_list <- infxn_line_list %>%
    dplyr::mutate(age = factor(age, levels = agelvls)) %>%  # need this for later summarize
    dplyr::mutate(event_obs_day = cut(todr, breaks = c((min_day-1), min_day:(curr_day-1)),
                                      labels = min_day:(curr_day-1))) %>%
    dplyr::filter(!is.na(event_obs_day)) %>%   # drop "future" deaths
    dplyr::select(c("age", "hosp", "day", "outcome", "event_obs_day")) %>%
    dplyr::mutate(day = as.numeric(as.character(day)),
                  event_obs_day = as.numeric(as.character(event_obs_day))) %>%  # protect against factor but out numeric
    dplyr::rename(AgeGroup = age,
                  OnsetDay = day,
                  Outcome = outcome,
                  EventDay = event_obs_day)

  #..................
  # run seroprev
  #..................
  #TODO -- update me

  #..................
  # out
  #..................
  if (simulate_seroprevalence) {
    ret <- list(
      linelist = infxn_line_list,
      seroprev = seroprev)
    return(ret)
  } else {
    return(infxn_line_list)
  }
}



