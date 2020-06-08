#' Simulate Seroprevalence Study
#' @inheritParams Aggsim_infxn_2_death
#' @param seroinfxns numeric vector; the expected number of infections from the infection curve from the poisson draw
#' @noRd

sim_seroprev <- function(seroinfxns,
                         spec,
                         sens,
                         underreport,
                         sero_delay_rate,
                         popN,
                         fatalitydata,
                         min_day,
                         curr_day) {
  # store prevalence which is function of infections before delay to seroconversion
  dz_pos <- cumsum(seroinfxns)
  dz_neg <- popN - dz_pos
  fps <- dz_neg - (dz_neg * spec)

  # recast infections for seroprev delay
  sero.df <- data.frame(day =  min_day:curr_day,
                        infxncount = seroinfxns)

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
  sero_line_list$tosc <- sapply(sero_line_list$day, draw_tosc, sero_delay_rate = 1/sero_delay_rate)

  # Tidy up so that we observe deaths on a daily time step
  sero_line_list.agg <- sero_line_list %>%
    dplyr::mutate(event_obs_day = cut(tosc, breaks = c((min_day-1), min_day:curr_day),
                                      labels = min_day:curr_day)) %>%
    dplyr::filter(!is.na(event_obs_day)) %>%   # drop "future" deaths
    dplyr::group_by(event_obs_day, .drop = F) %>%
    dplyr::summarise(day_seros = dplyr::n()) %>%
    dplyr::mutate(day_seros = cumsum(day_seros))

  # these are the proportion of infxns we will observe -- and different truths
  sero_line_list.agg$day_seros_fn <- rpois(n = nrow(sero_line_list.agg), lambda = (sero_line_list.agg$day_seros * sens))
  sero_line_list.agg$day_seros_fp <- sero_line_list.agg$day_seros + fps
  sero_line_list.agg$day_seros_fn_fp <- sero_line_list.agg$day_seros_fn + fps

  #..................
  # out
  #..................
  sero_line_list.agg %>%
    dplyr::mutate(event_obs_day = as.numeric(as.character(event_obs_day))) %>% # protect against factor and min_day > 1
    dplyr::mutate(
      TrueSeroRate = day_seros/popN,
      SeroRateFP = day_seros_fp/popN,
      SeroRateFN = day_seros_fn/popN,
      SeroRateFNFP = day_seros_fn_fp/popN) %>%
    dplyr::rename(
      ObsDay = event_obs_day,
      TrueSeroCount = day_seros,
      SeroCountFP = day_seros_fp,
      SeroCountFN = day_seros_fn,
      SeroCountFNFP = day_seros_fn_fp)


}


#' @title Simulate Aggregate Expected Deaths
#' @param infections integer vector; The infections for each day up to the current day.
#' @param m_od double; The mean of the onset of infection to death (gamma distribution).
#' @param s_od double; The coefficient of variation of the onset of infection to death (gamma distribution).
#' @param m_or double; The mean of the onset of infection to recovery (gamma distribution).
#' @param s_or double; The coefficient of variation of the onset of infection to recovery (gamma distribution).
#' @param fatalitydata dataframe; The column names: strata, ifr, and pa correspond to (patient) strata, infection-fatality ratio, and the attack rate, respectively.
#' @param min_day numeric; First day epidemic was observed.
#' @param curr_day numeric; Current day of epidemic (considered up to but not including this day).
#' @param level character; Must either be "Time-Series" or "Cumulative", indicating whether daily death counts or cumulative deaths to the current day should be returned, respectively.
#' @param simulate_seroprevalence logical; Whether or not seroprevalence data should also be simulated
#' @param spec double; Specificity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param sens double; Sensitivity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param underreport double; Fraction of deaths that are assumed to be accounted for (i.e. not missed)
#' @param popN double; Population Size (only considered if simulate_seroprevalence is set to TRUE)
#' @param sero_delay_rate double; Rate of time from infection to seroconversion, assumed to be exponentially distributed (only considered if simulate_seroprevalence is set to TRUE)
#' @importFrom magrittr %>%
#' @export

Aggsim_infxn_2_death <- function(fatalitydata, infections, m_od = 18.8, s_od = 0.45,
                                 min_day = 1, curr_day, level,
                                 simulate_seroprevalence = TRUE,
                                 spec, sens, underreport = 1, popN, sero_delay_rate){

  #..................
  # Assertions that are specific to this project
  #..................
  assert_single_numeric(m_od)
  assert_single_numeric(s_od)
  assert_dataframe(fatalitydata)
  assert_in(colnames(fatalitydata), c("strata", "ifr", "pa"))
  assert_single_int(min_day)
  assert_single_int(curr_day)
  assert_single_string(level)
  assert_in(x = level, y = c("Time-Series", "Cumulative"))
  assert_vector(infections)
  assert_same_length(infections, min_day:curr_day)
  assert_logical(simulate_seroprevalence)
  if (simulate_seroprevalence) {
    assert_numeric(spec)
    assert_bounded(spec, left = 0, right = 1)
    assert_numeric(sens)
    assert_bounded(sens, left = 0, right = 1)
    assert_numeric(sero_delay_rate)
    assert_pos_int(popN)
  }

  if (sum(fatalitydata$pa) != 1) {
    warning("Prob. of infection (pa) does not sum 1. Standardizing now.")
    fatalitydata$pa <- fatalitydata$pa/sum(fatalitydata$pa)
  }

  #..................
  # Run Infxns and Deaths
  #..................
  # get  number of infections for each day
  expected_inf.day <- rpois(length(infections), lambda = infections)

  # Split infxns by the Age Prop.
  expected_inf.age.day <- matrix(NA, nrow = length(fatalitydata$pa), ncol = length(expected_inf.day))
  for (i in 1:length(expected_inf.day)) {
    expected_inf.age.day[,i] <- rmultinom(n = 1, size = expected_inf.day[i], prob = fatalitydata$pa)
  }

  # draw deaths among Infected
  t_death <- matrix(NA, nrow = nrow(expected_inf.age.day), ncol = ncol(expected_inf.age.day))
  for (i in 1:nrow(expected_inf.age.day)) {
    for(j in 1:ncol(expected_inf.age.day)) {
      t_death[i,j] <- rbinom(n = 1, size = expected_inf.age.day[i,j], prob = fatalitydata$ifr[i])
    }
  }
  # account for underreporting
  t_death.ur <- matrix(rpois(n = length(t_death), lambda = t_death * underreport),
                       nrow = nrow(expected_inf.age.day), ncol = ncol(expected_inf.age.day))
  # expand out the death grid
  t_death.df <- cbind.data.frame(age = fatalitydata$strata, t_death.ur)
  colnames(t_death.df)[2:ncol(t_death.df)] <- min_day:(curr_day)
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
  agelvls <- unique(as.character(fatalitydata$strata)) # protect against data.frame, string as factor = F
  death_line_list <- death_line_list %>%
    dplyr::mutate(age = factor(age, levels = agelvls)) %>%  # need this for later summarize
    dplyr::select(c("age", "tod")) %>%
    dplyr::mutate(obs_day = cut(tod, breaks = c((min_day-1), min_day:curr_day),
                                labels = min_day:curr_day)) %>%
    dplyr::filter(!is.na(obs_day)) %>% # drop "future" deaths
    dplyr::group_by(obs_day, age, .drop = F) %>%
    dplyr::summarise(
      day_deaths = dplyr::n()) %>%
    dplyr::ungroup(obs_day, age)

  #..................
  # run seroprev
  #..................
  if (simulate_seroprevalence) {
    seroprev <- sim_seroprev(seroinfxns = expected_inf.day, spec = spec, sens = sens, sero_delay_rate = sero_delay_rate,
                             popN = popN, fatalitydata = fatalitydata, min_day = min_day, curr_day = curr_day)
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
        Strata = age,
        Deaths = deaths)
  } else {
    death_line_list <- death_line_list %>%
      dplyr::mutate(obs_day = as.numeric(as.character(obs_day))) %>% # protect against factor and min_day > 1
      dplyr::rename(
        ObsDay = obs_day,
        Strata = age,
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


