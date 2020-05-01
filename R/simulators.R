#' @title Simulate Expected Deaths for a Given Day Assuming Exponential Growth of the Incidence Curve
#' @param I0 integer; number of infected individuals at the beginning of the epidemic for exponential growth
#' @param r double; the expontential growth rate of the epidemic
#' @param m_od double; the mean of the onset of infection to death (gamma distribution)
#' @param s_od double; the coefficient of variation of the onset of infection to death (gamma distribution)
#' @param casefat dataframe; column names age, cfr, and pa that correspond to age-band, case-fatality ratio, and the attack rate in given age band, respectively
#' @param min_day numeric; first day epidemic was observed
#' @param curr_day numeric; current day of epidemic, considered up to but not including this day
#' @param level character; Either "Cumulative" or "Time-Series" if expected deaths should be returned as a total or for each day, respectively
#' @param expgrowth logical; Simulate exponential growth or not
#' @param infections integer vector; The infections for each day
#' @importFrom magrittr %>%
#' @export

sim_infxn_2_death <- function(casefat, I0, r, m_od = 18.8, s_od = 0.45,
                              min_day = 1, curr_day, level, expgrowth, infections){

  #..............................................................
  # Assertions that are specific to this project
  #..............................................................
  assert_single_numeric(m_od)
  assert_single_numeric(s_od)
  assert_dataframe(casefat)
  assert_single_int(min_day)
  assert_single_int(curr_day)
  assert_single_string(level)
  assert_single_logical(expgrowth)
  assert_in(x = level, y = c("Time-Series", "Cumulative"))
  if (!expgrowth) {
    assert_vector(infections)
    assert_same_length(infections, min_day:(curr_day - 1))
  }
  if (expgrowth) {
    assert_single_int(I0)
    assert_single_numeric(r)
  }
  #..............................................................
  # Run
  #..............................................................
  # get  number of infections for each day
  if (expgrowth){
    It <- sapply(min_day:(curr_day-1), function(t){I0*exp(r*t)})
  } else {
    It <- infections
  }

  expected_inf.day <- rpois(length(It), lambda = It)

  # Split infxns by the Age Prop.
  expected_inf.age.day <- matrix(NA, nrow = length(casefat$pa), ncol = length(expected_inf.day))
  for (i in 1:length(expected_inf.day)) {
    expected_inf.age.day[,i] <- rmultinom(n = 1, size = expected_inf.day[i], prob = casefat$pa)
  }

  # draw deaths among Infected
  t_death <- matrix(NA, nrow = nrow(expected_inf.age.day), ncol = ncol(expected_inf.age.day))
  for (i in 1:nrow(expected_inf.age.day)) {
    for(j in 1:ncol(expected_inf.age.day)) {
      t_death[i,j] <- rbinom(n = 1, size = expected_inf.age.day[i,j], prob = casefat$cfr[i])
    }
  }

  # expand out the death grid
  t_death.df <- cbind.data.frame(age = casefat$age, t_death)
  colnames(t_death.df)[2:ncol(t_death.df)] <- min_day:(curr_day - 1)
  t_death.df <- t_death.df %>%
    tidyr::gather(., key = "day", value = "deathcount", 2:ncol(.))

  casefat_expand <- function(datrow){
    datrow <- datrow[rep(1, times = datrow$deathcount), ]
    datrow
    return(datrow)
  }
  death_line_list <- split(t_death.df, 1:nrow(t_death.df))
  death_line_list <- lapply(death_line_list, casefat_expand) %>%
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

  # out
  if (level == "Cumulative") {
    death_line_list <- death_line_list %>%
      dplyr::group_by(age) %>%
      dplyr::summarise(deaths = sum(day_deaths))
    ret <- list(
      death_line_list = death_line_list,
      infxns = expected_inf.day
    )
    return(ret)

  } else {
    death_line_list <- death_line_list %>%
      dplyr::ungroup(age)
    ret <- list(
      death_line_list = death_line_list,
      infxns = expected_inf.day
    )
    return(ret)
  }
}

