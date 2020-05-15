#' @title Simulate Aggregate Expected Deaths for a Given Day Assuming Exponential Growth of the Incidence Curve
#' @param I0 integer; number of infected individuals at the beginning of the epidemic for exponential growth
#' @param r double; the expontential growth rate of the epidemic
#' @param m_od double; the mean of the onset of infection to death (gamma distribution)
#' @param s_od double; the coefficient of variation of the onset of infection to death (gamma distribution)
#' @param m_or double; the mean of the onset of infection to recovery (gamma distribution)
#' @param s_or double; the coefficient of variation of the onset of infection to recovery (gamma distribution)
#' @param casefat dataframe; column names age, ifr, and pa that correspond to age-band, infection-fatality ratio, and the attack rate in given age band, respectively
#' @param hospprob numeric; probability of hospitalization given infection
#' @param min_day numeric; first day epidemic was observed
#' @param curr_day numeric; current day of epidemic, considered up to but not including this day
#' @param expgrowth logical; Simulate exponential growth or not
#' @param infections integer vector; The infections for each day
#' @importFrom magrittr %>%
#' @export

LineListsim_infxn_2_death <- function(casefat, I0, r,
                                      m_od = 18.8, s_od = 0.45,
                                      m_or = 24.7, s_or = 0.45,
                                      hospprob,
                                      min_day = 1, curr_day, expgrowth, infections){

  #..............................................................
  # Assertions that are specific to this project
  #..............................................................
  assert_single_numeric(m_od)
  assert_single_numeric(s_od)
  assert_single_numeric(m_or)
  assert_single_numeric(s_or)
  assert_single_numeric(hospprob)
  assert_dataframe(casefat)
  assert_single_int(min_day)
  assert_single_int(curr_day)
  assert_single_logical(expgrowth)
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
  infxn_line_list$outcome <- unlist( purrr::map(infxn_line_list$out, "outcome") )
  infxn_line_list$todr <- unlist( purrr::map(infxn_line_list$out, "todr") )

  # draw hospitilization
  infxn_line_list$hosp <- rbinom(n = nrow(infxn_line_list), size = 1, prob = hospprob)


  # Tidy up so that we observe deaths on a daily time step
  agelvls <- unique(as.character(casefat$age)) # protect against data.frame, string as factor = F
  # out
  infxn_line_list %>%
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


}

