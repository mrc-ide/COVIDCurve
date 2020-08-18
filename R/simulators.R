#' Simulate Seroprevalence Study
#' @inheritParams Aggsim_infxn_2_death
#' @param seroinfxns numeric vector; the expected number of infections from the infection curve from the poisson draw
#' @importFrom magrittr %>%
#' @noRd

sim_seroprev <- function(seroinfxns,
                         spec,
                         sens,
                         sero_delay_rate,
                         demog,
                         fatalitydata,
                         min_day,
                         curr_day) {

  # seroinfxns matrix of stratified infections (rows) by day -- minday:currday -- of infection onset (columns)
  sero.df <- cbind.data.frame(fatalitydata$Strata, seroinfxns)
  colnames(sero.df) <- c("Strata", min_day:curr_day)
  sero.df <- sero.df %>%
    tidyr::pivot_longer(cols = -c("Strata"), names_to = "ObsDay", values_to = "infxncount") %>%
    dplyr::mutate(ObsDay = as.numeric(ObsDay)) # coercing character to numeric

  # expand this out to line-list
  df_expand <- function(datrow){
    datrow <- datrow[rep(1, times = datrow$infxncount), ]
    return(datrow)
  }
  sero_line_list <- split(sero.df, 1:nrow(sero.df))
  sero_line_list <- lapply(sero_line_list, df_expand) %>%
    dplyr::bind_rows(.) %>%
    dplyr::select(-c("infxncount")) %>%
    tibble::as_tibble(.)

  # draw time to seroconversion
  draw_tosc <- function(day, sero_delay_rate){
    as.numeric(day) + stats::rexp(n = 1, rate = 1/sero_delay_rate)
  }
  sero_line_list$tosc <- sapply(sero_line_list$ObsDay, draw_tosc, sero_delay_rate = sero_delay_rate)

  #..................
  # Tidy up so that we observe deaths on a daily time step
  #..................
  sero_line_list <- sero_line_list %>%
    dplyr::mutate(event_obs_day = cut(tosc, breaks = c((min_day-1), min_day:curr_day),
                                      labels = min_day:curr_day)) %>%
    dplyr::filter(!is.na(event_obs_day))   # drop "future" deaths

  sero_agg <- sero_line_list %>%
    dplyr::group_by(Strata, event_obs_day, .drop = F) %>%
    dplyr::summarise(day_seros = dplyr::n()) %>%
    dplyr::left_join(., demog, by = "Strata") %>%
    dplyr::mutate(event_obs_day = as.numeric(as.character(event_obs_day)), # protect against factor and min_day > 1
                  TrueSeroCount = cumsum(day_seros),
                  TruePrev = TrueSeroCount/popN,
                  ObsPrev = sens*TruePrev + (1-spec)*(1-TruePrev)) %>%
    dplyr::ungroup(.)

  # out
  out <- list(sero_line_list = sero_line_list,
              sero_agg = sero_agg)
  return(out)
}


#' @title Simulate Aggregate Expected Deaths
#' @param infections integer vector; The infections for each day up to the current day.
#' @param m_od double; The mean of the onset of infection to death (gamma distribution).
#' @param s_od double; The coefficient of variation of the onset of infection to death (gamma distribution).
#' @param fatalitydata dataframe; The strata-specific fatalities to simulate given a probability of infection and a noise effect. The column names: strata, ifr, rho, and Ne correspond to (patient) strata, infection-fatality ratio, the probability of infection (i.e. a probalistic attack rate), and a noise effect, respectively.
#' @param demog dataframe; Strata-specific population (demographic) counts. The columns names strata and popN correspond to (patient) strata and the number of individuals within that strata. The demography strata must match the fatalitydata strata. Only considered if \code{simulate_seroprevalence = TRUE}
#' @param min_day numeric; First day epidemic was observed.
#' @param curr_day numeric; Current day of epidemic (considered up to but not including this day).
#' @param level character; Must either be "Time-Series" or "Cumulative", indicating whether daily death counts or cumulative deaths to the current day should be returned, respectively.
#' @param simulate_seroprevalence logical; Whether or not seroprevalence data should also be simulated
#' @param spec double; Specificity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param sens double; Sensitivity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param sero_delay_rate double; Rate of time from infection to seroconversion, assumed to be exponentially distributed (only considered if simulate_seroprevalence is set to TRUE)
#' @importFrom magrittr %>%
#' @export

Aggsim_infxn_2_death <- function(fatalitydata, infections, m_od = 14.2, s_od = 0.79,
                                 min_day = 1, curr_day,
                                 simulate_seroprevalence = TRUE,
                                 spec, sens, demog, sero_delay_rate){

  #..................
  # Assertions that are specific to this project
  #..................
  assert_single_numeric(m_od)
  assert_single_numeric(s_od)
  assert_dataframe(fatalitydata)
  assert_eq(colnames(fatalitydata), c("Strata", "IFR", "Rho", "Ne"))
  assert_single_int(min_day)
  assert_single_int(curr_day)
  assert_vector(infections)
  assert_same_length(infections, min_day:curr_day)
  assert_logical(simulate_seroprevalence)
  if (simulate_seroprevalence) {
    assert_numeric(spec)
    assert_bounded(spec, left = 0, right = 1)
    assert_numeric(sens)
    assert_bounded(sens, left = 0, right = 1)
    assert_numeric(sero_delay_rate)
    assert_dataframe(demog)
    assert_eq(colnames(demog), c("Strata", "popN"))
    assert_eq(demog$Strata, fatalitydata$Strata,
              message = "%s must equal %s -- check that your strata are in the same order")
  }

  #..................
  # Run Infxns and Deaths
  #..................
  # get  number of infections for each day
  expected_inf.day <- stats::rpois(length(infections), lambda = infections)

  # Split infxns by strata prop. and noise effect (i.e. a "random" effect but calling noise effect as we are not using a traditional MLM)
  expected_inf.strt.day <- matrix(NA, nrow = length(fatalitydata$Rho), ncol = length(expected_inf.day))
  Pinfxnsero <- (fatalitydata$Rho * fatalitydata$Ne)/sum((fatalitydata$Rho * fatalitydata$Ne))
  for (i in 1:length(expected_inf.day)) {
    expected_inf.strt.day[,i] <- stats::rmultinom(n = 1, size = expected_inf.day[i],
                                                  prob = Pinfxnsero)
  }

  # Aside: tidy this for out
  tidy_expected_inf.strt.day <- expected_inf.strt.day %>%
    t(.) %>%
    tibble::as_tibble(., .name_repair = c("minimal")) %>%
    magrittr::set_colnames(paste0("infxns_", as.character(fatalitydata$Strata))) %>%
    dplyr::mutate(ObsDay = min_day:curr_day) %>%
    dplyr::select(c("ObsDay", dplyr::everything(.)))

  # draw deaths among Infected
  t_death <- matrix(NA, nrow = nrow(expected_inf.strt.day), ncol = ncol(expected_inf.strt.day))
  for (i in 1:nrow(expected_inf.strt.day)) {
    for(j in 1:ncol(expected_inf.strt.day)) {
      t_death[i,j] <- stats::rbinom(n = 1, size = expected_inf.strt.day[i,j], prob = fatalitydata$IFR[i])
    }
  }

  # expand out the death grid
  t_death.df <- cbind.data.frame(strata = fatalitydata$Strata, t_death)
  colnames(t_death.df)[2:ncol(t_death.df)] <- min_day:(curr_day)


  t_death.df <- t_death.df %>%
    tidyr::pivot_longer(cols = -c("strata"), names_to = "day", values_to = "deathcount")

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
  stratalvls <- unique(as.character(fatalitydata$Strata)) # protect against data.frame, string as factor = F
  death_line_list <- death_line_list %>%
    dplyr::mutate(strata = factor(strata, levels = stratalvls)) %>%  # need this for later summarize
    dplyr::select(c("strata", "tod")) %>%
    dplyr::mutate(obs_day = cut(tod, breaks = c((min_day-1), min_day:curr_day),
                                labels = min_day:curr_day)) %>%
    dplyr::filter(!is.na(obs_day)) %>%  # drop "future" deaths
    dplyr::rename(
      ObsDay = obs_day,
      Strata = strata)

  death_agg <- death_line_list %>%
    dplyr::group_by(ObsDay, Strata, .drop = F) %>%
    dplyr::summarise(
      day_deaths = dplyr::n()) %>%
    dplyr::ungroup(ObsDay, Strata) %>%
    dplyr::mutate(ObsDay = as.numeric(as.character(ObsDay))) %>% # protect against factor and min_day > 1
    dplyr::rename(Deaths = day_deaths)


  #..................
  # run seroprev
  #..................
  if (simulate_seroprevalence) {
    seroprev <- sim_seroprev(seroinfxns = expected_inf.strt.day, spec = spec, sens = sens, sero_delay_rate = sero_delay_rate,
                             demog = demog, fatalitydata = fatalitydata, min_day = min_day, curr_day = curr_day)
  }

  #..................
  # out
  #..................

  if (simulate_seroprevalence) {

    ret <- list(
      AggDeath = death_agg,
      AggSeroPrev = seroprev$sero_agg,
      AggInfxns = tidy_expected_inf.strt.day,
      LineListDeath = death_line_list,
      LineListSero = seroprev$sero_line_list)
    return(ret)

  } else {

    ret <- list(
      AggDeath = death_agg,
      AggInfxns = tidy_expected_inf.strt.day,
      LineListDeath = death_line_list)

    return(ret)
  }
}


