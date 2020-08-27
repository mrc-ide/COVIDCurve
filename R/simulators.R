#' Simulate Seroprevalence Study
#' @inheritParams Aggsim_infxn_2_death
#' @param serolin numeric vector; the expected number of infections from the infection curve from the poisson draw
#' @importFrom magrittr %>%
#' @noRd

sim_seroprev <- function(sero_line_list,
                         spec,
                         sens,
                         smplfrac,
                         sero_delay_rate,
                         demog,
                         fatalitydata,
                         curr_day) {


  # draw time from onset to seroconversion
  sero_line_list$otsc <- stats::rexp(n = nrow(sero_line_list), rate = 1/sero_delay_rate)
  # observed time of seroconversion
  sero_line_list$tosc <- as.numeric(sero_line_list$doi) + sero_line_list$otsc

  #..................
  # Tidy up so that we observe deaths on a daily time step
  #..................
  sero_line_list <- sero_line_list %>%
    dplyr::mutate(obs_day_of_serocon = cut(tosc, breaks = c(0, 1:curr_day),
                                           labels = 1:curr_day))


  #......................
  # observed seroprevalences taking into account sampling fractions
  #......................
  sero_line_list_sampl <- sero_line_list %>%
    dplyr::filter(!is.na(obs_day_of_serocon))   # drop "future" seroconversions
  keeprows <- as.logical(rbinom(n = nrow(sero_line_list_sampl), size = 1, prob = smplfrac))
  sero_line_list_sampl <- sero_line_list_sampl[keeprows, ]
  # get serotested after sampling fraction
  serotested <- tibble::tibble(strata = fatalitydata$strata,
                               n_tested = demog$popN * smplfrac)

  # get aggregate counts for model
  sero_strata_agg <- sero_line_list_sampl %>%
    dplyr::group_by(strata, obs_day_of_serocon, .drop = F) %>%
    dplyr::summarise(day_seros = dplyr::n()) %>%
    dplyr::left_join(., demog, by = "strata") %>%
    dplyr::left_join(., serotested, by = "strata") %>%
    dplyr::mutate(obsday = as.numeric(as.character(obs_day_of_serocon)), # protect against factor
                  true_serocount = cumsum(day_seros),
                  true_prev = true_serocount/n_tested,
                  obs_prev = sens*true_prev + (1-spec)*(1-true_prev)) %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-c("day_seros", "popN", "obs_day_of_serocon"))

  # out
  out <- list(sero_line_list = sero_line_list,
              sero_strata_agg = sero_strata_agg)
  return(out)
}


#' @title Simulate Aggregate Expected Deaths
#' @param infections integer vector; The infections for each day up to the current day.
#' @param m_od double; The mean of the onset of infection to death (gamma distribution).
#' @param s_od double; The coefficient of variation of the onset of infection to death (gamma distribution).
#' @param fatalitydata dataframe; The strata-specific fatalities to simulate given a probability of infection and a noise effect. The column names: strata, ifr, rho, and Ne correspond to (patient) strata, infection-fatality ratio, the probability of infection (i.e. a probalistic attack rate), and a noise effect, respectively.
#' @param demog dataframe; Strata-specific population (demographic) counts. The columns names strata and popN correspond to (patient) strata and the number of individuals within that strata. The demography strata must match the fatalitydata strata. Only considered if \code{simulate_seroprevalence = TRUE}
#' @param curr_day numeric; Current day of epidemic (considered up to but not including this day).
#' @param level character; Must either be "Time-Series" or "Cumulative", indicating whether daily death counts or cumulative deaths to the current day should be returned, respectively.
#' @param spec double; Specificity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param sens double; Sensitivity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param sero_delay_rate double; Rate of time from infection to seroconversion, assumed to be exponentially distributed (only considered if simulate_seroprevalence is set to TRUE)
#' @param smplfrac numeric; Sampling fraction for the observed seroprevalence study (assumed to be a simple random sample of all infected)
#' @importFrom magrittr %>%
#' @export

Aggsim_infxn_2_death <- function(fatalitydata, infections, m_od = 14.26, s_od = 0.79,
                                 curr_day,
                                 spec, sens, demog, sero_delay_rate, smplfrac = 1){

  #..................
  # Assertions that are specific to this project
  #..................
  assert_single_numeric(m_od)
  assert_single_numeric(s_od)
  assert_dataframe(fatalitydata)
  assert_eq(colnames(fatalitydata), c("strata", "IFR", "rho"))
  assert_single_int(curr_day)
  assert_vector(infections)
  assert_same_length(infections, 1:curr_day)
  assert_numeric(spec)
  assert_bounded(spec, left = 0, right = 1)
  assert_numeric(sens)
  assert_bounded(sens, left = 0, right = 1)
  assert_numeric(sero_delay_rate)
  assert_dataframe(demog)
  assert_eq(colnames(demog), c("strata", "popN"))
  assert_eq(demog$strata, fatalitydata$strata,
            message = "%s must equal %s -- check that your strata are in the same order")
  assert_bounded(smplfrac, left = 0, right = 1, inclusive_left = FALSE)

  #......................
  # misc helper function
  #......................
  df_expand <- function(datrow, col){
    assert_string(col)
    assert_length(nrow(datrow), 1, message = "df_expand internal function issue. Datrow must only be a single dataframe row")
    datrow <- datrow[rep(1, times = datrow[, col]), ]
    return(datrow)
  }

  #..................
  # Run Infxns and Deaths
  #..................
  # get  number of infections for each day
  expected_inf.day <- stats::rpois(length(infections), lambda = infections)

  # Split infxns by strata prop. and rho
  expected_inf.strt.day <- matrix(NA, nrow = length(fatalitydata$rho), ncol = length(expected_inf.day))
  Pinfxnsero <- fatalitydata$rho/sum(fatalitydata$rho)
  for (i in 1:length(expected_inf.day)) {
    expected_inf.strt.day[,i] <- stats::rmultinom(n = 1, size = expected_inf.day[i],
                                                  prob = Pinfxnsero)
  }

  # tidy infxn df
  # not needed here but for clarity
  colnames(expected_inf.strt.day) <- 1:curr_day
  infxn_df <- expected_inf.strt.day %>%
    cbind.data.frame(strata = fatalitydata$strata, .) %>%
    tidyr::pivot_longer(cols = -c("strata"), names_to = "doi", values_to = "infxncount")

  # expand out the infection linelist
  infxn_line_list <- split(infxn_df, 1:nrow(infxn_df))
  # tidy up and expand out death line list
  infxn_line_list <- lapply(infxn_line_list, df_expand, col = "infxncount") %>%
    do.call("rbind.data.frame", .) %>%
    dplyr::select(-c("infxncount")) %>%
    dplyr::mutate(id = 1:nrow(.)) %>%
    dplyr::select(c("id", "doi", "strata"))

  #..................
  # get seroprevalence
  #..................
  seroprev <- sim_seroprev(sero_line_list = infxn_line_list, spec = spec, sens = sens,
                           sero_delay_rate = sero_delay_rate, smplfrac = smplfrac,
                           demog = demog, fatalitydata = fatalitydata, curr_day = curr_day)

  #..................
  # get deaths
  #..................
  death_line_list <- dplyr::left_join(infxn_line_list, fatalitydata, by = "strata")
  # draw deaths among infected
  death_line_list <- death_line_list %>%
    dplyr::mutate(dies = purrr::map_int(IFR, function(x){rbinom(n = 1, size = 1, prob = x)})) %>%
    dplyr::filter(dies == 1) %>%
    dplyr::select(-c("dies"))

  # simulate onset of infection to death
  death_line_list$otd <- rgamma(nrow(death_line_list), shape = 1/s_od^2, scale = m_od*s_od^2)
  # simulate time of death -- note, we have a discrete day + a continuous time
  death_line_list$tod <- as.numeric(death_line_list$doi) + death_line_list$otd

  # Tidy up so that we observe deaths on a daily time step
  stratalvls <- unique(as.character(fatalitydata$strata)) # protect against data.frame, string as factor = F
  death_line_list <- death_line_list %>%
    dplyr::mutate(strata = factor(strata, levels = stratalvls)) %>%  # need this for later summarize
    dplyr::mutate(obs_day_of_death = cut(tod, breaks = c(0, 1:curr_day),
                                         labels = 1:curr_day))

  # tidy up for out
  death_strata_agg <- death_line_list %>%
    dplyr::filter(!is.na(obs_day_of_death)) %>% # drop "future" deaths
    dplyr::select(c("strata", "obs_day_of_death")) %>%
    dplyr::rename(
      obsday = obs_day_of_death,
      strata = strata) %>%
    dplyr::group_by(obsday, strata, .drop = F) %>%
    dplyr::summarise(
      day_deaths = dplyr::n()) %>%
    dplyr::ungroup(obsday, strata) %>%
    dplyr::mutate(obsday = as.numeric(as.character(obsday))) %>% # protect against factor
    dplyr::rename(deaths = day_deaths)
  # marginalize over for model
  death_agg <- death_strata_agg %>%
    dplyr::group_by(obsday) %>%
    dplyr::summarise(deaths = sum(deaths)) %>%
    dplyr::ungroup(.)

  # full linelist
  full_linelist <- dplyr::full_join(seroprev$sero_line_list, death_line_list, by = c("id", "doi", "strata")) %>%
    dplyr::select(-c("IFR", "rho")) # holdovers from merge with fatality data




  #..................
  # out
  #..................
  ret <- list(
    StrataAgg_TimeSeries_Death = death_strata_agg,
    Agg_TimeSeries_Death = death_agg,
    StrataAgg_Seroprev = seroprev$sero_strata_agg,
    full_linelist = full_linelist)

  return(ret)
}


