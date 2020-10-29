#' Simulate Seroprevalence Study
#' @inheritParams Agesim_infxn_2_death
#' @param sero_line_list dataframe; small line list generated from simulator
#' @importFrom magrittr %>%
#' @noRd

sim_seroprev <- function(sero_line_list,
                         spec,
                         sens,
                         smplfrac,
                         sero_delay_rate,
                         simulate_seroreversion,
                         sero_rev_shape,
                         sero_rev_scale,
                         demog,
                         fatalitydata,
                         curr_day) {


  # draw time from onset to seroconversion
  sero_line_list$otsc <- stats::rexp(n = nrow(sero_line_list), rate = 1/sero_delay_rate)
  # observed time of seroconversion
  sero_line_list$tosc <- (as.numeric(sero_line_list$doi)-1) + sero_line_list$otsc

  if (simulate_seroreversion) {
    # draw time from seroconversion to seroreversion
    sero_line_list$otsr <- stats::rweibull(n = nrow(sero_line_list),
                                           shape = sero_rev_shape, scale = sero_rev_scale)
    # observed time of seroreversion
    sero_line_list$tosr <- as.numeric(sero_line_list$tosc) + sero_line_list$otsr
  } else {
    sero_line_list$otsr <- NA
    sero_line_list$tosr <- Inf
  }

  #..................
  # Tidy up so that we observe deaths on a daily time step
  #..................
  sero_line_list <- sero_line_list %>%
    dplyr::mutate(ObsDaySeroCon = cut(tosc, breaks = c(0, 1:curr_day),
                                      labels = 1:curr_day),
                  ObsDaySeroRev = cut(tosr, breaks = c(0, 1:curr_day),
                                      labels = 1:curr_day),
                  ObsDaySeroRev = ifelse(is.na(ObsDaySeroRev), Inf, ObsDaySeroRev)) # if seroreversion is missing, the subject doesn't revert within the study period


  #......................
  # observed seroprevalences taking into account sampling fractions
  #......................
  sero_line_list_sampl <- sero_line_list %>%
    dplyr::filter(!is.na(ObsDaySeroCon))   # drop "future" seroconversions
  keeprows <- as.logical(stats::rbinom(n = nrow(sero_line_list_sampl), size = 1, prob = smplfrac))
  sero_line_list_sampl <- sero_line_list_sampl[keeprows, ]
  # get serotested after sampling fraction
  serotested <- tibble::tibble(Strata = fatalitydata$Strata,
                               testedN = demog$popN * smplfrac)

  # get count of seroreversions by day -- these are lost infections
  sero_revertsdf <- sero_line_list_sampl %>%
    dplyr::mutate(ObsDaySeroRev = factor(ObsDaySeroRev, levels = c(1:curr_day))) %>%
    dplyr::group_by(Strata, ObsDaySeroRev, .drop = F) %>%
    dplyr::summarise(daily_seroreverts = dplyr::n()) %>%
    dplyr::mutate(ObsDay = as.numeric(as.character(ObsDaySeroRev)), # protect against factor
                  daily_cum_seroreverts = cumsum(daily_seroreverts)) %>%
    dplyr::ungroup(.) %>%
    dplyr::filter(ObsDay <= curr_day) %>%
    dplyr::select(c("ObsDay", "Strata", "daily_cum_seroreverts"))


  # get aggregate counts for model
  sero_strata_agg <- sero_line_list_sampl %>%
    dplyr::group_by(Strata, ObsDaySeroCon, .drop = F) %>%
    dplyr::summarise(daily_seroconverts = dplyr::n()) %>%
    dplyr::mutate(ObsDay = as.numeric(as.character(ObsDaySeroCon)), # protect against factor
                  daily_cum_seroconverts = cumsum(daily_seroconverts)) %>%
    dplyr::left_join(., sero_revertsdf, by = c("ObsDay", "Strata")) %>%
    dplyr::left_join(., demog, by = "Strata") %>%
    dplyr::left_join(., serotested, by = "Strata") %>%
    dplyr::mutate(
      TrueSeroCount = daily_cum_seroconverts,
      RevertSeroCount = daily_cum_seroconverts - daily_cum_seroreverts,
      TruePrev = TrueSeroCount/testedN,
      ObsPrev = RevertSeroCount/testedN, # observed prev corrected for seroreverts
      ObsPrev = sens*ObsPrev + (1-spec)*(1-ObsPrev)) %>% # observed prev corrected for spec/sens
    dplyr::ungroup(.) %>%
    dplyr::select(-c("daily_seroconverts", "daily_cum_seroconverts", "daily_cum_seroreverts", "RevertSeroCount",
                     "popN", "ObsDaySeroCon"))

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
#' @param spec double; Specificity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param sens double; Sensitivity of the Seroprevalence Study (only considered if simulate_seroprevalence is set to TRUE)
#' @param sero_delay_rate double; Rate of time from infection to seroconversion, assumed to be exponentially distributed (only considered if simulate_seroprevalence is set to TRUE)
#' @param simulate_seroreversion logical; Whether seroreversion (due to waning of antibodies) should be simulated or not
#' @param sero_rev_shape double; The shape parameter of the Weibull seroreversion distribution
#' @param sero_rev_scale double; The scale parameter of the Weibull seroreversion distribution
#' @param smplfrac numeric; Sampling fraction for the observed seroprevalence study (assumed to be a simple random sample of all infected)
#' @param return_linelist logical; Whether or not the linelist that was used to create the observed marginal data should be returned. N.B. the linelist can be quite large/burdensome for memory depending on population size and number of days considered.
#' @importFrom magrittr %>%
#' @export

Agesim_infxn_2_death <- function(fatalitydata, infections, m_od = 14.26, s_od = 0.79,
                                 curr_day,
                                 spec, sens, demog, sero_delay_rate,
                                 simulate_seroreversion, sero_rev_shape = NULL, sero_rev_scale = NULL,
                                 smplfrac = 1, return_linelist = FALSE){

  #..................
  # Assertions that are specific to this project
  #..................
  assert_single_numeric(m_od)
  assert_single_numeric(s_od)
  assert_dataframe(fatalitydata)
  assert_eq(colnames(fatalitydata), c("Strata", "IFR", "Rho"))
  assert_single_int(curr_day)
  assert_vector(infections)
  assert_same_length(infections, 1:curr_day)
  assert_numeric(spec)
  assert_bounded(spec, left = 0, right = 1)
  assert_numeric(sens)
  assert_bounded(sens, left = 0, right = 1)
  assert_numeric(sero_delay_rate)
  assert_dataframe(demog)
  assert_eq(colnames(demog), c("Strata", "popN"))
  assert_eq(demog$Strata, fatalitydata$Strata,
            message = "%s must equal %s -- check that your strata are in the same order")
  assert_bounded(smplfrac, left = 0, right = 1, inclusive_left = FALSE)
  assert_logical(simulate_seroreversion)
  if (simulate_seroreversion) {
    assert_numeric(sero_rev_shape)
    assert_numeric(sero_rev_scale)
  }


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
  # get  number of total infections in the entire population for each day
  expected_inf.day <- stats::rpois(length(infections), lambda = infections)

  # Given the total number of infections in the population, split by population density
  # while accounting for some degree of differential attack rates through Rho
  # P(Infxn) = P(Infxn|age) * P(age)
  expected_inf.strt.day <- matrix(NA, nrow = nrow(fatalitydata), ncol = length(expected_inf.day))
  Pinfxnsero <- (demog$popN * fatalitydata$Rho) / sum((demog$popN * fatalitydata$Rho))
  for (i in 1:length(expected_inf.day)) {
    expected_inf.strt.day[,i] <- stats::rmultinom(n = 1, size = expected_inf.day[i],
                                                  prob = Pinfxnsero)
  }

  # tidy infxn df
  # not needed here but for clarity
  colnames(expected_inf.strt.day) <- 1:curr_day
  infxn_df <- expected_inf.strt.day %>%
    cbind.data.frame(Strata = fatalitydata$Strata, .) %>%
    tidyr::pivot_longer(cols = -c("Strata"), names_to = "doi", values_to = "infxncount")

  # expand out the infection linelist
  infxn_line_list <- split(infxn_df, 1:nrow(infxn_df))
  # tidy up and expand out death line list
  infxn_line_list <- lapply(infxn_line_list, df_expand, col = "infxncount") %>%
    do.call("rbind.data.frame", .) %>%
    dplyr::select(-c("infxncount")) %>%
    dplyr::mutate(id = 1:nrow(.)) %>%
    dplyr::select(c("id", "doi", "Strata"))

  #..................
  # get seroprevalence
  #..................
  seroprev <- sim_seroprev(sero_line_list = infxn_line_list, spec = spec, sens = sens,
                           sero_delay_rate = sero_delay_rate,
                           simulate_seroreversion = simulate_seroreversion,
                           sero_rev_shape = sero_rev_shape, sero_rev_scale = sero_rev_scale,
                           smplfrac = smplfrac,
                           demog = demog, fatalitydata = fatalitydata, curr_day = curr_day)


  #..................
  # get deaths
  #..................
  death_line_list <- dplyr::left_join(infxn_line_list, fatalitydata, by = "Strata")
  # draw deaths among infected
  death_line_list <- death_line_list %>%
    dplyr::mutate(dies = purrr::map_int(IFR, function(x){stats::rbinom(n = 1, size = 1, prob = x)})) %>%
    dplyr::filter(dies == 1) %>%
    dplyr::select(-c("dies"))

  # simulate onset of infection to death
  death_line_list$otd <- stats::rgamma(nrow(death_line_list), shape = 1/s_od^2, scale = m_od*s_od^2)
  # simulate time of death -- note, we have a discrete day + a continuous time -- but allow it to happen on the day it was observed
  death_line_list$tod <- (as.numeric(death_line_list$doi)-1) + death_line_list$otd

  # Tidy up so that we observe deaths on a daily time step
  stratalvls <- unique(as.character(fatalitydata$Strata)) # protect against data.frame, string as factor = F
  death_line_list <- death_line_list %>%
    dplyr::mutate(Strata = factor(Strata, levels = stratalvls)) %>%  # need this for later summarize
    dplyr::mutate(ObsDayDeath = cut(tod, breaks = c(0, 1:curr_day),
                                    labels = 1:curr_day))

  # tidy up for out
  death_strata_agg <- death_line_list %>%
    dplyr::filter(!is.na(ObsDayDeath)) %>% # drop "future" deaths
    dplyr::select(c("Strata", "ObsDayDeath")) %>%
    dplyr::rename( ObsDay = ObsDayDeath ) %>%
    dplyr::group_by(ObsDay, Strata, .drop = F) %>%
    dplyr::summarise( day_deaths = dplyr::n() ) %>%
    dplyr::ungroup(ObsDay, Strata) %>%
    dplyr::mutate(ObsDay = as.numeric(as.character(ObsDay))) %>% # protect against factor
    dplyr::rename(Deaths = day_deaths)
  # marginalize over for model
  death_agg <- death_strata_agg %>%
    dplyr::group_by(ObsDay) %>%
    dplyr::summarise(Deaths = sum(Deaths)) %>%
    dplyr::ungroup(.)

  # full linelist
  full_linelist <- dplyr::full_join(seroprev$sero_line_list, death_line_list, by = c("id", "doi", "Strata")) %>%
    dplyr::select(-c("IFR", "Rho")) # holdovers from merge with fatality data




  #..................
  # out
  #..................
  if (return_linelist) {
    ret <- list(
      StrataAgg_TimeSeries_Death = death_strata_agg,
      Agg_TimeSeries_Death = death_agg,
      StrataAgg_Seroprev = seroprev$sero_strata_agg,
      full_linelist = full_linelist)
  } else {
    ret <- list(
      StrataAgg_TimeSeries_Death = death_strata_agg,
      Agg_TimeSeries_Death = death_agg,
      StrataAgg_Seroprev = seroprev$sero_strata_agg)
  }


  return(ret)
}


