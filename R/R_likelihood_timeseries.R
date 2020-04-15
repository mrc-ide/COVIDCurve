#' @title Likelihood for Curve Aware
r_tod_log_like_timeseries <- function(params, param_i, data, misc) {

  # from lancet id paper
  r <- 0.14
  m_od <- 18.8
  s_od <- 0.45

  # free params
  ma2 <- params[2]
  r1 <- params[1]
  #I0 <- params[3]
  I0 <- 2

  # LEGACY FOR FULL MODEL
  #ma9 <- params[1]
  #r1 <- params[2]
  # r2 <- params[3]
  # r3 <- params[4]
  # r4 <- params[5]
  # r5 <- params[6]
  # r6 <- params[7]
  # r7 <- params[8]
  # r8 <- params[9]
  # I0 <- params[10]
  # scalars <- c(r1, r2, r3, r4, r5, r6, r7, r8)
  # ma <- sapply(scalars, function(x){x * ma9})
  # ma <- c(ma, ma9)
  ma <- c(r1*ma2, ma2)

  get_exp_deaths_by_day <- function(day){
    # integrate for expected incidence mapped onto onset-death time lag
    integrand <- function(t, gr = r){ return(
      I0 * exp(gr*t) * (pgamma(day-t+1, shape = 1/s_od^2, scale = m_od*s_od^2) -
                          pgamma(day-t, shape = 1/s_od^2, scale = m_od*s_od^2)))}
    integral <- integrate(integrand, lower = -Inf, upper = day)

    # expected deaths
    exp.deaths <- misc$pa * ma * integral$value
    return(exp.deaths)
  }
  # loop through deaths each day
  exp_deaths.day <- lapply(misc$min_day:misc$curr_day, get_exp_deaths_by_day)

  # poisson LL
  get_pois <- function(obs_death, exp_death){
    ret <- sum(dpois(x = obs_death, lambda = exp_death, log = T)) + log(ma2)
    #ret <- ret + length(scalars) * log(ma9) # account for reparameterization -- LEGACY FULL MODEL
    return(ret)
  }

  loglik <- mapply(get_pois, obs_death = data$obs_deaths, exp_death = exp_deaths.day)
  loglik <- sum(loglik)
  return(loglik)
}

#' @title Prior for Curve Aware
r_tod_log_prior_timeseries <- function(params, param_i, misc) {

  ma2 <- params[2]
  r1 <- params[1]
  #I0 <- params[3]


  # LEGACY FOR FULL MODEL
  #ma9 <- params[1]
  #r1 <- params[2]
  # r2 <- params[3]
  # r3 <- params[4]
  # r4 <- params[5]
  # r5 <- params[6]
  # r6 <- params[7]
  # r7 <- params[8]
  # r8 <- params[9]
  # I0 <- params[10]

  # set up strong prior on "average" M_a*P_a
  ma1 <- r1 * ma2
  ma <- c(ma1, ma2)
  mapa <- mean( ma * misc$pa )

  # priors
  ret <- dunif(ma2, min = 0, max = 1, log = TRUE) +
    dlnorm(r1, meanlog = 0, sdlog = 5, log = TRUE)

#  + dlnorm(I0, meanlog = 0.69, sdlog = 0.05, log = TRUE)

  return(ret)
}





