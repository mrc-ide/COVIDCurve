#' @title Likelihood for Curve Aware
r_timeseries_likelihood <- function(params, param_i, data, misc) {

  # from lancet id paper
  r <- 0.14
  m_od <- 18.8
  s_od <- 0.45

  # free params
  I0 <- params["I0"]
  r1 <- params["r1"]
  ma2 <- params["ma2"]

  # get age-specific mortality rates
  ma <- c(r1*ma2, ma2)

  get_exp_deaths_by_day <- function(day){
    # integrate for expected incidence mapped onto onset-death time lag
    integrand <- function(t, gr = r){
      I0 * exp(gr*t) * ( pgamma(day-t, shape = 1/s_od^2, scale = m_od*s_od^2) -
                         pgamma(day-t-1, shape = 1/s_od^2, scale = m_od*s_od^2) )
    }
    integral <- integrate(integrand, lower = -Inf, upper = day)

    # expected deaths
    exp.deaths <- misc$pa * ma * integral$value
    return(exp.deaths)
  }
  # loop through deaths each day
  exp_deaths.day <- lapply(misc$min_day:misc$curr_day, get_exp_deaths_by_day)

  # poisson LL
  get_pois <- function(obs_death, exp_death){
    ret <- sum(dpois(x = obs_death, lambda = exp_death, log = T))
    #ret <- ret + length(scalars) * log(ma9) # account for reparameterization -- LEGACY FULL MODEL
    return(ret)
  }

  loglik <- mapply(get_pois, obs_death = data$obs_deaths, exp_death = exp_deaths.day)
  loglik <- sum(loglik)
  return(loglik)
}


