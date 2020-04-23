#' @title Likelihood for Curve Aware
r_cumulative_likelihood <- function(params, param_i, data, misc) {

  # from lancet id paper
  m_od <- 18.8
  r <- 0.14
  s_od <- 0.45

  # store misc item for integral
  curr_day <- misc$curr_day

  # free params
  I0 <- params["I0"]
  ma1 <- params["ma1"]
  ma2 <- params["ma2"]

  # get age specific mortality
  ma <- c(ma1, ma2)

  # integrate for expected incidence mapped onto onset-death time lag
  integrand <- function(t, gr = r){ return(
    I0 * exp(gr*t) * pgamma(curr_day - t, shape = 1/s_od^2, scale = m_od*s_od^2)) }
  integral <- integrate(integrand, lower = -Inf, upper = curr_day)

  # total exp deaths
  exp.deaths <- misc$pa * ma * integral$value

  # poisson LL
  ret <- sum(dpois(x = data$obs_deaths, lambda = exp.deaths, log = T))
  return(ret)

}
