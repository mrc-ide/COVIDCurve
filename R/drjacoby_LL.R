#' @title Likelihood for Curve Aware
r_tod_log_like <- function(params, param_i, data, misc) {
  # assume r is fixed
  r <- 0.14
  curr_day <- misc$curr_day

  # alpha and beta for time from infection to time to death
  # from lancet id paper
  m_od <- 18.8
  s_od <- 0.45
  # free params
  I0 <- params[1]
  ma <- params[2]
  #ma1 <- params[2]
  #ma2 <- params[3]
  # ma3 <- params[4]
  # ma4 <- params[5]
  # ma5 <- params[6]
  # ma6 <- params[7]
  # ma7 <- params[8]
  # ma8 <- params[9]
  # ma9 <- params[10]
  # ma <- c(ma1, ma2, ma3, ma4, ma5, ma6, ma7, ma8, ma9)

  #..................
  # expected deaths
  #..................
  integrand <- function(t, gr = r){ return(
    I0 * exp(gr*t) * pgamma(curr_day - t, shape = 1/s_od^2, scale = m_od*s_od^2)) }
  integral <- integrate(integrand, lower = -Inf, curr_day)

  # total exp deaths
  exp.deaths <- misc$pa * ma * integral$value

  #..................
  # poisson
  #..................
  ret <- sum(dpois(x = data$obs_deaths, lambda = exp.deaths, log = T))
  return(ret)

}

#' @title Prior for Curve Aware
r_tod_log_prior <- function(params, param_i, misc) {
  I0 <- params[1]
  ma <- params[2]
  #I0 <- params[1]
  #ma1 <- params[2]
  #ma2 <- params[3]
  # ma3 <- params[4]
  # ma4 <- params[5]
  # ma5 <- params[6]
  # ma6 <- params[7]
  # ma7 <- params[8]
  # ma8 <- params[9]
  # ma9 <- params[10]
  # ma <- c(ma1, ma2, ma3, ma4, ma5, ma6, ma7, ma8, ma9)
  #ma <- c(ma1, ma2)

  # flat prior
  #ret <- dunif(I0, min = 0, max = 10, log = TRUE) +
   # sum( sapply(ma, function(x){dunif(x, min = 0, max = 1, log = TRUE)}) )
  ret <-  dunif(I0, min = 0, max = 10, log = TRUE) +
    dunif(ma, min = 0, max = 1, log = TRUE)

  return(ret)
}



