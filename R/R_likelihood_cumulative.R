#' @title Likelihood for Curve Aware
r_tod_log_like_cumulative <- function(params, param_i, data, misc) {

  # from lancet id paper
  m_od <- 18.8
  r <- 0.14
  s_od <- 0.45

  # free params
  curr_day <- misc$curr_day
  ma2 <- params[2]
  r1 <- params[1]
  #I0 <- params[3]
  I0 <- 2

  # LEGACY FOR FULL MODEL
  # ma9 <- params[1]
  # r1 <- params[2]
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


  # integrate for expected incidence mapped onto onset-death time lag
  integrand <- function(t, gr = r){ return(
    I0 * exp(gr*t) * pgamma(curr_day - t, shape = 1/s_od^2, scale = m_od*s_od^2)) }
  integral <- integrate(integrand, lower = -Inf, upper = curr_day)

  # total exp deaths
  exp.deaths <- misc$pa * ma * integral$value

  # poisson LL
  ret <- sum(dpois(x = data$obs_deaths, lambda = exp.deaths, log = T)) + log(ma2)
  #ret <- ret + length(scalars) * log(ma2) # account for reparameterization
  return(ret)

}

#' @title Prior for Curve Aware
r_tod_log_prior_cumulative <- function(params, param_i, misc) {
  # free params
  ma2 <- params[2]
  r1 <- params[1]
  #I0 <- params[3]

  # LEGACY FOR FULL PROGRAM
  # ma9 <- params[1]
  # r1 <- params[2]
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
  ma1 <- r1 * ma2
  ma <- c(ma1, ma2)
  mapa <- mean( ma * misc$pa )

  # # flat prior
  # ret <- dunif(ma9, min = 0, max = 1, log = TRUE) +
  #   sum( sapply(scalars, function(x){dlnorm(x, meanlog = 0, sdlog = 5, log = TRUE)}) )
  # + dunif(I0, min = 0, max = 10, log = TRUE)

  ret <- dunif(ma2, min = 0, max = 1, log = TRUE) +
    dlnorm(r1, meanlog = 0, sdlog = 5, log = TRUE)

  return(ret)
}
