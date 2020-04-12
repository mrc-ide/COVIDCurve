#' @title Likelihood for Curve Aware
r_tod_log_like <- function(params, param_i, data, misc) {

  # assume r is fixed
  r <- 0.14
  # alpha and beta for time from infection to time to death
  # from lancet id paper
  m_od <- 18.8
  s_od <- 0.45
  I0 <- 2
  # free params
  ma <- params[1]
  #r1 <- params[2]
  #I0 <- params[2]
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
  #ma <- sapply(scalars, function(x){x * ma9})

  get_exp_deaths_by_day <- function(bigT){

    # expected deaths
    integrand <- function(t, gr = r){ return(
      I0 * exp(gr*t) * pgamma(bigT - t, shape = 1/s_od^2, scale = m_od*s_od^2)) }
    integral <- integrate(integrand, lower = -Inf, bigT)

    # total exp deaths
    exp.deaths <- misc$pa * ma * integral$value
    # out
    return(exp.deaths)
  }

  exp_deaths.day <- lapply(misc$min_day:misc$curr_day, get_exp_deaths_by_day)


  #..................
  # poisson
  #..................
  get_pois <- function(obs_death, exp_death){
    ret <- sum(dpois(x = obs_death, lambda = exp_death, log = T))
    #ret <- ret + length(scalars) * log(ma9) # account for reparameterization
    #ret <- ret + log(ma2) # account for reparameterization
    return(ret)
  }

  loglik <- mapply(get_pois, obs_death = data$obs_deaths, exp_death = exp_deaths.day)
  loglik <- sum(loglik)
  return(loglik)

}

#' @title Prior for Curve Aware
r_tod_log_prior <- function(params, param_i, misc) {

  ma <- params[1]
  #r1 <- params[2]
  #I0 <- params[2]
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
  # flat prior
  dayprior <- sum(
                  sapply(misc$min_day:misc$curr_day, function(dummy){
                   # ret <- dunif(ma9, min = 0, max = 1, log = TRUE) +
                   #   sum( sapply(scalars, function(x){dlnorm(x, meanlog = 0, sdlog = 5, log = TRUE)}) )
                    #+ dunif(I0, min = 0, max = 10, log = TRUE)
                    ret <- dunif(ma, min = 0, max = 1, log = TRUE)
                      #dlnorm(r1, meanlog = 0, sdlog = 5, log = TRUE) +
                      #+ dlnorm(I0, meanlog = 0.7, sd = 0.05, log = TRUE)

                    return(ret)
                  })
              )
  return(dayprior)
}



