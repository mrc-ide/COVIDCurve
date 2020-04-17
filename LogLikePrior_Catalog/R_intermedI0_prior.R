#' @title Prior for Curve Aware
r_intermedI0_prior <- function(params, param_i, misc) {

  I0 <- params["I0"]
  r1 <- params["r1"]
  ma2 <- params["ma2"]

  # scale ma1
  ma1 <- r1 * ma2

  # get prior
  ret <- dlnorm(I0, meanlog = 0.69, sdlog = 0.5, log = TRUE) +
    dlnorm(r1, meanlog = 0, sdlog = 5, log = TRUE) +
    dunif(ma2, min = 0, max = 1, log = TRUE) +
    log(ma2) # account for reparameterization

  return(ret)
}

