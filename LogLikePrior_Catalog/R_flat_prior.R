#' @title Prior for Curve Aware
r_flat_prior <- function(params, param_i, misc) {

  # free params
  I0 <- params["I0"]
  r1 <- params["r1"]
  ma2 <- params["ma2"]

  # get prior
  ret <- dunif(I0, min = 1, max = 10, log = TRUE) +
         dlnorm(r1, meanlog = 0, sdlog = 5, log = TRUE) +
         dunif(ma2, min = 0, max = 1, log = TRUE)


  return(ret)
}





