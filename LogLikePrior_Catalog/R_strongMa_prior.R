#' @title Prior for Curve Aware
r_strongMa_prior <- function(params, param_i, misc) {

  I0 <- params["I0"]
  ma1 <- params["ma1"]
  ma2 <- params["ma2"]


  # get prior
  ret <- dunif(I0, min = 1, max = 10, log = TRUE) +
    dbeta(ma1, shape1 = 100, shape2 = 900, log = TRUE) +
    dbeta(ma2, shape1 = 500, shape2 = 500, log = TRUE)

  return(ret)
}

