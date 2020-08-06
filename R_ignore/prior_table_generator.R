# https://stackoverflow.com/questions/12088080/how-to-convert-integer-number-into-binary-vector
# https://stackoverflow.com/questions/54059917/generate-all-length-n-permutations-of-true-false

#......................
# get combination of T/F vectors
#......................
n <- 5
out <- list()
for(i in 0:(2^n-1)){
  out <- append(out, list(binaryLogic::as.binary(i, n = n)))
}
priortab <- out %>%
  do.call("rbind.data.frame", .) %>%
  magrittr::set_colnames(c("rIFR", "rKnots", "rInfxn", "rSeros", "rNes"))

#......................
# reparam lines for each param
#......................
priortab <- priortab[32:1, ]

