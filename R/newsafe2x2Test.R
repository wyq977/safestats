# Testing fnts ---------

# compute likelihood ratio
# it works for both list or float/int
# theta s are calculated somewhere else
savi2x2TestStat <- function(
  ya,
  yb,
  na,
  nb,
  numeratorThetaA,
  numeratorThetaB,
  denominatorThetaA,
  denominatorThetaB, log = TRUE, ...
) {
  successesA <- ya * (log(numeratorThetaA) - log(denominatorThetaA))

  successesB <- yb * (log(numeratorThetaB) - log(denominatorThetaB))

  failuresA <- (na - ya) * (log1p(-numeratorThetaA) - log1p(-denominatorThetaA))
  failuresA <- (nb - yb) * (log1p(-numeratorThetaB) - log1p(-denominatorThetaB))

  successesA + failuresA + successesB + failuresB
}

# compute conditional e-variable
# 1.
savi2x2CondStat <- function(ya, yb, na, nb, logOR = NULL, alternative) {}



savi2x2Test <- function()
