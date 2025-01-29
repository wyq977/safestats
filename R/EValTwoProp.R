likelihoodBernoulli <- function(n, x, theta) {
  # Use log and exp to avoid underflow/overflow
  return(x * log(theta) + (n - x) * log(1 - theta))
}


# TODO: update priors with none restrictions in k x 2 settings
# TODO: change Beta dist. parameters
# TODO: default two-sided tests only
# TODO: na, nb as vectors of length(ya)
#
#' Calculate e-value for 2 x 2 contingency table
#' @export
EValTwoProp <- function(ya, yb, alternative = c("twoSided", "greater", "less"),
                        na = 1, nb = 1, esType = c("logOddsRatio", "difference", "none"),
                        esMin = NULL, alpha = 0.05, prior = NULL) {
  betaA1 <- betaA2 <- betaB1 <- betaB2 <- 0.18
  eVal <- 1

  if (esType == "none") {
    # init vector
    totalSuccessA <- cumsum(ya)
    totalSuccessB <- cumsum(yb)
    groupSizeVecA <- seq_along(totalSuccessA) * na
    groupSizeVecB <- seq_along(totalSuccessB) * nb
    totalFailA <- groupSizeVecA - totalSuccessA
    totalFailB <- groupSizeVecB - totalSuccessB

    thetaA <- thetaB <- theta0 <- rep(0.5, length(ya))
    # case 1: No restriction on prior, just calculate

    # vectorize operation except 1st
    theta <- safestats:::updateETwoProportions(
      totalSuccessA[-length(totalSuccessA)], totalFailA[-length(totalFailA)],
      totalSuccessB[-length(totalSuccessB)], totalFailB[-length(totalFailB)],
      na, nb,
      betaA1, betaA2,
      betaB1, betaB2
    )

    # assign values to all elements except 1st
    thetaA[-1] <- theta$thetaA
    thetaB[-1] <- theta$thetaB
    theta0[-1] <- theta$theta0

    eValVec <- safestats:::calculateETwoProportions(
      na1 = ya, na = na, nb1 = yb, nb = nb,
      thetaA = thetaA, thetaB = thetaB, theta0 = theta0
    )

    eVal <- prod(eValVec)
  } else if (esType == "logOddsRatio") {
    # case 2: log-odds ratio difference
    if (!is.numeric(esMin)) {
      stop("esMin must be numeric for difference")
    }
  } else if (esType == "difference") {
    # case 3: absolute difference
    if (!is.numeric(esMin)) {
      stop("esMin must be numeric for difference")
    }
    # FIXME: plotting debugs
    gridSize <- 1e3
    K <- 1 / gridSize

    rhoGrid <- seq(K / (1 - esMin), 1 - K, length.out = gridSize)

    # prepare prob. grid = rho x (1 - esMin)
    thetaAgrid <- rhoGrid * (1 - esMin)
    thetaBgrid <- thetaAgrid + esMin

    # prior density for thetaA at each grid points with normalization
    # rho = thetaA / (1 - esMin) ~ Beta(shape1, shape2)
    thetaADensity <- stats::dbeta(rhoGrid, shape1 = betaA1, shape2 = betaA2)
    thetaADensity <- thetaADensity / sum(thetaADensity)

    thetaB <- thetaA <- theta0 <- rep(0.5, length(ya))

    for (i in 1:length(ya)) {
      # Compute log-likelihoods for thetaA and thetaB
      logLikelihoodA <- likelihoodBernoulli(na, ya[i], thetaAgrid)
      logLikelihoodB <- likelihoodBernoulli(nb, yb[i], thetaBgrid)

      # Compute posterior in log-space for stability
      logPosterior <- logLikelihoodA + logLikelihoodB + log(thetaADensity)

      # Convert back to probability space
      posteriorDensity <- exp(logPosterior)

      # Normalize to update thetaADensity
      thetaADensity <- posteriorDensity / sum(posteriorDensity)

      # Compute expectation of thetaA (posterior mean)
      thetaA[i] <- thetaAgrid %*% thetaADensity
    }

    thetaB <- thetaA + esMin
    theta0 <- (na * thetaA + nb * thetaB)/(na + nb)

    # calculate eValVec
    eValVec <- safestats:::calculateETwoProportions(
      na1 = ya, na = na, nb1 = yb, nb = nb,
      thetaA = thetaA, thetaB = thetaB, theta0 = theta0
    )

    eVal <- prod(eValVec)

  }

  return(eVal)
}
