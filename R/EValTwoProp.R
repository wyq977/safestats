likelihoodBernoulli <- function(n, x, theta) {
  # Use log and exp to avoid underflow/overflow
  return(x * log(theta) + (n - x) * log(1 - theta))
}


# TODO: update priors with none restrictions in k x 2 settings
# TODO: change Beta dist. parameters
# TODO: default two-sided tests only
#
#' Calculate e-value for 2 x 2 contingency table
#' @export
EValTwoProp <- function(ya, yb, alternative = c("twoSided", "greater", "less"),
                        na = 1, nb = 1, esType = c("logOddsRatio", "difference", "none"),
                        esMin = NULL, alpha = 0.05, prior = NULL) {

  if (is.null(prior)) {
    # print("Using the default parameters (alpha = beta = 0.18) for Beta distribution...")
    betaA1 <- betaA2 <- betaB1 <- betaB2 <- 0.18
  } else {
    # unpack the prior values
    betaA1 <- prior[["betaA1"]]
    betaA2 <- prior[["betaA2"]]
    betaB1 <- prior[["betaB1"]]
    betaB2 <- prior[["betaB2"]]
  }

  eVal <- 1

  # case 1: No restriction on prior --------------------------------------------
  if (esType == "none") {
    if (!is.null(esMin)) {
      stop("There should be no minimum difference in this case!")
    } else if (length(ya) != length(yb)) {
      stop("The numer of samples in group A and group B must be equal!")
    }

    # init vector
    totalSuccessA <- cumsum(ya)
    totalSuccessB <- cumsum(yb)
    groupSizeVecA <- seq_along(totalSuccessA) * na
    groupSizeVecB <- seq_along(totalSuccessB) * nb
    totalFailA <- groupSizeVecA - totalSuccessA
    totalFailB <- groupSizeVecB - totalSuccessB

    # # naive approach
    # thetaA <- thetaB <- theta0 <- 0.5
    #
    # for (i in seq_along(ya)) {
    #   newE <- safestats:::calculateETwoProportions(
    #     na1 = ya[i],
    #     na = na,
    #     nb1 = yb[i],
    #     nb = nb,
    #     thetaA = thetaA,
    #     thetaB = thetaB,
    #     theta0 = theta0
    #   )
    #
    #   eVal <- eVal * newE
    #
    #   thetaA <- safestats:::bernoulliMLTwoProportions(totalSuccessA[i], totalFailA[i], betaA1, betaA2)
    #
    #   thetaB <- safestats:::bernoulliMLTwoProportions(totalSuccessB[i], totalFailB[i], betaB1, betaB2)
    #
    #   theta0 <- (na * thetaA + nb * thetaB) / (na + nb)
    # }

    # vectorized approach
    thetaA <- thetaB <- theta0 <- rep(0.5, length(ya))

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

  } else if (esType == "difference" || esType == "logOddsRatio") {
    # case 2/3: restrict difference on (0, 1)^2 para. space --------------------
    # they all used Bayensian updating in each timesteps between blocks
    # the only difference are the calculation of thetaB and grid
    if (!is.numeric(esMin)) {
      stop("esMin must be numeric for difference")
    }
    gridSize <- 1e3
    K <- 1 / gridSize

    rhoGrid <- seq(K, 1 - K, length.out = gridSize)

    ## difference in probability grid ------------------------------------------
    if (esType == "difference") {
      # prepare prob. grid = rho x (1 - esMin)
      thetaAgrid <- rhoGrid * (1 - esMin)
      thetaBgrid <- thetaAgrid + esMin
    } else if (esType == "logOddsRatio") {
      thetaAgrid <- rhoGrid
      thetaBgrid <- sapply(thetaAgrid, calculateThetaBFromThetaAAndLOR, lOR = delta)
    }

    # prior density for thetaA at each grid points with normalization
    # rho = thetaA / (1 - esMin) ~ Beta(shape1, shape2)
    thetaADensity <- stats::dbeta(rhoGrid, shape1 = betaA1, shape2 = betaA2)
    thetaADensity <- thetaADensity / sum(thetaADensity)

    # TODO: Is the initial theta 1 / 2?
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

    ## difference in calculating thetaB ----------------------------------------
    if (esType == "difference") {
      thetaB <- thetaA + esMin
    } else if (esType == "logOddsRatio") {
      thetaB <- safestats:::calculateThetaBFromThetaAAndLOR(thetaA = thetaA, lOR = esMin)
    }

    theta0 <- (na * thetaA + nb * thetaB) / (na + nb)

    # calculate eValVec
    eValVec <- safestats:::calculateETwoProportions(
      na1 = ya, na = na, nb1 = yb, nb = nb,
      thetaA = thetaA, thetaB = thetaB, theta0 = theta0
    )

    eVal <- prod(eValVec)
  }

  return(eVal)
}
