#' Calculate the e-value for a Safe Experiment to Test Two Proportions in Stream Data
#'
#' @param ya positive observations/events per data block in group a: a numeric with integer values
#' between (and including) 0 and \code{na}, the number of observations in group a per block.
#' @param yb positive observations/events per data block in group b: a numeric with integer values
#' between (and including) 0 and \code{nb}, the number of observations in group b per block.
#' @param alternative comparision in the alternative hypothesis
#' @param na number of observations in group a per data block
#' @param nb number of observations in group b per data block.
#' The number of data blocks arrived in two streams (group a/b) must be the same
#' while the size of data block can vary (\code{na} != \code{nb}).
#' @param esType a character string specifying an optional restriction on the alternative hypothesis; must be one of "none" (default),
#' "difference" (difference group mean b minus group b) or "logOddsRatio" (the log odds ratio between group means b and a).
#' @param esMin a priori minimal relevant divergence between group means b and a, either a numeric between -1 and 1 for
#' no alternative restriction or a restriction on difference, or a real for a restriction on the log odds ratio.
#' @param alpha numeric in (0, 1)
#' @param prior a list of shape parameters for priori Beta distribution for group a, group b.
#'
#' @return e-value between 0 and \code{Inf}
#' @export
EValTwoProp <- function(ya, yb, alternative = c("twoSided", "greater", "less"),
                        na = 1, nb = 1, esType = c("none", "logOddsRatio", "difference"),
                        esMin = NULL, alpha = 0.05, prior = NULL) {
  # argument validation
  alternative <- match.arg(alternative)
  esType <- match.arg(esType)

  if (length(ya) != length(yb)) stop("The numer of samples in group A and group B must be equal!")

  if (max(ya) > na || min(ya) < 0 || max(yb) > nb || min(yb) < 0) {
    stop(paste0(
      "Number of observations cannot be negative or larger than the block size!\n",
      "    Expected:  ya in [0, ", na, "]\t", "Actual: ya in [", min(ya), ", ", max(ya), "]\n",
      "    Expected:  yb in [0, ", nb, "]\t", "Actual: yb in [", min(yb), ", ", max(yb), "]\n"
    ))
  }

  if (esType == "none" && !is.null(esMin)) {
    stop("There should be no minimum difference in this case!")
  } else if (esType == "difference") {
    if (!is.numeric(esMin)) {
      stop("For absolute difference cases, esMin must be numeric.")
    } else if (esMin < 0 || esMin > 1) {
      stop("For absolute difference cases, esMin must be between 0 and 1.")
    }
  } else if (esType == "logOddsRatio") {
    # esMin can be a float between -Inf and Inf
    if (!is.finite(esMin)) {
      stop("For logOddsRatio cases, esMin must be a finite value.\n
           Consider setting esType == \"none\" as there is no restriction placed\n
           on parameters when logOddsRatio == Inf or -Inf.")
    }
  }

  # unpack the prior values
  if (is.null(prior)) {
    # print("Using the default parameters (alpha = beta = 0.18) for Beta distribution...")
    betaA1 <- betaA2 <- betaB1 <- betaB2 <- 0.18
  } else {
    betaA1 <- prior[["betaA1"]]
    betaA2 <- prior[["betaA2"]]
    betaB1 <- prior[["betaB1"]]
    betaB2 <- prior[["betaB2"]]
  }


  # Case 1: No restriction on prior --------------------------------------------
  if (esType == "none") {
    eVal <- 1

    # naive approach
    thetaA <- thetaB <- theta0 <- 0.5

    # init vector
    totalSuccessA <- cumsum(ya)
    totalSuccessB <- cumsum(yb)
    groupSizeVecA <- seq_along(totalSuccessA) * na
    groupSizeVecB <- seq_along(totalSuccessB) * nb
    totalFailA <- groupSizeVecA - totalSuccessA
    totalFailB <- groupSizeVecB - totalSuccessB

    for (i in seq_along(ya)) {
      newE <- safestats:::calculateETwoProportions(
        na1 = ya[i], na = na, nb1 = yb[i], nb = nb, thetaA = thetaA,
        thetaB = thetaB, theta0 = theta0
      )

      eVal <- eVal * newE

      thetaA <- safestats:::bernoulliMLTwoProportions(totalSuccessA[i], totalFailA[i], betaA1, betaA2)

      thetaB <- safestats:::bernoulliMLTwoProportions(totalSuccessB[i], totalFailB[i], betaB1, betaB2)

      theta0 <- (na * thetaA + nb * thetaB) / (na + nb)
    }
  } else if (esType == "difference" || esType == "logOddsRatio") {
    # Case 2: Specified minimal effect difference (esMin) ----------------------
    # on (0, 1)^2 para. space
    # they all used Bayensian updating in each timesteps between blocks
    # the only difference are the calculation of thetaB and grid

    # avoid 0 and 1 to prevent numerical instability
    gridSize <- 1e3

    rhoGrid <- seq(1 / gridSize, 1 - 1 / gridSize, length.out = gridSize)

    ## difference in probability grid ------------------------------------------
    if (esType == "difference") {
      # prepare prob. grid = rho x (1 - esMin)
      thetaAgrid <- rhoGrid * (1 - esMin)
      thetaBgrid <- thetaAgrid + esMin
    } else if (esType == "logOddsRatio") {
      thetaAgrid <- rhoGrid
      thetaBgrid <- safestats:::calculateThetaBFromThetaAAndLOR(thetaAgrid, esMin)
    }

    # prepare log theta
    logThetaAgrid <- log(thetaAgrid)
    logThetaBgrid <- log(thetaBgrid)
    logOneMinusThetaAgrid <- log(1 - thetaAgrid)
    logOneMinusThetaBgrid <- log(1 - thetaBgrid)


    # prior density for thetaA at each grid points with normalization
    # rho = thetaA / (1 - esMin) ~ Beta(shape1, shape2)
    thetaADensity <- stats::dbeta(rhoGrid, shape1 = betaA1, shape2 = betaA2)
    thetaADensity <- thetaADensity / sum(thetaADensity)

    # TODO: Is the initial theta 1 / 2?
    thetaB <- thetaA <- theta0 <- numeric(length(ya))

    for (i in seq_along(ya)) {
      # Compute log-likelihoods for thetaA and thetaB
      logLikelihoodA <- ya[i] * logThetaAgrid + (na - ya[i]) * logOneMinusThetaAgrid
      logLikelihoodB <- yb[i] * logThetaBgrid + (nb - yb[i]) * logOneMinusThetaBgrid

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
    # TODO: add validation for esMin < 0
    if (esType == "difference") {
      thetaB <- thetaA + esMin
    } else if (esType == "logOddsRatio") {
      thetaB <- safestats:::calculateThetaBFromThetaAAndLOR(thetaA, esMin)
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
