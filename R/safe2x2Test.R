# Testing fnts ---------

## Oneshot ----
oneshotTurner <- function(ya, na, yb, nb, log = FALSE) {
  gamma <- 0.18
  # Formula: theta = shape * previousYa / (previousNa + shape + scale)
  # shape=gamma, scale=gamma, previousYa=1, previousNa=1
  thetaA <- thetaB <- gamma / (1 + 2 * gamma)
  thetaNull <- (na * thetaA + nb * thetaB) / (na + nb)

  logEValue <- dbinom(ya, na, thetaA, log = TRUE) +
    dbinom(yb, nb, thetaB, log = TRUE) -
    dbinom(ya, na, thetaNull, log = TRUE) -
    dbinom(yb, nb, thetaNull, log = TRUE)

  if (log) return(logEValue) else return(exp(logEValue))
}

oneshotNml <- function(ya, yb, na, nb, pseudo = FALSE) {
  n1 <- ya + yb

  deltaMle <- fnchMle(ya, yb, na, nb)

  ePseudo <- conditionalEValueFixedAlternative(ya, yb, na, nb, deltaMle)

  if (pseudo) {
    return(ePseudo)
  }

  # Support domain for ya
  kVals <- max(0, n1 - nb):min(na, n1)

  densities <- sapply(kVals, function(k) {
    # Calculate MLE for the hypothetical outcome k
    kMle <- fnchMle(k, n1 - k, na, nb)
    BiasedUrn::dFNCHypergeo(k, na, nb, n1, exp(kMle))
  })

  normConst <- sum(densities)

  return(ePseudo / normConst)
}

oneshotBayes <- function(ya, na, yb, nb, priorDist = stats::dnorm, ..., log = FALSE) {
  n1 <- ya + yb

  # Psi(0) is just the log of the central hypergeometric total sum
  # which is log(choose(na + nb, n1))
  psiZero <- lchoose(na + nb, n1)

  integrand <- function(delta) {
    # Vectorized calculation of the E-value kernel
    # exp( delta * ya - Psi(delta) + Psi(0) + logPrior )
    sapply(delta, function(d) {
      logWeight <- (d * ya) - fnchPsi(na, nb, n1, d) + psiZero
      exp(logWeight + priorDist(d, ..., log = TRUE))
    })
  }

  # Use a sensible finite range for the log-odds integration
  eValue <- stats::integrate(integrand, lower = -40, upper = 40)$value

  if (log) return(log(eValue))
  eValue
}


#' Conditional E-variable for fixed delta
conditionalEValueFixedAlternative <- function(ya, yb, na, nb, delta) {
  odds <- exp(delta)
  n1 <- ya + yb

  numerator <- BiasedUrn::dFNCHypergeo(
    x = ya,
    m1 = na,
    m2 = nb,
    n = n1,
    odds = odds
  )
  denominator <- stats::dhyper(x = ya, m = na, n = nb, k = n1)

  result <- numerator / denominator
  result[denominator == 0] <- 0

  return(result)
}

## Seq ----

#' Sequential Turner's method
seqTurner <- function(ya, yb, na, nb, gamma = 1.0, log = FALSE) {
  nSteps <- length(ya)

  # Previous cumulative sums
  prevYaCumsum <- c(0, cumsum(ya[-nSteps]))
  prevYbCumsum <- c(0, cumsum(yb[-nSteps]))
  prevNaCumsum <- c(0, cumsum(na[-nSteps]))
  prevNbCumsum <- c(0, cumsum(nb[-nSteps]))

  # Prior parameters (Symmetric Beta(gamma, gamma))
  alphaA <- betaA <- alphaB <- betaB <- gamma

  # Predictive means
  # Formula: theta = (alpha + sum_prev_y) / (sum_prev_n + alpha + beta)
  breveMeanA <- (alphaA + prevYaCumsum) / (prevNaCumsum + alphaA + betaA)
  breveMeanB <- (alphaB + prevYbCumsum) / (prevNbCumsum + alphaB + betaB)

  # Null predictive mean (weighted average)
  breveMeanNull <- (na * breveMeanA + nb * breveMeanB) / (na + nb)

  # Likelihood ratio (incremental)
  logNum <- stats::dbinom(ya, na, breveMeanA, log = TRUE) +
    stats::dbinom(yb, nb, breveMeanB, log = TRUE)
  logDen <- stats::dbinom(ya, na, breveMeanNull, log = TRUE) +
    stats::dbinom(yb, nb, breveMeanNull, log = TRUE)

  logEIncremental <- logNum - logDen
  logECumulative <- cumsum(logEIncremental)

  if (log) return(logECumulative) else return(exp(logECumulative))
}

#' Sequential NML method
seqNml <- function(ya, yb, na, nb, log = FALSE, pseudo = FALSE) {
  nSteps <- length(ya)

  yaCumsum <- cumsum(ya)
  ybCumsum <- cumsum(yb)
  naCumsum <- cumsum(na)
  nbCumsum <- cumsum(nb)

  prevYaCumsum <- c(0, yaCumsum[-nSteps])
  prevYbCumsum <- c(0, ybCumsum[-nSteps])

  # Global MLE sequence for the numerator
  deltaMle <- numeric(nSteps)
  for (t in 1:nSteps) {
    deltaMle[t] <- fnchMle(yaCumsum[t], ybCumsum[t], naCumsum[t], nbCumsum[t])
  }
  oddsMle <- exp(deltaMle)

  # Numerator: log P_delta(ya_t | n1_t)
  logNum <- numeric(nSteps)
  for (t in 1:nSteps) {
    logNum[t] <- log(BiasedUrn::dFNCHypergeo(ya[t], na[t], nb[t], ya[t] + yb[t], oddsMle[t]))
  }

  # Denominator: Central Hypergeometric
  logDenHg <- stats::dhyper(ya, na, nb, ya + yb, log = TRUE)

  logDiff <- logNum - logDenHg

  if (pseudo) {
    logECumulative <- cumsum(logDiff)
    if (log) return(logECumulative) else return(exp(logECumulative))
  }

  # Normalization constants (log-space)
  normConst <- numeric(nSteps)
  n1 <- ya + yb

  for (t in 1:nSteps) {
    if (n1[t] == 0) {
      normConst[t] <- 0
      next
    }

    lb <- max(0, n1[t] - nb[t])
    ub <- min(n1[t], na[t])
    if (lb == ub) {
      normConst[t] <- 0
      next
    }

    kVals <- lb:ub
    # MLE for all possible points at step t based on cumulative data
    logDensities <- sapply(kVals, function(k) {
      ya_up <- prevYaCumsum[t] + k
      yb_up <- prevYbCumsum[t] + (n1[t] - k)

      # Use tryCatch for fnchMle
      kDeltaMle <- tryOrFailWithNA(fnchMle(ya_up, yb_up, naCumsum[t], nbCumsum[t]))
      if (is.na(kDeltaMle)) return(-Inf)

      log(BiasedUrn::dFNCHypergeo(k, na[t], nb[t], n1[t], exp(kDeltaMle)))
    })

    normConst[t] <- logSumExp(logDensities)
  }

  logEIncremental <- logDiff - normConst
  logECumulative <- cumsum(logEIncremental)

  if (log) return(logECumulative) else return(exp(logECumulative))
}

seqBayes <- function(ya, yb, na, nb, priorDist = stats::dnorm, ..., log = FALSE, returnPosteriors = FALSE) {
  nSteps <- length(ya)
  n1 <- ya + yb
  deltaGrid <- seq(-10, 10, length.out = 1000)

  # Initial prior weights
  logPriorWeights <- priorDist(deltaGrid, ..., log = TRUE)
  logPriorWeights <- logPriorWeights - logSumExp(logPriorWeights)

  logECumulative <- numeric(nSteps)
  currentLogPostWeights <- logPriorWeights

  if (returnPosteriors) {
    posteriors <- matrix(NA, nrow = nSteps + 1, ncol = length(deltaGrid))
    posteriors[1, ] <- exp(logPriorWeights)
  }

  for (t in 1:nSteps) {
    # Incremental likelihood ratio log(f_delta / f_0)
    # logLR(delta) = delta * ya[t] - psi(delta) + psi(0)
    psiZero <- lchoose(na[t] + nb[t], n1[t])

    logLR <- sapply(deltaGrid, function(d) {
      (d * ya[t]) - fnchPsi(na[t], nb[t], n1[t], d) + psiZero
    })

    # Step t e-value: sum( w_{t-1} * LR_t )
    logEIncremental <- logSumExp(currentLogPostWeights + logLR)

    if (t == 1) {
      logECumulative[t] <- logEIncremental
    } else {
      logECumulative[t] <- logECumulative[t-1] + logEIncremental
    }

    # Update posterior for next step: w_t = w_{t-1} * LR_t / e_t
    currentLogPostWeights <- currentLogPostWeights + logLR - logEIncremental

    if (returnPosteriors) {
      posteriors[t+1, ] <- exp(currentLogPostWeights)
    }
  }

  res <- if (log) logECumulative else exp(logECumulative)

  if (returnPosteriors) {
    return(list(eValues = res, posteriors = posteriors, deltaGrid = deltaGrid))
  } else {
    return(res)
  }
}


# Helper ----

#' Log-sum-exp transformation
logSumExp <- function(x) {
  xMax <- max(x)
  xMax + log(sum(exp(x - xMax)))
}

#' MLE estimator for FNCH
fnchMle <- function(ya, yb, na, nb) {
  n1 <- ya + yb

  # Solve for delta: observed ya - expected value = 0
  objFunc <- function(delta) {
    ya - BiasedUrn::meanFNCHypergeo(na, nb, n1, exp(delta))
  }

  # interval c(-50, 50) is safe for log-odds; extendInt handles extremes
  result <- stats::uniroot(objFunc, interval = c(-50, 50), extendInt = "yes")
  result$root
}

#' Log Partition function for FNCH
fnchPsi <- function(na, nb, n1, delta) {
  kVals <- max(0, n1 - nb):min(na, n1)

  # Log of the base measure r(k)
  logH <- lchoose(na, kVals) + lchoose(nb, n1 - kVals)
  logTerms <- (delta * kVals) + logH

  logSumExp(logTerms)
}
