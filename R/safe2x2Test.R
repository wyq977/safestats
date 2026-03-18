# Testing fnts ---------

saviTwoPropCondStat <- function(
  ya,
  yb,
  na,
  nb,
  logOddsRatio = NULL,
  alternative = c("twoSided", "less", "greater"),
  eType = c("eGauss", "nml", "turner", "grow"),
  designObj = NULL,
  sequential = NULL,
  debug = FALSE,
  ...
) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  result <- constructSaviTestObj("Two Proportions")

  ### Check: designObj ----
  if (is.null(designObj)) {
    warningMessage <- paste(
      "No designObj given. Default pilot test computed",
      "based on default settings."
    )
    warning(warningMessage)

    designObj <- list(
      alpha = 0.05,
      alternative = alternative,
      h0 = 0,
      testType = "2x2",
      testName = "Two Proportions",
      parameter = list(logOddsRatio = logOddsRatio),
      pilot = TRUE,
      eType = eType
    )
  }

  alpha <- designObj[["alpha"]]
  alternative <- designObj[["alternative"]]
  h0 <- designObj[["h0"]]

  ### Check: Data -----
  if (!isValid2x2Vec(ya, yb, na, nb)) {
    stop("Input data wrong! Please check")
  }

  mIter <- length(ya)

  if (length(na) == 1 && mIter > 1) {
    na <- rep(na, mIter)
  }
  if (length(nb) == 1 && mIter > 1) {
    nb <- rep(nb, mIter)
  }

  if (is.null(sequential)) {
    sequential <- if (mIter > 1) TRUE else FALSE
  }

  # Aggregate the data if not sequential but multiple vectors provided
  if (!sequential && mIter > 1) {
    warningMessage <- paste(
      "Data has been aggregated into a single table!"
    )
    warning(warningMessage)

    ya <- sum(ya)
    yb <- sum(yb)
    na <- sum(na)
    nb <- sum(nb)
    mIter <- 1
  }

  dataName <- paste(deparse1(substitute(ya)), "and", deparse1(substitute(yb)))

  n <- c(sum(na), sum(nb))
  names(n) <- c("nObsA", "nObsB")

  ### Early Return: Debug -----
  if (debug) {
    result[["n"]] <- n
    result[["eValue"]] <- NA_real_
    result[["dataName"]] <- dataName
    result[["alternative"]] <- alternative
    result[["testType"]] <- "2x2"
    result[["testName"]] <- "Two Proportions"
    result[["designObj"]] <- designObj
    result[["h0"]] <- h0
    result[["eValueVec"]] <- if (sequential) numeric(mIter) else NULL
    result[["confSeqMatrix"]] <- if (sequential) matrix(NA_real_, nrow = mIter, ncol = 2) else NULL
    result[["eValueApproxError"]] <- NULL
    result[["call"]] <- sys.call()

    class(result) <- "saviTest"
    return(result)
  }

  # Vars for analysis
  eValueVec <- NULL
  confSeqMatrix <- NULL
  eValueApproxError <- NULL

  ### Compute: eValue ----
  if (eType == "eGauss") {
    if (mIter == 1) {
      bayesRes <- oneshotBayes(ya, yb, na, nb, log = FALSE, returnError = TRUE)
      eValueVec <- bayesRes[["eValue"]]
      eValueApproxError <- bayesRes[["abs.error"]]
    } else {
      eValueVec <- seqBayes(ya, yb, na, nb)
    }
  } else {
    eValueVec <- switch(eType,
      "turner" = computeTurner(ya, yb, na, nb),
      "nml"    = computeNml(ya, yb, na, nb),
      "grow"   = conditionalEValueFixedAlternative(ya, yb, na, nb, logOddsRatio),
      rep(NA_real_, mIter)
    )
  }

  # Final e-value is the last one in the sequence
  eValue <- eValueVec[mIter]

  ### Compute: Confidence Intervals ----
  if (sequential) {
    confSeqMatrix <- matrix(nrow = mIter, ncol = 2)
    for (i in seq_along(ya)) {
      # TODO: Implement sequence CI math
      confSeqMatrix[i, ] <- c(NA_real_, NA_real_)
    }
  }

  ### Fill: Result -----
  result[["n"]] <- n
  result[["eValue"]] <- eValue
  result[["dataName"]] <- dataName
  result[["alternative"]] <- alternative
  result[["testType"]] <- "2x2"
  result[["testName"]] <- "Two Proportions"
  result[["designObj"]] <- designObj
  result[["h0"]] <- h0
  result[["eValueVec"]] <- if (sequential) eValueVec else NULL
  result[["confSeqMatrix"]] <- confSeqMatrix
  result[["eValueApproxError"]] <- eValueApproxError
  result[["call"]] <- sys.call()

  class(result) <- "saviTest"
  return(result)
}

## Math Functions ----

oneshotBayes <- function(
  ya,
  yb,
  na,
  nb,
  priorDist = stats::dnorm,
  ...,
  log = FALSE,
  returnError = FALSE
) {
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
  res <- stats::integrate(integrand, lower = -40, upper = 40)
  eValue <- res$value

  if (log) {
    eValue <- log(eValue)
  }

  if (returnError) {
    return(list(eValue = eValue, abs.error = res$abs.error))
  }
  eValue
}

#' Conditional E-variable for fixed delta
conditionalEValueFixedAlternative <- function(ya, yb, na, nb, logOddsRatio) {
  if (is.null(logOddsRatio)) {
    logOddsRatio <- 0
  }
  odds <- exp(logOddsRatio)
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

#' Turner's method
computeTurner <- function(ya, yb, na, nb, priorValues = NULL, log = FALSE) {
  nSteps <- length(ya)

  if (is.null(priorValues)) {
    # Default to REGRET optimal based on first data block sample sizes
    priorValues <- list(
      betaA1 = 0.18, 
      betaA2 = 0.18, 
      betaB1 = (nb[1]/na[1]) * 0.18, 
      betaB2 = (nb[1]/na[1]) * 0.18
    )
  }

  # unpack the prior values
  betaA1 <- priorValues[["betaA1"]]
  betaA2 <- priorValues[["betaA2"]]
  betaB1 <- priorValues[["betaB1"]]
  betaB2 <- priorValues[["betaB2"]]

  # Previous cumulative sums
  if (nSteps == 1) {
    prevYaCumsum <- 0
    prevYbCumsum <- 0
    prevNaCumsum <- 0
    prevNbCumsum <- 0
  } else {
    prevYaCumsum <- c(0, cumsum(ya[-nSteps]))
    prevYbCumsum <- c(0, cumsum(yb[-nSteps]))
    prevNaCumsum <- c(0, cumsum(na[-nSteps]))
    prevNbCumsum <- c(0, cumsum(nb[-nSteps]))
  }

  # Predictive means
  breveMeanA <- (betaA1 + prevYaCumsum) / (prevNaCumsum + betaA1 + betaA2)
  breveMeanB <- (betaB1 + prevYbCumsum) / (prevNbCumsum + betaB1 + betaB2)

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

#' NML method
computeNml <- function(ya, yb, na, nb, log = FALSE, pseudo = FALSE) {
  nSteps <- length(ya)

  yaCumsum <- cumsum(ya)
  ybCumsum <- cumsum(yb)
  naCumsum <- cumsum(na)
  nbCumsum <- cumsum(nb)

  if (nSteps == 1) {
    prevYaCumsum <- 0
    prevYbCumsum <- 0
  } else {
    prevYaCumsum <- c(0, yaCumsum[-nSteps])
    prevYbCumsum <- c(0, ybCumsum[-nSteps])
  }

  # Global MLE sequence for the numerator
  deltaMle <- numeric(nSteps)
  for (t in 1:nSteps) {
    deltaMle[t] <- fnchMle(yaCumsum[t], ybCumsum[t], naCumsum[t], nbCumsum[t])
  }
  oddsMle <- exp(deltaMle)

  # Numerator: log P_delta(ya_t | n1_t)
  logNum <- numeric(nSteps)
  for (t in 1:nSteps) {
    logNum[t] <- log(BiasedUrn::dFNCHypergeo(
      ya[t],
      na[t],
      nb[t],
      ya[t] + yb[t],
      oddsMle[t]
    ))
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
      if (t == 1) {
        ya_up <- k
        yb_up <- (n1[t] - k)
      } else {
        ya_up <- prevYaCumsum[t] + k
        yb_up <- prevYbCumsum[t] + (n1[t] - k)
      }

      # Use tryCatch for fnchMle
      kDeltaMle <- tryOrFailWithNA(fnchMle(
        ya_up,
        yb_up,
        naCumsum[t],
        nbCumsum[t]
      ))
      if (is.na(kDeltaMle)) {
        return(-Inf)
      }

      log(BiasedUrn::dFNCHypergeo(k, na[t], nb[t], n1[t], exp(kDeltaMle)))
    })

    normConst[t] <- logSumExp(logDensities)
  }

  logEIncremental <- logDiff - normConst
  logECumulative <- cumsum(logEIncremental)

  if (log) return(logECumulative) else return(exp(logECumulative))
}


seqBayes <- function(
  ya,
  yb,
  na,
  nb,
  priorDist = stats::dnorm,
  ...,
  log = FALSE,
  returnPosteriors = FALSE
) {
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
      logECumulative[t] <- logECumulative[t - 1] + logEIncremental
    }

    # Update posterior for next step: w_t = w_{t-1} * LR_t / e_t
    currentLogPostWeights <- currentLogPostWeights + logLR - logEIncremental

    if (returnPosteriors) {
      posteriors[t + 1, ] <- exp(currentLogPostWeights)
    }
  }

  res <- if (log) logECumulative else exp(logECumulative)

  if (returnPosteriors) {
    return(list(eValues = res, posteriors = posteriors, deltaGrid = deltaGrid))
  } else {
    return(res)
  }
}

# Design fnts ----

#' Designs a Savi Experiment to Test Two Proportions
#'
#' @param na number of observations in group a per data block
#' @param nb number of observations in group b per data block
#' @param alpha numeric in (0, 1) that specifies the tolerable type I error control
#' @param beta numeric in (0, 1) that specifies the tolerable type II error control
#' @param logOddsRatio true log odds ratio to detect
#' @param nBlocksPlan planned number of data blocks collected
#' @param alternative a character string specifying the alternative hypothesis
#' @param eType character specifying the e-variable type (eGauss, nml, turner, grow)
#'
#' @return Returns a 'saviDesign' object
#' @export
designSaviTwoProportions <- function(
  na = 1, nb = 1,
  alpha = 0.05,
  beta = NULL,
  logOddsRatio = NULL,
  nBlocksPlan = NULL,
  alternative = c("twoSided", "less", "greater"),
  eType = c("eGauss", "nml", "turner", "grow"),
  pilot = FALSE
) {
  warning("TODO: none restriction on parameters now!")
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  result <- constructSaviDesignObj("Two Proportions")

  # Design scenarios logic (simplified for structural compliance)
  # Scenario 1a: logOddsRatio + beta known -> Find nBlocksPlan
  if (!is.null(logOddsRatio) && !is.null(beta) && is.null(nBlocksPlan)) {
    designScenario <- "1a"
    # TODO: run simulation/math to find worstCaseQuantile (nBlocksPlan)
    nBlocksPlan <- NA_real_

  # Scenario 1b: pilot (only nBlocksPlan known)
  } else if (!is.null(nBlocksPlan) && is.null(logOddsRatio) && is.null(beta)) {
    designScenario <- "1b"
    pilot <- TRUE

  # Scenario 2: logOddsRatio + nBlocksPlan known -> Find beta
  } else if (!is.null(nBlocksPlan) && !is.null(logOddsRatio) && is.null(beta)) {
    designScenario <- "2"
    # TODO: run simulation/math to find worstCasePower (beta)
    beta <- NA_real_

  # Scenario 3: beta + nBlocksPlan known -> Find logOddsRatio
  } else if (!is.null(nBlocksPlan) && is.null(logOddsRatio) && !is.null(beta)) {
    designScenario <- "3"
    # TODO: run simulation/math to find min effect size
    logOddsRatio <- NA_real_

  } else if (pilot) {
    designScenario <- "pilot"
  } else {
    stop("Provide two of: nBlocksPlan, logOddsRatio, and beta. Or set pilot = TRUE.")
  }

  nPlan <- c(na, nb, nBlocksPlan)
  names(nPlan) <- c("na", "nb", "nBlocksPlan")

  esMin <- logOddsRatio
  if (!is.null(esMin)) {
    names(esMin) <- "log odds ratio"
  }

  tempResult <- list(
    "nPlan" = nPlan,
    "alpha" = alpha,
    "beta" = beta,
    "esMin" = esMin,
    "h0" = 0,
    "testType" = "2x2",
    "alternative" = alternative,
    "pilot" = pilot,
    "designScenario" = designScenario,
    "eType" = eType,
    "call" = sys.call()
  )

  # Merge using modifyList and filter out NULLs (Linus/zTest style)
  result <- utils::modifyList(result, tempResult)
  result <- Filter(Negate(is.null), result)
  class(result) <- "saviDesign"

  return(result)
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

#' Validate data for 2 stream data
isValid2x2Vec <- function(ya, yb, na, nb) {
  # 1. Ensure all are finite and non-missing
  if (!all(is.finite(c(ya, yb, na, nb)))) return(FALSE)

  # 2. Vectorized check for:
  # - Non-negative integers (ya, yb)
  # - Positive integers (na, nb)
  # - Data does not exceed group size
  all(
    ya >= 0, yb >= 0,
    na > 0,  nb > 0,
    ya %% 1 == 0, yb %% 1 == 0,
    na %% 1 == 0, nb %% 1 == 0,
    ya <= na, yb <= nb
  )
}
