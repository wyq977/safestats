# Testing fnts ---------

#' A Wrapper for the futility
saviFutilityTwoPropCondStat <- function(
  ya,
  yb,
  na,
  nb,
  logOR = NULL,
  alternative = c("twoSided", "greater", "less"),
  designObj = NULL,
  ...
) {
  alternative <- match.arg(alternative)
  result <- list()
  result[["alternative"]] <- alternative

  if (alternative %in% c("twoSided", "greater")) {
    greaterEValues <- conditionalEValueFixedAlternative(
      ya = ya,
      yb = yb,
      na = na,
      nb = nb,
      logOR = logOR,
      alternative = "greater"
    )
  }

  if (alternative %in% c("twoSided", "less")) {
    lessEValues <- conditionalEValueFixedAlternative(
      ya = ya,
      yb = yb,
      na = na,
      nb = nb,
      logOR = logOR,
      alternative = "less"
    )
  }

  if (alternative == "greater") {
    result[["eValueVec"]] <- greaterEValues
  } else if (alternative == "less") {
    result[["eValueVec"]] <- lessEValues
  } else {
    result[["eValueVec"]] <- pmax(greaterEValues, lessEValues)
  }

  return(result)
}

saviTwoPropCondStat <- function(
  ya,
  yb,
  na,
  nb,
  logOR = NULL,
  alternative = c("twoSided", "greater", "less"),
  eType = c("eGauss", "nml", "turner", "grow", "eBayes"),
  designObj = NULL,
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
      parameter = list(logOR = logOR),
      pilot = TRUE,
      eType = eType
    )
  }

  alpha <- designObj[["alpha"]]
  alternative <- designObj[["alternative"]]
  h0 <- designObj[["h0"]]

  ### Check: Data -----
  # if (!isValid2x2Vec(ya = ya, yb = yb, na = na, nb = nb)) {
  #   stop("Input data wrong! Please check")
  # }

  nSteps <- length(ya)

  if (length(na) == 1 && nSteps > 1) {
    na <- rep(na, nSteps)
  }
  if (length(nb) == 1 && nSteps > 1) {
    nb <- rep(nb, nSteps)
  }

  dataName <- paste(deparse1(substitute(ya)), "and", deparse1(substitute(yb)))

  n <- c(sum(na), sum(nb))
  names(n) <- c("nObsA", "nObsB")

  # Vars for analysis
  eValueVec <- NULL

  ### Compute: eValue ----
  if (eType == "turner") {
    eValueVec <- turnerEProcess(ya, yb, na, nb)
    eValue <- eValueVec[nSteps]
  } else if (eType %in% c("eGauss", "eBayes")) {
    eValueVec <- computeEGaussGrid(ya, yb, na, nb)
  } else {
    eValueVec <- switch(
      eType,
      "nml" = computeNml(ya, yb, na, nb),
      "grow" = conditionalEValueFixedAlternative(
        ya,
        yb,
        na,
        nb,
        logOR,
        alternative
      ),
      rep(NA_real_, nSteps)
    )
  }

  # Other methods above return blockwise e-factors.
  if (eType != "turner") eValue <- prod(eValueVec)

  ### Fill: Result -----
  result[["n"]] <- n
  result[["eValue"]] <- eValue
  result[["dataName"]] <- dataName
  result[["alternative"]] <- alternative
  result[["testType"]] <- "2x2"
  result[["testName"]] <- "Two Proportions"
  result[["designObj"]] <- designObj
  result[["h0"]] <- h0
  result[["eValueVec"]] <- if (nSteps > 1) eValueVec else NULL
  result[["call"]] <- sys.call()

  class(result) <- "saviTest"
  return(result)
}

#' This calculate the bounds on log-odds ratio.
computeConfidenceInterval2x2 <- function(
  ya,
  yb,
  na,
  nb,
  ciValue = 0.95,
  alternative = c("twoSided", "greater", "less"),
  eType = c("eGauss", "nml", "turner", "grow", "eBayes")
) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  trivialConfidenceInterval <- switch(
    alternative,
    "twoSided" = c(-Inf, Inf),
    "greater" = c(0, Inf),
    "less" = c(-Inf, 0)
  )

  return(trivialConfidenceInterval)
}


## Math Functions ----

#' Sequential conditional Gaussian-mixture E-values
#'
#' The initial standard-normal prior on the log odds ratio is represented on
#' a fixed grid. Each block's conditional likelihood ratio is averaged under
#' the grid posterior based on preceding blocks, then used to update that
#' posterior. Positive log odds ratios favour group A.
computeEGaussGrid <- function(
  ya,
  yb,
  na,
  nb,
  log = FALSE,
  returnPosteriors = FALSE
) {
  nSteps <- length(ya)
  if (length(na) == 1L) {
    na <- rep(na, nSteps)
  }
  if (length(nb) == 1L) {
    nb <- rep(nb, nSteps)
  }
  if (length(yb) != nSteps || length(na) != nSteps || length(nb) != nSteps) {
    stop("ya, yb, na, and nb must have equal lengths after recycling")
  }
  totalSuccesses <- ya + yb
  gridSize <- 2000L
  logORGrid <- seq(-20, 20, length.out = gridSize)

  # Discrete approximation to the standard-normal prior
  priorLogWeights <- stats::dnorm(logORGrid, log = TRUE)
  priorLogWeights <- priorLogWeights - logSumExp(priorLogWeights)

  logE <- numeric(nSteps)
  logPosteriorWeights <- priorLogWeights

  if (returnPosteriors) {
    posteriors <- matrix(NA, nrow = nSteps + 1, ncol = gridSize)
    posteriors[1, ] <- exp(priorLogWeights)
  }

  for (block in seq_len(nSteps)) {
    # Incremental likelihood ratio log(f_logOR / f_0), which for the Fisher
    # noncentral hypergeometric is logOR * ya - psi(logOR) + psi(0).
    logPartitionAtZero <- lchoose(na[block] + nb[block], totalSuccesses[block])

    logLikelihoodRatio <- vapply(
      logORGrid,
      function(candidateLogOR) {
        candidateLogOR * ya[block] -
          fnchPsi(
            na[block], nb[block], totalSuccesses[block], candidateLogOR
          ) +
          logPartitionAtZero
      },
      numeric(1)
    )

    # This block's e-factor: the likelihood ratio averaged under the weights
    # from blocks 1 through block - 1 only.
    logE[block] <- logSumExp(logPosteriorWeights + logLikelihoodRatio)

    # Then fold this block in, ready to predict the next one.
    logPosteriorWeights <-
      logPosteriorWeights + logLikelihoodRatio - logE[block]

    if (returnPosteriors) {
      posteriors[block + 1, ] <- exp(logPosteriorWeights)
    }
  }

  result <- if (log) logE else exp(logE)

  if (returnPosteriors) {
    return(list(
      eValues = result,
      posteriors = posteriors,
      logORGrid = logORGrid
    ))
  } else {
    return(result)
  }
}

#' Conditional E-variable for fixed delta
conditionalEValueFixedAlternative <- function(
  ya,
  yb,
  na,
  nb,
  logOR = NULL,
  alternative = c("twoSided", "greater", "less")
) {
  if (is.null(logOR)) {
    stop("No logOR is given!")
  }

  alternative <- match.arg(alternative)

  oddsRatio <- switch(
    alternative,
    "twoSided" = exp(logOR),
    "greater" = exp(abs(logOR)),
    "less" = exp(-abs(logOR))
  )

  totalSuccesses <- ya + yb

  if (length(oddsRatio) == 1 && length(ya) > 1) {
    oddsRatio <- rep(oddsRatio, length(ya))
  }

  # BiasedUrn takes only scalar,
  numerator <- mapply(
    BiasedUrn::dFNCHypergeo,
    x = ya,
    m1 = na,
    m2 = nb,
    n = totalSuccesses,
    odds = oddsRatio
  )
  denominator <- stats::dhyper(x = ya, m = na, n = nb, k = totalSuccesses)

  result <- numerator / denominator
  result[denominator == 0] <- 0

  return(result)
}

#' Sequential conditional plug-in E-values
#'
#' For block one, no past estimate is available and the E-factor is one. For
#' each later block, the conditional maximum-likelihood estimate of the log
#' odds ratio is calculated from all preceding blocks and plugged into the
#' Fisher noncentral-hypergeometric versus central-hypergeometric likelihood
#' ratio for the current block. Positive log odds ratios favour group A.
seqCond <- function(ya, yb, na, nb, log = FALSE) {
  nSteps <- length(ya)
  if (length(na) == 1L) {
    na <- rep(na, nSteps)
  }
  if (length(nb) == 1L) {
    nb <- rep(nb, nSteps)
  }
  if (length(yb) != nSteps || length(na) != nSteps || length(nb) != nSteps) {
    stop("ya, yb, na, and nb must have equal lengths after recycling")
  }
  if (!all(is.finite(c(ya, yb, na, nb)))) {
    stop("ya, yb, na, and nb must be finite")
  }
  if (
    any(na <= 0) || any(nb <= 0) || any(ya < 0) || any(yb < 0) ||
      any(ya > na) || any(yb > nb) ||
      any(c(ya, yb, na, nb) %% 1 != 0)
  ) {
    stop("ya, yb, na, and nb must define valid 2x2 tables")
  }

  logE <- numeric(nSteps)
  if (nSteps < 2L) {
    return(if (log) logE else exp(logE))
  }

  totalSuccesses <- ya + yb
  cumulativeYa <- cumsum(ya)
  cumulativeYb <- cumsum(yb)
  cumulativeNa <- cumsum(na)
  cumulativeNb <- cumsum(nb)

  for (block in 2:nSteps) {
    previousBlock <- block - 1L
    previousSuccesses <- cumulativeYa[previousBlock] +
      cumulativeYb[previousBlock]
    previousLower <- max(0, previousSuccesses - cumulativeNb[previousBlock])
    previousUpper <- min(previousSuccesses, cumulativeNa[previousBlock])

    if (previousLower == previousUpper) {
      mleLogOR <- 0
    } else if (cumulativeYa[previousBlock] == previousLower) {
      mleLogOR <- -Inf
    } else if (cumulativeYa[previousBlock] == previousUpper) {
      mleLogOR <- Inf
    } else {
      mleLogOR <- fnchMle(
        cumulativeYa[previousBlock],
        cumulativeYb[previousBlock],
        cumulativeNa[previousBlock],
        cumulativeNb[previousBlock]
      )
    }

    currentLower <- max(0, totalSuccesses[block] - nb[block])
    currentUpper <- min(totalSuccesses[block], na[block])
    if (currentLower == currentUpper) {
      logE[block] <- 0
    } else if (is.infinite(mleLogOR)) {
      mleSupportPoint <- if (mleLogOR < 0) currentLower else currentUpper
      logE[block] <- if (ya[block] == mleSupportPoint) {
        lchoose(na[block] + nb[block], totalSuccesses[block]) -
          lchoose(na[block], ya[block]) -
          lchoose(nb[block], yb[block])
      } else {
        -Inf
      }
    } else {
      logE[block] <-
        mleLogOR * ya[block] -
        fnchPsi(na[block], nb[block], totalSuccesses[block], mleLogOR) +
        lchoose(na[block] + nb[block], totalSuccesses[block])
    }
  }

  if (log) logE else exp(logE)
}

#' NML method
computeNml <- function(ya, yb, na, nb, log = FALSE, pseudo = FALSE) {
  nSteps <- length(ya)
  if (length(na) == 1L) {
    na <- rep(na, nSteps)
  }
  if (length(nb) == 1L) {
    nb <- rep(nb, nSteps)
  }
  if (length(yb) != nSteps || length(na) != nSteps || length(nb) != nSteps) {
    stop("ya, yb, na, and nb must have equal lengths after recycling")
  }

  cumulativeYa <- cumsum(ya)
  cumulativeYb <- cumsum(yb)
  cumulativeNa <- cumsum(na)
  cumulativeNb <- cumsum(nb)

  if (nSteps == 1) {
    previousYa <- 0
    previousYb <- 0
  } else {
    previousYa <- c(0, cumulativeYa[-nSteps])
    previousYb <- c(0, cumulativeYb[-nSteps])
  }

  # Global MLE sequence for the numerator
  mleLogOR <- numeric(nSteps)
  for (block in 1:nSteps) {
    mleLogOR[block] <- fnchMle(
      cumulativeYa[block],
      cumulativeYb[block],
      cumulativeNa[block],
      cumulativeNb[block]
    )
  }
  mleOdds <- exp(mleLogOR)

  # Numerator: log P_logOR(ya[block] | totalSuccesses[block])
  logNumerator <- numeric(nSteps)
  for (block in 1:nSteps) {
    logNumerator[block] <- log(BiasedUrn::dFNCHypergeo(
      ya[block],
      na[block],
      nb[block],
      ya[block] + yb[block],
      mleOdds[block]
    ))
  }

  # Denominator: Central Hypergeometric
  logDenominator <- stats::dhyper(ya, na, nb, ya + yb, log = TRUE)

  logLikelihoodRatio <- logNumerator - logDenominator

  if (pseudo) {
    if (log) {
      return(logLikelihoodRatio)
    } else {
      return(exp(logLikelihoodRatio))
    }
  }

  # Normalisation constants (log-space)
  logNormalisingConstant <- numeric(nSteps)
  totalSuccesses <- ya + yb

  for (block in 1:nSteps) {
    if (totalSuccesses[block] == 0) {
      logNormalisingConstant[block] <- 0
      next
    }

    feasibleLower <- max(0, totalSuccesses[block] - nb[block])
    feasibleUpper <- min(totalSuccesses[block], na[block])
    if (feasibleLower == feasibleUpper) {
      logNormalisingConstant[block] <- 0
      next
    }

    feasibleSuccesses <- feasibleLower:feasibleUpper
    # MLE at every outcome this block could have taken, given the cumulative
    # data that preceded it.
    logDensities <- vapply(feasibleSuccesses, function(candidateYa) {
      candidateCumulativeYa <- previousYa[block] + candidateYa
      candidateCumulativeYb <- previousYb[block] +
        (totalSuccesses[block] - candidateYa)

      candidateMleLogOR <- tryOrFailWithNA(fnchMle(
        candidateCumulativeYa,
        candidateCumulativeYb,
        cumulativeNa[block],
        cumulativeNb[block]
      ))
      if (is.na(candidateMleLogOR)) {
        return(-Inf)
      }

      log(BiasedUrn::dFNCHypergeo(
        candidateYa,
        na[block],
        nb[block],
        totalSuccesses[block],
        exp(candidateMleLogOR)
      ))
    }, numeric(1))

    logNormalisingConstant[block] <- logSumExp(logDensities)
  }

  logE <- logLikelihoodRatio - logNormalisingConstant

  if (log) return(logE) else return(exp(logE))
}

# Helper ----

#' Log-sum-exp transformation
logSumExp <- function(x) {
  xMax <- max(x)
  xMax + log(sum(exp(x - xMax)))
}

#' MLE estimator for FNCH
fnchMle <- function(ya, yb, na, nb) {
  totalSuccesses <- ya + yb

  # At the boundary of the support the MLE is infinite and uniroot's extended
  # search runs into odds that BiasedUrn rejects. Return the ends of the search
  # interval instead: exp(50) is far enough out that the density at the
  # boundary point is one to machine precision.
  feasibleLower <- max(0, totalSuccesses - nb)
  feasibleUpper <- min(na, totalSuccesses)
  if (feasibleLower == feasibleUpper) {
    return(0)
  }
  if (ya <= feasibleLower) {
    return(-50)
  }
  if (ya >= feasibleUpper) {
    return(50)
  }

  # Solve for logOR: observed ya - expected value = 0
  objectiveFunction <- function(logOR) {
    ya - BiasedUrn::meanFNCHypergeo(na, nb, totalSuccesses, exp(logOR))
  }

  # interval c(-50, 50) is safe for log-odds; extendInt handles extremes
  result <- stats::uniroot(
    objectiveFunction,
    interval = c(-50, 50),
    extendInt = "yes"
  )
  result$root
}

#' Log Partition function for FNCH
fnchPsi <- function(na, nb, totalSuccesses, logOR) {
  feasibleSuccesses <-
    max(0, totalSuccesses - nb):min(na, totalSuccesses)

  # Log of the base measure
  logBaseMeasure <- lchoose(na, feasibleSuccesses) +
    lchoose(nb, totalSuccesses - feasibleSuccesses)
  logTerms <- (logOR * feasibleSuccesses) + logBaseMeasure

  logSumExp(logTerms)
}

#' Validate data for 2 stream data
isValid2x2Vec <- function(ya, yb, na, nb) {
  if (!all(is.finite(c(ya, yb, na, nb)))) {
    warning(
      "Data contains infinite value!"
    )
    return(FALSE)
  }

  if (!all(na > 0, nb > 0, ya >= 0, yb >= 0)) {
    warning(
      "Observations and group size must be greater than 0!"
    )
    return(FALSE)
  }

  if (!all(ya %% 1 == 0, yb %% 1 == 0, na %% 1 == 0, nb %% 1 == 0)) {
    warning(
      "Observations and group size must be integer!"
    )
    return(FALSE)
  }

  if (!all(ya <= na, yb <= nb)) {
    warning(
      "Observations must be smaller than group size!"
    )
    return(FALSE)
  }
}
