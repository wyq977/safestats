# Turner e-process ----

#' Cumulative log-likelihood-ratio process for two Bernoulli streams
#'
#' Calculates the cumulative log likelihood ratio for given predictable
#' numerator and denominator probabilities. The caller is responsible for
#' constructing all four probability vectors using only information available
#' before the corresponding data block.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param na,nb Number of observations in groups A and B in each data block.
#' @param numeratorThetaA,numeratorThetaB Predictable Bernoulli probabilities
#'   in the numerator.
#' @param denominatorThetaA,denominatorThetaB Predictable Bernoulli
#'   probabilities in the denominator.
#'
#' @return A numeric vector containing the cumulative log likelihood ratio
#'   after every block.
logLikelihoodRatioProcess <- function(
  ya,
  yb,
  na,
  nb,
  numeratorThetaA,
  numeratorThetaB,
  denominatorThetaA,
  denominatorThetaB
) {
  nSteps <- length(ya)
  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)
  inputLengths <- c(
    length(yb),
    length(na),
    length(nb),
    length(numeratorThetaA),
    length(numeratorThetaB),
    length(denominatorThetaA),
    length(denominatorThetaB)
  )
  if (!all(inputLengths == nSteps)) {
    stop("Data, group sizes, and all four theta vectors must have the same length!")
  }

  logNumerator <- stats::dbinom(
    ya,
    na,
    numeratorThetaA,
    log = TRUE
  ) + stats::dbinom(
    yb,
    nb,
    numeratorThetaB,
    log = TRUE
  )
  logDenominator <- stats::dbinom(
    ya,
    na,
    denominatorThetaA,
    log = TRUE
  ) + stats::dbinom(
    yb,
    nb,
    denominatorThetaB,
    log = TRUE
  )

  cumsum(logNumerator - logDenominator)
}

# Supply and validate the four Turner Beta prior parameters.
resolveTurnerPriorParameters <- function(
  na,
  nb,
  priorParameters = NULL
) {
  if (is.null(priorParameters)) {
    priorParameters <- list(
      betaA1 = 0.18,
      betaA2 = 0.18,
      betaB1 = (nb[1] / na[1]) * 0.18,
      betaB2 = (nb[1] / na[1]) * 0.18
    )
  }

  requiredNames <- c("betaA1", "betaA2", "betaB1", "betaB2")
  if (!all(requiredNames %in% names(priorParameters))) {
    stop(
      "priorParameters must contain betaA1, betaA2, betaB1, and betaB2."
    )
  }
  priorValues <- unlist(priorParameters[requiredNames], use.names = FALSE)
  if (!is.numeric(priorValues) || any(!is.finite(priorValues)) ||
      any(priorValues <= 0)) {
    stop("Turner prior parameters must be finite and greater than zero.")
  }

  priorParameters[requiredNames]
}

# Calculate the two predictable Beta posterior means before each data block.
betaPredictiveMeansTwoProportions <- function(
  ya,
  yb,
  na,
  nb,
  priorParameters
) {
  nSteps <- length(ya)
  prevYa <- c(0, cumsum(ya))[seq_len(nSteps)]
  prevYb <- c(0, cumsum(yb))[seq_len(nSteps)]
  prevNa <- c(0, cumsum(na))[seq_len(nSteps)]
  prevNb <- c(0, cumsum(nb))[seq_len(nSteps)]

  list(
    thetaA = (priorParameters[["betaA1"]] + prevYa) / (
      prevNa + priorParameters[["betaA1"]] + priorParameters[["betaA2"]]
    ),
    thetaB = (priorParameters[["betaB1"]] + prevYb) / (
      prevNb + priorParameters[["betaB1"]] + priorParameters[["betaB2"]]
    )
  )
}

#' Cumulative Turner e-process for two Bernoulli streams
#'
#' Computes the Turner e-process for `thetaA - thetaB = 0`. This is the
#' difference-zero specialization of the linear-difference likelihood-ratio
#' process. The returned vector is cumulative: element `i` uses blocks `1`
#' through `i`.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param na,nb Number of observations in groups A and B in each data block.
#'   A scalar is recycled over all blocks.
#' @param priorParameters Optional named list with `betaA1`, `betaA2`,
#'   `betaB1`, and `betaB2`.
#' @param log Return the cumulative log e-process instead of the e-process.
#'
#' @return A numeric vector containing the cumulative e-process after every
#'   block.
#' @export
turnerEProcess <- function(
  ya,
  yb,
  na,
  nb,
  priorParameters = NULL,
  log = FALSE
) {
  nSteps <- length(ya)

  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)
  if (!all(c(length(yb), length(na), length(nb)) == nSteps)) {
    stop("ya, yb, na, and nb must have the same length.")
  }

  priorParameters <- resolveTurnerPriorParameters(
    na = na,
    nb = nb,
    priorParameters = priorParameters
  )

  breveMean <- betaPredictiveMeansTwoProportions(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    priorParameters = priorParameters
  )
  breveMeanNull <- (
    na * breveMean[["thetaA"]] + nb * breveMean[["thetaB"]]
  ) / (na + nb)
  logEProcess <- logLikelihoodRatioProcess(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    numeratorThetaA = breveMean[["thetaA"]],
    numeratorThetaB = breveMean[["thetaB"]],
    denominatorThetaA = breveMeanNull,
    denominatorThetaB = breveMeanNull
  )

  if (log) logEProcess else exp(logEProcess)
}

# Linear-difference confidence sequence ----

# Solve for the thetaA coordinate of the reverse information projection onto
# the null curve thetaA - thetaB = difference at every data block.
solveLinearDifferenceRIPrThetaA <- function(
  numeratorThetaA,
  numeratorThetaB,
  na,
  nb,
  difference,
  bWeight = 1
) {
  nSteps <- length(numeratorThetaA)
  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)
  if (!all(c(length(numeratorThetaB), length(na), length(nb)) == nSteps)) {
    stop("Numerator theta vectors and group sizes must align by block.")
  }

  solveOne <- function(thetaStarA, thetaStarB, blockSizeA, blockSizeB) {
    A <- blockSizeA
    B <- blockSizeB * bWeight
    # The polynomial below is written for thetaB - thetaA, so its difference
    # is the negative of the public thetaA - thetaB convention.
    bMinusADifference <- -difference
    coefficients <- c(
      -A * thetaStarA * bMinusADifference * (1 - bMinusADifference),
      A * (
        bMinusADifference * (1 - bMinusADifference) -
          thetaStarA * (1 - 2 * bMinusADifference)
      ) + B * (bMinusADifference - thetaStarB),
      A * (1 - 2 * bMinusADifference + thetaStarA) +
        B * (1 - bMinusADifference + thetaStarB),
      -(A + B)
    )
    roots <- Re(base::polyroot(coefficients))

    # The feasible RIPr solution is the middle of the three real roots.
    sum(roots) - min(roots) - max(roots)
  }

  mapply(
    solveOne,
    numeratorThetaA,
    numeratorThetaB,
    na,
    nb,
    USE.NAMES = FALSE
  )
}

#' Calculate e-processes over a grid of linear differences
#'
#' For every candidate `delta = thetaA - thetaB`, computes the Turner
#' likelihood-ratio process whose denominator is the reverse information
#' projection onto that candidate null curve. Every process is calculated once
#' over the complete data sequence and returned on the log scale.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param na,nb Number of observations in groups A and B in each data block.
#'   A scalar is recycled over all blocks.
#' @param priorParameters Named list with `betaA1`, `betaA2`, `betaB1`, and
#'   `betaB2` for the Turner predictor.
#' @param gridSize Number of equally spaced candidate differences in the
#'   feasible open interval `(-1, 1)`.
#'
#' @return A list containing `delta` and `logEProcesses`. Rows of the log
#'   e-process matrix are data blocks and columns are candidate differences.
calculateEValuesForLinearDeltaGrid <- function(
  ya,
  yb,
  na,
  nb,
  priorParameters,
  gridSize = 100
) {
  nSteps <- length(ya)
  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)
  if (!all(c(length(yb), length(na), length(nb)) == nSteps)) {
    stop("ya, yb, na, and nb must have the same length.")
  }

  breveMean <- betaPredictiveMeansTwoProportions(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    priorParameters = priorParameters
  )
  deltaGrid <- seq_len(gridSize) * (2 / (gridSize + 1)) - 1

  logEProcesses <- matrix(
    NA_real_,
    nrow = nSteps,
    ncol = gridSize,
    dimnames = list(NULL, as.character(deltaGrid))
  )

  for (deltaIndex in seq_along(deltaGrid)) {
    delta <- deltaGrid[deltaIndex]
    denominatorThetaA <- solveLinearDifferenceRIPrThetaA(
      numeratorThetaA = breveMean[["thetaA"]],
      numeratorThetaB = breveMean[["thetaB"]],
      na = na,
      nb = nb,
      difference = delta
    )
    logEProcesses[, deltaIndex] <- logLikelihoodRatioProcess(
      ya = ya,
      yb = yb,
      na = na,
      nb = nb,
      numeratorThetaA = breveMean[["thetaA"]],
      numeratorThetaB = breveMean[["thetaB"]],
      denominatorThetaA = denominatorThetaA,
      denominatorThetaB = denominatorThetaA - delta
    )
  }

  list(delta = deltaGrid, logEProcesses = logEProcesses)
}

#' Confidence sequence for the difference between two proportions
#'
#' Computes every candidate difference's e-process once over the full data
#' sequence and inverts their running intersections at every block.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param confidenceBoundGridPrecision Number of candidate differences used for
#'   the grid approximation in the feasible open interval `(-1, 1)`.
#' @param saviDesign A `saviDesign` returned by
#'   [designSaviTwoProportions()].
#'
#' @return A data frame with `block`, `lowerBound`, and `upperBound`. If no
#'   grid point remains at a block, both bounds are `NA`.
computeConfidenceSequenceForDifferenceTwoProportions <- function(
  ya,
  yb,
  confidenceBoundGridPrecision,
  saviDesign
) {
  if (length(confidenceBoundGridPrecision) != 1L ||
      !is.finite(confidenceBoundGridPrecision) ||
      confidenceBoundGridPrecision < 1L ||
      confidenceBoundGridPrecision %% 1 != 0) {
    stop("confidenceBoundGridPrecision must be a positive integer.")
  }

  gridProcesses <- calculateEValuesForLinearDeltaGrid(
    ya = ya,
    yb = yb,
    na = saviDesign[["nPlan"]][["na"]],
    nb = saviDesign[["nPlan"]][["nb"]],
    priorParameters = saviDesign[["betaPriorParameterValues"]],
    gridSize = confidenceBoundGridPrecision
  )
  deltaGrid <- gridProcesses[["delta"]]
  logEProcesses <- gridProcesses[["logEProcesses"]]
  nSteps <- nrow(logEProcesses)
  logThreshold <- log(1 / saviDesign[["alpha"]])
  retained <- matrix(
    FALSE,
    nrow = confidenceBoundGridPrecision,
    ncol = nSteps
  )

  for (deltaIndex in seq_along(deltaGrid)) {
    retained[deltaIndex, ] <-
      cummax(logEProcesses[, deltaIndex]) < logThreshold
  }

  lowerBound <- upperBound <- rep(NA_real_, nSteps)
  for (block in seq_len(nSteps)) {
    retainedDelta <- deltaGrid[retained[, block]]
    if (length(retainedDelta) > 0L) {
      lowerBound[block] <- min(retainedDelta)
      upperBound[block] <- max(retainedDelta)
    }
  }

  data.frame(
    block = seq_len(nSteps),
    lowerBound = lowerBound,
    upperBound = upperBound
  )
}

# Precompute candidate-independent quantities for adaptive crossing searches.
prepareLinearDifferenceCrossingSearch <- function(
  ya,
  yb,
  na,
  nb,
  priorParameters
) {
  nSteps <- length(ya)
  if (nSteps < 1L) {
    stop("ya and yb must contain at least one data block.")
  }
  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)
  if (!all(c(length(yb), length(na), length(nb)) == nSteps)) {
    stop("ya, yb, na, and nb must have the same length.")
  }

  breveMean <- betaPredictiveMeansTwoProportions(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    priorParameters = priorParameters
  )
  logNumerator <- stats::dbinom(
    ya,
    na,
    breveMean[["thetaA"]],
    log = TRUE
  ) + stats::dbinom(
    yb,
    nb,
    breveMean[["thetaB"]],
    log = TRUE
  )

  list(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    breveMean = breveMean,
    logNumerator = logNumerator,
    nSteps = nSteps
  )
}

# Evaluate one grid candidate from block one up to a known upper crossing-time
# bound. The prefix cannot be skipped because its denominator depends on delta.
firstCrossingTimeLinearDelta <- function(
  preparedSearch,
  delta,
  logThreshold,
  lowerCrossingTime = 1L,
  upperCrossingTime = preparedSearch[["nSteps"]] + 1L
) {
  nSteps <- preparedSearch[["nSteps"]]
  upperBlock <- min(nSteps, upperCrossingTime)
  blockIndices <- seq_len(upperBlock)
  denominatorThetaA <- solveLinearDifferenceRIPrThetaA(
    numeratorThetaA = preparedSearch[["breveMean"]][["thetaA"]][blockIndices],
    numeratorThetaB = preparedSearch[["breveMean"]][["thetaB"]][blockIndices],
    na = preparedSearch[["na"]][blockIndices],
    nb = preparedSearch[["nb"]][blockIndices],
    difference = delta
  )
  logDenominator <- stats::dbinom(
    preparedSearch[["ya"]][blockIndices],
    preparedSearch[["na"]][blockIndices],
    denominatorThetaA,
    log = TRUE
  ) + stats::dbinom(
    preparedSearch[["yb"]][blockIndices],
    preparedSearch[["nb"]][blockIndices],
    denominatorThetaA - delta,
    log = TRUE
  )
  logEProcess <- cumsum(
    preparedSearch[["logNumerator"]][blockIndices] - logDenominator
  )
  crossingTime <- match(TRUE, logEProcess >= logThreshold, nomatch = 0L)
  if (crossingTime == 0L) crossingTime <- nSteps + 1L

  list(
    crossingTime = crossingTime,
    blocksEvaluated = upperBlock,
    respectsBounds = crossingTime >= lowerCrossingTime &&
      crossingTime <= upperCrossingTime
  )
}

# Reconstruct grid confidence bounds without a candidate-by-block matrix.
linearDifferenceSequenceFromCrossingTimes <- function(
  deltaGrid,
  crossingTimes,
  nSteps
) {
  lowerBound <- upperBound <- rep(NA_real_, nSteps)
  for (block in seq_len(nSteps)) {
    retainedIndices <- which(crossingTimes > block)
    if (length(retainedIndices) > 0L) {
      lowerBound[block] <- deltaGrid[min(retainedIndices)]
      upperBound[block] <- deltaGrid[max(retainedIndices)]
    }
  }

  data.frame(
    block = seq_len(nSteps),
    lowerBound = lowerBound,
    upperBound = upperBound
  )
}

# Check the unimodal crossing-time and contiguous-retention assumptions on a
# completely evaluated grid.
checkLinearDifferenceCrossingStructure <- function(crossingTimes, nSteps) {
  peakTime <- max(crossingTimes)
  peakIndices <- which(crossingTimes == peakTime)
  firstPeak <- peakIndices[1L]
  lastPeak <- peakIndices[length(peakIndices)]
  lowerMonotone <- firstPeak == 1L ||
    all(diff(crossingTimes[seq_len(firstPeak)]) >= 0L)
  upperIndices <- seq.int(lastPeak, length(crossingTimes))
  upperMonotone <- length(upperIndices) == 1L ||
    all(diff(crossingTimes[upperIndices]) <= 0L)
  contiguous <- all(vapply(seq_len(nSteps), function(block) {
    retainedIndices <- which(crossingTimes > block)
    length(retainedIndices) < 2L || all(diff(retainedIndices) == 1L)
  }, logical(1)))

  list(
    unimodal = lowerMonotone && upperMonotone,
    contiguous = contiguous,
    firstPeakIndex = firstPeak,
    lastPeakIndex = lastPeak
  )
}

# Internal implementation returning instrumentation for tests and benchmarks.
computeConfidenceSequenceForDifferenceTwoProportionsAdaptiveDetails <- function(
  ya,
  yb,
  confidenceBoundGridPrecision,
  saviDesign,
  coarseGridSize = 8L
) {
  if (length(confidenceBoundGridPrecision) != 1L ||
      !is.finite(confidenceBoundGridPrecision) ||
      confidenceBoundGridPrecision < 1L ||
      confidenceBoundGridPrecision %% 1 != 0) {
    stop("confidenceBoundGridPrecision must be a positive integer.")
  }
  if (length(coarseGridSize) != 1L || !is.finite(coarseGridSize) ||
      coarseGridSize < 2L || coarseGridSize %% 1 != 0) {
    stop("coarseGridSize must be an integer of at least two.")
  }

  gridSize <- as.integer(confidenceBoundGridPrecision)
  nSteps <- length(ya)
  deltaGrid <- seq_len(gridSize) * (2 / (gridSize + 1)) - 1
  logThreshold <- log(1 / saviDesign[["alpha"]])
  preparedSearch <- prepareLinearDifferenceCrossingSearch(
    ya = ya,
    yb = yb,
    na = saviDesign[["nPlan"]][["na"]],
    nb = saviDesign[["nPlan"]][["nb"]],
    priorParameters = saviDesign[["betaPriorParameterValues"]]
  )
  crossingTimes <- rep(NA_integer_, gridSize)
  evaluated <- rep(FALSE, gridSize)
  candidateBlockEvaluations <- 0L
  candidateEvaluations <- 0L
  finalBlockEvaluations <- 0L
  refinements <- 0L
  maximumRefinementDepth <- 0L
  fallbackReason <- NA_character_

  evaluateIndex <- function(
    index,
    lowerCrossingTime = 1L,
    upperCrossingTime = nSteps + 1L
  ) {
    if (evaluated[index]) {
      return(crossingTimes[index])
    }
    result <- firstCrossingTimeLinearDelta(
      preparedSearch = preparedSearch,
      delta = deltaGrid[index],
      logThreshold = logThreshold,
      lowerCrossingTime = lowerCrossingTime,
      upperCrossingTime = upperCrossingTime
    )
    candidateEvaluations <<- candidateEvaluations + 1L
    candidateBlockEvaluations <<-
      candidateBlockEvaluations + result[["blocksEvaluated"]]
    if (result[["blocksEvaluated"]] == nSteps) {
      finalBlockEvaluations <<- finalBlockEvaluations + 1L
    }
    if (!result[["respectsBounds"]]) {
      fallbackReason <<- "monotone crossing-time bound failed"
    }
    crossingTimes[index] <<- result[["crossingTime"]]
    evaluated[index] <<- TRUE
    crossingTimes[index]
  }

  fallbackToFullGrid <- function(reason) {
    fullGrid <- calculateEValuesForLinearDeltaGrid(
      ya = ya,
      yb = yb,
      na = saviDesign[["nPlan"]][["na"]],
      nb = saviDesign[["nPlan"]][["nb"]],
      priorParameters = saviDesign[["betaPriorParameterValues"]],
      gridSize = gridSize
    )
    fullCrossingTimes <- apply(
      fullGrid[["logEProcesses"]] >= logThreshold,
      2L,
      function(crossed) {
        crossingTime <- match(TRUE, crossed, nomatch = 0L)
        if (crossingTime == 0L) nSteps + 1L else crossingTime
      }
    )
    list(
      confidenceSequence = linearDifferenceSequenceFromCrossingTimes(
        deltaGrid = fullGrid[["delta"]],
        crossingTimes = fullCrossingTimes,
        nSteps = nSteps
      ),
      crossingTimes = fullCrossingTimes,
      diagnostics = list(
        fallback = TRUE,
        fallbackReason = reason,
        distinctDeltaCandidatesEvaluated = gridSize,
        deltaCandidateEvaluations = candidateEvaluations + gridSize,
        candidatePrefixReplays = candidateEvaluations + gridSize,
        candidateBlockEvaluations = candidateBlockEvaluations + gridSize * nSteps,
        candidatesEvaluatedToFinalBlock = finalBlockEvaluations + gridSize,
        inferredDeltaCandidates = 0L,
        adaptiveRefinements = refinements,
        maximumRefinementDepth = maximumRefinementDepth,
        assumptionStatus = "verified by full-grid fallback"
      )
    )
  }

  coarseIndices <- unique(as.integer(round(seq(
    1,
    gridSize,
    length.out = min(coarseGridSize, gridSize)
  ))))
  priorParameters <- saviDesign[["betaPriorParameterValues"]]
  posteriorThetaA <- (priorParameters[["betaA1"]] + sum(ya)) / (
    priorParameters[["betaA1"]] + priorParameters[["betaA2"]] +
      sum(preparedSearch[["na"]])
  )
  posteriorThetaB <- (priorParameters[["betaB1"]] + sum(yb)) / (
    priorParameters[["betaB1"]] + priorParameters[["betaB2"]] +
      sum(preparedSearch[["nb"]])
  )
  supportedIndex <- which.min(abs(
    deltaGrid - (posteriorThetaA - posteriorThetaB)
  ))
  coarseIndices <- sort(unique(c(coarseIndices, supportedIndex)))
  for (index in coarseIndices) evaluateIndex(index)

  survivingCoarse <- coarseIndices[
    crossingTimes[coarseIndices] == nSteps + 1L
  ]
  if (length(survivingCoarse) == 0L) {
    return(fallbackToFullGrid("coarse grid found no final retained candidate"))
  }
  centreIndex <- survivingCoarse[ceiling(length(survivingCoarse) / 2)]

  refineBranch <- function(knownIndices, direction) {
    if (length(knownIndices) < 2L) return(TRUE)
    knownTimes <- crossingTimes[knownIndices]
    monotone <- if (direction == "increasing") {
      all(diff(knownTimes) >= 0L)
    } else {
      all(diff(knownTimes) <= 0L)
    }
    if (!monotone) return(FALSE)

    stack <- lapply(seq_len(length(knownIndices) - 1L), function(position) {
      list(
        left = knownIndices[position],
        right = knownIndices[position + 1L],
        depth = 1L
      )
    })
    while (length(stack) > 0L) {
      interval <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      left <- interval[["left"]]
      right <- interval[["right"]]
      depth <- interval[["depth"]]
      maximumRefinementDepth <<- max(maximumRefinementDepth, depth)
      if (right - left <= 1L) next

      leftTime <- crossingTimes[left]
      rightTime <- crossingTimes[right]
      if (leftTime == rightTime) {
        crossingTimes[(left + 1L):(right - 1L)] <<- leftTime
        next
      }

      midpoint <- floor((left + right) / 2)
      lowerTime <- min(leftTime, rightTime)
      upperTime <- max(leftTime, rightTime)
      midpointTime <- evaluateIndex(
        midpoint,
        lowerCrossingTime = lowerTime,
        upperCrossingTime = upperTime
      )
      refinements <<- refinements + 1L
      if (!is.na(fallbackReason)) return(FALSE)

      stack[[length(stack) + 1L]] <- list(
        left = midpoint,
        right = right,
        depth = depth + 1L
      )
      stack[[length(stack) + 1L]] <- list(
        left = left,
        right = midpoint,
        depth = depth + 1L
      )
    }
    TRUE
  }

  leftKnown <- sort(unique(c(coarseIndices[coarseIndices < centreIndex], centreIndex)))
  rightKnown <- sort(unique(c(centreIndex, coarseIndices[coarseIndices > centreIndex])))
  if (!refineBranch(leftKnown, "increasing") ||
      !refineBranch(rightKnown, "decreasing") ||
      anyNA(crossingTimes)) {
    reason <- if (is.na(fallbackReason)) {
      "evaluated crossing times were not unimodal"
    } else {
      fallbackReason
    }
    return(fallbackToFullGrid(reason))
  }

  result <- linearDifferenceSequenceFromCrossingTimes(
    deltaGrid = deltaGrid,
    crossingTimes = crossingTimes,
    nSteps = nSteps
  )
  list(
    confidenceSequence = result,
    crossingTimes = crossingTimes,
    diagnostics = list(
      fallback = FALSE,
      fallbackReason = NA_character_,
      distinctDeltaCandidatesEvaluated = sum(evaluated),
      deltaCandidateEvaluations = candidateEvaluations,
      candidatePrefixReplays = candidateEvaluations,
      candidateBlockEvaluations = candidateBlockEvaluations,
      candidatesEvaluatedToFinalBlock = finalBlockEvaluations,
      inferredDeltaCandidates = sum(!evaluated),
      adaptiveRefinements = refinements,
      maximumRefinementDepth = maximumRefinementDepth,
      assumptionStatus = "assumed from monotonicity at evaluated grid indices"
    )
  )
}

#' Experimental adaptive confidence sequence for a linear difference
#'
#' Computes the same finite-grid running-intersection confidence sequence as
#' [computeConfidenceSequenceForDifferenceTwoProportions()] while adaptively
#' evaluating first threshold crossing times over the exact same candidate
#' grid. The method assumes crossing times are nondecreasing on the lower side
#' of a retained centre and nonincreasing on the upper side. It falls back to
#' the complete grid when evaluated candidates contradict that assumption or a
#' coarse retained centre cannot be found.
#'
#' Newly evaluated candidates must replay their prefix from block one because
#' the RIPr denominator depends nonlinearly on the candidate difference. This
#' experimental routine therefore saves candidate-block evaluations only when
#' monotone regions with equal endpoint crossing times can be inferred.
#'
#' @inheritParams computeConfidenceSequenceForDifferenceTwoProportions
#' @param coarseGridSize Number of candidate indices in the initial subset of
#'   the final grid.
#'
#' @return A data frame with `block`, `lowerBound`, and `upperBound`, matching
#'   [computeConfidenceSequenceForDifferenceTwoProportions()].
computeConfidenceSequenceForDifferenceTwoProportionsAdaptive <- function(
  ya,
  yb,
  confidenceBoundGridPrecision,
  saviDesign,
  coarseGridSize = 8L
) {
  computeConfidenceSequenceForDifferenceTwoProportionsAdaptiveDetails(
    ya = ya,
    yb = yb,
    confidenceBoundGridPrecision = confidenceBoundGridPrecision,
    saviDesign = saviDesign,
    coarseGridSize = coarseGridSize
  )[["confidenceSequence"]]
}

# Log-odds-ratio confidence sequence ----

# Solve for the thetaA coordinate of the reverse information projection onto
# the log-odds-ratio boundary logit(thetaA) - logit(thetaB) = logOddsRatio.
solveLogOddsRatioRIPrThetaA <- function(
  numeratorThetaA,
  numeratorThetaB,
  na,
  nb,
  logOddsRatio
) {
  nSteps <- length(numeratorThetaA)
  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)
  if (!all(c(length(numeratorThetaB), length(na), length(nb)) == nSteps)) {
    stop("Numerator theta vectors and group sizes must align by block.")
  }
  if (length(logOddsRatio) != 1L || !is.finite(logOddsRatio)) {
    stop("logOddsRatio must be a finite scalar.")
  }

  solveOne <- function(thetaStarA, thetaStarB, blockSizeA, blockSizeB) {
    targetMean <- blockSizeA * thetaStarA + blockSizeB * thetaStarB
    score <- function(thetaA) {
      thetaB <- stats::plogis(stats::qlogis(thetaA) - logOddsRatio)
      blockSizeA * thetaA + blockSizeB * thetaB - targetMean
    }

    # The score is strictly increasing from a negative to a positive value.
    stats::uniroot(
      score,
      interval = c(0, 1),
      tol = .Machine$double.eps^0.75
    )$root
  }

  mapply(
    solveOne,
    numeratorThetaA,
    numeratorThetaB,
    na,
    nb,
    USE.NAMES = FALSE
  )
}

#' Calculate e-processes over a grid of log odds ratios
#'
#' For every candidate log odds ratio, computes the Turner likelihood-ratio
#' process whose denominator is the reverse information projection onto the
#' corresponding one-sided null. The log odds ratio is
#' `logit(thetaA) - logit(thetaB)`. Lower-bound processes use the null
#' `logOddsRatio <= candidate`; upper-bound processes use the null
#' `logOddsRatio >= candidate`. Either direction accepts signed candidates.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param na,nb Number of observations in groups A and B in each data block.
#'   A scalar is recycled over all blocks.
#' @param priorParameters Named list with `betaA1`, `betaA2`, `betaB1`, and
#'   `betaB2` for the Turner predictor.
#' @param logOddsRatioGrid Finite candidate log odds ratios.
#' @param bound Whether to test the one-sided null `logOddsRatio <= candidate`
#'   for a lower bound or `logOddsRatio >= candidate` for an upper bound.
#'
#' @return A list containing `logOddsRatio` and `logEProcesses`. Rows of the
#'   log e-process matrix are data blocks and columns are candidate log odds
#'   ratios.
calculateEValuesForLogOddsRatioGrid <- function(
  ya,
  yb,
  na,
  nb,
  priorParameters,
  logOddsRatioGrid,
  bound = c("lower", "upper")
) {
  bound <- match.arg(bound)
  nSteps <- length(ya)
  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)
  if (!all(c(length(yb), length(na), length(nb)) == nSteps)) {
    stop("ya, yb, na, and nb must have the same length.")
  }
  if (length(logOddsRatioGrid) < 1L ||
      any(!is.finite(logOddsRatioGrid))) {
    stop("logOddsRatioGrid must contain at least one finite candidate.")
  }
  breveMean <- betaPredictiveMeansTwoProportions(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    priorParameters = priorParameters
  )
  predictorLogOddsRatio <- stats::qlogis(breveMean[["thetaA"]]) -
    stats::qlogis(breveMean[["thetaB"]])
  gridSize <- length(logOddsRatioGrid)
  logEProcesses <- matrix(
    NA_real_,
    nrow = nSteps,
    ncol = gridSize,
    dimnames = list(NULL, as.character(logOddsRatioGrid))
  )

  for (gridIndex in seq_along(logOddsRatioGrid)) {
    logOddsRatio <- logOddsRatioGrid[gridIndex]
    insideNull <- if (bound == "lower") {
      predictorLogOddsRatio <= logOddsRatio
    } else {
      predictorLogOddsRatio >= logOddsRatio
    }
    denominatorThetaA <- breveMean[["thetaA"]]
    denominatorThetaB <- breveMean[["thetaB"]]

    if (any(!insideNull)) {
      outsideNull <- !insideNull
      denominatorThetaA[outsideNull] <- solveLogOddsRatioRIPrThetaA(
        numeratorThetaA = breveMean[["thetaA"]][outsideNull],
        numeratorThetaB = breveMean[["thetaB"]][outsideNull],
        na = na[outsideNull],
        nb = nb[outsideNull],
        logOddsRatio = logOddsRatio
      )
      denominatorThetaB[outsideNull] <- stats::plogis(
        stats::qlogis(denominatorThetaA[outsideNull]) - logOddsRatio
      )
    }

    logEProcesses[, gridIndex] <- logLikelihoodRatioProcess(
      ya = ya,
      yb = yb,
      na = na,
      nb = nb,
      numeratorThetaA = breveMean[["thetaA"]],
      numeratorThetaB = breveMean[["thetaB"]],
      denominatorThetaA = denominatorThetaA,
      denominatorThetaB = denominatorThetaB
    )
  }

  list(logOddsRatio = logOddsRatioGrid, logEProcesses = logEProcesses)
}

#' Confidence sequence for the log odds ratio
#'
#' Constructs a symmetric candidate grid from the supplied resolution and
#' search bounds. Every candidate e-process is computed once over the full data
#' sequence in each direction. Both one-sided inversions use the full signed
#' grid, and each side uses alpha/2.
#' The log odds ratio is `logit(thetaA) - logit(thetaB)`.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param confidenceBoundGridPrecision Number of candidate values on each side
#'   of zero.
#' @param logOddsConfidenceSearchBounds Positive finite lower and upper bounds
#'   for the absolute candidate log odds ratios.
#' @param saviDesign A `saviDesign` returned by
#'   [designSaviTwoProportions()].
#'
#' @return A data frame with `block`, `lowerBound`, and `upperBound`. A bound is
#'   infinite until its one-sided sequence rejects a candidate and is `NA` if
#'   every candidate on that side has been rejected.
computeConfidenceSequenceForLogOddsRatioTwoProportions <- function(
  ya,
  yb,
  confidenceBoundGridPrecision,
  logOddsConfidenceSearchBounds,
  saviDesign
) {
  if (length(confidenceBoundGridPrecision) != 1L ||
      !is.finite(confidenceBoundGridPrecision) ||
      confidenceBoundGridPrecision < 1L ||
      confidenceBoundGridPrecision %% 1 != 0) {
    stop("confidenceBoundGridPrecision must be a positive integer.")
  }
  if (length(logOddsConfidenceSearchBounds) != 2L ||
      any(!is.finite(logOddsConfidenceSearchBounds)) ||
      logOddsConfidenceSearchBounds[1L] <= 0 ||
      logOddsConfidenceSearchBounds[1L] >=
        logOddsConfidenceSearchBounds[2L]) {
    stop(
      paste(
        "logOddsConfidenceSearchBounds must contain two positive,",
        "increasing finite values."
      )
    )
  }

  # use tanh to transform log uniform in log odds ratio
  positiveTransformedBounds <- tanh(logOddsConfidenceSearchBounds / 4)

  positiveGrid <- 4 * atanh(seq(
    positiveTransformedBounds[1L],
    positiveTransformedBounds[2L],
    length.out = confidenceBoundGridPrecision
  ))

  candidateGrid <- c(-rev(positiveGrid), 0, positiveGrid)

  lowerGridProcesses <- calculateEValuesForLogOddsRatioGrid(
    ya = ya,
    yb = yb,
    na = saviDesign[["nPlan"]][["na"]],
    nb = saviDesign[["nPlan"]][["nb"]],
    priorParameters = saviDesign[["betaPriorParameterValues"]],
    logOddsRatioGrid = candidateGrid,
    bound = "lower"
  )
  upperGridProcesses <- calculateEValuesForLogOddsRatioGrid(
    ya = ya,
    yb = yb,
    na = saviDesign[["nPlan"]][["na"]],
    nb = saviDesign[["nPlan"]][["nb"]],
    priorParameters = saviDesign[["betaPriorParameterValues"]],
    logOddsRatioGrid = candidateGrid,
    bound = "upper"
  )

  lowerLogEProcesses <- lowerGridProcesses[["logEProcesses"]]
  upperLogEProcesses <- upperGridProcesses[["logEProcesses"]]
  nSteps <- nrow(lowerLogEProcesses)
  logThreshold <- log(2 / saviDesign[["alpha"]])
  retainedLower <- matrix(
    FALSE,
    nrow = length(candidateGrid),
    ncol = nSteps
  )
  retainedUpper <- matrix(
    FALSE,
    nrow = length(candidateGrid),
    ncol = nSteps
  )

  for (gridIndex in seq_along(candidateGrid)) {
    retainedLower[gridIndex, ] <-
      cummax(lowerLogEProcesses[, gridIndex]) < logThreshold
    retainedUpper[gridIndex, ] <-
      cummax(upperLogEProcesses[, gridIndex]) < logThreshold
  }

  lowerBound <- rep(-Inf, nSteps)
  upperBound <- rep(Inf, nSteps)
  for (block in seq_len(nSteps)) {
    retainedCandidates <- candidateGrid[retainedLower[, block]]
    if (length(retainedCandidates) == 0L) {
      lowerBound[block] <- NA_real_
    } else if (!all(retainedLower[, block])) {
      lowerBound[block] <- min(retainedCandidates)
    }

    retainedCandidates <- candidateGrid[retainedUpper[, block]]
    if (length(retainedCandidates) == 0L) {
      upperBound[block] <- NA_real_
    } else if (!all(retainedUpper[, block])) {
      upperBound[block] <- max(retainedCandidates)
    }
  }

  data.frame(
    block = seq_len(nSteps),
    lowerBound = lowerBound,
    upperBound = upperBound
  )
}

# Stopping-time simulation ----

# Calculate cumulative log e-processes for one simulated data chunk. The input
# and output matrices have blocks as rows and active simulation paths as
# columns. Previous sufficient statistics have one value per active path;
# previousBlocks is shared because all active paths have reached the same time.
turnerLogEProcessChunk <- function(
  yaChunk,
  ybChunk,
  na,
  nb,
  previousYa,
  previousYb,
  previousBlocks,
  previousLogE,
  priorParameters
) {
  nBlocks <- nrow(yaChunk)
  nPaths <- ncol(yaChunk)
  if (!identical(dim(yaChunk), dim(ybChunk))) {
    stop("yaChunk and ybChunk must have the same dimensions.")
  }
  if (!all(c(length(previousYa), length(previousYb), length(previousLogE)) ==
      nPaths)) {
    stop("Previous path states must have one value per chunk column.")
  }

  # helper to calculate the cumulative sums in each column
  columnCumsums <- function(values) {
    matrix(
      apply(values, 2L, cumsum),
      nrow = nrow(values),
      ncol = ncol(values)
    )
  }

  cumulativeYaWithinChunk <- columnCumsums(yaChunk)
  cumulativeYbWithinChunk <- columnCumsums(ybChunk)
  laggedYaWithinChunk <- rbind(
    rep(0, nPaths),
    cumulativeYaWithinChunk[-nBlocks, , drop = FALSE]
  )
  laggedYbWithinChunk <- rbind(
    rep(0, nPaths),
    cumulativeYbWithinChunk[-nBlocks, , drop = FALSE]
  )
  predictableYa <- sweep(laggedYaWithinChunk, 2L, previousYa, "+")
  predictableYb <- sweep(laggedYbWithinChunk, 2L, previousYb, "+")
  predictableNa <- (previousBlocks + seq_len(nBlocks) - 1L) * na
  predictableNb <- (previousBlocks + seq_len(nBlocks) - 1L) * nb

  numeratorThetaA <- (priorParameters[["betaA1"]] + predictableYa) / (
    predictableNa + priorParameters[["betaA1"]] + priorParameters[["betaA2"]]
  )
  numeratorThetaB <- (priorParameters[["betaB1"]] + predictableYb) / (
    predictableNb + priorParameters[["betaB1"]] + priorParameters[["betaB2"]]
  )
  denominatorTheta <- (
    na * numeratorThetaA + nb * numeratorThetaB
  ) / (na + nb)

  logEIncrements <-
    stats::dbinom(yaChunk, na, numeratorThetaA, log = TRUE) +
    stats::dbinom(ybChunk, nb, numeratorThetaB, log = TRUE) -
    stats::dbinom(yaChunk, na, denominatorTheta, log = TRUE) -
    stats::dbinom(ybChunk, nb, denominatorTheta, log = TRUE)

  sweep(columnCumsums(logEIncrements), 2L, previousLogE, "+")
}

# Simulate paths for one fixed (thetaA, thetaB) pair. Each chunk has blocks as
# rows and active simulation paths as columns. Paths are removed after crossing
# so later chunks only simulate paths still below the threshold. The two output
# vectors have length nSimulations; a non-crossing path has Inf and NA.
simulateTurnerStoppingTimes <- function(
  na,
  nb,
  thetaA,
  thetaB,
  alpha = 0.05,
  priorParameters = NULL,
  nSimulations = 1e3,
  maxBlocks = 1e4,
  chunkSize = 50L
) {
  if (length(na) != 1L || length(nb) != 1L ||
      any(!is.finite(c(na, nb))) || any(c(na, nb) <= 0) ||
      any(c(na, nb) %% 1 != 0)) {
    stop("na and nb must be positive integer block sizes.")
  }
  if (length(thetaA) != 1L || length(thetaB) != 1L ||
      any(!is.finite(c(thetaA, thetaB))) ||
      any(c(thetaA, thetaB) < 0) || any(c(thetaA, thetaB) > 1)) {
    stop("thetaA and thetaB must be probabilities between 0 and 1.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }
  if (any(!is.finite(c(nSimulations, maxBlocks, chunkSize))) ||
      any(c(nSimulations, maxBlocks, chunkSize) < 1) ||
      any(c(nSimulations, maxBlocks, chunkSize) %% 1 != 0)) {
    stop("nSimulations, maxBlocks, and chunkSize must be positive integers.")
  }

  logThreshold <- log(1 / alpha)
  stoppingTimes <- rep(Inf, nSimulations)
  eValuesAtStopping <- rep(NA_real_, nSimulations)
  logEValues <- numeric(nSimulations)
  cumulativeYa <- numeric(nSimulations)
  cumulativeYb <- numeric(nSimulations)
  active <- seq_len(nSimulations)
  blocksCompleted <- 0L

  priorParameters <- resolveTurnerPriorParameters(
    na = na,
    nb = nb,
    priorParameters = priorParameters
  )

  while (length(active) > 0L && blocksCompleted < maxBlocks) {
    blocksInChunk <- min(chunkSize, maxBlocks - blocksCompleted)
    nActive <- length(active)
    yaChunk <- matrix(
      stats::rbinom(
        n = blocksInChunk * nActive,
        size = na,
        prob = thetaA
      ),
      nrow = blocksInChunk
    )
    ybChunk <- matrix(
      stats::rbinom(
        n = blocksInChunk * nActive,
        size = nb,
        prob = thetaB
      ),
      nrow = blocksInChunk
    )
    chunkLogE <- turnerLogEProcessChunk(
      yaChunk = yaChunk,
      ybChunk = ybChunk,
      na = na,
      nb = nb,
      previousYa = cumulativeYa[active],
      previousYb = cumulativeYb[active],
      previousBlocks = blocksCompleted,
      previousLogE = logEValues[active],
      priorParameters = priorParameters
    )
    crossingWithinChunk <- apply(
      chunkLogE >= logThreshold,
      2L,
      function(crossed) {
        crossing <- which(crossed)[1L]
        if (length(crossing) == 0L) NA_integer_ else crossing
      }
    )
    crossedColumns <- which(!is.na(crossingWithinChunk))

    if (length(crossedColumns) > 0L) {
      crossedPaths <- active[crossedColumns]
      crossingRows <- crossingWithinChunk[crossedColumns]
      stoppingTimes[crossedPaths] <- blocksCompleted + crossingRows
      eValuesAtStopping[crossedPaths] <- exp(
        chunkLogE[cbind(crossingRows, crossedColumns)]
      )
    }

    survivingColumns <- which(is.na(crossingWithinChunk))
    if (length(survivingColumns) > 0L) {
      survivingPaths <- active[survivingColumns]
      cumulativeYa[survivingPaths] <- cumulativeYa[survivingPaths] +
        colSums(yaChunk[, survivingColumns, drop = FALSE])
      cumulativeYb[survivingPaths] <- cumulativeYb[survivingPaths] +
        colSums(ybChunk[, survivingColumns, drop = FALSE])
      logEValues[survivingPaths] <-
        chunkLogE[blocksInChunk, survivingColumns]
      active <- survivingPaths
    } else {
      active <- integer()
    }
    blocksCompleted <- blocksCompleted + blocksInChunk
  }

  list(
    stoppingTimes = stoppingTimes,
    eValuesAtStopping = eValuesAtStopping
  )
}

# Simulate stopping times over a feasible baseline grid for one fixed effect.
# The returned stopping-time and e-value matrices have baseline-grid points as
# rows and independent simulation paths as columns.
simulateStoppingTimesForParameter <- function(
  na,
  nb,
  difference = NULL,
  logOddsRatio = NULL,
  alpha = 0.05,
  priorParameters = NULL,
  nSimulations = 1e3,
  maxBlocks = 1e4,
  baselineGridSize = 8L,
  chunkSize = 50L
) {
  if (is.null(difference) == is.null(logOddsRatio)) {
    stop("Supply exactly one of difference and logOddsRatio.")
  }
  if (!is.null(difference) &&
      (length(difference) != 1L || !is.finite(difference) ||
       difference == 0 || abs(difference) >= 1)) {
    stop("difference must be nonzero and strictly between -1 and 1.")
  }
  if (!is.null(logOddsRatio) &&
      (length(logOddsRatio) != 1L || !is.finite(logOddsRatio) ||
       logOddsRatio == 0)) {
    stop("logOddsRatio must be a nonzero finite scalar.")
  }
  if (length(baselineGridSize) != 1L || !is.finite(baselineGridSize) ||
      baselineGridSize < 1L || baselineGridSize %% 1 != 0) {
    stop("baselineGridSize must be a positive integer.")
  }

  rho <- seq_len(baselineGridSize) / (baselineGridSize + 1)
  if (!is.null(difference)) {
    thetaA <- rho * (1 - abs(difference)) + max(difference, 0)
    thetaB <- thetaA - difference
  } else {
    thetaA <- rho
    thetaB <- stats::plogis(stats::qlogis(thetaA) - logOddsRatio)
  }
  stoppingTimes <- matrix(
    Inf,
    nrow = baselineGridSize,
    ncol = nSimulations
  )
  eValuesAtStopping <- matrix(
    NA_real_,
    nrow = baselineGridSize,
    ncol = nSimulations
  )

  for (gridIndex in seq_len(baselineGridSize)) {
    simulation <- simulateTurnerStoppingTimes(
      na = na,
      nb = nb,
      thetaA = thetaA[gridIndex],
      thetaB = thetaB[gridIndex],
      alpha = alpha,
      priorParameters = priorParameters,
      nSimulations = nSimulations,
      maxBlocks = maxBlocks,
      chunkSize = chunkSize
    )
    stoppingTimes[gridIndex, ] <- simulation[["stoppingTimes"]]
    eValuesAtStopping[gridIndex, ] <- simulation[["eValuesAtStopping"]]
  }

  list(
    difference = difference,
    logOddsRatio = logOddsRatio,
    alpha = alpha,
    maxBlocks = maxBlocks,
    thetaA = thetaA,
    thetaB = thetaB,
    stoppingTimes = stoppingTimes,
    eValuesAtStopping = eValuesAtStopping
  )
}

#' Simulate worst-case stopping times for two proportions
#'
#' Simulates the Turner e-process under either a fixed linear difference or a
#' fixed log odds ratio over a grid of feasible baseline probabilities. At
#' every grid point, paths are generated and evaluated in fixed-size chunks. A
#' path that does not cross `1 / alpha` within `maxBlocks` has stopping time
#' `Inf` and stopping e-value `NA`.
#'
#' The baseline grid is a simulation device, not an alternative restriction.
#' It ranges over values of `thetaA` for which both probabilities are in
#' `(0, 1)`. The worst-case stopping time is the largest empirical
#' `(1 - beta)` quantile over this grid. Every bootstrap replicate resamples
#' all grid rows and recomputes that maximum, thereby including uncertainty
#' about which baseline probability is worst.
#'
#' @param na,nb Number of observations in groups A and B per block.
#' @param difference Data-generating difference `thetaA - thetaB`. Supply
#'   exactly one of `difference` and `logOddsRatio`.
#' @param logOddsRatio Data-generating value
#'   `logit(thetaA) - logit(thetaB)`. Supply exactly one of `difference` and
#'   `logOddsRatio`.
#' @param alpha E-process rejection threshold is `1 / alpha`.
#' @param beta Target type-II error used to select the stopping-time quantile.
#' @param priorParameters Optional Turner Beta prior parameters.
#' @param nSimulations Number of simulated paths at each baseline probability.
#' @param maxBlocks Maximum number of blocks simulated per path.
#' @param baselineGridSize Number of equally spaced feasible baseline
#'   probabilities.
#' @param nBoot Number of nonparametric bootstrap samples used to estimate the
#'   standard error of the worst-case stopping-time quantile; must be at least
#'   two.
#' @param chunkSize Maximum number of new blocks generated and evaluated at a
#'   time for each active path.
#'
#' @return A list with the worst-case baseline probabilities, stopping-time
#'   quantiles, and `baselineGridSize` by `nSimulations` matrices containing
#'   the per-path stopping results. Non-crossing paths have stopping time `Inf`
#'   and stopping e-value `NA`.
#'   `worstCaseStoppingTimeBootstrapSe` is the bootstrap standard error of the
#'   maximum grid quantile, and `worstCaseStoppingTimeTwoSe` is twice that
#'   standard error.
#' @export
simulateWorstCaseStoppingTimes <- function(
  na,
  nb,
  difference = NULL,
  logOddsRatio = NULL,
  alpha = 0.05,
  beta = 0.2,
  priorParameters = NULL,
  nSimulations = 1e3,
  maxBlocks = 1e4,
  baselineGridSize = 8L,
  nBoot = 1e3,
  chunkSize = 50L
) {
  if (length(beta) != 1L || !is.finite(beta) || beta <= 0 || beta >= 1) {
    stop("beta must be strictly between 0 and 1.")
  }
  if (length(nBoot) != 1L || !is.finite(nBoot) ||
      nBoot < 2L || nBoot %% 1 != 0) {
    stop("nBoot must be an integer of at least two.")
  }

  gridSimulation <- simulateStoppingTimesForParameter(
    na = na,
    nb = nb,
    difference = difference,
    logOddsRatio = logOddsRatio,
    alpha = alpha,
    priorParameters = priorParameters,
    nSimulations = nSimulations,
    maxBlocks = maxBlocks,
    baselineGridSize = baselineGridSize,
    chunkSize = chunkSize
  )

  stoppingTimes <- gridSimulation[["stoppingTimes"]]
  stoppingTimeQuantiles <- apply(
    stoppingTimes,
    1L,
    stats::quantile,
    probs = 1 - beta,
    names = FALSE
  )
  worstCaseIndex <- which.max(stoppingTimeQuantiles)
  bootstrapStoppingTimes <- replicate(nBoot, {
    max(vapply(seq_len(nrow(stoppingTimes)), function(gridIndex) {
      stats::quantile(
        sample(
          stoppingTimes[gridIndex, ],
          size = ncol(stoppingTimes),
          replace = TRUE
        ),
        probs = 1 - beta,
        names = FALSE
      )
    }, numeric(1)))
  })
  bootstrapSe <- if (any(!is.finite(bootstrapStoppingTimes))) {
    Inf
  } else {
    stats::sd(bootstrapStoppingTimes)
  }

  c(gridSimulation, list(
    beta = beta,
    stoppingTimeQuantiles = stoppingTimeQuantiles,
    worstCaseIndex = worstCaseIndex,
    worstCaseThetaA = gridSimulation[["thetaA"]][worstCaseIndex],
    worstCaseThetaB = gridSimulation[["thetaB"]][worstCaseIndex],
    worstCaseStoppingTime = stoppingTimeQuantiles[worstCaseIndex],
    worstCaseStoppingTimeBootstrapSe = bootstrapSe,
    worstCaseStoppingTimeTwoSe = 2 * bootstrapSe
  ))
}

# Simulate the lowest rejection probability over the feasible baseline grid at
# a fixed stopping horizon. Rows are baseline-grid points and columns are
# independent paths. The bootstrap repeats both the row-wise power estimates
# and the minimisation over the grid.
simulateWorstCasePower <- function(
  na,
  nb,
  difference = NULL,
  logOddsRatio = NULL,
  alpha = 0.05,
  priorParameters = NULL,
  nSimulations = 1e3,
  maxBlocks,
  baselineGridSize = 8L,
  nBoot = 1e3,
  chunkSize = 50L
) {
  if (!is.null(nBoot) &&
      (length(nBoot) != 1L || !is.finite(nBoot) ||
       nBoot < 2L || nBoot %% 1 != 0)) {
    stop("nBoot must be an integer of at least two.")
  }

  gridSimulation <- simulateStoppingTimesForParameter(
    na = na,
    nb = nb,
    difference = difference,
    logOddsRatio = logOddsRatio,
    alpha = alpha,
    priorParameters = priorParameters,
    nSimulations = nSimulations,
    maxBlocks = maxBlocks,
    baselineGridSize = baselineGridSize,
    chunkSize = chunkSize
  )
  rejected <- is.finite(gridSimulation[["stoppingTimes"]])
  power <- rowMeans(rejected)
  worstCaseIndex <- which.min(power)
  worstCasePowerTwoSe <- NULL
  if (!is.null(nBoot)) {
    bootstrapPower <- replicate(nBoot, {
      min(vapply(seq_len(nrow(rejected)), function(gridIndex) {
        mean(sample(
          rejected[gridIndex, ],
          size = ncol(rejected),
          replace = TRUE
        ))
      }, numeric(1)))
    })
    worstCasePowerTwoSe <- 2 * stats::sd(bootstrapPower)
  }

  c(gridSimulation[c("difference", "logOddsRatio", "thetaA", "thetaB")], list(
    power = power,
    worstCasePower = power[worstCaseIndex],
    worstCasePowerTwoSe = worstCasePowerTwoSe,
    worstCaseIndex = worstCaseIndex,
    worstCaseThetaA = gridSimulation[["thetaA"]][worstCaseIndex],
    worstCaseThetaB = gridSimulation[["thetaB"]][worstCaseIndex]
  ))
}

# Search a descending grid of positive effect sizes and return the smallest
# value whose worst-case power reaches 1 - beta.
simulateMinimumDetectableDifference <- function(
  na,
  nb,
  alpha,
  beta,
  effectMeasure = c("linearDifference", "logOddsRatio"),
  priorParameters = NULL,
  nSimulations = 1e3,
  maxBlocks,
  deltaGridSize = 10L,
  deltaMax = NULL,
  deltaMin = 0.01,
  baselineGridSize = 8L,
  chunkSize = 50L
) {
  effectMeasure <- match.arg(effectMeasure)
  if (length(beta) != 1L || !is.finite(beta) || beta <= 0 || beta >= 1) {
    stop("beta must be strictly between 0 and 1.")
  }
  if (length(deltaGridSize) != 1L || !is.finite(deltaGridSize) ||
      deltaGridSize < 1L || deltaGridSize %% 1 != 0) {
    stop("deltaGridSize must be a positive integer.")
  }
  if (is.null(deltaMax)) {
    deltaMax <- if (effectMeasure == "linearDifference") 0.99 else 5
  }
  if (any(!is.finite(c(deltaMin, deltaMax))) || deltaMin <= 0 ||
      deltaMin > deltaMax ||
      (effectMeasure == "linearDifference" && deltaMax >= 1)) {
    stop(paste(
      "deltaMin and deltaMax must be positive and ordered;",
      "linear differences must be less than one."
    ))
  }

  deltaGrid <- seq(deltaMax, deltaMin, length.out = deltaGridSize)
  effectArgument <- if (effectMeasure == "linearDifference") {
    "difference"
  } else {
    "logOddsRatio"
  }
  minimumDetectableEffect <- NULL
  for (delta in deltaGrid) {
    simulationArguments <- list(
      na = na,
      nb = nb,
      alpha = alpha,
      priorParameters = priorParameters,
      nSimulations = nSimulations,
      maxBlocks = maxBlocks,
      baselineGridSize = baselineGridSize,
      nBoot = NULL,
      chunkSize = chunkSize
    )
    simulationArguments[[effectArgument]] <- delta
    simulationResult <- do.call(simulateWorstCasePower, simulationArguments)
    if (simulationResult[["worstCasePower"]] < 1 - beta) {
      break
    }
    minimumDetectableEffect <- delta
  }
  if (is.null(minimumDetectableEffect)) {
    stop("No candidate effect reached the requested power at this horizon.")
  }

  minimumDetectableEffect
}

# Design ----

#' Designs a Savi Experiment to Test Two Proportions in Stream Data
#'
#' The Turner design uses independent Beta predictive distributions for the two
#' Bernoulli streams. The data-generating effect `delta` is either
#' `thetaA - thetaB` or `logit(thetaA) - logit(thetaB)`; it is not supplied to
#' the unrestricted Turner e-process.
#'
#' Supply `delta` and `beta` to estimate the required number of data blocks.
#' Supply `nBlocksPlan` and `delta` to estimate the worst-case type-II error at
#' that horizon, supply `nBlocksPlan` and `beta` to search a fixed grid for the
#' minimum detectable effect, or supply only `nBlocksPlan` to construct a pilot
#' design. Restricted alternatives are not yet implemented.
#'
#' @param na Number of observations in group A per data block.
#' @param nb Number of observations in group B per data block.
#' @param nBlocksPlan Planned number of data blocks collected.
#' @param beta Numeric in `(0, 1)` specifying the tolerable type II error.
#' @param delta Minimal relevant effect used to generate data in the design
#'   simulation.
#' @param effectMeasure Whether `delta` is a linear difference or log odds
#'   ratio.
#' @param alpha Numeric in `(0, 1)` specifying the tolerable type I error.
#' @param pilot Logical specifying whether this is a pilot design.
#' @param hyperParameterValues Named list containing positive `betaA1`,
#'   `betaA2`, `betaB1`, and `betaB2` values for the two Beta priors.
#' @param previousSaviTestResult Optional previous test result whose
#'   `posteriorHyperParameters` are used as the new prior parameters.
#' @param M Number of simulations at each baseline probability used to estimate
#'   `nBlocksPlan`, `beta`, or the minimum detectable `delta`.
#'
#' @return A `saviDesign` object. Its `nPlan` component contains `na`, `nb`, and
#'   the planned number of blocks; `nPlanTwoSe` contains twice the bootstrap
#'   standard error of these values.
#' @export
#'
#' @examples
#' designSaviTwoProportions(na = 1, nb = 1, nBlocksPlan = 20)
designSaviTwoProportions <- function(
  na,
  nb,
  nBlocksPlan = NULL,
  beta = NULL,
  delta = NULL,
  alpha = 0.05,
  pilot = FALSE,
  hyperParameterValues = NULL,
  previousSaviTestResult = NULL,
  M = 1e3,
  effectMeasure = c("linearDifference", "logOddsRatio")
) {
  effectMeasure <- match.arg(effectMeasure)
  # Validate inputs
  if (length(na) != 1L || length(nb) != 1L ||
      any(!is.finite(c(na, nb))) || any(c(na, nb) <= 0) ||
      any(c(na, nb) %% 1 != 0)) {
    stop("na and nb must be positive integer block sizes.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }
  if (length(M) != 1L || !is.finite(M) || M < 1L || M %% 1 != 0) {
    stop("M must be a positive integer.")
  }

  # Setup Beta prior
  note <- NULL
  if (!is.null(previousSaviTestResult)) {
    hyperParameterValues <- previousSaviTestResult[["posteriorHyperParameters"]]
    if (is.null(hyperParameterValues)) {
      stop("previousSaviTestResult has no posteriorHyperParameters.")
    }
    priorValuesForPrint <- paste(unlist(hyperParameterValues), collapse = " ")
    note <- c(
      note,
      "Hyperparameters set according to posterior values from previous test result"
    )
  } else if (is.null(hyperParameterValues)) {
    hyperParameterValues <- list(
      betaA1 = 0.18,
      betaB1 = (nb / na) * 0.18,
      betaA2 = 0.18,
      betaB2 = (nb / na) * 0.18
    )
    priorValuesForPrint <- "standard, REGRET optimal"
    note <- c(
      note,
      "Optimality of hyperparameters only verified for equal group sizes (na = nb = 1)"
    )
  } else {
    priorValuesForPrint <- paste(unlist(hyperParameterValues), collapse = " ")
  }

  hyperParameterValues <- resolveTurnerPriorParameters(
    na = na,
    nb = nb,
    priorParameters = hyperParameterValues
  )

  names(priorValuesForPrint) <- "Beta hyperparameters"
  nPlanTwoSe <- betaTwoSe <- NULL
  logImpliedTarget <- logImpliedTargetTwoSe <- NULL

  # Classify and validate the requested design case
  hasNBlocksPlan <- !is.null(nBlocksPlan)
  hasEffect <- !is.null(delta)
  hasBeta <- !is.null(beta)
  if (hasNBlocksPlan &&
      (length(nBlocksPlan) != 1L || !is.finite(nBlocksPlan) ||
       nBlocksPlan < 1L || nBlocksPlan %% 1 != 0)) {
    stop("nBlocksPlan must be a positive integer.")
  }
  if (hasEffect &&
      (length(delta) != 1L || !is.finite(delta) || delta == 0 ||
       (effectMeasure == "linearDifference" && abs(delta) >= 1))) {
    stop("delta must be finite and nonzero; linear differences must be between -1 and 1.")
  }
  effectArgument <- if (effectMeasure == "linearDifference") {
    "difference"
  } else {
    "logOddsRatio"
  }
  if (hasBeta &&
      (length(beta) != 1L || !is.finite(beta) || beta <= 0 || beta >= 1)) {
    stop("beta must be strictly between 0 and 1.")
  }

  # Case 1: nBlocksPlan only -> pilot design
  if (hasNBlocksPlan && !hasEffect && !hasBeta) {
    pilot <- TRUE
    warning("no simulation done!")


  # Case 2: delta and beta -> worst-case nBlocksPlan
  # bootstrapping to calculate the standard deviation
  } else if (!hasNBlocksPlan && hasEffect && hasBeta) {
    simulationArguments <- list(
      na = na,
      nb = nb,
      priorParameters = hyperParameterValues,
      alpha = alpha,
      beta = beta,
      nSimulations = M
    )
    simulationArguments[[effectArgument]] <- delta
    simulationResult <- do.call(
      simulateWorstCaseStoppingTimes,
      simulationArguments
    )
    nBlocksPlan <- simulationResult[["worstCaseStoppingTime"]]
    nPlanTwoSe <- c(0, 0, simulationResult[["worstCaseStoppingTimeTwoSe"]])


  # Case 3: nBlocksPlan and delta -> worst-case power
  } else if (hasNBlocksPlan && hasEffect && !hasBeta) {
    simulationArguments <- list(
      na = na,
      nb = nb,
      alpha = alpha,
      priorParameters = hyperParameterValues,
      nSimulations = M,
      maxBlocks = nBlocksPlan
    )
    simulationArguments[[effectArgument]] <- delta
    simulationResult <- do.call(simulateWorstCasePower, simulationArguments)
    beta <- 1 - simulationResult[["worstCasePower"]]
    betaTwoSe <- simulationResult[["worstCasePowerTwoSe"]]
    note <- c(
      note,
      "Implied target is not calculated for the stop-on-crossing simulator."
    )

  # Case 4: nBlocksPlan and beta -> minimum detectable delta
  } else if (hasNBlocksPlan && !hasEffect && hasBeta) {
    delta <- simulateMinimumDetectableDifference(
      na = na,
      nb = nb,
      alpha = alpha,
      beta = beta,
      effectMeasure = effectMeasure,
      priorParameters = hyperParameterValues,
      nSimulations = M,
      maxBlocks = nBlocksPlan
    )
    note <- c(
      note,
      "Effect selected from a descending simulation grid."
    )

  # Unsupported combinations
  } else {
    stop(
      paste(
        "Provide delta and beta to estimate nBlocksPlan; nBlocksPlan and",
        "delta to estimate beta; nBlocksPlan and beta to estimate delta;",
        "or only nBlocksPlan for a pilot design."
      )
    )
  }

  # Assemble the saviDesign result
  nPlan <- c(na, nb, nBlocksPlan)
  names(nPlan) <- c("na", "nb", "nBlocksPlan")

  if (!is.null(delta)) {
    names(delta) <- if (effectMeasure == "linearDifference") {
      "difference"
    } else {
      "logOddsRatio"
    }
  }

  result <- list(
    "nPlan" = nPlan,
    "nPlanTwoSe" = nPlanTwoSe,
    "parameter" = priorValuesForPrint,
    "betaPriorParameterValues" = hyperParameterValues,
    "alpha" = alpha,
    "beta" = beta,
    "betaTwoSe" = betaTwoSe,
    "logImpliedTarget" = logImpliedTarget,
    "logImpliedTargetTwoSe" = logImpliedTargetTwoSe,
    "esMin" = delta,
    "effectMeasure" = effectMeasure,
    "h0" = 0,
    "testType" = "2x2",
    "testName" = "Two Proportions",
    "alternativeRestriction" = "none",
    "alternative" = "twoSided",
    "pilot" = pilot,
    "lowN" = NULL,
    "highN" = NULL,
    "call" = sys.call(),
    "timeStamp" = Sys.time(),
    "note" = note
  )
  class(result) <- "saviDesign"

  return(result)
}


# Calculation ----

# saviTwoProportionsTest <- function(ya, yb, designObj = NULL, wantConfidenceSequence = FALSE, ciValue = NULL,
# confidenceBoundGridPrecision = 20, logOddsConfidenceSearchBounds = c(0.01, 5), pilot = FALSE)
