# Turner's method and design ----

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
#' Computes the Turner e-process for `thetaB - thetaA = 0`. This is the
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

  if (is.null(priorParameters)) {
    # Default to REGRET optimal based on first data block sample sizes
    priorParameters <- list(
      betaA1 = 0.18,
      betaA2 = 0.18,
      betaB1 = (nb[1] / na[1]) * 0.18,
      betaB2 = (nb[1] / na[1]) * 0.18
    )
  }

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
# the null curve thetaB - thetaA = difference at every data block.
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
    coefficients <- c(
      -A * thetaStarA * difference * (1 - difference),
      A * (
        difference * (1 - difference) -
          thetaStarA * (1 - 2 * difference)
      ) + B * (difference - thetaStarB),
      A * (1 - 2 * difference + thetaStarA) +
        B * (1 - difference + thetaStarB),
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
#' For every candidate `delta = thetaB - thetaA`, computes the Turner
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
      denominatorThetaB = denominatorThetaA + delta
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
#' @param gridSize Number of candidate differences used for the grid
#'   approximation in the feasible open interval `(-1, 1)`.
#' @param saviDesign A `saviDesign` returned by
#'   [designSaviTwoProportions()].
#'
#' @return A data frame with `block`, `lowerBound`, and `upperBound`. If no
#'   grid point remains at a block, both bounds are `NA`.
computeConfidenceSequenceForDifferenceTwoProportions <- function(
  ya,
  yb,
  gridSize,
  saviDesign
) {
  gridProcesses <- calculateEValuesForLinearDeltaGrid(
    ya = ya,
    yb = yb,
    na = saviDesign[["nPlan"]][["na"]],
    nb = saviDesign[["nPlan"]][["nb"]],
    priorParameters = saviDesign[["betaPriorParameterValues"]],
    gridSize = gridSize
  )
  deltaGrid <- gridProcesses[["delta"]]
  logEProcesses <- gridProcesses[["logEProcesses"]]
  nSteps <- nrow(logEProcesses)
  logThreshold <- log(1 / saviDesign[["alpha"]])
  retained <- matrix(FALSE, nrow = gridSize, ncol = nSteps)

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

# Prepare the data shared by candidate linear-difference e-processes. The
# numerator probabilities and log likelihood do not depend on the candidate
# difference, so they are calculated only once for adaptive searches.
prepareLinearDifferenceEProcess <- function(
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
  if (any(!is.finite(c(ya, yb, na, nb))) ||
      any(c(ya, yb, na, nb) %% 1 != 0) ||
      any(na <= 0) || any(nb <= 0) ||
      any(ya < 0) || any(yb < 0) ||
      any(ya > na) || any(yb > nb)) {
    stop("Success counts and group sizes must be valid integer-valued blocks.")
  }

  priorNames <- c("betaA1", "betaA2", "betaB1", "betaB2")
  if (!all(priorNames %in% names(priorParameters)) ||
      any(!is.finite(unlist(priorParameters[priorNames]))) ||
      any(unlist(priorParameters[priorNames]) <= 0)) {
    stop("priorParameters must contain positive betaA1, betaA2, betaB1, and betaB2 values.")
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
    priorParameters = priorParameters,
    breveMean = breveMean,
    logNumerator = logNumerator,
    nSteps = nSteps
  )
}

# Evaluate one candidate from block one, in chunks, until its running
# intersection excludes it or all blocks have been used. A chunk may evaluate
# up to chunkSize - 1 blocks beyond the crossing in exchange for vectorized
# RIPr and likelihood calculations.
evaluateLinearDifferenceUntilThreshold <- function(
  preparedProcess,
  difference,
  logThreshold,
  chunkSize = 50L
) {
  if (length(difference) != 1L || !is.finite(difference) ||
      difference <= -1 || difference >= 1) {
    stop("difference must be a finite scalar in the open interval (-1, 1).")
  }
  if (length(logThreshold) != 1L || !is.finite(logThreshold)) {
    stop("logThreshold must be a finite scalar.")
  }
  if (length(chunkSize) != 1L || !is.finite(chunkSize) ||
      chunkSize < 1L || chunkSize %% 1 != 0) {
    stop("chunkSize must be a positive integer.")
  }

  currentLogE <- 0
  maxLogE <- -Inf
  startBlock <- 1L
  evaluatedBlocks <- 0L
  nSteps <- preparedProcess[["nSteps"]]

  while (startBlock <= nSteps) {
    endBlock <- min(nSteps, startBlock + chunkSize - 1L)
    blockIndices <- startBlock:endBlock
    denominatorThetaA <- solveLinearDifferenceRIPrThetaA(
      numeratorThetaA = preparedProcess[["breveMean"]][["thetaA"]][blockIndices],
      numeratorThetaB = preparedProcess[["breveMean"]][["thetaB"]][blockIndices],
      na = preparedProcess[["na"]][blockIndices],
      nb = preparedProcess[["nb"]][blockIndices],
      difference = difference
    )
    logDenominator <- stats::dbinom(
      preparedProcess[["ya"]][blockIndices],
      preparedProcess[["na"]][blockIndices],
      denominatorThetaA,
      log = TRUE
    ) + stats::dbinom(
      preparedProcess[["yb"]][blockIndices],
      preparedProcess[["nb"]][blockIndices],
      denominatorThetaA + difference,
      log = TRUE
    )
    logEChunk <- currentLogE + cumsum(
      preparedProcess[["logNumerator"]][blockIndices] - logDenominator
    )
    maxLogE <- max(maxLogE, logEChunk)
    crossingInChunk <- match(TRUE, logEChunk >= logThreshold, nomatch = 0L)
    evaluatedBlocks <- endBlock

    if (crossingInChunk != 0L) {
      crossingBlock <- startBlock + crossingInChunk - 1L
      return(list(
        difference = difference,
        crossed = TRUE,
        crossingBlock = crossingBlock,
        logEAtExit = logEChunk[crossingInChunk],
        maxLogE = maxLogE,
        blocksEvaluated = evaluatedBlocks
      ))
    }

    currentLogE <- logEChunk[length(logEChunk)]
    startBlock <- endBlock + 1L
  }

  list(
    difference = difference,
    crossed = FALSE,
    crossingBlock = NA_integer_,
    logEAtExit = currentLogE,
    maxLogE = maxLogE,
    blocksEvaluated = evaluatedBlocks
  )
}

#' Heuristic adaptive final interval for a difference between two proportions
#'
#' Computes a final, running-intersection interval for `thetaB - thetaA` using
#' chunked e-process evaluations and adaptive endpoint bisection. Each queried
#' difference is evaluated from the first data block, so newly queried values
#' use their complete evidence history. This is a numerical heuristic: it
#' assumes the final retained set is one connected interval and does not prove
#' that property.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param precision Maximum width of each retained/rejected endpoint bracket.
#' @param saviDesign A `saviDesign` returned by
#'   [designSaviTwoProportions()].
#'
#' @return A list with conservative `lowerBound` and `upperBound`, the inner
#'   retained/rejected endpoint brackets, `achievedPrecision`, `status`, and
#'   diagnostics for every evaluated candidate. `heuristic` is always `TRUE`:
#'   the connected-interval assumption has not been established for this
#'   e-process.
#' @export
computeAdaptiveFinalIntervalForDifferenceTwoProportions <- function(
  ya,
  yb,
  precision = 1e-8,
  saviDesign
) {
  if (length(precision) != 1L || !is.finite(precision) ||
      precision <= 0 || precision >= 2) {
    stop("precision must be a finite scalar strictly between zero and two.")
  }
  if (is.null(saviDesign[["nPlan"]]) || is.null(saviDesign[["alpha"]]) ||
      is.null(saviDesign[["betaPriorParameterValues"]])) {
    stop("saviDesign must contain nPlan, alpha, and betaPriorParameterValues.")
  }
  if (length(saviDesign[["alpha"]]) != 1L ||
      !is.finite(saviDesign[["alpha"]]) ||
      saviDesign[["alpha"]] <= 0 || saviDesign[["alpha"]] >= 1) {
    stop("saviDesign alpha must be strictly between zero and one.")
  }

  preparedProcess <- prepareLinearDifferenceEProcess(
    ya = ya,
    yb = yb,
    na = saviDesign[["nPlan"]][["na"]],
    nb = saviDesign[["nPlan"]][["nb"]],
    priorParameters = saviDesign[["betaPriorParameterValues"]]
  )
  logThreshold <- log(1 / saviDesign[["alpha"]])
  domainMargin <- max(.Machine$double.eps * 8, min(1e-10, precision / 4))
  domainLower <- -1 + domainMargin
  domainUpper <- 1 - domainMargin
  candidateStates <- new.env(parent = emptyenv())

  evaluateCandidate <- function(difference) {
    difference <- min(domainUpper, max(domainLower, difference))
    key <- sprintf("%.17g", difference)
    if (!exists(key, envir = candidateStates, inherits = FALSE)) {
      assign(
        key,
        evaluateLinearDifferenceUntilThreshold(
          preparedProcess = preparedProcess,
          difference = difference,
          logThreshold = logThreshold
        ),
        envir = candidateStates
      )
    }
    get(key, envir = candidateStates, inherits = FALSE)
  }

  candidateDiagnostics <- function() {
    keys <- ls(candidateStates, all.names = TRUE)
    if (length(keys) == 0L) {
      return(data.frame(
        delta = numeric(),
        crossed = logical(),
        crossingBlock = integer(),
        logEAtExit = numeric(),
        maxLogE = numeric(),
        blocksEvaluated = integer()
      ))
    }
    states <- lapply(keys, get, envir = candidateStates, inherits = FALSE)
    diagnostics <- do.call(rbind, lapply(states, function(state) {
      data.frame(
        delta = state[["difference"]],
        crossed = state[["crossed"]],
        crossingBlock = state[["crossingBlock"]],
        logEAtExit = state[["logEAtExit"]],
        maxLogE = state[["maxLogE"]],
        blocksEvaluated = state[["blocksEvaluated"]]
      )
    }))
    diagnostics[order(diagnostics[["delta"]]), , drop = FALSE]
  }

  shapeWarning <- function() {
    diagnostics <- candidateDiagnostics()
    if (nrow(diagnostics) == 0L) return(FALSE)
    retainedRuns <- rle(!diagnostics[["crossed"]])
    sum(retainedRuns[["values"]]) > 1L
  }

  emptyResult <- function(status) {
    list(
      lowerBound = NA_real_,
      upperBound = NA_real_,
      lowerBracket = c(rejected = NA_real_, retained = NA_real_),
      upperBracket = c(retained = NA_real_, rejected = NA_real_),
      achievedPrecision = NA_real_,
      status = status,
      heuristic = TRUE,
      shapeWarning = shapeWarning(),
      candidateDiagnostics = candidateDiagnostics()
    )
  }

  priorParameters <- preparedProcess[["priorParameters"]]
  posteriorThetaA <- (priorParameters[["betaA1"]] + sum(preparedProcess[["ya"]])) / (
    sum(preparedProcess[["na"]]) + priorParameters[["betaA1"]] +
      priorParameters[["betaA2"]]
  )
  posteriorThetaB <- (priorParameters[["betaB1"]] + sum(preparedProcess[["yb"]])) / (
    sum(preparedProcess[["nb"]]) + priorParameters[["betaB1"]] +
      priorParameters[["betaB2"]]
  )
  centreDifference <- min(
    domainUpper,
    max(domainLower, posteriorThetaB - posteriorThetaA)
  )
  centreState <- evaluateCandidate(centreDifference)

  if (centreState[["crossed"]]) {
    pilotGrid <- seq(domainLower, domainUpper, length.out = 21L)
    pilotStates <- lapply(pilotGrid, evaluateCandidate)
    retainedPilot <- vapply(
      pilotStates,
      function(state) !state[["crossed"]],
      logical(1)
    )
    if (!any(retainedPilot)) {
      return(emptyResult("empty_or_undetermined"))
    }
    retainedDelta <- pilotGrid[retainedPilot]
    centreDifference <- retainedDelta[which.min(abs(retainedDelta - centreDifference))]
    centreState <- evaluateCandidate(centreDifference)
  }

  findEndpointBracket <- function(direction) {
    retainedDifference <- centreDifference
    step <- 0.05

    repeat {
      candidateDifference <- min(
        domainUpper,
        max(domainLower, centreDifference + direction * step)
      )
      if (identical(candidateDifference, retainedDifference)) {
        return(c(retained = retainedDifference, rejected = NA_real_))
      }
      candidateState <- evaluateCandidate(candidateDifference)
      if (candidateState[["crossed"]]) {
        return(c(retained = retainedDifference, rejected = candidateDifference))
      }
      retainedDifference <- candidateDifference
      if (candidateDifference %in% c(domainLower, domainUpper)) {
        return(c(retained = retainedDifference, rejected = NA_real_))
      }
      step <- step * 2
    }
  }

  refineEndpointBracket <- function(bracket) {
    if (is.na(bracket[["rejected"]])) return(bracket)

    while (abs(bracket[["retained"]] - bracket[["rejected"]]) > precision) {
      midpoint <- mean(bracket)
      midpointState <- evaluateCandidate(midpoint)
      if (midpointState[["crossed"]]) {
        bracket[["rejected"]] <- midpoint
      } else {
        bracket[["retained"]] <- midpoint
      }
    }
    bracket
  }

  lowerBracket <- refineEndpointBracket(findEndpointBracket(direction = -1))
  upperBracket <- refineEndpointBracket(findEndpointBracket(direction = 1))
  checkOuterNeighbour <- function(bracket, direction) {
    if (is.na(bracket[["rejected"]])) return(FALSE)

    bracketWidth <- abs(bracket[["retained"]] - bracket[["rejected"]])
    neighbourDifference <- min(
      domainUpper,
      max(
        domainLower,
        bracket[["rejected"]] + direction * bracketWidth
      )
    )
    if (identical(neighbourDifference, bracket[["rejected"]])) return(FALSE)

    !evaluateCandidate(neighbourDifference)[["crossed"]]
  }
  lowerNeighbourWarning <- checkOuterNeighbour(lowerBracket, direction = -1)
  upperNeighbourWarning <- checkOuterNeighbour(upperBracket, direction = 1)
  lowerWidth <- abs(lowerBracket[["retained"]] - lowerBracket[["rejected"]])
  upperWidth <- abs(upperBracket[["retained"]] - upperBracket[["rejected"]])
  bracketWidths <- c(lowerWidth, upperWidth)
  achievedPrecision <- if (all(is.na(bracketWidths))) {
    NA_real_
  } else {
    max(bracketWidths, na.rm = TRUE)
  }
  hasDomainBoundary <- any(is.na(c(
    lowerBracket[["rejected"]],
    upperBracket[["rejected"]]
  )))

  list(
    lowerBound = if (is.na(lowerBracket[["rejected"]])) -1 else lowerBracket[["rejected"]],
    upperBound = if (is.na(upperBracket[["rejected"]])) 1 else upperBracket[["rejected"]],
    lowerBracket = c(
      rejected = lowerBracket[["rejected"]],
      retained = lowerBracket[["retained"]]
    ),
    upperBracket = c(
      retained = upperBracket[["retained"]],
      rejected = upperBracket[["rejected"]]
    ),
    achievedPrecision = achievedPrecision,
    status = if (hasDomainBoundary) "heuristic_domain_boundary" else "heuristic",
    heuristic = TRUE,
    shapeWarning = shapeWarning() || lowerNeighbourWarning || upperNeighbourWarning,
    candidateDiagnostics = candidateDiagnostics()
  )
}

# Stopping Time Simulation ----

# Simulate complete Turner paths at fixed Bernoulli probabilities. A path that
# does not cross has no stopping e-value, so its time is Inf, its stopping
# e-value is NA, and its crossed indicator is FALSE.
simulateTurnerStoppingTimes <- function(
  na,
  nb,
  thetaA,
  thetaB,
  alpha = 0.05,
  priorParameters = NULL,
  nSimulations = 1e3,
  maxBlocks = 1e4
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
  if (any(!is.finite(c(nSimulations, maxBlocks))) ||
      any(c(nSimulations, maxBlocks) < 1) ||
      any(c(nSimulations, maxBlocks) %% 1 != 0)) {
    stop("nSimulations and maxBlocks must be positive integers.")
  }

  logThreshold <- log(1 / alpha)
  stoppingTimes <- rep(Inf, nSimulations)
  eValuesAtStopping <- rep(NA_real_, nSimulations)
  crossed <- rep(FALSE, nSimulations)

  for (simulation in seq_len(nSimulations)) {
    ya <- stats::rbinom(maxBlocks, size = na, prob = thetaA)
    yb <- stats::rbinom(maxBlocks, size = nb, prob = thetaB)
    logEProcess <- turnerEProcess(
      ya,
      yb,
      na,
      nb,
      priorParameters = priorParameters,
      log = TRUE
    )
    crossing <- match(TRUE, logEProcess >= logThreshold, nomatch = 0L)

    if (crossing != 0L) {
      stoppingTimes[simulation] <- crossing
      eValuesAtStopping[simulation] <- exp(logEProcess[crossing])
      crossed[simulation] <- TRUE
    }
  }

  list(
    stoppingTimes = stoppingTimes,
    eValuesAtStopping = eValuesAtStopping,
    crossed = crossed,
    blocksSimulated = rep(maxBlocks, nSimulations)
  )
}

#' Simulate the worst-case stopping time for a linear difference
#'
#' Simulates the Turner e-process under `thetaB - thetaA = difference` over a
#' grid of feasible baseline probabilities. A path that does not cross
#' `1 / alpha` within `maxBlocks` has stopping time `Inf`.
#'
#' The baseline grid is a simulation device, not an alternative restriction.
#' It ranges over values of `thetaA` for which both probabilities are in
#' `[0, 1]`.
#'
#' @param na,nb Number of observations in groups A and B per block.
#' @param difference Data-generating difference `thetaB - thetaA`.
#' @param alpha E-process rejection threshold is `1 / alpha`.
#' @param beta Target type-II error used to select the stopping-time quantile.
#' @param priorParameters Optional Turner Beta prior parameters.
#' @param nSimulations Number of simulated paths at each baseline probability.
#' @param maxBlocks Maximum number of blocks simulated per path.
#' @param gridSize Number of equally spaced feasible baseline probabilities.
#' @param nBoot Number of nonparametric bootstrap samples used to estimate the
#'   standard error of the worst-case stopping-time quantile; must be at least
#'   two.
#'
#' @return A list with the worst-case baseline probabilities, stopping-time
#'   quantiles, and per-path stopping results. `crossed` records whether a path
#'   crossed the threshold. Non-crossing paths have stopping time `Inf` and
#'   stopping e-value `NA` in `eValuesAtStopping`.
#'   `worstCaseStoppingTimeBootstrapSe` is the bootstrap standard error of the
#'   selected worst-case quantile, and `worstCaseStoppingTimeTwoSe` is twice
#'   that standard error.
#' @export
simulateWorstCaseStoppingTimeLinearDifference <- function(
  na,
  nb,
  difference,
  alpha = 0.05,
  beta = 0.2,
  priorParameters = NULL,
  nSimulations = 1e3,
  maxBlocks = 1e4,
  gridSize = 8,
  nBoot = 1e3
) {
  if (length(na) != 1L || length(nb) != 1L ||
      any(!is.finite(c(na, nb))) || any(c(na, nb) <= 0) ||
      any(c(na, nb) %% 1 != 0)) {
    stop("na and nb must be positive integer block sizes.")
  }
  if (length(difference) != 1L || !is.finite(difference) ||
      difference == 0 || abs(difference) >= 1) {
    stop("difference must be nonzero and strictly between -1 and 1.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(beta) != 1L || !is.finite(beta) || beta <= 0 || beta >= 1) {
    stop("alpha and beta must be strictly between 0 and 1.")
  }
  if (any(!is.finite(c(nSimulations, maxBlocks, gridSize, nBoot))) ||
      any(c(nSimulations, maxBlocks, gridSize) < 1) || nBoot < 2 ||
      any(c(nSimulations, maxBlocks, gridSize, nBoot) %% 1 != 0)) {
    stop(
      paste(
        "nSimulations, maxBlocks, and gridSize must be positive integers;",
        "nBoot must be an integer of at least two."
      )
    )
  }

  rhoGrid <- seq(1 / gridSize, 1 - 1 / gridSize, length.out = gridSize)
  thetaAValues <- rhoGrid * (1 - abs(difference)) -
    ifelse(difference < 0, difference, 0)

  simulationResults <- lapply(thetaAValues, function(thetaA) {
    simulateTurnerStoppingTimes(
      na = na,
      nb = nb,
      thetaA = thetaA,
      thetaB = thetaA + difference,
      alpha = alpha,
      priorParameters = priorParameters,
      nSimulations = nSimulations,
      maxBlocks = maxBlocks
    )
  })
  stoppingTimes <- lapply(simulationResults, `[[`, "stoppingTimes")
  stoppingTimeQuantiles <- vapply(
    stoppingTimes,
    stats::quantile,
    numeric(1),
    probs = 1 - beta,
    names = FALSE
  )
  empiricalPower <- vapply(simulationResults, function(result) {
    mean(result[["crossed"]])
  }, numeric(1))
  worstCaseIndex <- which.max(stoppingTimeQuantiles)
  worstCaseBootstrap <- computeBootObj(
    values = stoppingTimes[[worstCaseIndex]],
    beta = beta,
    nBoot = nBoot,
    objType = "nPlan"
  )
  bootstrapQuantiles <- worstCaseBootstrap[["t"]]
  worstCaseStoppingTimeBootstrapSe <- if (any(!is.finite(bootstrapQuantiles))) {
    Inf
  } else {
    worstCaseBootstrap[["bootSe"]]
  }

  list(
    difference = difference,
    alpha = alpha,
    beta = beta,
    maxBlocks = maxBlocks,
    thetaAValues = thetaAValues,
    thetaBValues = thetaAValues + difference,
    stoppingTimes = stoppingTimes,
    eValuesAtStopping = lapply(
      simulationResults,
      `[[`,
      "eValuesAtStopping"
    ),
    crossed = lapply(simulationResults, `[[`, "crossed"),
    blocksSimulated = lapply(simulationResults, `[[`, "blocksSimulated"),
    stoppingTimeQuantiles = stoppingTimeQuantiles,
    empiricalPower = empiricalPower,
    worstCasePower = min(empiricalPower),
    worstCaseThetaA = thetaAValues[worstCaseIndex],
    worstCaseThetaB = thetaAValues[worstCaseIndex] + difference,
    worstCaseStoppingTime = stoppingTimeQuantiles[worstCaseIndex],
    worstCaseStoppingTimeBootstrapSe = worstCaseStoppingTimeBootstrapSe,
    worstCaseStoppingTimeTwoSe = 2 * worstCaseStoppingTimeBootstrapSe
  )
}

# Design fnts ----

#' Designs a Savi Experiment to Test Two Proportions in Stream Data
#'
#' The Turner design uses independent Beta predictive distributions for the two
#' Bernoulli streams. For the currently supported design scenario, `delta`
#' describes the difference between the data-generating probabilities,
#' `thetaB - thetaA`; it is not supplied to the unrestricted Turner e-process.
#'
#' Supply `delta` and `beta` to estimate the required number of data blocks, or
#' supply only `nBlocksPlan` to construct a pilot design. The inverse design
#' scenarios and restricted alternatives are not yet implemented.
#'
#' @param na Number of observations in group A per data block.
#' @param nb Number of observations in group B per data block.
#' @param nBlocksPlan Planned number of data blocks collected.
#' @param beta Numeric in `(0, 1)` specifying the tolerable type II error.
#' @param delta Minimal relevant difference `thetaB - thetaA` used to generate
#'   data in the design simulation.
#' @param alpha Numeric in `(0, 1)` specifying the tolerable type I error.
#' @param pilot Logical specifying whether this is a pilot design.
#' @param hyperParameterValues Named list containing positive `betaA1`,
#'   `betaA2`, `betaB1`, and `betaB2` values for the two Beta priors.
#' @param previousSaviTestResult Optional previous test result whose
#'   `posteriorHyperParameters` are used as the new prior parameters.
#' @param M Number of simulations used to estimate `nBlocksPlan`.
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
  M = 1e3
) {
  if (length(na) != 1L || length(nb) != 1L ||
      any(!is.finite(c(na, nb))) || any(c(na, nb) <= 0) ||
      any(c(na, nb) %% 1 != 0)) {
    stop("na and nb must be positive integer block sizes.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }
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

  requiredPriorNames <- c("betaA1", "betaA2", "betaB1", "betaB2")
  if (!all(requiredPriorNames %in% names(hyperParameterValues))) {
    stop(
      paste(
        "Provide hyperparameters as a named list for betaA1, betaA2,",
        "betaB1 and betaB2."
      )
    )
  }
  priorValues <- unlist(hyperParameterValues[requiredPriorNames], use.names = FALSE)
  if (!is.numeric(priorValues) || any(!is.finite(priorValues)) ||
      any(priorValues <= 0)) {
    stop("Beta prior hyperparameters must be finite and greater than zero.")
  }

  names(priorValuesForPrint) <- "Beta hyperparameters"
  nPlanTwoSe <- NULL

  planningDesign <- is.null(nBlocksPlan) && !is.null(delta) && !is.null(beta)
  pilotDesign <- !is.null(nBlocksPlan) && is.null(delta) && is.null(beta)

  if (!planningDesign && !pilotDesign) {
    stop(
      paste(
        "Provide delta and beta to estimate nBlocksPlan, or provide only",
        "nBlocksPlan for a pilot design."
      )
    )
  }

  if (planningDesign) {
    if (length(delta) != 1L || !is.finite(delta) || delta == 0 || abs(delta) >= 1) {
      stop("delta must be a nonzero difference strictly between -1 and 1.")
    }
    if (length(beta) != 1L || !is.finite(beta) || beta <= 0 || beta >= 1) {
      stop("beta must be strictly between 0 and 1.")
    }

    simulationResult <- simulateWorstCaseQuantileTwoProportions(
      delta = delta,
      na = na,
      nb = nb,
      priorValues = hyperParameterValues,
      alpha = alpha,
      beta = beta,
      M = M
    )
    nBlocksPlan <- simulationResult[["nBlocksPlan"]]
    nPlanTwoSe <- c(0, 0, simulationResult[["nBlocksPlanTwoSe"]])
  } else {
    pilot <- TRUE
  }

  nPlan <- c(na, nb, nBlocksPlan)
  names(nPlan) <- c("na", "nb", "nBlocksPlan")

  if (!is.null(delta)) {
    names(delta) <- "difference"
  }

  result <- list(
    "nPlan" = nPlan,
    "nPlanTwoSe" = nPlanTwoSe,
    "parameter" = priorValuesForPrint,
    "betaPriorParameterValues" = hyperParameterValues,
    "alpha" = alpha,
    "beta" = beta,
    "betaTwoSe" = NULL,
    "logImpliedTarget" = NULL,
    "logImpliedTargetTwoSe" = NULL,
    "esMin" = delta,
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
