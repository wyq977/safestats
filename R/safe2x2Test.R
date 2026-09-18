# Likelihood-ratio calculations ----

# One term of a log likelihood ratio: count * (logP - logQ). A zero count
# contributes nothing and equal log probabilities cancel exactly, and both
# rules hold at a probability boundary, where the plain arithmetic would give
# 0 * Inf or (-Inf) - (-Inf), i.e. NaN.
countWeightedLogRatio <- function(count, logP, logQ) {
  ifelse(count == 0 | logP == logQ, 0, count * (logP - logQ))
}

# Calculate blockwise log-likelihood-ratio increments. Arguments may be
# vectors or conformable matrices; callers own predictability and accumulation.
logLikelihoodRatioIncrements <- function(
  ya,
  yb,
  na,
  nb,
  numeratorThetaA,
  numeratorThetaB,
  denominatorThetaA,
  denominatorThetaB
) {
  successesA <- countWeightedLogRatio(
    ya, log(numeratorThetaA), log(denominatorThetaA)
  )
  successesB <- countWeightedLogRatio(
    yb, log(numeratorThetaB), log(denominatorThetaB)
  )
  failuresA <- countWeightedLogRatio(
    na - ya, log1p(-numeratorThetaA), log1p(-denominatorThetaA)
  )
  failuresB <- countWeightedLogRatio(
    nb - yb, log1p(-numeratorThetaB), log1p(-denominatorThetaB)
  )

  successesA + failuresA + successesB + failuresB
}

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
  cumsum(logLikelihoodRatioIncrements(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    numeratorThetaA = numeratorThetaA,
    numeratorThetaB = numeratorThetaB,
    denominatorThetaA = denominatorThetaA,
    denominatorThetaB = denominatorThetaB
  ))
}

# Predictable theta calculations ----

# Supply and validate the four Beta prior parameters. The B-group defaults are
# scaled by nb / na so that both streams carry the same prior weight per
# observation.
resolveBetaPriorParameters <- function(
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
    stop("Beta prior parameters must be finite and greater than zero.")
  }

  priorParameters[requiredNames]
}

# The one-dimensional set of theta pairs with a fixed effect, discretised on
# `nWeight` points. The free coordinate rho is thetaA rescaled to its feasible
# interval under a propDiff restriction and thetaA itself under a logOR one;
# either way rho stays in (0, 1), so the same Beta prior applies to both.
restrictedThetaWeightGrid <- function(restriction, delta, nWeight = 1000L) {
  restriction <- match.arg(restriction, c("propDiff", "logOR"))
  if (length(delta) != 1L || !is.numeric(delta) || !is.finite(delta)) {
    stop("delta must be a finite numeric scalar for a restricted e-process.")
  }
  if (restriction == "propDiff" && abs(delta) >= 1) {
    stop("delta must lie strictly between -1 and 1 for a propDiff restriction.")
  }
  if (length(nWeight) != 1L || !is.finite(nWeight) ||
      nWeight < 2L || nWeight %% 1 != 0) {
    stop("nWeight must be an integer greater than or equal to 2.")
  }

  rhoGrid <- seq(1 / nWeight, 1 - 1 / nWeight, length.out = nWeight)
  if (restriction == "propDiff") {
    thetaA <- rhoGrid * (1 - abs(delta)) + max(-delta, 0)
  } else {
    thetaA <- rhoGrid
  }
  thetaB <- thetaBFromRestriction(thetaA, restriction, delta)

  # A support point at exactly 0 or 1 has an infinite log probability, which
  # turns a zero count into 0 * -Inf = NaN and then spreads through the
  # posterior normalisation. plogis() saturates to exactly 1 once its argument
  # reaches about 37, so a large enough logOR degenerates the whole curve.
  # There is no meaningful posterior there, so refuse to build it.
  if (any(thetaA <= 0 | thetaA >= 1 | thetaB <= 0 | thetaB >= 1)) {
    stop(paste0(
      "delta = ", format(delta), " is too extreme for a ", restriction,
      " restriction at nWeight = ", nWeight, ": it places a support point ",
      "at exactly 0 or 1, where the likelihood is degenerate. Use a delta ",
      "closer to zero."
    ))
  }

  list(rho = rhoGrid, thetaA = thetaA, thetaB = thetaB)
}

# thetaB implied by thetaA once the effect is fixed. `restriction` must already
# be normalised to "propDiff" or "logOR".
thetaBFromRestriction <- function(thetaA, restriction, delta) {
  if (restriction == "propDiff") {
    thetaA + delta
  } else {
    stats::plogis(stats::qlogis(thetaA) + delta)
  }
}

# Normalised log prior weights on restrictedThetaWeightGrid()'s free coordinate.
# Only the betaA* parameters are used: fixing the effect leaves one free
# probability, so the second Beta prior has nothing left to describe.
restrictedPriorLogWeights <- function(priorParameters, rhoGrid) {
  logWeights <-
    (priorParameters[["betaA1"]] - 1) * log(rhoGrid) +
    (priorParameters[["betaA2"]] - 1) * log1p(-rhoGrid)
  logWeights - logSumExp(logWeights)
}

# Predictable numerator probabilities for the Turner e-process. Unrestricted
# probabilities are independent Beta posterior means. Restricted probabilities
# use a grid posterior on the fixed-effect curve. In either case block i is
# predicted from blocks 1 through i - 1 only.
learnPredictiveThetas <- function(
  ya,
  yb,
  na,
  nb,
  priorParameters,
  restriction = c("none", "propDiff", "logOR"),
  delta = NULL,
  nWeight = 1000L
) {
  restriction <- match.arg(restriction)
  nSteps <- length(ya)

  if (restriction == "none") {
    previousYa <- c(0, cumsum(ya))[seq_len(nSteps)]
    previousYb <- c(0, cumsum(yb))[seq_len(nSteps)]
    previousNa <- c(0, cumsum(na))[seq_len(nSteps)]
    previousNb <- c(0, cumsum(nb))[seq_len(nSteps)]

    return(list(
      thetaA = (priorParameters[["betaA1"]] + previousYa) / (
        previousNa +
          priorParameters[["betaA1"]] + priorParameters[["betaA2"]]
      ),
      thetaB = (priorParameters[["betaB1"]] + previousYb) / (
        previousNb +
          priorParameters[["betaB1"]] + priorParameters[["betaB2"]]
      )
    ))
  }

  weightGrid <- restrictedThetaWeightGrid(restriction, delta, nWeight)
  logWeights <- restrictedPriorLogWeights(priorParameters, weightGrid[["rho"]])
  logThetaA <- log(weightGrid[["thetaA"]])
  logOneMinusThetaA <- log1p(-weightGrid[["thetaA"]])
  logThetaB <- log(weightGrid[["thetaB"]])
  logOneMinusThetaB <- log1p(-weightGrid[["thetaB"]])
  thetaA <- thetaB <- numeric(nSteps)

  for (block in seq_along(ya)) {
    thetaA[block] <- sum(weightGrid[["thetaA"]] * exp(logWeights))
    thetaB[block] <- thetaBFromRestriction(thetaA[block], restriction, delta)

    logWeights <- logWeights +
      ya[block] * logThetaA +
      (na[block] - ya[block]) * logOneMinusThetaA +
      yb[block] * logThetaB +
      (nb[block] - yb[block]) * logOneMinusThetaB
    logWeights <- logWeights - logSumExp(logWeights)
  }

  list(thetaA = thetaA, thetaB = thetaB)
}

# Simulation theta grids ----

# Construct feasible data-generating theta pairs for one fixed effect.
makeSimulationThetaGrid <- function(
  propDiff = NULL,
  logOR = NULL,
  nTheta = 8L
) {
  if (is.null(propDiff) == is.null(logOR)) {
    stop("Supply exactly one of propDiff and logOR.")
  }
  if (!is.null(propDiff) &&
      (length(propDiff) != 1L || !is.finite(propDiff) ||
       propDiff == 0 || abs(propDiff) >= 1)) {
    stop("propDiff must be nonzero and strictly between -1 and 1.")
  }
  if (!is.null(logOR) &&
      (length(logOR) != 1L || !is.finite(logOR) || logOR == 0)) {
    stop("logOR must be a nonzero finite scalar.")
  }
  if (length(nTheta) != 1L || !is.finite(nTheta) ||
      nTheta < 1L || nTheta %% 1 != 0) {
    stop("nTheta must be a positive integer.")
  }

  rho <- seq_len(nTheta) / (nTheta + 1)
  restriction <- if (is.null(propDiff)) "logOR" else "propDiff"
  delta <- if (is.null(propDiff)) logOR else propDiff
  thetaA <- if (restriction == "propDiff") {
    rho * (1 - abs(delta)) + max(-delta, 0)
  } else {
    rho
  }

  # The restriction and its value travel with the grid: the simulated
  # e-process has to place its numerator on this same effect curve, not on an
  # unrestricted alternative.
  list(
    thetaA = thetaA,
    thetaB = thetaBFromRestriction(thetaA, restriction, delta),
    restriction = restriction,
    delta = delta
  )
}

# Equality-null e-process: built only from the predictive-theta and prior
# helpers above. Unlike the propDiff and logOR sections
# below, it tests thetaB - thetaA = 0 directly and never solves a RIPr.
#' Cumulative Turner e-process for two Bernoulli streams
#'
#' Computes the Turner e-process for the equality null `thetaB - thetaA = 0`.
#' By default the numerator uses independent Beta predictors for the two
#' proportions. A restricted numerator can instead impose either
#' `propDiff = thetaB - thetaA = delta` or
#' `logOR = logit(thetaB) - logit(thetaA) = delta`; its
#' grid posterior is updated with each observed block before predicting the
#' next. The returned vector is cumulative: element `i` uses blocks `1`
#' through `i`.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param na,nb Number of observations in groups A and B in each data block.
#'   A scalar is recycled over all blocks.
#' @param priorParameters Optional named list with `betaA1`, `betaA2`,
#'   `betaB1`, and `betaB2`. Restricted predictors use the `betaA*` prior on
#'   their one-dimensional grid.
#' @param restriction Restriction placed on the numerator alternative: `"none"`
#'   (default), `"propDiff"`, or `"logOR"`.
#' @param delta The restricted numerator's value of `thetaB - thetaA` or
#'   `logit(thetaB) - logit(thetaA)`, according to `restriction`.
#' @param nWeight Number of points in the restricted-alternative posterior
#'   grid.
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
  log = FALSE,
  restriction = c("none", "propDiff", "logOR"),
  delta = NULL,
  nWeight = 1000L
) {
  restriction <- match.arg(restriction)
  nSteps <- length(ya)

  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)
  if (!all(c(length(yb), length(na), length(nb)) == nSteps)) {
    stop("ya, yb, na, and nb must have the same length.")
  }
  values <- c(ya, yb, na, nb)
  if (!is.numeric(values) || any(!is.finite(values)) ||
      any(values %% 1 != 0)) {
    stop("ya, yb, na, and nb must contain finite integer values.")
  }
  if (any(na <= 0) || any(nb <= 0)) {
    stop("na and nb must contain positive values.")
  }
  if (any(ya < 0) || any(yb < 0) || any(ya > na) || any(yb > nb)) {
    stop("Success counts must lie between zero and their group sizes.")
  }

  priorParameters <- resolveBetaPriorParameters(
    na = na,
    nb = nb,
    priorParameters = priorParameters
  )

  predictiveThetas <- learnPredictiveThetas(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    priorParameters = priorParameters,
    restriction = restriction,
    delta = delta,
    nWeight = nWeight
  )
  pooledPredictiveTheta <- (
    na * predictiveThetas[["thetaA"]] + nb * predictiveThetas[["thetaB"]]
  ) / (na + nb)
  logEProcess <- logLikelihoodRatioProcess(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    numeratorThetaA = predictiveThetas[["thetaA"]],
    numeratorThetaB = predictiveThetas[["thetaB"]],
    denominatorThetaA = pooledPredictiveTheta,
    denominatorThetaB = pooledPredictiveTheta
  )

  if (log) logEProcess else exp(logEProcess)
}

# Proportion difference: RIPr projection and confidence sequence ----

# Effect convention:
#   propDiff = thetaB - thetaA.
# A positive propDiff means thetaB > thetaA.

# Solve the reverse information projection onto thetaB - thetaA = propDiff
# for one data block. thetaA and thetaB are the predictor's Bernoulli
# probabilities being projected and na and nb the block's group sizes. The KL
# projection's first-order condition reduces to a cubic in the null thetaA, of
# which exactly one root lies in the feasible interval; that root is returned.
solvePropDiffRIPr <- function(
  thetaA,
  thetaB,
  na,
  nb,
  propDiff
) {
  A <- na
  B <- nb
  coefficients <- c(
    -A * thetaA * propDiff * (1 - propDiff),
    A * (
      propDiff * (1 - propDiff) -
        thetaA * (1 - 2 * propDiff)
    ) + B * (propDiff - thetaB),
    A * (1 - 2 * propDiff + thetaA) +
      B * (1 - propDiff + thetaB),
    -(A + B)
  )
  roots <- base::polyroot(coefficients)
  tolerance <- sqrt(.Machine$double.eps)
  realRoots <- Re(roots)[abs(Im(roots)) < tolerance]
  feasibleLower <- max(0, -propDiff)
  feasibleUpper <- min(1, 1 - propDiff)
  feasibleRoots <- realRoots[
    realRoots > feasibleLower + tolerance &
      realRoots < feasibleUpper - tolerance
  ]
  if (length(feasibleRoots) != 1L) {
    stop("Could not identify a unique feasible propDiff RIPr.")
  }
  feasibleRoots
}

#' Confidence sequence for the proportion difference
#'
#' Walks through the data one block at a time. At every block the predictor
#' learned from the earlier blocks is projected onto each candidate
#' `propDiff = thetaB - thetaA` still in the confidence set, and that
#' candidate's log e-process is advanced by one likelihood-ratio increment.
#' Each candidate is a two-sided point null tested at level `alpha`: it
#' leaves the set once its e-process reaches `1 / alpha`. The reported bounds
#' are the smallest and largest candidates remaining after each block.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param confidenceBoundGridPrecision Number of candidate proportion
#'   differences on each side of zero. Zero itself is never a candidate, so a
#'   bound of exactly zero is never reported; the grid holds
#'   `2 * confidenceBoundGridPrecision` values in the open interval `(-1, 1)`.
#' @param saviDesign A `saviDesign` returned by
#'   [designSaviTwoProportions()].
#' @param runningIntersection Logical. If `TRUE` (default), a candidate
#'   rejected at some block stays rejected at every later block, so the sets
#'   are nested and a rejected candidate's e-process is no longer updated. If
#'   `FALSE`, every candidate is updated at every block and each block's set
#'   is read from the cumulative e-process values at that block only, with no
#'   memory of earlier rejections, so a rejected candidate may re-enter
#'   later. Both versions have the same time-uniform coverage; the
#'   intersection is never wider.
#'
#' @return A data frame with `block`, `lowerBound`, and `upperBound`. Bounds
#'   remain at `-1` or `1` while their corresponding edge candidate remains;
#'   if no grid point remains at a block, both bounds are `NA`.
computeConfidenceSequenceForPropDiffTwoProportions <- function(
  ya,
  yb,
  confidenceBoundGridPrecision,
  saviDesign,
  runningIntersection = TRUE
) {
  if (length(confidenceBoundGridPrecision) != 1L ||
      !is.finite(confidenceBoundGridPrecision) ||
      confidenceBoundGridPrecision < 1L ||
      confidenceBoundGridPrecision %% 1 != 0) {
    stop("confidenceBoundGridPrecision must be a positive integer.")
  }
  if (!is.logical(runningIntersection) || length(runningIntersection) != 1L ||
      is.na(runningIntersection)) {
    stop("runningIntersection must be TRUE or FALSE.")
  }

  nBlocks <- length(ya)
  na <- rep_len(saviDesign[["nPlan"]][["na"]], nBlocks)
  nb <- rep_len(saviDesign[["nPlan"]][["nb"]], nBlocks)
  if (length(yb) != nBlocks) {
    stop("ya and yb must have the same length.")
  }
  logThreshold <- log(1 / saviDesign[["alpha"]])

  # The predictor for block t depends on the data before t only, so it can be
  # learned for every block up front.
  predictiveThetas <- learnPredictiveThetas(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    priorParameters = saviDesign[["betaPriorParameterValues"]],
    restriction = "none"
  )

  # Equally spaced candidates in (0, 1), mirrored about zero. Zero is
  # deliberately not a candidate, matching the log-odds-ratio grid: a bound of
  # exactly 0 has no unambiguous reading, so a finite bound always sits
  # strictly on one side of the null.
  positiveGrid <- seq_len(confidenceBoundGridPrecision) /
    (confidenceBoundGridPrecision + 1)
  propDiffGrid <- c(-rev(positiveGrid), positiveGrid)
  nCandidates <- length(propDiffGrid)
  logEProcess <- numeric(nCandidates)
  inSet <- rep(TRUE, nCandidates)
  lowerBound <- upperBound <- rep(NA_real_, nBlocks)

  for (block in seq_len(nBlocks)) {
    thetaStarA <- predictiveThetas[["thetaA"]][block]
    thetaStarB <- predictiveThetas[["thetaB"]][block]

    # With the running intersection a rejected candidate never returns, so
    # its e-process need not be advanced any further.
    toUpdate <- if (runningIntersection) which(inSet) else seq_len(nCandidates)
    for (candidate in toUpdate) {
      propDiff <- propDiffGrid[candidate]
      nullThetaA <- solvePropDiffRIPr(
        thetaA = thetaStarA,
        thetaB = thetaStarB,
        na = na[block],
        nb = nb[block],
        propDiff = propDiff
      )
      logEProcess[candidate] <- logEProcess[candidate] +
        logLikelihoodRatioIncrements(
          ya = ya[block],
          yb = yb[block],
          na = na[block],
          nb = nb[block],
          numeratorThetaA = thetaStarA,
          numeratorThetaB = thetaStarB,
          denominatorThetaA = nullThetaA,
          denominatorThetaB = nullThetaA + propDiff
        )
    }

    notRejected <- logEProcess < logThreshold
    inSet <- if (runningIntersection) inSet & notRejected else notRejected

    # If the outermost candidate is still in, the true value may lie beyond
    # the grid, so the bound stays at the parameter limit. If no candidate is
    # left, both bounds stay NA.
    if (any(inSet)) {
      lowerBound[block] <-
        if (inSet[1L]) -1 else min(propDiffGrid[inSet])
      upperBound[block] <-
        if (inSet[nCandidates]) 1 else max(propDiffGrid[inSet])
    }
  }

  data.frame(
    block = seq_len(nBlocks),
    lowerBound = lowerBound,
    upperBound = upperBound
  )
}

# Log odds ratio: RIPr projection and confidence sequence ----

# Effect convention:
#   logOR = logit(thetaB) - logit(thetaA)
#         = log(oddsB / oddsA).
# A positive logOR means oddsB > oddsA.

# Log of exp(x) + exp(y), computed pairwise without exponentiating either
# term on its own. Vectorised over both arguments.
logAddExp <- function(x, y) {
  pmax(x, y) + log1p(exp(-abs(x - y)))
}

# Positive root of a*x^2 + b*x + c = 0, returned as log(root). Assumes a > 0
# and c < 0, which guarantees exactly one positive root, and takes log(a) and
# log(|c|) rather than a and c so that callers whose coefficients underflow to
# zero in linear space can still pass them as ordinary finite logs.
#
# The two quadratic formulas below are algebraically identical, but each
# cancels catastrophically where the other does not: (-b + sqrt(disc)) / (2a)
# subtracts near-equal terms when b > 0, and 2c / (-b - sqrt(disc)) does the
# same when b < 0. Picking by the sign of b always takes the sum, never the
# difference.
logPositiveQuadraticRoot <- function(logA, b, logAbsC) {
  logAbsB <- log(abs(b))
  # sqrt(b^2 - 4ac), with -4ac = 4a|c| positive, so the discriminant is a sum.
  logSqrtDiscriminant <-
    0.5 * logAddExp(2 * logAbsB, log(4) + logA + logAbsC)
  logAbsBPlusSqrtDiscriminant <- logAddExp(logAbsB, logSqrtDiscriminant)

  ifelse(
    b >= 0,
    log(2) + logAbsC - logAbsBPlusSqrtDiscriminant,  # 2c / (-b - sqrt(disc))
    logAbsBPlusSqrtDiscriminant - log(2) - logA      # (-b + sqrt(disc)) / (2a)
  )
}

# Solve the reverse information projection onto
# logit(thetaB) - logit(thetaA) = logOR. As in solvePropDiffRIPr, thetaA and
# thetaB are the predictor's Bernoulli probabilities being projected and na and
# nb the group sizes.
#
# As with solvePropDiffRIPr, the KL projection's first-order condition reduces
# to matching the predictor's weighted mean of successes,
#
#   na*nullThetaA + nb*nullThetaB = na*thetaA + nb*thetaB =: successes
#
# subject to the constraint. Writing each theta as odds/(1 + odds) and clearing
# the two denominators turns that into a quadratic in one group's odds; the
# constraint supplies the other group's odds as a fixed multiple of it.
#
# Which group to solve for is a numerical choice. Parameterising by the group
# with the larger constrained logit (the "base" group) makes the multiplier for
# the other one r = exp(-abs(logOR)), which stays in (0, 1] and so keeps the
# quadratic's coefficients bounded for either sign of a large logOR. With
# failures := na + nb - successes, the quadratic in the base odds u is
#
#   r*failures*u^2 + (baseSize - successes + r*(shiftedSize - successes))*u
#     - successes = 0
#
# whose leading coefficient is positive and constant term negative whenever
# 0 < successes < na + nb, leaving exactly one positive root.
#
# Every argument is elementwise and recycled to a common length, so one call
# solves a single candidate over many blocks (theta vectors with a scalar
# logOR) or many candidates at a single block (scalar thetas with a logOR
# vector), which is how the confidence sequence uses it.
solveLogORRIPr <- function(
  thetaA,
  thetaB,
  na,
  nb,
  logOR
) {
  argumentLengths <- c(
    length(thetaA),
    length(thetaB),
    length(na),
    length(nb),
    length(logOR)
  )
  nSolves <- max(argumentLengths)
  if (!all(argumentLengths == 1L | argumentLengths == nSolves)) {
    stop(
      "Thetas, group sizes, and logOR must have length one or a common ",
      "length."
    )
  }
  if (any(!is.finite(logOR))) {
    stop("logOR must be finite.")
  }
  thetaA <- rep_len(thetaA, nSolves)
  thetaB <- rep_len(thetaB, nSolves)
  na <- rep_len(na, nSolves)
  nb <- rep_len(nb, nSolves)
  logOR <- rep_len(logOR, nSolves)

  totalSize <- na + nb
  successes <- na * thetaA + nb * thetaB
  nullThetaA <- nullThetaB <- numeric(nSolves)

  # At these boundaries the mean-matching condition only holds at theta = 0
  # or 1; theta = 0 is already the default value set above.
  atLowerBoundary <- successes == 0
  atUpperBoundary <- successes == totalSize
  nullThetaA[atUpperBoundary] <- nullThetaB[atUpperBoundary] <- 1
  interior <- !atLowerBoundary & !atUpperBoundary

  if (!any(interior)) {
    return(list(thetaA = nullThetaA, thetaB = nullThetaB))
  }

  interiorLogOR <- logOR[interior]
  interiorNa <- na[interior]
  interiorNb <- nb[interior]
  interiorSuccesses <- successes[interior]
  interiorFailures <- totalSize[interior] - interiorSuccesses

  # logOR = logit(thetaB) - logit(thetaA), so A is the base group exactly when
  # logOR is negative.
  aIsBase <- interiorLogOR < 0
  baseSize <- ifelse(aIsBase, interiorNa, interiorNb)
  shiftedSize <- ifelse(aIsBase, interiorNb, interiorNa)
  r <- exp(-abs(interiorLogOR))

  # coefA and coefC are passed as logs: r underflows to exactly 0 for a large
  # abs(logOR) even though log(r) = -abs(logOR) stays an ordinary finite
  # number, so log(r * failures) is formed as a sum rather than via log() of
  # the linear-space product.
  logCoefA <- log(interiorFailures) - abs(interiorLogOR)
  coefB <- baseSize - interiorSuccesses + r * (shiftedSize - interiorSuccesses)
  logAbsCoefC <- log(interiorSuccesses)  # coefC = -interiorSuccesses

  baseLogit <- logPositiveQuadraticRoot(logCoefA, coefB, logAbsCoefC)

  interiorThetaA <- ifelse(
    aIsBase,
    stats::plogis(baseLogit),
    stats::plogis(baseLogit - interiorLogOR)
  )
  interiorThetaB <- ifelse(
    aIsBase,
    stats::plogis(baseLogit + interiorLogOR),
    stats::plogis(baseLogit)
  )

  # The equality null: base and shifted groups coincide at the pooled
  # proportion. The quadratic has that root too; this only keeps it exact.
  pooled <- interiorLogOR == 0
  interiorThetaA[pooled] <- interiorThetaB[pooled] <-
    interiorSuccesses[pooled] /
      (interiorSuccesses[pooled] + interiorFailures[pooled])

  nullThetaA[interior] <- interiorThetaA
  nullThetaB[interior] <- interiorThetaB
  list(thetaA = nullThetaA, thetaB = nullThetaB)
}

#' Confidence sequence for the log odds ratio
#'
#' Constructs a symmetric candidate grid from the supplied resolution and
#' search bounds; the grid excludes zero. The log odds ratio is
#' `logit(thetaB) - logit(thetaA)`. Two one-sided families are inverted, each
#' at `alpha / 2`: the lower family tests `logOR <= candidate`, the upper
#' family `logOR >= candidate`, and both run over the full signed grid.
#'
#' The function walks through the data one block at a time. At every block
#' the predictor learned from the earlier blocks lies inside one family's
#' null for each candidate, where that family's likelihood-ratio increment is
#' zero, and outside the other family's null, where the predictor is
#' projected onto the candidate and the increment is added to that family's
#' running log e-process. Each family keeps a running intersection: a
#' candidate leaves a family's set once its log e-process reaches
#' `log(2 / alpha)` and is never projected again. The lower bound after a
#' block is the smallest candidate still in the lower family's set and the
#' upper bound the largest candidate still in the upper family's set.
#'
#' @param ya,yb Number of successes in groups A and B in each data block.
#' @param confidenceBoundGridPrecision Number of candidate values on each side
#'   of zero. Zero itself is not a candidate, so a bound of exactly zero is
#'   never reported; the grid holds `2 * confidenceBoundGridPrecision` values.
#' @param logORConfidenceSearchBounds Positive finite lower and upper bounds
#'   for the absolute candidate log odds ratios.
#' @param saviDesign A `saviDesign` returned by
#'   [designSaviTwoProportions()].
#'
#' @return A data frame with `block`, `lowerBound`, and `upperBound`. A bound is
#'   infinite until its one-sided sequence rejects a candidate and is `NA` if
#'   every candidate on that side has been rejected.
computeConfidenceSequenceForLogORTwoProportions <- function(
  ya,
  yb,
  confidenceBoundGridPrecision,
  logORConfidenceSearchBounds,
  saviDesign
) {
  if (length(confidenceBoundGridPrecision) != 1L ||
      !is.finite(confidenceBoundGridPrecision) ||
      confidenceBoundGridPrecision < 1L ||
      confidenceBoundGridPrecision %% 1 != 0) {
    stop("confidenceBoundGridPrecision must be a positive integer.")
  }
  if (length(logORConfidenceSearchBounds) != 2L ||
      any(!is.finite(logORConfidenceSearchBounds)) ||
      logORConfidenceSearchBounds[1L] <= 0 ||
      logORConfidenceSearchBounds[1L] >=
        logORConfidenceSearchBounds[2L]) {
    stop(
      paste(
        "logORConfidenceSearchBounds must contain two positive,",
        "increasing finite values."
      )
    )
  }

  nBlocks <- length(ya)
  if (length(yb) != nBlocks) {
    stop("ya and yb must have the same length.")
  }
  na <- rep_len(saviDesign[["nPlan"]][["na"]], nBlocks)
  nb <- rep_len(saviDesign[["nPlan"]][["nb"]], nBlocks)
  logThreshold <- log(2 / saviDesign[["alpha"]])

  # The predictor for block t depends on the data before t only, so it can be
  # learned for every block up front.
  predictiveThetas <- learnPredictiveThetas(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    priorParameters = saviDesign[["betaPriorParameterValues"]],
    restriction = "none"
  )
  predictorLogOR <- stats::qlogis(predictiveThetas[["thetaB"]]) -
    stats::qlogis(predictiveThetas[["thetaA"]])

  # use tanh to transform log uniform in log odds ratio
  positiveTransformedBounds <- tanh(logORConfidenceSearchBounds / 4)

  positiveGrid <- 4 * atanh(seq(
    positiveTransformedBounds[1L],
    positiveTransformedBounds[2L],
    length.out = confidenceBoundGridPrecision
  ))

  # Zero is deliberately not a candidate, matching the proportion-difference
  # grid: a bound of exactly 0 has no unambiguous reading, so a finite bound
  # always sits strictly on one side of the null.
  candidateGrid <- c(-rev(positiveGrid), positiveGrid)
  nCandidates <- length(candidateGrid)

  lowerLogEProcess <- upperLogEProcess <- numeric(nCandidates)
  lowerInSet <- upperInSet <- rep(TRUE, nCandidates)
  lowerBound <- upperBound <- rep(NA_real_, nBlocks)

  for (block in seq_len(nBlocks)) {
    thetaStarA <- predictiveThetas[["thetaA"]][block]
    thetaStarB <- predictiveThetas[["thetaB"]][block]

    # For a candidate below the predictor, the predictor lies outside the
    # lower null logOR <= candidate and inside the upper null, so only the
    # lower family gains a nonzero increment; above the predictor the roles
    # swap. The projection does not depend on the family, so each candidate
    # is solved at most once per block, and only while the family that needs
    # it still holds the candidate.
    lowerNeedsSolve <- lowerInSet & candidateGrid < predictorLogOR[block]
    upperNeedsSolve <- upperInSet & candidateGrid > predictorLogOR[block]
    toSolve <- which(lowerNeedsSolve | upperNeedsSolve)
    if (length(toSolve) > 0L) {
      nullTheta <- solveLogORRIPr(
        thetaA = thetaStarA,
        thetaB = thetaStarB,
        na = na[block],
        nb = nb[block],
        logOR = candidateGrid[toSolve]
      )
      increments <- logLikelihoodRatioIncrements(
        ya = ya[block],
        yb = yb[block],
        na = na[block],
        nb = nb[block],
        numeratorThetaA = thetaStarA,
        numeratorThetaB = thetaStarB,
        denominatorThetaA = nullTheta[["thetaA"]],
        denominatorThetaB = nullTheta[["thetaB"]]
      )
      forLower <- lowerNeedsSolve[toSolve]
      lowerLogEProcess[toSolve[forLower]] <-
        lowerLogEProcess[toSolve[forLower]] + increments[forLower]
      upperLogEProcess[toSolve[!forLower]] <-
        upperLogEProcess[toSolve[!forLower]] + increments[!forLower]
    }

    lowerInSet <- lowerInSet & lowerLogEProcess < logThreshold
    upperInSet <- upperInSet & upperLogEProcess < logThreshold

    # If the outermost candidate on a side is still in, the true value may lie
    # beyond the grid, so that bound stays infinite. If a family has rejected
    # every candidate, its bound stays NA.
    if (any(lowerInSet)) {
      lowerBound[block] <-
        if (lowerInSet[1L]) -Inf else min(candidateGrid[lowerInSet])
    }
    if (any(upperInSet)) {
      upperBound[block] <-
        if (upperInSet[nCandidates]) Inf else max(candidateGrid[upperInSet])
    }
  }

  data.frame(
    block = seq_len(nBlocks),
    lowerBound = lowerBound,
    upperBound = upperBound
  )
}

# Stopping-time simulation ----

# Simulate one path of the restricted Turner e-process until it crosses
# 1 / alpha or reaches maxBlocks. All maxBlocks blocks are drawn up front,
# which is cheap; the process is then evaluated block by block and the loop
# returns at the first crossing, so no block past the stopping time is
# evaluated. The numerator is the grid posterior on the effect curve, predicted
# from the blocks before the current one; the denominator is the pooled
# predictable probability, as in turnerEProcess(). Returns the stopping time
# (Inf if the path never crosses) and the e-value there (NA if it never
# crosses). Arguments are validated by the caller.
simulateStoppingTimeWithRestriction <- function(
  thetaA,
  thetaB,
  na,
  nb,
  restriction,
  delta,
  alpha,
  priorParameters,
  nWeight,
  maxBlocks
) {
  weightGrid <- restrictedThetaWeightGrid(restriction, delta, nWeight)
  weightGridThetaA <- weightGrid[["thetaA"]]
  logThetaA <- log(weightGridThetaA)
  logOneMinusThetaA <- log1p(-weightGridThetaA)
  logThetaB <- log(weightGrid[["thetaB"]])
  logOneMinusThetaB <- log1p(-weightGrid[["thetaB"]])
  logWeights <- restrictedPriorLogWeights(priorParameters, weightGrid[["rho"]])
  logThreshold <- log(1 / alpha)

  ya <- stats::rbinom(maxBlocks, na, thetaA)
  yb <- stats::rbinom(maxBlocks, nb, thetaB)
  logEValue <- 0

  for (block in seq_len(maxBlocks)) {
    # Predict this block from the posterior over the blocks before it.
    weights <- exp(logWeights)
    numeratorThetaA <- sum(weightGridThetaA * weights) / sum(weights)
    numeratorThetaB <- thetaBFromRestriction(
      numeratorThetaA,
      restriction,
      delta
    )
    pooledTheta <- (na * numeratorThetaA + nb * numeratorThetaB) / (na + nb)

    logEValue <- logEValue + logLikelihoodRatioIncrements(
      ya = ya[block],
      yb = yb[block],
      na = na,
      nb = nb,
      numeratorThetaA = numeratorThetaA,
      numeratorThetaB = numeratorThetaB,
      denominatorThetaA = pooledTheta,
      denominatorThetaB = pooledTheta
    )
    if (logEValue >= logThreshold) {
      return(list(stoppingTime = block, eValue = exp(logEValue)))
    }

    # Only now add the block to the posterior, so the next block is predicted
    # from the past alone. Only ratios of weights matter; pinning the maximum
    # at zero keeps the weights representable over arbitrarily many blocks.
    logWeights <- logWeights +
      ya[block] * logThetaA +
      (na - ya[block]) * logOneMinusThetaA +
      yb[block] * logThetaB +
      (nb - yb[block]) * logOneMinusThetaB
    logWeights <- logWeights - max(logWeights)
  }

  list(stoppingTime = Inf, eValue = NA_real_)
}

# Simulate stopping times for G fixed data-generating theta pairs: nSim
# independent calls to simulateStoppingTimeWithRestriction() per pair. Rows of
# each output matrix correspond to theta-pair indices and columns to paths.
# Non-crossing paths have stopping time Inf and e-value NA.
sampleStoppingTimesSaviTwoProportions <- function(
  thetaA,
  thetaB,
  na,
  nb,
  alpha = 0.05,
  priorParameters = NULL,
  restriction = c("propDiff", "logOR"),
  delta = NULL,
  nWeight = 1000L,
  nSim = 1e3,
  maxBlocks = 1e4
) {
  if (length(na) != 1L || length(nb) != 1L ||
      any(!is.finite(c(na, nb))) || any(c(na, nb) <= 0) ||
      any(c(na, nb) %% 1 != 0)) {
    stop("na and nb must be positive integer block sizes.")
  }
  if (!is.numeric(thetaA) || !is.numeric(thetaB) ||
      length(thetaA) < 1L || length(thetaA) != length(thetaB) ||
      any(!is.finite(c(thetaA, thetaB))) ||
      any(c(thetaA, thetaB) < 0) || any(c(thetaA, thetaB) > 1)) {
    stop(paste(
      "thetaA and thetaB must be equal-length, nonempty numeric vectors",
      "of probabilities between 0 and 1."
    ))
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }
  if (any(!is.finite(c(nSim, maxBlocks))) ||
      any(c(nSim, maxBlocks) < 1) ||
      any(c(nSim, maxBlocks) %% 1 != 0)) {
    stop("nSim and maxBlocks must be positive integers.")
  }
  restriction <- match.arg(restriction)
  priorParameters <- resolveBetaPriorParameters(
    na = na,
    nb = nb,
    priorParameters = priorParameters
  )

  nTheta <- length(thetaA)
  stoppingTimes <- matrix(Inf, nrow = nTheta, ncol = nSim)
  eValuesAtStopping <- matrix(NA_real_, nrow = nTheta, ncol = nSim)

  for (thetaIndex in seq_len(nTheta)) {
    for (path in seq_len(nSim)) {
      simulated <- simulateStoppingTimeWithRestriction(
        thetaA[thetaIndex], thetaB[thetaIndex], na, nb,
        restriction, delta, alpha, priorParameters, nWeight, maxBlocks
      )
      stoppingTimes[thetaIndex, path] <- simulated[["stoppingTime"]]
      eValuesAtStopping[thetaIndex, path] <- simulated[["eValue"]]
    }
  }

  list(
    stoppingTimes = stoppingTimes,
    eValuesAtStopping = eValuesAtStopping
  )
}

#' Simulate worst-case stopping times for two proportions
#'
#' Simulates the Turner e-process under either a fixed proportion difference or
#' a fixed log odds ratio over a grid of feasible baseline probabilities. A
#' path that does not cross `1 / alpha` within `maxBlocks` has stopping time
#' `Inf` and stopping e-value `NA`.
#'
#' The simulated numerator is restricted to the supplied effect, so the
#' stopping times describe the e-process that
#' [turnerEProcess()] runs with the same `restriction` and `delta` rather than
#' the unrestricted one.
#'
#' The baseline grid is a separate simulation device: it ranges over the values
#' of `thetaA` for which both probabilities implied by the effect are in
#' `(0, 1)`. The worst-case stopping time is the largest empirical
#' `(1 - beta)` quantile over this grid. The bootstrap resamples only the
#' already-identified worst-case row: with `nSim` paths already spent
#' on that theta pair, its own Monte Carlo noise is the relevant uncertainty.
#' To refine which baseline is worst, increase `nTheta` rather than
#' re-deriving the worst case inside the bootstrap.
#'
#' @param na,nb Number of observations in groups A and B per block.
#' @param propDiff Data-generating proportion difference `thetaB - thetaA`.
#'   Supply exactly one of `propDiff` and `logOR`.
#' @param logOR Data-generating value `logit(thetaB) - logit(thetaA)`.
#'   Supply exactly one of `propDiff` and `logOR`.
#' @param alpha E-process rejection threshold is `1 / alpha`.
#' @param beta Target type-II error used to select the stopping-time quantile.
#' @param priorParameters Optional Turner Beta prior parameters.
#' @param nSim Number of simulated paths at each baseline probability.
#' @param maxBlocks Maximum number of blocks simulated per path.
#' @param nTheta Number of equally spaced feasible baseline theta pairs
#'   generated for this effect.
#' @param nBoot Number of nonparametric bootstrap samples used to estimate the
#'   standard error of the worst-case stopping-time quantile; must be at least
#'   two.
#' @param nWeight Number of points in the restricted numerator's posterior
#'   grid. It should match the `nWeight` the test itself will use; lowering it
#'   trades numerator accuracy for simulation speed.
#'
#' @return A list with the worst-case baseline probabilities, stopping-time
#'   quantiles, and `nTheta` by `nSim` matrices containing
#'   the per-path stopping results. Non-crossing paths have stopping time `Inf`
#'   and stopping e-value `NA`.
#'   `worstCaseStoppingTimeBootstrapSe` is the bootstrap standard error of the
#'   worst-case row's quantile, and `worstCaseStoppingTimeTwoSe` is twice that
#'   standard error.
#' @export
computeNPlanSaviTwoProportions <- function(
  na,
  nb,
  propDiff = NULL,
  logOR = NULL,
  alpha = 0.05,
  beta = 0.2,
  priorParameters = NULL,
  nSim = 1e3,
  maxBlocks = 1e4,
  nTheta = 8L,
  nBoot = 1e3,
  nWeight = 1000L
) {
  if (length(beta) != 1L || !is.finite(beta) || beta <= 0 || beta >= 1) {
    stop("beta must be strictly between 0 and 1.")
  }
  if (length(nBoot) != 1L || !is.finite(nBoot) ||
      nBoot < 2L || nBoot %% 1 != 0) {
    stop("nBoot must be an integer of at least two.")
  }

  thetaGrid <- makeSimulationThetaGrid(
    propDiff = propDiff,
    logOR = logOR,
    nTheta = nTheta
  )
  gridSimulation <- sampleStoppingTimesSaviTwoProportions(
    thetaA = thetaGrid[["thetaA"]],
    thetaB = thetaGrid[["thetaB"]],
    na = na,
    nb = nb,
    alpha = alpha,
    priorParameters = priorParameters,
    restriction = thetaGrid[["restriction"]],
    delta = thetaGrid[["delta"]],
    nWeight = nWeight,
    nSim = nSim,
    maxBlocks = maxBlocks
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
  # Resample only the worst-case row identified above (matching the previous
  # implementation's approach), not the argmax over a resampled grid: with
  # nSim paths already spent there, that row's own Monte Carlo noise
  # is what we want the SE to reflect. A grid resampled every replicate would
  # also fold in "which baseline is worst" selection noise; if that's a
  # concern, increase nTheta instead of picking it up here.
  bootstrapStoppingTimes <- replicate(nBoot, {
    stats::quantile(
      sample(
        stoppingTimes[worstCaseIndex, ],
        size = ncol(stoppingTimes),
        replace = TRUE
      ),
      probs = 1 - beta,
      names = FALSE
    )
  })
  bootstrapSe <- if (any(!is.finite(bootstrapStoppingTimes))) {
    Inf
  } else {
    stats::sd(bootstrapStoppingTimes)
  }

  if (!is.finite(stoppingTimeQuantiles[worstCaseIndex])) {
    crossingFraction <- mean(is.finite(stoppingTimes[worstCaseIndex, ]))
    warning(sprintf(paste(
      "Cannot determine the worst-case stopping time with maxBlocks = %g:",
      "only %.1f%% of paths at the worst-case baseline (thetaA = %g,",
      "thetaB = %g) crossed 1 / alpha within that horizon, which is short of",
      "the %.1f%% needed for the 1 - beta quantile. Try increasing maxBlocks."
    ), maxBlocks, 100 * crossingFraction,
      thetaGrid[["thetaA"]][worstCaseIndex],
      thetaGrid[["thetaB"]][worstCaseIndex], 100 * (1 - beta)))
  }

  c(list(
    propDiff = propDiff,
    logOR = logOR,
    alpha = alpha,
    maxBlocks = maxBlocks,
    thetaA = thetaGrid[["thetaA"]],
    thetaB = thetaGrid[["thetaB"]]
  ), gridSimulation, list(
    beta = beta,
    stoppingTimeQuantiles = stoppingTimeQuantiles,
    worstCaseIndex = worstCaseIndex,
    worstCaseThetaA = thetaGrid[["thetaA"]][worstCaseIndex],
    worstCaseThetaB = thetaGrid[["thetaB"]][worstCaseIndex],
    worstCaseStoppingTime = stoppingTimeQuantiles[worstCaseIndex],
    worstCaseStoppingTimeBootstrapSe = bootstrapSe,
    worstCaseStoppingTimeTwoSe = 2 * bootstrapSe
  ))
}

# Simulate the lowest rejection probability over the feasible baseline grid at
# a fixed stopping horizon. Rows are baseline-grid points and columns are
# independent paths. The bootstrap resamples only the already-identified
# worst-case row; see the comment in computeNPlanSaviTwoProportions.
computePowerSaviTwoProportions <- function(
  na,
  nb,
  propDiff = NULL,
  logOR = NULL,
  alpha = 0.05,
  priorParameters = NULL,
  nSim = 1e3,
  maxBlocks,
  nTheta = 8L,
  nBoot = 1e3,
  nWeight = 1000L
) {
  if (!is.null(nBoot) &&
      (length(nBoot) != 1L || !is.finite(nBoot) ||
       nBoot < 2L || nBoot %% 1 != 0)) {
    stop("nBoot must be an integer of at least two.")
  }

  thetaGrid <- makeSimulationThetaGrid(
    propDiff = propDiff,
    logOR = logOR,
    nTheta = nTheta
  )
  gridSimulation <- sampleStoppingTimesSaviTwoProportions(
    thetaA = thetaGrid[["thetaA"]],
    thetaB = thetaGrid[["thetaB"]],
    na = na,
    nb = nb,
    alpha = alpha,
    priorParameters = priorParameters,
    restriction = thetaGrid[["restriction"]],
    delta = thetaGrid[["delta"]],
    nWeight = nWeight,
    nSim = nSim,
    maxBlocks = maxBlocks
  )
  rejected <- is.finite(gridSimulation[["stoppingTimes"]])
  power <- rowMeans(rejected)
  worstCaseIndex <- which.min(power)
  worstCasePowerTwoSe <- NULL
  if (!is.null(nBoot)) {
    bootstrapPower <- replicate(nBoot, {
      mean(sample(
        rejected[worstCaseIndex, ],
        size = ncol(rejected),
        replace = TRUE
      ))
    })
    worstCasePowerTwoSe <- 2 * stats::sd(bootstrapPower)
  }

  c(list(
    propDiff = propDiff,
    logOR = logOR,
    thetaA = thetaGrid[["thetaA"]],
    thetaB = thetaGrid[["thetaB"]]
  ), list(
    power = power,
    worstCasePower = power[worstCaseIndex],
    worstCasePowerTwoSe = worstCasePowerTwoSe,
    worstCaseIndex = worstCaseIndex,
    worstCaseThetaA = thetaGrid[["thetaA"]][worstCaseIndex],
    worstCaseThetaB = thetaGrid[["thetaB"]][worstCaseIndex]
  ))
}

# Use monotonicity to binary-search a descending positive-effect grid and
# return the smallest value whose worst-case power reaches 1 - beta.
#
# TODO: the search assumes worst-case power is monotone in the effect, but
# powerAtIndex() draws a fresh simulation at every candidate, so Monte Carlo
# noise can invert the comparison near the boundary and send the search down
# the wrong half. Options are a common random number stream across candidates,
# or requiring the power gap to exceed its own standard error before branching.
# Left as is for now; the returned effect carries no error estimate.
computeMinEsSaviTwoProportions <- function(
  na,
  nb,
  alpha,
  beta,
  effectMeasure = c("propDiff", "logOR"),
  priorParameters = NULL,
  nSim = 1e3,
  maxBlocks,
  effectGridSize = 10L,
  effectMax = NULL,
  effectMin = 0.01,
  nTheta = 8L,
  nWeight = 1000L
) {
  effectMeasure <- match.arg(effectMeasure)
  if (length(beta) != 1L || !is.finite(beta) || beta <= 0 || beta >= 1) {
    stop("beta must be strictly between 0 and 1.")
  }
  if (length(effectGridSize) != 1L || !is.finite(effectGridSize) ||
      effectGridSize < 1L || effectGridSize %% 1 != 0) {
    stop("effectGridSize must be a positive integer.")
  }
  if (is.null(effectMax)) {
    effectMax <- if (effectMeasure == "propDiff") 0.99 else 5
  }
  if (any(!is.finite(c(effectMin, effectMax))) || effectMin <= 0 ||
      effectMin > effectMax ||
      (effectMeasure == "propDiff" && effectMax >= 1)) {
    stop(paste(
      "effectMin and effectMax must be positive and ordered;",
      "proportion differences must be less than one."
    ))
  }

  effectGrid <- seq(effectMax, effectMin, length.out = effectGridSize)
  effectArgument <- if (effectMeasure == "propDiff") {
    "propDiff"
  } else {
    "logOR"
  }
  targetPower <- 1 - beta
  powerAtIndex <- function(effectIndex) {
    simulationArguments <- list(
      na = na,
      nb = nb,
      alpha = alpha,
      priorParameters = priorParameters,
      nSim = nSim,
      maxBlocks = maxBlocks,
      nTheta = nTheta,
      nBoot = NULL,
      nWeight = nWeight
    )
    simulationArguments[[effectArgument]] <- effectGrid[effectIndex]
    do.call(
      computePowerSaviTwoProportions,
      simulationArguments
    )[["worstCasePower"]]
  }

  largestEffectPower <- powerAtIndex(1L)
  if (largestEffectPower < targetPower) {
    stop("No candidate effect reached the requested power at this horizon.")
  }
  if (effectGridSize == 1L || effectMin == effectMax) {
    return(effectGrid[1L])
  }

  smallestEffectPower <- powerAtIndex(effectGridSize)
  if (smallestEffectPower >= targetPower) {
    return(effectGrid[effectGridSize])
  }

  passingIndex <- 1L
  failingIndex <- effectGridSize
  while (failingIndex - passingIndex > 1L) {
    middleIndex <- (passingIndex + failingIndex) %/% 2L
    if (powerAtIndex(middleIndex) >= targetPower) {
      passingIndex <- middleIndex
    } else {
      failingIndex <- middleIndex
    }
  }

  effectGrid[passingIndex]
}

# Design ----

#' Designs a Savi Experiment to Test Two Proportions in Stream Data
#'
#' The Turner design uses independent Beta predictive distributions for the two
#' Bernoulli streams. The data-generating effect `delta` is either
#' `propDiff = thetaB - thetaA` or
#' `logOR = logit(thetaB) - logit(thetaA)`; it is not supplied to the
#' unrestricted Turner e-process.
#'
#' Supply `delta` and `beta` to estimate the required number of data blocks.
#' Supply `nBlocksPlan` and `delta` to estimate the worst-case type-II error at
#' that horizon, supply `nBlocksPlan` and `beta` to search a fixed grid for the
#' minimum detectable effect, or supply only `nBlocksPlan` to construct a pilot
#' design.
#'
#' Whenever a `delta` is known or found, the design simulation restricts the
#' e-process numerator to that effect and records it as the design's
#' `alternativeRestriction`; a pilot design leaves the numerator unrestricted.
#'
#' @param na Number of observations in group A per data block.
#' @param nb Number of observations in group B per data block.
#' @param nBlocksPlan Planned number of data blocks collected.
#' @param beta Numeric in `(0, 1)` specifying the tolerable type II error.
#' @param delta Minimal relevant effect used to generate data in the design
#'   simulation.
#' @param effectMeasure Whether `delta` is a proportion difference (`"propDiff"`)
#'   or log odds ratio (`"logOR"`). This is the only entry point that still
#'   accepts the legacy spellings `"difference"`, `"linearDifference"`, and
#'   `"logOddsRatio"`; they are mapped on entry and never stored or passed on.
#' @param alpha Numeric in `(0, 1)` specifying the tolerable type I error.
#' @param pilot Logical specifying whether this is a pilot design.
#' @param hyperParameterValues Named list containing positive `betaA1`,
#'   `betaA2`, `betaB1`, and `betaB2` values for the two Beta priors.
#' @param previousSaviTestResult Optional previous test result whose
#'   `posteriorHyperParameters` are used as the new prior parameters.
#' @param nSim Number of simulated paths at each baseline probability used to
#'   estimate `nBlocksPlan`, `beta`, or the minimum detectable `delta`.
#' @param nBoot Number of nonparametric bootstrap samples used for the
#'   standard errors reported alongside `nBlocksPlan` and `beta`.
#' @param maxBlocks Simulation horizon used when estimating `nBlocksPlan` from
#'   `delta` and `beta`. Raise it if `nBlocksPlan` comes back infinite. The
#'   other design cases take their horizon from `nBlocksPlan` and ignore this.
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
  nSim = 1e3L,
  nBoot = 1e3L,
  effectMeasure = c(
    "propDiff", "logOR", "difference", "linearDifference", "logOddsRatio"
  ),
  maxBlocks = 1e4
) {
  # The only place legacy effect names are accepted. Everything downstream,
  # including the stored effectMeasure, speaks propDiff and logOR only.
  effectMeasure <- switch(
    match.arg(effectMeasure),
    difference = ,
    linearDifference = "propDiff",
    logOddsRatio = "logOR",
    match.arg(effectMeasure)
  )
  # Validate inputs
  if (length(na) != 1L || length(nb) != 1L ||
      any(!is.finite(c(na, nb))) || any(c(na, nb) <= 0) ||
      any(c(na, nb) %% 1 != 0)) {
    stop("na and nb must be positive integer block sizes.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }
  if (length(nSim) != 1L || !is.finite(nSim) || nSim < 1L || nSim %% 1 != 0) {
    stop("nSim must be a positive integer.")
  }
  if (length(nBoot) != 1L || !is.finite(nBoot) ||
      nBoot < 2L || nBoot %% 1 != 0) {
    stop("nBoot must be an integer of at least two.")
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

  hyperParameterValues <- resolveBetaPriorParameters(
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
       (effectMeasure == "propDiff" && abs(delta) >= 1))) {
    stop("delta must be finite and nonzero; propDiff must be between -1 and 1.")
  }
  effectArgument <- if (effectMeasure == "propDiff") {
    "propDiff"
  } else {
    "logOR"
  }
  if (hasBeta &&
      (length(beta) != 1L || !is.finite(beta) || beta <= 0 || beta >= 1)) {
    stop("beta must be strictly between 0 and 1.")
  }

  designScenario <- NULL

  # Case 1: nBlocksPlan only -> pilot design
  if (hasNBlocksPlan && !hasEffect && !hasBeta) {
    designScenario <- "1"
    pilot <- TRUE
    warning("no simulation done!")


  # Case 2: delta and beta -> worst-case nBlocksPlan
  # bootstrapping to calculate the standard deviation
  } else if (!hasNBlocksPlan && hasEffect && hasBeta) {
    designScenario <- "1a"
    simulationArguments <- list(
      na = na,
      nb = nb,
      priorParameters = hyperParameterValues,
      alpha = alpha,
      beta = beta,
      nSim = nSim,
      maxBlocks = maxBlocks
    )
    simulationArguments[[effectArgument]] <- delta
    simulationResult <- do.call(
      computeNPlanSaviTwoProportions,
      simulationArguments
    )
    # A quantile of integer stopping times is interpolated, so round up to a
    # whole number of blocks rather than planning a fractional horizon.
    nBlocksPlan <- ceiling(simulationResult[["worstCaseStoppingTime"]])
    nPlanTwoSe <- c(0, 0, simulationResult[["worstCaseStoppingTimeTwoSe"]])
    if (!is.finite(nBlocksPlan)) {
      note <- c(
        note,
        "nBlocksPlan is infinite: raise maxBlocks or accept a larger delta."
      )
    }


  # Case 3: nBlocksPlan and delta -> worst-case power
  } else if (hasNBlocksPlan && hasEffect && !hasBeta) {
    designScenario <- "2"
    simulationArguments <- list(
      na = na,
      nb = nb,
      alpha = alpha,
      priorParameters = hyperParameterValues,
      nSim = nSim,
      maxBlocks = nBlocksPlan
    )
    simulationArguments[[effectArgument]] <- delta
    simulationResult <- do.call(computePowerSaviTwoProportions, simulationArguments)
    beta <- 1 - simulationResult[["worstCasePower"]]
    betaTwoSe <- simulationResult[["worstCasePowerTwoSe"]]
    note <- c(
      note,
      "Implied target is not calculated for the stop-on-crossing simulator."
    )

  # Case 4: nBlocksPlan and beta -> minimum detectable delta
  } else if (hasNBlocksPlan && !hasEffect && hasBeta) {
    designScenario <- "3"
    delta <- computeMinEsSaviTwoProportions(
      na = na,
      nb = nb,
      alpha = alpha,
      beta = beta,
      effectMeasure = effectMeasure,
      priorParameters = hyperParameterValues,
      nSim = nSim,
      maxBlocks = nBlocksPlan
    )
    note <- c(
      note,
      "Effect selected by monotone binary search over a simulation grid."
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
    names(delta) <- if (effectMeasure == "propDiff") {
      "propDiff"
    } else {
      "logOR"
    }
  }

  result <- constructSaviDesignObj("Two Proportions")

  result[["nPlan"]] <- nPlan
  result[["nPlanTwoSe"]] <- nPlanTwoSe
  result[["parameter"]] <- priorValuesForPrint
  result[["betaPriorParameterValues"]] <- hyperParameterValues
  result[["alpha"]] <- alpha
  result[["beta"]] <- beta
  result[["betaTwoSe"]] <- betaTwoSe
  result[["logImpliedTarget"]] <- logImpliedTarget
  result[["logImpliedTargetTwoSe"]] <- logImpliedTargetTwoSe
  result[["esMin"]] <- delta
  result[["effectMeasure"]] <- effectMeasure
  result[["h0"]] <- c("propDiff" = 0)
  result[["testType"]] <- "2x2"
  result[["alternative"]] <- "twoSided"
  result[["alternativeRestriction"]] <-
    if (is.null(delta)) "none" else effectMeasure
  result[["eType"]] <- "turner"
  result[["designScenario"]] <- designScenario
  result[["pilot"]] <- pilot
  result[["call"]] <- sys.call()
  result[["timeStamp"]] <- Sys.time()
  result[["note"]] <- note

  # Drop the slots this design never fills, as the other design functions do,
  # so printing and plotting only see what was actually computed.
  result <- Filter(Negate(is.null), result)
  class(result) <- "saviDesign"

  return(result)
}


# Calculation ----

#' Savi Test of Two Proportions in Stream Data
#'
#' Computes the Turner e-process for two Bernoulli streams observed in blocks,
#' optionally together with an anytime-valid confidence sequence for the
#' effect. The e-process is cumulative: element `i` uses blocks `1` through
#' `i`, so it may be monitored continuously and stopped as soon as it reaches
#' `1 / alpha`.
#'
#' The `designObj` supplies the Beta prior, the type-I error, and whether the
#' numerator is restricted to a fixed effect. A design carrying an `esMin`
#' restricts the numerator to that effect; a pilot design leaves it
#' unrestricted.
#'
#' @param ya,yb Number of successes in groups A and B in each data block. For
#'   the formula method, `ya` is replaced by `formula` and `data`.
#' @param na,nb Number of observations in groups A and B in each data block.
#'   A scalar is recycled over all blocks.
#' @param designObj A `saviDesign` object returned by
#'   [designSaviTwoProportions()]. If omitted, a pilot design is constructed
#'   and a warning is issued.
#' @param ciValue Numeric in `(0, 1)`, the coverage of the reported confidence
#'   sequence. Defaults to `1 - alpha` of the design.
#' @param wantConfidenceSequence Logical. If `TRUE`, a confidence sequence for
#'   the design's effect measure is computed at every block.
#' @param confidenceBoundGridPrecision Number of confidence-sequence candidate
#'   values on each side of zero. See
#'   [computeConfidenceSequenceForPropDiffTwoProportions()] and
#'   [computeConfidenceSequenceForLogORTwoProportions()].
#' @param logORConfidenceSearchBounds Positive, increasing search bounds for
#'   the absolute candidate log odds ratios. Used only when the design's
#'   effect measure is `"logOR"`.
#' @param ... Further arguments passed to methods.
#'
#' @return A `saviTest` object. Alongside the usual fields it carries
#'   `eValueVec`, the cumulative e-process; `confSeqMatrix`, the running
#'   bounds when requested; and `posteriorHyperParameters`, the updated Beta
#'   parameters, which [designSaviTwoProportions()] accepts as
#'   `previousSaviTestResult` to carry the prior forward into a new design.
#'
#' @export
#'
#' @examples
#' designObj <- designSaviTwoProportions(na = 1, nb = 1, nBlocksPlan = 20)
#' ya <- c(0, 1, 0, 0, 1, 0, 0, 1, 0, 0)
#' yb <- c(1, 1, 0, 1, 1, 1, 0, 1, 1, 0)
#' saviTwoProportionsTest(ya, yb, designObj = designObj)
saviTwoProportionsTest <- function(ya, ...) {
  UseMethod("saviTwoProportionsTest")
}

#' @rdname saviTwoProportionsTest
#' @aliases saviTwoProportionsTest
#' @export
saviTwoProportionsTest.default <- function(
  ya,
  yb,
  na = 1,
  nb = 1,
  designObj = NULL,
  ciValue = NULL,
  wantConfidenceSequence = FALSE,
  confidenceBoundGridPrecision = 20,
  logORConfidenceSearchBounds = c(0.01, 5),
  ...
) {
  result <- constructSaviTestObj("Two Proportions")

  nSteps <- length(ya)
  if (length(na) == 1L) na <- rep(na, nSteps)
  if (length(nb) == 1L) nb <- rep(nb, nSteps)

  ### Check: designObj ----
  if (is.null(designObj)) {
    designObj <- suppressWarnings(
      designSaviTwoProportions(na = na[1], nb = nb[1], nBlocksPlan = nSteps)
    )
    designObj[["pilot"]] <- TRUE
    warning(paste(
      "No designObj given. Default pilot test computed with an unrestricted",
      "numerator and the standard Beta hyperparameters."
    ))
  }

  if (!identical(designObj[["testName"]], "Two Proportions")) {
    stop("designObj must be a two-proportion design.")
  }
  if (na[1] != designObj[["nPlan"]][["na"]] ||
      nb[1] != designObj[["nPlan"]][["nb"]]) {
    warning(paste(
      "Block sizes differ from the design's na and nb;",
      "the design's Beta prior may no longer be appropriate."
    ))
  }

  alpha <- designObj[["alpha"]]
  if (is.null(ciValue)) ciValue <- 1 - alpha

  ### Calculate: e-process ----
  restriction <- designObj[["alternativeRestriction"]]
  if (is.null(restriction)) restriction <- "none"
  delta <- if (restriction == "none") NULL else unname(designObj[["esMin"]])

  eValueVec <- turnerEProcess(
    ya = ya,
    yb = yb,
    na = na,
    nb = nb,
    priorParameters = designObj[["betaPriorParameterValues"]],
    restriction = restriction,
    delta = delta
  )

  ### Calculate: confidence sequence ----
  effectMeasure <- designObj[["effectMeasure"]]
  if (is.null(effectMeasure)) effectMeasure <- "propDiff"

  confSeqMatrix <- NULL
  confSeq <- NULL
  if (wantConfidenceSequence) {
    confidenceSequence <- if (effectMeasure == "logOR") {
      computeConfidenceSequenceForLogORTwoProportions(
        ya = ya,
        yb = yb,
        confidenceBoundGridPrecision = confidenceBoundGridPrecision,
        logORConfidenceSearchBounds = logORConfidenceSearchBounds,
        saviDesign = designObj
      )
    } else {
      computeConfidenceSequenceForPropDiffTwoProportions(
        ya = ya,
        yb = yb,
        confidenceBoundGridPrecision = confidenceBoundGridPrecision,
        saviDesign = designObj
      )
    }
    confSeqMatrix <- as.matrix(
      confidenceSequence[, c("lowerBound", "upperBound")]
    )
    dimnames(confSeqMatrix) <- list(NULL, c("lowerBound", "upperBound"))
    confSeq <- confSeqMatrix[nSteps, ]
  }

  ### Calculate: summary statistics ----
  totalYa <- sum(ya)
  totalYb <- sum(yb)
  totalNa <- sum(na)
  totalNb <- sum(nb)
  observedThetaA <- totalYa / totalNa
  observedThetaB <- totalYb / totalNb

  # Reported on the design's own scale, so the estimate, the confidence
  # sequence and esMin are all the same quantity.
  observedEffect <- if (effectMeasure == "logOR") {
    stats::qlogis(observedThetaB) - stats::qlogis(observedThetaA)
  } else {
    observedThetaB - observedThetaA
  }
  estimate <- c(observedThetaA, observedThetaB, observedEffect)
  names(estimate) <- c("prop in group A", "prop in group B", effectMeasure)

  nObs <- c(totalNa, totalNb, nSteps)
  names(nObs) <- c("na", "nb", "nBlocks")

  priorParameters <- designObj[["betaPriorParameterValues"]]

  ### Fill: Result -----
  result[["eValue"]] <- unname(eValueVec[nSteps])
  result[["eValueVec"]] <- unname(eValueVec)
  result[["estimate"]] <- estimate
  result[["n"]] <- nObs
  result[["confSeq"]] <- confSeq
  result[["confSeqMatrix"]] <- confSeqMatrix
  result[["ciValue"]] <- if (wantConfidenceSequence) ciValue else NULL
  result[["designObj"]] <- designObj
  result[["testType"]] <- "2x2"
  result[["alternative"]] <- designObj[["alternative"]]
  result[["h0"]] <- designObj[["h0"]]
  result[["dataName"]] <- "ya and yb"
  result[["call"]] <- sys.call()

  # The Beta posterior after all observed blocks, so a follow-up design can
  # carry this analysis' information forward through previousSaviTestResult.
  result[["posteriorHyperParameters"]] <- list(
    betaA1 = priorParameters[["betaA1"]] + totalYa,
    betaA2 = priorParameters[["betaA2"]] + totalNa - totalYa,
    betaB1 = priorParameters[["betaB1"]] + totalYb,
    betaB2 = priorParameters[["betaB2"]] + totalNb - totalYb
  )
  result[["sumStats"]] <- list(
    ya = ya, yb = yb, na = na, nb = nb,
    totalYa = totalYa, totalYb = totalYb,
    totalNa = totalNa, totalNb = totalNb
  )

  return(result)
}

#' @rdname saviTwoProportionsTest
#' @aliases saviTwoProportionsTest
#'
#' @param formula A formula of the form `success ~ group`, where `success` is
#'   a numeric 0/1 or logical outcome and `group` is a factor with two levels.
#'   The first level is group A. Observations are paired into blocks of one
#'   observation per group in the order they appear.
#' @param data A data frame containing the variables in `formula`.
#' @param subset,na.action Passed on to [stats::model.frame()].
#'
#' @export
saviTwoProportionsTest.formula <- function(
  formula,
  data,
  subset,
  na.action,
  ...
) {
  if (missing(formula) || length(formula) != 3L) {
    stop("formula must be of the form 'success ~ group'.")
  }
  modelFrameCall <- match.call(expand.dots = FALSE)
  matched <- match(
    c("formula", "data", "subset", "na.action"),
    names(modelFrameCall),
    0L
  )
  modelFrameCall <- modelFrameCall[c(1L, matched)]
  modelFrameCall[[1L]] <- quote(stats::model.frame)
  modelFrameCall[["na.action"]] <- quote(stats::na.omit)
  modelFrame <- eval(modelFrameCall, parent.frame())

  dataName <- paste(names(modelFrame), collapse = " by ")
  groupingFactor <- factor(modelFrame[[2L]])
  if (nlevels(groupingFactor) != 2L) {
    stop("The grouping factor must have exactly two levels.")
  }
  response <- modelFrame[[1L]]
  if (is.logical(response)) response <- as.integer(response)
  if (!all(response %in% c(0, 1))) {
    stop("The response must be a 0/1 or logical outcome.")
  }

  streams <- split(response, groupingFactor)
  # Blocks pair one observation from each group, so a trailing imbalance has
  # no block to belong to.
  nBlocks <- min(lengths(streams))
  if (any(lengths(streams) != nBlocks)) {
    warning(paste0(
      "Groups have unequal sizes; using the first ", nBlocks,
      " observations of each."
    ))
  }

  result <- saviTwoProportionsTest(
    ya = streams[[1L]][seq_len(nBlocks)],
    yb = streams[[2L]][seq_len(nBlocks)],
    na = 1,
    nb = 1,
    ...
  )
  effectName <- names(result[["estimate"]])[3L]
  names(result[["estimate"]]) <- c(
    paste("prop in group", levels(groupingFactor)),
    paste0(effectName, " (", levels(groupingFactor)[2L], " - ",
           levels(groupingFactor)[1L], ")")
  )
  result[["dataName"]] <- dataName
  return(result)
}

#' @rdname saviTwoProportionsTest
#' @aliases saviTwoProportionsTest
#' @export
savi.prop.test <- function(ya, yb, na = 1, nb = 1, designObj = NULL, ...) {
  result <- saviTwoProportionsTest(
    ya = ya, yb = yb, na = na, nb = nb, designObj = designObj, ...
  )
  argumentNames <- getArgs()
  result[["dataName"]] <- paste(
    extractNameFromArgs(argumentNames, "ya"),
    "and",
    extractNameFromArgs(argumentNames, "yb")
  )
  return(result)
}
