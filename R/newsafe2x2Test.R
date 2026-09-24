# Testing fnts ----

# Cumulative likelihood-ratio for two Bernoulli streams.
savi2x2TestStat <- function(ya, yb, na, nb,
                            numeratorThetaA, numeratorThetaB,
                            denominatorThetaA, denominatorThetaB,
                            log = FALSE, ...) {
  successesA <- ya * (log(numeratorThetaA) - log(denominatorThetaA))
  successesB <- yb * (log(numeratorThetaB) - log(denominatorThetaB))

  failuresA <- (na - ya) * (log1p(-numeratorThetaA) - log1p(-denominatorThetaA))
  failuresB <- (nb - yb) * (log1p(-numeratorThetaB) - log1p(-denominatorThetaB))

  logEValueVec <- cumsum(successesA + failuresA + successesB + failuresB)

  if (log) logEValueVec else exp(logEValueVec)
}

# Cumulative log e-value against thetaA = thetaB.
logEProcess2x2PlugIn <- function(ya, yb, na, nb, priorHyperParameters,
                                 propDiff = NULL) {
  thetas <- if (is.null(propDiff)) {
    predictiveThetas2x2(ya, yb, na, nb, priorHyperParameters)
  } else {
    predictiveThetas2x2PropDiff(ya, yb, na, nb, priorHyperParameters,
      propDiff = propDiff
    )
  }
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]
  thetaNull <- (na * thetaA + nb * thetaB) / (na + nb)

  savi2x2TestStat(
    ya = ya, yb = yb, na = na, nb = nb,
    numeratorThetaA = thetaA, numeratorThetaB = thetaB,
    denominatorThetaA = thetaNull, denominatorThetaB = thetaNull,
    log = TRUE
  )
}

# Conditional e-factor of ONE table at a fixed logOR (B minus A) against the
# null nullLogOR, given the block's total successes ya + yb: under either
# value yb is Fisher's noncentral hypergeometric with odds exp(logOR) on group
# B (central at 0), and the log ratio of the two is
#   yb * (logOR - nullLogOR) - fnchLogPartition(logOR) + fnchLogPartition(nullLogOR),
# the last term being lchoose(na + nb, ya + yb) when nullLogOR = 0. Weighting
# yb, not ya, is what makes logOR the log odds of B over A.
conditionalEValueFixedAlternative <- function(ya, yb, na, nb, logOR,
                                              nullLogOR = 0, log = FALSE) {
  totalSuccesses <- ya + yb
  logEFactor <- yb * (logOR - nullLogOR) -
    fnchLogPartition(nb, na, totalSuccesses, logOR) +
    fnchLogPartition(nb, na, totalSuccesses, nullLogOR)

  if (log) logEFactor else exp(logEFactor)
}

# Per-block conditional e-factors (not cumulated) at a plug-in logOR.
# 1. given a logOR, just plug in: GROW case or UMP case
# 2. given a prior weight, numerically integrate it (report posterior maybe later)
savi2x2CondStat <- function(ya, yb, na, nb, logOR = NULL, parameter,
                            alternative = c("twoSided", "less", "greater"),
                            eType = c("grow", "ump", "eGauss"),
                            log = FALSE, ...) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  if (is.null(logOR)) {
    stop("savi2x2CondStat() needs logOR")
  }
  if (length(logOR) != 1 || !is.finite(logOR)) {
    stop("logOR must be one finite number")
  }

  # TODO: dispatch on eType and alternative (not designed yet):
  # - grow + less/greater: plug in the given logOR, one e-factor per block;
  # - ump: only for one block, logOR from solveUmpLogOR();
  # - eGauss + twoSided: integrate the e-factor over a prior on logOR.
  # For now the plug-in case. fnchLogPartition() is scalar, hence mapply().
  logEFactorVec <- mapply(
    conditionalEValueFixedAlternative,
    ya = ya, yb = yb, na = na, nb = nb,
    MoreArgs = list(logOR = logOR, log = TRUE)
  )

  if (log) logEFactorVec else exp(logEFactorVec)
}


#' Safe anytime-valid 2x2 test
#'
#' Tests the equality null `thetaA = thetaB` on data that arrive in blocks of
#' `na[i]` observations from group A and `nb[i]` from group B. The numerator predicts each block with the Beta
#' posterior means of `thetaA` and `thetaB` given the earlier blocks only; the
#' denominator uses their size-weighted average, the projection of that
#' prediction onto the null.
#'
#' With `propDiffMin` set in the design, the numerator is restricted to the
#' curve `thetaB - thetaA = propDiffMin` and learns `thetaA` along it from
#' the earlier blocks. For `alternative = "greater"` that is the e-process;
#' for `"twoSided"` it is the average of the e-processes restricted at
#' `+propDiffMin` and `-propDiffMin`. The null is always `thetaA = thetaB`.
#'
#' @param ya,yb integer vectors, the number of successes in group A and group
#'   B in each block, in the order the blocks were observed.
#' @param na,nb the group sizes per block: `NULL` takes the planned sizes
#'   from the design, one positive integer is used for every block, and a
#'   vector of length `length(ya)` gives each block its own size.
#' @param designObj a `saviDesign` object from `designSavi2x2()`. `NULL`
#'   gives a pilot design with the default settings and a warning.
#' @param wantCi logical, whether to compute the anytime-valid confidence
#'   sequence for `propDiff`; its coverage is `1 - alpha` from the design.
#' @param runningIntersection logical, whether a candidate `propDiff` that is
#'   rejected once stays rejected in all later blocks (`TRUE`, nested sets),
#'   or each block reports the candidates its current e-value has not
#'   rejected (`FALSE`).
#'
#' @return A `saviTest` object. `eValueVec[i]` is the cumulative e-process
#'   after block `i`, using blocks `1` to `i` only, and `eValue` is its last
#'   element. `estimate` holds the observed proportions and their difference
#'   B minus A. `n` holds the total group sizes and the number of blocks.
#'   `posteriorHyperParameters`
#'   holds the Beta posterior after the last block, which is the prior the
#'   next block would use.
#' @noRd
savi2x2Test <- function(ya, yb, na = NULL, nb = NULL, designObj = NULL,
                        wantCi = TRUE, runningIntersection = TRUE) {
  result <- constructSaviTestObj("Two Proportions")

  if (is.null(designObj)) {
    designObj <- designSavi2x2()
    designObj[["pilot"]] <- TRUE
    warning("No designObj given. Default pilot design used.")
  }

  if (!identical(designObj[["testName"]], "Two Proportions")) {
    stop("designObj must be a design from designSavi2x2().")
  }

  propDiffMin <- designObj[["esMin"]]
  alternative <- designObj[["alternative"]]

  if (alternative == "less") {
    stop("alternative = \"less\" is not implemented yet.")
  }

  if (alternative == "greater" && is.null(propDiffMin)) {
    stop("alternative = \"greater\" needs a positive propDiffMin.")
  }

  # The e-variable below is built for thetaA = thetaB only; a shifted null
  # thetaB - thetaA = h0 needs a different projection.
  if (designObj[["h0"]] != 0) {
    stop("Only h0 = 0 is implemented yet.")
  }

  nBlocks <- length(ya)

  if (nBlocks < 1L || length(yb) != nBlocks) {
    stop("ya and yb must have the same, positive length.")
  }

  if (!is.numeric(ya) || !is.numeric(yb) ||
    any(!is.finite(c(ya, yb))) || any(c(ya, yb) %% 1 != 0)) {
    stop("ya and yb must contain finite integer counts.")
  }

  # Observed block sizes: the planned size from the design, one size for all
  # blocks, or one size per block. They may differ from the plan.
  if (is.null(na)) na <- designObj[["nPlan"]][["na"]]
  if (is.null(nb)) nb <- designObj[["nPlan"]][["nb"]]
  if (length(na) == 1L) na <- rep(na, nBlocks)
  if (length(nb) == 1L) nb <- rep(nb, nBlocks)

  if (length(na) != nBlocks || length(nb) != nBlocks) {
    stop("na and nb must be NULL, one number, or have one entry per block.")
  }

  if (!is.numeric(na) || !is.numeric(nb) ||
    any(!is.finite(c(na, nb))) || any(c(na, nb) < 1) ||
    any(c(na, nb) %% 1 != 0)) {
    stop("na and nb must contain positive integer block sizes.")
  }

  if (any(ya < 0) || any(yb < 0) || any(ya > na) || any(yb > nb)) {
    stop(
      "Success counts must lie between zero and the block sizes ",
      "na and nb of their block."
    )
  }

  prior <- designObj[["priorHyperParameters"]]

  if (is.null(propDiffMin)) {
    logEValueVec <- logEProcess2x2PlugIn(ya, yb, na, nb, prior)
  } else if (alternative == "greater") {
    logEValueVec <- logEProcess2x2PlugIn(ya, yb, na, nb, prior,
      propDiff = propDiffMin
    )
  } else {
    # Two-sided with a restriction: the equal-weight mixture of the two
    # cumulative e-processes at +propDiffMin and -propDiffMin, which is again
    # an e-process. Averaged on the log scale to avoid overflow.
    logEPlus <- logEProcess2x2PlugIn(ya, yb, na, nb, prior,
      propDiff = propDiffMin
    )
    logEMinus <- logEProcess2x2PlugIn(ya, yb, na, nb, prior,
      propDiff = -propDiffMin
    )
    logEValueVec <- pmax(logEPlus, logEMinus) +
      log1p(exp(-abs(logEPlus - logEMinus))) - log(2)
  }

  eValueVec <- exp(logEValueVec)

  ## Confidence Sequence ----
  if (wantCi) {
    alpha <- designObj[["alpha"]]
    confSetRuns <- computeConfidenceInterval2x2PropDiff(
      ya = ya, yb = yb, na = na, nb = nb,
      priorHyperParameters = prior, alpha = alpha,
      runningIntersection = runningIntersection
    )

    # One row per block, as for the other tests: the outermost bounds of that
    # block's union, NA when every candidate is rejected. The hull contains
    # the union, so coverage is kept; only the last block is kept exact.
    block <- factor(confSetRuns[, "block"], levels = seq_len(nBlocks))
    confSeqMatrix <- cbind(
      "lowerBound" = as.vector(tapply(confSetRuns[, "lowerBound"], block, min)),
      "upperBound" = as.vector(tapply(confSetRuns[, "upperBound"], block, max))
    )

    lastBlock <- confSetRuns[, "block"] == nBlocks

    result[["confSeqMatrix"]] <- confSeqMatrix
    result[["confSeq"]] <- confSetRuns[lastBlock, c("lowerBound", "upperBound"),
      drop = FALSE
    ]
    result[["ciValue"]] <- 1 - alpha
  }

  result[["eValue"]] <- eValueVec[nBlocks]
  result[["eValueVec"]] <- eValueVec
  result[["n"]] <- c("na" = sum(na), "nb" = sum(nb), "nBlocks" = nBlocks)
  result[["posteriorHyperParameters"]] <- list(
    "betaA1" = prior[["betaA1"]] + sum(ya),
    "betaA2" = prior[["betaA2"]] + sum(na) - sum(ya),
    "betaB1" = prior[["betaB1"]] + sum(yb),
    "betaB2" = prior[["betaB2"]] + sum(nb) - sum(yb)
  )
  result[["designObj"]] <- designObj
  result[["testType"]] <- "2x2"
  result[["alternative"]] <- designObj[["alternative"]]
  result[["h0"]] <- designObj[["h0"]]
  result[["dataName"]] <- paste(
    deparse1(substitute(ya)), "and",
    deparse1(substitute(yb))
  )
  result[["call"]] <- sys.call()

  # x-axis of plot.saviTest(): the block index.
  result[["n1Vec"]] <- seq_len(nBlocks)
  # TODO: decided later: MLE for propDiff and or MLE for logOR?
  result[["estimate"]] <- 0

  return(result)
}


# Design fnts ----

#' Design a safe anytime-valid 2x2 test
#'
#' Sets up the design object for the two-proportion test. Data arrive in
#' blocks of `na` observations from group A and `nb` from group B. No
#' sample-size planning is done yet: `nPlan[["nBlocks"]]` is `NA`.
#'
#' @param propDiffMin `NULL`, or a number strictly between 0 and 1: the
#'   minimal relevant difference `thetaB - thetaA`. When set, the e-variable's
#'   numerator is restricted to `thetaB - thetaA = propDiffMin` (and, for
#'   `"twoSided"`, also to `-propDiffMin`).
#' @param na,nb positive integers, the planned number of observations per block in
#'   group A and group B.
#' @param nPlan integer or `NULL`, the planned number of blocks. Not used yet.
#' @param alpha numeric in (0, 1), the tolerable type I error rate. The null
#'   is rejected once the e-process reaches `1/alpha`.
#' @param power numeric in (0, 1) or `NULL`, the target power. Not used yet.
#' @param h0 numeric, the difference under the null.
#' @param alternative one of `"twoSided"`, `"greater"`, `"less"`. The
#'   direction refers to the effect B minus A, so `"greater"` means group B
#'   has the larger success probability. `"greater"` requires `propDiffMin`;
#'   `"less"` is not implemented yet.
#' @param eType the e-variable type, currently only `"grow"`.
#' @param priorHyperParameters `NULL`, or a list named `betaA1`, `betaA2`,
#'   `betaB1`, `betaB2`: the success and failure shapes of the Beta priors on
#'   `thetaA` and `thetaB`. `NULL` keeps the default of 0.18 for all four,
#'   set in [constructSaviDesignObj()].
#'
#' @return A `saviDesign` object. `h0` is named `propDiff`, `nPlan` holds
#'   `nBlocks`, `na` and `nb`, and `priorHyperParameters` holds the Beta
#'   prior used by the e-variable.
#' @noRd
designSavi2x2 <- function(propDiffMin = NULL, alpha = 0.05, na = 1, nb = 1,
                          nPlan = NULL, power = NULL, h0 = 0,
                          alternative = c("twoSided", "greater", "less"),
                          eType = c("grow"),
                          priorHyperParameters = NULL) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  if (length(na) != 1L || length(nb) != 1L ||
    !is.finite(na) || !is.finite(nb) ||
    na < 1 || nb < 1 || na %% 1 != 0 || nb %% 1 != 0) {
    stop("na and nb must be positive integer block sizes.")
  }

  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }

  if (!is.null(propDiffMin) &&
    (length(propDiffMin) != 1L || !is.finite(propDiffMin) ||
      propDiffMin <= 0 || propDiffMin >= 1)) {
    stop("propDiffMin must be NULL or strictly between 0 and 1.")
  }

  if (alternative == "less") {
    stop("alternative = \"less\" is not implemented yet.")
  }

  if (alternative == "greater" && is.null(propDiffMin)) {
    stop("alternative = \"greater\" needs a positive propDiffMin.")
  }

  result <- constructSaviDesignObj("Two Proportions")

  if (!is.null(priorHyperParameters)) {
    requiredNames <- names(result[["priorHyperParameters"]])

    if (!is.list(priorHyperParameters) ||
      !setequal(names(priorHyperParameters), requiredNames)) {
      stop(
        "priorHyperParameters must be a list named ",
        paste(requiredNames, collapse = ", "), "."
      )
    }

    priorValues <- unlist(priorHyperParameters[requiredNames])

    if (!is.numeric(priorValues) || any(!is.finite(priorValues)) ||
      any(priorValues <= 0)) {
      stop("priorHyperParameters must be finite and positive.")
    }

    result[["priorHyperParameters"]] <- priorHyperParameters[requiredNames]
  }

  result[["esMin"]] <- propDiffMin
  result[["parameter"]] <- c(
    "Beta hyperparameters" =
      paste(unlist(result[["priorHyperParameters"]]), collapse = " ")
  )
  result[["eType"]] <- eType
  result[["alpha"]] <- alpha
  result[["alternative"]] <- alternative
  result[["h0"]] <- c("propDiff" = h0)
  result[["nPlan"]] <- c("nBlocks" = NA_real_, "na" = na, "nb" = nb)
  result[["testType"]] <- "2x2"
  result[["call"]] <- sys.call()
  result[["timeStamp"]] <- Sys.time()

  return(result)
}


# Confidence Interval ----


#' Anytime-valid confidence sequence for the proportion difference
#'
#' Inverts the test: each candidate `propDiff = thetaB - thetaA` on a grid is
#' a point null with its own e-process, whose numerator is the same
#' predictable Beta posterior mean as in `savi2x2Test()` and whose
#' denominator is that prediction projected onto the candidate's null line.
#' With `runningIntersection = TRUE` a candidate leaves the set for good once
#' its e-process reaches `1/alpha`, so the sets are nested over blocks; with
#' `FALSE` each block keeps the candidates whose current e-value is below
#' `1/alpha`.
#'
#' @param ya,yb integer vectors, the successes in group A and group B in each
#'   block.
#' @param na,nb integer vectors of length `length(ya)`, the block sizes.
#' @param priorHyperParameters list with `betaA1`, `betaA2`, `betaB1`,
#'   `betaB2`.
#' @param alpha numeric in (0, 1); the sequence has coverage `1 - alpha`.
#' @param precision positive integer, the number of equally spaced
#'   candidates strictly inside `(-1, 1)`.
#' @param runningIntersection logical, see above.
#'
#' @return A matrix with columns `block`, `lowerBound` and `upperBound`.
#'   Each run of consecutive non-rejected candidates after a block is one
#'   interval; the confidence set after that block is the union of its rows.
#'   Without holes there is one row per block, as for the z-test; a block
#'   whose candidates are all rejected has no row.
#' @noRd
computeConfidenceInterval2x2PropDiff <- function(ya, yb, na, nb,
                                                 priorHyperParameters,
                                                 alpha, precision = 100,
                                                 runningIntersection = TRUE) {
  nBlocks <- length(ya)
  thetas <- predictiveThetas2x2(ya, yb, na, nb, priorHyperParameters)

  propDiffGrid <- seq(-1, 1, length.out = precision + 2)[-c(1, precision + 2)]
  logEValues <- numeric(precision)
  inSet <- rep(TRUE, precision)

  confSeqMatrix <- matrix(numeric(0),
    ncol = 3,
    dimnames = list(NULL, c(
      "block", "lowerBound",
      "upperBound"
    ))
  )

  for (i in seq_len(nBlocks)) {
    # Under the running intersection a rejected candidate never returns, so
    # its e-process is not advanced; otherwise every candidate is followed.
    activeCandidates <- if (runningIntersection) which(inSet) else seq_len(precision)

    for (j in activeCandidates) {
      propDiff <- propDiffGrid[j]
      nullThetaA <- solveRIPr2x2PropDiff(
        thetaA = thetas[["thetaA"]][i], thetaB = thetas[["thetaB"]][i],
        na = na[i], nb = nb[i], propDiff = propDiff
      )

      logEValues[j] <- logEValues[j] + savi2x2TestStat(
        ya = ya[i], yb = yb[i], na = na[i], nb = nb[i],
        numeratorThetaA = thetas[["thetaA"]][i],
        numeratorThetaB = thetas[["thetaB"]][i],
        denominatorThetaA = nullThetaA,
        denominatorThetaB = nullThetaA + propDiff,
        log = TRUE
      )
    }

    notRejected <- logEValues < log(1 / alpha)
    inSet <- if (runningIntersection) inSet & notRejected else notRejected

    # Split the non-rejected candidates into runs of neighbours on the grid;
    # each run is one interval of the union.
    runs <- rle(inSet)
    runEnds <- cumsum(runs[["lengths"]])[runs[["values"]]]
    runStarts <- runEnds - runs[["lengths"]][runs[["values"]]] + 1

    confSeqMatrix <- rbind(
      confSeqMatrix,
      cbind(
        "block" = rep(i, length(runStarts)),
        "lowerBound" = propDiffGrid[runStarts],
        "upperBound" = propDiffGrid[runEnds]
      )
    )
  }

  return(confSeqMatrix)
}

#' Anytime-valid confidence sequence for the log odds ratio
#'
#' Inverts the conditional test: each candidate `logOR = logit(thetaB) -
#' logit(thetaA)` on a grid is the null of its own e-process, the product over
#' blocks of the conditional e-factors of
#' `conditionalEValueFixedAlternative()`. The plug-in alternative for block
#' `i` is the log odds ratio of the Beta posterior means given blocks `1` to
#' `i - 1`, so it is predictable and finite even with empty cells. A
#' candidate leaves the set for good once its e-process reaches `1/alpha`
#' when `runningIntersection = TRUE`, so the sets are nested over blocks;
#' with `FALSE` each block keeps the candidates whose current e-value is
#' below `1/alpha`.
#'
#' @param ya,yb integer vectors, the successes in group A and group B in each
#'   block.
#' @param na,nb integer vectors of length `length(ya)`, the block sizes.
#' @param priorHyperParameters list with `betaA1`, `betaA2`, `betaB1`,
#'   `betaB2`.
#' @param alpha numeric in (0, 1); the sequence has coverage `1 - alpha`.
#' @param precision positive integer, the number of equally spaced
#'   candidates strictly inside `(-logORBound, logORBound)`.
#' @param logORBound positive number, the half-width of the candidate grid.
#' @param runningIntersection logical, see above.
#'
#' @return A matrix with columns `block`, `lowerBound` and `upperBound`, one
#'   row per run of consecutive non-rejected candidates after each block, as
#'   for `computeConfidenceInterval2x2PropDiff()`.
#' @noRd
computeConfidenceInterval2x2LogOR <- function(ya, yb, na, nb,
                                              priorHyperParameters,
                                              alpha, precision = 100,
                                              logORBound = 40,
                                              runningIntersection = TRUE) {
  nBlocks <- length(ya)
  thetas <- predictiveThetas2x2(ya, yb, na, nb, priorHyperParameters)
  # Predictable plug-in alternative: B minus A on the logit scale.
  # TODO: to be decided here
  plugInLogOR <- stats::qlogis(thetas[["thetaB"]]) -
    stats::qlogis(thetas[["thetaA"]])

  logORGrid <- seq(-logORBound, logORBound,
    length.out = precision + 2
  )[-c(1, precision + 2)]
  logEValues <- numeric(precision)
  inSet <- rep(TRUE, precision)

  confSeqMatrix <- matrix(numeric(0),
    ncol = 3,
    dimnames = list(NULL, c("block", "lowerBound", "upperBound"))
  )

  for (i in seq_len(nBlocks)) {
    # Under the running intersection a rejected candidate never returns, so
    # its e-process is not advanced; otherwise every candidate is followed.
    activeCandidates <- if (runningIntersection) which(inSet) else seq_len(precision)

    for (j in activeCandidates) {
      logEValues[j] <- logEValues[j] + conditionalEValueFixedAlternative(
        ya = ya[i], yb = yb[i], na = na[i], nb = nb[i],
        logOR = plugInLogOR[i], nullLogOR = logORGrid[j], log = TRUE
      )
    }

    notRejected <- logEValues < log(1 / alpha)
    inSet <- if (runningIntersection) inSet & notRejected else notRejected

    runs <- rle(inSet)
    runEnds <- cumsum(runs[["lengths"]])[runs[["values"]]]
    runStarts <- runEnds - runs[["lengths"]][runs[["values"]]] + 1

    confSeqMatrix <- rbind(
      confSeqMatrix,
      cbind(
        "block" = rep(i, length(runStarts)),
        "lowerBound" = logORGrid[runStarts],
        "upperBound" = logORGrid[runEnds]
      )
    )
  }

  return(confSeqMatrix)
}

# Helpers ----

## propDiff ----

# Predictable plug-in for the numerator, for block i given the counts of
# blocks 1 to i - 1 only. predictiveThetas2x2(): the independent Beta
# posterior means of thetaA and thetaB. predictiveThetas2x2PropDiff(): the
# numerator lives on the curve thetaB - thetaA = propDiff, and thetaA is the
# posterior mean under a grid posterior on that curve.
predictiveThetas2x2 <- function(ya, yb, na, nb, priorHyperParameters) {
  nBlocks <- length(ya)
  betaA1 <- priorHyperParameters[["betaA1"]]
  betaA2 <- priorHyperParameters[["betaA2"]]
  betaB1 <- priorHyperParameters[["betaB1"]]
  betaB2 <- priorHyperParameters[["betaB2"]]

  # Successes and sizes accumulated over blocks 1 to i - 1.
  previousYa <- c(0, cumsum(ya))[seq_len(nBlocks)]
  previousYb <- c(0, cumsum(yb))[seq_len(nBlocks)]
  previousNa <- c(0, cumsum(na))[seq_len(nBlocks)]
  previousNb <- c(0, cumsum(nb))[seq_len(nBlocks)]

  return(list(
    "thetaA" = (betaA1 + previousYa) / (betaA1 + betaA2 + previousNa),
    "thetaB" = (betaB1 + previousYb) / (betaB1 + betaB2 + previousNb)
  ))
}

predictiveThetas2x2PropDiff <- function(ya, yb, na, nb, priorHyperParameters,
                                        propDiff, nWeight = 1000L) {
  nBlocks <- length(ya)

  # output placeholder
  thetaA <- numeric(nBlocks)

  # Only the thetaA prior is used
  betaA1 <- priorHyperParameters[["betaA1"]]
  betaA2 <- priorHyperParameters[["betaA2"]]

  if (length(propDiff) != 1L || !is.finite(propDiff) || abs(propDiff) >= 1) {
    stop("propDiff must lie strictly between -1 and 1.")
  }

  # thetaA are restricted by propDiff
  rho <- seq(1 / nWeight, 1 - 1 / nWeight, length.out = nWeight)
  thetaAGrid <- max(0, -propDiff) + rho * (1 - abs(propDiff))
  thetaBGrid <- thetaAGrid + propDiff

  logThetaA <- log(thetaAGrid)
  logOneMinusThetaA <- log1p(-thetaAGrid)
  logThetaB <- log(thetaBGrid)
  logOneMinusThetaB <- log1p(-thetaBGrid)

  # Un-normalised log posterior weights, shifted so their maximum is 0: the
  # largest weight is then exactly 1 and the sum can neither underflow nor
  # overflow, however many blocks have been seen.
  logWeights <- (betaA1 - 1) * log(rho) + (betaA2 - 1) * log1p(-rho)
  logWeights <- logWeights - max(logWeights)

  for (i in seq_len(nBlocks)) {
    # posterior mean for thetaA
    weights <- exp(logWeights)
    thetaA[i] <- sum(thetaAGrid * weights) / sum(weights)

    logWeights <- logWeights +
      ya[i] * logThetaA + (na[i] - ya[i]) * logOneMinusThetaA +
      yb[i] * logThetaB + (nb[i] - yb[i]) * logOneMinusThetaB
    logWeights <- logWeights - max(logWeights)
  }

  list("thetaA" = thetaA, "thetaB" = thetaA + propDiff)
}

# Find the means that minimize the KL between the alternative and null
# The null is H0: thetaB - thetaA = propDiff
# TODO: this can be a cubic function solver
solveRIPr2x2PropDiff <- function(thetaA, thetaB, na, nb, propDiff) {
  derivativeKL <- function(nullThetaA) {
    nullThetaB <- nullThetaA + propDiff
    na * ((1 - thetaA) / (1 - nullThetaA) - thetaA / nullThetaA) +
      nb * ((1 - thetaB) / (1 - nullThetaB) - thetaB / nullThetaB)
  }

  # The derivative is infinite at the edges, so search just inside them.
  stats::uniroot(derivativeKL,
    lower = max(0, -propDiff) + 1e-12, upper = min(1, 1 - propDiff) - 1e-12,
    tol = 1e-12
  )[["root"]]
}

## logOR ----
# Log partition function of Fisher's noncentral hypergeometric distribution
fnchLogPartition <- function(na, nb, totalSuccesses, logOR) {
  if (logOR == 0) {
    return(lchoose(na + nb, totalSuccesses))
  }

  feasibleSuccesses <-
    max(0, totalSuccesses - nb):min(na, totalSuccesses)
  logTerms <- lchoose(na, feasibleSuccesses) +
    lchoose(nb, totalSuccesses - feasibleSuccesses) +
    logOR * feasibleSuccesses

  # log-sum-exp, shifted by the largest term against overflow
  maxLogTerm <- max(logTerms)
  maxLogTerm + log(sum(exp(logTerms - maxLogTerm)))
}

# UMP plug-in for one table: the logOR on the side of the alternative at which
#   KL(FNCH(logOR) || FNCH(nullLogOR)) = log(1 / alpha),
# with KL = (logOR - nullLogOR) * E_logOR[yb]
#           - fnchLogPartition(logOR) + fnchLogPartition(nullLogOR),
# the mean E_logOR[yb] taken from BiasedUrn (odds = exp(logOR) on B). At
# nullLogOR = 0 the last term is lchoose(na + nb, ya + yb).
# The KL is 0 at nullLogOR and increases away from it, but is bounded by
# -log P0(ya at its feasible extreme), so the equation may have no root:
# then NULL is returned and the caller uses the trivial e-factor 1.
solveUmpLogOR <- function(na, nb, totalSuccesses, alpha,
                          alternative = c("greater", "less"),
                          nullLogOR = 0,
                          searchBound = 100) {
  alternative <- match.arg(alternative)

  # logOR is B minus A, so the weighted count is yb: group B goes first.
  klMinusTarget <- function(logOR) {
    (logOR - nullLogOR) *
      BiasedUrn::meanFNCHypergeo(nb, na, totalSuccesses, exp(logOR)) -
      fnchLogPartition(nb, na, totalSuccesses, logOR) +
      fnchLogPartition(nb, na, totalSuccesses, nullLogOR) + log(alpha)
  }

  # The KL is bounded, so the target may be unreachable within the search
  # interval; uniroot() would error on equal signs, hence the explicit check.
  bounds <- if (alternative == "greater") {
    c(nullLogOR, nullLogOR + searchBound)
  } else {
    c(nullLogOR - searchBound, nullLogOR)
  }
  if (klMinusTarget(bounds[1]) * klMinusTarget(bounds[2]) > 0) {
    return(NULL)
  }

  stats::uniroot(klMinusTarget,
    lower = bounds[1], upper = bounds[2],
    tol = 1e-10
  )[["root"]]
}
