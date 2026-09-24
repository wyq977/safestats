# Testing fnts ----

# Cumulative likelihood-ratio e-process for blocks of two Bernoulli streams.
# The four theta vectors are predictable: element i uses blocks 1 to i - 1
# only. Element i of the result is the e-process after block i, so a call on
# a single block returns that block's e-factor.
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

# Predictable plug-in for the numerator: the Beta posterior means of thetaA
# and thetaB for block i, given the counts of blocks 1 to i - 1 only.
predictiveThetas2x2 <- function(ya, yb, na, nb, priorHyperParameters) {
  nBlocks <- length(ya)
  previousYa <- c(0, cumsum(ya))[seq_len(nBlocks)]
  previousYb <- c(0, cumsum(yb))[seq_len(nBlocks)]
  previousBlocks <- seq_len(nBlocks) - 1

  prior <- priorHyperParameters
  list(
    "thetaA" = (prior[["betaA1"]] + previousYa) /
      (prior[["betaA1"]] + prior[["betaA2"]] + na * previousBlocks),
    "thetaB" = (prior[["betaB1"]] + previousYb) /
      (prior[["betaB1"]] + prior[["betaB2"]] + nb * previousBlocks))
}

# compute conditional e-variable
# 1. given a logOR, just plug in: GROW case or UMP case
# 2. given a prior weight, numerically integrate it (report posterior maybe later)
savi2x2CondStat <- function(ya, yb, na, nb, logOR = NULL, alternative=c("twoSided", "greater", "less")) {}



#' Safe anytime-valid 2x2 test
#'
#' Tests the equality null `thetaA = thetaB` on data that arrive in blocks of
#' `na` observations from group A and `nb` from group B, with `na` and `nb`
#' taken from the design. The numerator predicts each block with the Beta
#' posterior means of `thetaA` and `thetaB` given the earlier blocks only; the
#' denominator uses their size-weighted average, the projection of that
#' prediction onto the null. Only the unrestricted, two-sided case exists.
#'
#' @param ya,yb integer vectors, the number of successes in group A and group
#'   B in each block, in the order the blocks were observed.
#' @param designObj a `saviDesign` object from `designSavi2x2()`. `NULL`
#'   gives a pilot design with the default settings and a warning.
#'
#' @return A `saviTest` object. `eValueVec[i]` is the cumulative e-process
#'   after block `i`, using blocks `1` to `i` only, and `eValue` is its last
#'   element. `estimate` holds the observed proportions and their difference
#'   B minus A. `posteriorHyperParameters` holds the Beta posterior after the
#'   last block, which is the prior the next block would use.
#' @noRd
savi2x2Test <- function(ya, yb, designObj = NULL, wantCi = TRUE) {
  result <- constructSaviTestObj("Two Proportions")

  if (is.null(designObj)) {
    designObj <- designSavi2x2()
    designObj[["pilot"]] <- TRUE
    warning("No designObj given. Default pilot design used.")
  }

  if (!identical(designObj[["testName"]], "Two Proportions"))
    stop("designObj must be a design from designSavi2x2().")

  if (!is.null(designObj[["esMin"]]))
    stop("A restricted alternative (propDiffMin) is not implemented yet.")

  if (designObj[["alternative"]] != "twoSided")
    stop("Only alternative = \"twoSided\" is implemented yet.")

  # The e-variable below is built for thetaA = thetaB only; a shifted null
  # thetaB - thetaA = h0 needs a different projection.
  if (designObj[["h0"]] != 0)
    stop("Only h0 = 0 is implemented yet.")

  na <- designObj[["nPlan"]][["na"]]
  nb <- designObj[["nPlan"]][["nb"]]
  nBlocks <- length(ya)

  if (nBlocks < 1L || length(yb) != nBlocks)
    stop("ya and yb must have the same, positive length.")

  if (!is.numeric(ya) || !is.numeric(yb) ||
      any(!is.finite(c(ya, yb))) || any(c(ya, yb) %% 1 != 0))
    stop("ya and yb must contain finite integer counts.")

  if (any(ya < 0) || any(yb < 0) || any(ya > na) || any(yb > nb))
    stop("Success counts must lie between zero and the block sizes ",
         "na = ", na, " and nb = ", nb, ".")

  prior <- designObj[["priorHyperParameters"]]

  thetas <- predictiveThetas2x2(ya, yb, na, nb, prior)
  thetaA <- thetas[["thetaA"]]
  thetaB <- thetas[["thetaB"]]

  # Reverse information projection of the product Bernoulli prediction onto
  # thetaA = thetaB: the common theta is the size-weighted average.
  thetaNull <- (na * thetaA + nb * thetaB) / (na + nb)

  # Positive Beta shapes keep every theta strictly inside (0, 1), so no
  # 0 * log(0) term arises.
  eValueVec <- savi2x2TestStat(
    ya = ya, yb = yb, na = na, nb = nb,
    numeratorThetaA = thetaA, numeratorThetaB = thetaB,
    denominatorThetaA = thetaNull, denominatorThetaB = thetaNull)

  ## Confidence Sequence ----
  if (wantCi) {
    alpha <- designObj[["alpha"]]
    confSeqMatrix <- computeConfidenceInterval2x2PropDiff(
      ya = ya, yb = yb, na = na, nb = nb,
      priorHyperParameters = prior, alpha = alpha)

    lastBlock <- confSeqMatrix[, "block"] == nBlocks

    result[["confSeqMatrix"]] <- confSeqMatrix
    result[["confSeq"]] <- confSeqMatrix[lastBlock, c("lowerBound", "upperBound"),
                                         drop = FALSE]
    result[["ciValue"]] <- 1 - alpha
  }

  result[["eValue"]] <- eValueVec[nBlocks]
  result[["eValueVec"]] <- eValueVec
  result[["n"]] <- c("na" = na, "nb" = nb, "nBlocks" = nBlocks)
  # x-axis of plot.saviTest(): the block index.
  result[["n1Vec"]] <- seq_len(nBlocks)
  result[["estimate"]] <- 0 # TODO: decided later
  result[["posteriorHyperParameters"]] <- list(
    "betaA1" = prior[["betaA1"]] + sum(ya),
    "betaA2" = prior[["betaA2"]] + na * nBlocks - sum(ya),
    "betaB1" = prior[["betaB1"]] + sum(yb),
    "betaB2" = prior[["betaB2"]] + nb * nBlocks - sum(yb))
  result[["designObj"]] <- designObj
  result[["testType"]] <- "2x2"
  result[["alternative"]] <- designObj[["alternative"]]
  result[["h0"]] <- designObj[["h0"]]
  result[["dataName"]] <- paste(deparse1(substitute(ya)), "and",
                                deparse1(substitute(yb)))
  result[["call"]] <- sys.call()

  return(result)
}


# Design fnts ----

#' Design a safe anytime-valid 2x2 test
#'
#' Sets up the design object for the two-proportion test. Data arrive in
#' blocks of `na` observations from group A and `nb` from group B. No
#' sample-size planning is done yet: `nPlan[["nBlocks"]]` is `NA`.
#'
#' @param propDiffMin numeric or `NULL`, the minimal relevant difference
#'   `thetaB - thetaA`. Not used yet.
#' @param na,nb positive integers, the number of observations per block in
#'   group A and group B.
#' @param nPlan integer or `NULL`, the planned number of blocks. Not used yet.
#' @param alpha numeric in (0, 1), the tolerable type I error rate. The null
#'   is rejected once the e-process reaches `1/alpha`.
#' @param power numeric in (0, 1) or `NULL`, the target power. Not used yet.
#' @param h0 numeric, the difference under the null.
#' @param alternative one of `"twoSided"`, `"greater"`, `"less"`. The
#'   direction refers to the effect B minus A, so `"greater"` means group B
#'   has the larger success probability.
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
designSavi2x2 <- function(propDiffMin=NULL, na = 1, nb = 1, nPlan=NULL, alpha = 0.05, power=NULL,
                          h0=0, alternative=c("twoSided", "greater", "less"),
                          eType = c("grow"),
                          priorHyperParameters = NULL) {
  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  if (length(na) != 1L || length(nb) != 1L ||
      !is.finite(na) || !is.finite(nb) ||
      na < 1 || nb < 1 || na %% 1 != 0 || nb %% 1 != 0)
    stop("na and nb must be positive integer block sizes.")

  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("alpha must be strictly between 0 and 1.")

  result <- constructSaviDesignObj("Two Proportions")

  if (!is.null(priorHyperParameters)) {
    requiredNames <- names(result[["priorHyperParameters"]])

    if (!is.list(priorHyperParameters) ||
        !setequal(names(priorHyperParameters), requiredNames))
      stop("priorHyperParameters must be a list named ",
           paste(requiredNames, collapse = ", "), ".")

    priorValues <- unlist(priorHyperParameters[requiredNames])

    if (!is.numeric(priorValues) || any(!is.finite(priorValues)) ||
        any(priorValues <= 0))
      stop("priorHyperParameters must be finite and positive.")

    result[["priorHyperParameters"]] <- priorHyperParameters[requiredNames]
  }

  result[["esMin"]] <- propDiffMin
  result[["parameter"]] <- c(
    "Beta hyperparameters" =
      paste(unlist(result[["priorHyperParameters"]]), collapse = " "))
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

# Find the means that minimize the KL between the alternative and null
# The alternative is usually learnt and given
# The null is H0: thetaB - thetaA = propDiff
solveRIPr2x2PropDiff <- function(thetaA, thetaB, na, nb, propDiff) {
  derivativeKL <- function(nullThetaA) {
    nullThetaB <- nullThetaA + propDiff
    na * ((1 - thetaA) / (1 - nullThetaA) - thetaA / nullThetaA) +
      nb * ((1 - thetaB) / (1 - nullThetaB) - thetaB / nullThetaB)
  }

  # The derivative is infinite at the edges, so search just inside them.
  lower <- max(0, -propDiff)
  upper <- min(1, 1 - propDiff)
  edge <- 1e-12 * (upper - lower)

  stats::uniroot(derivativeKL, lower = lower + edge, upper = upper - edge,
                 tol = 1e-12)[["root"]]
}

#' Anytime-valid confidence sequence for the proportion difference
#'
#' Inverts the test: each candidate `propDiff = thetaB - thetaA` on a grid is
#' a point null with its own e-process, whose numerator is the same
#' predictable Beta posterior mean as in `savi2x2Test()` and whose
#' denominator is that prediction projected onto the candidate's null line.
#' A candidate leaves the set for good once its e-process reaches `1/alpha`
#' (running intersection), so the sets are nested over blocks.
#'
#' @param ya,yb integer vectors, the successes in group A and group B in each
#'   block.
#' @param na,nb positive integers, the block sizes.
#' @param priorHyperParameters list with `betaA1`, `betaA2`, `betaB1`,
#'   `betaB2`.
#' @param alpha numeric in (0, 1); the sequence has coverage `1 - alpha`.
#' @param precision positive integer, the number of equally spaced
#'   candidates strictly inside `(-1, 1)`.
#'
#' @return A matrix with columns `block`, `lowerBound` and `upperBound`.
#'   Each run of consecutive non-rejected candidates after a block is one
#'   interval; the confidence set after that block is the union of its rows.
#'   Without holes there is one row per block, as for the z-test; a block
#'   whose candidates are all rejected has no row.
#' @noRd
computeConfidenceInterval2x2PropDiff <- function(ya, yb, na, nb,
                                                 priorHyperParameters,
                                                 alpha, precision = 100) {
  nBlocks <- length(ya)
  thetas <- predictiveThetas2x2(ya, yb, na, nb, priorHyperParameters)

  propDiffGrid <- seq(-1, 1, length.out = precision + 2)[-c(1, precision + 2)]
  logEValues <- numeric(precision)
  inSet <- rep(TRUE, precision)
  logThreshold <- log(1 / alpha)

  confSeqMatrix <- matrix(numeric(0), ncol = 3,
                          dimnames = list(NULL, c("block", "lowerBound",
                                                  "upperBound")))

  for (i in seq_len(nBlocks)) {
    # A rejected candidate never returns, so its e-process is not advanced.
    for (j in which(inSet)) {
      propDiff <- propDiffGrid[j]
      nullThetaA <- solveRIPr2x2PropDiff(
        thetaA = thetas[["thetaA"]][i], thetaB = thetas[["thetaB"]][i],
        na = na, nb = nb, propDiff = propDiff)

      logEValues[j] <- logEValues[j] + savi2x2TestStat(
        ya = ya[i], yb = yb[i], na = na, nb = nb,
        numeratorThetaA = thetas[["thetaA"]][i],
        numeratorThetaB = thetas[["thetaB"]][i],
        denominatorThetaA = nullThetaA,
        denominatorThetaB = nullThetaA + propDiff,
        log = TRUE)
    }

    inSet <- inSet & logEValues < logThreshold

    # Split the non-rejected candidates into runs of neighbours on the grid;
    # each run is one interval of the union.
    runs <- rle(inSet)
    runEnds <- cumsum(runs[["lengths"]])[runs[["values"]]]
    runStarts <- runEnds - runs[["lengths"]][runs[["values"]]] + 1

    confSeqMatrix <- rbind(
      confSeqMatrix,
      cbind("block" = rep(i, length(runStarts)),
            "lowerBound" = propDiffGrid[runStarts],
            "upperBound" = propDiffGrid[runEnds]))
  }

  return(confSeqMatrix)
}
