# Testing fnts ----

# compute conditional e-variable
# 1. given a logOR, just plug in: GROW case or UMP case
# 2. given a prior weight, numerically integrate it (report posterior maybe later)
savi2x2CondStat <- function(ya, yb, na, nb, logOR = NULL, alternative=c("twoSided", "greater", "less")) {}


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
