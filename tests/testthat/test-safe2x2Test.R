testthat::test_that("Turner stopping-time simulation retains non-crossing paths", {
  result <- simulateTurnerStoppingTimes(
    na = 1,
    nb = 1,
    thetaA = 0,
    thetaB = 0,
    alpha = 0.05,
    nSimulations = 3,
    maxBlocks = 5
  )

  testthat::expect_identical(result$crossed, rep(FALSE, 3))
  testthat::expect_true(all(is.infinite(result$stoppingTimes)))
  testthat::expect_true(all(is.na(result$eValuesAtStopping)))
  testthat::expect_identical(result$blocksSimulated, rep(5, 3))
})

testthat::test_that("Turner stopping-time simulation records crossing e-values", {
  set.seed(99)
  result <- simulateTurnerStoppingTimes(
    na = 2,
    nb = 1,
    thetaA = 0.02,
    thetaB = 0.82,
    alpha = 0.5,
    nSimulations = 5,
    maxBlocks = 100
  )

  testthat::expect_true(all(result$crossed))
  testthat::expect_true(all(is.finite(result$stoppingTimes)))
  testthat::expect_true(all(result$eValuesAtStopping >= 1 / 0.5))
  testthat::expect_identical(result$blocksSimulated, rep(100, 5))
})

testthat::test_that("e-process grid produces the prefix confidence sequence", {
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 12,
    alpha = 0.05
  )
  set.seed(20260813)
  ya <- stats::rbinom(12, size = 1, prob = 0.2)
  yb <- stats::rbinom(12, size = 1, prob = 0.5)

  confidenceSequence <-
    computeConfidenceSequenceForDifferenceTwoProportions(
      ya = ya,
      yb = yb,
      gridSize = 21,
      saviDesign = design
    )
  completeGrid <- calculateEValuesForLinearDeltaGrid(
    ya = ya,
    yb = yb,
    na = design[["nPlan"]][["na"]],
    nb = design[["nPlan"]][["nb"]],
    priorParameters = design[["betaPriorParameterValues"]],
    gridSize = 21
  )
  prefixBounds <- t(vapply(seq_along(ya), function(block) {
    gridResult <- calculateEValuesForLinearDeltaGrid(
      ya = ya[seq_len(block)],
      yb = yb[seq_len(block)],
      na = design[["nPlan"]][["na"]],
      nb = design[["nPlan"]][["nb"]],
      priorParameters = design[["betaPriorParameterValues"]],
      gridSize = 21
    )
    retained <- colSums(
      gridResult[["logEProcesses"]] >= log(1 / design[["alpha"]])
    ) == 0L
    included <- gridResult[["delta"]][retained]
    if (length(included) == 0L) {
      c(NA_real_, NA_real_)
    } else {
      range(included)
    }
  }, numeric(2)))

  testthat::expect_length(completeGrid[["delta"]], 21)
  testthat::expect_identical(
    dim(completeGrid[["logEProcesses"]]),
    c(12L, 21L)
  )
  testthat::expect_identical(confidenceSequence[["block"]], seq_along(ya))
  testthat::expect_equal(
    as.matrix(confidenceSequence[c("lowerBound", "upperBound")]),
    prefixBounds,
    tolerance = 0,
    ignore_attr = TRUE
  )
})

testthat::test_that("linear-difference grid covers the feasible interior", {
  priorParameters <- list(
    betaA1 = 0.18,
    betaA2 = 0.18,
    betaB1 = 0.18,
    betaB2 = 0.18
  )
  gridResult <- calculateEValuesForLinearDeltaGrid(
    ya = c(0, 1, 0),
    yb = c(1, 1, 0),
    na = 1,
    nb = 1,
    priorParameters = priorParameters,
    gridSize = 7
  )

  testthat::expect_equal(
    gridResult[["delta"]],
    seq_len(7) * (2 / 8) - 1
  )
  testthat::expect_identical(dim(gridResult[["logEProcesses"]]), c(3L, 7L))
})
