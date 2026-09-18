testthat::test_that("Turner stopping-time simulation retains non-crossing paths", {
  result <- sampleStoppingTimesSaviTwoProportions(
    thetaA = 0,
    thetaB = 0,
    na = 1,
    nb = 1,
    alpha = 0.05,
    restriction = "propDiff",
    delta = 0,
    nSim = 3,
    maxBlocks = 5
  )

  testthat::expect_identical(dim(result$stoppingTimes), c(1L, 3L))
  testthat::expect_identical(dim(result$eValuesAtStopping), c(1L, 3L))
  testthat::expect_true(all(is.infinite(result$stoppingTimes)))
  testthat::expect_true(all(is.na(result$eValuesAtStopping)))
  testthat::expect_named(
    result,
    c("stoppingTimes", "eValuesAtStopping"),
    ignore.order = FALSE
  )
})

testthat::test_that("Turner stopping-time simulation records crossing e-values", {
  set.seed(99)
  result <- sampleStoppingTimesSaviTwoProportions(
    thetaA = 0.02,
    thetaB = 0.82,
    na = 2,
    nb = 1,
    alpha = 0.5,
    restriction = "propDiff",
    delta = 0.8,
    nSim = 5,
    maxBlocks = 100
  )

  testthat::expect_true(all(is.finite(result$stoppingTimes)))
  testthat::expect_true(all(result$eValuesAtStopping >= 1 / 0.5))
  testthat::expect_true(all(result$stoppingTimes < 50))
})

testthat::test_that("restricted simulations match turnerEProcess paths", {
  priorParameters <- list(
    betaA1 = 0.18, betaA2 = 0.18, betaB1 = 0.18, betaB2 = 0.18
  )
  settings <- list(
    list(restriction = "propDiff", delta = 0.6),
    list(restriction = "logOR", delta = log(16))
  )
  na <- 2
  nb <- 1

  for (setting in settings) {
    set.seed(4242)
    simulated <- sampleStoppingTimesSaviTwoProportions(
      thetaA = 0.2,
      thetaB = 0.8,
      na = na,
      nb = nb,
      alpha = 0.05,
      priorParameters = priorParameters,
      restriction = setting$restriction,
      delta = setting$delta,
      gridSize = 101,
      nSim = 1,
      maxBlocks = 200
    )
    stoppingTime <- simulated$stoppingTimes[1, 1]
    testthat::expect_true(is.finite(stoppingTime))

    set.seed(4242)
    ya <- yb <- numeric(stoppingTime)
    for (block in seq_len(stoppingTime)) {
      ya[block] <- stats::rbinom(1, na, 0.2)
      yb[block] <- stats::rbinom(1, nb, 0.8)
    }
    replayed <- turnerEProcess(
      ya = ya, yb = yb, na = na, nb = nb,
      priorParameters = priorParameters,
      restriction = setting$restriction,
      delta = setting$delta,
      gridSize = 101
    )

    testthat::expect_equal(
      replayed[stoppingTime],
      simulated$eValuesAtStopping[1, 1],
      tolerance = 1e-10
    )
    testthat::expect_true(all(replayed[-stoppingTime] < 1 / 0.05))
  }
})

testthat::test_that("shared log-likelihood-ratio increments reproduce the process", {
  ya <- c(0, 1, 1)
  yb <- c(1, 0, 1)
  numeratorThetaA <- c(0.4, 0.3, 0.5)
  numeratorThetaB <- c(0.6, 0.7, 0.5)
  denominatorThetaA <- denominatorThetaB <- c(0.5, 0.45, 0.5)

  increments <- logLikelihoodRatioIncrements(
    ya, yb, na = 1, nb = 1,
    numeratorThetaA, numeratorThetaB,
    denominatorThetaA, denominatorThetaB
  )
  process <- logLikelihoodRatioProcess(
    ya, yb, na = 1, nb = 1,
    numeratorThetaA, numeratorThetaB,
    denominatorThetaA, denominatorThetaB
  )

  testthat::expect_equal(process, cumsum(increments), tolerance = 0)
})

testthat::test_that("coefficient-free log ratios match binomial log ratios", {
  ya <- c(0, 1, 2, 3)
  yb <- c(2, 1, 0, 2)
  na <- c(3, 3, 3, 3)
  nb <- c(2, 2, 2, 2)
  numeratorThetaA <- c(0.2, 0.4, 0.6, 0.8)
  numeratorThetaB <- c(0.7, 0.5, 0.3, 0.6)
  denominatorThetaA <- c(0.3, 0.5, 0.7, 0.9)
  denominatorThetaB <- c(0.6, 0.4, 0.2, 0.5)

  expected <-
    stats::dbinom(ya, na, numeratorThetaA, log = TRUE) +
    stats::dbinom(yb, nb, numeratorThetaB, log = TRUE) -
    stats::dbinom(ya, na, denominatorThetaA, log = TRUE) -
    stats::dbinom(yb, nb, denominatorThetaB, log = TRUE)

  testthat::expect_equal(
    logLikelihoodRatioIncrements(
      ya, yb, na, nb,
      numeratorThetaA, numeratorThetaB,
      denominatorThetaA, denominatorThetaB
    ),
    expected,
    tolerance = 1e-14
  )
})

testthat::test_that("coefficient-free log ratios handle shared boundaries", {
  increments <- logLikelihoodRatioIncrements(
    ya = c(1, 0),
    yb = c(0, 1),
    na = 1,
    nb = 1,
    numeratorThetaA = 0,
    numeratorThetaB = 0,
    denominatorThetaA = 0,
    denominatorThetaB = 0
  )

  testthat::expect_equal(increments, c(0, 0))
})

testthat::test_that("conditional Gaussian mixture updates only after each block", {
  ya <- c(2, 0, 2)
  yb <- c(0, 1, 0)
  na <- c(2, 1, 2)
  nb <- c(2, 1, 2)

  result <- computeEGaussGrid(
    ya, yb, na, nb,
    log = TRUE,
    returnPosteriors = TRUE
  )
  prefix <- computeEGaussGrid(ya[1:2], yb[1:2], na[1:2], nb[1:2], log = TRUE)

  logORGrid <- result$logORGrid
  initialLogWeights <- stats::dnorm(logORGrid, log = TRUE)
  initialLogWeights <- initialLogWeights - logSumExp(initialLogWeights)
  firstLogLikelihoodRatio <- vapply(
    logORGrid,
    function(logOR) {
      logOR * ya[1] - fnchPsi(na[1], nb[1], ya[1] + yb[1], logOR) +
        lchoose(na[1] + nb[1], ya[1] + yb[1])
    },
    numeric(1)
  )
  continuousFirst <- stats::integrate(
    function(logOR) {
      logLikelihoodRatio <- vapply(
        logOR,
        function(value) {
          value * ya[1] -
            fnchPsi(na[1], nb[1], ya[1] + yb[1], value) +
            lchoose(na[1] + nb[1], ya[1] + yb[1])
        },
        numeric(1)
      )
      exp(logLikelihoodRatio + stats::dnorm(logOR, log = TRUE))
    },
    lower = -Inf,
    upper = Inf,
    rel.tol = 1e-12,
    abs.tol = 1e-12
  )$value

  testthat::expect_equal(
    result$eValues[1],
    logSumExp(initialLogWeights + firstLogLikelihoodRatio),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    exp(result$eValues[1]),
    continuousFirst,
    tolerance = 1e-10
  )
  testthat::expect_equal(result$eValues[1:2], prefix, tolerance = 0)
  testthat::expect_equal(
    rowSums(result$posteriors),
    rep(1, length(ya) + 1L),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    exp(result$eValues),
    computeEGaussGrid(ya, yb, na, nb),
    tolerance = 1e-12
  )
})

testthat::test_that("Turner e-process rejects invalid binomial block data", {
  testthat::expect_error(
    turnerEProcess(ya = -1, yb = 0, na = 1, nb = 1),
    "Success counts"
  )
  testthat::expect_error(
    turnerEProcess(ya = 0.5, yb = 0, na = 1, nb = 1),
    "integer values"
  )
  testthat::expect_error(
    turnerEProcess(ya = 2, yb = 0, na = 1, nb = 1),
    "Success counts"
  )
  testthat::expect_error(
    turnerEProcess(ya = 0, yb = 0, na = 0, nb = 1),
    "positive values"
  )
})

testthat::test_that("Turner predictors use past data and B-minus-A restrictions", {
  priorParameters <- list(
    betaA1 = 1,
    betaA2 = 1,
    betaB1 = 2,
    betaB2 = 2
  )
  ya <- c(1, 0)
  yb <- c(0, 1)
  unrestricted <- learnPredictiveThetas(
    ya, yb, na = c(1, 1), nb = c(1, 1),
    priorParameters = priorParameters
  )
  propDiff <- learnPredictiveThetas(
    ya, yb, na = c(1, 1), nb = c(1, 1),
    priorParameters = priorParameters,
    restriction = "propDiff", delta = 0.25, gridSize = 101
  )
  logOdds <- learnPredictiveThetas(
    ya, yb, na = c(1, 1), nb = c(1, 1),
    priorParameters = priorParameters,
    restriction = "logOR", delta = 0.7, gridSize = 101
  )

  testthat::expect_equal(unrestricted$thetaA, c(0.5, 2 / 3))
  testthat::expect_equal(unrestricted$thetaB, c(0.5, 0.4))
  testthat::expect_equal(propDiff$thetaB - propDiff$thetaA, rep(0.25, 2))
  testthat::expect_equal(
    stats::qlogis(logOdds$thetaB) - stats::qlogis(logOdds$thetaA),
    rep(0.7, 2),
    tolerance = 1e-12
  )
})

testthat::test_that("restricted Turner processes use likelihood-updated predictors", {
  priorParameters <- list(
    betaA1 = 0.18,
    betaA2 = 0.18,
    betaB1 = 0.18,
    betaB2 = 0.18
  )
  ya <- c(0, 1, 0, 1)
  yb <- c(1, 1, 0, 0)

  propDiffProcess <- turnerEProcess(
    ya, yb, na = 1, nb = 1,
    priorParameters = priorParameters,
    restriction = "propDiff", delta = 0.25, gridSize = 101
  )
  logORProcess <- turnerEProcess(
    ya, yb, na = 1, nb = 1,
    priorParameters = priorParameters,
    restriction = "logOR", delta = 0.7, gridSize = 101
  )

  testthat::expect_equal(
    propDiffProcess,
    c(1.5625, 1.46484375, 1.1474712811576, 0.645452595651148),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    logORProcess,
    c(1.37527821341647, 1.33960921187147, 1.25210874368758, 0.855839806438304),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    log(propDiffProcess),
    turnerEProcess(
      ya, yb, na = 1, nb = 1,
      priorParameters = priorParameters,
      restriction = "propDiff", delta = 0.25, gridSize = 101,
      log = TRUE
    ),
    tolerance = 1e-12
  )
})

testthat::test_that("legacy public effect names map to propDiff and logOR", {
  priorParameters <- list(
    betaA1 = 0.18,
    betaA2 = 0.18,
    betaB1 = 0.18,
    betaB2 = 0.18
  )
  process <- function(restriction) {
    turnerEProcess(
      ya = c(0, 1),
      yb = c(1, 0),
      na = 1,
      nb = 1,
      priorParameters = priorParameters,
      restriction = restriction,
      delta = 0.25,
      gridSize = 101
    )
  }

  # Legacy spellings are confined to the design entry point: everything else
  # speaks propDiff and logOR only, and rejects the old names outright.
  testthat::expect_error(process("difference"), "should be one of")
  testthat::expect_error(process("linearDifference"), "should be one of")
  testthat::expect_error(process("logOddsRatio"), "should be one of")
  testthat::expect_error(
    learnPredictiveThetas(
      ya = 0, yb = 1, na = 1, nb = 1,
      priorParameters = priorParameters,
      restriction = "logOddsRatio", delta = 0.5
    ),
    "should be one of"
  )
  testthat::expect_equal(
    suppressWarnings(designSaviTwoProportions(
      na = 1, nb = 1, nBlocksPlan = 1,
      effectMeasure = "linearDifference"
    ))$effectMeasure,
    "propDiff"
  )
  testthat::expect_equal(
    suppressWarnings(designSaviTwoProportions(
      na = 1, nb = 1, nBlocksPlan = 1,
      effectMeasure = "logOddsRatio"
    ))$effectMeasure,
    "logOR"
  )
  testthat::expect_equal(
    suppressWarnings(designSaviTwoProportions(
      na = 1, nb = 1, nBlocksPlan = 1,
      effectMeasure = "difference"
    ))$effectMeasure,
    "propDiff"
  )
})

testthat::test_that("batched Turner simulation is reproducible across a partial batch", {
  set.seed(20260817)
  first <- sampleStoppingTimesSaviTwoProportions(
    thetaA = 0.1,
    thetaB = 0.9,
    na = 1,
    nb = 1,
    alpha = 0.05,
    restriction = "propDiff",
    delta = 0.8,
    nSim = 6,
    maxBlocks = 51
  )
  set.seed(20260817)
  second <- sampleStoppingTimesSaviTwoProportions(
    thetaA = 0.1,
    thetaB = 0.9,
    na = 1,
    nb = 1,
    alpha = 0.05,
    restriction = "propDiff",
    delta = 0.8,
    nSim = 6,
    maxBlocks = 51
  )

  testthat::expect_identical(first, second)
  testthat::expect_true(all(first$stoppingTimes < 50))

  tailBatch <- sampleStoppingTimesSaviTwoProportions(
    thetaA = 0,
    thetaB = 0,
    na = 1,
    nb = 1,
    alpha = 0.05,
    restriction = "propDiff",
    delta = 0,
    nSim = 2,
    maxBlocks = 51
  )
  testthat::expect_true(all(is.infinite(tailBatch$stoppingTimes)))
})

testthat::test_that("fixed propDiff simulation returns baseline-grid matrices", {
  set.seed(20260817)
  thetaGrid <- makeSimulationThetaGrid(
    propDiff = 0.6,
    thetaGridSize = 3
  )
  result <- sampleStoppingTimesSaviTwoProportions(
    thetaA = thetaGrid$thetaA,
    thetaB = thetaGrid$thetaB,
    na = 1,
    nb = 1,
    alpha = 0.05,
    restriction = thetaGrid$restriction,
    delta = thetaGrid$delta,
    nSim = 4,
    maxBlocks = 51
  )

  testthat::expect_equal(thetaGrid$thetaA, c(0.1, 0.2, 0.3))
  testthat::expect_equal(thetaGrid$thetaB, c(0.7, 0.8, 0.9))
  testthat::expect_identical(dim(result$stoppingTimes), c(3L, 4L))
  testthat::expect_identical(dim(result$eValuesAtStopping), c(3L, 4L))
  testthat::expect_identical(
    is.na(result$eValuesAtStopping),
    is.infinite(result$stoppingTimes)
  )
})

testthat::test_that("simulation theta grids use the B-minus-A convention", {
  negativeGrid <- makeSimulationThetaGrid(
    propDiff = -0.6,
    thetaGridSize = 3
  )
  logOddsGrid <- makeSimulationThetaGrid(
    logOR = 0.8,
    thetaGridSize = 3
  )
  negativeLogOddsGrid <- makeSimulationThetaGrid(
    logOR = -0.8,
    thetaGridSize = 3
  )

  testthat::expect_equal(negativeGrid$thetaA, c(0.7, 0.8, 0.9))
  testthat::expect_equal(negativeGrid$thetaB, c(0.1, 0.2, 0.3))
  testthat::expect_equal(
    stats::qlogis(logOddsGrid$thetaB) - stats::qlogis(logOddsGrid$thetaA),
    rep(0.8, 3),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    stats::qlogis(negativeLogOddsGrid$thetaB) -
      stats::qlogis(negativeLogOddsGrid$thetaA),
    rep(-0.8, 3),
    tolerance = 1e-12
  )
})

testthat::test_that("Turner stopping grids preserve dimensions when N is one", {
  result <- sampleStoppingTimesSaviTwoProportions(
    thetaA = c(0, 0.1),
    thetaB = c(0, 0.9),
    na = 1,
    nb = 1,
    restriction = "propDiff",
    delta = 0.8,
    nSim = 1,
    maxBlocks = 3
  )

  testthat::expect_identical(dim(result$stoppingTimes), c(2L, 1L))
  testthat::expect_identical(dim(result$eValuesAtStopping), c(2L, 1L))
  testthat::expect_error(
    sampleStoppingTimesSaviTwoProportions(
      thetaA = c(0.1, 0.2), thetaB = 0.3,
      na = 1, nb = 1, restriction = "propDiff", delta = 0.1,
      nSim = 1, maxBlocks = 1
    ),
    "equal-length"
  )
})

testthat::test_that("Turner stopping grids require a restricted numerator", {
  testthat::expect_error(
    sampleStoppingTimesSaviTwoProportions(
      thetaA = 0.2,
      thetaB = 0.8,
      na = 1,
      nb = 1,
      restriction = "none",
      nSim = 1,
      maxBlocks = 1
    ),
    "should be one of"
  )
})

testthat::test_that("worst-case wrapper simulates and bootstraps the grid", {
  set.seed(20260817)
  result <- computeNPlanSaviTwoProportions(
    na = 1,
    nb = 1,
    propDiff = 0.6,
    alpha = 0.05,
    beta = 0.2,
    nSim = 4,
    maxBlocks = 51,
    thetaGridSize = 3,
    nBoot = 2
  )

  testthat::expect_identical(dim(result$stoppingTimes), c(3L, 4L))
  testthat::expect_identical(dim(result$eValuesAtStopping), c(3L, 4L))
  testthat::expect_length(result$stoppingTimeQuantiles, 3)
  testthat::expect_true(result$worstCaseIndex %in% seq_len(3))
  testthat::expect_equal(
    result$worstCaseThetaB - result$worstCaseThetaA,
    0.6
  )
})

testthat::test_that("two-proportion planning uses the worst-case simulator", {
  set.seed(7)
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    delta = 0.8,
    beta = 0.2,
    alpha = 0.5,
    nSim = 2
  )

  testthat::expect_true(is.finite(design$nPlan[["nBlocksPlan"]]))
  testthat::expect_length(design$nPlanTwoSe, 3)
  testthat::expect_equal(design$esMin[["propDiff"]], 0.8)
})

testthat::test_that("two-proportion design estimates beta at a fixed horizon", {
  set.seed(8)
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 20,
    delta = 0.8,
    alpha = 0.5,
    nSim = 10
  )

  testthat::expect_gte(design$beta, 0)
  testthat::expect_lte(design$beta, 1)
  testthat::expect_gte(design$betaTwoSe, 0)
  testthat::expect_null(design$logImpliedTarget)
  testthat::expect_true(any(grepl(
    "Implied target is not calculated",
    design$note,
    fixed = TRUE
  )))
})

testthat::test_that("two-proportion design searches for a detectable propDiff", {
  set.seed(2)
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 20,
    beta = 0.8,
    alpha = 0.5,
    nSim = 10
  )

  testthat::expect_gte(design$esMin[["propDiff"]], 0.01)
  testthat::expect_lte(design$esMin[["propDiff"]], 0.99)
  testthat::expect_equal(design$beta, 0.8)
  testthat::expect_true(any(grepl(
    "binary search over a simulation grid",
    design$note,
    fixed = TRUE
  )))
})

testthat::test_that("minimum detectable effect uses boundary-checked binary search", {
  calls <- numeric()
  powerRule <- function(effect) effect
  testthat::local_mocked_bindings(
    computePowerSaviTwoProportions = function(
      ...,
      propDiff = NULL,
      logOR = NULL
    ) {
      effect <- if (is.null(propDiff)) logOR else propDiff
      calls <<- c(calls, effect)
      list(worstCasePower = powerRule(effect))
    },
    .package = "safestats"
  )

  selected <- computeMinEsSaviTwoProportions(
    na = 1,
    nb = 1,
    alpha = 0.05,
    beta = 0.5,
    maxBlocks = 10,
    effectGridSize = 10,
    effectMax = 0.9,
    effectMin = 0.1
  )
  candidateGrid <- seq(0.9, 0.1, length.out = 10)
  testthat::expect_equal(selected, min(candidateGrid[candidateGrid >= 0.5]))
  testthat::expect_lte(length(calls), 6L)

  powerRule <- function(effect) 0
  testthat::expect_error(
    computeMinEsSaviTwoProportions(
      na = 1, nb = 1, alpha = 0.05, beta = 0.5,
      maxBlocks = 10, effectGridSize = 4
    ),
    "No candidate effect"
  )

  powerRule <- function(effect) 1
  testthat::expect_equal(
    computeMinEsSaviTwoProportions(
      na = 1, nb = 1, alpha = 0.05, beta = 0.5,
      maxBlocks = 10, effectGridSize = 4,
      effectMax = 0.9, effectMin = 0.1
    ),
    0.1
  )
  testthat::expect_equal(
    computeMinEsSaviTwoProportions(
      na = 1, nb = 1, alpha = 0.05, beta = 0.5,
      maxBlocks = 10, effectGridSize = 1,
      effectMax = 0.8, effectMin = 0.1
    ),
    0.8
  )
})

testthat::test_that("propDiff confidence sequence is a prefix property", {
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
    computeConfidenceSequenceForPropDiffTwoProportions(
      ya = ya,
      yb = yb,
      confidenceBoundGridPrecision = 21,
      saviDesign = design
    )
  # The bounds at block t may depend on the data up to t only, so running
  # the sequence on every prefix must reproduce the full sequence row by row.
  prefixBounds <- t(vapply(seq_along(ya), function(block) {
    prefixSequence <- computeConfidenceSequenceForPropDiffTwoProportions(
      ya = ya[seq_len(block)],
      yb = yb[seq_len(block)],
      confidenceBoundGridPrecision = 21,
      saviDesign = design
    )
    unlist(prefixSequence[block, c("lowerBound", "upperBound")])
  }, numeric(2)))

  testthat::expect_identical(confidenceSequence[["block"]], seq_along(ya))
  testthat::expect_equal(
    as.matrix(confidenceSequence[c("lowerBound", "upperBound")]),
    prefixBounds,
    tolerance = 0,
    ignore_attr = TRUE
  )
  testthat::expect_true(all(
    confidenceSequence[["lowerBound"]] <= confidenceSequence[["upperBound"]],
    na.rm = TRUE
  ))
})

testthat::test_that("propDiff grid covers the feasible interior", {
  testthat::expect_equal(
    propDiffCandidateGrid(7),
    c(-rev(seq_len(7) / 8), seq_len(7) / 8)
  )
  testthat::expect_length(propDiffCandidateGrid(21), 42)
  testthat::expect_error(propDiffCandidateGrid(0), "gridSize")
  testthat::expect_error(propDiffCandidateGrid(2.5), "gridSize")
  testthat::expect_error(propDiffCandidateGrid(c(3, 4)), "gridSize")
})

testthat::test_that("propDiff candidates exclude zero and stay symmetric", {
  for (gridSize in c(1, 2, 7, 8, 100, 101)) {
    candidates <- propDiffCandidateGrid(gridSize)

    testthat::expect_false(any(candidates == 0))
    testthat::expect_gt(min(abs(candidates)), 0)
    testthat::expect_equal(candidates, -rev(candidates))
    testthat::expect_false(is.unsorted(candidates, strictly = TRUE))
    testthat::expect_true(all(abs(candidates) < 1))
  }
})

testthat::test_that("propDiff RIPr is feasible and stationary", {
  numeratorThetaA <- c(0.2, 0.7, 0.4)
  numeratorThetaB <- c(0.7, 0.2, 0.6)
  na <- c(1, 3, 2)
  nb <- c(2, 1, 5)
  propDiff <- 0.2
  riprThetaA <- mapply(
    solveOnePropDiffRIPr,
    thetaStarA = numeratorThetaA,
    thetaStarB = numeratorThetaB,
    blockSizeA = na,
    blockSizeB = nb,
    MoreArgs = list(propDiff = propDiff)
  )
  ripr <- list(thetaA = riprThetaA, thetaB = riprThetaA + propDiff)
  score <-
    na * (ripr$thetaA - numeratorThetaA) /
      (ripr$thetaA * (1 - ripr$thetaA)) +
    nb * (ripr$thetaB - numeratorThetaB) /
      (ripr$thetaB * (1 - ripr$thetaB))

  testthat::expect_equal(ripr$thetaB - ripr$thetaA, rep(propDiff, 3))
  testthat::expect_true(all(ripr$thetaA > 0 & ripr$thetaA < 1))
  testthat::expect_true(all(ripr$thetaB > 0 & ripr$thetaB < 1))
  testthat::expect_equal(score, rep(0, 3), tolerance = 1e-10)
})

testthat::test_that("logOR inversion keeps sentinels and permanent rejection", {
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 4,
    alpha = 0.1
  )
  threshold <- log(2 / design[["alpha"]])
  # Four candidates c(-g2, -g1, g1, g2). Column j is candidate j; a family's
  # e-process reaching the threshold in a block rejects that candidate there.
  # Candidates 1 and 4 fall back below the threshold afterwards but must stay
  # rejected, and in the last block each family has rejected everything.
  lowerLogE <- rbind(
    c(0, 0, 0, 0),
    c(threshold, 0, 0, 0),
    c(0, threshold, 0, 0),
    c(0, 0, threshold, threshold)
  )
  upperLogE <- rbind(
    c(0, 0, 0, 0),
    c(0, 0, 0, threshold),
    c(0, 0, threshold, 0),
    c(threshold, threshold, 0, 0)
  )
  seenGrid <- NULL
  testthat::local_mocked_bindings(
    calculateEValuesForLogORGrid = function(logORGrid, bound, ...) {
      seenGrid <<- logORGrid
      list(
        logOR = logORGrid,
        logEProcesses = if (bound == "lower") lowerLogE else upperLogE
      )
    }
  )

  bounds <- computeConfidenceSequenceForLogORTwoProportions(
    ya = rep(0, 4),
    yb = rep(1, 4),
    confidenceBoundGridPrecision = 2,
    logORConfidenceSearchBounds = c(0.5, 2),
    saviDesign = design
  )

  testthat::expect_length(seenGrid, 4L)
  testthat::expect_equal(
    as.matrix(bounds[c("lowerBound", "upperBound")]),
    rbind(
      c(-Inf, Inf),
      c(seenGrid[2L], seenGrid[3L]),
      c(seenGrid[3L], seenGrid[2L]),
      c(NA_real_, NA_real_)
    ),
    tolerance = 0,
    ignore_attr = TRUE
  )
})

testthat::test_that("propDiff running intersection can be switched off", {
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 12,
    alpha = 0.05
  )
  set.seed(20260917)
  ya <- stats::rbinom(12, size = 1, prob = 0.2)
  yb <- stats::rbinom(12, size = 1, prob = 0.8)
  # Reference: every candidate's complete log e-process, built from the
  # same predictor, solver, and increments the sequence uses block by block.
  predictiveThetas <- learnPredictiveThetas(
    ya = ya,
    yb = yb,
    na = rep(1, length(ya)),
    nb = rep(1, length(ya)),
    priorParameters = design[["betaPriorParameterValues"]],
    restriction = "none"
  )
  grid <- propDiffCandidateGrid(10)
  logEProcesses <- vapply(grid, function(propDiff) {
    nullThetaA <- mapply(
      solveOnePropDiffRIPr,
      thetaStarA = predictiveThetas[["thetaA"]],
      thetaStarB = predictiveThetas[["thetaB"]],
      blockSizeA = 1,
      blockSizeB = 1,
      MoreArgs = list(propDiff = propDiff)
    )
    logLikelihoodRatioProcess(
      ya = ya,
      yb = yb,
      na = 1,
      nb = 1,
      numeratorThetaA = predictiveThetas[["thetaA"]],
      numeratorThetaB = predictiveThetas[["thetaB"]],
      denominatorThetaA = nullThetaA,
      denominatorThetaB = nullThetaA + propDiff
    )
  }, numeric(length(ya)))
  logThreshold <- log(1 / design[["alpha"]])
  boundsFrom <- function(inSetMatrix) {
    t(apply(inSetMatrix, 1L, function(inSet) {
      if (!any(inSet)) return(c(NA_real_, NA_real_))
      c(
        if (inSet[1L]) -1 else min(grid[inSet]),
        if (inSet[length(grid)]) 1 else max(grid[inSet])
      )
    }))
  }

  withIntersection <- computeConfidenceSequenceForPropDiffTwoProportions(
    ya = ya, yb = yb, confidenceBoundGridPrecision = 10, saviDesign = design
  )
  withoutIntersection <- computeConfidenceSequenceForPropDiffTwoProportions(
    ya = ya, yb = yb, confidenceBoundGridPrecision = 10, saviDesign = design,
    runningIntersection = FALSE
  )

  testthat::expect_equal(
    withIntersection,
    computeConfidenceSequenceForPropDiffTwoProportions(
      ya = ya, yb = yb, confidenceBoundGridPrecision = 10,
      saviDesign = design, runningIntersection = TRUE
    )
  )
  testthat::expect_equal(
    as.matrix(withoutIntersection[c("lowerBound", "upperBound")]),
    boundsFrom(logEProcesses < logThreshold),
    ignore_attr = TRUE
  )
  testthat::expect_equal(
    as.matrix(withIntersection[c("lowerBound", "upperBound")]),
    boundsFrom(
      apply(logEProcesses, 2L, cummax) < logThreshold
    ),
    ignore_attr = TRUE
  )
  testthat::expect_true(all(
    withIntersection[["lowerBound"]] >= withoutIntersection[["lowerBound"]],
    na.rm = TRUE
  ))
  testthat::expect_true(all(
    withIntersection[["upperBound"]] <= withoutIntersection[["upperBound"]],
    na.rm = TRUE
  ))
  testthat::expect_error(
    computeConfidenceSequenceForPropDiffTwoProportions(
      ya = ya, yb = yb, confidenceBoundGridPrecision = 10,
      saviDesign = design, runningIntersection = NA
    ),
    "runningIntersection"
  )
})

testthat::test_that("propDiff confidence sequence starts at the parameter limits", {
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 1,
    alpha = 0.05
  )
  confidenceSequence <- computeConfidenceSequenceForPropDiffTwoProportions(
    ya = 0,
    yb = 0,
    confidenceBoundGridPrecision = 3L,
    saviDesign = design
  )

  testthat::expect_equal(
    unname(unlist(confidenceSequence[1, c("lowerBound", "upperBound")])),
    c(-1, 1)
  )
})

testthat::test_that("log-odds-ratio RIPr satisfies the score equation", {
  thetaA <- c(0.2, 0.6, 0.8)
  thetaB <- c(0.7, 0.3, 0.4)
  na <- c(1, 3, 2)
  nb <- c(2, 1, 4)
  logOR <- 0.8
  ripr <- solveLogORRIPr(
    numeratorThetaA = thetaA,
    numeratorThetaB = thetaB,
    na = na,
    nb = nb,
    logOR = logOR
  )

  testthat::expect_equal(
    na * ripr$thetaA + nb * ripr$thetaB,
    na * thetaA + nb * thetaB,
    tolerance = 1e-10
  )
  testthat::expect_equal(
    stats::qlogis(ripr$thetaB) - stats::qlogis(ripr$thetaA),
    rep(logOR, 3),
    tolerance = 1e-10
  )
})

testthat::test_that("analytic log-odds-ratio RIPr matches numerical minimization", {
  cases <- expand.grid(
    thetaA = c(0.01, 0.2, 0.8, 0.99),
    thetaB = c(0.02, 0.4, 0.9),
    na = c(1, 7),
    nb = c(2, 5),
    logOR = c(-5, -0.5, 0, 0.5, 5)
  )

  for (case in seq_len(nrow(cases))) {
    values <- cases[case, ]
    ripr <- solveLogORRIPr(
      numeratorThetaA = values$thetaA,
      numeratorThetaB = values$thetaB,
      na = values$na,
      nb = values$nb,
      logOR = values$logOR
    )
    targetSuccesses <-
      values$na * values$thetaA + values$nb * values$thetaB
    numericalLogitA <- stats::uniroot(
      function(logitA) {
        values$na * stats::plogis(logitA) +
          values$nb * stats::plogis(logitA + values$logOR) -
          targetSuccesses
      },
      interval = c(-50, 50),
      tol = .Machine$double.eps^0.75
    )$root

    testthat::expect_equal(
      unname(c(ripr$thetaA, ripr$thetaB)),
      c(
        stats::plogis(numericalLogitA),
        stats::plogis(numericalLogitA + values$logOR)
      ),
      tolerance = 1e-10
    )
  }
})

testthat::test_that("log-odds-ratio RIPr handles pooled and boundary cases", {
  zeroEffect <- solveLogORRIPr(
    numeratorThetaA = 0.2,
    numeratorThetaB = 0.8,
    na = 1,
    nb = 3,
    logOR = 0
  )
  boundaries <- solveLogORRIPr(
    numeratorThetaA = c(0, 1),
    numeratorThetaB = c(0, 1),
    na = 1,
    nb = 2,
    logOR = 0.8
  )

  testthat::expect_equal(zeroEffect$thetaA, 0.65)
  testthat::expect_equal(zeroEffect$thetaB, 0.65)
  testthat::expect_equal(boundaries$thetaA, c(0, 1))
  testthat::expect_equal(boundaries$thetaB, c(0, 1))
})

testthat::test_that("log-odds-ratio grid produces the prefix confidence sequence", {
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 12,
    alpha = 0.05
  )
  set.seed(20260814)
  ya <- stats::rbinom(12, size = 1, prob = 0.2)
  yb <- stats::rbinom(12, size = 1, prob = 0.5)
  confidenceBoundGridPrecision <- 10L
  logORConfidenceSearchBounds <- c(0.01, 3)

  confidenceSequence <-
    computeConfidenceSequenceForLogORTwoProportions(
      ya = ya,
      yb = yb,
      confidenceBoundGridPrecision = confidenceBoundGridPrecision,
      logORConfidenceSearchBounds = logORConfidenceSearchBounds,
      saviDesign = design
  )
  prefixBounds <- t(vapply(seq_along(ya), function(block) {
    prefixSequence <- computeConfidenceSequenceForLogORTwoProportions(
      ya = ya[seq_len(block)],
      yb = yb[seq_len(block)],
      confidenceBoundGridPrecision = confidenceBoundGridPrecision,
      logORConfidenceSearchBounds = logORConfidenceSearchBounds,
      saviDesign = design
    )
    unname(unlist(prefixSequence[block, c("lowerBound", "upperBound")]))
  }, numeric(2)))

  testthat::expect_identical(
    confidenceSequence[["block"]],
    seq_along(ya)
  )
  testthat::expect_equal(
    as.matrix(confidenceSequence[c("lowerBound", "upperBound")]),
    prefixBounds,
    tolerance = 0,
    ignore_attr = TRUE
  )
})

testthat::test_that("log-odds-ratio grid respects bound direction and group swaps", {
  priorParameters <- list(
    betaA1 = 0.18,
    betaA2 = 0.18,
    betaB1 = 0.18,
    betaB2 = 0.18
  )
  ya <- c(0, 1, 0, 1)
  yb <- c(1, 0, 1, 1)
  positiveGrid <- c(0, 0.5, 1)
  lowerGrid <- calculateEValuesForLogORGrid(
    ya = ya,
    yb = yb,
    na = 1,
    nb = 1,
    priorParameters = priorParameters,
    logORGrid = positiveGrid,
    bound = "lower"
  )
  upperGrid <- calculateEValuesForLogORGrid(
    ya = yb,
    yb = ya,
    na = 1,
    nb = 1,
    priorParameters = priorParameters,
    logORGrid = -positiveGrid,
    bound = "upper"
  )

  testthat::expect_equal(
    unname(lowerGrid[["logEProcesses"]]),
    unname(upperGrid[["logEProcesses"]]),
    tolerance = 1e-10
  )
  signedGrid <- calculateEValuesForLogORGrid(
    ya = ya,
    yb = yb,
    na = 1,
    nb = 1,
    priorParameters = priorParameters,
    logORGrid = c(-1, 0, 1),
    bound = "lower"
  )
  testthat::expect_identical(
    dim(signedGrid[["logEProcesses"]]),
    c(4L, 3L)
  )
})

testthat::test_that("log-odds-ratio bounds include zero before any grid rejection", {
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 1,
    alpha = 0.05
  )

  confidenceSequence <- computeConfidenceSequenceForLogORTwoProportions(
    ya = 0,
    yb = 0,
    confidenceBoundGridPrecision = 2L,
    logORConfidenceSearchBounds = c(0.1, 0.5),
    saviDesign = design
  )

  testthat::expect_equal(
    unname(unlist(confidenceSequence[1, c("lowerBound", "upperBound")])),
    c(-Inf, Inf)
  )
})

testthat::test_that("logOR candidates exclude zero and stay symmetric", {
  searchBounds <- c(0.01, 3)

  for (precision in c(1L, 2L, 7L, 20L)) {
    # Rebuild the grid the wrapper uses, then confirm the wrapper agrees with
    # it by checking that no reported bound ever lands on zero.
    positiveTransformedBounds <- tanh(searchBounds / 4)
    positiveGrid <- 4 * atanh(seq(
      positiveTransformedBounds[1L], positiveTransformedBounds[2L],
      length.out = precision
    ))
    candidateGrid <- c(-rev(positiveGrid), positiveGrid)

    testthat::expect_length(candidateGrid, 2 * precision)
    testthat::expect_false(any(candidateGrid == 0))
    testthat::expect_equal(candidateGrid, -rev(candidateGrid))
    testthat::expect_false(is.unsorted(candidateGrid, strictly = TRUE))
  }

  set.seed(20260826)
  ya <- stats::rbinom(60, size = 1, prob = 0.15)
  yb <- stats::rbinom(60, size = 1, prob = 0.6)
  confidenceSequence <- computeConfidenceSequenceForLogORTwoProportions(
    ya = ya, yb = yb,
    confidenceBoundGridPrecision = 40L,
    logORConfidenceSearchBounds = searchBounds,
    saviDesign = suppressWarnings(designSaviTwoProportions(
      na = 1, nb = 1, nBlocksPlan = 60, alpha = 0.05
    ))
  )
  bounds <- unlist(confidenceSequence[c("lowerBound", "upperBound")])
  testthat::expect_false(any(bounds == 0, na.rm = TRUE))
})

testthat::test_that("positive log-odds data can produce finite bounds", {
  design <- designSaviTwoProportions(
    na = 1,
    nb = 1,
    nBlocksPlan = 95,
    alpha = 0.05
  )
  set.seed(19012022)
  ya <- stats::rbinom(95, size = 1, prob = 0.2)
  yb <- stats::rbinom(95, size = 1, prob = 0.5)

  confidenceSequence <- computeConfidenceSequenceForLogORTwoProportions(
    ya = ya,
    yb = yb,
    confidenceBoundGridPrecision = 100L,
    logORConfidenceSearchBounds = c(0.001, 5),
    saviDesign = design
  )

  testthat::expect_true(any(is.finite(confidenceSequence[["lowerBound"]])))
  testthat::expect_gt(confidenceSequence[["lowerBound"]][95], 0)
  testthat::expect_gt(confidenceSequence[["upperBound"]][95], 0)
})

testthat::test_that("saviTwoProportionsTest reproduces turnerEProcess", {
  design <- suppressWarnings(
    designSaviTwoProportions(na = 1, nb = 1, nBlocksPlan = 20)
  )
  set.seed(707)
  ya <- stats::rbinom(20, 1, 0.2)
  yb <- stats::rbinom(20, 1, 0.7)

  result <- saviTwoProportionsTest(ya, yb, designObj = design)

  testthat::expect_s3_class(result, "saviTest")
  testthat::expect_equal(
    result$eValueVec,
    unname(turnerEProcess(
      ya, yb, na = 1, nb = 1,
      priorParameters = design$betaPriorParameterValues
    ))
  )
  testthat::expect_equal(result$eValue, result$eValueVec[20])
  testthat::expect_equal(unname(result$n), c(20, 20, 20))
  testthat::expect_equal(
    unname(result$estimate),
    c(mean(ya), mean(yb), mean(yb) - mean(ya))
  )
})

testthat::test_that("saviTwoProportionsTest honours the design's restriction", {
  set.seed(708)
  design <- designSaviTwoProportions(
    na = 1, nb = 1, nBlocksPlan = 30, delta = 0.4, nSim = 20, nBoot = 20
  )
  ya <- stats::rbinom(30, 1, 0.2)
  yb <- stats::rbinom(30, 1, 0.7)

  testthat::expect_equal(design$alternativeRestriction, "propDiff")
  testthat::expect_equal(
    saviTwoProportionsTest(ya, yb, designObj = design)$eValueVec,
    unname(turnerEProcess(
      ya, yb, na = 1, nb = 1,
      priorParameters = design$betaPriorParameterValues,
      restriction = "propDiff", delta = 0.4
    ))
  )
})

testthat::test_that("the test result carries the Beta posterior forward", {
  design <- suppressWarnings(
    designSaviTwoProportions(na = 1, nb = 1, nBlocksPlan = 10)
  )
  ya <- c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1)
  yb <- c(0, 0, 1, 0, 1, 0, 0, 1, 0, 0)
  result <- saviTwoProportionsTest(ya, yb, designObj = design)
  prior <- design$betaPriorParameterValues

  testthat::expect_equal(
    result$posteriorHyperParameters$betaA1, prior$betaA1 + sum(ya)
  )
  testthat::expect_equal(
    result$posteriorHyperParameters$betaA2, prior$betaA2 + 10 - sum(ya)
  )

  # designSaviTwoProportions accepts this result directly, which was
  # previously unreachable: nothing produced posteriorHyperParameters.
  followUp <- suppressWarnings(designSaviTwoProportions(
    na = 1, nb = 1, nBlocksPlan = 5, previousSaviTestResult = result
  ))
  testthat::expect_equal(
    followUp$betaPriorParameterValues[["betaA1"]], prior$betaA1 + sum(ya)
  )
})

testthat::test_that("the confidence sequence matches the standalone wrapper", {
  design <- suppressWarnings(
    designSaviTwoProportions(na = 1, nb = 1, nBlocksPlan = 25)
  )
  set.seed(709)
  ya <- stats::rbinom(25, 1, 0.2)
  yb <- stats::rbinom(25, 1, 0.7)
  result <- saviTwoProportionsTest(
    ya, yb, designObj = design,
    wantConfidenceSequence = TRUE, confidenceBoundGridPrecision = 30
  )
  expected <- computeConfidenceSequenceForPropDiffTwoProportions(
    ya = ya, yb = yb, confidenceBoundGridPrecision = 30, saviDesign = design
  )

  testthat::expect_equal(
    result$confSeqMatrix,
    as.matrix(expected[, c("lowerBound", "upperBound")]),
    ignore_attr = TRUE
  )
  testthat::expect_equal(unname(result$confSeq), unname(result$confSeqMatrix[25, ]))
  testthat::expect_equal(result$ciValue, 1 - design$alpha)
})

testthat::test_that("the formula method splits by group and warns on imbalance", {
  design <- suppressWarnings(
    designSaviTwoProportions(na = 1, nb = 1, nBlocksPlan = 6)
  )
  frame <- data.frame(
    success = c(1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0),
    group = factor(rep(c("a", "b"), each = 6))
  )
  fromFormula <- saviTwoProportionsTest(
    success ~ group, data = frame, designObj = design
  )
  fromVectors <- saviTwoProportionsTest(
    ya = frame$success[1:6], yb = frame$success[7:12], designObj = design
  )

  testthat::expect_equal(fromFormula$eValueVec, fromVectors$eValueVec)
  testthat::expect_equal(fromFormula$dataName, "success by group")
  testthat::expect_match(names(fromFormula$estimate)[3L], "^propDiff \\(b - a\\)$")

  unbalanced <- frame[-12L, ]
  testthat::expect_warning(
    saviTwoProportionsTest(success ~ group, data = unbalanced, designObj = design),
    "unequal sizes"
  )
  testthat::expect_error(
    saviTwoProportionsTest(
      success ~ group,
      data = data.frame(success = c(1, 0, 1), group = factor(c("a", "b", "c"))),
      designObj = design
    ),
    "exactly two levels"
  )
})

testthat::test_that("a missing design falls back to a pilot with a warning", {
  testthat::expect_warning(
    result <- saviTwoProportionsTest(c(0, 1, 1), c(1, 1, 0)),
    "No designObj given"
  )
  testthat::expect_true(result$designObj$pilot)
  testthat::expect_length(result$eValueVec, 3)
})

testthat::test_that("restricted supports reject a degenerate delta", {
  priorParameters <- list(
    betaA1 = 0.18, betaA2 = 0.18, betaB1 = 0.18, betaB2 = 0.18
  )
  # plogis() saturates to exactly 1 near 37, which used to make 0 * -Inf = NaN
  # and spread through the posterior normalisation.
  testthat::expect_error(
    turnerEProcess(
      c(0, 1), c(1, 0), na = 1, nb = 1, priorParameters = priorParameters,
      restriction = "logOR", delta = 60, gridSize = 101
    ),
    "too extreme"
  )
  testthat::expect_error(
    restrictedThetaSupport("logOR", delta = 800, gridSize = 101),
    "too extreme"
  )
  moderate <- turnerEProcess(
    c(0, 1), c(1, 0), na = 1, nb = 1, priorParameters = priorParameters,
    restriction = "logOR", delta = 25, gridSize = 101
  )
  testthat::expect_true(all(is.finite(moderate)))
})

testthat::test_that("conditional block functions recycle scalar group sizes", {
  ya <- c(2, 1, 3, 2)
  yb <- c(1, 2, 1, 2)
  na <- 4
  nb <- 4

  for (compute in list(computeNml, seqCond, computeEGaussGrid)) {
    testthat::expect_equal(
      compute(ya, yb, na, nb),
      compute(ya, yb, rep(na, 4), rep(nb, 4))
    )
  }
  testthat::expect_error(
    computeNml(ya, yb, rep(na, 3), nb),
    "equal lengths"
  )
})

testthat::test_that("the conditional entry point speaks propDiff/logOR", {
  ya <- c(2, 1, 3, 2)
  yb <- c(1, 2, 1, 2)

  fixed <- conditionalEValueFixedAlternative(
    ya, yb, na = 4, nb = 4, logOR = 1, alternative = "greater"
  )
  testthat::expect_length(fixed, 4)
  testthat::expect_error(
    conditionalEValueFixedAlternative(ya, yb, na = 4, nb = 4),
    "No logOR is given"
  )

  futility <- saviFutilityTwoPropCondStat(
    ya, yb, na = 4, nb = 4, logOR = 1, alternative = "twoSided"
  )
  testthat::expect_equal(
    futility$eValueVec,
    pmax(
      conditionalEValueFixedAlternative(
        ya, yb, na = 4, nb = 4, logOR = 1, alternative = "greater"
      ),
      conditionalEValueFixedAlternative(
        ya, yb, na = 4, nb = 4, logOR = 1, alternative = "less"
      )
    )
  )

  # computeConfidenceInterval2x2 used to error on its own default, because it
  # switched on the unmatched choices vector.
  testthat::expect_equal(
    computeConfidenceInterval2x2(ya, yb, na = 4, nb = 4),
    c(-Inf, Inf)
  )
  testthat::expect_equal(
    computeConfidenceInterval2x2(
      ya, yb, na = 4, nb = 4, alternative = "greater"
    ),
    c(0, Inf)
  )
})

testthat::test_that("deprecated 2x2 wrappers forward to the new signatures", {
  # designSafeTwoProportions used to pass M= to a function without that formal
  # (hard error); the two test wrappers passed logOddsConfidenceSearchBounds and
  # pilot, which the new signatures silently swallowed via `...`.
  design <- suppressWarnings(designSafeTwoProportions(
    na = 1, nb = 1, nBlocksPlan = 1,
    alternativeRestriction = "logOddsRatio", M = 2
  ))
  testthat::expect_s3_class(design, "saviDesign")
  testthat::expect_equal(design$effectMeasure, "logOR")
  testthat::expect_true(design$pilot)
  testthat::expect_warning(
    designSafeTwoProportions(na = 1, nb = 1, nBlocksPlan = 1, M = 2),
    "deprecated"
  )

  ya <- c(1, 0, 1)
  yb <- c(0, 0, 1)
  testthat::expect_warning(
    safeTwoProportionsTest(ya, yb, designObj = design), "deprecated"
  )
  testthat::expect_warning(
    safe.prop.test(ya, yb, designObj = design), "deprecated"
  )
  viaWrapper <- suppressWarnings(
    safeTwoProportionsTest(ya, yb, designObj = design,
                           logOddsConfidenceSearchBounds = c(0.1, 3))
  )
  viaAlias <- suppressWarnings(safe.prop.test(ya, yb, designObj = design))
  direct <- suppressWarnings(
    saviTwoProportionsTest(ya, yb, designObj = design)
  )
  testthat::expect_s3_class(viaWrapper, "saviTest")
  testthat::expect_equal(viaWrapper$eValue, direct$eValue)
  testthat::expect_equal(viaAlias$eValue, direct$eValue)
})
