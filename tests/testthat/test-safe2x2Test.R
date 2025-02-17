test_that("bernoulliMLTwoProportions", {
  totalSuccess <- 15
  totalFail <- 26
  priorSuccess <- 0.45555
  priorFail <- 0.23232

  expectedResult <- (totalSuccess + priorSuccess) / (totalSuccess + totalFail + priorSuccess + priorFail)
  result <- safestats:::bernoulliMLTwoProportions(totalSuccess, totalFail, priorSuccess, priorFail)
  expect_equal(result, expectedResult)
})

test_that("no difference, two-sided test, na = nb", {
  expectedResult <- c(1, 1, 1)
  result <- vector(mode = "numeric", length = length(expectedResult))

  # all 1s, e = 1
  ya <- rep(1, 10)
  yb <- rep(1, 10)
  result[1] <- safestats::EValTwoProp(ya, yb,
    alternative = "twoSided",
    na = 1, nb = 1, esType = "none"
  )

  # all 0s, e = 1
  ya <- rep(0, 10)
  yb <- rep(0, 10)
  result[2] <- safestats::EValTwoProp(ya, yb,
    alternative = "twoSided",
    na = 1, nb = 1, esType = "none"
  )

  # all 0s, e = 1
  ya <- rep(10, 10)
  yb <- rep(100, 10)
  result[3] <- safestats::EValTwoProp(ya, yb,
    alternative = "twoSided",
    na = 10, nb = 100, esType = "none"
  )

  expect_equal(result, expectedResult)
})

test_that("absolute difference, greater, na = nb", {
  expectedResult <- c(1, 1)
  result <- vector(mode = "numeric", length = length(expectedResult))

  # all 1s, e = 1
  ya <- rep(1, 10)
  yb <- rep(1, 10)
  result[1] <- safestats::EValTwoProp(ya, yb,
    alternative = "greater",
    na = 1, nb = 1, esType = "difference", esMin = 0.2
  )

  # all 0s, e = 1
  ya <- rep(0, 10)
  yb <- rep(0, 10)
  result[1] <- safestats::EValTwoProp(ya, yb,
    alternative = "greater",
    na = 1, nb = 1, esType = "difference", esMin = 0.2
  )

  expect_equal(result, expectedResult)
})

test_that("no difference, two-sided test, na != nb", {
  expectedResult <- c(1, 1)
  result <- vector(mode = "numeric", length = length(expectedResult))

  # all 1s, e = 1
  ya <- rep(2, 10)
  yb <- rep(20, 10)
  result[1] <- safestats::EValTwoProp(ya, yb,
    alternative = "twoSided",
    na = 2, nb = 20, esType = "none"
  )

  # all 0s, e = 1
  ya <- rep(0, 10)
  yb <- rep(0, 10)
  result[2] <- safestats::EValTwoProp(ya, yb,
    alternative = "twoSided",
    na = 20, nb = 2, esType = "none"
  )

  expect_equal(result, expectedResult)
})
