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
  expected <- c(1, 1, 0.9929555)
  actual <- vector(mode = "numeric", length = length(expected))

  # all 1s, e = 1
  ya <- rep(1, 10)
  yb <- rep(1, 10)
  actual[1] <- safestats::EValTwoProp(ya, yb,
    alternative = "twoSided",
    na = 1, nb = 1, esType = "none"
  )

  # all 0s, e = 1
  ya <- rep(0, 10)
  yb <- rep(0, 10)
  actual[2] <- safestats::EValTwoProp(ya, yb,
    alternative = "twoSided",
    na = 1, nb = 1, esType = "none"
  )

  # unbalanced, na = 1, nb = 100
  actual[3] <- safestats::EValTwoProp(ya, yb,
    alternative = "twoSided",
    na = 1, nb = 4, esType = "none"
  )

  expect_equal(actual, expected)
})


test_that("no difference, two-sided test, na = nb original", {
  expected <- c(1, 1, 0.9929555)
  actual <- vector(mode = "numeric", length = length(expected))

  designObj <- safestats::designSafeTwoProportions(
    na = 1, nb = 1, alternativeRestriction = "none", alpha = 0.05, M = 1, nBlocksPlan = 10
  )

  # all 1s, e = 1
  ya <- rep(1, 10)
  yb <- rep(1, 10)
  actual[1] <- safestats::safeTwoProportionsTest(ya, yb, designObj = designObj)$eValue

  # all 0s, e = 1
  ya <- rep(0, 10)
  yb <- rep(0, 10)
  actual[2] <- safestats::safeTwoProportionsTest(ya, yb, designObj = designObj)$eValue

  # unbalanced, na = 1, nb = 100
  designObj <- safestats::designSafeTwoProportions(
    na = 1, nb = 4, alternativeRestriction = "none", alpha = 0.05, M = 1, nBlocksPlan = 10
  )
  actual[3] <- safestats::safeTwoProportionsTest(ya, yb, designObj = designObj)$eValue

  expect_equal(actual, expected)
})
