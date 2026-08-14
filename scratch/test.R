# Compare a full fixed grid with the adaptive final-interval search.
# Run from the safestats repository root:
#
# Rscript scratch/test.R

source("R/safe2x2Test.R")

nBlocks <- 50L
gridSize <- 10000L
repetitions <- 5L
recordPath <- "scratch/adaptive-final-interval-benchmark.csv"
plotPath <- "scratch/adaptive-final-interval-benchmark.png"

# A grid of this size has this spacing over the feasible interval (-1, 1).
# Use the same value as the adaptive endpoint-bracket tolerance.
comparisonPrecision <- 2 / (gridSize + 1L)

set.seed(20260814)
ya <- stats::rbinom(nBlocks, size = 10, prob = 0.9)
yb <- stats::rbinom(nBlocks, size = 1, prob = 0.1)
design <- designSaviTwoProportions(
  na = 10,
  nb = 1,
  nBlocksPlan = nBlocks,
  alpha = 0.05
)

runFixedGrid <- function() {
  computeConfidenceSequenceForDifferenceTwoProportions(
    ya = ya,
    yb = yb,
    gridSize = gridSize,
    saviDesign = design
  )
}

runAdaptive <- function() {
  computeAdaptiveFinalIntervalForDifferenceTwoProportions(
    ya = ya,
    yb = yb,
    precision = comparisonPrecision,
    saviDesign = design
  )
}

# Warm both methods before timing so source loading and byte compilation do not
# affect the comparison.
fixedGridSequence <- runFixedGrid()
adaptiveInterval <- runAdaptive()

timeMethod <- function(method) {
  vapply(seq_len(repetitions), function(repetition) {
    system.time(invisible(method()))[["elapsed"]]
  }, numeric(1))
}

fixedGridTimes <- timeMethod(runFixedGrid)
adaptiveTimes <- timeMethod(runAdaptive)
fixedGridFinalInterval <- tail(
  fixedGridSequence[c("lowerBound", "upperBound")],
  n = 1L
)
adaptiveCandidateDiagnostics <- adaptiveInterval[["candidateDiagnostics"]]

benchmarkRecord <- rbind(
  data.frame(
    method = "Fixed grid",
    repetition = seq_len(repetitions),
    elapsedSeconds = fixedGridTimes,
    candidateProcesses = gridSize,
    candidateBlocks = gridSize * nBlocks,
    lowerBound = fixedGridFinalInterval[["lowerBound"]],
    upperBound = fixedGridFinalInterval[["upperBound"]],
    status = "grid"
  ),
  data.frame(
    method = "Adaptive final interval",
    repetition = seq_len(repetitions),
    elapsedSeconds = adaptiveTimes,
    candidateProcesses = nrow(adaptiveCandidateDiagnostics),
    candidateBlocks = sum(adaptiveCandidateDiagnostics[["blocksEvaluated"]]),
    lowerBound = adaptiveInterval[["lowerBound"]],
    upperBound = adaptiveInterval[["upperBound"]],
    status = adaptiveInterval[["status"]]
  )
)
utils::write.csv(benchmarkRecord, recordPath, row.names = FALSE)

methodOrder <- c("Fixed grid", "Adaptive final interval")
medianSeconds <- vapply(methodOrder, function(methodName) {
  stats::median(benchmarkRecord[["elapsedSeconds"]][
    benchmarkRecord[["method"]] == methodName
  ])
}, numeric(1))
candidateBlocks <- vapply(methodOrder, function(methodName) {
  benchmarkRecord[["candidateBlocks"]][
    benchmarkRecord[["method"]] == methodName
  ][1L]
}, numeric(1))
speedup <- medianSeconds[["Fixed grid"]] /
  medianSeconds[["Adaptive final interval"]]

print(benchmarkRecord)
cat(
  "\nMedian speedup (fixed grid / adaptive): ",
  sprintf("%.2fx", speedup),
  "\nEquivalent fixed-grid spacing: ",
  format(comparisonPrecision, scientific = TRUE),
  "\nRecord written to ", recordPath,
  "\nPlot written to ", plotPath,
  "\n",
  sep = ""
)

grDevices::png(plotPath, width = 1800, height = 650, res = 150)
graphics::par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))

speedBars <- graphics::barplot(
  medianSeconds,
  names.arg = c("Fixed\ngrid", "Adaptive\nfinal"),
  col = c("#999999", "#2166AC"),
  ylim = c(0, max(medianSeconds) * 1.15),
  ylab = "Median elapsed seconds",
  main = sprintf("Speed: %.2fx faster", speedup),
  sub = paste(repetitions, "repetitions")
)
graphics::text(
  speedBars,
  medianSeconds,
  labels = sprintf("%.3fs", medianSeconds),
  pos = 3
)

workBars <- graphics::barplot(
  candidateBlocks,
  names.arg = c("Fixed\ngrid", "Adaptive\nfinal"),
  col = c("#999999", "#2166AC"),
  ylim = c(0, max(candidateBlocks) * 1.15),
  ylab = "Candidate-block evaluations",
  main = "RIPr work"
)
graphics::text(
  workBars,
  candidateBlocks,
  labels = format(candidateBlocks, big.mark = ",", trim = TRUE),
  pos = 3
)

intervals <- rbind(
  c(
    lower = fixedGridFinalInterval[["lowerBound"]],
    upper = fixedGridFinalInterval[["upperBound"]]
  ),
  c(
    lower = adaptiveInterval[["lowerBound"]],
    upper = adaptiveInterval[["upperBound"]]
  )
)
graphics::plot(
  NA_real_,
  xlim = c(-1, 1),
  ylim = c(0.5, 2.5),
  yaxt = "n",
  xlab = expression(theta[B] - theta[A]),
  ylab = "",
  main = "Final interval"
)
graphics::axis(2, at = c(1, 2), labels = c("Fixed grid", "Adaptive"), las = 1)
graphics::segments(intervals[, "lower"], c(1, 2), intervals[, "upper"], c(1, 2),
  lwd = 4, col = c("#666666", "#2166AC")
)
graphics::points(intervals[, "lower"], c(1, 2), pch = 16, col = c("#666666", "#2166AC"))
graphics::points(intervals[, "upper"], c(1, 2), pch = 16, col = c("#666666", "#2166AC"))
graphics::abline(v = 0, lty = 2, col = "#BBBBBB")

grDevices::dev.off()
