# Two-proportion e-process: requirements

**This document is the reference for `R/safe2x2Test.R`.** It is the centre of
the implementation: every change to the code is reflected here in the same
commit, and any disagreement between the two is a bug in the code, not here.

Statements marked **[stated]** are what was asked for directly; **[inferred]**
are my reading of the intent and should be corrected if wrong.

Status legend: **done** · **deferred** (agreed to leave for now) · **open**
(identified, no decision yet).

---

## 0. Notation and conventions

| Symbol | Meaning |
| --- | --- |
| `thetaA`, `thetaB` | Success probabilities, floating point in `(0, 1)`. Boundaries are ignored in most cases. |
| `na`, `nb` | Integer group sizes per data block. |
| `L` | `length(ya)`, the number of observed blocks. |
| `ya`, `yb` | Number of ones in group A / B in each block. |
| `delta` | The effect. Either `propDiff` or `logOR`. |

Effect measures, both signed as **B minus A**:

- `propDiff = thetaB - thetaA`, in `(-1, 1)`.
- `logOR = logit(thetaB) - logit(thetaA)`, on the whole real line.

**R0.2 — No zero in any candidate grid.** **[stated]** Neither confidence
sequence offers `0` as a candidate effect, so a reported bound of exactly `0`
is impossible and a finite bound always falls strictly on one side of the
null. See R2.4 and R3.2 for how each grid achieves this.

`gridSize` is the resolution of the discrete grid used to learn a restricted
theta; `nSimulations` and `maxBlocks` size the simulation.

**R0.1 — One vocabulary.** **[stated]** `propDiff` and `logOR` are used
throughout. The legacy spellings `difference`, `linearDifference` and
`logOddsRatio` are accepted **only** by `designSaviTwoProportions()`, which
maps them on entry; they are never stored or passed on. Every other function
validates with `match.arg()` against the canonical names and rejects the old
ones. The `normalizeTwoProportionEffect()` translator this used to require has
been deleted. **done**

**R0.3 — One naming scheme across both files.** **[stated]**
`R/safe2x2TestCond.R` follows `R/safe2x2Test.R`, not its own conventions:

| Was | Now |
| --- | --- |
| `logOddsRatio` (argument), `delta`, `lOR` | `logOR` |
| `mIter` | `nSteps` |
| `t` (loop index) | `block` |
| `n1` | `totalSuccesses` |
| `deltaGrid` | `logORGrid` (also the name in the returned list) |
| `logPriorWeights` / `currentLogPostWeights` | `priorLogWeights` / `logPosteriorWeights` |
| `logLR` | `logLikelihoodRatio` |
| `yaCumsum`, `prevYaCumsum` | `cumulativeYa`, `previousYa` |
| `ya_up`, `yb_up` (snake_case) | `candidateCumulativeYa`, `candidateCumulativeYb` |
| `lb`, `ub`, `kVals` | `feasibleLower`, `feasibleUpper`, `feasibleSuccesses` |
| `deltaMle`, `oddsMle`, `kDeltaMle` | `mleLogOR`, `mleOdds`, `candidateMleLogOR` |
| `logNum`, `logDenHg`, `logDiff` | `logNumerator`, `logDenominator`, `logLikelihoodRatio` |
| `normConst` | `logNormalisingConstant` |
| `sPlus0`, `sMin0` | `greaterEValues`, `lessEValues` |
| `objFunc`, `res`, `previous`, `logH`, `trivialConfInt` | `objectiveFunction`, `result`, `previousBlock`, `logBaseMeasure`, `trivialConfidenceInterval` |

`alternative` choices are ordered `c("twoSided", "greater", "less")`
everywhere in the file, and `sapply()` is replaced by `vapply()` on the two
grid loops. **done**

Two crashes surfaced while doing this and were fixed:

- `computeNml()` did not recycle a scalar `na` / `nb`, unlike every other
  function in both files, so `cumsum(na)` stayed length one and every block
  after the first indexed `NA`. It now recycles and length-checks with the
  same idiom as `seqCond()` and `computeEGaussGrid()`.
- `computeConfidenceInterval2x2()` never called `match.arg(alternative)`, so
  calling it with its own default switched on a length-three vector and
  errored.

---

## 1. Numerical helpers

**R1.1 — Convert between the two thetas under a fixed effect.** **[stated]**
Given `thetaA`, the effect type and its value, return `thetaB`.
→ `thetaBFromRestriction()` handles both effects directly; the one-caller
`thetaBFromLogOR()` wrapper was removed. **done**

**R1.2 — Learn `thetaA`, `thetaB` from `ya, yb, na, nb`.** **[stated]**

- **R1.2.0 Unrestricted:** two independent Beta priors, one per stream.
- **R1.2.1 `propDiff`:** start from a Beta prior on `thetaA`, do Bayesian
  updating on a discrete grid; the prior is updated *after* evaluating the
  log likelihood of the block just seen.
- **R1.2.2 `logOR`:** identical, except for how `thetaB` is derived from
  `thetaA`.
- **R1.2.3 Naming and decomposition:** rename `turnerPredictiveThetas` to
  something informative (suggested `learnThetaWithRestriction`) and split it
  into two or three functions — one unrestricted, one per effect measure.

→ **done**, consolidated as:

| Function | Role |
| --- | --- |
| `restrictedThetaSupport()` | The one-dimensional support curve for a fixed effect. |
| `restrictedPriorLogWeights()` | Normalised Beta log prior on that curve. |
| `learnPredictiveThetas()` | Independent Beta posterior means when unrestricted; grid posterior when restricted. |

The one-caller `learnUnrestrictedThetas()` and `learnRestrictedThetas()`
helpers were removed. Their formulas remain visibly separated by a guard
clause inside `learnPredictiveThetas()`. **done**

**Predictability is a hard requirement.** The pair used for block `i` must be
computed from blocks `1 … i-1` only.

**R1.4 — Restricted supports must stay strictly interior.** **[stated]**
`plogis()` saturates to exactly `1` once its argument reaches about 37, so a
large enough `logOR` puts a support point at exactly `0` or `1`, where a zero
count gives `0 * -Inf = NaN` that then spreads through the posterior
normalisation. `restrictedThetaSupport()` now refuses to build such a curve
and says which `delta` and `gridSize` caused it, rather than returning `NaN`
(`turnerEProcess()`) or dying at an unrelated `if` (the sampler). **done**

**R1.3 — `resolveTurnerPriorParameters`.** **[stated]** Keep or drop; if kept,
rename to `resolveBetaPriorParameters`. → kept and renamed. **done**

**R1.4 — Coefficient-free likelihood calculations.** **[stated]** Binomial
coefficients cancel from every likelihood ratio and from normalized restricted
posterior weights. `logLikelihoodRatioIncrements()` therefore computes

$$
y\{\log p-\log q\}
+(n-y)\{\log(1-p)-\log(1-q)\}
$$

for each group directly, without `dbinom()`. A term with zero count is defined
as zero, including at probability boundaries. When numerator and denominator
probabilities are equal, their group log ratio is defined as zero, including
on a shared probability-zero event. Restricted posterior and Beta-prior grid
weights likewise omit normalizing constants that cancel when the weights are
normalized. **done**

---

## 2. Proportion difference

**R2.1 — Cubic solver.** **[stated]** Given `theta`, `na`, `nb` and the
difference, solve the third-degree polynomial and return the solution.
→ `solveOnePropDiffRIPr()` solves one block; `solvePropDiffRIPr()` is the
vectorised wrapper over blocks (`mapply`). The scalar solver exists so the
confidence sequence (R2.3) can call it directly inside its block loop without
the wrapper's validation and recycling overhead, which cost ~3x when the
wrapper was called once per block per candidate. **done** 2026-09-18

**R2.2 — Grid of e-processes.** **[stated]** Learn the unrestricted theta, find
the null thetas for every candidate via the solver, return a `gridSize × L`
matrix. → `calculateEValuesForPropDiffGrid()`, returning `L × gridSize`
(blocks as rows). Since 2026-09-18 the confidence sequence (R2.3) no longer
calls it; it remains as the diagnostic view of the complete e-process matrix
and as the reference the tests compare the sequence against. The candidate
grid itself comes from `propDiffCandidateGrid()`, shared with R2.3. **done**
(pre-existing)

**R2.3 — Confidence interval at each time.** **[stated]** A wrapper returning a
bound at every one of the `L` blocks.
→ `computeConfidenceSequenceForPropDiffTwoProportions()`. Each candidate is a
two-sided point null tested at `alpha`. The function walks the data block by
block: at block `t` it takes the predictor learned from blocks before `t`,
and for every candidate still in the confidence set solves the RIPr
(R2.1, scalar), adds one log-likelihood-ratio increment to that candidate's
running log e-process, and drops the candidate once the total reaches
`log(1/alpha)`. The bounds at block `t` are the smallest and largest
candidates remaining, `-1`/`1` while an edge candidate remains, `NA` when
none does. This is a direct transcription of the running-intersection
definition, and it never revisits a rejected candidate, so it skips 60–90 %
of the `polyroot()` solves the former matrix version performed (P1).
`runningIntersection = FALSE` uses the same loop but updates every candidate
at every block and reads each block's set from the cumulative values at that
block only, with no memory of earlier rejections, as the old
`calculateEValuesForLinearDeltaGrid()` flag did; it is kept for inspection
only and is `TRUE` by default. The former matrix version (all `L × gridSize`
e-processes via R2.2, then first rejection blocks per column) gave identical
bounds and is gone; R2.2 remains as the reference the tests compare against.
The propDiff construction (one two-sided family) and the log-odds-ratio
construction (two one-sided families, see §3) are different and are kept
separate on purpose. **done** 2026-09-18

**R2.4 — Candidates must exclude zero.** **[stated]** The old grid contained `0`
only when `gridSize` was odd, so whether the null value was testable depended
on the parity of a resolution argument. Resolved by *excluding* zero at every
parity: `gridSize` equally spaced candidates in `(0, 1)` are mirrored about
zero, so the grid holds `2 * gridSize` values and a bound of exactly `0` can
never be reported. `gridSize` therefore counts candidates per side, exactly as
`confidenceBoundGridPrecision` does for the log-odds-ratio grid (R3.2). **done**

---

## 3. Log odds ratio

**R3.1 — Same overall structure as §2.** **[stated]** Solver → grid of
e-processes → confidence-sequence wrapper.
→ `solveLogORRIPr()` (quadratic in the base group's odds, solved in log space
via `logPositiveQuadraticRoot()`), `calculateEValuesForLogORGrid()`,
`computeConfidenceSequenceForLogORTwoProportions()`. **done** (pre-existing)

**R3.1a — The log-space quadratic solve must read as one.** The mean-matching
condition `na*thetaA + nb*thetaB = na*numeratorThetaA + nb*numeratorThetaB`
under the log-odds-ratio constraint is a plain quadratic, and the code should
look like it. `logPositiveQuadraticRoot()` spelled its two log-sum-exp steps
out inline with `max`/`exp`/`log` arithmetic, which buried that. The repeated
step is now a named pairwise helper, `logAddExp(x, y)`, and the choice between
the two quadratic formulas is a single `ifelse()` on the sign of `b` rather
than assignment into a preallocated vector. `solveLogORRIPr()` keeps its
base/shifted parameterisation — it is what bounds the coefficients for a large
`abs(logOR)`, not a stylistic choice — but its derivation comment now states
the mean-matching equation first and the numerical reasoning second.

`base::polyroot()` was considered and rejected: it takes linear-space
coefficients, and `coefA` and `coefC` here are exactly the ones that underflow
to zero in linear space while staying finite as logs.

Numerically the refactor is a no-op except for one deliberate change:
`log(1 + z)` became `log1p(z)`, the correctly-rounded primitive. Roots are
bit-identical across 20k extreme coefficient draws; RIPr thetas move by at
most one ULP (2.2e-16), and the score equation holds to 2.3e-15 relative.
**done**

**R3.1b — `logAddExp()` is a sibling of `logSumExp()` (A11).** The pairwise
helper sits beside its caller in `R/safe2x2Test.R` while the reducing
`logSumExp()` sits in `R/safe2x2TestCond.R`. Both belong in the shared helper
file that A11 asks for; whoever does that move should take both. **open**

**R3.2 — Candidates must exclude zero.** **[stated]** R2.4 applies here too.
The grid used to insert `0` explicitly between the mirrored halves; that
insertion is removed, so the grid is now
`c(-rev(positiveGrid), positiveGrid)` and holds
`2 * confidenceBoundGridPrecision` candidates. Both search bounds are already
validated positive, so the positive half never reaches zero on its own and a
finite bound always sits strictly on one side of the null. **done**

**R3.2a — Confidence bounds are read off inline.**
`computeConfidenceSequenceForLogORTwoProportions()` inverts two one-sided
families, each at `alpha / 2`: the lower family tests `logOR <= candidate`,
the upper family `logOR >= candidate`, and both run over the full signed
grid. As in R2.3, a candidate is rejected by a family from the first block at
which its log e-process reaches `log(2/alpha)` onward (running intersection).
At each block the lower bound is the smallest candidate the lower family has
not yet rejected and the upper bound the largest candidate the upper family
has not yet rejected; a bound is `-Inf`/`Inf` while the outermost candidate
on its side remains and `NA` once a family has rejected every candidate. The
wrapper used to hand its two matrices to `confidenceBoundsFromLogEProcesses()`,
whose input validation, running-maximum matrices and sentinel handling made
the inversion harder to read than the rule above; that helper is removed and
the inversion sits in the wrapper next to the grid construction, mirroring the
propDiff wrapper. Outputs are identical. **done** 2026-09-18

**R3.3 — Sequential conditional Gaussian mixture.** `computeEGaussGrid()`
represents a standard-normal prior on the A-minus-B conditional log odds ratio
on a fixed grid. At block `t`, it averages the fixed-parameter conditional
likelihood ratio under weights based only on blocks `1, ..., t - 1`, returns
that blockwise e-factor, and then updates the weights with block `t`. The old
`computeEGauss()` calculation, which restarted from the same prior at every
block, is removed. **done**

---

## 4. Stopping-time simulation

This is the part flagged as both **the biggest time sink** and **incorrect**.

**R4.1 — Theta grid per effect.** **[stated]** For a given `propDiff` or
`logOR`, generate a grid of feasible `(thetaA, thetaB)` pairs.
→ `makeSimulationThetaGrid()`, which now also returns the `restriction` and
`delta` so they travel with the grid. **done**

**R4.2 — Simulate per pair.** **[stated]** For each theta pair, simulate
`nSimulations` trajectories up to `maxBlocks`. **done**

**R4.3 — Record stopping times and e-values.** **[stated]** The stopping time is
the index at which the process crosses `1 / alpha`; record the e-value there.
**A path that never crosses gets stopping time `Inf`** (not the horizon).
**done**

**R4.4 — The numerator must be restricted.** **[stated]** The core correctness
requirement. The process was built on the *unrestricted* case, with thetas
plugged in by `turnerLogEProcessChunk`. It must instead:

1. learn theta **under the restriction** given by the effect type and value,
2. compute `denominatorTheta` as the pooled predictable probability (as before),
3. form the log-likelihood ratio from those.

`simulateTurnerStoppingGrid()` is therefore restricted-only: its `restriction`
argument accepts exactly `"propDiff"` or `"logOR"`, never `"none"`. Its
unrestricted code path has been removed. `turnerEProcess()` continues to
support `"none"` outside this worst-case restricted simulation. **done**

Chunk-by-chunk evaluation was judged not possible, or too complicated, under
this requirement. **[stated]** → confirmed: `turnerLogEProcessChunk` and its
`chunkSize` argument were removed, and the simulator now advances one block at
a time. **done**

Verified against the legacy `upstream/futility88:R/safe2x2Test.R`'s `calculateSequential2x2E()` — final
e-values agree to ~1e-15 for both restricted effect measures.

**R4.5 — Do not simulate past the stopping time.** **[stated]** Out of one theta
pair, one effect and `maxBlocks`, it is very likely the process ends well
before the horizon. Neither the data generation nor the e-process evaluation
should run past a path's own crossing. → paths are dropped from the active set
the moment they cross. **done**

**R4.6 — Accuracy before speed.** **[stated]** A `for` loop is acceptable for
now; correctness first. → the block loop is sequential; paths are vectorised
within a block. **done**

**R4.7 — Direct sufficient-statistic posterior update; no learner object.**
**[stated]** The prior update must remain visible in the simulation rather than
being hidden behind a stateful learner abstraction. → `newTurnerThetaLearner()`
was deleted and the calculation is written directly in
`simulateTurnerStoppingGrid()`. **done**

The restricted posterior depends on the data **only through
`(blocks, sum ya, sum yb)`** — the binomial coefficients cancel in the
normalisation, so block order is irrelevant. With `t` completed blocks,
`Sa = sum(ya)` and `Sb = sum(yb)`, the support-point weights are

$$
\begin{aligned}
\log W_g ={}& \log \pi_g
  + S_A \log \theta_{A,g}
  + (t n_A-S_A)\log(1-\theta_{A,g}) \\
 &+ S_B \log \theta_{B,g}
  + (t n_B-S_B)\log(1-\theta_{B,g}).
\end{aligned}
$$

The simulation stores only cumulative successes for active paths. It
constructs and normalises this log likelihood immediately before generating
the next block, so the numerator for block `t + 1` uses blocks `1, ..., t`
only. The new observations are added to the cumulative counts only after that
block's e-factor is evaluated. There is no persistent `gridSize × nPaths`
posterior matrix.

Paths sharing a state are evaluated once:

- `propDiff` keys on the pair `(sum ya, sum yb)`.
- `logOR` keys on `sum ya + sum yb` alone. The extra collapse is valid because
  the term that distinguishes them, `sum yb · delta`, is constant across the
  support grid and cancels in the normalisation.

**R4.7.1 — Keep the likelihood update readable.** **[stated]** Both effects
must use the same four-term log-likelihood formula displayed above. The
`logOR` total-success optimization applies only when constructing the state
key; one representative `(sum ya, sum yb)` pair from each state is evaluated
with the common formula. Do not introduce a separate algebraic likelihood
branch for `logOR`, and do not hide the update behind a stateful object.
**done**

The restricted branch of `learnPredictiveThetas()` first calculates the
posterior mean of `thetaA`, then obtains `thetaB` with
`thetaBFromRestriction()`. For `propDiff` this is also the posterior mean of
`thetaB`. For `logOR`, it is a restricted plug-in value and is generally not
the separate marginal posterior mean of `thetaB`; changing that would be a
statistical change and requires an explicit decision here first.

Verified independently against `turnerEProcess()` at both signs of `delta`,
both effect measures and unequal block sizes: relative differences 0 to 7e-15.

**R4.8 — Explicit failure when the horizon is too short.** **[stated]** When too
few paths cross for the `1 - beta` quantile to exist, say so plainly: the
worst-case stopping time cannot be found at this `maxBlocks`, try increasing
it. → warning raised from `simulateWorstCaseStoppingTimes()`, naming the
horizon, the observed crossing fraction, the worst-case baseline and the
fraction needed. `maxBlocks` was added to `designSaviTwoProportions()` so the
advice is actionable from there. **done**

---

## 5. Testing functions

**R5.1 — `saviTwoProportionsTest()`.** **[stated]** A test entry point matching
the reference's `# Testing fnts` section, with an S3 generic, a `.default`
method taking `ya, yb, na, nb`, a `.formula` method for `success ~ group`, and
a `savi.prop.test` alias. Built on `constructSaviTestObj("Two Proportions")`.
**done**

- The e-process comes from `turnerEProcess()`, restricted whenever the design
  carries an `esMin`, so the test runs exactly the process the design planned.
- `wantConfidenceSequence = TRUE` dispatches to the propDiff or logOR
  confidence sequence according to the design's `effectMeasure`; both were
  previously unreachable from any entry point.
- `estimate` and the confidence sequence are reported on the design's own
  scale, so they and `esMin` are the same quantity.
- No `statistic` is set: the e-value is the test's summary and is printed on
  its own line.
- The result carries `posteriorHyperParameters`, the Beta posterior after all
  observed blocks. `designSaviTwoProportions()` has always accepted a
  `previousSaviTestResult` and read that field, but nothing produced it; the
  loop now closes.

---

## 6. Deferred and open items

| Ref | Item | Decision |
| --- | --- | --- |
| S1 | `gridSize` names three different things (posterior resolution, candidate count, baseline count) alongside `thetaGridSize`, `confidenceBoundGridPrecision` and `effectGridSize`. | **open** |
| S2 | `designSaviTwoProportions()` is a ~220-line four-case branch; the package already splits these (`designSaviT1aWantNPlan`, `designSaviT2WantBeta`, …). | **open** |
| S3 | Two-stream input validation (`na`/`nb` recycling plus length checks) is repeated in six places. | **open** |
| S4 | `simulateWorstCaseStoppingTimes()` and `simulateWorstCasePower()` are near-duplicates: same grid, same simulator call, same worst-row bootstrap. | **open** |
| P1 | `solvePropDiffRIPr()` calls `polyroot()` once per block per candidate (120k calls for a 600-block, 200-candidate sequence), making the propDiff sequence ~5x slower than the logOR one. | **partly done** 2026-09-18 — the block-by-block sequence (R2.3) stops solving for a candidate once it is rejected, skipping 60–90 % of the solves; measured 2.6x (100 blocks) to 8x (2000 blocks) faster at 40 candidates, 5x at 500 blocks × 200 candidates, with `runningIntersection = FALSE` at parity (0.8–1.0x). One `polyroot()` per surviving (block, candidate) pair remains. |
| P2 | `confidenceBoundsFromLogEProcesses()` computed the running maximum twice when `upperLogEProcesses` defaulted to `lowerLogEProcesses`, as it did for propDiff. | **done** 2026-09-18 — helper removed; both sequences read their bounds off first rejection blocks inline (R2.3, R3.2a) |
| P3 | The propDiff simulation state key is built with `paste()` every block; an integer key is ~4x faster. | **open** |
| C1 | `simulateTurnerStoppingGrid()` no longer accepts `restriction = "none"`, so the unrestricted process can no longer be simulated for comparison. | **open** — deliberate? |
| — | `propDiff` inverts at `alpha`, `logOR` at `alpha/2` (Bonferroni). Different coverage semantics between the two sequences. | **deferred** — "not the major concern now" |
| C2 | `vignettes/contingency-tables-vignette.Rmd` still calls six removed functions (`simulateTwoProportions`, `simulateOptionalStoppingScenarioTwoProportions`, `simulateIncorrectStoppingTimesFisher`, `plotConfidenceSequenceTwoProportions`, `simulateCoverageDifferenceTwoProportions`, `computeConfidenceBoundForLogOddsTwoProportions`) and loads `savi2x2Sim` objects whose print/plot methods are gone, so `R CMD build` fails on vignettes. | **open** — rewrite or retire the vignette (known gap as of 2026-09-17) |
| C3 | Ten roxygen-documented helpers are not exported (all of `R/safe2x2TestCond.R` plus `logLikelihoodRatioProcess`, `calculateEValuesFor*Grid`, `computeConfidenceSequenceFor*TwoProportions`), which `R CMD check` flags. | **open** — decide per function between `@export` and `@noRd` |
| C4 | Deprecated wrappers in `R/deprecate.R` (`designSafeTwoProportions`, `safeTwoProportionsTest`, `safe.prop.test`) forward to the new signatures: `M` → `nSim`, `alternativeRestriction` → `effectMeasure` (`"none"` → `"propDiff"`, restriction then comes only from `delta`), `logOddsConfidenceSearchBounds` → `logORConfidenceSearchBounds`; `pilot`, `simThetaAMin`, `simThetaAMax` are dropped on the test side / ignored with a warning. Covered by a test. | **done** 2026-09-17 |
| — | `simulateMinimumDetectableEffect()` binary-searches on noisy evaluations; Monte Carlo noise can break the monotonicity the search assumes, and the returned effect carries no error estimate. | **deferred** — marked with a `TODO` in the source, left as is |

## 7. Alignment with upstream `R/tTest.R`

Upstream's `R/tTest.R` (`upstream/futility88:R/tTest.R`, byte-identical to the former `scratch/Reference_tTest.R`) is the structural template, minus its relevance
testing (`saviRelevanceTStatNEffNu` and the `relevanceTest` / `relevanceSize`
/ `alphaRelevance` arguments). Its section order is:

```
# Testing fnts ----                    stat, S3 generic/.default/.formula,
                                       savi.t.test alias, confidence interval
# Design fnts ----                     designSaviT + 1aWantNPlan / 2WantPower
                                       / 3WantEsMin / 3bWantParameter
# Batch design fnts ----               computeNPlanBatch..., computeMinEsBatch...
# Sampling functions for design ----   sampleStoppingTimes..., computePower...,
                                       computeNPlan...
# Helper fnts ----
# Data generating fnt ----
```

| Ref | Gap against the template | Status |
| --- | --- | --- |
| A1 | Testing section. | **done** — `saviTwoProportionsTest()` with `.default` and `.formula` methods plus a `savi.prop.test` alias, built on `constructSaviTestObj("Two Proportions")`. See R5. |
| A2 | Design object built through `constructSaviDesignObj()`. | **done** — added a `"Two Proportions"` branch carrying `betaPriorParameterValues`, `effectMeasure` and `alternativeRestriction`; the design now sets `eType = "turner"` and `designScenario`, and drops unfilled slots with `Filter(Negate(is.null), ...)`. |
| A3 | The four-case branch should split into `designSaviTwoProportions1aWantNPlan` / `2WantBeta` / `3WantEsMin`, matching `designSaviT1aWantNPlan` etc., and set `designScenario`. (Same as S2.) | **open** |
| A4 | Bootstrapping is reimplemented with `replicate(nBoot, quantile(sample(...)))` instead of `computeBootObj()` / `computeNPlanBootstrapper()` / `computeBetaBootstrapper()`. | **open** — proposal to review before implementing; depends on A5, since the shared bootstrappers read `breakVector` and `eValuesAtNMax`. |
| A5 | The sampler returns only `stoppingTimes` and `eValuesAtStopping`. The template's `constructSampleStoppingTimesList()` also carries `breakVector`, `eValuesStopped`, `eValuesAtNMax`, `samplePaths` and `stoppedVector`. `computeBetaBootstrapper()` expects `breakVector`; `plot.saviDesign()` expects `samplePaths`. | **open** |
| A6 | Sampling-function naming. | **done** — `sampleStoppingTimesSaviTwoProportions()`, `computePowerSaviTwoProportions()`, `computeNPlanSaviTwoProportions()`, `computeMinEsSaviTwoProportions()`. |
| A7 | Argument vocabulary. | **done** — `M` and `nSimulations` are now `nSim`, and the design function takes `nBoot`. `maxBlocks`, `nBlocksPlan`, `thetaGridSize` and `gridSize` are kept: they carry more information in this scope than `nMax`/`nPlan` would. |
| A8 | `alternative` is hardcoded to `"twoSided"`. | **declined** |
| A9 | Roxygen/export hygiene, `addCite()` references. | **declined** |
| A10 | No `generateTwoProportionData()` to match `generateNormalData()`; data is generated inline with `rbinom()` inside the sampler. No `pb`, `seed`, `wantSamplePaths` or `wantSimData` arguments. | **open** |

**Section order.** The template's flat section order is *not* adopted wholesale.
The existing split by effect measure — `# Proportion difference: RIPr
projection and confidence sequence ----` and `# Log odds ratio: RIPr
projection and confidence sequence ----` — is kept deliberately: each holds a
solver, a grid function and a confidence-sequence wrapper that belong
together. **[stated]** Only the template's section *names* and function naming
are borrowed.
| A11 | `logSumExp()` lives in the untracked `R/safe2x2TestCond.R` while `R/safe2x2Test.R` depends on it. It belongs in a shared helper file — together with the pairwise `logAddExp()` added in R3.1a. | **agreed, open** |
| A12 | Function names in `R/safe2x2TestCond.R` not yet aligned: `seqCond()` (its own title says "Sequential conditional plug-in E-values"), `computeConfidenceInterval2x2()`, `saviTwoPropCondStat()`. Variables and arguments are aligned (R0.3); function names were left alone. | **open** |

---

## 8. Features present in the legacy `upstream/futility88:R/safe2x2Test.R` but not in the current code

Not requirements yet — listed so the omissions are deliberate rather than lost.

- `thetaAMin` / `thetaAMax`: narrow the baseline grid with prior knowledge.
  (Its branch in scratch is inverted and never worked as intended.)
- `estimateImpliedTarget`: keep collecting past the stopping time to report
  `logImpliedTarget` at `nPlan`.
- `deltaDesign`: generate data at one effect and test at another, for
  misspecification studies. `simulateTurnerStoppingGrid()` can already do this
  — `restriction`/`delta` are separate arguments from `thetaA`/`thetaB` — but
  `simulateWorstCaseStoppingTimes()` ties them together and does not expose it.
- `expectedStopTime`: mean instead of quantile.
