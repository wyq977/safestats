# Agent instructions

## Objective

Rewrite the safe anytime-valid 2x2 (two-proportion) test as a small, direct
implementation. Priorities, in order: statistical validity, explicit
conventions, numerically stable code, focused tests, then everything else.

Behaviour of the older code is not a requirement. Keep it only when it is
known to be intentional and correct. The `cond` branch holds the previous
attempt and may be consulted with `git show cond:<path>`, never merged.

## Scope

In scope: `R/newsafe2x2Test.R` (all new code goes here; the legacy
`R/safe2x2Test.R` is left untouched until the rewrite replaces it), its
tests under `tests/testthat/`, its roxygen and generated `man/` pages, and
the smallest 2x2-related edits to `NAMESPACE`, `DESCRIPTION`,
`R/safeS3Methods.R` and `R/deprecate.R`.

Out of scope: the t-test, z-test, log-rank test, the general design
framework, unrelated S3 methods, vignettes unless asked, package-wide style.

## Workflow

Each feature is agreed before it is coded:

1. the request is described in plain words;
2. it is reviewed and questions are settled;
3. the agreed contract is recorded under **Decisions** below;
4. only then is the R code written (tests only when asked).

Do not implement anything not yet recorded under Decisions. If code and a
decision disagree, the decision is the reference; change the code or reopen
the decision, never silently diverge.

## Conventions

- Groups are `A` and `B`. `ya`, `yb` are success counts, `na`, `nb` group
  sizes. Counts are integer-valued, non-negative and at most the group size.
- Two effect measures, named exactly `propDiff` (difference of
  proportions) and `logOR` (log odds ratio). Both are signed **B minus A**
  and both are anchored on `thetaA`: group B is always derived from group A
  and the effect, never the other way round.
  - `thetaB = thetaA + propDiff`, so `propDiff = thetaB - thetaA`.
  - `thetaB = plogis(qlogis(thetaA) + logOR)`, so
    `logOR = logit(thetaB) - logit(thetaA)`.
  The legacy spellings `difference`, `linearDifference` and `logOddsRatio`
  are not used anywhere in new code.
- `alternative` takes `c("twoSided", "greater", "less")`, in that order.
- Distinguish a blockwise e-factor from the cumulative e-process in names,
  documentation and tests. Step `i` uses only data through step `i`.
- Work on the log scale wherever underflow or overflow is plausible. Handle
  zero-probability cases explicitly; never let `NaN` or `Inf - Inf` decide.
- Mirror the layout and naming of `R/tTest.R`: a pure `...Stat` function,
  then the S3 generic with `.default` and `.formula`, then the `savi.*.test`
  alias, then design and sampling functions.

## Testing

Write tests only when asked. Do not test basic scaffolding such as object
fields, default values or name assignment. When asked, test against the
public API, not the internals. Prefer small exact tests: a hand-computable
table, one block versus the first step of a multi-block call, `log = TRUE`
versus the log of the plain output, a group swap with the mirrored
alternative, boundary and impossible inputs, and seeded simulation output. No Cartesian test grids.
Never loosen a tolerance before the discrepancy is understood.

Run 2x2 tests with `Rscript -e 'testthat::test_local(".", filter = "2x2")'`.
Regenerate docs with `LANG=en_US.UTF-8 Rscript -e 'devtools::document()'`;
the C locale corrupts unrelated `.Rd` files.

## Usage script

`~/Downloads/local-run.R` (outside the repository) mimics how users will
call the 2x2 design and test. Keep it working when the public API changes,
and run it from the package root after `devtools::load_all()` as a quick
end-to-end check.

## R style

2-space indent, `<-`, spaces around operators, `camelCase`, `match.arg()`
for enumerated arguments, guard clauses over nesting, roxygen2 Markdown on
exported functions with effect direction and return semantics stated.
Comments explain statistical intent, not syntax.

## Decisions

Agreed contracts, one entry per feature, newest last. Each entry states the
function, its inputs, its output, and the convention it relies on.

### 1. Sample-size vocabulary

`nSim` is the number of simulated paths, `nBoot` the number of bootstrap
resamples. `nPlan = c(nBlocks, na, nb)`, with `nBlocks` first because
`plot.saviDesign` reads `nPlan[1]`. The 2x2 design has no `nMax`,
`highN` or `nPlanBatch` (left `NULL`), and no `thetaA` scenario.

### 2. Design and test object fields

Only `eType = "grow"` for now, and the effect is always `propDiff`; there
is no `effectMeasure` field. `designSavi2x2(propDiffmin, na, nb, nPlan,
alpha, power, h0, alternative, eType, priorHyperParameters)` returns a
`saviDesign` with `testName = "Two Proportions"`, `testType = "2x2"`,
`h0 = c(propDiff = h0)`, and:

- `priorHyperParameters`: `list(betaA1, betaA2, betaB1, betaB2)`, the
  success and failure shapes of the Beta priors on `thetaA` and `thetaB`.
  The default, all four `0.18`, lives in `constructSaviDesignObj("Two
  Proportions")`; the design function replaces it only when the argument is
  not `NULL`, and errors on other names.
- `parameter`: a length-one named string summarising the prior, for
  printing.

Reserved for when the conditional e-variable returns: its prior on `logOR`
defaults to mean `0` and sd `1`.

`constructSaviTestObj("Two Proportions")` gets no `sumStats` or
`eFactorVec`. Note that `modifyList()` drops `NULL` placeholders, so a
field declared `NULL` in a constructor is absent until a function sets it.

### 3. Test function, unrestricted two-sided case

`savi2x2Test(ya, yb, designObj = NULL)`: `ya`, `yb` are per-block success
counts in observation order; `na`, `nb`, `alternative`, `h0` and the prior
all come from `designObj` (`NULL` gives a pilot `designSavi2x2()` with a
warning). No `ciValue` or confidence sequence yet.

- Numerator for block `i`: Beta posterior means of `thetaA`, `thetaB` given
  blocks `1..i-1` only. Denominator: the common
  `(na * thetaA + nb * thetaB) / (na + nb)`, the projection onto
  `thetaA = thetaB`.
- `eValueVec` is the cumulative e-process, `eValue` its last element.
- `n = c(na, nb, nBlocks)` in totals; `estimate` holds both observed
  proportions and `propDiff` (B minus A); `posteriorHyperParameters` is the
  Beta posterior after the last block.
- Errors, until each is designed: non-`NULL` `esMin`, `alternative` other
  than `"twoSided"`, `h0 != 0`. The design rejects non-positive Beta shapes.

### 4. Confidence sequence for propDiff

`computeConfidenceInterval2x2PropDiff(ya, yb, na, nb, priorHyperParameters,
alpha, precision = 100)` inverts the test on `precision` equally spaced
candidates strictly inside `(-1, 1)`. Each candidate `propDiff` is a point
null with its own e-process: numerator the predictable Beta posterior mean
(shared helper `predictiveThetas2x2`), denominator its reverse information
projection onto `thetaB - thetaA = propDiff`, found by `uniroot` on the KL
derivative (`solveRIPr2x2PropDiff`). Running intersection: a candidate
leaves for good at `1/alpha`. Each run of consecutive non-rejected
candidates is one interval, and the confidence set is the union of the
intervals; min and max over the whole set are not taken, since that would
fill holes. Returns a matrix with columns `block`, `lowerBound`,
`upperBound`: one row per block without holes, as for the z-test, several
rows for a block with holes, none for a block with everything rejected.

`savi2x2Test(..., wantCi = TRUE)` stores, as the other tests do, a two-column
`nBlocks x 2` `confSeqMatrix` (`lowerBound`, `upperBound`; other code reads it
by position): row `i` is the outermost bounds of block `i`'s union, `NA` for a
fully rejected block. This hull contains the union, so coverage holds, and
the rows stay nested. The exact union is kept only for the last block, as
`confSeq` (a `k x 2` matrix of `lowerBound`, `upperBound`; `k = 1` without
holes), with `ciValue = 1 - alpha` (no separate `ciValue` argument).
`savi2x2TestStat` returns the cumulative e-process, on the log scale when
`log = TRUE`.

### 5. Plotting with plot.saviTest

`plot.saviTest` works on 2x2 results unchanged in its general logic:

- `savi2x2Test` sets `n1Vec = seq_len(nBlocks)`, the block index, as the
  x-axis (legacy name the plot reads).
- `constructSaviDesignObj("Two Proportions")` defaults `relevanceTest =
  FALSE`; `NULL` makes the plot's `&&` error.
- The label switches map `"Two Proportions"` to `"Number of blocks"` (x) and
  `"propDiff"` (confidence-sequence y).
- For `wantConfSeqPlot = TRUE` the plot reads `confSeqMatrix` like any
  other test's; no 2x2 branch. The hull per block (Decision 4) fills holes
  in the picture; `confSeq` stays the exact union.

### 6. Restricted alternative on propDiff

Only `propDiff` restrictions; `logOR` is out of scope. `designSavi2x2`
accepts `propDiffMin` as `NULL` or one number strictly inside `(0, 1)`,
stored as `esMin`. Allowed combinations; everything else errors:

- `propDiffMin = NULL`, `"twoSided"`: the unrestricted test of Decision 3.
- `propDiffMin > 0`, `"greater"`: the numerator is restricted to the
  curve `thetaB - thetaA = propDiffMin`.
- `propDiffMin > 0`, `"twoSided"`: the e-process is the average of the
  two cumulative e-processes restricted at `+propDiffMin` and
  `-propDiffMin` (averaged as processes, not per block).
- `"less"`, or `"greater"` without `propDiffMin`, is not designed yet.

The restricted numerator builds on `learnPredictiveThetas` from `cond`:
`predictiveThetas2x2PropDiff(..., propDiff, nWeight = 1000)`, next to
`predictiveThetas2x2(...)` for the Beta posterior means of Decision 3;
`logEProcess2x2PlugIn(..., propDiff = NULL)` picks between them. The free
coordinate `rho` is `thetaA` rescaled to its feasible interval
`(max(0, -propDiff), min(1, 1 - propDiff))`, on `nWeight` equally spaced
grid points strictly inside `(0, 1)`, with the prior `Beta(betaA1, betaA2)`;
`betaB*` are not used. Block `i` uses the posterior mean of `thetaA` given blocks `1..i-1`,
and `thetaB = thetaA + propDiff`. The weights are updated on the log scale.
The denominator stays the pooled projection onto `thetaA = thetaB`, so the
null is always the point `thetaA = thetaB`; `"greater"` names the direction
of the alternative, not a composite null. The confidence sequence keeps the
unrestricted numerator. `nWeight` is not a design field.

### 7. Plug-in conditional e-factor on logOR

`savi2x2CondStat(ya, yb, na, nb, logOR, weightGrid = NULL, log = FALSE)`
returns the conditional e-factor of **one** block (scalar counts), not a
cumulative e-process. It conditions on the block's total `ya + yb`: under
the null `ya` is hypergeometric, under `logOR` (B minus A, anchored on
`thetaA`) it is Fisher's noncentral hypergeometric, and the log e-factor is
`ya * logOR - fnchLogPartition(na, nb, ya + yb, logOR) + lchoose(na + nb,
ya + yb)`, on the log scale when `log = TRUE`. `logOR` is one finite number
supplied by the caller (plug-in, e.g. GROW or UMP); `weightGrid` (a prior on
`logOR`) is reserved and errors for now.

`fnchLogPartition(na, nb, totalSuccesses, logOR)` is the log of
`sum_k choose(na, k) choose(nb, totalSuccesses - k) exp(logOR * k)` over the
feasible `k`, computed by a max-shifted log-sum-exp; at `logOR = 0` it
returns `lchoose(na + nb, totalSuccesses)` exactly.
