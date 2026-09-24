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
