# Repository agent instructions

## Current objective

Work only on the safe anytime-valid 2x2 test and its direct supporting code.
The priorities, in order, are:

1. statistical validity;
2. explicit assumptions and parameter orientation;
3. numerically stable, direct implementations;
4. focused tests of statistical invariants;
5. clean, maintainable R code and accurate documentation.

This is an active port and implementation effort. Preserve behaviour only when
that behaviour is known to be intentional and statistically correct. Do not
retain a questionable result merely because it is present in older code.

## Scope boundary

The normal working scope is:

- `R/safe2x2Test.R`;
- `R/safe2x2TestCond.R`;
- 2x2-specific tests under `tests/testthat/`;
- 2x2-specific roxygen comments and generated `man/` pages;
- `vignettes/savi2x2.Rmd` and
  `vignettes/contingency-tables-vignette.Rmd` when the task explicitly requires
  user-facing documentation;
- the smallest necessary 2x2-related changes to `NAMESPACE`, `DESCRIPTION`,
  `R/safeS3Methods.R`, or compatibility wrappers in `R/deprecate.R`.

Everything else is out of scope. In particular, do not modify the t-test,
z-test, log-rank test, general design framework, unrelated S3 behaviour,
package-wide style, or unrelated documentation. If a 2x2 change appears to
require broader work, stop and explain the dependency before editing it.

Treat `scratch/` as exploratory material, not as a source of truth or a package
dependency. Do not promote scratch code without checking its derivation,
licence/provenance, numerical behaviour, and fit with the package API.

## Model roles

These are responsibility and routing rules, not a requirement to spawn three
agents for every task. If the active model is known, follow its role. If model
selection or delegation is available, route work as follows.

### GPT-5.6 Sol — statistical lead and final reviewer

Use Sol for work where a subtle error could invalidate inference:

- define the null, alternative, conditioning argument, and parameter direction;
- assess whether a quantity is an e-variable, an e-factor, or a cumulative
  e-process;
- review optional-stopping and optional-continuation claims;
- derive or review likelihood ratios, mixtures, NML normalisers, confidence-set
  inversion, and design criteria;
- resolve discrepancies between the R port, a paper, and a reference
  implementation;
- approve changes to statistical formulas, public semantics, or exported APIs;
- perform the final statistical review of consequential changes.

Sol must state the assumptions and a checkable correctness criterion before a
statistically consequential edit. A simulation may diagnose behaviour, but it
does not by itself prove validity.

### GPT-5.6 Terra — implementation and test maintainer

Use Terra for the main engineering work once the statistical contract is clear:

- implement focused R changes and ports;
- translate defining formulas directly into readable code;
- add high-value unit tests and deterministic numerical regressions;
- improve numerical stability without changing statistical meaning;
- maintain roxygen, package exports, and 2x2-facing documentation;
- run focused tests and package checks, then report exact results.

Terra must escalate to Sol when the formula, parameter orientation, cumulative
semantics, edge-case convention, or expected result is ambiguous. It must not
silently choose among statistically different interpretations.

### GPT-5.6 Luna — bounded mechanical support

Use Luna for low-risk, precisely specified work:

- inventory 2x2 symbols, call sites, tests, and documentation;
- compare signatures, names, and fixed reference outputs;
- locate stale generated documentation or exports;
- make narrow spelling, formatting, and roxygen consistency edits;
- run prescribed checks and summarise failures without diagnosing new theory.

Luna must not decide statistical formulas, infer missing method definitions,
change public semantics, loosen tolerances to make tests pass, or declare a
method statistically valid. Route those decisions to Sol; route non-trivial
implementation to Terra.

## Required workflow

Before editing:

1. Read `git status` and preserve all unrelated work in the dirty worktree.
2. Read the complete functions, tests, and roxygen blocks involved.
3. Identify the source being ported: paper/equation, earlier implementation,
   issue, or agreed specification. Do not guess when sources disagree.
4. State the statistical contract: data layout, conditioning, null and
   alternative, direction of each effect parameter, and whether outputs are
   blockwise or cumulative.
5. Define the narrowest observable success criterion.

During implementation:

- Make the smallest coherent change that satisfies the task.
- Keep formulas visibly comparable with their definitions.
- Separate statistical logic from simulation, printing, and plotting.
- Preserve public APIs unless the user explicitly authorises a change.
- Do not add dependencies or abstractions without a demonstrated need.
- Do not rewrite, rename, reformat, stage, or commit unrelated code.
- Do not hand-edit generated `.Rd` files. Change roxygen source and regenerate
  documentation when the toolchain is available.

After implementation:

1. Inspect the diff for accidental out-of-scope edits.
2. Run the narrowest relevant deterministic tests.
3. Run broader package checks in proportion to the change.
4. Report what was verified, what was not, and any remaining statistical or
   numerical uncertainty.

## Statistical correctness rules

- Keep the meanings of `ya`, `yb`, `na`, and `nb` explicit. Success counts must
  be integer-valued, non-negative, and no larger than their group sizes.
- Define effect orientation at every API boundary. In particular, do not assume
  that a probability difference and a conditional log-odds ratio use the same
  A/B direction merely because their signs look similar.
- Distinguish a blockwise e-factor from a cumulative e-process. Name, document,
  combine, and test them according to their actual contract.
- For sequential methods, verify that step `i` uses exactly the information
  available through step `i`; prevent look-ahead through full-sample estimates.
- Optional-stopping validity requires the appropriate null expectation or
  conditional expectation property. Rejection-rate simulations are supporting
  diagnostics, not a substitute for that argument.
- Compute products, likelihood ratios, mixtures, and normalising constants on
  the log scale when underflow or overflow is plausible. Handle zero-probability
  cases deliberately rather than allowing `NaN`, `Inf - Inf`, or silent zeros to
  choose the result.
- For conditional 2x2 methods, verify the support of the (noncentral)
  hypergeometric distribution and behaviour at degenerate margins.
- For one-sided methods, verify the sign convention and group-swap relationship.
  Do not obtain a two-sided method by an ad hoc maximum, sum, or factor of two
  without a documented validity argument.
- Confidence intervals or sequences must be obtained by valid inversion of the
  intended test/e-process. Label placeholders as unimplemented; never present a
  trivial bound as a computed interval.
- Simulations must expose or respect reproducibility via the caller's RNG state,
  retain non-crossing paths explicitly, and report Monte Carlo uncertainty when
  used for design decisions.

## Testing expectations

For a legacy 2x2 migration, begin with an end-to-end comparison in
`scratch/debug-2x2.R`: load the legacy code in a separate environment, use
explicit data and prior parameters, print old and new outputs side by side, and
fail on a justified numerical tolerance. Only add narrower regression tests
after this comparison establishes the intended contract.

Prefer small, exact, high-information tests. As applicable, cover:

- a hand-computable 2x2 table or comparison with an independent implementation;
- one block versus the corresponding first step of a multi-block call;
- cumulative output versus the product or sum of its documented increments;
- `log = TRUE` versus the logarithm of ordinary-scale output;
- scalar block-size recycling versus explicit vectors;
- swapping groups together with the appropriate alternative transformation;
- boundary tables, degenerate margins, and impossible observations;
- seeded simulation output, non-crossing paths, and stopping-time threshold
  conventions;
- regression values for a previously identified statistical or numerical bug.

Do not create broad Cartesian test matrices for routine validation. Never relax
a numerical tolerance until the reason for the discrepancy is understood.

For focused 2x2 tests, prefer:

```sh
Rscript -e 'testthat::test_local(".", filter = "2x2")'
```

When roxygen or exports change, regenerate documentation and inspect the
resulting `NAMESPACE` and relevant `.Rd` diff. For a completed package-facing
change, also run an appropriate package build/check, normally:

```sh
R CMD build .
R CMD check --no-manual safestats_*.tar.gz
```

If a command is unavailable, too expensive, or blocked by unrelated failures,
say so explicitly and run the strongest narrower check available.

## R style

- Use 2-space indentation, `<-` for assignment, and spaces around operators.
- Use `camelCase` for functions, arguments, and local variables, matching the
  established package API.
- Prefer guard clauses and a flat happy path over deeply nested branches.
- Use `match.arg()` for enumerated character arguments.
- Validate inputs that could otherwise yield a plausible but statistically
  misleading answer; avoid defensive machinery that obscures the method.
- Keep comments concise and explain statistical intent or a non-obvious
  numerical choice, not syntax already clear from the code.
- Use roxygen2 Markdown for exported functions and document return semantics,
  effect direction, and sequential accumulation precisely.
- Use structured section comments such as `# Section ----` and
  `### Action: Target ----` where they genuinely improve navigation.
- Use `TODO(Name):` for a concrete pending task with an identifiable owner.

## Completion standard

A 2x2 task is complete only when the implementation, tests, documentation, and
exports agree on the same statistical contract; the relevant checks pass; and
the final report identifies the source or invariant used to establish
correctness. Clean code without a sound statistical argument is not complete,
and statistically plausible code without focused verification is not complete.
