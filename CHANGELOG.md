# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-07-18

### Changed
- **Breaking**: All compound function names converted to camelCase (e.g.
  `lowerIncompleteGamma`, `inverseLowerIncompleteGamma`, `inverseErf`).
  Math names (`beta`, `digamma`, `trigamma`, `erfcx`, `erf`, `erfc`, `lgamma`)
  retain their standard lowercase form.
- **Breaking**: All fallible functions now return `NaN` on out-of-domain inputs
  or when iterative algorithms fail to converge, instead of raising
  `ValueError` or silently producing garbage.
- **Breaking**: Require Nim >= 2.0.0.
- **Breaking**: Inverse incomplete gamma/beta functions return `Inf` instead of
  `NaN` when the inverse does not have a finite solution (e.g.
  `inverseLowerIncompleteGamma(a, gamma(a))` → `Inf`).
- Added `{.experimental: "strictFuncs".}` and `{.push raises: [].}` to all
  source modules.
- Added `## ...` docstrings on every exported proc and func.
- Domain checks use `a <= 0.0` / `x < 0.0` (replaced `a < EPS` / `x < EPS`),
  and all exported functions propagate `NaN` explicitly via `isNaN` checks.
- `EPS` constant is full float64 machine epsilon (`2.220446049250313e-16`).
- `mathutils.nim`: Lentz continued-fraction floor uses `LENTZ_TINY` (`1e-30`);
  dead `ctdfrac_eval` removed.
- `trigamma` asymptotic expansion extended to 7 Bernoulli terms (was 5),
  shift threshold raised to 16 (was 12).  `digamma` shift raised to 20.
- `erfcx` asymptotic series includes a divergence guard.
- Removed deprecated `nanChecks` / `infChecks` and `locks: 0` pragmas.
- Replaced `round`-based test comparisons with relative-tolerance `approxEq`.
- Performance: Newton loops for inverse functions no longer recompute
  `gamma(a)` / `beta(a, b)` on every iteration; Lentz loops inlined in
  `beta.nim` and `gamma.nim` (zero closure allocation).
- **`.float` → `float(...)`**: all member-call syntax replaced with
  function-call style for Nim 2.2+ compatibility.
- **Acklam inverse normal CDF**: replaced `as26_2_22` (A&S 26.2.22,
  ~4.5e-4 absolute error) with Acklam's rational approximation (~1.15e-9
  relative error), accelerating Newton convergence in inverses.
- **README accuracy**: softened to "~1e-13 where convergent" rather than
  "14 decimal places"; removed property-test TODO.

### Fixed
- **x=1 bug in `beta.nim`**: short-circuit for `abs(xn) < EPS` corrected.
- **Incomplete-gamma overflow**: `exp(a·ln(x) − x)` now normalized by
  `lgamma(a)` in the exponent, preventing `Inf/Inf → NaN` for large `a, x`
  (e.g. `a=500, x=600`).
- **NaN guards in Newton loops**: inverse gamma/beta now detect `isNaN(xn)`
  and enforce a halving cap (`≤ 20` consecutive halvings → `NaN`) to prevent
  `MAX_ITER` runaway.

### Added
- **Generic over `SomeFloat`**: all public functions are now generic,
  accepting both `float32` and `float64`.  Literals default to `float64`
  via Nim type inference; use `'f32` suffix for `float32`.
- **Native `float32` implementations**: `lgamma` (Lanczos approximation),
  `gamma`, `erf` (A&S 7.1.26 rational), and `erfc` all have pure-Nim
  `float32` paths, activated automatically when `T = float32`.  On
  `float64`, the library still delegates to the C standard library via
  `std/math`.
- **Type-tuned tolerances**: `EPS`, `TOL`, and `LENTZ_TINY` are now
  per-type generic helpers (`epsT`, `tolT`, `lentzTinyT`), ensuring
  Newton and Lentz continued-fraction loops converge properly for both
  precisions without oscillating at the lower float32 ULP.
- New internal module `special_functions/realmath.nim` housing the
  generic dispatch, native `float32` primitives, and per-type
  `infT`/`nanT`/`piT` helpers.
- **New functions**: `digamma`, `trigamma`, `erfcx`.
- **Re-exports**: `erf`, `erfc`, `lgamma` are now generic functions defined in
  `realmath.nim` rather than re-exports from `std/math`.  For `float64` callers,
  behavior is unchanged (they delegate to `std/math` internally).
- Precomputed `G`/`B` overloads for all incomplete and inverse functions,
  eliminating redundant `gamma(a)` / `beta(a,b)` calls in batch workloads.
- Explicit symbol-level exports in `special_functions.nim` with a module
  docstring.
- CI: `windows-latest` in OS matrix; `nim check` compilation step.
- Benchmarks: `benchmarks/bench.nim` (Nim) and `benchmarks/scipy_bench.py`
  (Python/SciPy).
- Tests (55 total): NaN propagation, `x=0`/`x=1`/`Inf` edge cases,
  precomputed-overload equivalence, round-trip invariants, `lbeta` symmetry,
  regularized-upper round-trips, digamma/trigamma large-x (`100`, `1000`) and
  small-x (`0.01`) accuracy, erfcx crossover (`5.9`, `6.0`, `6.1`, `15`, `25`),
  and gamma overflow regression.

### Removed
- All Julia references from README and benchmarks.
- `src/htmldocs/` (generated artifact).

## [0.2.0] - 2022-07-15

### Removed
- Modernized lib format

## [0.1.1] - 2021-06-06

### Removed
- Removed binomial submodule

## [0.1.0] - 2021-06-05

### Added
- Created special_functions library
