![Linux Build Status (Github Actions)](https://github.com/ayman-albaz/special-functions/actions/workflows/install_and_test.yml/badge.svg) [![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# Special Functions

Special functions is a Nim library for calculating special mathematical functions
that are not found in Nim's standard `math` module.

All functions are **generic over ``SomeFloat``**, supporting both ``float32`` and
``float64``.  On ``float64``, ``lgamma``, ``gamma``, ``erf``, and ``erfc``
delegate directly to the C standard library.  On ``float32``, the library
provides native Nim implementations (Lanczos approximation for ``lgamma``,
A&S 7.1.26 for ``erf``).

Functions return ``NaN`` when given out-of-domain inputs or when iterative
algorithms fail to converge. Use ``isNaN`` from ``math`` to detect failures.

## Requirements

- Nim >= 2.0.0

## Supported Functions

All functions are generic over ``SomeFloat``.  Literals default to
``float64``; use ``'f32`` suffix for ``float32``.

### Gamma Functions

```nim
import special_functions

let x1 = lowerIncompleteGamma(2.0, 0.25)         # γ(a, x)  — float64
let x1f = lowerIncompleteGamma(2.0'f32, 0.25'f32)  # γ(a, x)  — float32
let x2 = upperIncompleteGamma(2.0, 0.25)         # Γ(a, x)
let x3 = regularizedLowerIncompleteGamma(2.0, 0.25)  # P(a, x)
let x4 = regularizedUpperIncompleteGamma(2.0, 0.25)  # Q(a, x)

let x5 = inverseLowerIncompleteGamma(2.0, 0.25)
let x6 = inverseUpperIncompleteGamma(2.0, 0.25)
let x7 = inverseRegularizedLowerIncompleteGamma(2.0, 0.25)
let x8 = inverseRegularizedUpperIncompleteGamma(2.0, 0.25)
```

For batch operations with the same shape parameter, use the precomputed ``G``
overloads to avoid recomputing ``gamma(a)``:

```nim
let G = gamma(3.5)
for x in values:
  let y = lowerIncompleteGamma(3.5, x, G)
```

### Beta Functions

```nim
import special_functions

let b  = beta(2.0, 2.0)          # B(a, b)
let lb = lbeta(2.0, 2.0)        # ln B(a, b)

let y1 = lowerIncompleteBeta(2.0, 2.0, 0.25)
let y2 = upperIncompleteBeta(2.0, 2.0, 0.25)
let y3 = regularizedLowerIncompleteBeta(2.0, 2.0, 0.25)
let y4 = regularizedUpperIncompleteBeta(2.0, 2.0, 0.25)

let y5 = inverseLowerIncompleteBeta(2.0, 2.0, 0.25)
let y6 = inverseUpperIncompleteBeta(2.0, 2.0, 0.25)
let y7 = inverseRegularizedLowerIncompleteBeta(2.0, 2.0, 0.25)
let y8 = inverseRegularizedUpperIncompleteBeta(2.0, 2.0, 0.25)
```

Precomputed ``B`` overloads are also available:

```nim
let B = beta(3.0, 1.5)
for x in values:
  let y = lowerIncompleteBeta(3.0, 1.5, x, B)
```

### Error Functions

```nim
import special_functions

let e   = erf(1.0)          # generic — float64 → C stdlib
let ef  = erf(1.0'f32)       # generic — float32 → native A&S 7.1.26
let ec  = erfc(1.0)         # generic
let ecx = erfcx(1.0)        # exp(x²)·erfc(x)
let z   = inverseErf(0.25)
```

### Polygamma Functions

```nim
import special_functions

let d = digamma(1.0)        # ψ(x) = Γ'(x)/Γ(x)
let t = trigamma(1.0)       # ψ₁(x) = ψ'(x)
```

### Log Gamma / Gamma

```nim
import special_functions

let lg = lgamma(2.0)        # generic — float64 → C stdlib
let lg32 = lgamma(2.0'f32)   # generic — float32 → native Lanczos
let g  = gamma(5.0)         # generic
```

## Accuracy

All functions achieve ~1e-13 relative accuracy where convergent (float64).
On float32, accuracy is ~1e-5 relative (approximately 1 ULP for most
inputs).  The asymptotic branches (``erfcx`` for ``x >= 10``, inverse
Newton iterations) may reach ~1e-12 in adverse regimes; consult the test
suite for exact tolerances per function.  Reference values are validated
against SciPy (Python).

## Performance

Representative timings for 100,000 iterations each (--d:release, Nim 2.x;
SciPy 1.9, Python 3.11).  Lower is better.

### Nim

```
lowerIncompleteGamma(2, 0.5):          5.6 ms
regularizedLowerIncompleteGamma(2, 0.5):  5.6 ms
inverseLowerIncompleteGamma(2, 0.5):     42.2 ms
lgamma(5):                              1.8 ms
beta(2, 2):                             3.5 ms
regularizedLowerIncompleteBeta(2, 2, 0.5): 15.0 ms
inverseRegLowerIncompleteBeta(2, 2, 0.5):  33.8 ms
erf(1):                               < 0.1 ms
erfcx(3):                               2.7 ms
inverseErf(0.5):                        5.5 ms
digamma(1):                             3.6 ms
trigamma(1):                            2.1 ms
```

Precomputed `G`/`B` overloads cut cost for batch use:

```
lowerIncompleteGamma(2, 0.5, G):       3.5 ms  (40% faster)
lowerIncompleteBeta(2, 2, 0.5, B):    11.3 ms  (25% faster)
```

### SciPy (Python 3.11)

```
gammainc(0.5, 2)           [lower reg]:  91.9 ms
gammaincc(0.5, 2)          [upper reg]:  91.5 ms
gammaincinv(0.5, 2)       [lower inv reg]: 849 ms
gamma(2)                   [scalar gamma]:  55.1 ms
gammaln(5)                [log gamma]:      50.2 ms
beta(2, 2)                [beta]:           72.0 ms
betainc(2, 2, 0.5)        [lower reg]:     104 ms
betaincinv(2, 2, 0.5)     [lower inv reg]: 105 ms
erf(1):                    52.2 ms
erfcx(3):                  52.8 ms
erfinv(0.5)               [inverse erf]:    54.1 ms
digamma(1):                59.5 ms
polygamma(1, 1)           [trigamma]:      979 ms
```

To reproduce:

```bash
nim r -d:release --path:src benchmarks/bench.nim
python3 benchmarks/scipy_bench.py
```

## Returning NaN

All fallible functions in this library return `NaN` when given out-of-domain
inputs or when iterative algorithms fail to converge. Examples:

- `digamma(0)` — pole at zero → `NaN`
- `inverseErf(2.0)` — out of domain → `NaN`
- Non-convergence of an iterative solver → `NaN`

Use `isNaN` from `math` to detect failures at the call site.

## TODO

- Add more special functions on an as-needed basis.

Performance, feature, and documentation PRs are always welcome.

## Acknowledgments

Many of the algorithms used here are based off of different peoples' code.
Thank you to:

- https://github.com/jkovacic/math/
- http://libit.sourceforge.net/math_8c-source.html

## Contact

I can be reached at aymanalbaz98@gmail.com
