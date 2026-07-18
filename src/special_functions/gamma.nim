{.experimental: "strictFuncs".}
{.push raises: [].}

import std/math except lgamma, gamma, erf, erfc
import realmath
import mathutils

## Incomplete gamma function (lower/upper), regularized variants
## ``P(a,x)`` / ``Q(a,x)``, and their inverses.  Every function is generic
## over ``SomeFloat``.  Precomputed-``G`` overloads accept ``G = gamma(a)``
## to avoid recomputing the normalization constant in batch workloads.
##
## The forward incomplete gamma uses a continued fraction for ``x > a+1``
## and a power series otherwise, both computed in regularized form to
## prevent overflow for large ``a, x``.  The inverse uses Newton's method
## with a log-space step to guard against overflow.
##
## All fallible functions return ``NaN`` for out-of-domain inputs
## (``a ≤ 0``, ``x < 0``) or when iterative algorithms fail to converge.

func incompleteGammaInternal[T: SomeFloat](a, x, G: T, lower, reg: bool): T =
  if isNaN(a) or isNaN(x) or a <= T(0.0) or x < T(0.0):
    return nanT(T)

  if x == T(0.0):
    if lower:
      return T(0.0)
    else:
      return if reg: T(1.0) else: G

  if x == infT(T):
    if lower:
      return if reg: T(1.0) else: G
    else:
      return T(0.0)

  let logGinc = a * ln(x) - x - lgamma(a)
  var ginc = exp(logGinc)

  if x > (a + T(1.0)):
    var
      f = x - a + T(1.0)
      c = f
      d = T(0.0)
      delta = T(0.0)
      j = 1

    while abs(delta - T(1.0)) > tolT(T) and j < MAX_ITER:
      let fa = -T(j) * (T(j) - a)
      let fb = x - a + T(1.0) + T(2.0) * T(j)
      d = fb + fa * d
      if abs(d) < lentzTinyT(T):
        d = lentzTinyT(T)
      c = fb + fa / c
      if abs(c) < lentzTinyT(T):
        c = lentzTinyT(T)
      d = T(1.0) / d
      delta = c * d
      f *= delta
      j += 1

    if j >= MAX_ITER:
      return nanT(T)

    ginc /= f
    if lower:
      ginc = T(1.0) - ginc
  else:
    ginc /= a
    var term = ginc
    var at = a
    var i = 1
    while abs(term) > tolT(T) and i < MAX_ITER:
      at += T(1.0)
      term *= x / at
      ginc += term
      i += 1
    if i >= MAX_ITER:
      return nanT(T)
    if not lower:
      ginc = T(1.0) - ginc

  if reg:
    ginc
  else:
    ginc * G

func lowerIncompleteGamma*[T: SomeFloat](a, x: T): T =
  ## Lower incomplete gamma function
  ## :math:`\gamma(a,x) = \int_0^x t^{a-1} e^{-t} dt`.
  incompleteGammaInternal(a, x, gamma(a), lower = true, reg = false)

func lowerIncompleteGamma*[T: SomeFloat](a, x, G: T): T =
  ## Precomputed ``G = gamma(a)`` overload for batch workloads.
  incompleteGammaInternal(a, x, G, lower = true, reg = false)

func upperIncompleteGamma*[T: SomeFloat](a, x: T): T =
  ## Upper incomplete gamma function
  ## :math:`\Gamma(a,x) = \int_x^\infty t^{a-1} e^{-t} dt`.
  incompleteGammaInternal(a, x, gamma(a), lower = false, reg = false)

func upperIncompleteGamma*[T: SomeFloat](a, x, G: T): T =
  ## Precomputed ``G = gamma(a)`` overload for batch workloads.
  incompleteGammaInternal(a, x, G, lower = false, reg = false)

func regularizedLowerIncompleteGamma*[T: SomeFloat](a, x: T): T =
  ## Regularized lower incomplete gamma function
  ## :math:`P(a,x) = \gamma(a,x) / \Gamma(a)`.
  incompleteGammaInternal(a, x, gamma(a), lower = true, reg = true)

func regularizedLowerIncompleteGamma*[T: SomeFloat](a, x, G: T): T =
  ## Precomputed ``G = gamma(a)`` overload for batch workloads.
  incompleteGammaInternal(a, x, G, lower = true, reg = true)

func regularizedUpperIncompleteGamma*[T: SomeFloat](a, x: T): T =
  ## Regularized upper incomplete gamma function
  ## :math:`Q(a,x) = \Gamma(a,x) / \Gamma(a)`.
  incompleteGammaInternal(a, x, gamma(a), lower = false, reg = true)

func regularizedUpperIncompleteGamma*[T: SomeFloat](a, x, G: T): T =
  ## Precomputed ``G = gamma(a)`` overload for batch workloads.
  incompleteGammaInternal(a, x, G, lower = false, reg = true)

func inverseIncompleteGammaInternal[T: SomeFloat](a, y, G: T, lower, reg: bool): T =
  if isNaN(a) or isNaN(y) or a <= T(0.0) or y < T(0.0):
    return nanT(T)

  var p = y
  if not reg:
    p /= G
  if not lower:
    p = T(1.0) - p

  if abs(p) < epsT(T):
    return T(0.0)

  if p >= T(1.0) - epsT(T):
    return infT(T)

  let
    c1 = T(0.253)
    c2 = T(0.12)

  var x: T
  if a <= T(1.0):
    let pa = a * (c1 + c2 * a)
    let cpa = T(1.0) - pa
    if p < cpa:
      x = pow(p / cpa, T(1.0) / a)
    else:
      x = T(1.0) - ln((T(1.0) - p) / pa)
  else:
    x = -as26_2_22(p)
    x = T(1.0) - (T(1.0) / (T(9.0) * a)) + (x / (T(3.0) * sqrt(a)))
    x = a * x * x * x

  var i = 0
  var halvings = 0
  let logG = lgamma(a)
  var f = incompleteGammaInternal(a, x, G, lower = true, reg = true) - p

  while abs(f) > tolT(T) and i < MAX_ITER:
    let logFactor = logG + x - (a - T(1.0)) * ln(x)
    var step: T
    if logFactor > maxLogT(T):
      step = infT(T)
    else:
      step = f * exp(logFactor)
    if step == infT(T) or isNaN(step):
      x *= T(0.5)
      inc halvings
      if halvings > 20:
        return nanT(T)
    else:
      let xn = x - step
      if isNaN(xn) or xn < epsT(T):
        x *= T(0.5)
        inc halvings
        if halvings > 20:
          return nanT(T)
      else:
        x = xn
        halvings = 0
    f = incompleteGammaInternal(a, x, G, lower = true, reg = true) - p
    i += 1

  if i >= MAX_ITER:
    return nanT(T)

  x

func inverseLowerIncompleteGamma*[T: SomeFloat](a, y: T): T =
  ## Inverse of the lower incomplete gamma function.
  ## Returns ``x`` such that ``lowerIncompleteGamma(a, x) = y``.
  inverseIncompleteGammaInternal(a, y, gamma(a), lower = true, reg = false)

func inverseLowerIncompleteGamma*[T: SomeFloat](a, y, G: T): T =
  ## Precomputed ``G = gamma(a)`` overload for batch workloads.
  inverseIncompleteGammaInternal(a, y, G, lower = true, reg = false)

func inverseUpperIncompleteGamma*[T: SomeFloat](a, y: T): T =
  ## Inverse of the upper incomplete gamma function.
  ## Returns ``x`` such that ``upperIncompleteGamma(a, x) = y``.
  inverseIncompleteGammaInternal(a, y, gamma(a), lower = false, reg = false)

func inverseUpperIncompleteGamma*[T: SomeFloat](a, y, G: T): T =
  ## Precomputed ``G = gamma(a)`` overload for batch workloads.
  inverseIncompleteGammaInternal(a, y, G, lower = false, reg = false)

func inverseRegularizedLowerIncompleteGamma*[T: SomeFloat](a, y: T): T =
  ## Inverse of the regularized lower incomplete gamma function.
  ## Returns ``x`` such that ``regularizedLowerIncompleteGamma(a, x) = y``.
  inverseIncompleteGammaInternal(a, y, gamma(a), lower = true, reg = true)

func inverseRegularizedLowerIncompleteGamma*[T: SomeFloat](a, y, G: T): T =
  ## Precomputed ``G = gamma(a)`` overload for batch workloads.
  inverseIncompleteGammaInternal(a, y, G, lower = true, reg = true)

func inverseRegularizedUpperIncompleteGamma*[T: SomeFloat](a, y: T): T =
  ## Inverse of the regularized upper incomplete gamma function.
  ## Returns ``x`` such that ``regularizedUpperIncompleteGamma(a, x) = y``.
  inverseIncompleteGammaInternal(a, y, gamma(a), lower = false, reg = true)

func inverseRegularizedUpperIncompleteGamma*[T: SomeFloat](a, y, G: T): T =
  ## Precomputed ``G = gamma(a)`` overload for batch workloads.
  inverseIncompleteGammaInternal(a, y, G, lower = false, reg = true)

{.pop.}
