{.experimental: "strictFuncs".}
{.push raises: [].}

import std/math except lgamma, gamma, erf, erfc
import realmath
import mathutils

## Beta function, its logarithm, and all incomplete / inverse / regularized
## variants.  Every function is generic over ``SomeFloat``.  Precomputed-``B``
## overloads accept ``B = beta(a, b)`` to avoid recomputing the normalization
## constant in batch workloads.
##
## All fallible functions return ``NaN`` for out-of-domain inputs
## (``a ≤ 0``, ``b ≤ 0``, ``x ∉ [0, 1]``) or when iterative algorithms
## fail to converge.

func lbeta*[T: SomeFloat](a, b: T): T =
  ## Natural log of the beta function :math:`\ln B(a,b)`.
  lgamma(a) + lgamma(b) - lgamma(a + b)

func beta*[T: SomeFloat](a, b: T): T =
  ## Beta function :math:`B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}`.
  exp(lbeta(a, b))

func incompleteBetaInternal[T: SomeFloat](a, b, x, B: T, lower, reg: bool): T =
  if isNaN(a) or isNaN(b) or isNaN(x) or a <= T(0.0) or b <= T(0.0) or x < T(0.0) or x > T(1.0):
    return nanT(T)

  let xlarge = x > ((a + T(1.0)) / (a + b + T(2.0)))
  let xn = if xlarge: T(1.0) - x else: x
  let an = if xlarge: b else: a
  let bn = if xlarge: a else: b

  var binc: T

  if abs(xn) < epsT(T):
    binc = T(0.0)
  else:
    binc = pow(xn, an) * pow(T(1.0) - xn, bn) / an

    var
      f = T(1.0)
      c = f
      d = T(0.0)
      delta = T(0.0)
      j = 1

    while abs(delta - T(1.0)) > tolT(T) and j < MAX_ITER:
      let m: T = T(j div 2)
      let fa = if floorMod(j, 2) == 1:
                 -((an + m) * (an + bn + m) * xn) / ((an + T(j) - T(1.0)) * (an + T(j)))
               else:
                 (m * (bn - m) * xn) / ((an + T(j) - T(1.0)) * (an + T(j)))
      let fb = T(1.0)
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

    binc /= f

  let complement = (xlarge == lower)
  if complement:
    binc = B - binc

  if reg:
    binc /= B

  binc

func lowerIncompleteBeta*[T: SomeFloat](a, b, x: T): T =
  ## Lower incomplete beta function
  ## :math:`B_x(a,b) = \int_0^x t^{a-1}(1-t)^{b-1} dt`.
  incompleteBetaInternal(a, b, x, beta(a, b), lower = true, reg = false)

func lowerIncompleteBeta*[T: SomeFloat](a, b, x, B: T): T =
  ## Precomputed ``B = beta(a, b)`` overload for batch workloads.
  incompleteBetaInternal(a, b, x, B, lower = true, reg = false)

func upperIncompleteBeta*[T: SomeFloat](a, b, x: T): T =
  ## Upper incomplete beta function :math:`B(a,b) - B_x(a,b)`.
  incompleteBetaInternal(a, b, x, beta(a, b), lower = false, reg = false)

func upperIncompleteBeta*[T: SomeFloat](a, b, x, B: T): T =
  ## Precomputed ``B = beta(a, b)`` overload for batch workloads.
  incompleteBetaInternal(a, b, x, B, lower = false, reg = false)

func regularizedLowerIncompleteBeta*[T: SomeFloat](a, b, x: T): T =
  ## Regularized lower incomplete beta function
  ## :math:`I_x(a,b) = B_x(a,b) / B(a,b)`.
  incompleteBetaInternal(a, b, x, beta(a, b), lower = true, reg = true)

func regularizedLowerIncompleteBeta*[T: SomeFloat](a, b, x, B: T): T =
  ## Precomputed ``B = beta(a, b)`` overload for batch workloads.
  incompleteBetaInternal(a, b, x, B, lower = true, reg = true)

func regularizedUpperIncompleteBeta*[T: SomeFloat](a, b, x: T): T =
  ## Regularized upper incomplete beta function
  ## :math:`1 - I_x(a,b)`.
  incompleteBetaInternal(a, b, x, beta(a, b), lower = false, reg = true)

func regularizedUpperIncompleteBeta*[T: SomeFloat](a, b, x, B: T): T =
  ## Precomputed ``B = beta(a, b)`` overload for batch workloads.
  incompleteBetaInternal(a, b, x, B, lower = false, reg = true)

func inverseIncompleteBetaInternal[T: SomeFloat](a, b, y, B: T, lower, reg: bool): T =
  if isNaN(a) or isNaN(b) or isNaN(y) or a <= T(0.0) or b <= T(0.0) or y < T(0.0):
    return nanT(T)

  var p = y
  if not reg:
    p /= B
  if not lower:
    p = T(1.0) - p

  if abs(p) < epsT(T):
    return T(0.0)
  if p >= T(1.0) - epsT(T):
    return infT(T)

  var x: T
  if a >= T(1.0) and b >= T(1.0):
    x = as26_2_22(p)
    let lambda = (x * x - T(3.0)) / T(6.0)
    let ta = T(1.0) / (T(2.0) * a - T(1.0))
    let tb = T(1.0) / (T(2.0) * b - T(1.0))
    let h = T(2.0) / (ta + tb)
    let w = x * sqrt(h + lambda) / h - (tb - ta) * (lambda + T(5.0) / T(6.0) - T(2.0) / (T(3.0) * h))
    x = a / (a + b * exp(T(2.0) * w))
  else:
    let ta = pow(a / (a + b), a) / a
    let tb = pow(b / (a + b), b) / b
    let S = ta + tb
    if p < ta / S:
      x = pow(p * S * a, T(1.0) / a)
    else:
      x = T(1.0) - pow((T(1.0) - p) * S * b, T(1.0) / b)

  var i = 0
  var halvings = 0
  let logB = lbeta(a, b)
  var f = incompleteBetaInternal(a, b, x, B, lower = true, reg = true) - p
  while abs(f) > tolT(T) and i < MAX_ITER:
    var logFactor: T
    if x <= epsT(T) or x >= T(1.0) - epsT(T):
      logFactor = infT(T)
    else:
      logFactor = logB - (a - T(1.0)) * ln(x) - (b - T(1.0)) * ln(T(1.0) - x)
    var step: T
    if logFactor > maxLogT(T):
      step = infT(T)
    else:
      step = f * exp(logFactor)
    if step == infT(T) or isNaN(step):
      x *= T(0.5)
      inc halvings
    elif isNaN(x - step) or x - step < epsT(T):
      x *= T(0.5)
      inc halvings
    elif x - step > T(1.0) - epsT(T):
      x = (T(1.0) + x) * T(0.5)
      inc halvings
    else:
      x = x - step
      halvings = 0
    if halvings > 20:
      return nanT(T)
    f = incompleteBetaInternal(a, b, x, B, lower = true, reg = true) - p
    i += 1

  if i >= MAX_ITER:
    return nanT(T)

  x

func inverseLowerIncompleteBeta*[T: SomeFloat](a, b, y: T): T =
  ## Inverse of the lower incomplete beta function.
  ## Returns ``x`` such that ``lowerIncompleteBeta(a, b, x) = y``.
  inverseIncompleteBetaInternal(a, b, y, beta(a, b), lower = true, reg = false)

func inverseLowerIncompleteBeta*[T: SomeFloat](a, b, y, B: T): T =
  ## Precomputed ``B = beta(a, b)`` overload for batch workloads.
  inverseIncompleteBetaInternal(a, b, y, B, lower = true, reg = false)

func inverseUpperIncompleteBeta*[T: SomeFloat](a, b, y: T): T =
  ## Inverse of the upper incomplete beta function.
  ## Returns ``x`` such that ``upperIncompleteBeta(a, b, x) = y``.
  inverseIncompleteBetaInternal(a, b, y, beta(a, b), lower = false, reg = false)

func inverseUpperIncompleteBeta*[T: SomeFloat](a, b, y, B: T): T =
  ## Precomputed ``B = beta(a, b)`` overload for batch workloads.
  inverseIncompleteBetaInternal(a, b, y, B, lower = false, reg = false)

func inverseRegularizedLowerIncompleteBeta*[T: SomeFloat](a, b, y: T): T =
  ## Inverse of the regularized lower incomplete beta function.
  ## Returns ``x`` such that ``regularizedLowerIncompleteBeta(a, b, x) = y``.
  inverseIncompleteBetaInternal(a, b, y, beta(a, b), lower = true, reg = true)

func inverseRegularizedLowerIncompleteBeta*[T: SomeFloat](a, b, y, B: T): T =
  ## Precomputed ``B = beta(a, b)`` overload for batch workloads.
  inverseIncompleteBetaInternal(a, b, y, B, lower = true, reg = true)

func inverseRegularizedUpperIncompleteBeta*[T: SomeFloat](a, b, y: T): T =
  ## Inverse of the regularized upper incomplete beta function.
  ## Returns ``x`` such that ``regularizedUpperIncompleteBeta(a, b, x) = y``.
  inverseIncompleteBetaInternal(a, b, y, beta(a, b), lower = false, reg = true)

func inverseRegularizedUpperIncompleteBeta*[T: SomeFloat](a, b, y, B: T): T =
  ## Precomputed ``B = beta(a, b)`` overload for batch workloads.
  inverseIncompleteBetaInternal(a, b, y, B, lower = false, reg = true)

{.pop.}
