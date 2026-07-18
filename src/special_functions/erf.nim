{.experimental: "strictFuncs".}
{.push raises: [].}

import std/math except erf, erfc
import realmath
import mathutils

## Inverse error function (``inverseErf``) and scaled complementary error
## function (``erfcx``).  Both are generic over ``SomeFloat``.
##
## ``inverseErf`` uses a rational approximation followed by Newton
## refinement.  ``erfcx`` uses the product ``exp(x²)·erfc(x)`` for
## ``x < 6`` and an asymptotic series for ``x ≥ 6``.
##
## Returns ``NaN`` for out-of-domain inputs or when iterative algorithms
## fail to converge.

func inverseErf*[T: SomeFloat](x: T): T =
  ## Inverse error function :math:`\operatorname{erf}^{-1}(x)`.
  ##
  ## Uses a rational approximation followed by Newton refinement.
  ## Returns ``Inf`` for ``x = 1``, ``-Inf`` for ``x = -1``, and ``NaN``
  ## for ``|x| > 1`` or when the iteration fails to converge.
  if isNaN(x) or x > T(1.0) or x < T(-1.0):
    return nanT(T)

  if x == T(0.0):
    return T(0.0)

  if x == T(1.0):
    return infT(T)
  if x == T(-1.0):
    return -infT(T)

  let
    invErfA3 = T(-0.140543331)
    invErfA2 = T(0.914624893)
    invErfA1 = T(-1.645349621)
    invErfA0 = T(0.886226899)

    invErfB4 = T(0.012229801)
    invErfB3 = T(-0.329097515)
    invErfB2 = T(1.442710462)
    invErfB1 = T(-2.118377725)
    invErfB0 = T(1.0)

    invErfC3 = T(1.641345311)
    invErfC2 = T(3.429567803)
    invErfC1 = T(-1.62490649)
    invErfC0 = T(-1.970840454)

    invErfD2 = T(1.637067800)
    invErfD1 = T(3.543889200)
    invErfD0 = T(1.0)

  let xSign = if x > T(0.0): T(1.0) else: T(-1.0)
  let xAbs = abs(x)

  var r: T

  if xAbs <= T(0.7):
    let x2 = xAbs * xAbs
    r = xAbs * (((invErfA3 * x2 + invErfA2) * x2 + invErfA1) * x2 + invErfA0)
    r /= (((invErfB4 * x2 + invErfB3) * x2 + invErfB2) * x2 + invErfB1) * x2 + invErfB0
  else:
    let y = sqrt(-ln((T(1.0) - xAbs) / T(2.0)))
    r = (((invErfC3 * y + invErfC2) * y + invErfC1) * y + invErfC0)
    r /= ((invErfD2 * y + invErfD1) * y + invErfD0)

  r = r * xSign

  let twoOverSqrtPi = T(2.0) / sqrt(piT(T))

  var iter = 0
  var f = erf(r) - x
  while abs(f) > tolT(T) and iter < MAX_ITER:
    r -= f / (twoOverSqrtPi * exp(-r * r))
    if isNaN(r):
      return nanT(T)
    f = erf(r) - x
    iter += 1

  if iter >= MAX_ITER:
    return nanT(T)

  r

func erfcx*[T: SomeFloat](x: T): T =
  ## Scaled complementary error function
  ## :math:`\operatorname{erfcx}(x) = e^{x^2} \operatorname{erfc}(x)`.
  ##
  ## For ``x < 6`` computes the product directly; for ``x ≥ 6`` uses an
  ## asymptotic series.  For ``x`` so negative that ``exp(x²)`` overflows,
  ## returns ``Inf``.
  if isNaN(x):
    return nanT(T)

  if x < T(0.0):
    if x * x > maxLogT(T):
      return infT(T)
    return T(2.0) * exp(x * x) - erfcx(-x)

  if x < T(6.0):
    return exp(x * x) * erfc(x)

  let oneOverXSq = T(1.0) / (x * x)
  var term = T(1.0) / x
  var s = term
  var prevTerm = term
  var k = 1

  while abs(term) > tolT(T) * abs(s) and k < MAX_ITER:
    prevTerm = term
    term *= -oneOverXSq * (T(2.0) * T(k) - T(1.0)) / T(2.0)
    if abs(term) > abs(prevTerm):
      break
    s += term
    k += 1

  if k >= MAX_ITER:
    return nanT(T)

  s / sqrt(piT(T))

{.pop.}
