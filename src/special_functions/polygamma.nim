{.experimental: "strictFuncs".}
{.push raises: [].}

import std/math
import realmath

## Digamma ``ψ(x)`` and trigamma ``ψ₁(x)`` functions, generic over
## ``SomeFloat``.  Both use recurrence to shift ``x`` above a threshold
## (20 for digamma, 16 for trigamma) followed by an asymptotic expansion
## with 7 Bernoulli terms.
##
## Returns ``NaN`` for ``x = 0, -1, -2, …`` (poles).

func digamma*[T: SomeFloat](x: T): T =
  ## Digamma function :math:`\psi(x) = \Gamma'(x) / \Gamma(x)`.
  ##
  ## Uses recurrence to shift ``x`` above 20 followed by an asymptotic
  ## expansion with 7 Bernoulli terms.  Returns ``NaN`` for ``x = 0, -1, -2, …``.
  if isNaN(x):
    return nanT(T)

  if x <= T(0.0) and floor(x) == x:
    return nanT(T)

  var s = x
  var w = T(0.0)

  while s < T(20.0):
    w += T(1.0) / s
    s += T(1.0)

  let r = T(1.0) / s
  var z = r * r
  var y = ln(s) - T(0.5) * r

  z = -z / T(12.0) +
      z * z / T(120.0) -
      z * z * z / T(252.0) +
      z * z * z * z / T(240.0) -
      z * z * z * z * z / T(132.0) +
      z * z * z * z * z * z * T(691.0) / T(32760.0) -
      z * z * z * z * z * z * z / T(12.0)

  y - w + z

func trigamma*[T: SomeFloat](x: T): T =
  ## Trigamma function :math:`\psi_1(x) = \psi'(x)`.
  ##
  ## Uses recurrence to shift ``x`` above 16 followed by an asymptotic
  ## expansion with 7 Bernoulli terms.  Returns ``NaN`` for ``x = 0, -1, -2, …``.
  if isNaN(x):
    return nanT(T)

  if x <= T(0.0) and floor(x) == x:
    return nanT(T)

  var s = x
  var w = T(0.0)

  while s < T(16.0):
    w += T(1.0) / (s * s)
    s += T(1.0)

  let r = T(1.0) / s
  let w1 = r * r
  let r3 = r * w1

  let z = T(1.0) / T(6.0) +
    w1 * (-T(1.0) / T(30.0) +
    w1 * (T(1.0) / T(42.0) +
    w1 * (-T(1.0) / T(30.0) +
    w1 * (T(5.0) / T(66.0) +
    w1 * (-T(691.0) / T(2730.0) +
    w1 * (T(7.0) / T(6.0)))))))

  let y = r + T(0.5) * w1 + r3 * z

  w + y

{.pop.}
