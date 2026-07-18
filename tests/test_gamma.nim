import std/math
import std/unittest

import special_functions

func approxEq(a, b: float): bool =
  abs(a - b) < 1e-13 * max(1.0, max(abs(a), abs(b)))

suite "Gamma":

  test "lowerIncompleteGamma(2.0, 0.25)":
    check approxEq(lowerIncompleteGamma(2.0, 0.25), 0.02649902116074391)

  test "lowerIncompleteGamma(2.0, 0.50)":
    check approxEq(lowerIncompleteGamma(2.0, 0.50), 0.09020401043104986)

  test "lowerIncompleteGamma(2.0, 0.75)":
    check approxEq(lowerIncompleteGamma(2.0, 0.75), 0.17335853270322426)

  test "lowerIncompleteGamma(2.0, 0.0) == 0":
    check lowerIncompleteGamma(2.0, 0.0) == 0.0

  test "upperIncompleteGamma(2.0, 0.25)":
    check approxEq(upperIncompleteGamma(2.0, 0.25), 0.9735009788392561)

  test "upperIncompleteGamma(2.0, 0.50)":
    check approxEq(upperIncompleteGamma(2.0, 0.50), 0.9097959895689501)

  test "upperIncompleteGamma(2.0, 0.75)":
    check approxEq(upperIncompleteGamma(2.0, 0.75), 0.8266414672967757)

  test "upperIncompleteGamma(2.0, 0.0) == gamma(2.0)":
    check approxEq(upperIncompleteGamma(2.0, 0.0), gamma(2.0))

  test "lower + upper incomplete gamma sums to full gamma":
    check approxEq(
      lowerIncompleteGamma(2.0, 0.5) + upperIncompleteGamma(2.0, 0.5),
      gamma(2.0))

  test "regularizedLowerIncompleteGamma(2.0, 0.5)":
    check approxEq(regularizedLowerIncompleteGamma(2.0, 0.5), 0.09020401043104986)

  test "regularizedUpperIncompleteGamma(2.0, 0.5)":
    check approxEq(regularizedUpperIncompleteGamma(2.0, 0.5), 0.9097959895689501)

  test "regularized lower + upper == 1":
    check approxEq(
      regularizedLowerIncompleteGamma(3.5, 2.0) +
      regularizedUpperIncompleteGamma(3.5, 2.0), 1.0)

  test "inverseLowerIncompleteGamma(2.0, 0.25)":
    check approxEq(inverseLowerIncompleteGamma(2.0, 0.25), 0.9612787631147771)

  test "inverseLowerIncompleteGamma(2.0, 0.50)":
    check approxEq(inverseLowerIncompleteGamma(2.0, 0.50), 1.6783469900166612)

  test "inverseLowerIncompleteGamma(2.0, 0.75)":
    check approxEq(inverseLowerIncompleteGamma(2.0, 0.75), 2.692634528889695)

  test "inverseLowerIncompleteGamma(2.0, 0.0) == 0":
    check inverseLowerIncompleteGamma(2.0, 0.0) == 0.0

  test "inverseUpperIncompleteGamma(2.0, 0.0) == Inf":
    check inverseUpperIncompleteGamma(2.0, 0.0) == Inf

  test "inverseUpperIncompleteGamma(2.0, gamma(2.0)) == 0":
    check approxEq(inverseUpperIncompleteGamma(2.0, gamma(2.0)), 0.0)

  test "inverseRegularizedLowerIncompleteGamma round-trip":
    check approxEq(
      inverseRegularizedLowerIncompleteGamma(2.0, 0.09020401043104986),
      0.5)

  test "inverseRegularizedUpperIncompleteGamma round-trip":
    check approxEq(
      inverseRegularizedUpperIncompleteGamma(2.0, 0.9097959895689501),
      0.5)

  test "inverse round-trip with a > 1":
    let p = regularizedLowerIncompleteGamma(3.5, 2.0)
    check approxEq(inverseRegularizedLowerIncompleteGamma(3.5, p), 2.0)

  test "inverse round-trip with a < 1":
    let p = regularizedLowerIncompleteGamma(0.5, 0.3)
    check approxEq(inverseRegularizedLowerIncompleteGamma(0.5, p), 0.3)

  test "inverseRegularizedLowerIncompleteGamma p close to 1 returns Inf":
    check inverseRegularizedLowerIncompleteGamma(2.0, 0.9999999999999999) == Inf

  test "precomputed G overload matches direct call":
    let G = gamma(2.0)
    check approxEq(lowerIncompleteGamma(2.0, 0.5, G),
                   lowerIncompleteGamma(2.0, 0.5))
    check approxEq(upperIncompleteGamma(2.0, 0.5, G),
                   upperIncompleteGamma(2.0, 0.5))
    check approxEq(inverseLowerIncompleteGamma(2.0, 0.5, G),
                   inverseLowerIncompleteGamma(2.0, 0.5))

  test "incompleteGamma a <= 0 returns NaN":
    check isNaN(lowerIncompleteGamma(0.0, 0.5))
    check isNaN(lowerIncompleteGamma(-1.0, 0.5))

  test "incompleteGamma x < 0 returns NaN":
    check isNaN(lowerIncompleteGamma(2.0, -0.1))

  test "incompleteGamma NaN input returns NaN":
    check isNaN(lowerIncompleteGamma(NaN, 0.5))
    check isNaN(lowerIncompleteGamma(2.0, NaN))

  test "inverseIncompleteGamma y < 0 returns NaN":
    check isNaN(inverseLowerIncompleteGamma(2.0, -0.1))

  test "inverseIncompleteGamma NaN input returns NaN":
    check isNaN(inverseLowerIncompleteGamma(NaN, 0.5))
    check isNaN(inverseLowerIncompleteGamma(2.0, NaN))

  test "inverseRegularizedUpper p ≈ 1 returns small value":
    let p = regularizedUpperIncompleteGamma(2.0, 0.01)
    check approxEq(inverseRegularizedUpperIncompleteGamma(2.0, p), 0.01)

  test "lgamma re-export":
    check approxEq(lgamma(5.0), 3.1780538303479458)

  test "regularizedLowerIncompleteGamma large a,x (overflow guard)":
    let r = regularizedLowerIncompleteGamma(170.0, 200.0)
    check not isNaN(r)
    check r >= 0.0 and r <= 1.0

  test "regularizedUpperIncompleteGamma round-trip a > 1":
    let q = regularizedUpperIncompleteGamma(3.5, 2.0)
    check approxEq(inverseRegularizedUpperIncompleteGamma(3.5, q), 2.0)

  test "regularizedUpperIncompleteGamma round-trip a < 1":
    let q = regularizedUpperIncompleteGamma(0.7, 0.5)
    check approxEq(inverseRegularizedUpperIncompleteGamma(0.7, q), 0.5)

  test "upperIncompleteGamma large x returns small value":
    check upperIncompleteGamma(2.0, 20.0) < 1e-7

  test "incompleteGamma at x=Inf":
    check regularizedLowerIncompleteGamma(2.0, Inf) == 1.0
    check regularizedUpperIncompleteGamma(2.0, Inf) == 0.0
    check lowerIncompleteGamma(2.0, Inf) == gamma(2.0)
    check upperIncompleteGamma(2.0, Inf) == 0.0

  test "inverseIncompleteGamma Inf input — p at boundary":
    check inverseRegularizedLowerIncompleteGamma(2.0, 1.0) == Inf
    check inverseRegularizedUpperIncompleteGamma(2.0, 0.0) == Inf

  test "gamma(Inf)":
    check gamma(Inf) == Inf

  test "lgamma(Inf)":
    check lgamma(Inf) == Inf
