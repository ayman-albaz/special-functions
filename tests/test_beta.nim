import std/math
import std/unittest

import special_functions

func approxEq(a, b: float): bool =
  abs(a - b) < 1e-13 * max(1.0, max(abs(a), abs(b)))

suite "Beta":

  test "beta(1.0, 2.0)":
    check approxEq(beta(1.0, 2.0), 0.5)

  test "beta(2.0, 1.0)":
    check approxEq(beta(2.0, 1.0), 0.5)

  test "beta(2.0, 2.0)":
    check approxEq(beta(2.0, 2.0), 0.16666666666666666)

  test "lbeta(2.0, 2.0) == ln(beta(2.0, 2.0))":
    check approxEq(lbeta(2.0, 2.0), ln(beta(2.0, 2.0)))

  test "lbeta(a, b) == lbeta(b, a)":
    check approxEq(lbeta(2.0, 5.0), lbeta(5.0, 2.0))

  test "lowerIncompleteBeta(2.0, 2.0, 0.25)":
    check approxEq(lowerIncompleteBeta(2.0, 2.0, 0.25), 0.02604166666666667)

  test "lowerIncompleteBeta(2.0, 2.0, 0.5)":
    check approxEq(lowerIncompleteBeta(2.0, 2.0, 0.5), 0.08333333333333333)

  test "lowerIncompleteBeta(2.0, 2.0, 0.75)":
    check approxEq(lowerIncompleteBeta(2.0, 2.0, 0.75), 0.140625)

  test "lowerIncompleteBeta(2.0, 2.0, 0.0) == 0":
    check approxEq(lowerIncompleteBeta(2.0, 2.0, 0.0), 0.0)

  test "lowerIncompleteBeta(2.0, 2.0, 1.0) == beta(2.0, 2.0)":
    check approxEq(lowerIncompleteBeta(2.0, 2.0, 1.0), beta(2.0, 2.0))

  test "upperIncompleteBeta(2.0, 2.0, 0.0) == beta(2.0, 2.0)":
    check approxEq(upperIncompleteBeta(2.0, 2.0, 0.0), beta(2.0, 2.0))

  test "upperIncompleteBeta(2.0, 2.0, 1.0) == 0":
    check approxEq(upperIncompleteBeta(2.0, 2.0, 1.0), 0.0)

  test "lower + upper incomplete beta sums to full beta":
    check approxEq(
      lowerIncompleteBeta(2.0, 2.0, 0.5) + upperIncompleteBeta(2.0, 2.0, 0.5),
      beta(2.0, 2.0))

  test "upperIncompleteBeta(2.0, 2.0, 0.5)":
    check approxEq(upperIncompleteBeta(2.0, 2.0, 0.5), 0.08333333333333333)

  test "regularizedLowerIncompleteBeta(2.0, 2.0, 0.25)":
    check approxEq(regularizedLowerIncompleteBeta(2.0, 2.0, 0.25), 0.15625)

  test "regularizedLowerIncompleteBeta(2.0, 2.0, 0.5)":
    check approxEq(regularizedLowerIncompleteBeta(2.0, 2.0, 0.5), 0.5)

  test "regularizedLowerIncompleteBeta(2.0, 2.0, 0.75)":
    check approxEq(regularizedLowerIncompleteBeta(2.0, 2.0, 0.75), 0.84375)

  test "regularizedUpperIncompleteBeta(2.0, 2.0, 0.5)":
    check approxEq(regularizedUpperIncompleteBeta(2.0, 2.0, 0.5), 0.5)

  test "regularized lower + upper == 1":
    check approxEq(
      regularizedLowerIncompleteBeta(3.0, 2.0, 0.4) +
      regularizedUpperIncompleteBeta(3.0, 2.0, 0.4), 1.0)

  test "inverseLowerIncompleteBeta(0.5, 0.5, 0.01)":
    check approxEq(inverseLowerIncompleteBeta(0.5, 0.5, 0.01), 2.4999791667e-5)

  test "inverseLowerIncompleteBeta(1.0, 1.0, 0.01)":
    check approxEq(inverseLowerIncompleteBeta(1.0, 1.0, 0.01), 0.01)

  test "inverseLowerIncompleteBeta(2.0, 2.0, 0.01)":
    check approxEq(inverseLowerIncompleteBeta(2.0, 2.0, 0.01), 0.149016952545326)

  test "inverseLowerIncompleteBeta(2.0, 2.0, 0.05)":
    check approxEq(inverseLowerIncompleteBeta(2.0, 2.0, 0.05), 0.363257491090568)

  test "inverseLowerIncompleteBeta(2.0, 2.0, 0.10)":
    check approxEq(inverseLowerIncompleteBeta(2.0, 2.0, 0.10), 0.567068922852268)

  test "inverseLowerIncompleteBeta(2.0, 2.0, 0.0) == 0":
    check inverseLowerIncompleteBeta(2.0, 2.0, 0.0) == 0.0

  test "inverseUpperIncompleteBeta(2.0, 2.0, 0.0) == Inf":
    check inverseUpperIncompleteBeta(2.0, 2.0, 0.0) == Inf

  test "inverseUpperIncompleteBeta(2.0, 2.0, beta(2,2)) == 0":
    check approxEq(inverseUpperIncompleteBeta(2.0, 2.0, beta(2.0, 2.0)), 0.0)

  test "inverseRegularizedLowerIncompleteBeta(2.0, 2.0, 0.15625) ≈ 0.25":
    check approxEq(inverseRegularizedLowerIncompleteBeta(2.0, 2.0, 0.15625), 0.25)

  test "inverseRegularizedUpperIncompleteBeta(2.0, 2.0, 0.5) ≈ 0.5":
    check approxEq(inverseRegularizedUpperIncompleteBeta(2.0, 2.0, 0.5), 0.5)

  test "inverse round-trip":
    let p = regularizedLowerIncompleteBeta(3.0, 1.5, 0.3)
    check approxEq(inverseRegularizedLowerIncompleteBeta(3.0, 1.5, p), 0.3)

  test "inverse round-trip with a < 1":
    let p = regularizedLowerIncompleteBeta(0.5, 0.8, 0.3)
    check approxEq(inverseRegularizedLowerIncompleteBeta(0.5, 0.8, p), 0.3)

  test "precomputed B overload matches direct call":
    let B = beta(2.0, 2.0)
    check approxEq(lowerIncompleteBeta(2.0, 2.0, 0.5, B),
                   lowerIncompleteBeta(2.0, 2.0, 0.5))
    check approxEq(upperIncompleteBeta(2.0, 2.0, 0.5, B),
                   upperIncompleteBeta(2.0, 2.0, 0.5))
    check approxEq(inverseLowerIncompleteBeta(2.0, 2.0, 0.01, B),
                   inverseLowerIncompleteBeta(2.0, 2.0, 0.01))

  test "incompleteBeta a <= 0 returns NaN":
    check isNaN(lowerIncompleteBeta(0.0, 2.0, 0.5))
    check isNaN(lowerIncompleteBeta(-1.0, 2.0, 0.5))

  test "incompleteBeta b <= 0 returns NaN":
    check isNaN(lowerIncompleteBeta(2.0, 0.0, 0.5))

  test "incompleteBeta x < 0 returns NaN":
    check isNaN(lowerIncompleteBeta(2.0, 2.0, -0.1))

  test "incompleteBeta x > 1 returns NaN":
    check isNaN(lowerIncompleteBeta(2.0, 2.0, 1.5))

  test "incompleteBeta NaN input returns NaN":
    check isNaN(lowerIncompleteBeta(NaN, 2.0, 0.5))
    check isNaN(lowerIncompleteBeta(2.0, NaN, 0.5))
    check isNaN(lowerIncompleteBeta(2.0, 2.0, NaN))

  test "inverseIncompleteBeta y < 0 returns NaN":
    check isNaN(inverseLowerIncompleteBeta(2.0, 2.0, -0.1))

  test "inverseIncompleteBeta NaN input returns NaN":
    check isNaN(inverseLowerIncompleteBeta(NaN, 2.0, 0.5))
    check isNaN(inverseLowerIncompleteBeta(2.0, NaN, 0.5))
    check isNaN(inverseLowerIncompleteBeta(2.0, 2.0, NaN))

  test "inverseIncompleteBeta y >= 1-eps returns Inf":
    check inverseLowerIncompleteBeta(2.0, 2.0, 2.0) == Inf

  test "inverseRegularized with a < 1, b < 1":
    check approxEq(inverseRegularizedLowerIncompleteBeta(0.5, 0.5, 0.3), 0.2061073738537633)

  test "regularizedUpperIncompleteBeta round-trip":
    let q = regularizedUpperIncompleteBeta(2.0, 2.0, 0.3)
    check approxEq(inverseRegularizedUpperIncompleteBeta(2.0, 2.0, q), 0.3)

  test "regularizedUpperIncompleteBeta round-trip a<1,b<1":
    let q = regularizedUpperIncompleteBeta(0.4, 0.7, 0.3)
    check approxEq(inverseRegularizedUpperIncompleteBeta(0.4, 0.7, q), 0.3)

  test "incompleteBeta at boundary x=0, x=1":
    check regularizedLowerIncompleteBeta(2.0, 2.0, 0.0) == 0.0
    check regularizedLowerIncompleteBeta(2.0, 2.0, 1.0) == 1.0
    check regularizedUpperIncompleteBeta(2.0, 2.0, 0.0) == 1.0
    check regularizedUpperIncompleteBeta(2.0, 2.0, 1.0) == 0.0
