import std/math
import std/unittest
import special_functions

func approxEq32(a, b: float32): bool =
  abs(a - b) < 1e-5'f32 * max(1.0'f32, max(abs(a), abs(b)))

suite "Generic functions with float32":

  test "lgamma float32 via generic":
    check approxEq32(lgamma(1.0'f32), float32(lgamma(1.0)))
    check approxEq32(lgamma(5.0'f32), float32(lgamma(5.0)))
    check approxEq32(lgamma(10.0'f32), float32(lgamma(10.0)))

  test "gamma float32 via generic":
    check approxEq32(gamma(1.0'f32), float32(gamma(1.0)))
    check approxEq32(gamma(5.0'f32), float32(gamma(5.0)))

  test "erf float32 via generic":
    check approxEq32(erf(0.5'f32), float32(erf(0.5)))
    check approxEq32(erf(1.0'f32), float32(erf(1.0)))
    check approxEq32(erf(-0.5'f32), float32(erf(-0.5)))

  test "erfc float32 via generic":
    check approxEq32(erfc(0.5'f32), float32(erfc(0.5)))
    check approxEq32(erfc(1.0'f32), float32(erfc(1.0)))
    check approxEq32(erf(0.5'f32) + erfc(0.5'f32), 1.0'f32)

  test "inverseErf float32":
    check approxEq32(inverseErf(0.5'f32), float32(inverseErf(0.5)))
    check inverseErf(0.0'f32) == 0.0'f32
    check inverseErf(1.0'f32) == float32(Inf)
    check inverseErf(-1.0'f32) == -float32(Inf)

  test "erfcx float32":
    check approxEq32(erfcx(0.0'f32), 1.0'f32)
    check approxEq32(erfcx(1.0'f32), float32(erfcx(1.0)))
    check approxEq32(erfcx(10.0'f32), float32(erfcx(10.0)))

  test "erfcx float32 negative":
    check approxEq32(erfcx(-1.0'f32), float32(erfcx(-1.0)))

  test "erfc float32 asymptotic branch (x > 2)":
    check approxEq32(erfc(3.0'f32), float32(erfc(3.0)))
    check approxEq32(erfc(5.0'f32), float32(erfc(5.0)))

  test "digamma float32":
    check approxEq32(digamma(1.0'f32), float32(digamma(1.0)))
    check approxEq32(digamma(2.0'f32), float32(digamma(2.0)))
    check approxEq32(digamma(0.5'f32), float32(digamma(0.5)))
    check isNaN(digamma(0'f32))

  test "trigamma float32":
    check approxEq32(trigamma(1.0'f32), float32(trigamma(1.0)))
    check approxEq32(trigamma(2.0'f32), float32(trigamma(2.0)))
    check isNaN(trigamma(0'f32))

  test "lowerIncompleteGamma float32":
    let a = 2.0'f32
    let x = 0.5'f32
    check approxEq32(lowerIncompleteGamma(a, x),
                     float32(lowerIncompleteGamma(float(a), float(x))))

  test "regularizedLowerIncompleteGamma float32":
    let a = 2.0'f32
    let x = 0.5'f32
    check approxEq32(regularizedLowerIncompleteGamma(a, x),
                     float32(regularizedLowerIncompleteGamma(float(a), float(x))))

  test "regularized lower + upper incomplete gamma == 1 (float32)":
    let p = regularizedLowerIncompleteGamma(3.5'f32, 2.0'f32)
    let q = regularizedUpperIncompleteGamma(3.5'f32, 2.0'f32)
    check approxEq32(p + q, 1.0'f32)

  test "beta float32":
    check approxEq32(beta(2.0'f32, 2.0'f32), float32(beta(2.0, 2.0)))

  test "lbeta float32":
    check approxEq32(lbeta(2.0'f32, 2.0'f32), float32(lbeta(2.0, 2.0)))

  test "lowerIncompleteBeta float32":
    check approxEq32(
      lowerIncompleteBeta(2.0'f32, 2.0'f32, 0.5'f32),
      float32(lowerIncompleteBeta(2.0, 2.0, 0.5)))

  test "regularizedLowerIncompleteBeta float32":
    check approxEq32(
      regularizedLowerIncompleteBeta(2.0'f32, 2.0'f32, 0.5'f32),
      float32(regularizedLowerIncompleteBeta(2.0, 2.0, 0.5)))

  test "regularized lower + upper incomplete beta == 1 (float32)":
    let p = regularizedLowerIncompleteBeta(2.0'f32, 2.0'f32, 0.5'f32)
    let q = regularizedUpperIncompleteBeta(2.0'f32, 2.0'f32, 0.5'f32)
    check approxEq32(p + q, 1.0'f32)

  test "inverse regularized lower incomplete gamma round-trip (float32)":
    let a = 2.0'f32
    let x = 0.5'f32
    let p = regularizedLowerIncompleteGamma(a, x)
    let x2 = inverseRegularizedLowerIncompleteGamma(a, p)
    check approxEq32(x2, x)

  test "inverse regularized lower incomplete beta round-trip (float32)":
    let a = 2.0'f32
    let b = 3.0'f32
    let x = 0.5'f32
    let p = regularizedLowerIncompleteBeta(a, b, x)
    let x2 = inverseRegularizedLowerIncompleteBeta(a, b, p)
    check approxEq32(x2, x)

  test "NaN propagation (float32)":
    let nan32 = float32(NaN)
    check isNaN(lowerIncompleteGamma(2.0'f32, nan32))
    check isNaN(lowerIncompleteBeta(2.0'f32, 2.0'f32, nan32))
    check isNaN(inverseErf(nan32))
    check isNaN(digamma(nan32))
