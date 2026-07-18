import std/math
import std/unittest
import special_functions/realmath

func approxEq(a, b: float32): bool =
  abs(a - b) < 1e-5'f32 * max(1.0'f32, max(abs(a), abs(b)))

suite "realmath — native float32 primitives":

  test "lgamma float32 positive values":
    let expected = float32(lgamma(1.0))
    check approxEq(lgamma(1.0'f32), expected)
    check approxEq(lgamma(2.0'f32), float32(lgamma(2.0)))
    check approxEq(lgamma(3.0'f32), float32(lgamma(3.0)))
    check approxEq(lgamma(5.0'f32), float32(lgamma(5.0)))
    check approxEq(lgamma(10.0'f32), float32(lgamma(10.0)))
    check approxEq(lgamma(20.0'f32), float32(lgamma(20.0)))

  test "lgamma float32 near zero":
    let f32 = lgamma(0.001'f32)
    let f64 = float32(lgamma(0.001))
    check approxEq(f32, f64)

  test "lgamma float32 reflection for negative values":
    let f32 = lgamma(-0.5'f32)
    let f64 = float32(lgamma(-0.5))
    check approxEq(f32, f64)

  test "lgamma float32 at pole returns Inf":
    let r0 = lgamma(0'f32)
    check r0 == float32(Inf)
    let rn1 = lgamma(-1'f32)
    check rn1 == float32(Inf)

  test "gamma float32 positive values":
    check approxEq(gamma(1.0'f32), float32(gamma(1.0)))
    check approxEq(gamma(2.0'f32), float32(gamma(2.0)))
    check approxEq(gamma(3.0'f32), float32(gamma(3.0)))
    check approxEq(gamma(4.0'f32), float32(gamma(4.0)))
    check approxEq(gamma(5.0'f32), float32(gamma(5.0)))

  test "gamma float32 at pole returns Inf":
    check gamma(0'f32) == float32(Inf)

  test "erf float32":
    check approxEq(erf(0.0'f32), float32(erf(0.0)))
    check approxEq(erf(0.5'f32), float32(erf(0.5)))
    check approxEq(erf(1.0'f32), float32(erf(1.0)))
    check approxEq(erf(2.0'f32), float32(erf(2.0)))
    check approxEq(erf(-1.0'f32), float32(erf(-1.0)))
    check approxEq(erf(-2.0'f32), float32(erf(-2.0)))

  test "erfc float32":
    check approxEq(erfc(0.0'f32), float32(erfc(0.0)))
    check approxEq(erfc(0.5'f32), float32(erfc(0.5)))
    check approxEq(erfc(1.0'f32), float32(erfc(1.0)))
    check approxEq(erfc(2.0'f32), float32(erfc(2.0)))
    check approxEq(erfc(3.0'f32), float32(erfc(3.0)))

  test "erfc float32 + erf float32 == 1":
    check approxEq(erf(0.5'f32) + erfc(0.5'f32), 1.0'f32)

  test "erfc float32 for negative x":
    check approxEq(erfc(-1.0'f32), float32(erfc(-1.0)))
    check approxEq(erfc(-2.0'f32), float32(erfc(-2.0)))

  test "infT produces correct float32 Inf":
    let fi = infT(float32)
    check fi == float32(Inf)
    check fi > 1e30'f32

  test "infT produces correct float64 Inf":
    let fi = infT(float64)
    check fi == Inf
    check fi > 1e300

  test "nanT produces correct float32 NaN":
    check isNaN(nanT(float32))

  test "nanT produces correct float64 NaN":
    check isNaN(nanT(float64))

  test "piT produces correct float32 PI":
    check abs(piT(float32) - 3.141592'f32) < 1e-5'f32

  test "piT produces correct float64 PI":
    check abs(piT(float64) - PI) < 1e-15
