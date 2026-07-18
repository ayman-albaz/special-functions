import std/math
import std/unittest

import special_functions

func approxEq(a, b: float): bool =
  abs(a - b) < 1e-13 * max(1.0, max(abs(a), abs(b)))

suite "Erf":

  test "erf is available":
    check approxEq(erf(1.0), 0.8427007929497149)

  test "erfc is available":
    check approxEq(erfc(1.0), 0.1572992070502851)

  test "erf + erfc == 1":
    check approxEq(erf(1.5) + erfc(1.5), 1.0)

  test "inverseErf(0.0) == 0":
    check inverseErf(0.0) == 0.0

  test "inverseErf(0.25)":
    check approxEq(inverseErf(0.25), 0.2253120550121781)

  test "inverseErf(0.50)":
    check approxEq(inverseErf(0.50), 0.4769362762044699)

  test "inverseErf(0.75)":
    check approxEq(inverseErf(0.75), 0.8134198475976185)

  test "inverseErf(-0.25)":
    check approxEq(inverseErf(-0.25), -0.2253120550121781)

  test "inverseErf(-0.50)":
    check approxEq(inverseErf(-0.50), -0.4769362762044699)

  test "inverseErf(-0.75)":
    check approxEq(inverseErf(-0.75), -0.8134198475976185)

  test "inverseErf(0.99)":
    check approxEq(inverseErf(0.99), 1.821386367718449)

  test "inverseErf(-0.99)":
    check approxEq(inverseErf(-0.99), -1.821386367718449)

  test "inverseErf(0.999)":
    check approxEq(inverseErf(0.999), 2.326753765513524)

  test "inverseErf(1.0) == Inf":
    check inverseErf(1.0) == Inf

  test "inverseErf(-1.0) == -Inf":
    check inverseErf(-1.0) == -Inf

  test "inverseErf(x > 1) returns NaN":
    check isNaN(inverseErf(1.5))

  test "inverseErf(x < -1) returns NaN":
    check isNaN(inverseErf(-1.5))

  test "inverseErf(NaN) returns NaN":
    check isNaN(inverseErf(NaN))

  test "erfcx(0.0) == 1":
    check approxEq(erfcx(0.0), 1.0)

  test "erfcx(1.0)":
    check approxEq(erfcx(1.0), 0.427583576155807)

  test "erfcx(-1.0)":
    check approxEq(erfcx(-1.0), 5.008980080762283)

  test "erfcx(5.0)":
    check approxEq(erfcx(5.0), 0.1107046377330686)

  test "erfcx(10.0)":
    check approxEq(erfcx(10.0), 0.0561409927438226)

  test "erfcx(-5.0)":
    check approxEq(erfcx(-5.0), 1.4400979867466104e11)

  test "erfcx(20.0) — asymptotic branch":
    check approxEq(erfcx(20.0), 0.02817434874105132)

  test "erfcx(NaN) returns NaN":
    check isNaN(erfcx(NaN))

  test "erfcx round-trip for x > 0":
    let x = 3.0
    check approxEq(erfcx(x) / exp(x * x), erfc(x))

  test "erfcx at crossover (5.9, 6.0, 6.1)":
    check approxEq(erfcx(5.9), 0.09430713614832685)
    check approxEq(erfcx(6.0), 0.09277656780053835)
    check approxEq(erfcx(6.1), 0.09129430036868366)

  test "erfcx(15.0)":
    check approxEq(erfcx(15.0), 0.03752960638850578)

  test "erfcx(25.0)":
    check approxEq(erfcx(25.0), 0.02254957243264136)

  test "erfcx(50.0)":
    check approxEq(erfcx(50.0), 0.011281536265323772)

  test "erfcx(100.0)":
    check approxEq(erfcx(100.0), 0.005641613782989433)

  test "erfcx(Inf)":
    check erfcx(Inf) == 0.0

  test "erfcx large negative returns Inf":
    check erfcx(-30.0) == Inf

  test "inverseErf near 1.0":
    check approxEq(inverseErf(0.999999), 3.4589107372728956)

  test "inverseErf near -1.0":
    check approxEq(inverseErf(-0.999999), -3.4589107372728956)

  test "erf(Inf)":
    check erf(Inf) == 1.0

  test "erf(-Inf)":
    check erf(-Inf) == -1.0
