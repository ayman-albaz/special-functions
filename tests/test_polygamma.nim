import std/math
import std/unittest

import special_functions

func approxEq(a, b: float): bool =
  abs(a - b) < 1e-13 * max(1.0, max(abs(a), abs(b)))

suite "Polygamma":

  test "digamma(1.0) == -γ":
    check approxEq(digamma(1.0), -0.5772156649015329)

  test "digamma(2.0) == 1 - γ":
    check approxEq(digamma(2.0), 0.4227843350984671)

  test "digamma(0.5)":
    check approxEq(digamma(0.5), -1.9635100260214235)

  test "digamma(10.0)":
    check approxEq(digamma(10.0), 2.251752589066721)

  test "digamma recurrence ψ(x+1) == ψ(x) + 1/x":
    check approxEq(digamma(3.5), digamma(2.5) + 1.0 / 2.5)

  test "digamma at negative integer returns NaN":
    check isNaN(digamma(-1.0))
    check isNaN(digamma(-2.0))
    check isNaN(digamma(-3.0))

  test "digamma at zero returns NaN":
    check isNaN(digamma(0.0))

  test "digamma at negative non-integer":
    check approxEq(digamma(-0.5), 0.03648997397857652)

  test "digamma at negative non-integer recurrence":
    let x = -1.3
    check approxEq(digamma(x), digamma(x + 1.0) - 1.0 / x)

  test "digamma(NaN) returns NaN":
    check isNaN(digamma(NaN))

  test "digamma at large x":
    check approxEq(digamma(50.0), 3.901989673427892)

  test "trigamma(1.0) == π²/6":
    check approxEq(trigamma(1.0), 1.6449340668482264)

  test "trigamma(2.0) == π²/6 - 1":
    check approxEq(trigamma(2.0), 0.6449340668482264)

  test "trigamma(0.5) == π²/2":
    check approxEq(trigamma(0.5), 4.934802200544679)

  test "trigamma recurrence ψ₁(x+1) == ψ₁(x) - 1/x²":
    check approxEq(trigamma(3.5), trigamma(2.5) - 1.0 / (2.5 * 2.5))

  test "trigamma at negative integer returns NaN":
    check isNaN(trigamma(-1.0))
    check isNaN(trigamma(-2.0))

  test "trigamma at zero returns NaN":
    check isNaN(trigamma(0.0))

  test "trigamma(NaN) returns NaN":
    check isNaN(trigamma(NaN))

  test "trigamma at negative non-integer":
    check approxEq(trigamma(-0.5), 8.934802200544679)
    let x = -1.3
    check approxEq(trigamma(x), trigamma(x + 1.0) + 1.0 / (x * x))

  test "digamma(Inf)":
    check digamma(Inf) == Inf

  test "trigamma(Inf)":
    check trigamma(Inf) == 0.0

  test "trigamma at large x — verify asymptotic accuracy":
    check approxEq(trigamma(50.0), 0.02020133322669713)

  test "trigamma(20.0)":
    check approxEq(trigamma(20.0), 0.05127082293520312)

  test "digamma(100)":
    check approxEq(digamma(100.0), 4.600161852738087)

  test "digamma(1000)":
    check approxEq(digamma(1000.0), 6.907255195648812)

  test "trigamma(100)":
    check approxEq(trigamma(100.0), 0.010050166663334)

  test "trigamma(1000)":
    check approxEq(trigamma(1000.0), 0.0010005001666666)

  test "digamma near pole (small positive x)":
    check approxEq(digamma(0.01), -100.5608854578687)

  test "trigamma near pole (small positive x)":
    check approxEq(trigamma(0.01), 10001.621213528315)
