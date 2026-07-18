{.experimental: "strictFuncs".}
{.push raises: [].}

import std/math

## Generic ``float32`` / ``float64`` dispatch for basic special functions
## (``lgamma``, ``gamma``, ``erf``, ``erfc``) and per-type ``infT``,
## ``nanT``, ``piT`` helpers.
##
## On ``float64`` delegates to the C standard library via ``std/math``.
## On ``float32`` uses native Nim implementations (Lanczos approximation
## for ``lgamma``, A&S 7.1.26 for ``erf``).

func infT*(T: typedesc): T {.inline.} =
  ## Returns positive infinity for type ``T``.
  ##
  ## Internal helper; not exported through `special_functions.nim`.
  when T is float32: float32(Inf)
  elif T is float64: float64(Inf)
  else: T(Inf)

func nanT*(T: typedesc): T {.inline.} =
  ## Returns a quiet NaN for type ``T``.
  ##
  ## Internal helper; not exported through `special_functions.nim`.
  when T is float32: float32(NaN)
  elif T is float64: float64(NaN)
  else: T(NaN)

func piT*(T: typedesc): T {.inline.} =
  ## Returns π for type ``T``.
  ##
  ## Internal helper; not exported through `special_functions.nim`.
  when T is float32: float32(3.14159265358979)
  elif T is float64: float64(PI)
  else: T(PI)

func lgammaNative32(x: float32): float32 =
  if x <= 0.0'f32:
    if x == floor(x):
      return float32(Inf)
    let s = sin(float32(PI) * x)
    if s == 0.0'f32:
      return float32(Inf)
    return ln(float32(PI)) - ln(abs(s)) - lgammaNative32(1.0'f32 - x)

  if x == float32(Inf):
    return float32(Inf)

  # Lanczos approximation (g=7, n=8), accurate for float32.
  # lgamma(z) = 0.5·ln(2π) + (z-0.5)·ln(z+g-0.5) - (z+g-0.5) + ln(A(z-1))
  # A(t) = c0 + c1/(t+1) + c2/(t+2) + ... + c8/(t+8)
  const
    g = 7.0'f32
    halfLn2Pi = 0.9189385332046727'f32
    c0 =  0.9999999'f32
    c1 =  676.5204'f32
    c2 = -1259.139'f32
    c3 =  771.3234'f32
    c4 = -176.6150'f32
    c5 =  12.50734'f32
    c6 = -0.1385711'f32
    c7 =  9.98437e-6'f32
    c8 =  1.505633e-7'f32

  var z = x
  if z < 0.5'f32:
    return lgammaNative32(z + 1.0'f32) - ln(z)

  let t = z - 1.0'f32
  let s = c0 + (c1 / (t + 1.0'f32) +
                c2 / (t + 2.0'f32) +
                c3 / (t + 3.0'f32) +
                c4 / (t + 4.0'f32) +
                c5 / (t + 5.0'f32) +
                c6 / (t + 6.0'f32) +
                c7 / (t + 7.0'f32) +
                c8 / (t + 8.0'f32))
  halfLn2Pi + (z - 0.5'f32) * ln(z + g - 0.5'f32) - (z + g - 0.5'f32) + ln(s)

func erfNative32(x: float32): float32 =
  let sign = if x < 0.0'f32: -1.0'f32 else: 1.0'f32
  let z = abs(x)
  const
    p = 0.3275911'f32
    a1 =  0.254829592'f32
    a2 = -0.284496736'f32
    a3 =  1.421413741'f32
    a4 = -1.453152027'f32
    a5 =  1.061405429'f32
  let t = 1.0'f32 / (1.0'f32 + p * z)
  let y = ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t * exp(-z*z)
  sign * (1.0'f32 - y)

func erfcNative32(x: float32): float32 =
  if x < 0.0'f32:
    return 2.0'f32 - erfcNative32(-x)
  if x < 2.0'f32:
    return 1.0'f32 - erfNative32(x)
  let x2 = x * x
  let invX = 1.0'f32 / x
  var term = invX
  var s = term
  var k = 1.0'f32
  var prev: float32
  while k < 20.0'f32:
    prev = term
    term *= -(2.0'f32*k - 1.0'f32) / (2.0'f32 * x2)
    if abs(term) > abs(prev):
      break
    s += term
    k += 1.0'f32
  exp(-x2) * s / sqrt(float32(PI))

func lgamma*[T: SomeFloat](x: T): T =
  ## Natural log of the absolute value of the gamma function.
  ##
  ## On ``float64`` delegates to the C standard library; on ``float32``
  ## uses a Lanczos approximation (g=7, n=8).
  when T is float64:
    math.lgamma(x)
  elif T is float32:
    lgammaNative32(x)

func gamma*[T: SomeFloat](x: T): T =
  ## Gamma function :math:`\Gamma(x)`.
  ##
  ## On ``float64`` delegates to the C standard library; on ``float32``
  ## computes ``exp(lgamma(x))`` with an overflow guard.
  when T is float64:
    math.gamma(x)
  elif T is float32:
    let lg = lgammaNative32(x)
    if lg == float32(Inf):
      float32(Inf)
    else:
      exp(lg)

func erf*[T: SomeFloat](x: T): T =
  ## Error function :math:`\operatorname{erf}(x)`.
  ##
  ## On ``float64`` delegates to the C standard library; on ``float32``
  ## uses the A&S 7.1.26 rational approximation (error < 1.5e-7).
  when T is float64:
    math.erf(x)
  elif T is float32:
    erfNative32(x)

func erfc*[T: SomeFloat](x: T): T =
  ## Complementary error function :math:`\operatorname{erfc}(x) = 1 - \operatorname{erf}(x)`.
  ##
  ## On ``float64`` delegates to the C standard library; on ``float32``
  ## uses A&S 7.1.26 for ``|x| < 2`` and an asymptotic series for ``|x| ≥ 2``.
  when T is float64:
    math.erfc(x)
  elif T is float32:
    erfcNative32(x)

{.pop.}
