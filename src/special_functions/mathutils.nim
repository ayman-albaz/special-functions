{.experimental: "strictFuncs".}
{.push raises: [].}

import std/math
import realmath

## Internal utilities shared across the library.
##
## Provides type-tuned constants (``epsT``, ``tolT``, ``lentzTinyT``,
## ``maxLogT``), the iteration cap ``MAX_ITER``, and Acklam's inverse
## normal CDF approximation (``as26_2_22``).  These symbols carry ``*``
## for sibling-module access but are **not** re-exported through
## ``special_functions.nim``.

const
  MAX_ITER* = 10000
    ## Maximum number of iterations for Newton and continued-fraction loops.

func epsT*(T: typedesc): T {.inline.} =
  ## Machine epsilon for type ``T``.
  when T is float32: 1.19209290e-7'f32
  elif T is float64: 2.220446049250313e-16
  else: T(2.220446049250313e-16)

func tolT*(T: typedesc): T {.inline.} =
  ## Convergence tolerance for type ``T`` (Newton / continued-fraction loops).
  when T is float32: 1e-6'f32
  elif T is float64: 1e-15
  else: T(1e-15)

func lentzTinyT*(T: typedesc): T {.inline.} =
  ## Floor value for Lentz continued-fraction denominators (prevents underflow).
  when T is float32: 1e-19'f32
  elif T is float64: 1e-30
  else: T(1e-30)

func maxLogT*(T: typedesc): T {.inline.} =
  ## Safe exponent ceiling; ``exp(maxLogT(T))`` fits in ``T`` without overflow.
  when T is float32: 85.0'f32
  elif T is float64: 700.0
  else: T(700.0)

func as26_2_22*[T: SomeFloat](p: T): T =
  ## Inverse normal CDF (Acklam's rational approximation).
  ## Relative error < 1.15e-9.  Used internally for initial guesses
  ## in inverse incomplete gamma/beta.
  ##
  ## Sign convention: returns positive for ``p < 0.5``, negative for
  ## ``p > 0.5`` (inverse of the standard probit convention), matching
  ## the prior A&S 26.2.22 implementation used by existing callers.

  if p <= T(0.0):
    return -infT(T)
  if p >= T(1.0):
    return infT(T)

  let
    a1 = T(-3.969683028665376e+01)
    a2 = T(2.209460984245205e+02)
    a3 = T(-2.759285104469687e+02)
    a4 = T(1.383577518672690e+02)
    a5 = T(-3.066479806614716e+01)
    a6 = T(2.506628277459239e+00)

    b1 = T(-5.447609879822406e+01)
    b2 = T(1.615858368580409e+02)
    b3 = T(-1.556989798598866e+02)
    b4 = T(6.680131188771972e+01)
    b5 = T(-1.328068155288572e+01)

    c1 = T(-7.784894002430293e-03)
    c2 = T(-3.223964580411365e-01)
    c3 = T(-2.400758277161838e+00)
    c4 = T(-2.549732539343734e+00)
    c5 = T(4.374664141464968e+00)
    c6 = T(2.938163982698783e+00)

    d1 = T(7.784695709041462e-03)
    d2 = T(3.224671290700398e-01)
    d3 = T(2.445134137142996e+00)
    d4 = T(3.754408661907416e+00)

    thresh = T(0.02425)
    two = T(2.0)
    one = T(1.0)

  let lowTail = p < thresh
  let highTail = p > one - thresh

  if lowTail:
    let q = sqrt(-two * ln(p))
    let num = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6)
    let den = ((((d1 * q + d2) * q + d3) * q + d4) * q + one)
    return num / den

  if highTail:
    let q = sqrt(-two * ln(one - p))
    let num = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6)
    let den = ((((d1 * q + d2) * q + d3) * q + d4) * q + one)
    return -(num / den)

  let q = p - T(0.5)
  let r = q * q
  let num = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q
  let den = ((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + one
  -(num / den)

{.pop.}
