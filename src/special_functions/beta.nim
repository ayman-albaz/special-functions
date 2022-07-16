import math
import mathutils


{.nanChecks: on, infChecks: on.}


# GLOBAL CONSTANTS
const
  EPS = 2.2204e-16
  TOL = 1e-15
  MAX_ITER = 10000


proc lbeta*(a, b: float): float =
  return lgamma(a) + lgamma(b) - lgamma(a + b)

proc beta*(a, b: float): float =
  return exp(lbeta(a, b))

proc incomplete_beta_ctdfrac_closure(a, b, x: float): ((proc(i: int): float {.closure, noSideEffect, gcsafe, locks: 0.}), (proc(i: int): float {.noSideEffect, gcsafe, locks: 0.})) =
  # Based off of Jernej Kovačič's implementation: https://github.com/jkovacic/math/
  let 
    fa =  proc(i: int): float = 
      let m: float = (i div 2).float
      if floorMod(i, 2) == 1:
        return -((a + m) * (a + b + m) * x) / (( a + i.float - 1) * (a + i.float))
      else:
        return (m * (b - m) * x) / ((a + i.float - 1) * (a + i.float))

    fb = proc(i: int): float {.closure.} = 
      return 1.0

  return (fa, fb)

proc incomplete_beta(a, b, x: float, lower, reg: bool): float =
  # Based off of Jernej Kovačič's implementation: https://github.com/jkovacic/math/
  if a < EPS or b < EPS or x < 0.0 or x > 1:
    raise newException(ValueError, "Must not be: a < EPS or b < EPS or x < 0.0 or x > 1")

  var
    xlarge: bool = x > ((a + 1) / (a + b + 2))
    binc: float = 0.0
    xn: float
    an: float
    bn: float
    B: float
    fa, fb: ((proc(i: int): float {.closure, noSideEffect, gcsafe, locks: 0.}))

  if xlarge:
    xn = 1.0 - x
    an = b
    bn = a
  else:
    xn = x
    an = a
    bn = b

  (fa, fb) = incomplete_beta_ctdfrac_closure(an, bn, xn)

  if reg == true or (xlarge == false and lower == false) or (xlarge == true and lower == true):
    B = beta(a, b)
  else:
    B = 0.0

  if abs(xn) < EPS:
    if xlarge == true:
      binc = B
    else:
      binc = 0.0

  else:
    binc = pow(xn, an) * pow(1.0 - xn, bn) / an
    binc /= ctdfrac_eval(fa, fb)

  if (xlarge == false and lower == false) or (xlarge == true and lower == true):
    binc = B - binc

  if reg == true:
    binc /= B

  return binc


proc lower_incomplete_beta*(a, b, x: float): float =
  return incomplete_beta(a, b, x, true, false)


proc upper_incomplete_beta*(a, b, x: float): float =
  return incomplete_beta(a, b, x, false, false)


proc regularized_lower_incomplete_beta*(a, b, x: float): float =
  return incomplete_beta(a, b, x, true, true)


proc regularized_upper_incomplete_beta*(a, b, x: float): float =
  return incomplete_beta(a, b, x, false, true)


proc inverse_incomplete_beta(a, b, y: float, lower, reg: bool): float =
  # Based off of Jernej Kovačič's implementation: https://github.com/jkovacic/math/
  const
    c1 = 0.253
    c2 = 0.12

  # sanity check
  if a < EPS or b < EPS or y < 0:
    raise newException(ValueError, "a < EPS or b < EPS or y < 0")

  let B = beta(a, b)
  var
    i = 0
    p = y
    x = 0.0
    xn = 0.0
    pa: float
    cpa: float
    lambda: float
    ta: float
    tb: float
    h: float
    w: float
    S: float
    f: float

  if reg == false:
    p /= B

  if lower == false:
    p = 1 - p

  if abs(p) < EPS:
    return 0.0

  if p >= (1.0 - EPS):
    raise newException(ValueError, "p >= (1.0 - EPS)")

  if a >= 1.0 and b >= 1.0:
    x = as26_2_22(p)
    lambda = (x * x - 3.0) / 6.0
    ta = 1 / (2 * a - 1)
    tb = 1 / (2 * b - 1)
    h = 2 / (ta + tb)
    w = x * sqrt(h + lambda) / h - (tb - ta) * (lambda + 5.0 / 6.0 - 2.0 / (3.0 * h))
    x = a / (a + b * exp(2.0 * w))
  else:
    ta = pow(a / (a + b), a) / a
    tb = pow(b / (a + b), b) / b
    S = ta + tb
    if p < ta / S:
      x = pow(p * S * a, 1.0 / a)
    else:
      x = pow(1.0 - p * S * b, 1.0 / b)

  f = regularized_lower_incomplete_beta(a, b, x) - p
  while abs(f) > TOL and i < MAX_ITER:
    xn = x - f * B / (pow(x, (a - 1)) * pow((1 - x), (b - 1)))
    if xn < EPS:
      x *= 0.5
    elif xn > 1.0 - EPS:
      x = (1.0 + x) * 0.5
    else:
      x = xn
    f = regularized_lower_incomplete_beta(a, b, x) - p
    i += 1

  # Has the algorithm converged?
  if i >= MAX_ITER:
    raise newException(ValueError, "Algorithm has not converged")

  return x

proc inverse_lower_incomplete_beta*(a, b, x: float): float =
  return inverse_incomplete_beta(a, b, x, true, false)

proc inverse_upper_incomplete_beta*(a, b, x: float): float =
  return inverse_incomplete_beta(a, b, x, false, false)

proc inverse_regularized_lower_incomplete_beta*(a, b, x: float): float =
  return inverse_incomplete_beta(a, b, x, true, true)

proc inverse_regularized_upper_incomplete_beta*(a, b, x: float): float =
  return inverse_incomplete_beta(a, b, x, false, true)
