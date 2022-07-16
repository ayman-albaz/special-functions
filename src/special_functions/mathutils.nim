import math


# GLOBAL CONSTANTS
const
  EPS* = 2.2204e-16
  TOL* = 1e-15
  MAX_ITER* = 10000

proc as26_2_22*(p: float): float =
  ##[ 
    Based off of Jernej Kovačič's implementation: https://github.com/jkovacic/math/
    Implement the formula 26.2.22 in [Abramowitz & Stegun], i.e
    return such `x` that approximately satisfies Q(x) ~= p.
  ]##
  # [Abramowitz & Stegun], section 26.2.22:
  const
    a0 = 2.30753
    a1 = 0.27061
    b0 = 0.99229
    b1 = 0.04481
  var
    pp: float
    t: float
    x: float

  if p >= 0.5:
    pp = 1.0 - p
  else:
    pp = p

  t = sqrt(-2.0 * ln(pp))
  x = t - (a0 + a1 * t) / (1.0 + t * (b0 + b1 * t))

  if p > 0.5:
    return -x
  else:
    return x


proc ctdfrac_eval*(fa: (proc(i: int): float {.closure, noSideEffect, gcsafe, locks: 0.}), fb: (proc(i: int): float {.closure, noSideEffect, gcsafe, locks: 0.})): float = 
  ##[
    Based off of Jernej Kovačič's implementation: https://github.com/jkovacic/math/
  ]##
  var
    f: float = fb(0)
    c: float = f
    d: float = 0.0
    Delta: float = 0.0
    j: int = 1
    a: float
    b: float

  while abs(Delta - 1) > TOL and j < MAX_ITER:
    # obtain 'aj' and 'bj'
    a = fa(j)
    b = fb(j)

    # dj = bj + aj * d_j-1
    d = b + a * d

    if abs(d) < EPS:
      d = EPS

    # cj = bj + aj/c_j-1
    c = b + a / c

    if abs(c) < EPS:
      c = EPS

    # dj = 1 / dj
    d = 1.0 / d

    # Delta_j = cj * dj
    Delta = c * d

    # fj = f_j-1 * Delta_j
    f *= Delta

    j += 1

  # check if the algorithm has converged:
  if j >= MAX_ITER:
    raise newException(ValueError, "Algorithm has not converged")

  return f
