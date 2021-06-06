type
    Submodule* = object
        name*: string

proc initSubmodule*(): Submodule =
      ## Initialises a new ``Submodule`` object.
      Submodule(name: "Anonymous")


import math
import math_utils


{.nanChecks: on, infChecks: on.}


proc incomplete_gamma_ctdfrac_closure(a, x: float): ((proc(i: int): float {.closure, noSideEffect, gcsafe, locks: 0.}), (proc(i: int): float {.noSideEffect, gcsafe, locks: 0.})) =
    let 
        fa = proc(i: int): float = 
            return -i.float * (i.float - a)

        fb = proc(i: int): float = 
            return x - a + 1.0 + 2.0 * i.float

    return (fa, fb)


proc incomplete_gamma(a, x: float, lower, reg: bool): float = 
    #[
        Based off of Jernej Kovačič's implementation: https://github.com/jkovacic/math/
    ]#
    if a < EPS or x < EPS:
        raise newException(ValueError, "Must not be: a < EPS or b < EPS or x < 0.0 or x > 1")

    var
        ginc: float = exp(-x) * pow(x, a)
        i: int = 1
        G: float
        term: float
        at: float
        fa, fb: ((proc(i: int): float {.closure, noSideEffect, gcsafe, locks: 0.}))
    (fa, fb) = incomplete_gamma_ctdfrac_closure(a, x)

    if x > (a + 1):
        if lower == false and reg == false:
            G = 0.0
        else:
            G = gamma(a)
        ginc /=  ctdfrac_eval(fa, fb)
        if lower == true:
            ginc = G - ginc
        if reg == true:
            ginc /= G

    else:
        if lower == true and reg == false:
            G = 0.0
        else:
            G = gamma(a)
        ginc /= a
        term = ginc
        at = a
        while abs(term) > TOL and i < MAX_ITER:
            at += 1.0
            term *= x / at
            ginc += term
            i += 1

        # check if the algorithm has converged:
        if i >= MAX_ITER:
            raise newException(ValueError, "Algorithm has not converged")

        if lower == false:
            ginc = G - ginc
        if reg == true:
            ginc /= G
    
    return ginc


proc lower_incomplete_gamma*(a, x: float): float =
    return incomplete_gamma(a, x, true, false)


proc upper_incomplete_gamma*(a, x: float): float =
    return incomplete_gamma(a, x, false, false)


proc regularized_lower_incomplete_gamma*(a, x: float): float =
    return incomplete_gamma(a, x, true, true)


proc regularized_upper_incomplete_gamma*(a, x: float): float =
    return incomplete_gamma(a, x, false, true)


proc inverse_incomplete_gamma(a, g: float, lower, reg: bool): float =
    const
        c1 = 0.253
        c2 = 0.12

    # sanity check
    if a < EPS or g < 0.0:
        raise newException(ValueError, "a must be < 2.2204e-16 and g must be less than 0")

    let G = gamma(a)
    var
        p = g
        x = 1.0
        pa: float
        cpa: float
        xn: float = 0.0
        i: int = 0
        f: float

    if reg == false:
        p /= G

    if lower == false:
        p = 1.0 - p

    if abs(p) < EPS:
        return 0.0

    if p >= 1 - EPS:
        raise newException(ValueError, "p must be p >= 1 - 2.2204e-16")

    if a <= 1.0:
        pa = a * (c1 + c2 * a)
        cpa = 1.0 - pa

        if p < cpa:
            x = pow(p / cpa, 1.0 / a)
        else:
            x = 1.0 - ln((1.0 - p) / pa)

    else:
        x = -as26_2_22(p)
        # [Abramowitz & Stegun], section 26.4.17:
        x = 1.0 - (1.0 / (9.0 * a)) + (x / (3.0 * sqrt(a)))
        x  = a * x * x * x

    f = regularized_lower_incomplete_gamma(a, x) - p

    while abs(f) > TOL and i < MAX_ITER:
        xn = x - f * G * exp(x) / pow(x, (a - 1.0))

        # x must not go negative!
        if xn > EPS:
            x = xn
        else:
            x = x * 0.5

        f = regularized_lower_incomplete_gamma(a, x) - p
        i += 1

    # Has the algorithm converged?
    if i >= MAX_ITER:
        raise newException(ValueError, "Algorithm has not converged")

    return x


proc inverse_lower_incomplete_gamma*(a, x: float): float =
    return inverse_incomplete_gamma(a, x, true, false)


proc inverse_upper_incomplete_gamma*(a, x: float): float =
    return inverse_incomplete_gamma(a, x, false, false)


proc inverse_regularized_lower_incomplete_gamma*(a, x: float): float =
    return inverse_incomplete_gamma(a, x, true, true)


proc inverse_regularized_upper_incomplete_gamma*(a, x: float): float =
    return inverse_incomplete_gamma(a, x, false, true)
