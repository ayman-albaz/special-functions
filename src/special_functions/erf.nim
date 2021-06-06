type
    Submodule* = object
        name*: string

proc initSubmodule*(): Submodule =
      ## Initialises a new ``Submodule`` object.
      Submodule(name: "Anonymous")


import math


{.nanChecks: on, infChecks: on.}


proc inverse_erf*(x: float): float =
    #[
        Inverse error function.
        Based off of: http://libit.sourceforge.net/math_8c-source.html
    ]#
    const
        inv_erf_a3 = -0.140543331
        inv_erf_a2 = 0.914624893
        inv_erf_a1 = -1.645349621
        inv_erf_a0 = 0.886226899

        inv_erf_b4 = 0.012229801
        inv_erf_b3 = -0.329097515
        inv_erf_b2 = 1.442710462
        inv_erf_b1 = -2.118377725
        inv_erf_b0 = 1

        inv_erf_c3 = 1.641345311
        inv_erf_c2 = 3.429567803
        inv_erf_c1 = -1.62490649
        inv_erf_c0 = -1.970840454

        inv_erf_d2 = 1.637067800
        inv_erf_d1 = 3.543889200
        inv_erf_d0 = 1
    var 
        x_sign: float
        x_copy: float
        x2: float
        r: float
        y: float

    if x > 1:
        return NaN

    if x == 0:
        return 0.0

    elif x > 0:
        x_sign = 1.0
        x_copy = x

    elif x < 0:
        x_sign = -1.0
        x_copy = -x

    if x_copy <= 0.7:
        x2 = x_copy * x_copy
        r = x_copy * (((inv_erf_a3 * x2 + inv_erf_a2) * x2 + inv_erf_a1) * x2 + inv_erf_a0)
        r /= (((inv_erf_b4 * x2 + inv_erf_b3) * x2 + inv_erf_b2) * x2 + inv_erf_b1) * x2 + inv_erf_b0

    elif x_copy > 0.7:
        y = sqrt(-ln((1 - x_copy) / 2))
        r = (((inv_erf_c3 * y + inv_erf_c2) * y + inv_erf_c1) * y + inv_erf_c0)
        r /= ((inv_erf_d2 * y + inv_erf_d1) * y + inv_erf_d0)

    r = r * x_sign
    x_copy = x_copy * x_sign

    r -= (erf(r) - x_copy) / (2 / sqrt(PI) * exp(-r * r))
    r -= (erf(r) - x_copy) / (2 / sqrt(PI) * exp(-r * r))

    return r
