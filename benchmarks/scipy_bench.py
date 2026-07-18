#!/usr/bin/env python3
"""Performance comparison: special_functions (Nim) vs SciPy (Python).

Requires: pip install scipy
"""

import time
import scipy.special as sp

N = 100_000


def bench(name, fn, *args):
    start = time.monotonic()
    for _ in range(N):
        fn(*args)
    elapsed = time.monotonic() - start
    ms = elapsed * 1000.0
    print(f"{name}: {ms:.3f} ms")


print(f"SciPy — Python benchmark ({N} iterations each)")
print("-" * 50)

# ---- Gamma ----
bench("gammainc(0.5, 2)          [lower reg]", sp.gammainc, 0.5, 2.0)
bench("gammaincc(0.5, 2)         [upper reg]", sp.gammaincc, 0.5, 2.0)

# scipy has gammaincinv / gammainccinv — inverse of regularized
bench("gammaincinv(0.5, 2)       [lower inv reg]", sp.gammaincinv, 0.5, 2.0)

# scipy lower_incomplete_gamma → gammainc * gamma(a), not a direct name
# regularized versions are the primary API
bench("gamma(2)                   [scalar gamma]", sp.gamma, 2.0)
bench("gammaln(5)                [log gamma]", sp.gammaln, 5.0)

# ---- Beta ----
bench("beta(2, 2)                [beta]", sp.beta, 2.0, 2.0)
bench("betainc(2, 2, 0.5)        [lower reg]", sp.betainc, 2.0, 2.0, 0.5)
bench("betaincinv(2, 2, 0.5)     [lower inv reg]", sp.betaincinv, 2.0, 2.0, 0.5)

# ---- Erf ----
bench("erf(1)", sp.erf, 1.0)
bench("erfcx(3)", sp.erfcx, 3.0)
bench("erfinv(0.5)               [inverse erf]", sp.erfinv, 0.5)

# ---- Polygamma ----
bench("digamma(1)", sp.digamma, 1.0)
bench("polygamma(1, 1)           [trigamma]", sp.polygamma, 1, 1.0)

print("-" * 50)
print("Note: SciPy gammainc* returns regularized values.")
print("      multiply by gamma(a) for unregularized comparison.")
