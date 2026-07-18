## Special mathematical functions not found in Nim's standard ``math``
## module.
##
## All fallible functions return ``NaN`` for out-of-domain inputs or when
## iterative algorithms fail to converge.  See the `README`_ for usage
## examples.
##
## All functions are generic over ``SomeFloat``, supporting both
## ``float32`` and ``float64``.  On ``float64``, ``lgamma``, ``gamma``,
## ``erf``, and ``erfc`` delegate directly to the C standard library
## implementations in ``std/math``.
##
## .. _README: https://github.com/ayman-albaz/special-functions#readme

import special_functions/[beta, erf, gamma, polygamma]
import special_functions/realmath
export
  realmath.lgamma,  realmath.gamma,
  realmath.erf,     realmath.erfc

export
  beta,        lbeta,
  lowerIncompleteBeta,     upperIncompleteBeta,
  regularizedLowerIncompleteBeta,  regularizedUpperIncompleteBeta,
  inverseLowerIncompleteBeta,      inverseUpperIncompleteBeta,
  inverseRegularizedLowerIncompleteBeta, inverseRegularizedUpperIncompleteBeta,
  inverseErf,  erfcx,
  lowerIncompleteGamma,    upperIncompleteGamma,
  regularizedLowerIncompleteGamma,  regularizedUpperIncompleteGamma,
  inverseLowerIncompleteGamma,      inverseUpperIncompleteGamma,
  inverseRegularizedLowerIncompleteGamma,  inverseRegularizedUpperIncompleteGamma,
  digamma,     trigamma
