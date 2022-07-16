import math
import times
import unittest

import special_functions


suite "Gamma":
  
  #[
    Tested against https://keisan.casio.com/exec/system/1180573444
    All implementations are accurate to 14 digits
  ]#
  
  setup:
    let t0 = getTime()

  teardown:
    echo "\n  RUNTIME: ", getTime() - t0

  test "LowerIncompleteGamma(2.0, 0.25)":
    check lower_incomplete_gamma(2.0, 0.25).round(14) == 0.02649902116074391469354.round(14)

  test "LowerIncompleteGamma(2.0, 0.50)":
    check lower_incomplete_gamma(2.0, 0.50).round(14) == 0.0902040104310498645943.round(14)

  test "LowerIncompleteGamma(2.0, 0.75)":
    check lower_incomplete_gamma(2.0, 0.75).round(14) == 0.1733585327032242625084.round(14)

  test "UpperIncompleteGamma(2.0, 0.25)":
    check upper_incomplete_gamma(2.0, 0.25).round(14) == 0.9735009788392560853065.round(14)

  test "UpperIncompleteGamma(2.0, 0.50)":
    check upper_incomplete_gamma(2.0, 0.50).round(14) == 0.9097959895689501354057.round(14)

  test "UpperIncompleteGamma(2.0, 0.75)":
    check upper_incomplete_gamma(2.0, 0.75).round(14) == 0.8266414672967757374916.round(14)

  test "InverseLowerIncompleteGamma(2.0, 0.25)":
    check inverse_lower_incomplete_gamma(2.0, 0.25).round(14) == 0.9612787631147771.round(14)
  
  test "InverseLowerIncompleteGamma(2.0, 0.50)":
    check inverse_lower_incomplete_gamma(2.0, 0.50).round(14) == 1.6783469900166612.round(14)

  test "InverseLowerIncompleteGamma(2.0, 0.75)":
    check inverse_lower_incomplete_gamma(2.0, 0.75).round(14) == 2.692634528889695.round(14)

