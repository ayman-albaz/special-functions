import math
import times
import unittest

import special_functions


suite "Beta":
  
  #[
    Tested against https://keisan.casio.com/exec/system/1180573395
    All implementations are accurate to 14 digits
  ]#
  
  setup:
    let t0 = getTime()

  teardown:
    echo "\n  RUNTIME: ", getTime() - t0

  test "Beta(1.0, 2.0)":
    check beta(1.0, 2.0).round(14) == 0.5.round(14)

  test "Beta(2.0, 1.0)":
    check beta(2.0, 1.0).round(14) == 0.5.round(14)

  test "Beta(2.0, 2.0)":
    check beta(2.0, 2.0).round(14) == 0.1666666666666666666667.round(14)

  test "LowerIncompleteBeta(2.0, 2.0, 0.25)":
    check lower_incomplete_beta(2.0, 2.0, 0.25).round(14) == 0.02604166666666666666667.round(14)

  test "LowerIncompleteBeta(2.0, 2.0, 0.5)":
    check lower_incomplete_beta(2.0, 2.0, 0.5).round(14) == 0.08333333333333333333333.round(14)

  test "LowerIncompleteBeta(2.0, 2.0, 0.75)":
    check lower_incomplete_beta(2.0, 2.0, 0.75).round(14) == 0.140625.round(14)

  test "RegularizedLowerIncompleteBeta(2.0, 2.0, 0.25)":
    check regularized_lower_incomplete_beta(2.0, 2.0, 0.25).round(14) == 0.15625.round(14)

  test "RegularizedLowerIncompleteBeta(2.0, 2.0, 0.5)":
    check regularized_lower_incomplete_beta(2.0, 2.0, 0.5).round(14) == 0.5.round(14)

  test "RegularizedLowerIncompleteBeta(2.0, 2.0, 0.75)":
    check regularized_lower_incomplete_beta(2.0, 2.0, 0.75).round(14) == 0.84375.round(14)

  test "InverseLowerIncompleteBeta(0.5, 0.5, 0.01)":
    check inverse_lower_incomplete_beta(0.5, 0.5, 0.01).round(14) == 2.499979166736110987103e-5.round(14)

  test "InverseLowerIncompleteBeta(1.0, 1.0, 0.01)":
    check inverse_lower_incomplete_beta(1.0, 1.0, 0.01).round(14) == 0.01.round(14)

  test "InverseLowerIncompleteBeta(2.0, 2.0, 0.01)":
    check inverse_lower_incomplete_beta(2.0, 2.0, 0.01).round(14) == 0.1490169525453262588864.round(14)

  test "InverseLowerIncompleteBeta(2.0, 2.0, 0.05)":
    check inverse_lower_incomplete_beta(2.0, 2.0, 0.05).round(14) == 0.3632574910905676135767.round(14)

  test "InverseLowerIncompleteBeta(2.0, 2.0, 0.10)":
    check inverse_lower_incomplete_beta(2.0, 2.0, 0.10).round(14) == 0.5670689228522682362543.round(14)
