import math
import times
import unittest

import special_functions


suite "Beta":
  
  #[
    Tested against https://keisan.casio.com/exec/system/1180573395
    All implementations are accurate to 14 digits
  ]#
  
  const r1 = 14

  setup:
    let t0 = getTime()

  teardown:
    echo "\n  RUNTIME: ", getTime() - t0

  test "Beta(1.0, 2.0)":
    check beta(1.0, 2.0).round(r1) == 0.5.round(r1)

  test "Beta(2.0, 1.0)":
    check beta(2.0, 1.0).round(r1) == 0.5.round(r1)

  test "Beta(2.0, 2.0)":
    check beta(2.0, 2.0).round(r1) == 0.1666666666666666666667.round(r1)

  test "LowerIncompleteBeta(2.0, 2.0, 0.25)":
    check lower_incomplete_beta(2.0, 2.0, 0.25).round(r1) == 0.02604166666666666666667.round(r1)

  test "LowerIncompleteBeta(2.0, 2.0, 0.5)":
    check lower_incomplete_beta(2.0, 2.0, 0.5).round(r1) == 0.08333333333333333333333.round(r1)

  test "LowerIncompleteBeta(2.0, 2.0, 0.75)":
    check lower_incomplete_beta(2.0, 2.0, 0.75).round(r1) == 0.140625.round(r1)

  test "RegularizedLowerIncompleteBeta(2.0, 2.0, 0.25)":
    check regularized_lower_incomplete_beta(2.0, 2.0, 0.25).round(r1) == 0.15625.round(r1)

  test "RegularizedLowerIncompleteBeta(2.0, 2.0, 0.5)":
    check regularized_lower_incomplete_beta(2.0, 2.0, 0.5).round(r1) == 0.5.round(r1)

  test "RegularizedLowerIncompleteBeta(2.0, 2.0, 0.75)":
    check regularized_lower_incomplete_beta(2.0, 2.0, 0.75).round(r1) == 0.84375.round(r1)

  test "InverseLowerIncompleteBeta(0.5, 0.5, 0.01)":
    check inverse_lower_incomplete_beta(0.5, 0.5, 0.01).round(r1) == 2.499979166736110987103e-5.round(r1)

  test "InverseLowerIncompleteBeta(1.0, 1.0, 0.01)":
    check inverse_lower_incomplete_beta(1.0, 1.0, 0.01).round(r1) == 0.01.round(r1)

  test "InverseLowerIncompleteBeta(2.0, 2.0, 0.01)":
    check inverse_lower_incomplete_beta(2.0, 2.0, 0.01).round(r1) == 0.1490169525453262588864.round(r1)

  test "InverseLowerIncompleteBeta(2.0, 2.0, 0.05)":
    check inverse_lower_incomplete_beta(2.0, 2.0, 0.05).round(r1) == 0.3632574910905676135767.round(r1)

  test "InverseLowerIncompleteBeta(2.0, 2.0, 0.10)":
    check inverse_lower_incomplete_beta(2.0, 2.0, 0.10).round(r1) == 0.5670689228522682362543.round(r1)
