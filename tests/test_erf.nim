import math
import times
import unittest

import special_functions


suite "Erf":
  
  #[
    Tested against https://keisan.casio.com/exec/system/1180573448
    All implementations are accurate to 15
  ]#
  
  const r1 = 15

  setup:
    let t0 = getTime()

  teardown:
    echo "\n  RUNTIME: ", getTime() - t0

  test "inverse_erf(0.25)":
    check inverse_erf(0.25).round(r1) == 0.225312055012178104725.round(r1)

  test "inverse_erf(0.50)":
    check inverse_erf(0.50).round(r1) == 0.4769362762044698733814.round(r1)

  test "inverse_erf(0.75)":
    check inverse_erf(0.75).round(r1) == 0.8134198475976185416903.round(r1)
