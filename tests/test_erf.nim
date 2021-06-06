import math
import times
import unittest

import special_functions/erf


suite "Erf":
    
    #[
        Tested against https://keisan.casio.com/exec/system/1180573448
        All implementations are accurate to 15
    ]#
    
    setup:
        let t0 = getTime()

    teardown:
        echo "\n  RUNTIME: ", getTime() - t0

    test "inverse_erf(0.25)":
        check inverse_erf(0.25).round(15) == 0.225312055012178104725.round(15)

    test "inverse_erf(0.50)":
        check inverse_erf(0.50).round(15) == 0.4769362762044698733814.round(15)

    test "inverse_erf(0.75)":
        check inverse_erf(0.75).round(15) == 0.8134198475976185416903.round(15)
