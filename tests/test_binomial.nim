import times
import unittest

import special_functions/binomial


suite "Binomial":
    
    #[
        Tested against https://www.omnicalculator.com/math/binomial-coefficient
    ]#
    
    setup:
        let t0 = getTime()

    teardown:
        echo "\n  RUNTIME: ", getTime() - t0

    test "binomial_coefficient(20, 10)":
        check binomial_coefficient(20, 10) == 184756

    test "binomial_coefficient(20, 15)":
        check binomial_coefficient(20, 15) == 15504

    test "binomial_coefficient(20, 20)":
        check binomial_coefficient(20, 20) == 1
