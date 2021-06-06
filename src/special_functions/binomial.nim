type
    Submodule* = object
        name*: string

proc initSubmodule*(): Submodule =
      ## Initialises a new ``Submodule`` object.
      Submodule(name: "Anonymous")


{.nanChecks: on, infChecks: on.}


proc binomial_coefficient*(n, k: int): int {.inline.} =
    #[ 
        Based off of: https://www.geeksforgeeks.org/space-and-time-efficient-binomial-coefficient/
    ]#
    var 
        res: int = 1
        new_k: int = k

    if k > n - k:
        new_k = n - k
    
    for i in 0..<new_k:
        res = res * (n - i)
        res = res div (i + 1)
    return res
