# Special Functions
Special functions is a Nim library for calculating 'special' mathematical functions that are not found in Nim math.


## Supported Functions
Below are the supported gamma functions:
```Nim
import special_functions/gamma

var x1: float

# Incomplete gamma functions
x1 = lower_incomplete_gamma(2.0, 0.25)
x1 = upper_incomplete_gamma(2.0, 0.25)
x1 = regularized_lower_incomplete_gamma(2.0, 0.25)
x1 = regularized_upper_incomplete_gamma(2.0, 0.25)

# Inverse Incomplete gamma functions
x1 = inverse_lower_incomplete_gamma(2.0, 0.25)
x1 = inverse_upper_incomplete_gamma(2.0, 0.25)
x1 = inverse_regularized_lower_incomplete_gamma(2.0, 0.25)
x1 = inverse_regularized_upper_incomplete_gamma(2.0, 0.25)
```

Below are the supported beta functions:
```Nim
import special_functions/beta

var x1: float

# Incomplete beta functions
x1 = lower_incomplete_beta(2.0, 2.0, 0.25)
x1 = upper_incomplete_beta(2.0, 2.0, 0.25)
x1 = regularized_lower_incomplete_beta(2.0, 2.0, 0.25)
x1 = regularized_upper_incomplete_beta(2.0, 2.0, 0.25)

# Inverse Incomplete beta functions
x1 = inverse_lower_incomplete_beta(2.0, 2.0, 0.25)
x1 = inverse_upper_incomplete_beta(2.0, 2.0, 0.25)
x1 = inverse_regularized_lower_incomplete_beta(2.0, 2.0, 0.25)
x1 = inverse_regularized_upper_incomplete_beta(2.0, 2.0, 0.25)
```

Below are the supported erf functions:
```Nim
import special_functions/erf

var x1: float

# Erf functions
x1 = inverse_erf(0.25)
```


## Future Directions:
There are really two directions this library can grow.
1. Write more special functions in Nim.
	- Pros: Library is in Nim.
	- Cons: Longer dev time, calculations are more error prone.
2. Implement all special functions using a C-wrapper around the Rmath library.
	- Pros: Reduced dev time, calculations are less error prone.
	- Cons: Library is written in C and hard to read (its filled with macros).

I don't have a preferred direction and I prefer to leave the answer for the Nim scientific community to answer.


## Accuracy
All functions in this library are accurate up-to 14 decimal places (float64).


## Performance
This library was written with accuracy as a top priority as opposed to performance, however almost all implementations here are faster than SciPy's implementations and equal to, slower, or faster than Julia's Special Functions implementations. 


## TODO
List is organized from most important to least important:
- Add more special functions on an as-need-bases.
- Some functions do not have a unittest because I couldn't find a source of truth to compare them against. Finding a source of truth and adding a unittest is the next step.

Performance, feature, and documentation PR's are always welcome.


## Acknowledgments
Many of the algorithms used here are based off of different peoples' code.
Thank you to:
- https://github.com/jkovacic/math/
- http://libit.sourceforge.net/math_8c-source.html


## Contact
I can be reached at aymanalbaz98@gmail.com
