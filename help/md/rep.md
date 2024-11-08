## name
rep
## title
Replicate a value
## description
'rep' creates a vector of 'n' copies of the value 'x'.
## details
'rep' creates a vector of 'n' elements, each with value 'x', preserving the type of 'x' in the returned vector.

Note that the documentation describes the case where `x` is an integer; other variable types (e.g. Real, String, Bool) are also supported.

`n` will be rounded down to the nearest integer.

## authors
Sebastian Hoehna
## see_also
simplex
v
## example
	rep(0.1, 3)
	
## references
