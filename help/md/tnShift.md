## name
tnShift
## title
A shifted distribution
## description
Shifts a given distribution.
## details
If X ~ dist then tnShift(dist, d) is the distribution of X + d

## authors
Ben Redelings
## see_also
tnScale
## example
    x ~ tnShift(dExponential(1),2)       # An exponential variable starting at 2.
    x ~ dnExponential(1) |> tnShift(2)   # Expressed using pipes.

## references
