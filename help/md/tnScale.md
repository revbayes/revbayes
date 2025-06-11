## name
tnScale
## title
A scaled distribution
## description
Rescales a given distribution.
## details
If X ~ dist then tnScale(dist, lambda) is the distribution of X * lambda

## authors
Ben Redelings
## see_also
tnScale
## example
    x ~ tnScale(dExponential(1),2)       # An Exponential(rate=0.5) random variable.
    x ~ dnExponential(1) |> tnScale(2)   # Expressed using pipes.

## references
