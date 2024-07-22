## name
tnLog
## title
Log-transformed distribution
## description
Log-transforms a given distribution.
## details
If X ~ dist then tnLog(dist) is the distribution of log(X).

The distribution `dist` can be either univariate (dnExponential) or
multivariate (dnDirichlet).

This is NOT the same as dnLog(dist), which provides a distribution
that has distribution `dist` on the log-scale.

## authors
Ben Redelings
## see_also
tnExp, tnLogit, tnInvlogit
## example
    x ~ tnLog(dnExponential(1))       # The log of an Exponential random variable.
    x ~ dnExponential(1) |> tnLog()   # Expressed using pipes.

    y ~ dnExponential(1)
    x := log(y)                       # This is also equivalent.

    x ~ dnDirichlet([1,1,1,1]) |> tnLog()
## references
