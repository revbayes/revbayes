## name
tnExp
## title
Exp-transformed distribution
## description
Exp-transforms a given distribution.
## details
If X ~ dist then tnExp(dist) is the distribution of exp(X).

The distribution `dist` can be either univariate (dnNormal) or
multivariate (dnMultivariateNormal).

This turns out to be the same as dnLog(dist), which provides a distribution
that has distribution `dist` on the log-scale.

## authors
Ben Redelings
## see_also
tnLog, tnLogit, tnInvlogit
## example
    x ~ tnExp(dnNormal(0,1))          # Draw from the log-Normal distribution
    x ~ dnNormal(0,1) |> tnExp()      # Expressed using pipes.
    x ~ dnLognormal(0,1)              # This is equivalent.
    y ~ dnNormal(0,1)
    x := exp(y)                       # This is also equivalent.

    x ~ tnExp(dnGamma(2,3))           # There is no equivalent for this.
    x ~ dnIID(10,tnExp(dnGamma(2,3))) # Draw 10 log-Gamma(2,3) random variables.

    mu = [1.0, 2.0, 3.0, 4.0]
    Sigma ~ dnWishart(df=4, kappa=2, dim=4)
    x ~ dnMultivariateNormal(mu,Sigma) |> tnExp()
## references
