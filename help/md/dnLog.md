## name
dnLog
## title
Log-scaled distribution
## description
Log-scales a given distribution.
## details
If X ~ dist then exp(X) ~ dnLog(dist)

This provides a way to construct distributions like dnLognormal and
dnLoguniform directly from the underlying distribution in log-space.
It can therefore express distributions that are not directly implemented.

The distribution `dist` can be either univariate (dnNormal) or
multivariate (dnMultivariateNormal).

## authors
Ben Redelings
## see_also
dnLognormal
dnLoguniform
## example
    x ~ dnLog(dnNormal(0,1))          # Draw from the log-Normal distribution
    x ~ dnNormal(0,1) |> dnLog()      # Expressed using pipes.
    x ~ dnLognormal(0,1)              # This is equivalent.
    y ~ dnNormal(0,1)
    x := exp(y)                       # This is also equivalent.

    x ~ dnLog(dnGamma(2,3))           # There is no equivalent for this.
    x ~ dnIID(10,dnLog(dnGamma(2,3))) # Draw 10 log-Gamma(2,3) random variables.

    mu = [1.0, 2.0, 3.0, 4.0]
    Sigma ~ dnWishart(df=4, kappa=2, dim=4)
    x ~ dnLog(dnMultivariateNormal(mu,Sigma))
## references
