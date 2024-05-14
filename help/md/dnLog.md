## name
dnLog
## title
A log-transformed distribution
## description
Computes a log-transformed distribution.  If
        X ~ dist
then
        exp(X) ~ dnLog(dist)

This provides a way to construct distributions like dnLognormal and
dnLoguniform directly from the underlying distribution in log-space.
It can therefore express distributions that are not directly implemented.

## details
## authors
## see_also
dnLognormal
dnLoguniform
## example
    x ~ dnLog(dnNormal(0,1))          # Draw from the log-Normal distribution
    x ~ dnNormal(0,1) |> dnLog()      # Expressed using pipes.
    x ~ dnLognormal(0,1)              # This is equivalent.

    x ~ dnLog(dnGamma(2,3))           # There is no equivalent for this.
    x ~ dnIID(10,dnLog(dnGamma(2,3))) # Draw 10 log-Gamma(2,3) random variables.
## references
