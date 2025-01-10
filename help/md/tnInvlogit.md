## name
tnInvlogit
## title
Invlogit-transformed distribution
## description
Invlogit-transforms a given distribution.
## details
If X ~ dist then tnInvlogit(dist) is the distribution of exp(X)/(1+exp(X)).
The inverse logit function is also called the logistic function.

The distribution `dist` can be either univariate (dnNormal) or
multivariate (MultivariateNormal).

## authors
Ben Redelings
## see_also
logistic, tnExp, tnLog, tnLogit
## example
    p ~ tnInvlogit(dnNormal(0,1))      # The inverse-logit of a Normal random variable.
    p ~ dnNormal(0,1) |> tnInvlogit()  # Expressed using pipes.

    x ~ dnNormal(0,1)
    p := invlogit(x)                   # Expressed as a deterministic function of the log-odds.

    ps ~ dnIID(4,dnNormal(0,1)) |> tnInvlogit()

    mu = [1.0, 2.0, 3.0, 4.0]
    Sigma ~ dnWishart(df=4, kappa=2, dim=4)
    x ~ dnMultivariateNormal(mu,Sigma) |> tnInvlogit()
## references
