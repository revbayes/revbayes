## name
tnLogit
## title
Logit-transformed distribution
## description
Logit-transforms a given distribution.
## details
If P ~ dist then tnLogit(dist) is the distribution of log(P/(1-P)).

The distribution `dist` can be either univariate (dnBeta) or
multivariate (dnDirichlet).

## authors
Ben Redelings
## see_also
logit, tnExp, tnLog, tnInvlogitit
## example
    x ~ tnLogit(dnBeta(1,2))         # The log-odds of an Beta random variable.
    x ~ dnBeta(1,2)|> tnLogit()      # Expressed using pipes.

    p ~ dnBeta(1,2)
    x := logit(p)                    # Expressed as a deterministic function of the probability.

    xs ~ dnDirichlet([1,1,1,1]) |> tnLogit() 
## references
