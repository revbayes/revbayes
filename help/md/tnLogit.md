## name
tnLogit
## title
Logit-transformed distribution
## description
Logit-transforms a given distribution.
## details
If P ~ dist then tnLogit(dist) is the distribution of log(P/(1-P)).

## authors
Ben Redelings
## see_also
logit, tnExp, tnLog, tnInvLogitit
## example
    x ~ tnLogit(dnBeta(1,2))         # The log-odds of an Beta random variable.
    x ~ dnBeta(1,2)|> tnLogit()      # Expressed using pipes.

    p ~ dnBeta(1,2)
    x := logit(p)                    # Expressed as a deterministic function of the probability.
## references
