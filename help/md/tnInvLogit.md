## name
tnInvLogit
## title
InvLogit-transformed distribution
## description
InvLogit-transforms a given distribution.
## details
If X ~ dist then tnInvLogit(dist) is the distribution of exp(X)/(1+exp(X)).
The inverse logit function is also called the logistic function.

## authors
Ben Redelings
## see_also
logistic, tnExp, tnLog, tnLogit
## example
    p ~ tnInvLogit(dnNormal(0,1))      # The inverse-logit of a Normal random variable.
    p ~ dnNormal(0,1) |> tnInvLogit()  # Expressed using pipes.

    x ~ dnNormal(0,1)
    p := invlogit(x)                   # Expressed as a deterministic function of the log-odds.
## references
