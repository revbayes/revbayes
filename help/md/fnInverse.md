## name
fnInverse

## title
Inverse distribution 

## description
`fnInverse()` inverts a probability distribution.

`fnInverse(x).probability()` returns `1 / x.probability()`; 
`fnInverse(x).lnProbability()` returns `-x.lnProbability()`.

This provides a way to perform inference using conditional probabilities,
for example where 
Pr(x | Model, Condition) = Pr(x | Model) / Pr(x satisfies Condition)

## details
## authors
Martin R. Smith

## see_also
## example
x ~ dnNormal(mean = 0, sd = 1)
invX := fnInverse(x)
x.probability()
invX.probability()

## references
