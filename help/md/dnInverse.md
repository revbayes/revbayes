## name
dnInverse

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
```
# Draw X from an exponential distribution
x ~ dnExp(lambda = 1)
x.clamp(42)
x.lnProbability()

# Now calculate the inverse
invX ~ dnInverse(dnExp(lambda = 1))
# Clamp the value of the draw from the exponential distribution
invX.clamp(42)
invX.lnProbability()
```

## references
