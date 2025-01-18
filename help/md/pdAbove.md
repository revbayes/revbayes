## name
pdAbove

## title
PseudoDataLikelihood Above

## description
Pseudodata that x >= lower

## details
This represents a pseudodata likelihood that is
* 1 if x in [lower, infinity)
* 0 otherwise

Note the pseudodata likelihood is NOT a density, and is not
normalized to sum to 1.  Instead it is normalized so that the
maximum value is 1.

Specifying decayRate allows a "soft bound".  In this case
if x is below lower, the pseudodata likelihood is:
* exp(-D*decay)
where D = x - lower.

See dnPseudo for more information.

## authors
Benjamin Redelings

## see_also
dnPseudo
pdBelow
pdBetween
pdExponential

## example

# Evidence that x > 1
x ~ dnNormal(0,1)
f1 ~ dnPseudo( pdAbove(x,1) )
f1.clamp( pseudoObservation() )

# Evidence that t1 > t2
t1 ~ dnNormal(0,1)
t2 ~ dnNormal(0,1)
f2 ~ dnPseudo( pdAbove(t1,t2) )
f2.clamp( pseudoObservation() )

