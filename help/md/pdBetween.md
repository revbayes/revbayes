## name
pdBetween

## title
PseudoDataLikelihood Between

## description
PseudoData that lower <= x <= upper

## details
This represents a pseudodata likelihood that is
* 1 if x is in [lower, upper]
* 0 otherwise

Note the pseudodata likelihood is NOT a density, and is not
normalized to sum to 1.  Instead it is normalized so that the
maximum value is 1.

Specifying decayRate allows a "soft bound".  In this case
if x is not in the interval then the pseudodata likelihood is:
* exp(-D*decay)
where D is the  distance from x to the interval.

See dnPseudo for more information.

## authors
Benjamin Redelings

## see_also
dnPseudo
pdAbove
pdBelow
pdLogNormal

## example

# Evidence that x is in [1,2]
x ~ dnNormal(0,1)
f ~ dnPseudo( pdBetween(x,1,2) )
f.clamp( pseudoObservation() )

