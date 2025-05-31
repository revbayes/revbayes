## name
pdLogNormal

## title
LogNormal PseudoDataLikelihood

## description
A pseudodata likelihood that is shaped like a LogNormal.

## details
This represents a pseudodata likelihood that is given by

* exp( -(log(x-shift)-mean)^2 / (2*sd^2) ) if x > shift
* 0 if x <= shift

Note the pseudodata likelihood is NOT a density, and is not
normalized to sum to 1. Instead it is normalized so that the
maximum value is 1.

This is why the term (1/x)*(1/sqrt(2*pi*sigma^2)) is not present.

See dnPseudo for more information.

## authors
Benjamin Redelings

## see_also
dnPseudo
pdBetween
pdNormal

## example

# Evidence that x is near 2
x ~ dnNormal(0,10)
f ~ dnPseudo( pdLogNormal(x,2,1) )
f.clamp( pseudoObservation() )

## references
