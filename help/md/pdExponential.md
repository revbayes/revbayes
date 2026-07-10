## name
pdExponential

## title
Exponential PseudoDataLikelihood

## description
A pseudodata likelihood that is shaped like an Exponential

## details
This represents a pseudodata likelihood that is given by

* exp( -lambda * (x-shift) ) if x > shift
* 0 if x <= shift

Note the pseudodata likelihood is NOT a density, and is not
normalized to sum to 1. Instead it is normalized so that the
maximum value is 1.

This is why the additional lambda term is not present.

See dnPseudo for more information.

## authors
Benjamin Redelings

## see_also
dnPseudo
pdBetween
pdLogNormal

## example

# Evidence that x above 1
x ~ dnNormal(0,10)
f ~ dnPseudo( pdExponential(x,2,1) )
f.clamp( pseudoObservation() )

## references
