## name
pdConstant

## title
Constant PseudoDataLikelihood

## description
Pseudodata that gives no information.

## details
This represents a pseudodata likelihood that is always 1.

See the help for dnPseudo for more explanation.

## authors
Benjamin Redelings

## see_also
dnPseudo

## example
f ~ dnPseudo( pdConstant() )
f.clamp( pseudoObservation() )

## references
