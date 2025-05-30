## name
pseudoObservation

## title
A pseudo-observation

## description
An observation that is unspecified

## details
The function pseudoObservation() represents the observed
but unspecified value of a variable sampled from dnPseudo.

For example, if a fossil observation F=f provides evidence
that a taxon age t is between age_min and age_max, we could
write:

    F ~ dnPseudo( pdAbove(t, age_min, age_max) )
    F.clamp( pseudoObservation() )

Here pseudoObservation() represents the unspecified value f,
similar to how dnPseudo represents the unspecified distribution
D.

See the help for dnPseudo for a more complete explanation.

## authors
Benjamin Redelings

## see_also
dnPseudo
pdRequire

## example
x ~ dnNormal(0,1)
f ~ dnPseudo( pdAbove(x,1) )
f.clamp( pseudoObservation() )

## references
