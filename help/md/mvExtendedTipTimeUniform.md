## name
mvExtendedTipTimeUniform

## title
Extended tip extinction time move

## description
Draws a new extinction time for a random extinct tip of an extended tree, uniformly between the present and the taxon's youngest occurrence.

## details
The tips of an extended tree are extinction events rather than occurrences, so a tip age may fall below its taxon's fossil age range. This move therefore draws on the data-fixed window between the present and the youngest occurrence, rather than within the age range as `mvFossilTipTimeUniform` does. The proposal is symmetric and configurations that place the extinction above the augmented oldest age are rejected by the distribution.

Extant tips are pinned at the present and are never proposed.

## authors
June Walker

## see_also
dnFossilizedBirthDeathSpeciation
mvResampleAugmentedAges

## example
tr ~ dnFBDSP(origin=origin, lambda=lambda, mu=mu, psi=psi, rho=1, timeline=timeline, taxa=taxa)
moves.append( mvExtendedTipTimeUniform(tr, weight=taxa.size()) )

## references
