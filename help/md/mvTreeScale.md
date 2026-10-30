## name
mvTreeScale
## title
Scaling Move for Time Tree Node Ages
## description
Scales internal ages in a `TimeTree` by the same factor while leaving the topology unchanged.
Optionally scales an independent origin age with the tree.
## details
The `mvTreeScale` move scales internal ages by a random factor of exp(delta * (u - 0.5)),
where delta is a tuning parameter and u is a random draw from the uniform distribution on [0, 1].
Scaling stops at tips and sampled-ancestor parents, leaving their entire subtrees unchanged.
The Hastings factor counts only the internal ages being scaled.

Supply `rootAge` when the tree's root age is represented by a separate stochastic variable.
Its value must equal the tree's root age; it tracks that same coordinate during the move.
For an origin-conditioned tree, use `originAge` instead to scale the independent origin
from its own current value. The origin contributes one additional coordinate to the Hastings factor.
Supplying both arguments is an error. With neither argument, only the tree is scaled.
## authors
## see_also
mvSubtreeScale
mvNodeTimeSlide
## example
    taxa <- v(taxon("A"), taxon("B"), taxon("C"), taxon("D"), taxon("E"), taxon("F"))
    height ~ dnUniform(0, 10)
    moves = VectorMoves()
    
    # Simulate a simple TimeTree
    tree ~ dnBDP(lambda=1.0, mu=0.2, rootAge=height, taxa=taxa)
    
    # Assign it a mvTreeScale move
    moves.append( mvTreeScale(tree=tree, rootAge=height, delta=1, tune=TRUE, weight=5) )

    # Jointly scale an independent origin and an origin-conditioned tree
    origin ~ dnExponential(0.1)
    fossil_tree ~ dnFBDP(origin=origin, lambda=1.0, mu=0.2, psi=0.1, rho=1.0, taxa=taxa)
    moves.append( mvTreeScale(tree=fossil_tree, originAge=origin, weight=5) )

## references
- citation: Yang Z (2014). Molecular Evolution: A Statistical Approach. Oxford, UK: Oxford University Press.
  doi: 10.1093/acprof:oso/9780199602605.001.0001
  url: https://academic.oup.com/book/26340
