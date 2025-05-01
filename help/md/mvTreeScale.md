## name
mvTreeScale
## title
Scaling Move for Time Tree Node Ages
## description
Scales the ages of all internal nodes in a `TimeTree` by the same factor while
leaving the topology unchanged.
## details
The `mvTreeScale` move scales the ages of all internal nodes (including
the root) by a random factor of exp(delta * (u - 0.5)), where delta is a tuning
parameter and u is a random draw from the uniform distribution on [0, 1].
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

## references
- citation: Yang Z (2014). Molecular Evolution: A Statistical Approach. Oxford, UK: Oxford University Press.
  doi: 10.1093/acprof:oso/9780199602605.001.0001
  url: https://academic.oup.com/book/26340
