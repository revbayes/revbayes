## name
mvRandomCategoryWalk
## title
Random walk on vector of naturals (e.g., from multinomial distribution) keeping
the sum of the vector constant.
## description
This random walk proposal picks a random index of a vector. Then, it picks
a random neighbor, either one index left or one index right. Then it decreases
the current value at the chosen index by one, and increases the value of the
neighbor by one.

Such a move is useful for values drawn from a multinomial distribution because
it keeps the total constant.
## details
## authors
Sebastian Höhna
## see_also
mvRandomGeometricWalk
dnMultinomial
## example
    # create a vector of moves
    moves = VectorMoves()
    
    # draw from a multinomial distribution
    num_e_prior <- simplex(rep(1, 6))
    number_events_pi ~ dnMultinomial(p=num_e_prior, size=20)

    # place a move on the draw
    moves.append( mvRandomCategoryWalk(x=number_events_pi, weight=5) )
## references
