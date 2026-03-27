## name
mvRandomCategoryWalk
## title
Random walk on vector of naturals (e.g., from multinomial distribution) keeping the sum of the vector constant.
## description
This random walk proposal picks a random index of a vector. Then, it picks a random neighbor, either one index left or one index right. Then it decreases the current value at the chosen index by one, and increases the value of the neighbor by one.

Such a move is important for a value from a multinomial distribution because it keeps the total equal.
## details
## authors
Sebastian Höhna
## see_also
mvRandomGeometricWalk
dnMultinomial
## example
num_e_prior <- simplex(rep(1, 6))
number_events_pi ~ dnMultinomial(p=num_e_prior, size=20)

moves.append( mvRandomCategoryWalk(x=number_events_pi, weight=5) )
## references
