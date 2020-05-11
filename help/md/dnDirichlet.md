## name
dnDirichlet
## title
Dirichlet Distribution
## description
Dirichlet probability distribution on a simplex.
## details
The Dirichlet probability distribution is the generalization of the Beta distribution. A random variable from a Dirichlet distribution is a simplex, i.e., a vector of probabilities that sum to 1.
## authors
Sebastian Hoehna
## see_also
simplex
## example
	# lets get a draw from a Dirichlet distribution
	a <- [1,1,1,1]   # we could also use rep(1,4)
	b ~ dnDirichlet(a)
	b
	# let check if b really sums to 1
	sum(b)
	
## references
