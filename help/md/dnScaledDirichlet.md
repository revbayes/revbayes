## name
dnScaledDirichlet
## title
Scaled Dirichlet Distribution
## description
Scaled Dirichlet probability distribution on a simplex.
## details
The scaled Dirichlet probability distribution is the generalization of the dirichlet distribution. A random variable from a scaled Dirichlet distribution is a simplex, i.e., a vector of probabilities that sum to 1. If b[1]=b[2]=...=b[n], then the scaledDirichlet(alpha,beta) collapses to the Dirichlet with the same alphas.
## authors
Andrew Magee
## see_also
dnDirichlet
simplex
## example
	# lets get a draw from a Dirichlet distribution
	a <- [1,1,1,1]   # we could also use rep(1,4)
	b <- [1,2,3,4]   # if these are all equal, the scaled Dirichlet is equivilent to the Dirichlet(a)x ~ dnScaledDirichlet(a,b)
	x
	# let check if b really sums to 1
	sum(x)
	
## references
