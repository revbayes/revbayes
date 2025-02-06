## name
dnDirichlet
## title
Dirichlet Distribution
## description
The Dirichlet probability distribution is a generalization of the Beta distribution for multiple variables.
## details
The Dirichlet distribution takes one parameter, alpha, a vector of numbers representing the concentration of the distribution on each variable. It then returns a simplex (i.e. a vector whose elements sum to 1) representing the relative probability of each variable. 
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
