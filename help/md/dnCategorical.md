## name
dnCategorical
## title
The Categorical Distribution
## description
The categorical distribution, sometimes referred to as the generalized Bernoulli distribution. It describes the probability of one of K different outcomes, labeled from 1 to K, with each outcome probability separately specified.
## details
The argument to the constructor is a simplex containing the probabilities of the outcomes. The outcomes are labeled from 1 to K, where K is the number of elements in the simplex. Outcome i has probability specified by component i in the simplex.

A typical scenario where a categorical variable is used is in the definition of a variable drawn from a mixture. A vector of mixture components is set up first, and then a stochastic variable drawn from a categorical distribution is used as an index in a deterministic assignment that points to a component in the mixture. See example below.

## authors
Fredrik Ronquist
## see_also
dnBinomial
## example
	# Define a stochastic variable x that is drawn from
	# a categorical distribution with 4 categories, each
	# category having the same probability, then examine
	# the value of x.
	x ~ dnCat( simplex(1,1,1,1) )
	x
	
	# Draw 10 values from the distribution and place them
	# in a vector a, then examine a.
	for ( i in 1:10 ) {
	    a[i] <- x
	    x.redraw()
	}
	a
	
	# Use x in defining a deterministic variable y taking
	# on values from a mixture of RealPos values representing
	# rates from a discretized scaled gamma distribution
	# with four categories.
	shape ~ dnExp( 10.0 )
	rates := fnDiscretizeGamma( shape, shape, 4 )
	y := rates[x]
	
## references
