## name
dnCategorical
## title
The Categorical Distribution
## description
The Categorical distribution generalizes the Bernoulli distribution, describing the probability of choosing from a number of outcomes, each with their own probability. 
## details
The categorical distribution takes on a parameter p, a simplex (i.e. vector, the elements of which sum to 1). It returns outcome i with probability p[i].

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
