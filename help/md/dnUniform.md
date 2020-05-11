## name
dnUniform
## title
Uniform Distribution
## description
Uniform distribution with lower and uppper bounds.
## details
The uniform distribution has density, f(x) = 1/(b-a), where b is the upper bound and a is the lower bound.
## authors
Michael Landis
## see_also
## example
	# set the lower bound
	a <- -2.5
	# set the upper bound
	b <- -3.9
	# create a stochastic node with a uniform prior
	u ~ dnUniform(a, b)
	
## references
