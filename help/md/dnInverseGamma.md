## name
dnInverseGamma
## title
inverseGamma Distribution
## description
inverse-gamma probability distribution for positive real numbers.
## details
The inverse Gamma distribution is the probability of the sum of exponentially distributed variables. Thus, it provides a natural prior distribution for parameters that could be considered as sums of exponential variables.
## authors
Sebastian Hoehna
## see_also
dnExponential
## example
	# lets simulate
	a <- rinverseGamma(1000,shape=4,rate=4)
	# we expect a mean of 1
	mean(a)
	
	# create a random variable
	x ~ dnInverseGamma(shape=4,rate=1)
	x
	
## references
