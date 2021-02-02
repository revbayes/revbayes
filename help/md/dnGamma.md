## name
dnGamma
## title
Gamma Distribution
## description
Gamma probability distribution for positive real numbers.
## details
The gamma distribution is the probability of the sum of exponentially distributed variables. Thus, it provides a natural prior distribution for parameters that could be considered as sums of exponential variables.
## authors
Sebastian Hoehna
## see_also
dnExponential
## example
	# lets simulate
	a <- rgamma(1000,shape=4,rate=4)
	# we expect a mean of 1
	mean(a)
	
	# create a random variable
	x ~ dnGamma(shape=4,rate=1)
	x
	
## references
