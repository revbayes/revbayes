## name
dnGamma
## title
Gamma Distribution
## description
The Gamma distribution describes the probability of the sum of exponentially distributed variables. 
## details
The gamma distribution takes two parameters, shape and rate. Similar to how 1/rate gives the mean of the exponential, shape/rate gives the mean of the gamma. It provides a natural prior distribution for parameters that could be considered as sums of exponential variables.

The probability density function is f(x) = rate^shape * x^(shape - 1) * e^(-rate * x) / Gamma(shape), where Gamma is the gamma function. Note that for shape = 1, the gamma distribution reduces to an exponential distribution.
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
