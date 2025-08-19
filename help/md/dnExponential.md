## name
dnExponential
## title
Exponential Distribution
## description
The Exponential distribution describes the distribution of the times between events in a Poisson point process. 
## details
The exponential distribution takes one parameter, lambda, describing the rate (i.e. 1/mean). The probability density function is f(x) = lambda * exp(-lambda*x).

## authors
Michael Landis
## see_also
## example
	# we set a rate parameter
	rate <- 10.0
        # we create an exponentially distributed random variable
	x ~ dnExponential(lambda=rate)
	# compute the probability of the variable
	x.probability()
	
## references
