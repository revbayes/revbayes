## name
dnExponential
## title
Exponential Distribution
## description
Exponential distribution with rate equal to ‘lambda’.
## details
The exponential distribution has density:

f(x) = lambda * exp(-lambda*x)

where lambda is the rate parameter.

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
