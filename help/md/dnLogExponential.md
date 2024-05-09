## name
dnLogExponential
## title
Log-Exponential Distribution
## description
A real number x has a log-Exponential distribution if y = exp(x) has Exponential distribution.
## details
The log-Exponential distribution is defined over positive real numbers. Saying that x is log-Exponential is equivalent to saying that log(x) is Exponential.

The density is p(x) = lambda*exp(-lambda*log(x))/x = lambda * x**(lambda-1).

## authors
Sebastian Hoehna
## see_also
dnExponential
## example
	# a log-Exponential prior over the rate of change of a Brownian trait (or a Brownian relaxed clock)
	trueTree = readTrees("data/primates.tree")[1]
	log_sigma ~ dnLogExponential(lambda=1)
	sigma := exp(log_sigma)
	X ~ dnBrownian(trueTree,sigma)
	# ...
	
## references
