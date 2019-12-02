## name
dnLogExponential
## title
Log-Exponential Distribution
## description
A real number x has a log-Exponential distribution if y = exp(x) has Exponential distribution.
## details
The log-Exponential distribution is defined over real numbers. Saying that x is log-Exponential is equivalent to saying that y = exp(x) is Exponential. The log-Exponential distribution therefore expresses lack of information about the order of magnitude of a scale parameter:  if x has a log-Exponential distribution, then it has equal chance to be contained by any of the intervals of the form (10^k, 10^(k+1)) within the allowed range.

The density is p(x) = ???, which can be seen by defining x = ln(y) where y has Exponential distribution and apply the change-of-variable formula.

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
