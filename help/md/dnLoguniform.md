## name
dnLoguniform
## title
Log-Uniform Distribution
## description
A strictly positive real number x has a log-uniform distribution over interval (min,max) if its logarithm y = ln(x) has uniform distribution over interval (ln(min),ln(max)).
## details
The log-uniform distribution is defined over strictly positive real numbers. Saying that x is log-uniform is equivalent to saying that y = ln(x) is uniform. The log-uniform distribution therefore expresses lack of information about the order of magnitude of a scale parameter:  if x has a log-uniform distribution, then it has equal chance to be contained by any of the intervals of the form (10^k, 10^(k+1)) within the allowed range.

The density is p(x) = 1/x, which can be seen by defining x = exp(y) where y has uniform distribution and apply the change-of-variable formula.

The log-uniform distribution is improper when defined over the entire positive real line. To always make it proper, in RevBayes, a min and a max should always be specified.

## authors
Nicolas Lartillot
## see_also
dnUniform
## example
	# a log-uniform prior over the rate of change of a Brownian trait (or a Brownian relaxed clock)
	trueTree = readTrees("data/primates.tree")[1]
	sigma ~ dnLogUniform(min=0.001, max=1000)
	X ~ dnBrownian(trueTree,sigma)
	# ...
	
## references
