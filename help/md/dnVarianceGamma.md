## name
dnVarianceGamma
## title
Variance-gamma Distribution
## description
Variance-gamma distribution with location ‘mu’.
## details
The variance-gamma distribution has density:

f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))

where mu is the mean of the distribution and sigma the standard deviation.

## authors
Michael Landis
## see_also
## example
	# we simulate some oversations
	x <- rnorm(n=10,mean=5,sd=10)
	# let's see what the minum is (you could do the max too)
	min(x)
	# let's also see what the mean and the variance are
	mean(x)
	var(x)
	sd(x)
	
## references
