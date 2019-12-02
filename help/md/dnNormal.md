## name
dnNormal
## title
Normal Distribution
## description
Normal (gaussian) distribution with mean equal to ‘mean’ and standard deviation equal to ‘sd’.
## details
The normal distribution has density:

f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))

where mu is the mean of the distribution and sigma the standard deviation.

## authors
Sebastian Hoehna
## see_also
dnLognormal
## example
	# we simulate some observations
	x <- rnorm(n=10,mean=5,sd=10)
	# let's see what the minimum is (you could do the max too)
	min(x)
	# let's also see what the mean and the variance are
	mean(x)
	var(x)
	sd(x)
	
## references
