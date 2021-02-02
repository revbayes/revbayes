## name
dnHalfNormal
## title
half-Normal Distribution
## description
half-normal (gaussian) distribution with offset equal to ‘offset’ and standard deviation equal to ‘sd’.
## details
The half-normal distribution has density:

 f(x) = 2/(sqrt(2 pi) sigma) e^-((x - offset)^2/(2 sigma^2)) where offset is the offset of the distribution and sigma the standard deviation.

f(x) = 2/(sqrt(2 pi) sigma) e^-((x - offset)^2/(2 sigma^2))

where offset is the offset of the distribution and sigma the standard deviation.

## authors
Andrew Magee
## see_also
dnNormal
dnLognormal
## example
	# we simulate some oversations
	x <- rhalfNormal(n=10,offset=0,sd=10)
	# let's see what the minimum is (you could do the max too)
	# the minimum should never be less than the offset
	min(x)
	# let's also see what the mean and the variance are
	mean(x)
	var(x)
	sd(x)
	
## references
