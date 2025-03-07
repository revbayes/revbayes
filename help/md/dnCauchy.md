## name
dnCauchy
## title
Cauchy Distribution
## description
The Cauchy distribution describes the distribution of the ratio of two independent normal variables with mean 0 and variance 1. 
## details
The Cauchy distribution takes two parameters, location and scale. It is a symmetric distribution, but its tails are broad enough that it has no defined mean or variance. The probability density function is f(x) = 1/(pi * scale) * 1 / (1 + x - (location/scale)^2)
## authors
Andrew Magee
## see_also
dnNormal
dnChisq
## example
	# we simulate some obversations
	x <- rCauchy(n=10,location=0,scale=1)
	# let's see what the mean and the variance are.
	The mean will not converge with more samples, the Cauchy family has no moments.
	mean(x)
	var(x)
	sd(x)
	
## references
