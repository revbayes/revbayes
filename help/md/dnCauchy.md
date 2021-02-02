## name
dnCauchy
## title
Cauchy Distribution
## description
Cauchy distribution with location equal to ‘location’ and scale equal to ‘scale’.
## details
The Cauchy distribution has density:

f(x) = 1/(pi * scale) * 1/(1 + x-(location/scale)^2)

## authors
Andrew Magee
## see_also
## example
	# we simulate some obversations
	x <- rCauchy(n=10,location=0,scale=1)
	# let's see what the mean and the variance are.
	The mean will not converge with more samples, the Cauchy family has no moments.
	mean(x)
	var(x)
	sd(x)
	
## references
