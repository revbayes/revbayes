## name
dnHalfCauchy
## title
half-Cauchy Distribution
## description
Half-Cauchy distribution with location equal to ‘location’ and scale equal to ‘scale’.
## details
The half-Cauchy distribution has density:

f(x) = 2/(pi * sigma) * 1/(1 + x-(location/scale)^2)

## authors
Andrew Magee
## see_also
## example
	# we simulate some obversations
	x <- rHalfCauchy(n=10,location=0,scale=1)
	# let's see what the minimum is (you could do the max too). If this is not ‘location’, we're in trouble!
	min(x)
	# let's also see what the mean and the variance are.
	The mean will not converge with more samples, the Cauchy family has no moments.
	mean(x)
	var(x)
	sd(x)
	
## references
