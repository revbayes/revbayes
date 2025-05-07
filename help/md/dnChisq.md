## name
dnChisq
## title
Chi-Square Distribution
## description
The chi-square distribution with df degrees of freedom describes the distribution of the sum of the squares of df independent normal variables with mean 0 and variance 1. 
## details
The chi-square distribution takes one parameter, df, the number of degrees of freedom. The probability density function is f(x) = x^(df/2 - 1) * e^(-x/2) / (2^(df/2) * Gamma(df/2)), where Gamma is the gamma function.
## authors
Sebastian Hoehna
## see_also
## example
	# The most important use of the chi-square distribution
	# is arguable the quantile function.
	# You can access it the following way:
	df <- 10
	a := qchisq(0.025, df)
	a
	
## references
