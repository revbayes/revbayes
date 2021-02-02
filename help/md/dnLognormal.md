## name
dnLognormal
## title
Lognormal Distribution
## description
Lognormal distribution is the distribution for a log-transformed normally distributed random variable with mean 'mu' and standard deviation 'sigma'.
## details
The lognormal random variable is defined as

:X = exp(mu + sigma Z)

where mu is the mean parameter, sigma is the standard deviation, and Z is a standard normal random variable. Note, in effect, the mean and standard deviation provide the location and scale of the exponentiated normal variate, mu + sigma Z.The lognormal distribution has density:

f(x) = 1/(x sigma sqrt(2 pi)) e^-((ln x - mu)^2/(2 sigma^2))

where mu is the mean of the distribution and sigma the standard deviation.

## authors
Michael Landis
## see_also
## example
	# set an expected value for x
	expectation_of_x <- 1
	# set a mean and sd parameter
	sd <- 0.5
	mean <- ln(expectation_of_x) - 0.5 * sd^2
	# create a lognormal distribution with expected value of 1
	x ~ dnLognormal(mean=mean, sd=sd)
	
## references
