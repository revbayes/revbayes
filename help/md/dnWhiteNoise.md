## name
dnWhiteNoise
## title
White-Noise Process
## description
White-Noise process for positive real numbers.
## details
The white-noise process is a process of a positive continuous variable similar to Brownian motion and the Ornstein-Uhlenbeck process. However, the white-noise process has a large variance when the time is small, and has small variance if the time is large.
## authors
Sebastian Hoehna
## see_also
dnGamma
## example
	# lets simulate
	a <- rWhiteNoise(1000,mu=1,sigma=4,time=4)
	# we expect a mean of 1
	mean(a)
	
	# create a random variable
	x ~ dnWhiteNoise(mu=1.0,sigma=4,time=1)
	x
	
## references
