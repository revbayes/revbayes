## name
mvGammaScale
## title
## description
A move to scale a single continuous value by multiplying by a value drawn from a Gamma(lambda,1) distribution. Lambda is the tuning parameter that controls the size of the proposals.
## details
## authors
Jeremy M. Brown
## see_also
mvScale
## example
	# Here is a simple example for conducting MCMC on the mean and sd of a Normal distribution.
	
	# Uniform(0,1) priors on the mean and sd
	mean ~ dnUnif(0,1)
	sd ~ dnUnif(0,1)
	
	# Dummy data (will not actually be analyzed)
	data <- v(0.4,0.5,0.6)
	
	# Clamping data
	for (i in 1:data.size()){ outcomes[i] ~ dnNorm(mean,sd); outcomes[i].clamp(data[i]) }
	
	# Initializing move and monitor counters
	mvi = 1
	mni = 1
	
	# Adding Gamma scale moves for the mean and sd (THIS MOVE IS HERE)
	moves[mvi++] = mvGammaScale(mean)
	moves[mvi++] = mvGammaScale(sd)
	
	# Instantiating the model
	mymodel = model(outcomes)
	
	# Adding screen monitor for the mean
	monitors[mni++] = mnScreen(mean, printgen=1000)
	
	# Creating MCMC object
	mymcmc = mcmc(mymodel, moves, monitors)
	
	# Running MCMC under the prior
	mymcmc.run(30000,underPrior=TRUE);
	
## references
