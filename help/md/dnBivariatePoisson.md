## name
dnBivariatePoisson
## title
Bivariate Poisson Distribution
## description
A Bivariate Poisson distribution defines probabilities for pairs of natural numbers.
## details
## authors
Alexander Zarebski
Sebastian Hoehna
## see_also
dnPoisson
## example
	th1 ~ dnUniform(0.0,10.0)
	th2 ~ dnUniform(0.0,10.0)
	th0 ~ dnUniform(0.0,10.0)
	 ~ dnBivariatePoisson(th1, th2, th0)
	x.clamp([3, 3, 3])
	oves[1] = mvSlide(th1, delta=0.01, weight=1.0)
	moves[2] = mvSlide(th2, delta=0.01, weight=1.0)
	oves[3] = mvSlide(th0, delta=0.01, weight=1.0)
	monitors[1] = mnScreen(printgen=1000,  th0)
	ymodel = model(th1)
	mymcmc = mcmc(mymodel, monitors, moves)
	ymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references
- citation: Karlis D, Ntzoufras J (2003). Bayesian and Non-Bayesian Analysis of Soccer
    Data using Bivariate Poisson Regression Models.
  doi: null
  url: null
