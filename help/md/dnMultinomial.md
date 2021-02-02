## name
dnMultinomial
## title
Multinomial Distribution
## description
A multinomial distribution defines a probability distribution on a vector of natural numbers. It is understood as randomly picking n times from the k categories with replacement where each catefory has its own probability p[i].
## details
## authors
Sebastian Hoehna
## see_also
dnDirichlet
## example
	p <- simplex(1,1,1,1)
	x ~ dnMultinomial(10, p)
	y ~ dnDirichlet(x)
	y.clamp( simplex(1,2,3,4) )
	moves[1] = mvSlide(x, delta=0.1, weight=1.0)
	monitors[1] = screenmonitor(printgen=1000, separator = "	", x)
	mymodel = model(p)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references
