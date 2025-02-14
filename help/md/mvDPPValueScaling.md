## name
mvDPPValueScaling
## title
Scaling move applied to individual categories of a Dirichlet process mixture
## description
Operates on draws from a Dirichlet process prior (DPP) on mixtures of [RealPos](https://revbayes.github.io/documentation/RealPos.html) distributions, i.e., distributions defined over non-negative real numbers.
## details
In Dirichlet process mixtures, the number of categories (= clusters) is not specified beforehand but inferred from the data, and can range anywhere from 1 to the total number of elements (= observations). The move takes the current number of categories and simultaneously updates the value of every category using the scaling move with a scaling factor (lambda) of 1.
## authors
Sebastian Hoehna
## see_also
- dnDPP
- mvScale
- mvDPPValueBetaSimplex
- mvDPPValueSliding
## example
	# Here, we draw from a DP mixture for 3 elements, where every element
	# is a non-negative real number drawn from an exponential distribution
	x ~ dnDPP(1, dnExp(1), 3)
	
	# Next, we add the move. Note that without moves other than
	# mvDPPValueScaling, only the values of the categories will be updated:
        # the total number of categories and the assignment of elements to
	# categories will be determined by the initial draw.
	moves[1] = mvDPPValueScaling(x, weight=1)
	
	monitors[1] = mnScreen(x, printgen=1)
	mymodel = model(x)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.run(generations=50)
## references
