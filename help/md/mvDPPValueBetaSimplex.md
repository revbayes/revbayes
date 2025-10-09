## name
mvDPPValueBetaSimplex
## title
Beta simplex move applied to individual categories of a Dirichlet process mixture
## description
Operates on draws from a Dirichlet process prior (DPP) on mixtures of [Simplex](https://revbayes.github.io/documentation/Simplex.html) distributions, i.e., distributions defined over vectors whose elements are positive and sum to 1.
## details
In Dirichlet process mixtures, the number of categories (= clusters) is not specified beforehand but inferred from the data, and can range anywhere from 1 to the total number of elements (= observations). The move takes the current number of categories and simultaneously updates the value of every category using the beta simplex move with a concentration parameter (alpha) of 10.
## authors
Sebastian Hoehna
## see_also
- dnDPP
- mvBetaSimplex
- mvDPPValueScaling
- mvDPPValueSliding
## example
	# Here, we draw from a DP mixture for 3 elements, where every element
	# is itself a 2-element simplex drawn from a flat Dirichlet distribution
	x ~ dnDPP(1, dnDirichlet( [1, 1] ), 3)
	
	# Next, we add the move. Note that without moves other than
	# mvDPPValueBetaSimplex, only the values of the categories will be
        # updated: the total number of categories and the assignment of elements
	# to categories will be determined by the initial draw.
	moves[1] = mvDPPValueBetaSimplex(x, weight=1)
	
	monitors[1] = mnScreen(x, printgen=1)
	mymodel = model(x)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.run(generations=50)
## references
