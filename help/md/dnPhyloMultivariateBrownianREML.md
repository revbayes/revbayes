## name
dnPhyloMultivariateBrownianREML
## title
Phylogenetic Multivariate Brownian Motion
## description
Multivariate Brownian motion over a phylogeny with variance-covariance matrix rateMatrix.
## details
## authors
Michael R. May
Nicolai Vetr
## see_also
dnPhyloBrownianREML
dnPhyloBrownianMVN
## example
	
	# generate a tree and variance-covariance matrix.
	psi ~ dnUniformTimeTree(1, [taxon("A"),taxon("B"),taxon("C")])
	Sigma <- diagonalMatrix(5)
	
	# generate the multivariate data.
	x ~ dnPhyloMultivariateBrownianREML(tree=psi, rateMatrix=Sigma)
	
	# print the simulated data.
	x
	
## references
- citation: Huelsenbeck JP, B Rannala (2003). Detecting correlation between characters
    in a comparative analysis with uncertain phylogeny. Evolution, 57(6):1237-1247.
  doi: null
  url: null
