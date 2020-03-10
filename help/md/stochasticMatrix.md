## name
stochasticMatrix
## title
Building a stochastic matrix.
## description
A stochastic matrix is a matrix (not necessarily square) with rows that sum to 1.
## details
## authors
Michael R. May
## see_also
## example
	vec[1] ~ dnDirichlet( [1,1,1,1] )
	vec[2] ~ dnDirichlet( [1,1,1,1] )
	vec[3] ~ dnDirichlet( [1,1,1,1] )
	vec[4] ~ dnDirichlet( [1,1,1,1] )

	m := stochasticMatrix(vec)	

## references
