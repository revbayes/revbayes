## name
mvCorrelationMatrixElementSwap
## title
Correlation Matrix element swap move.
## description
Swaps elements i and j of the correlation matrix (i != j).
## details
## authors
Michael R. May
## see_also
## example
	
	# we draw a correlation matrix from an LKJ distribution
	R ~ dnLKJ(eta=1, dim=5)
	
	# we specify an element swap move
	moves[1] = mvCorrelationMatrixElementSwap(R)
	
## references
