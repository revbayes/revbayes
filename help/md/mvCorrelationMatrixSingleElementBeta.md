## name
mvCorrelationMatrixSingleElementBeta
## title
Correlation Matrix Beta proposal.
## description
Beta proposal on a random element of a correlation matrix.
## details
This move chooses a single element of the correlation matrix at random, and draws a proposed value from a Beta distribution centered on the current value (and stretched to range from -1 to 1).
## authors
Michael R. May
## see_also
mvCorrelationMatrixSpecificElementBeta
mvCorrelationMatrixRandomWalk
## example
	
	# we draw a correlation matrix from an LKJ distribution
	R ~ dnLKJ(eta=1, dim=5)
	
	# we specify a beta move on the correlation matrix
	moves[1] = mvCorrelationMatrixSingleElementBeta(R, alpha=10.0)
	
## references
