## name
dnEmpiricalTree
## title
## description
Creates a distribution of trees from a trace of trees.
## details
## authors
Will Freyman
Sebastian Hoehna
Will Pett
## see_also
mvEmpiricalTree
treeTrace
readTreeTrace
## example
	# Read in tree trace
	tree_trace = readTreeTrace("output/my.trees", burnin=0.25)
	
	# Create a distribution of trees
	tree ~ dnEmpiricalTree(tree_trace)
	
	# Add an MCMC move
	moves[1] = mvEmpiricalTree(tree)
	
## references
