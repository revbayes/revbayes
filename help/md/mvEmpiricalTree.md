## name
mvEmpiricalTree
## title
Move on an empirical tree distribution
## description
An MCMC move that operates on empirical tree distributions.
## details
## authors
Will Freyman
Sebastian Hoehna
Will Pett
Jiansi Gao
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
