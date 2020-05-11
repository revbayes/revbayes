## name
treeTrace
## title
## description
Creates a tree trace object from a vector of trees.
## details
## authors
Will Freyman
## see_also
mvEmpiricalTree
treeTrace
readTreeTrace
readTrees
## example
	# Read in a vector of trees
	trees = readTrees("trees.nex")
	
	# Create a tree trace
	tree_trace = treeTrace(trees, burnin=0.25)
	
	# Create a distribution of trees from the tree trace
	tree ~ dnEmpiricalTree(tree_trace)
	
	# Add an MCMC move
	moves[1] = mvEmpiricalTree(tree)
	
## references
