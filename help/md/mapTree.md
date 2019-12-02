## name
mapTree
## title
## description
Finds the maximum a posteriori (MAP) topology from a trace of trees and summarizes branch lengths.
## details
## authors
Will Freyman
Sebastian Hoehna
## see_also
consensusTree
mccTree
treeTrace
readTreeTrace
## example
	# Read in tree trace
	tree_trace = readTreeTrace("output/my.trees", burnin=0.25)
	
	# Generate the MAP tree
	map_tree = mapTree(trace=tree_trace, file="map.tree")
	
## references
