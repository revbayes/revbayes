## name
mccTree
## title
## description
Finds the maximum clade credibility (MCC) topology from a trace of trees and summarizes branch lengths.
## details
## authors
Will Pett
Sebastian Hoehna
## see_also
consensusTree
mapTree
treeTrace
readTreeTrace
## example
	# Read in tree trace
	tree_trace = readTreeTrace("output/my.trees", burnin=0.25)
	
	# Generate the MCC tree
	map_tree = mccTree(trace=tree_trace, file="mcc.tree")
	
## references
